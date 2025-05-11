#ifndef MSAGPU_HPP

#include <hip/hip_runtime.h>
#include "msa-gpu.hip.hpp"
#endif

void getGpuInfo(po::variables_map &vm, msa::option *option) {
    int maxGpuNum;
    hipGetDeviceCount(&maxGpuNum);
    int gpuNum = (vm.count("gpu")) ? vm["gpu"].as<int>() : maxGpuNum;
    if (gpuNum < 0)
    {
        std::cerr << "ERROR: Requested number of GPU <= 0.\n";
        exit(1);
    }
    if (gpuNum > maxGpuNum)
    {
        std::cerr << "ERROR: Requested number of GPU more than available GPUs.\n";
        exit(1);
    }
    if (option->cpuOnly)
        gpuNum = 0;
    else if (gpuNum > option->cpuNum)
    {
        if (option->cpuNum == 1)
            std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << option->cpuNum << " GPU.\n";
        else
            std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << option->cpuNum << " GPUs.\n";
        gpuNum = option->cpuNum;
    }
    std::vector<int> gpuIdx;
    if (vm.count("gpu-index"))
    {
        std::string gpuIdxString = vm["gpu-index"].as<std::string>();
        std::string id = "";
        for (int i = 0; i < gpuIdxString.size(); ++i)
        {
            if (gpuIdxString[i] == ',')
            {
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else if (i == gpuIdxString.size() - 1)
            {
                id += gpuIdxString[i];
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else
                id += gpuIdxString[i];
        }
    }
    else
    {
        for (int i = 0; i < gpuNum; ++i)
            gpuIdx.push_back(i);
    }
    if (gpuIdx.size() != gpuNum)
    {
        std::cerr << "ERROR: the number of requested GPUs does not match the number of specified gpu indexes.\n";
        exit(1);
    }
    for (auto id : gpuIdx)
    {
        if (id >= maxGpuNum)
        {
            std::cerr << "ERROR: specified gpu index >= the number of GPUs\n";
            exit(1);
        }
    }
    printf("Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);
    if (gpuNum == 0)
    {
        std::cout << "WARNING: Requested 0 GPU. CPU version is used.\n";
        option->cpuOnly = true;
    }

    option->gpuNum = gpuNum;
    option->gpuIdx = gpuIdx;
    return;
}

void msaOnSubtreeGpu(Tree *T, msa::utility *util, msa::option *option, partitionInfo_t *partition, Params &param)
{

    auto progressiveStart = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Node *, Node *>>> hier;
    for (auto &p : partition->partitionsRoot)
    {
        std::stack<Node *> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<std::pair<std::pair<Node *, Node *>, int>> subhier;
        int grpID = p.second.first->grpID;
        getMsaHierachy(subhier, msaStack, grpID);
        for (auto h : subhier)
        {
            while (hier.size() < h.second + 1)
            {
                std::vector<std::pair<Node *, Node *>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
    }
    auto scheduleEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scheduleTime = scheduleEnd - progressiveStart;
    if (option->printDetail)
        std::cout << "Scheduling in " << scheduleTime.count() / 1000 << " us\n";

    std::unordered_map<std::string, std::string> beforeAln;
    int level = 0;
    int cpuThres = option->cpuNum * 3;
    if (option->printDetail)
        std::cout << "Total " << hier.size() << " levels.\n";
    for (auto m : hier)
    {
        auto alnStart = std::chrono::high_resolution_clock::now();
        option->calSim = (option->psgopAuto && level >= 5 && !option->psgop);
        if (option->cpuOnly || m.size() < cpuThres || util->nowProcess == 2) msaCpu(T, m, util, option, param);
        else                                                                 msaGpu(T, m, util, option, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (option->printDetail)
        {
            if (m.size() > 1)
                std::cout << "Level " << level << ", aligned " << m.size() << " pairs in " << alnTime.count() / 1000000 << " ms\n";
            else
                std::cout << "Level " << level << ", aligned " << m.size() << " pair in " << alnTime.count() / 1000000 << " ms\n";
        }
        if (option->redo)
        {
            option->psgop = true;
            return;
        }
        ++level;
    }
    // Push msa results to roots of sub-subtrees
    for (auto p : partition->partitionsRoot)
    {
        Node *current = T->allNodes[p.first];
        while (true)
        {
            int chIdx = 0;
            for (int i = 0; i < current->children.size(); ++i)
            {
                if (current->children[i]->grpID == T->allNodes[p.first]->grpID)
                {
                    chIdx = i;
                    break;
                }
            }
            if (!current->children[chIdx]->msaIdx.empty())
            {
                T->allNodes[p.first]->msaIdx = current->children[chIdx]->msaIdx;
                util->seqsLen[p.first] = util->seqsLen[current->children[chIdx]->identifier];
                break;
            }
            current = current->children[chIdx];
        }
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progessiveTime = progressiveEnd - progressiveStart;
    std::cout << "Progressive alignment (length: " << util->seqsLen[T->root->identifier] << ") in " << progessiveTime.count() / 1000000000 << " s\n";
    if (util->badSequences.empty())
        return;
    // Adding bad sequences back
    util->nowProcess = 1;
    auto badStart = std::chrono::high_resolution_clock::now();
    int badSeqBefore = 0;
    int badProfileBefore = 0;
    std::unordered_map<std::string, int> badSeqs;
    hier.clear();
    hier = std::vector<std::vector<std::pair<Node *, Node *>>>(1);
    for (auto p : partition->partitionsRoot) {
        std::vector<Node*> badNodes;
        if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end()) {
            auto badSeqName = util->badSequences[p.second.first->grpID];
            for (auto n: badSeqName) {
                badNodes.push_back(T->allNodes[n]);
                badProfileBefore += 1;
                for (auto idx: T->allNodes[n]->msaIdx) badSeqs[n] = idx;
                badSeqBefore += T->allNodes[n]->msaIdx.size();
            }   
        }
        std::sort(badNodes.begin(), badNodes.end(), comp);
        for (int i = 0; i < badNodes.size(); ++i) {
            hier[0].push_back(std::make_pair(T->allNodes[p.second.first->identifier], T->allNodes[badNodes[i]->identifier]));
        }   
    }
    util->badSequences.clear();
    std::cout << "Adding bad profiles back. Total profiles/sequences: " << badProfileBefore << " / " << badSeqBefore << '\n';
    level = 0;
    if (!hier.empty()) {
        for (auto m : hier) {
            auto alnStart = std::chrono::high_resolution_clock::now();
            if (option->cpuOnly || m.size() < cpuThres) msaCpu(T, m, util, option, param);
            else                 msaGpu(T, m, util, option, param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            if (option->printDetail) {
                if (m.size() > 1)
                    std::cout << "Level " << level << ", aligned " << m.size() << " pairs in " << alnTime.count() / 1000000 << " ms\n";
                else
                    std::cout << "Level " << level << ", aligned " << m.size() << " pair in " << alnTime.count() / 1000000 << " ms\n";
            }
            ++level;
        }
    }
    util->nowProcess = 0;
    auto badEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds badTime = badEnd - badStart;
    std::cout << "Added bad profiles/sequences in " << badTime.count() / 1000000000 << " s\n";
    return;
}

void msaGpu(Tree *tree, std::vector<std::pair<Node *, Node *>> &nodes, msa::utility *util, msa::option *option, Params &param)
{
    auto startT = std::chrono::high_resolution_clock::now();
    
    updateNode(tree, nodes, util);
    int numBlocks = _BLOCKSIZE;
    int blockSize = _THREAD_NUM;
    int gpuNum = option->gpuNum;
    // get maximum sequence/profile length
    int32_t seqLen = 0;
    for (auto n : tree->allNodes)
        seqLen = (util->seqsLen[n.first] > seqLen) ? util->seqsLen[n.first] : seqLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size() % numBlocks == 0)
        roundGPU -= 1;
    if (roundGPU < gpuNum)
        gpuNum = roundGPU;

    int32_t maxProfileLen = MAX_PROFILE_LEN;
    maxProfileLen = maxProfileLen > seqLen ? seqLen : maxProfileLen;
    int paramSize = param.matrixSize * param.matrixSize + 4;
    float *hostParam = (float *)malloc(paramSize * sizeof(float));
    for (int i = 0; i < param.matrixSize; ++i)
        for (int j = 0; j < param.matrixSize; ++j)
            hostParam[i * param.matrixSize + j] = param.scoringMatrix[i][j];
    hostParam[paramSize-4] = param.gapOpen;
    hostParam[paramSize-3] = param.gapExtend;
    hostParam[paramSize-2] = param.gapBoundary;
    hostParam[paramSize-1] = param.xdrop;
    
    int profileSize = param.matrixSize + 1;
    
    

    // allocate memory on host and device
    int8_t ***wholeAln = new int8_t **[gpuNum];
    for (int gn = 0; gn < gpuNum; ++gn)
        wholeAln[gn] = new int8_t *[numBlocks];
    for (int gn = 0; gn < gpuNum; ++gn)
        for (int nb = 0; nb < numBlocks; ++nb)
            wholeAln[gn][nb] = new int8_t[2 * seqLen];

    float **hostFreq = new float *[gpuNum];
    float **hostGapOp = new float *[gpuNum]; // gap open
    float **hostGapEx = new float *[gpuNum]; // gap extend
    int8_t **hostAln = new int8_t *[gpuNum];
    int32_t **hostLen = new int32_t *[gpuNum];
    int32_t **hostNum = new int32_t *[gpuNum];
    int32_t **hostAlnLen = new int32_t *[gpuNum];
    int32_t **hostSeqInfo = new int32_t *[gpuNum];

    float **deviceFreq = new float *[gpuNum];
    float **deviceGapOp = new float *[gpuNum];
    float **deviceGapEx = new float *[gpuNum];
    float **deviceParam = new float *[gpuNum]; // parameters
    int8_t **deviceAln = new int8_t *[gpuNum];
    int32_t **deviceLen = new int32_t *[gpuNum];
    int32_t **deviceNum = new int32_t *[gpuNum];
    int32_t **deviceAlnLen = new int32_t *[gpuNum];
    int32_t **deviceSeqInfo = new int32_t *[gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::spin_rw_mutex fallbackMutex;

    std::atomic<uint64_t> kernelTime, copyTime;
    kernelTime.store(0);
    copyTime.store(0);
    std::vector<int> fallbackPairs;

    std::vector<std::vector<int8_t>> alnBad;
    if (util->nowProcess == 1) alnBad = std::vector<std::vector<int8_t>>(nodes.size());

    if (util->nowProcess == 1) {
        float refWeight = 0.0;
        int32_t refLen = util->seqsLen[nodes[0].first->identifier];
        int32_t refNum = tree->allNodes[nodes[0].first->identifier]->msaIdx.size();
        for (auto sIdx: tree->allNodes[nodes[0].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        if (tree->allNodes[nodes[0].first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes[0].first->identifier]->msaFreq = std::vector<std::vector<float>> (refLen, std::vector<float>(profileSize,0.0));
            for (auto sIdx: tree->allNodes[nodes[0].first->identifier]->msaIdx) { 
                int storage = util->seqsStorage[sIdx];
                std::string name = util->seqsName[sIdx];
                float w = tree->allNodes[name]->weight / refWeight * refNum;
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    int letterIndex = letterIdx(option->type, toupper(util->alnStorage[storage][sIdx][s]));
                    tree->allNodes[nodes[0].first->identifier]->msaFreq[s][letterIndex] += 1.0 * w;
                }
                });
                });
            }
            for (int s = 0; s < refLen; ++s) for (int t = 0; t < profileSize; ++t) tree->allNodes[nodes[0].first->identifier]->msaFreq[s][t] /= (refNum / refWeight);
        }
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range) { 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn]   =  (float*)  malloc(profileSize * 2 * maxProfileLen * numBlocks * sizeof(float));
            hostGapOp[gn]  =  (float*)  malloc(              2 * maxProfileLen * numBlocks * sizeof(float));
            hostGapEx[gn]  =  (float*)  malloc(              2 * maxProfileLen * numBlocks * sizeof(float));
            hostAln[gn]    =  (int8_t*) malloc(              2 * maxProfileLen * numBlocks * sizeof(int8_t));
            hostLen[gn]    =  (int32_t*)malloc(              2                 * numBlocks * sizeof(int32_t));
            hostNum[gn]    =  (int32_t*)malloc(              2                 * numBlocks * sizeof(int32_t));
            hostAlnLen[gn] =  (int32_t*)malloc(                                  numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(              3                             * sizeof(int32_t));
            
            hipSetDevice(option->gpuIdx[gn]);
            hipMalloc((void**)&deviceFreq[gn], profileSize * 2 * maxProfileLen * numBlocks * sizeof(float));
            hipMalloc((void**)&deviceGapOp[gn],              2 * maxProfileLen * numBlocks * sizeof(float));
            hipMalloc((void**)&deviceGapEx[gn],              2 * maxProfileLen * numBlocks * sizeof(float));
            hipMalloc((void**)&deviceAln[gn],                2 * maxProfileLen * numBlocks * sizeof(int8_t));
            hipMalloc((void**)&deviceLen[gn],                2 *                 numBlocks * sizeof(int32_t));
            hipMalloc((void**)&deviceNum[gn],                2 *                 numBlocks * sizeof(int32_t));
            hipMalloc((void**)&deviceAlnLen[gn],                                 numBlocks * sizeof(int32_t));
            hipMalloc((void**)&deviceSeqInfo[gn],            3 *                             sizeof(int32_t));
            hipMalloc((void**)&deviceParam[gn],  paramSize     *                             sizeof(float));
            std::string error = hipGetErrorString(hipGetLastError()); if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            
            hipMemcpy(deviceParam[gn], hostParam, paramSize * sizeof(float), hipMemcpyHostToDevice);
            error = hipGetErrorString(hipGetLastError()); if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());

            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            std::vector<std::vector<int8_t>> aln;
            std::vector<std::pair<int, int>> startPos;
            std::vector<std::pair<int, int>> profileLen;
            std::vector<bool> endAln;
            float globalRefWeight = 0.0;
            if (util->nowProcess == 1) {
                int32_t refLen = util->seqsLen[nodes[0].first->identifier];
                int32_t refNum = tree->allNodes[nodes[0].first->identifier]->msaIdx.size();
                for (auto sIdx: tree->allNodes[nodes[0].first->identifier]->msaIdx)  globalRefWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
            }

            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                if (util->nowProcess == 1 && option->printDetail) std::cout << "Round " << rn+1 << '/' << roundGPU << '\n';
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                aln.clear(); startPos.clear();
                aln = std::vector<std::vector<int8_t>> (alnPairs);
                startPos = std::vector<std::pair<int, int>> (alnPairs, std::make_pair(0,0));
                profileLen = std::vector<std::pair<int, int>> (alnPairs, std::make_pair(0,0));
                int32_t profileRound = 0;
                while (true) {
                    // Initialization
                    gappyColumns.clear();
                    endAln.clear();
                    for (int i = 0; i < profileSize * 2 * maxProfileLen * numBlocks; ++i) hostFreq[gn][i] = 0;
                    for (int i = 0; i <               2 * maxProfileLen * numBlocks; ++i) hostAln[gn][i] = 0;
                    for (int i = 0; i <               2 *                 numBlocks; ++i) hostLen[gn][i] = 0;
                    for (int i = 0; i <               2 *                 numBlocks; ++i) hostNum[gn][i] = 0;
                    for (int i = 0; i <                                   numBlocks; ++i) hostAlnLen[gn][i] = 0;
                    for (int i = 0; i <               2 * maxProfileLen * numBlocks; ++i) hostGapOp[gn][i] = 0;
                    for (int i = 0; i <               2 * maxProfileLen * numBlocks; ++i) hostGapEx[gn][i] = 0;
                    for (int i = 0; i < alnPairs; ++i) {
                        std::queue<std::pair<int,int>> gappyRef, gappyQry;
                        gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                        endAln.push_back(true);
                    }
                   
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> r) {
                    for (int n = r.begin(); n < r.end(); ++n) {
                        int32_t nIdx = n + rn*numBlocks;
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        if (startPos[n].first < refLen && startPos[n].second < qryLen) {
                            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                            int32_t offsetf = profileSize * 2 * maxProfileLen * n; 
                            int32_t offsetg =               2 * maxProfileLen * n;
                            int32_t maxLenLeft = std::max(refLen-startPos[n].first, qryLen-startPos[n].second);
                            float* rawProfile = new float [profileSize * 2 * maxLenLeft];
                            for (int i = 0; i < profileSize * 2 * maxLenLeft; ++i) rawProfile[i] = 0;
                            calculateProfileFreq(rawProfile, tree, nodes[nIdx], util, option->type, maxLenLeft, profileSize, startPos[n]);
                            if (util->nowProcess == 1) {
                                for (int s = startPos[n].first; s < std::min(startPos[n].first+maxLenLeft, refLen); ++s) {
                                    for (int v = 0; v < profileSize; ++v)  {
                                        rawProfile[profileSize*(s-startPos[n].first)+v] = tree->allNodes[nodes[0].first->identifier]->msaFreq[s][v]  / globalRefWeight * refNum;
                                    }
                                }  
                            }
                            std::pair<int,int> lens = std::make_pair(refLen-startPos[n].first, qryLen-startPos[n].second);
                            removeGappyColumns(rawProfile, tree, nodes[nIdx], util, option, gappyColumns[n], maxLenLeft, maxProfileLen, lens, profileLen[n]);
                            for (int i = 0; i < profileSize * maxProfileLen; ++i)  {
                                hostFreq[gn][offsetf+i] = (i < profileSize * lens.first) ? rawProfile[i] : 0.0;
                            }
                            for (int i = 0; i < profileSize * maxProfileLen; ++i) {
                                hostFreq[gn][offsetf+profileSize*maxProfileLen+i] = (i < profileSize * lens.second) ? rawProfile[profileSize * maxLenLeft+i] : 0.0;
                            }
                            if (startPos[n].first + profileLen[n].first < refLen || startPos[n].second + profileLen[n].second < qryLen ) {
                                endAln[n] = false;
                                if (startPos[n].first == 0 && startPos[n].second == 0) hostAlnLen[gn][n] = 2;
                                else hostAlnLen[gn][n] = 3;
                            }
                            else {
                                if (startPos[n].first == 0 && startPos[n].second == 0) hostAlnLen[gn][n] = 0;
                                else hostAlnLen[gn][n] = 1;
                            }
                            std::pair<int32_t, int32_t> offset = std::make_pair(offsetf, offsetg);
                            calculatePSGOP(hostFreq[gn], hostGapOp[gn], hostGapEx[gn], tree, nodes[nIdx], util, option, maxProfileLen, offset, lens, param);
                            hostLen[gn][2*n] = lens.first; hostLen[gn][2*n+1] = lens.second;
                            hostNum[gn][2*n] = refNum;        hostNum[gn][2*n+1] = qryNum;
                            
                            delete [] rawProfile;
                        }
                        else {
                            hostAlnLen[gn][n] = -2;
                        }
                        
                    }
                    });
                    });
                    // Break condition: all pairs of alignment reach the end
                    int endNum = 0;
                    for (auto end: endAln) {
                        if (!end) break;
                        else endNum++;
                    }
                    hostSeqInfo[gn][0] = alnPairs;
                    hostSeqInfo[gn][1] = maxProfileLen;
                    hostSeqInfo[gn][2] = profileSize;

                    
                    auto copyStart = std::chrono::high_resolution_clock::now();
                    hipMemcpy(deviceFreq[gn],    hostFreq[gn],   profileSize * 2 * maxProfileLen * alnPairs * sizeof(float),   hipMemcpyHostToDevice);
                    hipMemcpy(deviceGapOp[gn],   hostGapOp[gn],                2 * maxProfileLen * alnPairs * sizeof(float),   hipMemcpyHostToDevice);
                    hipMemcpy(deviceGapEx[gn],   hostGapEx[gn],                2 * maxProfileLen * alnPairs * sizeof(float),   hipMemcpyHostToDevice);
                    hipMemcpy(deviceAln[gn],     hostAln[gn],                  2 * maxProfileLen * alnPairs * sizeof(int8_t),  hipMemcpyHostToDevice);
                    hipMemcpy(deviceLen[gn],     hostLen[gn],                  2 *                 alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
                    hipMemcpy(deviceNum[gn],     hostNum[gn],                  2 *                 alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
                    hipMemcpy(deviceAlnLen[gn],  hostAlnLen[gn],                                   alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
                    hipMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn],              3 *                            sizeof(int32_t), hipMemcpyHostToDevice);
                    auto copyEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds cTime = copyEnd - copyStart;
                    int ct = copyTime.fetch_add(cTime.count());
                    std::string berr = hipGetErrorString(hipGetLastError());
                    if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                    auto kernelStart = std::chrono::high_resolution_clock::now();
                    alignGrpToGrp_freq<<<numBlocks, blockSize>>>(
                        deviceFreq[gn],
                        deviceAln[gn], 
                        deviceLen[gn],
                        deviceNum[gn],
                        deviceAlnLen[gn],
                        deviceSeqInfo[gn], 
                        deviceGapOp[gn],
                        deviceGapEx[gn],
                        deviceParam[gn]
                    );
                    hipDeviceSynchronize();
                    auto kernelEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds kTime = kernelEnd - kernelStart;
                    int kt = kernelTime.fetch_add(kTime.count());
                    hipMemcpy(hostAln[gn],    deviceAln[gn],    2 * maxProfileLen * alnPairs * sizeof(int8_t), hipMemcpyDeviceToHost);
                    hipMemcpy(hostAlnLen[gn], deviceAlnLen[gn],                     alnPairs * sizeof(int32_t), hipMemcpyDeviceToHost);
                    std::string aerr = hipGetErrorString(hipGetLastError());
                    if (aerr != "no error") {
                        printf("ERROR: After kernel %s!\n", aerr.c_str());
                        exit(1);
                    }
                    
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                    for (int n = range.begin(); n < range.end(); ++n) {
                        int32_t nIdx = n + rn*numBlocks;
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                            
                        if ((util->nowProcess == 0) && hostAlnLen[gn][n] == -1 && (refLen > 0 && qryLen > 0)) {
                            {
                                tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                                fallbackPairs.push_back(nIdx);
                            }
                            aln[n].clear();
                            if (util->nowProcess == 1) alnBad[nIdx].clear();
                            startPos[n].first = refLen;
                            startPos[n].second = qryLen;
                        }
                        else {
                            if (util->nowProcess == 0 && (refNum == 1 || qryNum == 1) && 
                               (util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]] || 
                                util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]])) {
                                aln[n].clear();
                                {
                                    tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                                    fallbackPairs.push_back(nIdx);
                                }
                            }
                            else {
                                std::vector<int8_t> aln_reduced;
                                int alnRef = 0, alnQry = 0;
                                // if (n == 0) std::cout << "HostAlnLen: " << hostAlnLen[gn][n] << '\n';
                                if (hostAlnLen[gn][n] > 0) {
                                    for (int j = 0; j < hostAlnLen[gn][n]; ++j) aln_reduced.push_back(hostAln[gn][n*2*maxProfileLen+j]);
                                }
                                else if (refLen == 0 || qryLen == 0) {
                                    if (refLen == 0) for (int j = 0; j < qryLen; ++j) aln_reduced.push_back(1);
                                    if (qryLen == 0) for (int j = 0; j < refLen; ++j) aln_reduced.push_back(2);
                                }
                                else {
                                    std::cout << "CPU on No. " << nIdx << " ( " << nodes[nIdx].second->identifier << " )\n";
                                    int32_t offsetf = profileSize * 2 * maxProfileLen * n; 
                                    int32_t offsetg =               2 * maxProfileLen * n;
                                    int newRef = hostLen[gn][2*n], newQry = hostLen[gn][2*n+1];
                                    std::vector<std::vector<float>> freqRef (newRef, std::vector<float>(profileSize, 0.0));
                                    std::vector<std::vector<float>> freqQry (newQry, std::vector<float>(profileSize, 0.0));
                                    std::vector<std::vector<float>> gapOp (2), gapEx (2);

                                    for (int s = 0; s < newRef; s++) for (int t = 0; t < profileSize; ++t) freqRef[s][t] = hostFreq[gn][offsetf+profileSize*s+t];
                                    for (int s = 0; s < newQry; s++) for (int t = 0; t < profileSize; ++t) freqQry[s][t] = hostFreq[gn][offsetf+profileSize*(seqLen+s)+t];     
                                    for (int r = 0; r < newRef; ++r) {
                                        gapOp[0].push_back(hostGapOp[gn][offsetg+r]);
                                        gapEx[0].push_back(hostGapEx[gn][offsetg+r]);
                                    }
                                    for (int q = 0; q < newQry; ++q) {
                                        gapOp[1].push_back(hostGapOp[gn][offsetg+seqLen+q]);
                                        gapEx[1].push_back(hostGapEx[gn][offsetg+seqLen+q]);
                                    }                 
                                
                                    std::pair<float, float> num = std::make_pair(static_cast<float>(refNum), static_cast<float>(qryNum));
                                    Talco_xdrop::Params* talco_params = new Talco_xdrop::Params(hostParam, param.matrixSize);
                                    if (refLen == 0) for (int j = 0; j < qryLen; ++j) aln_reduced.push_back(1);
                                    if (qryLen == 0) for (int j = 0; j < refLen; ++j) aln_reduced.push_back(2);
                                    while (aln_reduced.empty()) {
                                        int16_t errorType = 0;
                                        aln_reduced.clear();
                                        Talco_xdrop::Align_freq (
                                            talco_params,
                                            freqRef,
                                            freqQry,
                                            gapOp,
                                            gapEx,
                                            num,
                                            aln_reduced,
                                            errorType
                                        );
                                        if (errorType == 2) {
                                            if (option->printDetail) std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                            talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->fLen * 1.2) << 1, std::min(newRef, newQry)));
                                        }
                                        else if (errorType == 3) {
                                            std::cout << "There might be some bugs in the code!\n";
                                            exit(1);
                                        }
                                        else if (errorType == 1) {
                                            if (option->printDetail) std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                            talco_params->updateXDrop(static_cast<int32_t>(talco_params->xdrop * 1.2));
                                            talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->xdrop * 4) << 1, std::min(newRef, newQry)));
                                        }
                                    }
                                    delete talco_params;
                                
                                }
                                for (auto a: aln_reduced) {
                                    // if (n == 92) std::cout << (a & 0xFFF);
                                    if (a == 0) {alnRef += 1; alnQry += 1;}
                                    if (a == 1) {alnQry += 1;}
                                    if (a == 2) {alnRef += 1;}
                                }
                                if (endAln[n]) alnRef *= -1; // Pass end condition to the below function
                                if (!option->alignGappy) alnQry *= -1;
                                std::pair<int, int> debugIdx = std::make_pair(alnRef,alnQry);
                                addGappyColumnsBack(aln_reduced, aln[n], gappyColumns[n], debugIdx);
                                
                                if ((debugIdx.first != refLen || debugIdx.second != qryLen)) {
                                    std::cout << "Name (" << nIdx << "): " << nodes[nIdx].first->identifier << '-' << nodes[nIdx].second->identifier << '\n';
                                    std::cout << "Len: " << debugIdx.first << '/' << refLen << '-' << debugIdx.second << '/' << qryLen << '\n';
                                    std::cout << "Num: " << refNum << '-' << qryNum << '\n';
                                }
                                if (util->nowProcess == 1) alnBad[nIdx] = aln[n];
                                if (endAln[n] && (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty())) std::cout << "Gappy at " << n << ':' << gappyColumns[n].first.size() << '\t' << gappyColumns[n].second.size() << '\n';
                                startPos[n].first += (profileRound == 0) ? (debugIdx.first-1) : debugIdx.first;
                                startPos[n].second += (profileRound == 0) ? (debugIdx.second-1) : debugIdx.second;
                            }
                        }
                    }
                    });
                    }); 
                    if (endNum == alnPairs) break;
                    ++profileRound;
                }

                if (util->nowProcess == 0) {                
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                    for (int n = range.begin(); n < range.end(); ++n) {
                    // for (int n = 0; n < alnPairs; ++n) {
                        int32_t nIdx = n + rn*numBlocks;
                        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        std::pair<int, int> debugIdx = std::make_pair(0,0);
                        float refWeight = 0, qryWeight = 0; // , totalWeight = 0;
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        // Update alignment & frequency
                        if (!aln[n].empty()) {
                            updateAlignment(tree, nodes[nIdx], util, aln[n]);
                            updateFrequency(tree, nodes[nIdx], util, aln[n], refWeight, qryWeight, debugIdx);
                            if (option->calSim) {
                                double colSim = calColumnSimilarity(tree, nodes[nIdx].first, util, param);
                                if (colSim < option->divTH) {
                                    option->redo = true;
                                    return;
                                }
                            }
                        }
                    }
                    });
                    }); 
                }
                
            }
        } 
    });

    if (util->nowProcess == 1) {
        std::vector<Node*> badNodes;
        for (auto n: nodes) badNodes.push_back(n.second);
        mergeInsertions (tree, nodes[0].first, badNodes, util, alnBad);
    }

    // free memory
    for (int gn = 0; gn < gpuNum; ++gn)
    {
        hipSetDevice(option->gpuIdx[gn]);
        hipFree(deviceFreq[gn]);
        hipFree(deviceAln[gn]);
        hipFree(deviceLen[gn]);
        hipFree(deviceNum[gn]);
        hipFree(deviceAlnLen[gn]);
        hipFree(deviceSeqInfo[gn]);
        hipFree(deviceParam[gn]);
        hipFree(deviceGapOp[gn]);
        hipFree(deviceGapEx[gn]);
        std::string freeErr = hipGetErrorString(hipGetLastError());
        if (freeErr != "no error")
        {
            printf("CUDA_ERROR: Free memory %s!\n", freeErr.c_str());
            exit(1);
        }
        free(hostFreq[gn]);
        free(hostAln[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAlnLen[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEx[gn]);
    }

    delete[] deviceFreq;
    delete[] deviceAlnLen;
    delete[] deviceAln;
    delete[] deviceParam;
    delete[] deviceSeqInfo;
    delete[] deviceLen;
    delete[] deviceNum;
    delete[] deviceGapOp;
    delete[] deviceGapEx;
    delete[] hostFreq;
    delete[] hostAlnLen;
    delete[] hostLen;
    delete[] hostNum;
    delete[] hostAln;
    delete[] hostSeqInfo;
    delete[] hostGapOp;
    delete[] hostGapEx;

    for (auto n : tree->allNodes)
    {
        if (n.second->is_leaf())
        {
            if (util->memLen < util->seqMemLen[util->seqsIdx[n.first]])
                util->memLen = util->seqMemLen[util->seqsIdx[n.first]];
        }
    }

    free(hostParam);
    if (fallbackPairs.empty())
        return;
    fallback2cpu(fallbackPairs, tree, nodes, util, option);
    return;
}

bool comp(Node *A, Node *B)
{
    return (A->msaIdx.size() > B->msaIdx.size());
}
