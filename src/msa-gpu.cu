#ifndef MSAGPU_HPP
#include "msa-gpu.cuh"
#endif

void getGpuInfo(po::variables_map &vm, msa::option *option) {
    int maxGpuNum;
    cudaGetDeviceCount(&maxGpuNum);
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
    std::cout << "Progressive alignment in " << progessiveTime.count() / 1000000000 << " s\n";
    if (util->badSequences.empty())
        return;
    // Adding bad sequences back
    util->nowProcess = 1;
    auto badStart = std::chrono::high_resolution_clock::now();
    std::map<std::string, int> NodeAlnOrder;
    std::vector<std::pair<std::pair<Node *, Node *>, int>> alnOrder;
    hier.clear();
    int badSeqBefore = 0;
    int badProfileBefore = 0;
    for (auto p : partition->partitionsRoot)
    {
        if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end())
        {
            auto badSeqName = util->badSequences[p.second.first->grpID];
            std::vector<Node *> badSeqNode;
            for (auto name : badSeqName)
                badSeqNode.push_back(T->allNodes[name]);
            std::sort(badSeqNode.begin(), badSeqNode.end(), comp);
            badSeqName.clear();
            for (auto node : badSeqNode)
                badSeqName.push_back(node->identifier);
            badProfileBefore += badSeqNode.size();
            for (auto n : badSeqName)
                badSeqBefore += T->allNodes[n]->msaIdx.size();
            std::vector<std::string> nodeLeft;
            while (badSeqName.size() > 1)
            {
                nodeLeft.clear();
                for (int i = 0; i < badSeqName.size() - 1; i += 2)
                {
                    int firstIdx = (NodeAlnOrder.find(badSeqName[i]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i]] + 1 : 0;
                    int secondIdx = (NodeAlnOrder.find(badSeqName[i + 1]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i + 1]] + 1 : 0;
                    int maxIdx = max(firstIdx, secondIdx);
                    NodeAlnOrder[badSeqName[i]] = maxIdx;
                    NodeAlnOrder[badSeqName[i + 1]] = maxIdx;
                    alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[badSeqName[i]], T->allNodes[badSeqName[i + 1]]), maxIdx));
                    nodeLeft.push_back(badSeqName[i]);
                }
                if (badSeqName.size() % 2 == 1)
                    nodeLeft.push_back(badSeqName.back());
                badSeqName = nodeLeft;
            }
            assert(badSeqName.size() == 1);
            int idx = (NodeAlnOrder.find(badSeqName[0]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[0]] + 1 : 0;
            alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[p.second.first->identifier], T->allNodes[badSeqName[0]]), idx));
        }
    }
    std::cout << "Adding bad profiles back. Total profiles/sequences: " << badProfileBefore << " / " << badSeqBefore << '\n';
    for (auto h : alnOrder)
    {
        while (hier.size() < h.second + 1)
        {
            std::vector<std::pair<Node *, Node *>> temp;
            hier.push_back(temp);
        }
        hier[h.second].push_back(h.first);
    }
    util->badSequences.clear();
    level = 0;
    if (!hier.empty())
    {
        for (auto m : hier)
        {
            auto alnStart = std::chrono::high_resolution_clock::now();
            msaCpu(T, m, util, option, param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            if (option->printDetail)
            {
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
    updateNode(tree, nodes, util);
    int numBlocks = BLOCKSIZE;
    int blockSize = THREAD_NUM;
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
    // std::cout << "MaxProfile: " << maxProfileLen << '\n';

    float *hostParam = (float *)malloc(29 * sizeof(float));
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            hostParam[i * 5 + j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;

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

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range)
                      { 
        for (int gn = range.begin(); gn < range.end(); ++gn) {

            hostFreq[gn]   =  (float*)  malloc(12 * maxProfileLen * numBlocks * sizeof(float));
            hostGapOp[gn]  =  (float*)  malloc( 2 * maxProfileLen * numBlocks * sizeof(float));
            hostGapEx[gn]  =  (float*)  malloc( 2 * maxProfileLen * numBlocks * sizeof(float));
            hostAln[gn]    =  (int8_t*) malloc( 2 * maxProfileLen * numBlocks * sizeof(int8_t));
            hostLen[gn]    =  (int32_t*)malloc( 2                 * numBlocks * sizeof(int32_t));
            hostNum[gn]    =  (int32_t*)malloc( 2                 * numBlocks * sizeof(int32_t));
            hostAlnLen[gn] =  (int32_t*)malloc(                     numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc( 2                             * sizeof(int32_t));
            
            cudaSetDevice(option->gpuIdx[gn]);
            cudaMalloc((void**)&deviceFreq[gn],   12 * maxProfileLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapOp[gn],   2 * maxProfileLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],   2 * maxProfileLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],     2 * maxProfileLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],     2 *                 numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],     2 *                 numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],                      numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 *                             sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn],  29 *                             sizeof(float));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(float), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());

            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            std::vector<std::vector<int8_t>> aln;
            std::vector<std::pair<int, int>> startPos;
            std::vector<std::pair<int, int>> profileLen;
            std::vector<bool> endAln;
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                aln.clear(); startPos.clear();
                aln = std::vector<std::vector<int8_t>> (alnPairs);
                startPos = std::vector<std::pair<int, int>> (alnPairs, std::make_pair(0,0));
                profileLen = std::vector<std::pair<int, int>> (alnPairs, std::make_pair(0,0));
                int32_t profileRound = 0;
                while (true) {
                    gappyColumns.clear();
                    endAln.clear();
                    for (int i = 0; i < 12 * maxProfileLen * numBlocks; ++i) hostFreq[gn][i] = 0;
                    for (int i = 0; i <  2 * maxProfileLen * numBlocks; ++i) hostAln[gn][i] = 0;
                    for (int i = 0; i <  2 *                 numBlocks; ++i) hostLen[gn][i] = 0;
                    for (int i = 0; i <  2 *                 numBlocks; ++i) hostNum[gn][i] = 0;
                    for (int i = 0; i <                      numBlocks; ++i) hostAlnLen[gn][i] = 0;
                    for (int i = 0; i <  2 * maxProfileLen * numBlocks; ++i) hostGapOp[gn][i] = 0;
                    for (int i = 0; i <  2 * maxProfileLen * numBlocks; ++i) hostGapEx[gn][i] = 0;
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
                            int32_t offsetf = 12*maxProfileLen*n, offsetg = 2*maxProfileLen*n;
                            int32_t maxLenLeft = std::max(refLen-startPos[n].first, qryLen-startPos[n].second);
                            float* rawProfile = new float [12 * maxLenLeft];
                            for (int i = 0; i < 12 * maxLenLeft; ++i) rawProfile[i] = 0;
                            calculateProfileFreq(rawProfile, tree, nodes[nIdx], util, maxLenLeft, startPos[n]);
                            std::pair<int,int> lens = std::make_pair(refLen-startPos[n].first, qryLen-startPos[n].second);
                            removeGappyColumns(rawProfile, tree, nodes[nIdx], util, option, gappyColumns[n], maxLenLeft, maxProfileLen, lens, profileLen[n]);
                            for (int i = 0; i < 6*maxProfileLen; ++i)  {
                                hostFreq[gn][offsetf+i] = (i < 6*lens.first) ? rawProfile[i] : 0.0;
                            }
                            for (int i = 0; i < 6*maxProfileLen; ++i) {
                                hostFreq[gn][offsetf+6*maxProfileLen+i] = (i < 6*lens.second) ? rawProfile[6*maxLenLeft+i] : 0.0;
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
                    auto copyStart = std::chrono::high_resolution_clock::now();
                    cudaMemcpy(deviceFreq[gn],    hostFreq[gn],   12 * maxProfileLen * alnPairs * sizeof(float),   cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceGapOp[gn],   hostGapOp[gn],   2 * maxProfileLen * alnPairs * sizeof(float),   cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceGapEx[gn],   hostGapEx[gn],   2 * maxProfileLen * alnPairs * sizeof(float),   cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceAln[gn],     hostAln[gn],     2 * maxProfileLen * alnPairs * sizeof(int8_t),  cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceLen[gn],     hostLen[gn],     2 *                 alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceNum[gn],     hostNum[gn],     2 *                 alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceAlnLen[gn],  hostAlnLen[gn],                      alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                    cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 *                            sizeof(int32_t), cudaMemcpyHostToDevice);
                    auto copyEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds cTime = copyEnd - copyStart;
                    int ct = copyTime.fetch_add(cTime.count());
                    std::string berr = cudaGetErrorString(cudaGetLastError());
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
                    cudaDeviceSynchronize();
                    auto kernelEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds kTime = kernelEnd - kernelStart;
                    int kt = kernelTime.fetch_add(kTime.count());
                    cudaMemcpy(hostAln[gn],    deviceAln[gn],    2 * maxProfileLen * alnPairs * sizeof(int8_t), cudaMemcpyDeviceToHost);
                    cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn],                     alnPairs * sizeof(int32_t), cudaMemcpyDeviceToHost);
                    std::string aerr = cudaGetErrorString(cudaGetLastError());
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
                        if (hostAlnLen[gn][n] == -1) {
                            {
                                tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                                fallbackPairs.push_back(nIdx);
                            }
                            aln[n].clear();
                            startPos[n].first = refLen;
                            startPos[n].second = qryLen;
                        }
                        else if (hostAlnLen[gn][n] > 0) {
                            std::vector<int8_t> aln_reduced;
                            // if (n == 0) std::cout << "HostAlnLen: " << hostAlnLen[gn][n] << '\n';
                            for (int j = 0; j < hostAlnLen[gn][n]; ++j) aln_reduced.push_back(hostAln[gn][n*2*maxProfileLen+j]);
                            int alnRef = 0, alnQry = 0;
                            for (auto a: aln_reduced) {
                                // if (n == 0) std::cout << (a & 0xFFF);
                                if (a == 0) {alnRef += 1; alnQry += 1;}
                                if (a == 1) {alnQry += 1;}
                                if (a == 2) {alnRef += 1;}
                            }
                            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                            if (endAln[n]) alnRef *= -1; // Pass end condition to the below function
                            if (!option->alignGappy) alnQry *= -1;
                            std::pair<int, int> debugIdx = std::make_pair(alnRef,alnQry);
                            addGappyColumnsBack(aln_reduced, aln[n], gappyColumns[n], debugIdx);
                            if (endAln[n] && (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty())) std::cout << "Gappy at " << n << ':' << gappyColumns[n].first.size() << '\t' << gappyColumns[n].second.size() << '\n';
                            startPos[n].first += (profileRound == 0) ? (debugIdx.first-1) : debugIdx.first;
                            startPos[n].second += (profileRound == 0) ? (debugIdx.second-1) : debugIdx.second;
                        }
                    }
                    });
                    }); 
                    if (endNum == alnPairs) break;
                    ++profileRound;
                }
                
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
        } });

    // free memory
    for (int gn = 0; gn < gpuNum; ++gn)
    {
        cudaSetDevice(option->gpuIdx[gn]);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEx[gn]);
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
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

