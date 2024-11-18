#ifndef MSAGPU_HPP
#include "msa-gpu.cuh"
#endif

void getGpuInfo (po::variables_map& vm, msa::option* option) {
    int maxGpuNum;
    cudaGetDeviceCount(&maxGpuNum);
    int gpuNum = (vm.count("gpu-num")) ? vm["gpu-num"].as<int>() : maxGpuNum;
    if (gpuNum <= 0) {
        std::cerr << "ERROR: Requested number of GPU <= 0.\n";
        exit(1);
    }
    if (gpuNum > maxGpuNum) {
        std::cerr << "ERROR: Requested number of GPU more than available GPUs.\n";
        exit(1);
    }
    if (option->cpuOnly) gpuNum = 0;
    else if (gpuNum > option->cpuNum) {
        if (option->cpuNum == 1) std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << option->cpuNum << " GPU.\n";
        else                     std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << option->cpuNum << " GPUs.\n";
        gpuNum = option->cpuNum;
    }
    std::vector<int> gpuIdx;
    if (vm.count("gpu-index")) {
        std::string gpuIdxString = vm["gpu-index"].as<std::string>();
        std::string id = "";
        for (int i = 0; i < gpuIdxString.size(); ++i) {
            if (gpuIdxString[i] == ',') {
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else if (i == gpuIdxString.size()-1) {
                id += gpuIdxString[i];
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else id += gpuIdxString[i];
        }
    }
    else {
        for (int i = 0; i < gpuNum; ++i) gpuIdx.push_back(i);
    }
    if (gpuIdx.size() != gpuNum) {
        std::cerr << "ERROR: the number of requested GPUs does not match the number of specified gpu indexes.\n";
        exit(1);
    }
    for (auto id: gpuIdx) {
        if (id >= maxGpuNum) {
            std::cerr << "ERROR: specified gpu index >= the number of GPUs\n";
            exit(1);
        }
    }
    printf("Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);
    option->gpuNum = gpuNum;
    option->gpuIdx = gpuIdx;
    return;
}


void msaOnSubtreeGpu (Tree* T, msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param) {
    
    auto progressiveStart = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    for (auto &p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
        int grpID = p.second.first->grpID;
        getMsaHierachy(subhier, msaStack, grpID);
        for (auto h: subhier) {
            while (hier.size() < h.second+1) {
                std::vector<std::pair<Node*, Node*>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
    }
    std::unordered_map<std::string, std::string> beforeAln;
    int level = 0;
    for (auto m: hier) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        option->calSim = (option->psgopAuto && level >= 5 && !option->psgop);
        if (option->cpuOnly || m.size() < 300 || util->nowProcess == 2) msaCpu(T, m, util, option, param);
        else if (level < 5)                                             msaGpu_s(T, m, util, option, param);
        else                                                            msaGpu(T, m, util, option, param);
        // msaGpu(T, m, util, option, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (option->printDetail) {
            if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
            else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
        }
        if (option->redo) {
            option->psgop = true;
            return;
        }
        ++level;
    }
    // Push msa results to roots of sub-subtrees
    // for (auto p: partition->partitionsRoot) {
    //     std::stack<Node*> msaStack;
    //     getPostOrderList(p.second.first, msaStack);
    //     std::vector<Node*> msaArray;
    //     while (!msaStack.empty()) {
    //         msaArray.push_back(msaStack.top());
    //         msaStack.pop();
    //     }
    //     if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
    //         if (msaArray.size() == 2) {
    //             T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
    //             util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
    //             break;
    //         }
    //         for (int m = msaArray.size()-2; m >=0; --m) {
    //             if (msaArray[m]->msaIdx.size()>0) {
    //                 T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
    //                 util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
    //                 break;
    //             }
    //         }
    //     }
    // }
    for (auto p: partition->partitionsRoot) {
        Node* current = T->allNodes[p.first];
        while (true) {
            int chIdx = 0;
            for (int i = 0; i < current->children.size(); ++i) {
                if (current->children[i]->grpID == T->allNodes[p.first]->grpID) {
                    chIdx = i;
                    break;
                }
            }
            if (!current->children[chIdx]->msaIdx.empty()) {
                T->allNodes[p.first]->msaIdx = current->children[chIdx]->msaIdx;
                util->seqsLen[p.first] = util->seqsLen[current->children[chIdx]->identifier];
                break;
            }
            current = current->children[chIdx];
        }
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progessiveTime = progressiveEnd - progressiveStart;
    std::cout << "Progressive alignment in " <<  progessiveTime.count() / 1000000000 << " s\n";
    if (util->badSequences.empty()) return;
    // Adding bad sequences back
    util->nowProcess = 1;
    auto badStart = std::chrono::high_resolution_clock::now();
    std::map<std::string, int> NodeAlnOrder;
    std::vector<std::pair<std::pair<Node*, Node*>, int>> alnOrder;
    hier.clear();
    int badSeqBefore = 0;
    int badProfileBefore = 0;
    for (auto p: partition->partitionsRoot) {
        if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end()) {
            auto badSeqName = util->badSequences[p.second.first->grpID];
            std::vector<Node*> badSeqNode;
            for (auto name: badSeqName) badSeqNode.push_back(T->allNodes[name]);
            badSeqName.clear();
            for (auto node: badSeqNode) badSeqName.push_back(node->identifier);
            badProfileBefore += badSeqNode.size();
            for (auto n: badSeqName) badSeqBefore += T->allNodes[n]->msaIdx.size();
            std::vector<std::string> nodeLeft;
            while (badSeqName.size() > 1) {
                nodeLeft.clear();
                for (int i = 0; i < badSeqName.size()-1; i+=2) {
                    int firstIdx  = (NodeAlnOrder.find(badSeqName[i]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i]]+1 : 0;
                    int secondIdx = (NodeAlnOrder.find(badSeqName[i+1]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i+1]]+1 : 0;
                    int maxIdx = max(firstIdx, secondIdx);
                    NodeAlnOrder[badSeqName[i]] = maxIdx;
                    NodeAlnOrder[badSeqName[i+1]] = maxIdx;
                    alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[badSeqName[i]], T->allNodes[badSeqName[i+1]]), maxIdx));
                    nodeLeft.push_back(badSeqName[i]);
                }
                if (badSeqName.size()%2 == 1) nodeLeft.push_back(badSeqName.back());
                badSeqName = nodeLeft;
            }
            assert(badSeqName.size() == 1);
            int idx  = (NodeAlnOrder.find(badSeqName[0]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[0]]+1 : 0;
            alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[p.second.first->identifier], T->allNodes[badSeqName[0]]), idx));
        }
    }
    std::cout << "Adding bad profiles back. Total profiles/sequences: " << badProfileBefore << " / " << badSeqBefore << '\n';
    for (auto h: alnOrder) {
        while (hier.size() < h.second+1) {
            std::vector<std::pair<Node*, Node*>> temp;
            hier.push_back(temp);
        }
        hier[h.second].push_back(h.first);
    }
    util->badSequences.clear();
    level = 0;
    if (!hier.empty()) {
        for (auto m: hier) {
            auto alnStart = std::chrono::high_resolution_clock::now();
            msaCpu(T, m, util, option, param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            if (option->printDetail) {
                if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
                else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
            }
            ++level;
        }
    }
    util->nowProcess = 0;
    auto badEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds badTime = badEnd - badStart;
    std::cout << "Added bad profiles/sequences in " <<  badTime.count() / 1000000000 << " s\n";
    return;
}

void msaGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int numBlocks = BLOCKSIZE; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum; 
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: tree->allNodes) seqLen = (util->seqsLen[n.first] > seqLen) ? util->seqsLen[n.first] : seqLen;
    seqLen *= 2;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
            
    
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
    
    // allocate memory on host and device
    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    
    float** hostGapOp = new float* [gpuNum]; // gap open
    float** hostGapEx = new float* [gpuNum]; // gap extend

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float**  deviceParam = new float* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::spin_rw_mutex fallbackMutex;
    
    std::atomic<uint64_t> kernelTime, copyTime;
    kernelTime.store(0);
    copyTime.store(0);
    std::vector<int> fallbackPairs;
    
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] =   (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] =   (int8_t*)malloc( 2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] =  (int32_t*)malloc( 2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] =  (int32_t*)malloc( 2 *          numBlocks * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(  2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(  2 * seqLen * numBlocks * sizeof(float));
            hostAlnLen[gn] = (int32_t*)malloc(            numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2                    * sizeof(int32_t));
            
            
            cudaSetDevice(option->gpuIdx[gn]);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],  2 * seqLen * numBlocks * sizeof(float));

            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                gappyColumns.clear();
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEx[gn][n] = 0;
                // Calculate Frequency
                for (int n = 0; n < alnPairs; ++n) {
                    std::queue<std::pair<int,int>> gappyRef, gappyQry;
                    gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                }
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> r) {
                for (int n = r.begin(); n < r.end(); ++n) {  
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    int32_t offsetf = 12*seqLen*n, offsetg = 2*seqLen*n;
                    calculateProfileFreq(hostFreq[gn], tree, nodes[nIdx], util, seqLen, offsetf, offsetg);
                    calculatePSGOP(hostFreq[gn], hostGapOp[gn], hostGapEx[gn], tree, nodes[nIdx], util, option, seqLen, offsetf, offsetg, param);
                    int32_t newRef = refLen, newQry = qryLen;
                    removeGappyColumns(hostFreq[gn], hostGapOp[gn], hostGapEx[gn], tree, nodes[nIdx], util, option, gappyColumns[n], newRef, newQry, seqLen, offsetf, offsetg);
                    hostLen[gn][2*n] = newRef; hostLen[gn][2*n+1] = newQry;
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                }
                });
                });
                // for (int i = 0; i < 12; ++i) {
                // std::cout << hostLen[gn][2*i] << ':' << hostLen[gn][2*i+1] << '\n';
                // }
                // for (int i = 0; i < 1000; ++i) std::cout << hostFreq[gn][12*seqLen*11+6*(i)+2] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 1000; ++i) std::cout << hostFreq[gn][12*seqLen*11+6*(seqLen+i)+2] << ',';
                // std::cout << '\n';
                
                // for (int i = 0; i < 1000; ++i) std::cout << hostGapEx[gn][2*seqLen*11+i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 1000; ++i) std::cout << hostGapOp[gn][2*seqLen*11+seqLen+i] << ',';
                // std::cout << '\n';
                
                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                auto copyStart = std::chrono::high_resolution_clock::now();
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
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
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyDeviceToHost);
                

                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                if (rn % 10 == 0 && rn > 0 && option->printDetail) std::cout << rn*numBlocks << " pairs processed.\n";
                
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    std::vector<int8_t> aln_old, aln;
                    std::pair<int, int> debugIdx;
                    float refWeight = 0, qryWeight = 0; // , totalWeight = 0;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    if (hostAlnLen[gn][n] <= 0) {
                        {
                            tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                            fallbackPairs.push_back(nIdx);
                        }
                    }
                    else {
                        for (int j = 0; j < hostAlnLen[gn][n]; ++j) aln_old.push_back(hostAln[gn][n*2*seqLen+j]);
                        addGappyColumnsBack(aln_old, aln, gappyColumns[n], debugIdx, option);
                        // if (nIdx == 11 && (refNum > 1 || qryNum > 1)) {
                        //     for (int lll = 0; lll < aln.size(); lll++) {
                        //         std::cout << (aln[lll] & 0xFFFF) << ',';
                        //     }
                        //     std::cout << '\n';
                        // }
                     
                        if (debugIdx.first != refLen || debugIdx.second != qryLen) {
                            std::cout << "Name (" << nIdx << "): " << nodes[nIdx].first->identifier << '-' << nodes[nIdx].second->identifier << '\n';
                            std::cout << "Len: " << debugIdx.first << '/' << refLen << '-' << debugIdx.second << '/' << qryLen << '\n';
                            std::cout << "Num: " << refNum << '-' << qryNum << '\n';
                        }
                        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
                    }
                    // Update alignment & frequency
                    if (!aln.empty()) {
                        updateAlignment(tree, nodes[nIdx], util, aln);
                        updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
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
    });
    // });
    
    // if (option->printDetail) std::cout << "Average kernel Time: " << kernelTime/1000000/gpuNum << "ms\n";
    // if (option->printDetail) std::cout << "Average copy Time: " << copyTime/1000000/gpuNum << "ms\n";
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
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
        if (freeErr != "no error") {
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
    
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceGapOp;
    delete [] deviceGapEx;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;

    for (auto n: tree->allNodes) {
        if (n.second->is_leaf()) {
            if (util->memLen < util->seqMemLen[util->seqsIdx[n.first]]) util->memLen = util->seqMemLen[util->seqsIdx[n.first]];
        }
    }
    
    free(hostParam);   
    if (fallbackPairs.empty()) return;
    fallback2cpu(fallbackPairs, tree, nodes, util, option);
    return;
}

void msaGpu_s(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int numBlocks = BLOCKSIZE; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
  
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) seqLen = std::max(seqLen, std::max(util->seqsLen[n.first->identifier], util->seqsLen[n.second->identifier]));
    seqLen *= 2;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;

    int32_t totalSeqNum = 0;
    int32_t* avgSeqNum = new int32_t [gpuNum];
    for (auto n: nodes) totalSeqNum += (tree->allNodes[n.first->identifier]->msaIdx.size() + tree->allNodes[n.second->identifier]->msaIdx.size());
    for (int g = 0; g < gpuNum; ++g) avgSeqNum[g] = (totalSeqNum % (nodes.size()) == 0) ? totalSeqNum / (nodes.size()) : totalSeqNum / (nodes.size()) + 1;
    
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
    
    // allocate memory on host and device
    char** hostSeqs = new char* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    
    float** hostGapOp = new float* [gpuNum]; // gap open
    float** hostGapEx = new float* [gpuNum]; // gap open
    float** hostWeight = new float* [gpuNum];

    char** deviceSeqs = new char* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceWeight = new float* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::spin_rw_mutex memMutex;
    tbb::spin_rw_mutex fallbackMutex;
    std::atomic<uint64_t> kernelTime, copyTime;
    kernelTime.store(0);
    copyTime.store(0);
    std::vector<int> fallbackPairs;
    
    
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostSeqs[gn] =  (char*)malloc(avgSeqNum[gn] * seqLen * numBlocks * sizeof(char));
            hostAln[gn] = (int8_t*)malloc(        2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(       2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(       2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(                 numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2 * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostWeight[gn] = (float*)malloc(avgSeqNum[gn]        * numBlocks * sizeof(float));
            
            cudaSetDevice(option->gpuIdx[gn]);
            // cudaSetDevice(1);
            cudaMalloc((void**)&deviceSeqs[gn],     avgSeqNum[gn] * seqLen * numBlocks * sizeof(char));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],                 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],                 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceWeight[gn], avgSeqNum[gn]              * numBlocks * sizeof(float));
            
            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                int32_t roundSeqNum = 0;
                int32_t* startNum = new int32_t [2*alnPairs];
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    startNum[2*n] = roundSeqNum;
                    startNum[2*n+1] = roundSeqNum + refNum;
                    roundSeqNum += (refNum + qryNum);
                }
                
                int32_t roundAvgNum = (roundSeqNum % alnPairs == 0) ? roundSeqNum/alnPairs : roundSeqNum/alnPairs + 1;
                if (roundAvgNum > avgSeqNum[gn]) {
                    free(hostSeqs[gn]);
                    free(hostWeight[gn]);
                    cudaFree(deviceSeqs[gn]);
                    cudaFree(deviceWeight[gn]);
                    avgSeqNum[gn] = roundAvgNum;
                    hostSeqs[gn] =    (char*)malloc(avgSeqNum[gn] * seqLen * numBlocks * sizeof(char));
                    hostWeight[gn] = (float*)malloc(avgSeqNum[gn] *          numBlocks * sizeof(float));
                    cudaMalloc((void**)&deviceSeqs[gn],   avgSeqNum[gn] * seqLen * numBlocks * sizeof(char));
                    cudaMalloc((void**)&deviceWeight[gn], avgSeqNum[gn] *          numBlocks * sizeof(float));
                }
                // Initailize 
                for (int n = 0; n < avgSeqNum[gn] * seqLen * numBlocks; ++n) hostSeqs[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEx[gn][n] = 0;
                for (int n = 0; n <  avgSeqNum[gn] * numBlocks; ++n) hostWeight[gn][n] = 0;
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> r) {
                for (int n = r.begin(); n < r.end(); ++n) {  
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    int32_t refStart = startNum[2*n];
                    int32_t qryStart = startNum[2*n+1];
                    float refWeight = 0.0, qryWeight = 0.0;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    float* gapRatioRef = new float [refLen];
                    float* gapRatioQry = new float [qryLen];
                    for (int m = 0; m < refLen; ++m) gapRatioRef[m] = 0;
                    for (int m = 0; m < qryLen; ++m) gapRatioQry[m] = 0;
                    
                    for (int idx = 0; idx < refNum; ++idx) { 
                        int sIdx = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[idx];
                        int storage = util->seqsStorage[sIdx];
                        float w = tree->allNodes[util->seqsName[sIdx]]->weight / refWeight * refNum;  
                        hostWeight[gn][refStart+idx] = w;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'A';
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'C';
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'G';
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'T';
                            // else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'N';
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                             {hostSeqs[gn][(refStart+idx)*seqLen+s] = '-'; gapRatioRef[s] += w;}
                            else                                                                                             hostSeqs[gn][(refStart+idx)*seqLen+s] = 'N';
                        }
                        });
                        });
                    }
                    for (int idx = 0; idx < qryNum; ++idx) { 
                        int sIdx = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[idx];
                        int storage = util->seqsStorage[sIdx];
                        float w = tree->allNodes[util->seqsName[sIdx]]->weight / qryWeight * qryNum;  
                        hostWeight[gn][qryStart+idx] = w;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'A';
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'C';
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'G';
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'T';
                            // else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'N';
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                             {hostSeqs[gn][(qryStart+idx)*seqLen+s] = '-'; gapRatioQry[s] += w;}
                            else                                                                                             hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'N';
                        }
                        });
                        });
                    }
                    // ClustalW's method
                    if (option->psgop) {
                        for (int s = 0; s < refLen; ++s) {
                            if (gapRatioRef[s] > 0) {
                                hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-gapRatioRef[s]) * 1.0 / refNum);
                            }
                            else {
                                int backSites = 8;
                                int distance_from_gap = 0;
                                int increPenalty = false;
                                for (int d = s-1; d > max(s-backSites, 0); --d) {
                                    ++distance_from_gap;
                                    if (gapRatioRef[d] > 0) {
                                        increPenalty = true;
                                        break;
                                    }
                                }
                                hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                            }
                            hostGapEx[gn][2*seqLen*n+s] = param.gapExtend * (1-(gapRatioRef[s]/refNum));
                        }
                        for (int s = 0; s < qryLen; ++s) {
                            if (gapRatioQry[s] > 0) {
                                hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-gapRatioQry[s]) * 1.0 / qryNum);
                            }
                            else {
                                int backSites = 8;
                                int distance_from_gap = 0;
                                int increPenalty = false;
                                for (int d = s-1; d > max(s-backSites, 0); --d) {
                                    ++distance_from_gap;
                                    if (gapRatioQry[d] > 0) {
                                        increPenalty = true;
                                        break;
                                    }
                                }
                                hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                            }
                            hostGapEx[gn][2*seqLen*n+seqLen+s] = param.gapExtend * (1-(gapRatioQry[s]/qryNum));
                        }
                    }
                    else {
                        for (int s = 0; s < refLen; ++s) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen;
                            hostGapEx[gn][2*seqLen*n+s] = param.gapExtend;
                        }
                        for (int s = 0; s < qryLen; ++s) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen;
                            hostGapEx[gn][2*seqLen*n+seqLen+s] = param.gapExtend;
                        }
                    }
                    
                    hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                    delete [] gapRatioRef;
                    delete [] gapRatioQry;
                }
                });
                });
                delete [] startNum;

                // for (int i = 0; i < 1000; ++i) std::cout << hostGapEx[gn][2*seqLen*16+i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 1000; ++i) std::cout << hostGapOp[gn][2*seqLen*16+i] << ',';
                // std::cout << '\n';

                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
                auto copyStart = std::chrono::high_resolution_clock::now();
                cudaMemcpy(deviceSeqs[gn], hostSeqs[gn], roundSeqNum * seqLen * sizeof(char), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceWeight[gn], hostWeight[gn], roundSeqNum * sizeof(float), cudaMemcpyHostToDevice);
                auto copyEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds cTime = copyEnd - copyStart;
                int ct = copyTime.fetch_add(cTime.count());

                
                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                auto kernelStart = std::chrono::high_resolution_clock::now();
                alignGrpToGrp_seq<<<numBlocks, blockSize>>>(
                    deviceSeqs[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceWeight[gn],
                    deviceParam[gn]
                );
                cudaDeviceSynchronize();
                auto kernelEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds kTime = kernelEnd - kernelStart;
                int kt = kernelTime.fetch_add(kTime.count());
                
                // cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                // cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyDeviceToHost);
                
                
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                if (rn % 10 == 0 && rn > 0 && option->printDetail) std::cout << rn*numBlocks << " pairs processed.\n";
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    std::vector<int8_t> aln_old, aln;
                    std::pair<int, int> debugIdx;
                    float refWeight = 0, qryWeight = 0; // , totalWeight = 0;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    if (hostAlnLen[gn][n] <= 0) {
                        {
                            tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                            fallbackPairs.push_back(nIdx);
                        }
                    }
                    else {
                        for (int j = 0; j < hostAlnLen[gn][n]; ++j) aln.push_back(hostAln[gn][n*2*seqLen+j]);
                    }   
                    // Update alignment & frequency
                    if (!aln.empty()) {
                        updateAlignment(tree, nodes[nIdx], util, aln);
                        if (!tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
                            updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
                        }
                    }
                }
                });
                });
               
            }  

            
        }
    });
    // });

    // if (option->printDetail) std::cout << "Average kernel Time: " << kernelTime/1000000/gpuNum << "ms\n";
    // if (option->printDetail) std::cout << "Average copy Time: " << copyTime/1000000/gpuNum << "ms\n";
    
    
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(option->gpuIdx[gn]);
        cudaFree(deviceSeqs[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEx[gn]);
        cudaFree(deviceWeight[gn]);
        // cudaDeviceSynchronize();  
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
        if (freeErr != "no error") {
            printf("CUDA_ERROR: Free memory %s!\n", freeErr.c_str());
            exit(1);
        }
        free(hostSeqs[gn]);
        free(hostAln[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAlnLen[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEx[gn]);
        free(hostWeight[gn]);
    }
    delete [] avgSeqNum;
    delete [] deviceSeqs;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceGapOp;
    delete [] deviceGapEx;
    delete [] deviceWeight;
    delete [] hostSeqs;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostWeight;

    for (auto n: tree->allNodes) {
        if (n.second->is_leaf()) {
            if (util->memLen < util->seqMemLen[util->seqsIdx[n.first]]) util->memLen = util->seqMemLen[util->seqsIdx[n.first]];
        }
    } 
    free(hostParam);
    // CPU Fallback
    if (fallbackPairs.empty()) {
        return;
    }
    fallback2cpu(fallbackPairs, tree, nodes, util, option);
    return;
}
