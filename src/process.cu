#ifndef PROCESS_HPP
#include "process.cuh"
#endif

bool cmp_msaIdx(Node* a, Node* b) {
    return a->msaIdx.size() > b->msaIdx.size();
}

void msaOnSubtree (Tree* T, msa::utility* util, msa::option* option, paritionInfo_t* partition, Params& param) {
    
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
        if (option->cpuOnly || m.size() < 1000) msaCpu(T, m, util, option, param);
        else                                    msaGpu(T, m, util, option, param);
        // if (option->cpuOnly) msaCpu(T, m, util, option, param);
        // else                 msaGpu(T, m, util, option, param);
        
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
        else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
        ++level;
    }
    // Push msa results to roots of sub-subtrees
    for (auto p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<Node*> msaArray;
        while (!msaStack.empty()) {
            msaArray.push_back(msaStack.top());
            msaStack.pop();
        }
        if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
            if (msaArray.size() == 2) {
                T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
                util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
                break;
            }
            for (int m = msaArray.size()-2; m >=0; --m) {
                if (msaArray[m]->msaIdx.size()>0) {
                    T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
                    util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
                    break;
                }
            }
        }
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progessiveTime = progressiveEnd - progressiveStart;
    std::cout << "Progressive alignment in " <<  progessiveTime.count() / 1000000000 << " s\n";
    
    if (util->badSequences.empty()) return;
    auto badStart = std::chrono::high_resolution_clock::now();
    std::cout << "Adding bad profiles back.\n";
    
    int levelThreshold = 0;
    int maxIteration = 2;
    int iteration = 0;
    bool lastIter = false;
    std::map<std::string, int> NodeAlnOrder;
    std::vector<std::pair<std::pair<Node*, Node*>, int>> alnOrder;
        
    
    while (!util->badSequences.empty() && iteration < maxIteration) {
        if (lastIter || iteration == maxIteration-1) util->nowProcess = 1;
        ++iteration;
        hier.clear();
        NodeAlnOrder.clear();
        alnOrder.clear();
        int badSeqBefore = 0, badSeqAfter = 0;
        int badProfileBefore = 0, badProfileAfter = 0;
        for (auto p: partition->partitionsRoot) {
            if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end()) {
                auto badSeqName = util->badSequences[p.second.first->grpID];
                std::vector<Node*> badSeqNode;
                for (auto name: badSeqName) badSeqNode.push_back(T->allNodes[name]);
                // Sort bad profiles based on number of sequences
                std::sort(badSeqNode.begin(), badSeqNode.end(), cmp_msaIdx);
                if (badSeqNode.size() > levelThreshold && levelThreshold != 0) {
                    int tempSize = badSeqNode.size();
                    for (int i = 0; i < tempSize-levelThreshold; ++i) {
                        badSeqNode.pop_back();
                    }
                }
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
            std::cout << "Iteraton " << iteration-1 << ". Total " << hier.size() << " levels.\n";
            for (auto m: hier) {
                auto alnStart = std::chrono::high_resolution_clock::now();
                if (!option->cpuOnly) msaGpu(T, m, util, option, param);
                else                  msaCpu(T, m, util, option, param);
                auto alnEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
                else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
                ++level;
            }
        }

        for (auto bad: util->badSequences) badProfileAfter += bad.second.size();
        for (auto bad: util->badSequences) for (auto name: bad.second) badSeqAfter += T->allNodes[name]->msaIdx.size();
        if (badSeqBefore == badSeqAfter) lastIter = true;
        std::cout << "== The number of Bad profiles/sequences ==\n";
        std::cout << "Before: " << std::setw(5) << badProfileBefore << " / " << badSeqBefore << '\n';
        std::cout << "After:  " << std::setw(5) << badProfileAfter << " / " << badSeqAfter << '\n';
        std::cout << "==========================================\n";
    }

    auto badEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds badTime = badEnd - badStart;
    std::cout << "Added bad profiles in " <<  badTime.count() / 1000000000 << " s\n";
    return;
}

void alignSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param) {
    std::vector<std::pair<Node*, Node*>> type1Aln;
    for (auto n: newT->allNodes) {
        for (auto m: n.second->children) {
            if (newT->allNodes[m->identifier]->grpID == newT->allNodes[n.second->identifier]->grpID) {
                type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
            }
        }
    }
    // for (auto n: type1Aln) std::cout << n.first->identifier << ':' << util->seqsLen[n.first->identifier] << 
    //                              ',' << n.second->identifier << ':' << util->seqsLen[n.second->identifier] <<'\n';
    std::cout << "Align sub-subtrees, total " << type1Aln.size() << " pairs.\n"; 
    if (option->cpuOnly || type1Aln.size() < 1000) createOverlapAlnCpu(T, type1Aln, util, option, param);
    else                                           createOverlapAlnGpu(T, type1Aln, util, option, param);
    return;
}

void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util) {
    for (auto n: newT->allNodes) {
        T->allNodes[n.first]->msa.clear();
        T->allNodes[n.first]->msa.push_back(n.first);
    }
    
    std::vector<std::pair<Node*, Node*>> mergePairs;
    for (auto n: newT->allNodes) {
        if (n.second->children.size() > 1) {
            mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
            for (int i = 1; i < n.second->children.size(); ++i) {
                mergePairs.push_back(std::make_pair(n.second->children[0], n.second->children[i]));
            }
        }
        else if (n.second->children.size() == 1) {
            mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
        }
    }
    
    std::vector<std::pair<Node*, Node*>> singleLevel;
    std::map<std::string, char> addedNodes;
    
    int totalEdges = 0;
    int totalLevels = 0;
    while (true) {
        auto roundStart = std::chrono::high_resolution_clock::now();
        ++totalLevels;
        addedNodes.clear();
        singleLevel.clear();
        for (auto it = mergePairs.begin(); it != mergePairs.end();) {
            Node* a = it->first;
            Node* b = it->second;
            if ((a->parent != nullptr && b->parent != nullptr) &&  
                (addedNodes.find(a->identifier) == addedNodes.end() && addedNodes.find(b->identifier) == addedNodes.end())) {
                singleLevel.push_back(std::make_pair(a, b));
                for (auto id: T->allNodes[a->identifier]->msa) addedNodes[id] = 0;
                for (auto id: T->allNodes[b->identifier]->msa) addedNodes[id] = 0;
                mergePairs.erase(it);
            }
            else {
                ++it;
            }
        }
        bool breakLoop = false;
        if (singleLevel.empty()) {
            for (auto mp: mergePairs) {
                singleLevel.push_back(mp);
            }
            breakLoop = true;
        }
        transitivityMerge(T, newT, singleLevel, util);
        for (auto n: singleLevel) {
            // std::cout << T->allNodes[n.first->identifier]->msa.size() << '/' << T->allNodes[n.second->identifier]->msa.size() << '\n';
            // std::cout << n.first->identifier << '\n';
            // for (auto s: T->allNodes[n.first->identifier]->msa) std::cout << s << ',';
            // std::cout << '\n';
            // std::cout << n.second->identifier << '\n';
            // for (auto s: T->allNodes[n.second->identifier]->msa) std::cout << s << ',';
            // std::cout << '\n';
            if (util->nowProcess == 2 && n.first->parent == nullptr) T->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
            std::vector<std::string> temp = T->allNodes[n.first->identifier]->msa;
            
            for (int r = 1; r < T->allNodes[n.first->identifier]->msa.size(); ++r) {
                std::string grpNode = T->allNodes[n.first->identifier]->msa[r];
                for (auto id: T->allNodes[n.second->identifier]->msa) T->allNodes[grpNode]->msa.push_back(id);
            }
            for (auto id: T->allNodes[n.second->identifier]->msa) T->allNodes[n.first->identifier]->msa.push_back(id);
            for (int r = 1; r < T->allNodes[n.second->identifier]->msa.size(); ++r) {
                std::string grpNode = T->allNodes[n.second->identifier]->msa[r];
                for (auto id: temp) T->allNodes[grpNode]->msa.push_back(id);
            }
            for (auto id: temp) T->allNodes[n.second->identifier]->msa.push_back(id);
            // std::cout << T->allNodes[n.first->identifier]->msa.size() << '/' << T->allNodes[n.second->identifier]->msa.size() << '\n';
            // std::cout << n.first->identifier << '\n';
            // for (auto s: T->allNodes[n.first->identifier]->msa) std::cout << s << ',';
            // std::cout << '\n';
            // std::cout << n.second->identifier << '\n';
            // for (auto s: T->allNodes[n.second->identifier]->msa) std::cout << s << ',';
            // std::cout << '\n';
        }
        auto roundEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds roundTime = roundEnd - roundStart;
        if (singleLevel.size() > 1) {
            std::cout << "Merged "<< singleLevel.size() << " edges in " << roundTime.count() / 1000000 << " ms\n";
        }
        else {
            std::cout << "Merged "<< singleLevel.size() << " edge in " << roundTime.count() / 1000000 << " ms\n";
        }
        totalEdges += singleLevel.size();
        if (breakLoop) break;
        if (totalLevels % 100 == 0) std::cout << "=============" << totalLevels << "=================\n";
    }
    std::cout << "Total Edges/Levels: " << totalEdges << '/' << totalLevels << '\n';
    return;
}

void createOverlapAlnGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    int32_t seqLen = 0;
    if (util->nowProcess < 2) seqLen = util->memLen;
    else {
        for (auto pf: util->profileFreq) if (pf.second[0].size() > seqLen) seqLen = pf.second[0].size();
    }
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;

    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 

    
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;

    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    float** hostGapOp = new float* [gpuNum];
    float** hostGapEx = new float* [gpuNum];
    float** hostGapCl = new float* [gpuNum];

    

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceGapCl = new float* [gpuNum];
  
    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] = (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2 * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapCl[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            
            
            cudaSetDevice(gn);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            cudaMalloc((void**)&deviceGapOp[gn], 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn], 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapCl[gn], 2 * seqLen * numBlocks * sizeof(float));


            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            // error = cudaGetLastError(); printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            std::vector<std::pair<int, int>> seqIdx;
            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            
            
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEx[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCl[gn][n] = 0;
                gappyColumns.clear();
                for (int n = 0; n < alnPairs; ++n) {
                    std::queue<std::pair<int,int>> gappyRef, gappyQry;
                    gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                }
                // Calculate Frequency
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    
                    if (util->nowProcess < 2) {
                        float refWeight = 0.0; float qryWeight = 0.0;
                        // get sum of weights
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        if (tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
                            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                                int storage = util->seqsStorage[sIdx];
                                std::string name = util->seqsName[sIdx];
                                float w = tree->allNodes[name]->weight / refWeight * refNum;
                                tbb::this_task_arena::isolate( [&]{
                                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                                for (int s = r.begin(); s < r.end(); ++s) {
                                    if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                             util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                                }
                                });
                                });
                            }
                        }
                        else {
                            for (int s = 0; s < tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.size(); ++s) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*s+t] = tree->allNodes[nodes[nIdx].first->identifier]->msaFreq[s][t] / refWeight * refNum;
                            }
                        }
                        if (tree->allNodes[nodes[nIdx].second->identifier]->msaFreq.empty()) { 
                            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                                int storage = util->seqsStorage[sIdx];
                                std::string name = util->seqsName[sIdx];
                                float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                                tbb::this_task_arena::isolate( [&]{
                                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                                for (int s = r.begin(); s < r.end(); ++s) {
                                    if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                             util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                                    else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                                }
                                });
                                });
                            }
                        }
                        else {
                            for (int s = 0; s < tree->allNodes[nodes[nIdx].second->identifier]->msaFreq.size(); ++s) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+s)+t] = tree->allNodes[nodes[nIdx].second->identifier]->msaFreq[s][t] / qryWeight * qryNum;
                            }
                        }
                        for (int s = 0; s < refLen; ++s) { hostGapOp[gn][2*seqLen*n+s] /= refNum; hostGapEx[gn][2*seqLen*n+s] /= refNum; hostGapCl[gn][2*seqLen*n+s] /= refNum;} 
                        for (int s = 0; s < qryLen; ++s) { hostGapOp[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapEx[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapCl[gn][2*seqLen*n+seqLen+s] /= qryNum;}
                    }
                    else {
                        int subtreeRef = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
                        for (int i = 0; i < 6; ++i) {
                            for (int s = 0; s < refLen; ++s) {
                                hostFreq[gn][12*seqLen*n+6*s+i] = util->profileFreq[subtreeRef][i][s]; 
                            }
                        }
                        for (int s = 0; s < refLen; ++s) hostGapEx[gn][2*seqLen*n+s] = util->profileFreq[subtreeRef][5][s]; 
                        int subtreeQry = tree->allNodes[nodes[nIdx].second->identifier]->grpID;
                        for (int i = 0; i < 6; ++i) {
                            for (int s = 0; s < qryLen; ++s) {
                                hostFreq[gn][12*seqLen*n+6*(seqLen+s)+i] = util->profileFreq[subtreeQry][i][s]; 
                            }
                        }
                        for (int s = 0; s < refLen; ++s) hostGapEx[gn][2*seqLen*n+seqLen+s] = util->profileFreq[subtreeQry][5][s]; 
                    }
                    // Clustalw's method
                    for (int s = 0; s < refLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5])*1.0 / refNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    for (int s = 0; s < qryLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    
                    if (option->gappyHorizon > 0) {
                        float gappyVertical = option->gappyVertical;
                        int gappyHorizon = option->gappyHorizon, gappyLength;
                        int rawIdx = 0, newIdx = 0;
                        int gapRef = 0, gapQry = 0;
                        while (true) {
                            if (rawIdx >= refLen) {
                                for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*i+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*gappyLength+5]/refNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].first.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHR: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapRef += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, refLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*newIdx+t] = hostFreq[gn][12*seqLen*n+6*rawIdx+t];
                                hostGapOp[gn][2*seqLen*n+newIdx] = hostGapOp[gn][2*seqLen*n+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+newIdx] = hostGapEx[gn][2*seqLen*n+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+newIdx] = hostGapCl[gn][2*seqLen*n+rawIdx];
                                ++newIdx;
                            } 
                        }
                        hostLen[gn][2*n] = newIdx; 
                        // std::cout << n << '#' << rawIdx << '/' << newIdx << '\t';
                        rawIdx = 0; newIdx = 0;
                        while (true) {
                            if (rawIdx >= qryLen) {
                                for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*(seqLen+i)+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].second.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHQ: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapQry += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, qryLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+newIdx)+t] = hostFreq[gn][12*seqLen*n+6*(seqLen+rawIdx)+t];
                                hostGapOp[gn][2*seqLen*n+seqLen+newIdx] = hostGapOp[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+seqLen+newIdx] = hostGapEx[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+seqLen+newIdx] = hostGapCl[gn][2*seqLen*n+seqLen+rawIdx];
                                ++newIdx;
                            } 
                        }
                        // if (!gappyColumns[n].first.empty()) std::cout << "No. " << n << ": Removing " << gapRef << " columns from refernce.\n";
                        // if (!gappyColumns[n].second.empty()) std::cout << "No. " << n << ": Removing " << gapQry << " columns from query.\n";
                        hostLen[gn][2*n+1] = newIdx;
                        if (hostLen[gn][2*n] + gapRef != refLen) std::cout << "REF:" << hostLen[gn][2*n] << '+' << gapRef << " != " << refLen << '\n';
                        if (hostLen[gn][2*n+1] + gapQry != qryLen) std::cout << "QRY:" << hostLen[gn][2*n+1] << '+' << gapQry << " != " << qryLen << '\n';
                        
                        assert(hostLen[gn][2*n] + gapRef == refLen);
                        assert(hostLen[gn][2*n+1] + gapQry == qryLen);
                    }
                    else {
                        hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    }
                    
                    // hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;    
                }

                // for (int i = 0; i < 50; ++i) std::cout << hostFreq[gn][12*seqLen*0+6*i+0] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostFreq[gn][12*seqLen*0+6*i+5] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostFreq[gn][12*seqLen*0+6*i+3] << ',';
                // std::cout << '\n';


                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);


                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceGapCl[gn],
                    deviceParam[gn]
                );
                // cudaDeviceSynchronize();
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
                
                tbb::this_task_arena::isolate([&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                    // if (nIdx % 100 == 0) {
                    if (hostAlnLen[gn][n] <= 0) {
                        std::vector<int8_t> aln, aln_old;
                        std::vector<std::vector<float>> freqRef;
                        std::vector<std::vector<float>> freqQry;
                        std::vector<std::vector<float>> gapOp, gapEx, gapCl;
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    
            
                        for (int r = 0; r < refLen; r++) {
                            std::vector<float> temp;
                            for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*r+f]);
                            freqRef.push_back(temp);
                        }
                        for (int q = 0; q < qryLen; q++) {
                            std::vector<float> temp;
                            for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*(seqLen+q)+f]);
                            freqQry.push_back(temp);
                        }
                        for (int i = 0; i < 2; ++i) {
                            std::vector<float> temp;
                            gapOp.push_back(temp); gapEx.push_back(temp); gapCl.push_back(temp);
                        }
                        for (int r = 0; r < hostLen[gn][2*n]; ++r) {
                            gapOp[0].push_back(hostGapOp[gn][2*n*seqLen+r]);
                            gapEx[0].push_back(hostGapEx[gn][2*n*seqLen+r]);
                            gapCl[0].push_back(hostGapCl[gn][2*n*seqLen+r]);
                        }
                        for (int q = 0; q < hostLen[gn][2*n+1]; ++q) {
                            gapOp[1].push_back(hostGapOp[gn][2*n*seqLen+seqLen+q]);
                            gapEx[1].push_back(hostGapEx[gn][2*n*seqLen+seqLen+q]);
                            gapCl[1].push_back(hostGapCl[gn][2*n*seqLen+seqLen+q]);
                        }
                        // for (int s = 0; s < refLen; ++s) for (int j = 0; j < 6; ++j) freqRef[s][j] = hostFreq[gn][12*seqLen*n+6*s+j];
                        // for (int s = 0; s < qryLen; ++s) for (int j = 0; j < 6; ++j) freqQry[s][j] = hostFreq[gn][12*seqLen*n+6*(seqLen+s)+j];
                        std::pair<int32_t, int32_t> num = std::make_pair(refNum, qryNum);
                        Talco_xdrop::Params talco_params(hostParam);
                        while (aln.empty()) {
                            int16_t errorType = 0;
                            Talco_xdrop::Align_freq (
                                talco_params,
                                freqRef,
                                freqQry,
                                gapOp,
                                gapEx,
                                gapCl,
                                num,
                                aln_old,
                                errorType
                            );
                            // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                            if (errorType == 1) {
                                std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                talco_params.updateXDrop(2*talco_params.xdrop);
                            }
                            if (errorType == 2) {
                                std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                talco_params.updateFLen(talco_params.fLen << 1);
                            }
                        }
                        int j = 0, rIdx = 0, qIdx = 0;
                        while (j < aln_old.size() || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty())) {
                        // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                            bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                            int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                            int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                            if (gapR || gapQ) {
                                // std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                if (gapRLen >= gapQLen) {
                                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                }
                                else {
                                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                }
                                if (gapR) gappyColumns[n].first.pop();
                                if (gapQ) gappyColumns[n].second.pop();
                            }
                            else {
                                switch ((aln_old[j] & 0xFFFF)) {
                                    case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                    case 1: ++qIdx;         aln.push_back(1); break;
                                    case 2: ++rIdx;         aln.push_back(2); break;
                                }
                                ++j;
                            }
                            if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                for (j = j; j < aln_old.size(); ++j) {
                                    switch ((aln_old[j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                }
                            }
                        }
                        assert(rIdx == refLen); assert(qIdx == qryLen);
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                        // std::cout << "CPU fallback on No. " << n << " (" << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ")\n";
                        printf("CPU fallback on No. %d (%s), Alignment Length: %lu\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
                    }
                    else {
                        std::vector<int8_t> aln;
                        int j = 0, rIdx = 0, qIdx = 0; 
                        int oneInAln = 0;
                        while (j < hostAlnLen[gn][n] || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty()) ) {
                        // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                            bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                            int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                            int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                            if (gapR || gapQ) {
                                // if (n == 0 ) std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                if (gapRLen >= gapQLen) {
                                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                }
                                else {
                                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                }
                                if (gapR) gappyColumns[n].first.pop();
                                if (gapQ) gappyColumns[n].second.pop();
                            }
                            else {
                                switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                    case 0: ++rIdx; ++qIdx; aln.push_back(0); ++oneInAln; break;
                                    case 1: ++qIdx;         aln.push_back(1); ++oneInAln; break;
                                    case 2: ++rIdx;         aln.push_back(2); break;
                                }
                                ++j;
                            }
                            if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                for (j = j; j < hostAlnLen[gn][n]; ++j) {
                                    switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                }
                            }
                        }
                        assert(rIdx == refLen); assert(qIdx == qryLen);
                        // std::cout << "ALN:" << aln.size() << '\n';
                        // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                        //     aln.push_back(hostAln[gn][n*2*seqLen+j]);
                        // }
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                    }
                } 
                });   
                });
            }  
        }
    });
    
    // for (auto n: nodes) {
    //     tree->allNodes[n.first->identifier]->msa.clear();
    //     tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
    //     tree->allNodes[n.second->identifier]->msa.clear();
    //     tree->allNodes[n.second->identifier]->msa.push_back(n.second->identifier);    
    // }
    
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        // cudaSetDevice(1);
        cudaSetDevice(gn);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEx[gn]);
        cudaFree(deviceGapCl[gn]);
        cudaDeviceSynchronize();  
        free(hostFreq[gn]);
        free(hostAlnLen[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAln[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEx[gn]);
        free(hostGapCl[gn]);
    }
    
    free(hostParam);

    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceGapOp;
    delete [] deviceGapEx;
    delete [] deviceGapCl;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostGapCl;
    return;
}

void createOverlapAlnCpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    int32_t seqLen = 0;
    if (util->nowProcess < 2) seqLen = util->memLen;
    else {
        for (auto pf: util->profileFreq) if (pf.second.size() > seqLen) seqLen = pf.second.size();
    }
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
  
    // tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    // for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
    for (int nIdx = 0; nIdx < nodes.size(); ++nIdx) {
        float* hostFreq  = (float*)malloc(12 * seqLen * sizeof(float));
        float* hostGapOp = (float*)malloc(2 * seqLen * sizeof(float));
        float* hostGapEx = (float*)malloc(2 * seqLen * sizeof(float));
        float* hostGapCl = (float*)malloc(2 * seqLen * sizeof(float));
        for (int n = 0; n < 12*seqLen; ++n) hostFreq[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapOp[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapEx[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapCl[n] = 0;
        std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>> gappyColumns;
        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
        int32_t newRef = refLen, newQry = qryLen;
        float refNum = 0, qryNum = 0;
        if (util->nowProcess < 2) {
            calculateProfileFreq(hostFreq, hostGapOp, tree, nodes[nIdx], util, seqLen, param);
            refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        
        }
        else {
            int subtreeRef = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            for (int t = 0; t < 6; ++t) refNum += util->profileFreq[subtreeRef][0][t]; 
            for (int s = 0; s < refLen; ++s) for (int t = 0; t < 6; ++t) hostFreq[6*s+t] = util->profileFreq[subtreeRef][s][t]; 
            for (int s = 0; s < refLen; ++s) hostGapEx[s] = util->profileFreq[subtreeRef][s][5] / refNum; 
            int subtreeQry = tree->allNodes[nodes[nIdx].second->identifier]->grpID;
            for (int t = 0; t < 6; ++t) qryNum += util->profileFreq[subtreeQry][0][t]; 
            for (int s = 0; s < qryLen; ++s) for (int t = 0; t < 6; ++t) hostFreq[6*(seqLen+s)+t] = util->profileFreq[subtreeQry][s][t]; 
            for (int s = 0; s < qryLen; ++s) hostGapEx[seqLen+s] = util->profileFreq[subtreeQry][s][5] / qryNum; 
            // Clustalw's method
            for (int s = 0; s < refLen; ++s) {
                if (hostFreq[6*s+5] > 0) {
                    hostGapOp[s] = param.gapOpen * 0.3 * ((refNum-hostFreq[6*s+5])*1.0 / refNum);
                }
                else {
                    int backSites = 8;
                    int distance_from_gap = 0;
                    int increPenalty = false;
                    for (int d = s-1; d > max(s-backSites, 0); --d) {
                        ++distance_from_gap;
                        if (hostFreq[6*d+5] > 0) {
                            increPenalty = true;
                            break;
                        }
                    }
                    hostGapOp[s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                }
            }
            for (int s = 0; s < qryLen; ++s) {
                if (hostFreq[6*(seqLen+s)+5] > 0) {
                    hostGapOp[seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[6*(seqLen+s)+5]) * 1.0 / qryNum);
                }
                else {
                    int backSites = 8;
                    int distance_from_gap = 0;
                    int increPenalty = false;
                    for (int d = s-1; d > max(s-backSites, 0); --d) {
                        ++distance_from_gap;
                        if (hostFreq[6*(seqLen+d)+5] > 0) {
                            increPenalty = true;
                            break;
                        }
                    }
                    hostGapOp[seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                }
            }
        }
        if (option->gappyHorizon > 0) removeGappyColumns(hostFreq, hostGapOp, tree, nodes[nIdx], util, option, gappyColumns, newRef, newQry, seqLen);
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+0] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+1] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+2] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+3] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+5] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostGapEx[seqLen+i] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[seqLen+i] << ',';
        // std::cout << '\n';
        // exit(1);
        std::vector<int8_t> aln_old, aln;
        std::vector<std::vector<float>> freqRef (newRef, std::vector<float>(6, 0.0));
        std::vector<std::vector<float>> freqQry (newQry, std::vector<float>(6, 0.0));
        std::vector<std::vector<float>> gapOp (2), gapEx (2), gapCl (2);    
        for (int s = 0; s < newRef; s++) for (int t = 0; t < 6; ++t) freqRef[s][t] = hostFreq[6*s+t];
        for (int s = 0; s < newQry; s++) for (int t = 0; t < 6; ++t) freqQry[s][t] = hostFreq[6*(seqLen+s)+t];     
        for (int r = 0; r < newRef; ++r) {
            gapOp[0].push_back(hostGapOp[r]);
            gapEx[0].push_back(freqRef[r][5] / refNum);
            gapCl[0].push_back(hostGapCl[r]);
        }
        for (int q = 0; q < newQry; ++q) {
            gapOp[1].push_back(hostGapOp[seqLen+q]);
            gapEx[1].push_back(freqQry[q][5] / qryNum);
            gapCl[1].push_back(hostGapCl[seqLen+q]);
        }
        Talco_xdrop::Params talco_params(hostParam);
        std::pair<int32_t, int32_t> num = std::make_pair(static_cast<int32_t>(refNum), static_cast<int32_t>(qryNum));
        while (aln_old.empty()) {
            int16_t errorType = 0;
            Talco_xdrop::Align_freq (
                talco_params,
                freqRef,
                freqQry,
                gapOp,
                gapEx,
                gapCl,
                num,
                aln_old,
                errorType
            );
            // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
            if (errorType == 1) {
                std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                talco_params.updateXDrop(2*talco_params.xdrop);
            }
            if (errorType == 2) {
                std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                talco_params.updateFLen(talco_params.fLen << 1);
            }
        }
        std::pair<int, int> debugIdx;
        addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx);
        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
        free(hostFreq);
        free(hostGapCl);
        free(hostGapEx);
        free(hostGapOp);        
    }
    // });
    
    free(hostParam);
    return;
}

void transitivityMerge(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util) {
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range) {
    for (int s = range.begin(); s < range.end(); ++s) {
        auto n = nodes[s];
    // for (auto n: nodes) {
        if (newtree->allNodes[n.first->identifier]->parent == nullptr) {
            for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
            std::vector<int8_t> rootAln;
            for (int i = 0; i < tree->allNodes[n.second->identifier]->msaAln.size(); ++i) {
                if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 2) rootAln.push_back(1);
                else if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 1) rootAln.push_back(2);
                else rootAln.push_back(tree->allNodes[n.second->identifier]->msaAln[i]);
            }
            tree->allNodes[n.first->identifier]->msaAln = rootAln;

            int seqLen = tree->allNodes[n.first->identifier]->msaAln.size();
            util->seqsLen[n.first->identifier] = seqLen;
            if (util->nowProcess < 2) {
                util->memCheck(seqLen);
                for (auto id: tree->allNodes[n.first->identifier]->msa) {
                    // auto id = tree->root->msa[k];
                    std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
                    for (auto sIdx: tree->allNodes[id]->msaIdx) {
                        int orgIdx = 0;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        for (int j = 0; j < aln.size(); ++j) {
                            if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                                util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->alnStorage[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[id] = aln.size();
                        util->changeStorage(sIdx);
                    }
                }
                // });
            }
            for (auto id: tree->allNodes[n.first->identifier]->msa) {
                if (id != tree->allNodes[n.first->identifier]->identifier) {
                    for (auto Idx: tree->allNodes[id]->msaIdx) {
                        tree->allNodes[n.first->identifier]->msaIdx.push_back(Idx);
                    }
                }
            }
            tree->allNodes[n.first->identifier]->msa.clear();
            return;
            // else {
            //     util->seqLen = seqLen;
            //     if (util->seqNum != 0) {
            //         for (auto id: tree->allNodes[n.first->identifier]->msa) {
            //             std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
            //             for (auto sIdx: tree->allNodes[id]->msaIdx) {
            //                 char* seq = new char[seqLen+1];
            //                 int storage = util->seqsStorage[sIdx];
            //                 int orgIdx = 0;
            //                 for (int j = 0; j < seqLen+1; ++j) {
            //                     if (j < aln.size()) {
            //                         if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
            //                             seq[j] = util->seqs[sIdx][orgIdx];
            //                             orgIdx++;
            //                         }
            //                         else {
            //                             seq[j] = '-';
            //                         }
            //                     }
            //                     else {
            //                         seq[j] = 0;
            //                     }

            //                 }
            //                 util->seqsLen[id] = aln.size();
            //                 // util->changeStorage(sIdx);
            //                 delete [] util->seqs[sIdx];
            //                 util->seqs[sIdx] = seq;
            //             }
            //         }
            //         // });
            //         for (auto id: tree->allNodes[n.first->identifier]->msa) {
            //             if (id != tree->allNodes[n.first->identifier]->identifier) {
            //                 for (auto Idx: tree->allNodes[id]->msaIdx) {
            //                     tree->allNodes[n.first->identifier]->msaIdx.push_back(Idx);
            //                 }
            //             }
            //         }
            //         tree->allNodes[n.first->identifier]->msa.clear();
                
            //     }
            //     continue;
            // }
            
        }
        int8_t refGap = (newtree->allNodes[n.first->identifier]->level == newtree->allNodes[n.second->identifier]->level) ? 2 : 1; 
        std::vector<std::vector<int8_t>> refAln, qryAln;
        std::vector<std::vector<int8_t>> refNewAln, qryNewAln;
        // refAln.push_back(tree->allNodes[n.first->identifier]->msaAln);
        // qryAln.push_back(tree->allNodes[n.second->identifier]->msaAln);
        for (auto id: tree->allNodes[n.first->identifier]->msa)  refAln.push_back(tree->allNodes[id]->msaAln);
        for (auto id: tree->allNodes[n.second->identifier]->msa) qryAln.push_back(tree->allNodes[id]->msaAln);
        
        std::vector<int8_t> refOpAln = tree->allNodes[n.first->identifier]->msaAln;
        std::vector<int8_t> qryOpAln = tree->allNodes[n.second->identifier]->msaAln;
        int32_t refLen = refOpAln.size();
        int32_t qryLen = qryOpAln.size();
        int32_t refNum = refAln.size();
        int32_t qryNum = qryAln.size();
        int32_t seqLen = max(refLen, qryLen);
        for (int i = 0; i < refNum; ++i) {
            std::vector<int8_t> temp;
            refNewAln.push_back(temp);
        }
        for (int i = 0; i < qryNum; ++i) {
            std::vector<int8_t> temp;
            qryNewAln.push_back(temp);
        }
        // std::cout << refLen << ':';
        // for (int i = 0; i < refLen; ++i) {
        //     if ((refOpAln[i] & 0XFFFF) == refGap || (refOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // std::cout << qryLen << ':';
        // for (int i = 0; i < qryLen; ++i) {
        //     if ((qryOpAln[i] & 0XFFFF) == 2 || (qryOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // aln = 0: match
        // aln = 1: ref gap
        // aln = 2: qry gap
        // aln = 3: permenant gap
        int32_t rIdx = 0, qIdx = 0;
        while (rIdx < refLen && qIdx < qryLen) {
            if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                qIdx++;rIdx++;
            }
            else if (((refOpAln[rIdx] & 0xFFFF) == refGap || (refOpAln[rIdx] & 0xFFFF) == 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                int consecGap = 0;
                int k = rIdx;
                while (((refOpAln[k] & 0xFFFF) == refGap || (refOpAln[k] & 0xFFFF) == 3) && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
            }
            else if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) == 2 || (qryOpAln[qIdx] & 0xFFFF) == 3)) {
                int consecGap = 0;
                int k = qIdx;
                while (((qryOpAln[k] & 0xFFFF) == 2 || (qryOpAln[k] & 0xFFFF) == 3) && k < qryLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (((refOpAln[kr] & 0xFFFF) == refGap || (refOpAln[kr] & 0xFFFF) == 3) && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
                consecGap = 0;
                while (((qryOpAln[kq] & 0xFFFF) == 2 || (qryOpAln[kq] & 0xFFFF) == 3) && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
        }
        // printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][g]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);    
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][g]);
            }
        }
        assert (refNewAln[0].size() == qryNewAln[0].size());
        for (int i = 0; i < tree->allNodes[n.first->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.first->identifier]->msa[i];
            tree->allNodes[id]->msaAln = refNewAln[i];
        } 
        for (int i = 0; i < tree->allNodes[n.second->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.second->identifier]->msa[i];
            tree->allNodes[id]->msaAln = qryNewAln[i];
        }
    }
    });
    
    return;
}

/*
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{

    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        // is leaf
        if (n.first->is_leaf() && n.first->msaIdx.empty()) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.first->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->msaIdx.empty()) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.second->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        // if (nodes.size() < 5) {
        //     std::cout << n.first->identifier << '(' << util->seqsLen[n.first->identifier] << ')' <<  tree->allNodes[n.first->identifier]->msaIdx.size() << '\t';
        //     std::cout << n.second->identifier << '(' << util->seqsLen[n.second->identifier] << ')'<< tree->allNodes[n.second->identifier]->msaIdx.size() << '\n';
        // }
    }

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
            
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
            
    
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
    
    

    // std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    // allocate memory on host and device
    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    
    float** hostGapOp = new float* [gpuNum]; // gap open
    float** hostGapEx = new float* [gpuNum]; // gap extend
    float** hostGapCl = new float* [gpuNum]; // gap close

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceGapCl = new float* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];

    // float** rawFreq =  new float* [gpuNum];
    // float** rawGapOp = new float* [gpuNum];
    // float** rawGapCt = new float* [gpuNum];
    // float** rawGapEn = new float* [gpuNum];
    // for (int gn = 0; gn < gpuNum; ++gn) {
    //     rawFreq[gn] = (float*)malloc(12 * seqLen * sizeof(float));    
    //     rawGapOp[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    //     rawGapCt[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    //     rawGapEn[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    // }

    std::atomic<int> nowRound;
    tbb::mutex memMutex;
    nowRound.store(0);
    // nowRound.store(roundGPU-1);

    // int ThreadsPerGPU = maxThreads / gpuNum;
    bool* cpuFallback = new bool[nodes.size()];
    std::vector<std::vector<int8_t>> alnCpu;
    for (int i = 0; i < nodes.size(); ++i) {
        cpuFallback[i] = false;
        if (util->nowProcess == 1) {
            std::vector<int8_t> aln;
            alnCpu.push_back(aln);
        }
        
    }
    
            

    // tbb::task_arena limited_arena(2);
    // limited_arena.execute([&]{

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] =  (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2                     * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapCl[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            
            cudaSetDevice(option->gpuIdx[gn]);
            // cudaSetDevice(1);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapCl[gn],  2 * seqLen * numBlocks * sizeof(float));

            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            
            while (nowRound < roundGPU) {
            // while (nowRound >= 0) {
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
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCl[gn][n] = 0;

                // Calculate Frequency
                // tbb::task_arena ta {maxThreads - gpuNum};
                // tbb::task_group tg;
                // tg.run ([&] () {
                // tbb::this_task_arena::isolate ([&] () {
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> aln_range) { 
                // for (int n = aln_range.begin(); n < aln_range.end(); ++n) {
                
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
                    float refWeight = 0.0; float qryWeight = 0.0;

                    // get sum of weights
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        float w = tree->allNodes[name]->weight / refWeight * refNum;
                        // float w = 1.0;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < refLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                            if (s > 0) {
                                // if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+s] += 1.0*w;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+s] += 1.0*w;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+s] += 1.0*w;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+s] += 1.0*w;
                            }
                            
                        }
                        });
                        });
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        // float w = 1.0;
                        float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < qryLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                            if (s > 0) {
                                // if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+seqLen+s] += 1.0*w;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+seqLen+s] += 1.0*w;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+seqLen+s] += 1.0*w;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+seqLen+s] += 1.0*w;
                            }
                        }
                        });
                        });
                    }
                    for (int s = 0; s < refLen; ++s) { hostGapOp[gn][2*seqLen*n+s] /= refNum; hostGapEx[gn][2*seqLen*n+s] /= refNum; hostGapCl[gn][2*seqLen*n+s] /= refNum;} 
                    for (int s = 0; s < qryLen; ++s) { hostGapOp[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapEx[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapCl[gn][2*seqLen*n+seqLen+s] /= qryNum;}
                    
                    // ClustalW's method
                    
                    for (int s = 0; s < refLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5]) * 1.0 / refNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    for (int s = 0; s < qryLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    
                    if (option->gappyHorizon > 0) {
                        float gappyVertical = option->gappyVertical;
                        int gappyHorizon = option->gappyHorizon, gappyLength;
                        int rawIdx = 0, newIdx = 0;
                        int gapRef = 0, gapQry = 0;
                        while (true) {
                            if (rawIdx >= refLen) {
                                for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*i+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*gappyLength+5]/refNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].first.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHR: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapRef += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, refLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*newIdx+t] = hostFreq[gn][12*seqLen*n+6*rawIdx+t];
                                hostGapOp[gn][2*seqLen*n+newIdx] = hostGapOp[gn][2*seqLen*n+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+newIdx] = hostGapEx[gn][2*seqLen*n+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+newIdx] = hostGapCl[gn][2*seqLen*n+rawIdx];
                                ++newIdx;
                            } 
                        }
                        hostLen[gn][2*n] = newIdx; 
                        // std::cout << n << '#' << rawIdx << '/' << newIdx << '\t';
                        rawIdx = 0; newIdx = 0;
                        while (true) {
                            if (rawIdx >= qryLen) {
                                for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*(seqLen+i)+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].second.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHQ: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapQry += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, qryLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+newIdx)+t] = hostFreq[gn][12*seqLen*n+6*(seqLen+rawIdx)+t];
                                hostGapOp[gn][2*seqLen*n+seqLen+newIdx] = hostGapOp[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+seqLen+newIdx] = hostGapEx[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+seqLen+newIdx] = hostGapCl[gn][2*seqLen*n+seqLen+rawIdx];
                                ++newIdx;
                            } 
                        }
                        // if (!gappyColumns[n].first.empty()) std::cout << "No. " << n << ": Removing " << gapRef << " columns from refernce.\n";
                        // if (!gappyColumns[n].second.empty()) std::cout << "No. " << n << ": Removing " << gapQry << " columns from query.\n";
                    
                        hostLen[gn][2*n+1] = newIdx;
                        if (hostLen[gn][2*n] + gapRef != refLen) std::cout << "REF:" << hostLen[gn][2*n] << '+' << gapRef << " != " << refLen << '\n';
                        if (hostLen[gn][2*n+1] + gapQry != qryLen) std::cout << "QRY:" << hostLen[gn][2*n+1] << '+' << gapQry << " != " << qryLen << '\n';
                        
                        assert(hostLen[gn][2*n] + gapRef == refLen);
                        assert(hostLen[gn][2*n+1] + gapQry == qryLen);
                    }
                    else {
                        hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    }
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                    // gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                }
                });
                });
                    
                // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[gn][i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapEx[gn][i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapCl[gn][i] << ',';
                // std::cout << '\n';
                
                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
        
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);

                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceGapCl[gn],
                    deviceParam[gn]
                );
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();

                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                int maxAlnLen = 0;
                for (int n = 0; n <  alnPairs; ++n) {
                    if (hostAlnLen[gn][n] > maxAlnLen) maxAlnLen = hostAlnLen[gn][n];
                    // std::cout << n << ':' <<  hostAlnLen[gn][n] << '\n'; 
                }
                {
                    tbb::mutex::scoped_lock lock(memMutex);
                    util->memCheck(maxAlnLen);
                }
                if (rn % 10 == 0 && rn > 0) std::cout << rn*numBlocks << " pairs have been processed.\n";
                
                // tbb::this_task_arena::isolate( [&]{
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                // for (int n = range.begin(); n < range.end(); ++n) {
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                    if (hostAlnLen[gn][n] <= 0) {
                        cpuFallback[nIdx] = true;
                        if (util->nowProcess == 1) {
                            std::vector<int8_t> aln_old, aln;
                            std::vector<std::vector<float>> freqRef;
                            std::vector<std::vector<float>> freqQry;
                            std::vector<std::vector<float>> gapOp, gapEx, gapCl;
                            // int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                            // int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                            for (int r = 0; r < hostLen[gn][2*n]; r++) {
                                std::vector<float> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*r+f]);
                                freqRef.push_back(temp);
                            }
                            for (int q = 0; q < hostLen[gn][2*n+1]; q++) {
                                std::vector<float> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*(seqLen+q)+f]);
                                freqQry.push_back(temp);
                            }
                            for (int i = 0; i < 2; ++i) {
                                std::vector<float> temp;
                                gapOp.push_back(temp); gapEx.push_back(temp); gapCl.push_back(temp);
                            }
                            for (int r = 0; r < hostLen[gn][2*n]; ++r) {
                                gapOp[0].push_back(hostGapOp[gn][2*n*seqLen+r]);
                                gapEx[0].push_back(hostGapEx[gn][2*n*seqLen+r]);
                                gapCl[0].push_back(hostGapCl[gn][2*n*seqLen+r]);
                            }
                            for (int q = 0; q < hostLen[gn][2*n+1]; ++q) {
                                gapOp[1].push_back(hostGapOp[gn][2*n*seqLen+seqLen+q]);
                                gapEx[1].push_back(hostGapEx[gn][2*n*seqLen+seqLen+q]);
                                gapCl[1].push_back(hostGapCl[gn][2*n*seqLen+seqLen+q]);
                            }
                            
                            Talco_xdrop::Params talco_params(hostParam);
                            while (aln_old.empty()) {
                                int16_t errorType = 0;
                                Talco_xdrop::Align_freq (
                                    talco_params,
                                    freqRef,
                                    freqQry,
                                    gapOp,
                                    gapEx,
                                    gapCl,
                                    aln_old,
                                    errorType
                                );
                                // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                                if (errorType == 1) {
                                    std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                    talco_params.updateXDrop(talco_params.xdrop << 1);
                                }
                                if (errorType == 2) {
                                    std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                    talco_params.updateFLen(talco_params.fLen << 1);
                                }
                            }
                            int rIdx = 0, qIdx = 0, j = 0; 
                            // for (int j = 0; j < aln_old.size(); ++j) {
                            while (j < aln_old.size() || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty())) {
                            // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                                bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                                int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                                int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                                if (gapR || gapQ) {
                                    // std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                    if (gapRLen >= gapQLen) {
                                        for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                        for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                    }
                                    else {
                                        for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                        for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                    }
                                    if (gapR) gappyColumns[n].first.pop();
                                    if (gapQ) gappyColumns[n].second.pop();
                                }
                                else {
                                    switch ((aln_old[j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                    ++j;
                                }
                                if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                    for (j = j; j < aln_old.size(); ++j) {
                                        switch ((aln_old[j] & 0xFFFF)) {
                                            case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                            case 1: ++qIdx;         aln.push_back(1); break;
                                            case 2: ++rIdx;         aln.push_back(2); break;
                                        }
                                    }
                                }
                            }
                            // std::cout << n << '#' << rIdx << '-' << refLen << '^' << qIdx << '-' << qryLen << '\n';
                            assert(rIdx == refLen); assert(qIdx == qryLen);
                            alnCpu[nIdx] = aln;
                        }
                    }
                    else {
                        std::vector<int8_t> aln;
                        int rIdx = 0, qIdx = 0; 
                        // int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        // int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        int j = 0;
                        int oneInAln = 0;
                        while (j < hostAlnLen[gn][n] || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty()) ) {
                        // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                            bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                            int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                            int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                            if (gapR || gapQ) {
                                // if (n == 0 ) std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                if (gapRLen >= gapQLen) {
                                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                }
                                else {
                                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                }
                                if (gapR) gappyColumns[n].first.pop();
                                if (gapQ) gappyColumns[n].second.pop();
                            }
                            else {
                                switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                    case 0: ++rIdx; ++qIdx; aln.push_back(0); ++oneInAln; break;
                                    case 1: ++qIdx;         aln.push_back(1); ++oneInAln; break;
                                    case 2: ++rIdx;         aln.push_back(2); break;
                                }
                                ++j;
                            }
                            if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                for (j = j; j < hostAlnLen[gn][n]; ++j) {
                                    switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                }
                            }
                        }
                        {
                            tbb::mutex::scoped_lock lock(memMutex);
                            util->memCheck(aln.size());
                        }
                        // util->memCheck(aln.size());
                        if(rIdx != refLen || qIdx != qryLen) std::cout << n << '#' << rIdx << '-' << refLen << '^' << qIdx <<':' << oneInAln << '-' << qryLen << '\n';
                        // if (n == 0) std::cout << (aln[0] & 0xFFFF) << '\n';
                        // if (n == 0 && alnPairs == 3) {
                        //     for (int i = 0; i < aln.size(); ++i) {
                        //         if (i % 10 == 0) std::cout << '\n' << i/10 << '\t';
                        //         std::cout << (aln[i] & 0xFFFF);
                        //         
                        //     }
                        // }
                        assert(rIdx == refLen); assert(qIdx == qryLen);
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
                        for (int idx = range.begin(); idx < range.end(); ++idx) {
                            int sIdx = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[idx];
                        // for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                            int orgIdx = 0;
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {   
                            for (int k = 0; k < aln.size(); ++k) {
                                // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 2) {
                                if ((aln[k] & 0xFFFF) == 0 || (aln[k] & 0xFFFF) == 2) {
                                    util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][k] = '-';
                                }
                            }
                            // util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][n];
                            util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                            util->changeStorage(sIdx);
                        }
                        });
                        });
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
                        for (int idx = range.begin(); idx < range.end(); ++idx) {
                            int sIdx = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[idx];
                        // for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            int orgIdx = 0;
                            for (int k = 0; k < aln.size(); ++k) {
                                // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
                                if ((aln[k] & 0xFFFF) == 0 || (aln[k] & 0xFFFF) == 1) {
                                    util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                    
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][k] = '-';
                                }
                                
                            }
                            // util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
                            util->seqsLen[nodes[nIdx].second->identifier] = aln.size();  
                            util->changeStorage(sIdx);
                        }
                        });
                        });
                        for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                        }
                        tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                    }
                    
                }
                // });
                // });
            }  

            
        }
    });
    // });
    
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
        cudaFree(deviceGapCl[gn]);
        cudaDeviceSynchronize();  
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
        if (freeErr != "no error") {
            printf("ERROR: Free memory Last%s!\n", freeErr.c_str());
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
        free(hostGapCl[gn]);
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
    delete [] deviceGapCl;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostGapCl;
    // delete [] rawFreq;
    // delete [] rawGapOp;  
    // delete [] rawGapCt;
    // delete [] rawGapEn;    
    // CPU Fallback
    std::vector<int> fallbackPairs;
    for (int i = 0; i < nodes.size(); ++i) if (cpuFallback[i]) fallbackPairs.push_back(i);
    delete [] cpuFallback;
    
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    if (util->nowProcess == 0) {
        std::cout << "Bad alignments. Num of pairs: " << fallbackPairs.size() << '\n';
        for (int i = 0; i < fallbackPairs.size(); ++i) {
            int nIdx = fallbackPairs[i];
            int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
            if (util->badSequences.find(grpID) == util->badSequences.end()) {
                std::vector<std::string> temp;
                util->badSequences[grpID] = temp;
            }
            if (refNum < qryNum) {
                // proceed with qry and store ref as bad sequences
                int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                util->seqsLen[nodes[nIdx].second->identifier] = refLen;
                util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
                auto temp = nodes[nIdx].second->msaIdx;
                tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << refNum <<" seqeuences).\n";
            }
            else {
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << qryNum <<" seqeuences).\n";
            }
        }
    }
    else {
        std::cout << "CPU Fallback. Num of pairs: " << fallbackPairs.size() << '\n';
        std::vector<std::vector<int8_t>> alns;
        for (int i = 0; i < fallbackPairs.size(); ++i) alns.push_back(alnCpu[fallbackPairs[i]]);
        alnCpu.clear();
        int maxAlnLen = 0;
        for (auto aln: alns) maxAlnLen = (aln.size() > maxAlnLen) ? aln.size() : maxAlnLen;
        util->memCheck(maxAlnLen);

        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, fallbackPairs.size()), [&](tbb::blocked_range<int> range) {
        for (int n = range.begin(); n < range.end(); ++n) {
            int nIdx = fallbackPairs[n];
            auto aln = alns[n];
            // if (nodes[nIdx].first->identifier == "node_2") for (auto a: aln) std::cout << (a&0xFFFF);
            // if (nodes[nIdx].first->identifier == "node_2") std::cout << '\n';
            
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
            }
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
            // if (nodes[nIdx].first->identifier == "node_2") printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
            printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
        }
        });
        });
    }
    free(hostParam);
            
    return;
}
*/

void msaGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
      
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
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
    float** hostGapCl = new float* [gpuNum]; // gap close

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceGapCl = new float* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::mutex memMutex;
    // nowRound.store(roundGPU-1);

    // int ThreadsPerGPU = maxThreads / gpuNum;
    bool* cpuFallback = new bool[nodes.size()];
    std::vector<std::vector<int8_t>> alnCpu;
    for (int i = 0; i < nodes.size(); ++i) {
        cpuFallback[i] = false;
        if (util->nowProcess == 1) {
            std::vector<int8_t> aln;
            alnCpu.push_back(aln);
        }
        
    }
    
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] =  (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2                     * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapCl[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            
            cudaSetDevice(option->gpuIdx[gn]);
            // cudaSetDevice(1);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapCl[gn],  2 * seqLen * numBlocks * sizeof(float));

            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            
            while (nowRound < roundGPU) {
            // while (nowRound >= 0) {
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
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCl[gn][n] = 0;

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
                    int32_t numThreshold = 1000;
                    // if (nIdx == 0) std::cout << refLen << ':' << qryLen << ':' << refNum << ':' << qryNum << '\n';
                
                    float refWeight = 0.0, qryWeight = 0.0;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    if (tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                            int storage = util->seqsStorage[sIdx];
                            std::string name = util->seqsName[sIdx];
                            float w = tree->allNodes[name]->weight / refWeight * refNum;
                            tbb::this_task_arena::isolate( [&]{
                            tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                            for (int s = r.begin(); s < r.end(); ++s) {
                                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                            }
                            });
                            });
                        }
                        if (refNum >= numThreshold || qryNum > numThreshold) {
                            std::vector<std::vector<float>> freq (refLen, std::vector<float>(6,0.0));
                            tree->allNodes[nodes[nIdx].first->identifier]->msaFreq = freq;
                            for (int s = 0; s < refLen; ++s) for (int t = 0; t < 6; ++t) tree->allNodes[nodes[nIdx].first->identifier]->msaFreq[s][t] = hostFreq[gn][12*seqLen*n+6*s+t] / refNum * refWeight;
                        }
                    }
                    else {
                        for (int s = 0; s < tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.size(); ++s) {
                            for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*s+t] = tree->allNodes[nodes[nIdx].first->identifier]->msaFreq[s][t] / refWeight * refNum;
                        }
                    }
                    if (tree->allNodes[nodes[nIdx].second->identifier]->msaFreq.empty()) { 
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                            int storage = util->seqsStorage[sIdx];
                            std::string name = util->seqsName[sIdx];
                            float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                            tbb::this_task_arena::isolate( [&]{
                            tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                            for (int s = r.begin(); s < r.end(); ++s) {
                                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                            }
                            });
                            });
                        }
                        if (qryNum >= numThreshold || refNum >= numThreshold) {
                            std::vector<std::vector<float>> freq (qryLen, std::vector<float>(6,0.0));
                            tree->allNodes[nodes[nIdx].second->identifier]->msaFreq = freq;
                            for (int i = 0; i < qryLen; ++i) for (int j = 0; j < 6; ++j) tree->allNodes[nodes[nIdx].second->identifier]->msaFreq[i][j] = hostFreq[gn][12*seqLen*n+6*(seqLen+i)+j] / qryNum * qryWeight;
                        }
                    }
                    else {
                        for (int s = 0; s < tree->allNodes[nodes[nIdx].second->identifier]->msaFreq.size(); ++s) {
                            for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+s)+t] = tree->allNodes[nodes[nIdx].second->identifier]->msaFreq[s][t] / qryWeight * qryNum;
                        }
                    }
                    
                    for (int s = 0; s < tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.size(); ++s) {
                        hostGapEx[gn][2*seqLen*n+s] = hostFreq[gn][12*seqLen*n+6*s+5] / refNum;
                    }
                    for (int s = 0; s < tree->allNodes[nodes[nIdx].second->identifier]->msaFreq.size(); ++s) {
                        hostGapEx[gn][2*seqLen*n+seqLen+s] = hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] / qryNum;
                    }
                    // ClustalW's method
                    for (int s = 0; s < refLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5]) * 1.0 / refNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    for (int s = 0; s < qryLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    if (option->gappyHorizon > 0) {
                        float gappyVertical = option->gappyVertical;
                        int gappyHorizon = option->gappyHorizon, gappyLength;
                        int rawIdx = 0, newIdx = 0;
                        int gapRef = 0, gapQry = 0;
                        while (true) {
                            if (rawIdx >= refLen) {
                                for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*i+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*gappyLength+5]/refNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].first.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHR: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapRef += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, refLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*newIdx+t] = hostFreq[gn][12*seqLen*n+6*rawIdx+t];
                                hostGapOp[gn][2*seqLen*n+newIdx] = hostGapOp[gn][2*seqLen*n+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+newIdx] = hostGapEx[gn][2*seqLen*n+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+newIdx] = hostGapCl[gn][2*seqLen*n+rawIdx];
                                ++newIdx;
                            } 
                        }
                        hostLen[gn][2*n] = newIdx; 
                        // std::cout << n << '#' << rawIdx << '/' << newIdx << '\t';
                        rawIdx = 0; newIdx = 0;
                        while (true) {
                            if (rawIdx >= qryLen) {
                                for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*(seqLen+i)+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyColumns[n].second.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHQ: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapQry += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, qryLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+newIdx)+t] = hostFreq[gn][12*seqLen*n+6*(seqLen+rawIdx)+t];
                                hostGapOp[gn][2*seqLen*n+seqLen+newIdx] = hostGapOp[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+seqLen+newIdx] = hostGapEx[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+seqLen+newIdx] = hostGapCl[gn][2*seqLen*n+seqLen+rawIdx];
                                ++newIdx;
                            } 
                        }
                        
                        hostLen[gn][2*n+1] = newIdx;
                        if (hostLen[gn][2*n] + gapRef != refLen) std::cout << "REF:" << hostLen[gn][2*n] << '+' << gapRef << " != " << refLen << '\n';
                        if (hostLen[gn][2*n+1] + gapQry != qryLen) std::cout << "QRY:" << hostLen[gn][2*n+1] << '+' << gapQry << " != " << qryLen << '\n';
                        assert(hostLen[gn][2*n] + gapRef == refLen);
                        assert(hostLen[gn][2*n+1] + gapQry == qryLen);
                    }
                    else {
                        hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    }
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                    // gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                }
                });
                });
                // for (int i = 0; i < 50; ++i) std::cout << hostFreq[gn][6*(seqLen+i)+0] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapEx[gn][seqLen+i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[gn][seqLen+i] << ',';
                // std::cout << '\n';
                
                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
                
                // cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);

                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                
                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceGapCl[gn],
                    deviceParam[gn]
                );
                
                // cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                // cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyDeviceToHost);
                
                cudaDeviceSynchronize();

                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                if (rn % 10 == 0 && rn > 0) std::cout << rn*numBlocks << " pairs have been processed.\n";
                
                // tbb::this_task_arena::isolate( [&]{
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                // for (int n = range.begin(); n < range.end(); ++n) {
                for (int n = 0; n < alnPairs; ++n) {
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
                        cpuFallback[nIdx] = true;
                        if (util->nowProcess == 1) {
                            std::vector<std::vector<float>> freqRef (hostLen[gn][2*n], std::vector<float> (6, 0.0));
                            std::vector<std::vector<float>> freqQry (hostLen[gn][2*n+1], std::vector<float> (6, 0.0));
                            std::vector<std::vector<float>> gapOp (2), gapEx (2), gapCl (2);
                            
                            for (int s = 0; s < hostLen[gn][2*n]; s++)   for (int t = 0; t < 6; ++t) freqRef[s][t] = hostFreq[gn][12*seqLen*n+6*s+t];
                            for (int s = 0; s < hostLen[gn][2*n+1]; s++) for (int t = 0; t < 6; ++t) freqQry[s][t] = hostFreq[gn][12*seqLen*n+6*(seqLen+s)+t];
                            for (int r = 0; r < hostLen[gn][2*n]; ++r) {
                                gapOp[0].push_back(hostGapOp[gn][2*n*seqLen+r]);
                                gapEx[0].push_back(hostGapEx[gn][2*n*seqLen+r]);
                                gapCl[0].push_back(hostGapCl[gn][2*n*seqLen+r]);
                            }
                            for (int q = 0; q < hostLen[gn][2*n+1]; ++q) {
                                gapOp[1].push_back(hostGapOp[gn][2*n*seqLen+seqLen+q]);
                                gapEx[1].push_back(hostGapEx[gn][2*n*seqLen+seqLen+q]);
                                gapCl[1].push_back(hostGapCl[gn][2*n*seqLen+seqLen+q]);
                            }
                            std::pair<int32_t, int32_t> num = std::make_pair(refNum, qryNum);
                            Talco_xdrop::Params talco_params(hostParam);
                            while (aln_old.empty()) {
                                int16_t errorType = 0;
                                Talco_xdrop::Align_freq (
                                    talco_params,
                                    freqRef,
                                    freqQry,
                                    gapOp,
                                    gapEx,
                                    gapCl,
                                    num,
                                    aln_old,
                                    errorType
                                );
                                // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                                if (errorType == 1) {
                                    std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                    talco_params.updateXDrop(talco_params.xdrop << 1);
                                }
                                if (errorType == 2) {
                                    std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                    talco_params.updateFLen(talco_params.fLen << 1);
                                }
                            }
                            addGappyColumnsBack(aln_old, aln, gappyColumns[n], debugIdx);
                            assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
                            printf("CPU fallback on No. %d (%s), Alignment Length: %lu\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
                        }
                    }
                    else {
                        for (int j = 0; j < hostAlnLen[gn][n]; ++j) aln_old.push_back(hostAln[gn][n*2*seqLen+j]);
                        addGappyColumnsBack(aln_old, aln, gappyColumns[n], debugIdx);
                        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
                    }
                    // Update alignment & frequency
                    if (!aln.empty()) {
                        {
                            tbb::mutex::scoped_lock lock(memMutex);
                            util->memCheck(aln.size());
                        }
                        updateAlignment(tree, nodes[nIdx], util, aln);
                        if (!tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
                            updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
                        }
                    }
                }
                // });
                // });
               
            }  

            
        }
    });
    // });
    
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
        cudaFree(deviceGapCl[gn]);
        cudaDeviceSynchronize();  
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
        if (freeErr != "no error") {
            printf("ERROR: Free memory Last%s!\n", freeErr.c_str());
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
        free(hostGapCl[gn]);
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
    delete [] deviceGapCl;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostGapCl;
    // CPU Fallback
    std::vector<int> fallbackPairs;
    for (int i = 0; i < nodes.size(); ++i) if (cpuFallback[i]) fallbackPairs.push_back(i);
    delete [] cpuFallback;
    
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    if (util->nowProcess == 0) {
        std::cout << "Bad alignments. Num of pairs: " << fallbackPairs.size() << '\n';
        for (int i = 0; i < fallbackPairs.size(); ++i) {
            int nIdx = fallbackPairs[i];
            int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
            if (util->badSequences.find(grpID) == util->badSequences.end()) {
                std::vector<std::string> temp;
                util->badSequences[grpID] = temp;
            }
            if (refNum < qryNum) {
                // proceed with qry and store ref as bad sequences
                int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                util->seqsLen[nodes[nIdx].second->identifier] = refLen;
                util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
                auto temp = nodes[nIdx].second->msaIdx;
                tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << refNum <<" seqeuences).\n";
            }
            else {
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << qryNum <<" seqeuences).\n";
            }
        }
    }
    free(hostParam);
            
    return;
}

void msaCpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int32_t seqLen = util->memLen;
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;

    tbb::mutex memMutex;
    
    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
        // allocate memory
        float* hostFreq =  (float*)malloc(12 * seqLen * sizeof(float));
        float* hostGapOp = (float*)malloc( 2 * seqLen * sizeof(float));
        float* hostGapEx = (float*)malloc( 2 * seqLen * sizeof(float));
        float* hostGapCl = (float*)malloc( 2 * seqLen * sizeof(float));
        std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>> gappyColumns;
        // initialize
        for (int n = 0; n < 12*seqLen; ++n) hostFreq[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapOp[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapEx[n] = 0;
        for (int n = 0; n <  2*seqLen; ++n) hostGapCl[n] = 0;
        
        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        int32_t newRef = refLen, newQry = qryLen;
        float refWeight = 0.0, qryWeight = 0.0;
        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        
        calculateProfileFreq(hostFreq, hostGapOp, tree, nodes[nIdx], util, seqLen, param);
        if (option->gappyHorizon > 0) removeGappyColumns(hostFreq, hostGapOp, tree, nodes[nIdx], util, option, gappyColumns, newRef, newQry, seqLen);
        
        // Start alignment
        std::vector<int8_t> aln_old, aln;
        std::vector<std::vector<float>> freqRef (newRef, std::vector<float>(6, 0.0));
        std::vector<std::vector<float>> freqQry (newQry, std::vector<float>(6, 0.0));
        std::vector<std::vector<float>> gapOp (2), gapEx (2), gapCl (2);    
        for (int s = 0; s < newRef; s++) for (int t = 0; t < 6; ++t) freqRef[s][t] = hostFreq[6*s+t];
        for (int s = 0; s < newQry; s++) for (int t = 0; t < 6; ++t) freqQry[s][t] = hostFreq[6*(seqLen+s)+t];     
        for (int r = 0; r < newRef; ++r) {
            gapOp[0].push_back(hostGapOp[r]);
            gapEx[0].push_back(freqRef[r][5] / refNum);
            gapCl[0].push_back(hostGapCl[r]);
        }
        for (int q = 0; q < newQry; ++q) {
            gapOp[1].push_back(hostGapOp[seqLen+q]);
            gapEx[1].push_back(freqQry[q][5] / qryNum);
            gapCl[1].push_back(hostGapCl[seqLen+q]);
        }                 
        std::pair<int32_t, int32_t> num = std::make_pair(refNum, qryNum);
        Talco_xdrop::Params talco_params(hostParam);
        while (aln_old.empty()) {
            int16_t errorType = 0;
            Talco_xdrop::Align_freq (
                talco_params,
                freqRef,
                freqQry,
                gapOp,
                gapEx,
                gapCl,
                num,
                aln_old,
                errorType
            );
            if (errorType == 1) {
                std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                talco_params.updateXDrop(talco_params.xdrop << 1);
            }
            if (errorType == 2) {
                std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                talco_params.updateFLen(talco_params.fLen << 1);
            }
        }
        std::pair<int, int> debugIdx;
        addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx);
        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
        {
            tbb::mutex::scoped_lock lock(memMutex);
            util->memCheck(aln.size());
        }
        {
            tbb::mutex::scoped_lock lock(memMutex);
            updateAlignment(tree, nodes[nIdx], util, aln);
        }
        if (!tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
            updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
        }   
        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
        free(hostFreq);
        free(hostGapOp);
        free(hostGapEx);
        free(hostGapCl);
    }    
    });
    });
    free(hostParam);
    return;
}

void updateNode(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util) {
    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        if (n.first->is_leaf() && n.first->msaIdx.empty()) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.first->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.first->identifier]->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->msaIdx.empty()) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.second->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.second->identifier]->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        // if (nodes.size() < 5) {
        //     std::cout << n.first->identifier << '(' << util->seqsLen[n.first->identifier] << ')' <<  tree->allNodes[n.first->identifier]->msaIdx.size() << '\t';
        //     std::cout << n.second->identifier << '(' << util->seqsLen[n.second->identifier] << ')'<< tree->allNodes[n.second->identifier]->msaIdx.size() << '\n';
        // }
    }
}

void calculateProfileFreq(float* hostFreq, float* hostGapOp, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, int32_t seqLen, Params& param) {
    int32_t refLen = util->seqsLen[nodes.first->identifier];
    int32_t qryLen = util->seqsLen[nodes.second->identifier];
    int32_t refNum = tree->allNodes[nodes.first->identifier]->msaIdx.size();
    int32_t qryNum = tree->allNodes[nodes.second->identifier]->msaIdx.size();
    int32_t numThreshold = 1000;
    float refWeight = 0.0, qryWeight = 0.0;
    for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
    for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
    if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
        for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            std::string name = util->seqsName[sIdx];
            float w = tree->allNodes[name]->weight / refWeight * refNum;
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
            for (int s = r.begin(); s < r.end(); ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[6*s+0]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[6*s+1]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[6*s+2]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[6*s+3]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[6*s+4]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[6*s+5]+=1.0*w;
            }
            });
            });
        }
        if (refNum >= numThreshold || qryNum > numThreshold) {
            std::vector<std::vector<float>> freq (refLen, std::vector<float>(6,0.0));
            tree->allNodes[nodes.first->identifier]->msaFreq = freq;
            for (int s = 0; s < refLen; ++s) for (int t = 0; t < 6; ++t) tree->allNodes[nodes.first->identifier]->msaFreq[s][t] = hostFreq[6*s+t] / refNum * refWeight;
        }
    }
    else {
        for (int s = 0; s < tree->allNodes[nodes.first->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < 6; ++t) hostFreq[6*s+t] = tree->allNodes[nodes.first->identifier]->msaFreq[s][t] / refWeight * refNum;
        }
    }
    if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) { 
        for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            std::string name = util->seqsName[sIdx];
            float w = tree->allNodes[name]->weight / qryWeight * qryNum;
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
            for (int s = r.begin(); s < r.end(); ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[6*(seqLen+s)+0]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[6*(seqLen+s)+1]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[6*(seqLen+s)+2]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[6*(seqLen+s)+3]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[6*(seqLen+s)+4]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[6*(seqLen+s)+5]+=1.0*w;
            }
            });
            });
        }
        if (qryNum >= numThreshold || refNum >= numThreshold) {
            std::vector<std::vector<float>> freq (qryLen, std::vector<float>(6,0.0));
            tree->allNodes[nodes.second->identifier]->msaFreq = freq;
            for (int i = 0; i < qryLen; ++i) for (int j = 0; j < 6; ++j) tree->allNodes[nodes.second->identifier]->msaFreq[i][j] = hostFreq[6*(seqLen+i)+j] / qryNum * qryWeight;
        }
    }
    else {
        for (int s = 0; s < tree->allNodes[nodes.second->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < 6; ++t) hostFreq[6*(seqLen+s)+t] = tree->allNodes[nodes.second->identifier]->msaFreq[s][t] / qryWeight * qryNum;
        }
    } 
    // calculate gap-opne penalty (Clustalw's method)
    for (int s = 0; s < refLen; ++s) {
        if (hostFreq[6*s+5] > 0) {
            hostGapOp[s] = param.gapOpen * 0.3 * ((refNum-hostFreq[6*s+5]) * 1.0 / refNum);
        }
        else {
            int backSites = 8;
            int distance_from_gap = 0;
            int increPenalty = false;
            for (int d = s-1; d > max(s-backSites, 0); --d) {
                ++distance_from_gap;
                if (hostFreq[6*d+5] > 0) {
                    increPenalty = true;
                    break;
                }
            }
            hostGapOp[s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
        }
    }
    for (int s = 0; s < qryLen; ++s) {
        if (hostFreq[6*(seqLen+s)+5] > 0) {
            hostGapOp[seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[6*(seqLen+s)+5]) * 1.0 / qryNum);
        }
        else {
            int backSites = 8;
            int distance_from_gap = 0;
            int increPenalty = false;
            for (int d = s-1; d > max(s-backSites, 0); --d) {
                ++distance_from_gap;
                if (hostFreq[6*(seqLen+d)+5] > 0) {
                    increPenalty = true;
                    break;
                }
            }
            hostGapOp[seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
        }
    }
    return;
}

void removeGappyColumns(float* hostFreq, float* hostGapOp, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, int32_t& newRef, int32_t& newQry, int32_t seqLen) {
    float gappyVertical = option->gappyVertical;
    int gappyHorizon = option->gappyHorizon, gappyLength;
    int rawIdx = 0, newIdx = 0;
    int gapRef = 0, gapQry = 0;
    int32_t refLen = util->seqsLen[nodes.first->identifier];
    int32_t qryLen = util->seqsLen[nodes.second->identifier];
    float refNum_f = 0, qryNum_f = 0;
    for (int t = 0; t < 6; ++t) refNum_f += hostFreq[t]; 
    for (int t = 0; t < 6; ++t) qryNum_f += hostFreq[6*seqLen+t];       
    int32_t refNum = (util->nowProcess < 2) ? tree->allNodes[nodes.first->identifier]->msaIdx.size() : static_cast<int32_t>(round(refNum_f));
    int32_t qryNum = (util->nowProcess < 2) ? tree->allNodes[nodes.second->identifier]->msaIdx.size(): static_cast<int32_t>(round(qryNum_f));
    
    while (true) {
        if (rawIdx >= refLen) {
            for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[6*i+j] = 0;
            break;
        }
        for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) if (hostFreq[6*gappyLength+5]/refNum <= gappyVertical) break;
        if (gappyLength - rawIdx >= gappyHorizon) {
            gappyColumns.first.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
            gapRef += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, refLen); ++rawIdx) {
            for (int t = 0; t < 6; ++t) hostFreq[6*newIdx+t] = hostFreq[6*rawIdx+t];
            hostGapOp[newIdx] = hostGapOp[rawIdx]; 
            // hostGapEx[newIdx] = hostGapEx[rawIdx]; 
            // hostGapCl[newIdx] = hostGapCl[rawIdx];
            ++newIdx;
        } 
    }
    newRef = newIdx;
    rawIdx = 0; newIdx = 0;
    while (true) {
        if (rawIdx >= qryLen) {
            for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[6*(seqLen+i)+j] = 0;
            break;
        }
        for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) if (hostFreq[6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
        if (gappyLength - rawIdx >= gappyHorizon) {
            gappyColumns.second.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
            // if (alnPairs < 2) std::cout << "PUSHQ: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
            gapQry += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
        for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, qryLen); ++rawIdx) {
            for (int t = 0; t < 6; ++t) hostFreq[6*(seqLen+newIdx)+t] = hostFreq[6*(seqLen+rawIdx)+t];
            hostGapOp[seqLen+newIdx] = hostGapOp[seqLen+rawIdx]; 
            // hostGapEx[seqLen+newIdx] = hostGapEx[seqLen+rawIdx]; 
            // hostGapCl[seqLen+newIdx] = hostGapCl[seqLen+rawIdx];
            ++newIdx;
        } 
    }
    newQry = newIdx;
    if (newRef + gapRef != refLen) std::cout << "REF:" << newRef << '+' << gapRef << " != " << refLen << '\n';
    if (newQry + gapQry != qryLen) std::cout << "QRY:" << newQry << '+' << gapQry << " != " << qryLen << '\n';
    assert(newRef + gapRef == refLen); assert(newQry + gapQry == qryLen);
    return;
}

void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx) {
    int rIdx = 0, qIdx = 0, j = 0; 
    while (j < aln_old.size() || (!gappyColumns.first.empty() || !gappyColumns.second.empty())) {
        bool gapR = (gappyColumns.first.empty())  ? false : (rIdx == gappyColumns.first.front().first);
        bool gapQ = (gappyColumns.second.empty()) ? false : (qIdx == gappyColumns.second.front().first);
        int gapRLen = (!gapR) ? 0 : gappyColumns.first.front().second;
        int gapQLen = (!gapQ) ? 0 : gappyColumns.second.front().second;
        if (gapR || gapQ) {
            if (gapRLen >= gapQLen) {
                for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
            }
            else {
                for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
            }
            if (gapR) gappyColumns.first.pop();
            if (gapQ) gappyColumns.second.pop();
        }
        else {
            switch ((aln_old[j] & 0xFFFF)) {
                case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                case 1: ++qIdx;         aln.push_back(1); break;
                case 2: ++rIdx;         aln.push_back(2); break;
            }
            ++j;
        }
        if (gappyColumns.first.empty() && gappyColumns.second.empty()) {
            for (j = j; j < aln_old.size(); ++j) {
                switch ((aln_old[j] & 0xFFFF)) {
                    case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                    case 1: ++qIdx;         aln.push_back(1); break;
                    case 2: ++rIdx;         aln.push_back(2); break;
                }
            }
        }
    }
    debugIdx.first = rIdx; debugIdx.second = qIdx;
    return;
}

void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln) {
    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.first->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = tree->allNodes[nodes.first->identifier]->msaIdx[idx];  
    // for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
        int orgIdx = 0;
        int storeFrom = util->seqsStorage[sIdx];
        int storeTo = 1 - storeFrom;
        for (int k = 0; k < aln.size(); ++k) {
            if (aln[k] == 0 || aln[k] == 2) {
                util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                orgIdx++;
            }
            else {
                util->alnStorage[storeTo][sIdx][k] = '-';
            }
        }
        util->seqsLen[nodes.first->identifier] = aln.size();
        util->changeStorage(sIdx);
    }
    });
    });
    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.second->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
       int sIdx = tree->allNodes[nodes.second->identifier]->msaIdx[idx];
    // for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
        int storeFrom = util->seqsStorage[sIdx];
        int storeTo = 1 - util->seqsStorage[sIdx];
        int orgIdx = 0;
        for (int k = 0; k < aln.size(); ++k) {
            // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
            if ((aln[k] & 0xFFFF) == 0 || (aln[k] & 0xFFFF) == 1) {
                util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                orgIdx++;
            }
            else {
                util->alnStorage[storeTo][sIdx][k] = '-';
            }
        }
        // util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
        util->seqsLen[nodes.second->identifier] = aln.size();  
        util->changeStorage(sIdx);
    }
    });
    });
    for (auto q: tree->allNodes[nodes.second->identifier]->msaIdx) {
        tree->allNodes[nodes.first->identifier]->msaIdx.push_back(q);
    }
    tree->allNodes[nodes.second->identifier]->msaIdx.clear();
    return;
}

void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx) {
    int rIdx = 0, qIdx = 0;
    std::vector<std::vector<float>> mergeFreq(aln.size(), std::vector<float>(6,0.0));
    for (int j = 0; j < aln.size(); ++j) {
        if (aln[j] == 0) {
            for (int k = 0; k < 6; ++k) mergeFreq[j][k] = 
            (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k] + 
            tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k]);
            ++rIdx; ++qIdx;
        }
        else if (aln[j] == 1) {
            for (int k = 0; k < 5; ++k) mergeFreq[j][k] = 
            tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k];
            mergeFreq[j][5] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][5] + 1.0 * refWeight);
            ++qIdx;
        }
        else if (aln[j] == 2) {
            for (int k = 0; k < 5; ++k) mergeFreq[j][k] = 
            tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
            mergeFreq[j][5] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][5] + 1.0 * qryWeight);
            ++rIdx;
        }
    }  
    debugIdx.first = rIdx; debugIdx.second = qIdx;
    tree->allNodes[nodes.first->identifier]->msaFreq.clear();
    tree->allNodes[nodes.second->identifier]->msaFreq.clear();
    tree->allNodes[nodes.first->identifier]->msaFreq = mergeFreq;
    return;
}


__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result) {
    int32_t seqNum     = seqInfo[0]; 
    int32_t seqLen     = seqInfo[1]; 
    int32_t numBlocks  = seqInfo[2]; 
    int32_t numThreads = seqInfo[3]; 
    int32_t match      = seqInfo[4]; 
    int32_t mismatch   = seqInfo[5]; 
    int32_t gap        = seqInfo[6]; 

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;
    const int threadNum = 512;
    // int tidx = bx*bs+tx;
    if (bx < numBlocks) {
        __shared__ int64_t columnScore [threadNum];
        if (tx < numThreads) {
            columnScore[tx] = 0;
        }
        __syncthreads();
        for (int l = bx; l < seqLen; l=l+gs) {
            for (int i = 0; i < seqNum-1; ++i) {
                if (tx < numThreads) {
                    for (int j = i + 1 + tx; j < seqNum; j=j+bs) {
                        if (seqs[i*seqLen+l] == 'N' || seqs[j*seqLen+l] == 'N' || 
                            seqs[i*seqLen+l] == 'n' || seqs[j*seqLen+l] == 'n' ||
                           (seqs[i*seqLen+l] == '-' && seqs[j*seqLen+l] == '-' )) {
                            columnScore[tx] += 0;
                        }
                        else if (seqs[i*seqLen+l] == '-' || seqs[j*seqLen+l] == '-') {
                            columnScore[tx] += gap;
                        }
                        else if (seqs[i*seqLen+l] == seqs[j*seqLen+l]) {
                            columnScore[tx] += match;
                        }
                        else {
                            columnScore[tx] += mismatch;
                        }
                    }
                }
                __syncthreads();
            }
        }
        for (uint32_t r = threadNum/2; r > 0; r >>= 1) {
            if (tx < r) {
                columnScore[tx]   = columnScore[tx] + columnScore[tx+r];
            }
            __syncthreads();
        }
        if (tx == 0) result[bx] = columnScore[0];
    }

}

double getSPScore_gpu(msa::utility* util, Params& param) {
    // auto scoreStart = std::chrono::high_resolution_clock::now();
    int numBlocks = 1024;
    int blockSize = 512;
    size_t seqNum = 0;
    size_t seqLen = 0;
    if (!util->alnSeqs.empty()) {
        seqLen = util->alnSeqs.begin()->second.size();
        seqNum = util->alnSeqs.size();
        
    }
    else {
        while (util->alnStorage[util->seqsStorage[0]][0][seqLen] != 0) {
            ++seqLen;
        }
        seqNum = util->memNum;
    }
    char*    hostSeqs = (char*)malloc(seqLen * seqNum * sizeof(char));
    int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
    int64_t* hostResult = (int64_t*)malloc(numBlocks * sizeof(int64_t));
    

    
    
    // int seqCount = 0;
    if (util->alnSeqs.empty()) {
        printf("(Num, Len) - (%lu, %lu)\n", seqNum, seqLen);
        for (int i = 0; i < seqNum; ++i) {
            int storage = util->seqsStorage[i];
            for (int j = 0; j < seqLen; ++j) {
                hostSeqs[i*seqLen+j] = util->alnStorage[storage][i][j];
            }
        }
    }
    else {
        seqLen = util->alnSeqs.begin()->second.size();
        printf("(Num, Len) - (%lu, %lu)\n", seqNum, seqLen);
        int i = 0;
        for (auto s: util->alnSeqs) {
            // std::cout << i << ':' << s.first << '\n';
            // std::cout << s.second << '\n'; 
            for (int j = 0; j < seqLen; ++j) {
                hostSeqs[i*seqLen+j] = s.second[j];
            }
            ++i;
        }
    }
    // for (int j = 0; j < seqLen*seqNum; ++j) { 
    //     if (j%seqLen < alignment[seqCount].size()) {
    //         hostSeqs[j] = util->seqs[seqCount][j%seqLen];
    //     }
    //     else hostSeqs[j] = 0;
    //     if (j%seqLen == seqLen-1) ++seqCount;
    // }
    

    hostSeqInfo[0] = seqNum;
    hostSeqInfo[1] = seqLen;
    hostSeqInfo[2] = numBlocks;
    hostSeqInfo[3] = blockSize;
    hostSeqInfo[4] = param.scoringMatrix[0][0];
    hostSeqInfo[5] = param.scoringMatrix[0][1];
    hostSeqInfo[6] = param.gapExtend;
    for (int i = 0; i < numBlocks; ++i) hostResult[i] = 0;

    char*    deviceSeqs;
    int32_t* deviceSeqInfo;
    int64_t* deviceResult;

    cudaMalloc((void**)&deviceSeqs, seqLen * seqNum * sizeof(char));
    cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
    cudaMalloc((void**)&deviceResult, numBlocks * sizeof(int64_t));

    cudaMemcpy(deviceSeqs, hostSeqs, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceResult, hostResult, numBlocks * sizeof(int64_t), cudaMemcpyHostToDevice);

    calSPScore<<<numBlocks, blockSize>>>(
        deviceSeqs, 
        deviceSeqInfo,
        deviceResult
    );
    cudaDeviceSynchronize();
    cudaMemcpy(hostResult, deviceResult, numBlocks * sizeof(int64_t), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    std::string aerr = cudaGetErrorString(cudaGetLastError());
    if (aerr != "no error") {
        printf("ERROR: After kernel %s!\n", aerr.c_str());
        exit(1);
    }
    int64_t gpuScore = 0;
    for (int i = 0; i < numBlocks; ++i) {
        // std::cout << "Block: "<<i<<", score: " << hostResult[i] << '\n';
        gpuScore += hostResult[i];
    }
    double score = static_cast<double>(gpuScore)/static_cast<double>(seqNum*(seqNum-1)/2);
    // std::cout << "GPU score: " << score << '\n';
    // auto scoreEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds scoreTime = scoreEnd - scoreStart;
    
    // printf("GPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
    
    cudaFree(deviceSeqs);
    cudaFree(deviceSeqInfo);
    cudaFree(deviceResult);
    free(hostSeqs);
    free(hostSeqInfo);
    free(hostResult);
    return score;
} 
