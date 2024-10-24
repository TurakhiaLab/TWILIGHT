#ifndef MSAGPU_HPP
#include "msa-gpu.cuh"
#endif

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
        if (option->cpuOnly || m.size() < 1000 || util->nowProcess == 2) msaCpu(T, m, util, option, param);
        else if (level < 5)                                              msaGpu_s(T, m, util, option, param);
        else                                                             msaGpu(T, m, util, option, param);
        // msaGpu(T, m, util, option, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (option->printDetail) {
            if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
            else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
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
            if (!current->children[0]->msaIdx.empty()) {
                T->allNodes[p.first]->msaIdx = current->children[0]->msaIdx;
                util->seqsLen[p.first] = util->seqsLen[current->children[0]->identifier];
                break;
            }
            current = current->children[0];
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

void alignSubtreesGpu (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param) {
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
    if (option->printDetail) std::cout << "Align sub-subtrees, total " << type1Aln.size() << " pairs.\n"; 
    if (option->cpuOnly || type1Aln.size() < 1000) createOverlapAlnCpu(T, type1Aln, util, option, param);
    else                                           createOverlapAlnGpu(T, type1Aln, util, option, param);
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
                alignGrpToGrp_freq<<<numBlocks, blockSize>>>(
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

void msaGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
      
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
    float**  deviceParam = new float* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::spin_rw_mutex memMutex;
    tbb::spin_rw_mutex fallbackMutex;
    
    std::atomic<uint64_t> kernelTime, copyTime;
    kernelTime.store(0);
    copyTime.store(0);

    // int ThreadsPerGPU = maxThreads / gpuNum;
    std::vector<int> fallbackPairs;
    
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
                    if (option->gappyHorizon > 0 && option->gappyVertical < 1) {
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
                auto copyStart = std::chrono::high_resolution_clock::now();
                
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * alnPairs * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
                // cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * alnPairs * sizeof(float), cudaMemcpyHostToDevice);
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
                    deviceGapCl[gn],
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
                

                // for (int n = 0; n < alnPairs; ++n) std::cout << hostAlnLen[gn][n] << ',';
                // std::cout << '\n';


                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                if (rn % 10 == 0 && rn > 0 && option->printDetail) std::cout << rn*numBlocks << " pairs processed.\n";
                
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
                        {
                            tbb::spin_rw_mutex::scoped_lock fallbackLock;
                            fallbackLock.acquire(fallbackMutex);
                            fallbackPairs.push_back(nIdx);
                            fallbackLock.release();
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
                            tbb::spin_rw_mutex::scoped_lock memLock;
                            memLock.acquire(memMutex);
                            util->memCheck(aln.size(), option);
                            memLock.release();
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
    
    if (option->printDetail) std::cout << "Average kernel Time: " << kernelTime/1000000/gpuNum << "ms\n";
    if (option->printDetail) std::cout << "Average copy Time: " << copyTime/1000000/gpuNum << "ms\n";
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
        // cudaDeviceSynchronize();  
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
    
    
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    int totalSeqs = 0;
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        if (util->badSequences.find(grpID) == util->badSequences.end()) {
            util->badSequences[grpID] = std::vector<std::string> (0);
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
            totalSeqs += refNum;
        }
        else {
            util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
            // std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << qryNum <<" seqeuences).\n";
            totalSeqs += qryNum;
        }
    }
    if (option->printDetail) {
        if (fallbackPairs.size() == 1 && totalSeqs == 1) std::cout << "Deferring 1 pair (1 sequence).\n";
        else if (fallbackPairs.size() == 1 && totalSeqs > 1) printf("Deferring 1 pair (%d sequences).\n", totalSeqs); 
        else printf("Deferring %lu pair (%d sequences).\n", fallbackPairs.size(), totalSeqs); 
    }

    free(hostParam);
            
    return;
}

void msaGpu_s(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    updateNode(tree, nodes, util);
    int numBlocks = 1024; 
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
    float** hostWeight = new float* [gpuNum];

    char** deviceSeqs = new char* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
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
                for (int n = 0; n <  avgSeqNum[gn] * numBlocks; ++n) hostWeight[gn][n] = 0;
                // Calculate Frequency
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> r) {
                for (int n = r.begin(); n < r.end(); ++n) {  
                // for (int n = 0; n < alnPairs; ++n) {  
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
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostSeqs[gn][(refStart+idx)*seqLen+s] = 'N';
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                             {hostSeqs[gn][(refStart+idx)*seqLen+s] = '-'; gapRatioRef[s] += w;}
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
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostSeqs[gn][(qryStart+idx)*seqLen+s] = 'N';
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                             {hostSeqs[gn][(qryStart+idx)*seqLen+s] = '-'; gapRatioQry[s] += w;}
                        }
                        });
                        });
                    }
                    // ClustalW's method
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
                                if (gapRatioRef[s] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
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
                                if (gapRatioQry[s] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
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

                // for (int i = 0; i < 50; ++i) std::cout << hostSeqs[gn][i];
                // std::cout << '\n';    
                // for (int i = 0; i < 50; ++i) std::cout << hostSeqs[gn][seqLen+i];
                // std::cout << '\n';
                
                // for (int i = 0; i < 50; ++i) std::cout << hostWeight[gn][i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[gn][seqLen+i] << ',';
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
                        {
                            tbb::spin_rw_mutex::scoped_lock lock(memMutex);
                            util->memCheck(aln.size(), option);
                            updateAlignment(tree, nodes[nIdx], util, aln);
                        }
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

    if (option->printDetail) std::cout << "Average kernel Time: " << kernelTime/1000000/gpuNum << "ms\n";
    if (option->printDetail) std::cout << "Average copy Time: " << copyTime/1000000/gpuNum << "ms\n";
    
    
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
    delete [] deviceWeight;
    delete [] hostSeqs;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostWeight;

    // CPU Fallback
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    int totalSeqs = 0;
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        if (util->badSequences.find(grpID) == util->badSequences.end()) {
            util->badSequences[grpID] = std::vector<std::string> (0);
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
            totalSeqs += refNum;
        }
        else {
            util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
            totalSeqs += qryNum;
        }
    }
    if (option->printDetail) {
        if (fallbackPairs.size() == 1 && totalSeqs == 1) std::cout << "Deferring 1 pair (1 sequence).\n";
        else if (fallbackPairs.size() == 1 && totalSeqs > 1) printf("Deferring 1 pair (%d sequences).\n", totalSeqs); 
        else printf("Deferring %lu pair (%d sequences).\n", fallbackPairs.size(), totalSeqs); 
    }
    free(hostParam);
    return;
}
