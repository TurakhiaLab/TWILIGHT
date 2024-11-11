#ifndef MSA_HPP
#include "msa-cpu.hpp"
#endif

void msaOnSubtree (Tree* T, msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param) {
    
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
    auto scheduleEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scheduleTime = scheduleEnd - progressiveStart;
    if (option->printDetail) std::cout << "Scheduled alignment in " <<  scheduleTime.count() / 1000 << " us\n";
    
    std::unordered_map<std::string, std::string> beforeAln;
    int level = 0;
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
    // Push msa results to roots of sub-subtrees
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
                    int maxIdx = std::max(firstIdx, secondIdx);
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
    if (option->printDetail) std::cout << "Align sub-subtrees, total " << type1Aln.size() << " pairs.\n"; 
    createOverlapAlnCpu(T, type1Aln, util, option, param);
    return;
}

void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param) {
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
            if (option->alnMode == 1) break;
            for (auto mp: mergePairs) {
                singleLevel.push_back(mp);
            }
            breakLoop = true;
        }
        if (option->merger == "transitivity" || util->nowProcess < 2) {
            transitivityMerge(T, newT, singleLevel, util, option);
        }
        else if (option->merger == "profile" && util->nowProcess == 2) {
            msaCpu(T, singleLevel, util, option, param);
        }
       
        auto roundEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds roundTime = roundEnd - roundStart;
        if (option->printDetail) {
            if (singleLevel.size() > 1) {
                std::cout << "Merged "<< singleLevel.size() << " edges in " << roundTime.count() / 1000 << " us\n";
            }
            else {
                std::cout << "Merged "<< singleLevel.size() << " edge in " << roundTime.count() / 1000 << " us\n";
            }
        }
        totalEdges += singleLevel.size();
        if (breakLoop) break;
    }
    if (option->printDetail) std::cout << "Total Edges/Levels: " << totalEdges << '/' << totalLevels << '\n';
    return;
}

void createOverlapAlnCpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    int32_t seqLen = 0;
    if (util->nowProcess < 2) {
        for (auto n: tree->allNodes) seqLen = std::max(util->seqsLen[n.first], seqLen);
    }
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
                    for (int d = s-1; d > std::max(s-backSites, 0); --d) {
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
                    for (int d = s-1; d > std::max(s-backSites, 0); --d) {
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
        addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx, option);
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

void transitivityMerge(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option) {
    
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
                for (auto id: tree->allNodes[n.first->identifier]->msa) {
                    std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[id]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
                    for (int idx = range.begin(); idx < range.end(); ++idx) {
                        int sIdx = tree->allNodes[id]->msaIdx[idx];
                    // for (auto sIdx: tree->allNodes[id]->msaIdx) {
                        util->memCheck(seqLen, sIdx);
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
                    });
                    });
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
        // int32_t seqLen = std::max(refLen, qryLen);
        for (int i = 0; i < refNum; ++i) {
            std::vector<int8_t> temp;
            refNewAln.push_back(temp);
        }
        for (int i = 0; i < qryNum; ++i) {
            std::vector<int8_t> temp;
            qryNewAln.push_back(temp);
        }
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
    for (auto n: nodes) {
        if (util->nowProcess == 2 && n.first->parent == nullptr) tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
        std::vector<std::string> tempTotal = tree->allNodes[n.first->identifier]->msa;
        std::vector<std::string> tempR = tree->allNodes[n.first->identifier]->msa;
        std::vector<std::string> tempQ = tree->allNodes[n.second->identifier]->msa;
        for (auto id: tree->allNodes[n.second->identifier]->msa) tempTotal.push_back(id);
        for (auto id: tempR) tree->allNodes[id]->msa = tempTotal;
        for (auto id: tempQ) tree->allNodes[id]->msa = tempTotal;
    }
    return;
}

void msaCpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{
    if (util->nowProcess < 2) updateNode(tree, nodes, util);
    int32_t seqLen = 0;
    if (util->nowProcess < 2) {
        for (auto n: tree->allNodes) seqLen = (util->seqsLen[n.first] > seqLen) ? util->seqsLen[n.first] : seqLen;
        seqLen *= 2;
    }
    else {
        for (auto n: nodes) {
            if (tree->allNodes[n.first->identifier]->msaAln.size() > seqLen)  seqLen = tree->allNodes[n.first->identifier]->msaAln.size();
            if (tree->allNodes[n.second->identifier]->msaAln.size() > seqLen) seqLen = tree->allNodes[n.second->identifier]->msaAln.size();
        }
    }
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;

    // tbb::mutex memMutex;

    tbb::spin_rw_mutex fallbackMutex;
    std::vector<int> fallbackPairs;

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
        int32_t newRef = refLen, newQry = qryLen;
        float refWeight = 0.0, qryWeight = 0.0;
        if (util->nowProcess < 2) for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        if (util->nowProcess < 2) for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        calculateProfileFreq(hostFreq, hostGapOp, tree, nodes[nIdx], util, seqLen, param);
        float refNum_f = 0, qryNum_f = 0;
        for (int t = 0; t < 6; ++t) refNum_f += hostFreq[t]; 
        for (int t = 0; t < 6; ++t) qryNum_f += hostFreq[6*seqLen+t];       
        int32_t refNum = (util->nowProcess < 2) ? tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() : static_cast<int32_t>(round(refNum_f));
        int32_t qryNum = (util->nowProcess < 2) ? tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size(): static_cast<int32_t>(round(qryNum_f));
        
        if (option->gappyHorizon > 0) removeGappyColumns(hostFreq, hostGapOp, tree, nodes[nIdx], util, option, gappyColumns, newRef, newQry, seqLen);
        // if (nodes.size() == 1) {
        // for (int i = 0; i < 50; ++i) std::cout << hostFreq[6*(seqLen+i)+0] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostGapEx[seqLen+i] << ',';
        // std::cout << '\n';
        // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[seqLen+i] << ',';
        // std::cout << '\n';
        // }
        
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
        std::pair<float, float> num = std::make_pair(static_cast<float>(refNum), static_cast<float>(qryNum));
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
            if (util->nowProcess == 0 && errorType != 0) {
                aln_old.clear();
                {
                    tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                    fallbackPairs.push_back(nIdx);
                }
                break;
            }
            if (errorType == 1) {
                if (option->printDetail) std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                talco_params.updateXDrop(talco_params.xdrop << 1);
            }
            if (errorType == 2) {
                if (option->printDetail) std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                talco_params.updateFLen(talco_params.fLen << 1);
            }
        }
        if (!aln_old.empty()) {
            std::pair<int, int> debugIdx;
            addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx, option);
            if (debugIdx.first != refLen || debugIdx.second != qryLen) {
                std::cout << "Name (" << nIdx << "): " << nodes[nIdx].first->identifier << '-' << nodes[nIdx].second->identifier << '\n';
                std::cout << "Len: " << debugIdx.first << '/' << refLen << '-' << debugIdx.second << '/' << qryLen << '\n';
                std::cout << "Num: " << refNum << '-' << qryNum << '\n';
            }
            assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
            updateAlignment(tree, nodes[nIdx], util, aln); 
            if (!tree->allNodes[nodes[nIdx].first->identifier]->msaFreq.empty()) {
                updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
            }
            assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
        }
        free(hostFreq);
        free(hostGapOp);
        free(hostGapEx);
        free(hostGapCl);
    }    
    });
    });
    if (util->nowProcess < 2) {
        for (auto n: tree->allNodes) {
            if (n.second->is_leaf()) {
                if (util->memLen < util->seqMemLen[util->seqsIdx[n.first]]) util->memLen = util->seqMemLen[util->seqsIdx[n.first]];
            }
        }
    }
    
    free(hostParam);
    if (fallbackPairs.empty()) {
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
                    // if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
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
                    // if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
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
    if (util->nowProcess < 2) {
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
    }
    else {
        float refNum_f = 0, qryNum_f = 0;
        int subtreeRef = tree->allNodes[nodes.first->identifier]->grpID;
        int subtreeQry = tree->allNodes[nodes.second->identifier]->grpID; 
        if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.first->identifier]->msaFreq = util->profileFreq[subtreeRef];
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.second->identifier]->msaFreq = util->profileFreq[subtreeQry];
        }
        for (int t = 0; t < 6; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < 6; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
        
        for (int s = 0; s < tree->allNodes[nodes.first->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < 6; ++t) hostFreq[6*s+t] = tree->allNodes[nodes.first->identifier]->msaFreq[s][t];
        }
        for (int s = 0; s < tree->allNodes[nodes.second->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < 6; ++t) hostFreq[6*(seqLen+s)+t] = tree->allNodes[nodes.second->identifier]->msaFreq[s][t];
        }  
             
    }
    // Clustalw's method
    for (int s = 0; s < refLen; ++s) {
        if (hostFreq[6*s+5] > 0) {
            hostGapOp[s] = param.gapOpen * 0.3 * ((refNum-hostFreq[6*s+5])*1.0 / refNum);
        }
        else {
            int backSites = 8;
            int distance_from_gap = 0;
            int increPenalty = false;
            for (int d = s-1; d > std::max(s-backSites, 0); --d) {
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
            for (int d = s-1; d > std::max(s-backSites, 0); --d) {
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

/*
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
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, refLen); ++rawIdx) {
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
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, qryLen); ++rawIdx) {
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
*/

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
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) {
            if (hostFreq[6*gappyLength+5]/refNum <= gappyVertical) break;
            if (gappyLength == rawIdx) {
                if (hostFreq[6*gappyLength+0] == 0 && 
                    hostFreq[6*gappyLength+1] == 0 &&
                    hostFreq[6*gappyLength+2] == 0 &&
                    hostFreq[6*gappyLength+3] == 0) onlyN = true;
                else onlyN = false;
            }
            if (onlyN) {
                if (hostFreq[6*gappyLength+0] > 0 || 
                    hostFreq[6*gappyLength+1] > 0 ||
                    hostFreq[6*gappyLength+2] > 0 ||
                    hostFreq[6*gappyLength+3] > 0) {
                    gappyColumns.first.push(std::make_pair(tempStart, -1*(gappyLength-tempStart))); 
                    onlyN = false;
                    tempStart = gappyLength;
                }
            }
            else {
                if (hostFreq[6*gappyLength+0] == 0 && 
                    hostFreq[6*gappyLength+1] == 0 &&
                    hostFreq[6*gappyLength+2] == 0 &&
                    hostFreq[6*gappyLength+3] == 0) {
                    gappyColumns.first.push(std::make_pair(tempStart, (gappyLength-tempStart))); 
                    onlyN = true;
                    tempStart = gappyLength;
                }
            }

        }
        if (gappyLength - tempStart >= gappyHorizon) {
            if (onlyN) gappyColumns.first.push(std::make_pair(tempStart, -1*(gappyLength-tempStart))); 
            else       gappyColumns.first.push(std::make_pair(tempStart, gappyLength-tempStart)); 
            gapRef += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, refLen); ++rawIdx) {
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
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) {
            if (hostFreq[6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
            if (gappyLength == rawIdx) {
                if (hostFreq[6*(seqLen+gappyLength)+0] == 0 && 
                    hostFreq[6*(seqLen+gappyLength)+1] == 0 &&
                    hostFreq[6*(seqLen+gappyLength)+2] == 0 &&
                    hostFreq[6*(seqLen+gappyLength)+3] == 0) onlyN = true;
                else onlyN = false;
            }
            if (onlyN) {
                if (hostFreq[6*(seqLen+gappyLength)+0] > 0 || 
                    hostFreq[6*(seqLen+gappyLength)+1] > 0 ||
                    hostFreq[6*(seqLen+gappyLength)+2] > 0 ||
                    hostFreq[6*(seqLen+gappyLength)+3] > 0) {
                    gappyColumns.second.push(std::make_pair(tempStart, -1*(gappyLength-tempStart))); 
                    onlyN = false;
                    tempStart = gappyLength;
                }
            }
            else {
                if (hostFreq[6*(seqLen+gappyLength)+0] == 0 && 
                    hostFreq[6*(seqLen+gappyLength)+1] == 0 &&
                    hostFreq[6*(seqLen+gappyLength)+2] == 0 &&
                    hostFreq[6*(seqLen+gappyLength)+3] == 0) {
                    gappyColumns.second.push(std::make_pair(tempStart, (gappyLength-tempStart))); 
                    onlyN = true;
                    tempStart = gappyLength;
                }
            }
        }
        if (gappyLength - tempStart >= gappyHorizon) {
            if (onlyN) gappyColumns.second.push(std::make_pair(tempStart, -1*(gappyLength-tempStart))); 
            else       gappyColumns.second.push(std::make_pair(tempStart, gappyLength-tempStart)); 
            gapQry += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, qryLen); ++rawIdx) {
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

/*
void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx, msa::option* option) {
    int rIdx = 0, qIdx = 0, j = 0; 
    while (j < aln_old.size() || (!gappyColumns.first.empty() || !gappyColumns.second.empty())) {
        bool gapR = (gappyColumns.first.empty())  ? false : (rIdx == gappyColumns.first.front().first);
        bool gapQ = (gappyColumns.second.empty()) ? false : (qIdx == gappyColumns.second.front().first);
        int gapRLen = (!gapR) ? 0 : gappyColumns.first.front().second;
        int gapQLen = (!gapQ) ? 0 : gappyColumns.second.front().second;
        bool allNR = gapRLen < 0;
        bool allNQ = gapQLen < 0;
        gapRLen = std::abs(gapRLen);
        gapQLen = std::abs(gapQLen);
        if (gapR || gapQ) {
            // if (option->alignGappy) {
            //     if (gapRLen >= gapQLen) {
            //         for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
            //         for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
            //     }
            //     else {
            //         for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
            //         for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
            //     }
            // }
            // else {
            for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
            for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            // }
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
*/

void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx, msa::option* option) {
    int rIdx = 0, qIdx = 0, j = 0; 
    while (j < aln_old.size() || (!gappyColumns.first.empty() || !gappyColumns.second.empty())) {
        bool gapR = (gappyColumns.first.empty())  ? false : (rIdx == gappyColumns.first.front().first);
        bool gapQ = (gappyColumns.second.empty()) ? false : (qIdx == gappyColumns.second.front().first);
        int gapRLen = (!gapR) ? 0 : gappyColumns.first.front().second;
        int gapQLen = (!gapQ) ? 0 : gappyColumns.second.front().second;
        bool allNR = gapRLen < 0;
        bool allNQ = gapQLen < 0;
        gapRLen = std::abs(gapRLen);
        gapQLen = std::abs(gapQLen);
        if (gapR && !gapQ) {
            for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
            gappyColumns.first.pop();
        }
        else if (!gapR && gapQ) {
            for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            gappyColumns.second.pop();
        }
        else if (gapR && gapQ) {
            if (!allNR && !allNQ) {
                for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                gappyColumns.first.pop();
                gappyColumns.second.pop();
            }
            else if (allNR && !allNQ) {
                for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                gappyColumns.second.pop();
                // std::cout << "Q: " << gapQLen << '\n';
            }
            else if (allNR && !allNQ) {
                for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                gappyColumns.first.pop();
                // std::cout << "R: " << gapRLen << '\n';
            }
            else {
                // std::cout << "RQRQRQ: " << std::min(gapRLen, gapQLen) << '\n';
                if (gapRLen > gapQLen) {
                    for (int g = 0; g < gapQLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                    gappyColumns.first.front().first = rIdx;
                    gappyColumns.first.front().second = gapRLen - gapQLen;
                    gappyColumns.second.pop();
                }
                else {
                    for (int g = 0; g < gapRLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                    gappyColumns.second.front().first = qIdx;
                    gappyColumns.second.front().second = gapQLen - gapRLen;
                    gappyColumns.first.pop();
                }
            }
            // for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
            // for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            // if (gapR) gappyColumns.first.pop();
            // if (gapQ) gappyColumns.second.pop();
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
    if (util->nowProcess < 2) {
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.first->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
        for (int idx = range.begin(); idx < range.end(); ++idx) {
            int sIdx = tree->allNodes[nodes.first->identifier]->msaIdx[idx];  
            util->memCheck(aln.size(), sIdx);
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
            util->changeStorage(sIdx);
        }
        });
        });
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.second->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
        for (int idx = range.begin(); idx < range.end(); ++idx) {
            int sIdx = tree->allNodes[nodes.second->identifier]->msaIdx[idx];
            util->memCheck(aln.size(), sIdx);
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
            util->changeStorage(sIdx);
        }
        });
        });
        util->seqsLen[nodes.first->identifier] = aln.size();
        util->seqsLen[nodes.second->identifier] = aln.size();  
        for (auto q: tree->allNodes[nodes.second->identifier]->msaIdx) {
            tree->allNodes[nodes.first->identifier]->msaIdx.push_back(q);
        }
        tree->allNodes[nodes.second->identifier]->msaIdx.clear();
    }
    else {
        // int seqLenR = util->seqsLen[nodes.first->identifier];
        // int seqLenQ = util->seqsLen[nodes.second->identifier];
        if (tree->allNodes[nodes.first->identifier]->msa.empty())  tree->allNodes[nodes.first->identifier]->msa.push_back(nodes.first->identifier);
        if (tree->allNodes[nodes.second->identifier]->msa.empty()) tree->allNodes[nodes.second->identifier]->msa.push_back(nodes.second->identifier);
        for (auto seqName: tree->allNodes[nodes.first->identifier]->msa) {
            std::vector<int8_t> seqAln;
            int orgIdx = 0;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 2) {
                    seqAln.push_back(tree->allNodes[seqName]->msaAln[orgIdx]);
                    orgIdx++;
                }
                else {
                    seqAln.push_back(3);
                }
            }
            assert(orgIdx == tree->allNodes[seqName]->msaAln.size());
            tree->allNodes[seqName]->msaAln = seqAln;
            util->seqsLen[seqName] = aln.size();
        }
        for (auto seqName: tree->allNodes[nodes.second->identifier]->msa) {
            std::vector<int8_t> seqAln;
            int orgIdx = 0;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 1) {
                    seqAln.push_back(tree->allNodes[seqName]->msaAln[orgIdx]);
                    orgIdx++;
                }
                else {
                    seqAln.push_back(3);
                }
            }
            assert(orgIdx == tree->allNodes[seqName]->msaAln.size());
            tree->allNodes[seqName]->msaAln = seqAln;
            util->seqsLen[seqName] = aln.size();
        }
        for (auto q: tree->allNodes[nodes.second->identifier]->msa) {
            tree->allNodes[nodes.first->identifier]->msa.push_back(q);
        }
        for (auto q: tree->allNodes[nodes.first->identifier]->msa) {
            if (q != nodes.first->identifier) tree->allNodes[q]->msa = tree->allNodes[nodes.first->identifier]->msa;
        }
    }
    return;
}

void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx) {
    int rIdx = 0, qIdx = 0;
    std::vector<std::vector<float>> mergeFreq(aln.size(), std::vector<float>(6,0.0));
    if (util->nowProcess < 2) {
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
    }
    else {
        float refNum_f = 0, qryNum_f = 0;
        // int subtreeRef = tree->allNodes[nodes.first->identifier]->grpID;
        // int subtreeQry = tree->allNodes[nodes.second->identifier]->grpID; 
        for (int t = 0; t < 6; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < 6; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        int32_t refNum = static_cast<int32_t>(round(refNum_f));
        int32_t qryNum = static_cast<int32_t>(round(qryNum_f));
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
                mergeFreq[j][5] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][5] + 1.0 * refNum);
                ++qIdx;
            }
            else if (aln[j] == 2) {
                for (int k = 0; k < 5; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
                mergeFreq[j][5] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][5] + 1.0 * qryNum);
                ++rIdx;
            }
        }  
        debugIdx.first = rIdx; debugIdx.second = qIdx;
        for (auto node: tree->allNodes[nodes.first->identifier]->msa) {
            tree->allNodes[node]->msaFreq = mergeFreq;
            assert(tree->allNodes[node]->msaFreq.size() == tree->allNodes[node]->msaAln.size());
        }
    }
    return;
}

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID) {
    
    std::map<std::string, int> NodeAlnOrder;

    while(!postOrder.empty()) {
        Node* node = postOrder.top();
        // Leaf or not belongs to the subtree
        if (!(node->grpID==-1 || node->grpID==grpID) || node->is_leaf()) {
            postOrder.pop();
            continue;
        }
        // Collect children that belong to the subtree
        std::vector<Node*> children;
        for (int i = 0; i < node->children.size(); ++i) {
            if (node->children[i]->grpID == grpID) children.push_back(node->children[i]);
        }
        // Useless node, remove from the subtree
        if (children.empty()) {
            node->grpID = -2;
            postOrder.pop();
            for (auto it = node->parent->children.begin(); it != node->parent->children.end(); ++it) {
                if ((*it)->identifier == node->identifier) {
                    node->parent->children.erase(it);
                    break;
                }
            }
            continue;
        }
        // Only one child, merge the child and the parent 
        if (children.size() == 1 && node->parent != nullptr) {
            if (node->parent->grpID == grpID) {
                for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                    if (node->parent->children[chIdx]->identifier == node->identifier) {
                        node->parent->children[chIdx] = children[0];
                        children[0]->parent = node->parent;
                        break;
                    }
                }
                postOrder.pop();
                continue;
            }
        }
        // Pair children
        while (children.size() > 1) {
            std::vector<Node*> nodeLeft;
            for (int i = 0; i < children.size()-1; i+=2) {
                int firstIdx  = (NodeAlnOrder.find(children[i]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i]->identifier]+1 : 0;
                int secondIdx = (NodeAlnOrder.find(children[i+1]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i+1]->identifier]+1 : 0;
                int maxIdx = std::max(firstIdx, secondIdx);
                NodeAlnOrder[children[i]->identifier] = maxIdx;
                NodeAlnOrder[children[i+1]->identifier] = maxIdx;
                alnOrder.push_back(std::make_pair(std::make_pair(children[i], children[i+1]), maxIdx));
                nodeLeft.push_back(children[i]);
            }
            if (children.size()%2 == 1) nodeLeft.push_back(children.back());
            children = nodeLeft;
        }
        NodeAlnOrder[node->identifier] = NodeAlnOrder[children[0]->identifier];
        // children[0] = node;
        // std::map<int, std::vector<Node*>> childrenIdx;
        // for (int i = 0; i < children.size(); ++i) {
        //     int idx  = (NodeAlnOrder.find(children[i]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i]->identifier]+1 : 0;
        //     if (childrenIdx.find(idx) == childrenIdx.end()) childrenIdx[idx] = std::vector<Node*> (0);
        //     childrenIdx[idx].push_back(children[i]);
        // }
        // int start = childrenIdx.begin()->first, end = childrenIdx.rbegin()->first + 1;
        // for (int idx = start; idx < end; ++idx) {
        //     if (childrenIdx.find(idx+1) == childrenIdx.end()) childrenIdx[idx+1] = std::vector<Node*> (0);
        //     for (int i = 0; i < childrenIdx[idx].size()-1; i+=2) {
        //         alnOrder.push_back(std::make_pair(std::make_pair(children[i], children[i+1]), idx));
        //         childrenIdx[idx+1].push_back(children[i]);
        //     }
        //     if (childrenIdx[idx].size()%2 == 1) childrenIdx[idx+1].push_back(children.back());
        // }
        // int idx = end - 1;
        // while (true) {
        //     if (childrenIdx[idx].size() < 2) break;
        //     if (childrenIdx.find(idx+1) == childrenIdx.end()) childrenIdx[idx+1] = std::vector<Node*> (0);
        //     for (int i = 0; i < childrenIdx[idx].size()-1; i+=2) {
        //         alnOrder.push_back(std::make_pair(std::make_pair(children[i], children[i+1]), idx));
        //         childrenIdx[idx+1].push_back(children[i]);
        //     }
        //     if (childrenIdx[idx].size()%2 == 1) childrenIdx[idx+1].push_back(children.back());
        //     ++idx;
        // }
        // NodeAlnOrder[node->identifier] = idx - 1;
        postOrder.pop();
    }
}

void getPostOrderList(Node* node, std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(node); 
    Node* current; 
  
    while (!s1.empty()) {
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == current->grpID) s1.push(current->children[i]);     
        }
    } 
    return;
}