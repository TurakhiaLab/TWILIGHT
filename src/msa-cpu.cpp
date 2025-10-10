#ifndef MSA_HPP
#include "msa-cpu.hpp"
#endif

void msaOnSubtree (Tree* T, msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param) {
    
    auto progressiveStart = std::chrono::high_resolution_clock::now();
    // Schedule the alignment pairs
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
    if (option->printDetail) std::cout << "Scheduling in " <<  scheduleTime.count() / 1000 << " us\n";
    
    std::unordered_map<std::string, std::string> beforeAln;
    int level = 0;
    // Perform MSA
    if (option->printDetail) std::cout << "Total " << hier.size() << " levels.\n";
    for (auto m: hier) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        option->calSim = (option->psgopAuto && level >= 5 && !option->psgop);
        msaCpu(T, m, util, option, param);
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
    if (option->alnMode == 0) {
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
                    if (!current->children[chIdx]->msaFreq.empty()) T->allNodes[p.first]->msaFreq = std::move(current->children[chIdx]->msaFreq);
                    util->seqsLen[p.first] = util->seqsLen[current->children[chIdx]->identifier];
                    break;
                }
                current = current->children[chIdx];
            }
        }
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progessiveTime = progressiveEnd - progressiveStart;
    if (option->alnMode != 2) std::cout << "Progressive alignment in " <<  progessiveTime.count() / 1000000000 << " s\n";
    else std::cout << "Alignment in " <<  progessiveTime.count() / 1000000000 << " s\n";
    if (util->badSequences.empty()) return;
    // Adding bad sequences back
    util->nowProcess = 1;
    auto badStart = std::chrono::high_resolution_clock::now();
    int badSeqBefore = 0;
    int badProfileBefore = 0;
    hier.clear();
    hier = std::vector<std::vector<std::pair<Node *, Node *>>>(1);
    for (auto p : partition->partitionsRoot) {
        std::set<std::string> badNodes;
        if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end()) {
            auto badSeqName = util->badSequences[p.second.first->grpID];
            for (auto n: badSeqName) {
                badNodes.insert(n);
                badProfileBefore += 1;
                badSeqBefore += T->allNodes[n]->msaIdx.size();
            }   
        }
        for (auto name: badNodes) {
            hier[0].push_back(std::make_pair(T->allNodes[p.second.first->identifier], T->allNodes[name]));
        }   
    }
    util->badSequences.clear();
    std::cout << "Adding bad profiles back. Total profiles/sequences: " << badProfileBefore << " / " << badSeqBefore << '\n';
    level = 0;
    if (!hier.empty()) {
        for (auto m : hier) {
            auto alnStart = std::chrono::high_resolution_clock::now();
            msaCpu(T, m, util, option, param);
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
    std::cout << "Added bad profiles/sequences in " <<  badTime.count() / 1000000000 << " s\n";
    return;
}

void placement(Tree *T, msa::utility *util, msa::option *option, partitionInfo_t *partition, Params &param)
{

    auto progressiveStart = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<Node *, Node *>> alnPairs;

    createAlnPairs(T, util, option, alnPairs);
    auto scheduleEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scheduleTime = scheduleEnd - progressiveStart;
    
    if (option->printDetail)
        std::cout << "Create alignment pairs in " << scheduleTime.count() / 1000 << " us\n";

    
    auto alnStart = std::chrono::high_resolution_clock::now();
    option->calSim = false;
    msaCpu(T, alnPairs, util, option, param);
    mergedAlignedSeqs(T, util, option, alnPairs);
    auto alnEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    if (alnPairs.size() > 1)
        std::cout << "Aligned " << alnPairs.size() << " new sequences in " << alnTime.count() / 1000000 << " ms\n";
    else
            std::cout << "Aligned 1 new sequence in " << alnTime.count() / 1000000 << " ms\n";

    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progressiveTime = progressiveEnd - progressiveStart;
    if (option->alnMode != 2) std::cout << "Progressive alignment (length: " << util->seqsLen[T->root->identifier] << ") in " << progressiveTime.count() / 1000000000 << " s\n";
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
    // if (option->printDetail) std::cout << "Align sub-subtrees, total " << type1Aln.size() << " pairs.\n"; 
    createOverlapAlnCpu(T, type1Aln, util, option, param);
    return;
}

void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param) {
    for (auto n: newT->allNodes) {
        T->allNodes[n.first]->msa.clear();
        T->allNodes[n.first]->msa.push_back(n.first);
    }
    // Get all the edges to be merged
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
        // A greeady approach to schedule merged pairs
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
    // if (option->printDetail) std::cout << "Total Edges/Levels: " << totalEdges << '/' << totalLevels << '\n';
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
  
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
        float* hostFreq  = (float*)malloc(profileSize * 2 * seqLen * sizeof(float));
        float* hostGapOp = (float*)malloc(              2 * seqLen * sizeof(float));
        float* hostGapEx = (float*)malloc(              2 * seqLen * sizeof(float));
        float* hostGapCl = (float*)malloc(              2 * seqLen * sizeof(float));
        for (int n = 0; n < profileSize * 2 * seqLen; ++n) hostFreq[n] = 0;
        for (int n = 0; n <               2 * seqLen; ++n) hostGapOp[n] = 0;
        for (int n = 0; n <               2 * seqLen; ++n) hostGapEx[n] = 0;
        for (int n = 0; n <               2 * seqLen; ++n) hostGapCl[n] = 0;
        gappyColumnQueue gappyColumns;
        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
        int32_t newRef = refLen, newQry = qryLen;
        float refNum = 0, qryNum = 0;
        if (util->nowProcess < 2) {
            std::pair<int,int> startPos (std::make_pair(0,0));
            calculateProfileFreq(hostFreq, tree, nodes[nIdx], util, option, seqLen, profileSize, startPos);
            std::pair<int,int> lens = std::make_pair(refLen, qryLen), rawLens (std::make_pair(0, 0)), offset (std::make_pair(0, 0));
            removeGappyColumns(hostFreq, tree, nodes[nIdx], util, option, gappyColumns, seqLen, seqLen, lens, rawLens);
            newRef = lens.first, newQry = lens.second;
            calculatePSGOP(hostFreq, hostGapOp, hostGapEx, tree, nodes[nIdx], util, option, seqLen, offset, lens, param);                  
            refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        }
        else {
            int subtreeRef = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            int subtreeQry = tree->allNodes[nodes[nIdx].second->identifier]->grpID;
            for (int t = 0; t < profileSize; ++t) refNum += util->profileFreq[subtreeRef][0][t]; 
            for (int t = 0; t < profileSize; ++t) qryNum += util->profileFreq[subtreeQry][0][t]; 
            for (int s = 0; s < refLen; ++s) for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*s+t]          = util->profileFreq[subtreeRef][s][t]; 
            for (int s = 0; s < qryLen; ++s) for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*(seqLen+s)+t] = util->profileFreq[subtreeQry][s][t]; 
            for (int s = 0; s < refLen; ++s) hostGapEx[s]        = util->profileFreq[subtreeRef][s][profileSize-1] / refNum; 
            for (int s = 0; s < qryLen; ++s) hostGapEx[seqLen+s] = util->profileFreq[subtreeQry][s][profileSize-1] / qryNum; 
            std::pair<int,int> lens = std::make_pair(refLen, qryLen), rawLens (std::make_pair(0, 0)), offset (std::make_pair(0, 0));
            removeGappyColumns(hostFreq, tree, nodes[nIdx], util, option, gappyColumns, seqLen, seqLen, lens, rawLens);
            int32_t newRef = lens.first, newQry = lens.second;
            calculatePSGOP(hostFreq, hostGapOp, hostGapEx, tree, nodes[nIdx], util, option, seqLen, offset, lens, param);                  
        }
        std::vector<int8_t> aln_old, aln;
        std::vector<std::vector<float>> freqRef (newRef, std::vector<float>(profileSize, 0.0));
        std::vector<std::vector<float>> freqQry (newQry, std::vector<float>(profileSize, 0.0));
        std::vector<std::vector<float>> gapOp (2), gapEx (2);    
        for (int s = 0; s < newRef; s++) for (int t = 0; t < profileSize; ++t) freqRef[s][t] = hostFreq[profileSize*s+t];
        for (int s = 0; s < newQry; s++) for (int t = 0; t < profileSize; ++t) freqQry[s][t] = hostFreq[profileSize*(seqLen+s)+t];     
        for (int r = 0; r < newRef; ++r) {
            gapOp[0].push_back(hostGapOp[r]);
            gapEx[0].push_back(hostGapEx[r]);
        }
        for (int q = 0; q < newQry; ++q) {
            gapOp[1].push_back(hostGapOp[seqLen+q]);
            gapEx[1].push_back(hostGapEx[seqLen+q]);
        }   
        Talco_xdrop::Params* talco_params = new Talco_xdrop::Params(hostParam, param.matrixSize);
        std::pair<int32_t, int32_t> num = std::make_pair(static_cast<int32_t>(refNum), static_cast<int32_t>(qryNum));
        while (aln_old.empty()) {
            int16_t errorType = 0;
            Talco_xdrop::Align_freq (
                talco_params,
                freqRef,
                freqQry,
                gapOp,
                gapEx,
                num,
                aln_old,
                errorType
            );
            // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
            if (errorType == 1) {
                std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                talco_params->updateXDrop(2*talco_params->xdrop);
            }
            if (errorType == 2) {
                std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                talco_params->updateFLen(talco_params->fLen << 1);
            }
        }
        delete talco_params;
        int alnRef = 0, alnQry = 0;
        for (auto a: aln_old) {
            if (a == 0) {alnRef += 1; alnQry += 1;}
            if (a == 1) {alnQry += 1;}
            if (a == 2) {alnRef += 1;}
        }
        if (!option->alignGappy) alnQry *= -1;
        std::pair<int, int> debugIdx = std::make_pair(-1*alnRef,alnQry);   
        bool tempAlnGappy = option->alignGappy;
        option->alignGappy = false; 
        addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx, hostParam, param);
        option->alignGappy = tempAlnGappy;
        assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
        free(hostFreq);
        free(hostGapCl);
        free(hostGapEx);
        free(hostGapOp);        
    }
    });
   
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
        }
        int8_t refGap = (newtree->allNodes[n.first->identifier]->level == newtree->allNodes[n.second->identifier]->level) ? 2 : 1; 
        std::vector<std::vector<int8_t>> refAln, qryAln;
        std::vector<std::vector<int8_t>> refNewAln, qryNewAln;
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
    
    bool placement = (option->alnMode == 2 && option->nowProcess == 0);

    if (util->nowProcess < 2 && !placement) updateNode(tree, nodes, util);

    int32_t seqLen = 0;
    if (placement) {
        for (auto n : nodes) seqLen = std::max(seqLen, static_cast<int32_t>(std::max(n.first->msaFreq.size(), n.second->msaFreq.size())));
    }
    else if (util->nowProcess < 2) {
        for (auto n: nodes) seqLen = std::max(seqLen, std::max(util->seqsLen[n.first->identifier], util->seqsLen[n.second->identifier]));
    }
    else {
        for (auto n: nodes) seqLen = std::max(seqLen, static_cast<int32_t>(std::max(n.first->msaAln.size(), n.second->msaAln.size())));
        if (seqLen == 0) {
            for (auto p: util->profileFreq) seqLen = std::max(seqLen, static_cast<int32_t>(p.second.size()));   
            for (auto n: nodes) seqLen = std::max(seqLen, std::max(static_cast<int32_t>(tree->allNodes[n.first->identifier]->msaFreq.size()), static_cast<int32_t>(tree->allNodes[n.second->identifier]->msaFreq.size())));   
        }
    }

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
    

    tbb::spin_rw_mutex fallbackMutex;
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

    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
        // allocate memory
        float* hostFreq =  (float*)malloc(profileSize * 2 * seqLen * sizeof(float));
        float* hostGapOp = (float*)malloc(              2 * seqLen * sizeof(float));
        float* hostGapEx = (float*)malloc(              2 * seqLen * sizeof(float));
        gappyColumnQueue gappyColumns;
        // initialize
        for (int n = 0; n < profileSize * 2 * seqLen; ++n) hostFreq[n] = 0;
        for (int n = 0; n <               2 * seqLen; ++n) hostGapOp[n] = 0;
        for (int n = 0; n <               2 * seqLen; ++n) hostGapEx[n] = 0;
        int32_t refLen = (!placement) ? util->seqsLen[nodes[nIdx].first->identifier]  : nodes[nIdx].first->msaFreq.size();
        int32_t qryLen = (!placement) ? util->seqsLen[nodes[nIdx].second->identifier] : nodes[nIdx].second->msaFreq.size();
                        
        
        float refWeight = 0.0, qryWeight = 0.0;
        if (!placement) {
            if (util->nowProcess < 2) {
                if (util->nowProcess == 0) for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                else                       for (auto weight: tree->allNodes[nodes[nIdx].first->identifier]->msaFreq[0]) refWeight += weight;
            }
            if (util->nowProcess < 2) for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        }

        std::pair<int,int> startPos (std::make_pair(0,0));
        calculateProfileFreq(hostFreq, tree, nodes[nIdx], util, option, seqLen, profileSize, startPos);
        if (util->nowProcess == 1) {
            int32_t refNum = tree->allNodes[nodes[0].first->identifier]->msaIdx.size();
            for (int s = 0; s < util->seqsLen[nodes[0].first->identifier]; ++s) {
                for (int v = 0; v < profileSize; ++v)  {
                    hostFreq[profileSize*s+v] = tree->allNodes[nodes[0].first->identifier]->msaFreq[s][v]  / refWeight * refNum;
                }
            } 
        }
        std::pair<int,int> lens = std::make_pair(refLen, qryLen), rawLens (std::make_pair(0, 0)), offset (std::make_pair(0, 0));
        removeGappyColumns(hostFreq, tree, nodes[nIdx], util, option, gappyColumns, seqLen, seqLen, lens, rawLens);
        int32_t newRef = lens.first, newQry = lens.second;
        calculatePSGOP(hostFreq, hostGapOp, hostGapEx, tree, nodes[nIdx], util, option, seqLen, offset, lens, param);                  
        float refNum_f = 0, qryNum_f = 0;
        for (int t = 0; t < profileSize; ++t) refNum_f += hostFreq[t]; 
        for (int t = 0; t < profileSize; ++t) qryNum_f += hostFreq[profileSize*seqLen+t];       
        int32_t refNum = (placement) ? nodes[nIdx].first->msa.size() : ((util->nowProcess < 2) ? tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() : static_cast<int32_t>(round(refNum_f)));
        int32_t qryNum = (placement) ? 1                             : ((util->nowProcess < 2) ? tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size(): static_cast<int32_t>(round(qryNum_f)));
        
        // Start alignment
        std::vector<int8_t> aln_old, aln;
        std::vector<std::vector<float>> freqRef (newRef, std::vector<float>(profileSize, 0.0));
        std::vector<std::vector<float>> freqQry (newQry, std::vector<float>(profileSize, 0.0));
        std::vector<std::vector<float>> gapOp (2), gapEx (2);    
        for (int s = 0; s < newRef; s++) for (int t = 0; t < profileSize; ++t) freqRef[s][t] = hostFreq[profileSize*s+t];
        for (int s = 0; s < newQry; s++) for (int t = 0; t < profileSize; ++t) freqQry[s][t] = hostFreq[profileSize*(seqLen+s)+t];     
        for (int r = 0; r < newRef; ++r) {
            gapOp[0].push_back(hostGapOp[r]);
            gapEx[0].push_back(hostGapEx[r]);
        }
        for (int q = 0; q < newQry; ++q) {
            gapOp[1].push_back(hostGapOp[seqLen+q]);
            gapEx[1].push_back(hostGapEx[seqLen+q]);
        }                 
        free(hostFreq);
        free(hostGapOp);
        free(hostGapEx);

        std::pair<float, float> num = std::make_pair(static_cast<float>(refNum), static_cast<float>(qryNum));
        Talco_xdrop::Params* talco_params = new Talco_xdrop::Params(hostParam, param.matrixSize);
        if (util->nowProcess == 1) talco_params->gapBoundary = 0; 
        if (refLen == 0) for (int j = 0; j < qryLen; ++j) aln_old.push_back(1);
        if (qryLen == 0) for (int j = 0; j < refLen; ++j) aln_old.push_back(2);
        // if (util->nowProcess == 1) std::cout << nIdx << ':' << newRef << '\t' << newQry << std::endl;
        if (option->noFilter || (!util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]] && !util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]])) {
            while (aln_old.empty()) {
                int16_t errorType = 0;
                aln_old.clear();
                Talco_xdrop::Align_freq (
                    talco_params,
                    freqRef,
                    freqQry,
                    gapOp,
                    gapEx,
                    num,
                    aln_old,
                    errorType
                );
                if (util->nowProcess == 0 && errorType != 0 && !placement) {
                    aln_old.clear();
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                        fallbackPairs.push_back(nIdx);
                    }
                    break;
                }
                if (errorType == 2) {
                    if (option->printDetail) std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                    talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->fLen * 1.2) << 1, std::min(newRef, newQry)));
                }
                else if (errorType == 3) {
                    std::cout << "There might be some bugs in the code!\n";
                    exit(1);
                }
                else if (errorType == 1) {
                    talco_params->updateXDrop(static_cast<int32_t>(talco_params->xdrop * 1.2));
                    // talco_params->updateXDrop(static_cast<int32_t>(std::min(newRef, newQry) * -1 * param.gapExtend));
                    talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->xdrop * 4) << 1, std::min(newRef, newQry)));
                    if (option->printDetail) std::cout << "Updated x-drop value on No. " << nIdx << "\tNew Xdrop: " << talco_params->xdrop << '\n';

                }
            }
        }
        delete talco_params;
        if (util->nowProcess == 0 && (refNum == 1 || qryNum == 1) && !placement) {
            if (util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]] || 
                util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]]) {
                aln_old.clear();
                {
                    tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                    fallbackPairs.push_back(nIdx);
                }
            }
        }
        
        if (!aln_old.empty()) {
            int alnRef = 0, alnQry = 0;
            for (auto a: aln_old) {
                // if (nIdx == 0) std::cout << (a & 0xFFFF);
                if (a == 0) {alnRef += 1; alnQry += 1;}
                if (a == 1) {alnQry += 1;}
                if (a == 2) {alnRef += 1;}
            }
            // if (nIdx == 0) std::cout << '\n';
            if (!option->alignGappy) alnQry *= -1;
            std::pair<int, int> debugIdx = std::make_pair(-1*alnRef,alnQry);
            addGappyColumnsBack(aln_old, aln, gappyColumns, debugIdx, hostParam, param);
            if (debugIdx.first != refLen || debugIdx.second != qryLen) {
                std::cout << "Name (" << nIdx << "): " << nodes[nIdx].first->identifier << '-' << nodes[nIdx].second->identifier << '\n';
                std::cout << "Len: " << debugIdx.first << '/' << refLen << '-' << debugIdx.second << '/' << qryLen << '\n';
                std::cout << "Num: " << refNum << '-' << qryNum << '\n';
            }
            // assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
            if (util->nowProcess != 1) {
                if (!placement) {
                    updateAlignment(tree, nodes[nIdx], util, aln); 
                    updateFrequency(tree, nodes[nIdx], util, aln, refWeight, qryWeight, debugIdx);
                }
                else {
                    nodes[nIdx].second->msaAln = std::move(aln);
                }
            }
            else {
                alnBad[nIdx] = std::move(aln);
            }
            // assert(debugIdx.first == refLen); assert(debugIdx.second == qryLen);
            if (option->calSim && util->nowProcess < 2) {
                double colSim = calColumnSimilarity(tree, nodes[nIdx].first, util, param);
                if (colSim < option->divTH) {
                    option->redo = true;
                    return;
                }
            }
        }
        // }

    }    
    });
    if (util->nowProcess == 1) {
        std::vector<Node*> badNodes;
        for (auto n: nodes) badNodes.push_back(n.second);
        mergeInsertions (tree, nodes[0].first, badNodes, util, alnBad);
    }
    if (util->nowProcess < 2 && !placement) {
        for (auto n: tree->allNodes) {
            if (n.second->is_leaf()) {
                if (util->memLen < util->seqMemLen[util->seqsIdx[n.first]]) util->memLen = util->seqMemLen[util->seqsIdx[n.first]];
            }
        }
    }
    
    free(hostParam);
    if (fallbackPairs.empty()) return;
    fallback2cpu(fallbackPairs, tree, nodes, util, option);
    return;
}
