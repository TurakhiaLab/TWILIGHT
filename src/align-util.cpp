#ifndef ALNUTIL_HPP
#include "align-util.hpp"
#endif


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
    }
}

void fallback2cpu(std::vector<int>& fallbackPairs, Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option) {
    int totalSeqs = 0;
    std::sort(fallbackPairs.begin(), fallbackPairs.end());
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        // if (refNum < qryNum) {
        if ((refNum < qryNum) || (refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]])) {
            bool fallback = !(refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]]) || option->noFilter;
            if (fallback) {
                if (util->badSequences.find(grpID) == util->badSequences.end()) {
                    util->badSequences[grpID] = std::vector<std::string> (0);
                }
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                if (refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]]) {
                    util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]] = false;
                }
            }
            
            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            util->seqsLen[nodes[nIdx].second->identifier] = refLen;
            util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
            auto temp = nodes[nIdx].second->msaIdx;
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
            auto tempFreq = nodes[nIdx].second->msaFreq;
            tree->allNodes[nodes[nIdx].second->identifier]->msaFreq = nodes[nIdx].first->msaFreq;
            tree->allNodes[nodes[nIdx].first->identifier]->msaFreq = tempFreq;
            totalSeqs += refNum;
        }
        else {
            // std::cout << "Deferring " << nodes[nIdx].second->identifier << '\n';
            bool fallback = !(qryNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]]) || option->noFilter;
            if (fallback) {
                if (util->badSequences.find(grpID) == util->badSequences.end()) {
                    util->badSequences[grpID] = std::vector<std::string> (0);
                }
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                if (qryNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]]) {
                    util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]] = false;
                }
            }
            totalSeqs += qryNum;
        }
    }
    if (option->printDetail) {
        if (fallbackPairs.size() == 1 && totalSeqs == 1) std::cout << "Deferring/excluding 1 pair (1 sequence).\n";
        else if (fallbackPairs.size() == 1 && totalSeqs > 1) printf("Deferring/excluding 1 pair (%d sequences).\n", totalSeqs); 
        else printf("Deferring/excluding %lu pair (%d sequences).\n", fallbackPairs.size(), totalSeqs); 
    }
    return;
}

void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln) {
    if (util->nowProcess == 0) {
        if (aln.size() != util->seqsLen[nodes.first->identifier]) {
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
        }
        if (aln.size() != util->seqsLen[nodes.second->identifier]) {
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
        }
        util->seqsLen[nodes.first->identifier] = aln.size();
        util->seqsLen[nodes.second->identifier] = aln.size();  
        for (auto q: tree->allNodes[nodes.second->identifier]->msaIdx) {
            tree->allNodes[nodes.first->identifier]->msaIdx.push_back(q);
        }
        tree->allNodes[nodes.second->identifier]->msaIdx.clear();
    }
    else if (util->nowProcess == 2) {
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
    if (tree->allNodes[nodes.first->identifier]->msaFreq.empty() || tree->allNodes[nodes.second->identifier]->msaFreq.empty()) return;
    int m_size = tree->allNodes[nodes.first->identifier]->msaFreq[0].size();
    std::vector<std::vector<float>> mergeFreq(aln.size(), std::vector<float>(m_size,0.0));
    if (util->nowProcess < 2) {
        for (int j = 0; j < aln.size(); ++j) {
            if (aln[j] == 0) {
                for (int k = 0; k < m_size; ++k) mergeFreq[j][k] = 
                (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k] + 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k]);
                ++rIdx; ++qIdx;
            }
            else if (aln[j] == 1) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k];
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][m_size-1] + 1.0 * refWeight);
                ++qIdx;
            }
            else if (aln[j] == 2) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][m_size-1] + 1.0 * qryWeight);
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
        for (int t = 0; t < m_size; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < m_size; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        int32_t refNum = static_cast<int32_t>(round(refNum_f));
        int32_t qryNum = static_cast<int32_t>(round(qryNum_f));
        for (int j = 0; j < aln.size(); ++j) {
            if (aln[j] == 0) {
                for (int k = 0; k < m_size; ++k) mergeFreq[j][k] = 
                (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k] + 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k]);
                ++rIdx; ++qIdx;
            }
            else if (aln[j] == 1) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k];
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][m_size-1] + 1.0 * refNum);
                ++qIdx;
            }
            else if (aln[j] == 2) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][m_size-1] + 1.0 * qryNum);
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
                        children[0]->branchLength += node->branchLength;
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

double calColumnSimilarity(Tree* tree, Node* node, msa::utility* util, Params& param) {
    int len = util->seqsLen[node->identifier];
    std::vector<std::vector<float>> freq;
    int seqNum = tree->allNodes[node->identifier]->msaIdx.size();
    float totalWeight = 0.0;
    for (auto sIdx: tree->allNodes[node->identifier]->msaIdx)  totalWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
    if (tree->allNodes[node->identifier]->msaFreq.empty()) {
        freq = std::vector<std::vector<float>> (len, std::vector<float>(6,0)); 
        for (auto sIdx: tree->allNodes[node->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            std::string name = util->seqsName[sIdx];
            float w = tree->allNodes[name]->weight / totalWeight * seqNum;
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
            for (int s = r.begin(); s < r.end(); ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freq[s][0]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freq[s][1]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freq[s][2]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freq[s][3]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == '-')                                              freq[s][5]+=1.0*w;
                else                                                                                             freq[s][4]+=1.0*w;
            }
            });
            });
        }
    }
    else {
        freq = tree->allNodes[node->identifier]->msaFreq;
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            for (int t = 0; t < 6; ++t) freq[s][t] = freq[s][t] / totalWeight * seqNum;
        }
        });
        });
    }

    std::vector<double> spscore (len, 0);
    tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        for (int l = 0; l < 5; l++) spscore[s] += param.scoringMatrix[l][l] * (freq[s][l] * (freq[s][l] - 1) / 2);
        for (int l = 0; l < 5; l++) spscore[s] += param.gapExtend * (freq[s][l] * (freq[s][5]));
        for (int l = 0; l < 5; l++) {
            for (int m = 0; m < 5; m++) {
                if (m != l) spscore[s] += param.scoringMatrix[l][m] * (freq[s][m] * freq[s][l]);
            }
        }  
    }
    });
    double normScore = 0;
    for (int i = 0; i < len; i++) normScore += spscore[i];
    float matchAvg = 0;
    for (int i = 0; i < 4; ++i) matchAvg += param.scoringMatrix[i][i];
    matchAvg /= 4.0;
    normScore /= ((seqNum * (seqNum-1)) / 2);
    normScore /= len;
    normScore /= matchAvg;
    return normScore;
}   

void calculateProfileFreq(float* profile, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t profileLen, int32_t profileSize, std::pair<int, int> startPos) {
    bool placement = (option->alnMode == 2 && option->nowProcess == 0);
    int32_t refLen = (!placement) ? util->seqsLen[nodes.first->identifier] : nodes.first->msaFreq.size();
    int32_t qryLen = (!placement) ? util->seqsLen[nodes.second->identifier] : nodes.second->msaFreq.size();    
    int32_t refNum = (!placement) ? tree->allNodes[nodes.first->identifier]->msaIdx.size()  : nodes.first->msa.size();
    int32_t qryNum = (!placement) ? tree->allNodes[nodes.second->identifier]->msaIdx.size() : 1;
    
    if (util->nowProcess < 2) {
        float refWeight = 0.0, qryWeight = 0.0;
        int32_t numThreshold = SEQNUM_TH;
        // for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        bool storeFreq = (((refNum >= numThreshold || qryNum > numThreshold) || 
                           (!tree->allNodes[nodes.first->identifier]->msaFreq.empty() || !tree->allNodes[nodes.second->identifier]->msaFreq.empty())) &&
                           (startPos.first == 0 && startPos.second == 0));
        storeFreq = (storeFreq & (util->nowProcess != 1));

        if (util->nowProcess != 1) {
            if (!placement)                  for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
            else if (option->treeFile != "") for (auto sName: tree->allNodes[nodes.first->identifier]->msa) refWeight += tree->allNodes[sName]->weight;
            else refWeight = static_cast<float>(refNum);
            if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
                for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx) { 
                    int storage = util->seqsStorage[sIdx];
                    std::string name = util->seqsName[sIdx];
                    float w = tree->allNodes[name]->weight / refWeight * refNum;
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(startPos.first, refLen), [&](tbb::blocked_range<int> r) {
                    for (int s = r.begin(); s < r.end(); ++s) {
                        int t = s - startPos.first;
                        int letterIndex = letterIdx(option->type, toupper(util->alnStorage[storage][sIdx][s]));
                        profile[profileSize*t+letterIndex] += 1.0 * w;
                    }
                    });
                    });
                }
                if (storeFreq) {
                    std::vector<std::vector<float>> freq (refLen, std::vector<float>(profileSize,0.0));
                    tree->allNodes[nodes.first->identifier]->msaFreq = freq;
                    for (int s = 0; s < refLen; ++s) for (int t = 0; t < profileSize; ++t) tree->allNodes[nodes.first->identifier]->msaFreq[s][t] = profile[profileSize*s+t] / refNum * refWeight;
                }
            }
            else {
                for (int s = startPos.first; s < int(tree->allNodes[nodes.first->identifier]->msaFreq.size()); ++s) {
                    int t = s - startPos.first;
                    for (int v = 0; v < profileSize; ++v) profile[profileSize*t+v] = tree->allNodes[nodes.first->identifier]->msaFreq[s][v] / refWeight * refNum;
                }
            }
        }

        if (!placement)                  for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        else if (option->treeFile != "") qryWeight = nodes.second->weight;
        else                             qryWeight = static_cast<float>(qryNum);
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) { 
            for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) { 
                int storage = util->seqsStorage[sIdx];
                std::string name = util->seqsName[sIdx];
                float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(startPos.second, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    int t = s - startPos.second;
                    int letterIndex = letterIdx(option->type, toupper(util->alnStorage[storage][sIdx][s]));
                    profile[profileSize*(profileLen+t)+letterIndex]+=1.0*w;
                }
                });
                });
            }
            if (storeFreq) {
                std::vector<std::vector<float>> freq (qryLen, std::vector<float>(profileSize,0.0));
                tree->allNodes[nodes.second->identifier]->msaFreq = freq;
                for (int s = 0; s < qryLen; ++s) for (int t = 0; t < profileSize; ++t) tree->allNodes[nodes.second->identifier]->msaFreq[s][t] = profile[profileSize*(profileLen+s)+t] / qryNum * qryWeight;
            }
        }
        else {
            for (int s = startPos.second; s < int(tree->allNodes[nodes.second->identifier]->msaFreq.size()); ++s) {
                int t = s - startPos.first;
                for (int v = 0; v < profileSize; ++v) profile[profileSize*(profileLen+t)+v] = tree->allNodes[nodes.second->identifier]->msaFreq[s][v] / qryWeight * qryNum;
            }
        }
        
    }
    else { 
        float refNum_f = 0, qryNum_f = 0;
        int subtreeRef = util->seqsIdx[nodes.first->identifier];
        int subtreeQry = util->seqsIdx[nodes.second->identifier]; 
        
        if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.first->identifier]->msaFreq = util->profileFreq[subtreeRef];
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.second->identifier]->msaFreq = util->profileFreq[subtreeQry];
        }
        for (int t = 0; t < profileSize; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < profileSize; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
        
        for (int s = 0; s < tree->allNodes[nodes.first->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < profileSize; ++t) profile[profileSize*s+t] = tree->allNodes[nodes.first->identifier]->msaFreq[s][t];
        }
        for (int s = 0; s < tree->allNodes[nodes.second->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < profileSize; ++t) profile[profileSize*(profileLen+s)+t] = tree->allNodes[nodes.second->identifier]->msaFreq[s][t];
        } 
    }
    return;
}


void calculatePSGOP(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t seqLen, std::pair<int32_t,int32_t> offset, std::pair<int32_t,int32_t> lens, Params& param) {
    int32_t refLen = lens.first;
    int32_t qryLen = lens.second;
    int32_t offsetf = offset.first;
    int32_t offsetg = offset.second; 
    int32_t refNum = (!(option->alnMode == 2) && (util->nowProcess == 0)) ? tree->allNodes[nodes.first->identifier]->msaIdx.size()  : nodes.first->msa.size();
    int32_t qryNum = (!(option->alnMode == 2) && (util->nowProcess == 0)) ? tree->allNodes[nodes.second->identifier]->msaIdx.size() : 1;
    int32_t profileSize = (option->type == 'n') ? 6 : 22;

    if (util->nowProcess == 1) {
        float refNum_f = 0.0, qryNum_f = 0.0;
        for (int t = 0; t < profileSize; ++t) refNum_f += hostFreq[t]; 
        for (int t = 0; t < profileSize; ++t) qryNum_f += hostFreq[profileSize*seqLen+t];
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
    }

    if (util->nowProcess == 2) {
        float refNum_f = 0, qryNum_f = 0;
        int subtreeRef = util->seqsIdx[nodes.first->identifier];
        int subtreeQry = util->seqsIdx[nodes.second->identifier]; 
        if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.first->identifier]->msaFreq = util->profileFreq[subtreeRef];
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.second->identifier]->msaFreq = util->profileFreq[subtreeQry];
        }
        for (int t = 0; t < profileSize; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < profileSize; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
    }
    // Clustalw's method
    float scale = (option->type == 'n') ? 0.5 : 1.0;
    if (option->psgop) {
        for (int s = 0; s < seqLen; ++s) {
            if (s < refLen) {
                float gapRatio = hostFreq[offsetf+profileSize*s+profileSize-1];
                if (gapRatio > 0) {
                    hostGapOp[offsetg+s] = param.gapOpen * scale * ((refNum-gapRatio)*1.0 / refNum);
                    hostGapEx[offsetg+s] = param.gapExtend * ((refNum-gapRatio)*1.0 / refNum);
                }
                else {
                    hostGapOp[offsetg+s] = param.gapOpen;
                    hostGapEx[offsetg+s] = param.gapExtend;
                }
                // if (refLen == 3304 && qryLen == 1948) std::cout << s << ':' << gapRatio << '/' << refNum << '=' << hostGapOp[offsetg+s] << '\n';
                
            }
            else {
                hostGapOp[offsetg+s] = 0.0;
                hostGapEx[offsetg+s] = 0.0;
            }
        }
        for (int s = 0; s < seqLen; ++s) {
            if (s < qryLen) {
                float gapRatio = hostFreq[offsetf+profileSize*(seqLen+s)+profileSize-1];
                if (gapRatio > 0) {
                    hostGapOp[offsetg+seqLen+s] = param.gapOpen * scale * ((qryNum-gapRatio) * 1.0 / qryNum);
                    hostGapEx[offsetg+seqLen+s] = param.gapExtend * ((qryNum-gapRatio) * 1.0 / qryNum);
                }
                else {
                    hostGapOp[offsetg+seqLen+s] = param.gapOpen;
                    hostGapEx[offsetg+seqLen+s] = param.gapExtend;
                }
            }
            else {
                hostGapOp[offsetg+seqLen+s] = 0.0;
                hostGapEx[offsetg+seqLen+s] = 0.0;
            }
        }   
    }
    else {
        for (int s = 0; s < seqLen; ++s) {
            hostGapOp[offsetg+s] = (s < refLen) ? param.gapOpen   : 0;
            hostGapEx[offsetg+s] = (s < refLen) ? param.gapExtend : 0;
        }
        for (int s = 0; s < seqLen; ++s) {
            hostGapOp[offsetg+seqLen+s] = (s < qryLen) ? param.gapOpen   : 0;
            hostGapEx[offsetg+seqLen+s] = (s < qryLen) ? param.gapExtend : 0;
        }
    }   
    return;
}


void removeGappyColumns(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, gappyColumnQueue& gappyColumns, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<int,int>& rawLens) {
    float gappyVertical = option->gappyVertical;
    int gappyHorizon = option->gappyHorizon, gappyLength;
    if ((gappyVertical == 1 || gappyHorizon < 1)) {
        lens.first = std::min(lens.first, maxProfileLen);
        lens.second = std::min(lens.second, maxProfileLen);
        rawLens.first = lens.first;
        rawLens.second = lens.second;
        return;
    }
    int32_t profileSize = (option->type == 'n') ? 6 : 22;
    int32_t rawIdx = 0, newIdx = 0;
    int32_t gapRef = 0, gapQry = 0;
    int32_t refLen = lens.first, qryLen = lens.second;
    float refNum_f = 0, qryNum_f = 0;
    for (int t = 0; t < profileSize; ++t) refNum_f += hostFreq[t]; 
    for (int t = 0; t < profileSize; ++t) qryNum_f += hostFreq[profileSize*profileLen+t];       
    int32_t refNum = (util->nowProcess < 2) ? (!((option->alnMode == 2) && (util->nowProcess == 0)) ? nodes.first->msaIdx.size()  : nodes.first->msa.size() ) : static_cast<int32_t>(round(refNum_f));
    int32_t qryNum = (util->nowProcess < 2) ? (!((option->alnMode == 2) && (util->nowProcess == 0)) ? nodes.second->msaIdx.size() : 1)                        : static_cast<int32_t>(round(qryNum_f));
    int32_t numTH = static_cast<int32_t>(1/(1-gappyVertical));
    while (true) {
        if (rawIdx >= refLen) {
            for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < profileSize; ++j) hostFreq[profileSize*i+j] = 0;
            break;
        }
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) {
            if ((hostFreq[profileSize*gappyLength+profileSize-1]+hostFreq[profileSize*gappyLength+profileSize-2])/refNum <= gappyVertical) break;
            if ((hostFreq[profileSize*gappyLength+profileSize-1]+hostFreq[profileSize*gappyLength+profileSize-2])/refNum > gappyVertical && refNum < numTH) break;
        }
        if (gappyLength - tempStart >= gappyHorizon) {
            float* gapProfile = new float[profileSize * (gappyLength - tempStart)];
            for (int i = tempStart; i < gappyLength; ++i) {
                for (int j = 0; j < profileSize; ++j) gapProfile[profileSize * (i-tempStart) + j] = hostFreq[profileSize*i+j];
            }
            gappyColumns.first.push({tempStart, gappyLength-tempStart, gapProfile}); 
            gapRef += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, refLen); ++rawIdx) {
            if (newIdx == maxProfileLen) break;
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*newIdx+t] = hostFreq[profileSize*rawIdx+t];
            ++newIdx;
            
        } 
        if (newIdx == maxProfileLen) break;
    }
    lens.first = newIdx;
    rawLens.first = rawIdx;
    rawIdx = 0; newIdx = 0;
    while (true) {
        if (rawIdx >= qryLen) {
            for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < profileSize; ++j) hostFreq[profileSize*(profileLen+i)+j] = 0;
            break;
        }
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) {
            if ((hostFreq[profileSize*(profileLen+gappyLength)+profileSize-1]+hostFreq[profileSize*(profileLen+gappyLength)+profileSize-2])/qryNum <= gappyVertical) break;
            if ((hostFreq[profileSize*(profileLen+gappyLength)+profileSize-1]+hostFreq[profileSize*(profileLen+gappyLength)+profileSize-2])/qryNum > gappyVertical && qryNum < numTH) break;
        }
        if (gappyLength - tempStart >= gappyHorizon) {
            float* gapProfile = new float[profileSize * (gappyLength - tempStart)];
            for (int i = tempStart; i < gappyLength; ++i) {
                for (int j = 0; j < profileSize; ++j) gapProfile[profileSize * (i-tempStart) + j] = hostFreq[profileSize*(profileLen+i)+j];
            }
            gappyColumns.second.push({tempStart, gappyLength-tempStart, gapProfile}); 
            // gappyColumns.second.push(std::make_pair(tempStart, gappyLength-tempStart)); 
            gapQry += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, qryLen); ++rawIdx) {
            if (newIdx == maxProfileLen) break;
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*(profileLen+newIdx)+t] = hostFreq[profileSize*(profileLen+rawIdx)+t];
            ++newIdx;
            
        } 
        if (newIdx == maxProfileLen) break;
    }
    lens.second = newIdx;
    rawLens.second = rawIdx;
    
    return;
}

void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, gappyColumnQueue& gappyColumns, std::pair<int, int>& debugIdx, float* hostParam, Params& param) {
    int rIdx = 0, qIdx = 0, j = 0; 
    int arIdx = 0, aqIdx = 0;
    bool endAln = (debugIdx.first <= 0);
    bool mergeCols = (debugIdx.second > 0);
    int maxrIdx = std::abs(debugIdx.first), maxqIdx = std::abs(debugIdx.second);
    bool preN = false;
    int offset = (aln.empty()) ? 0 : 1;
    int lastDir = aln_old.back();

    while ((j < aln_old.size() || (!gappyColumns.first.empty() || !gappyColumns.second.empty())) && (arIdx < maxrIdx || aqIdx < maxqIdx)) {
        bool gapR = (gappyColumns.first.empty())  ? false : ((rIdx == std::get<0>(gappyColumns.first.front()) - offset) &&  (arIdx < maxrIdx));
        bool gapQ = (gappyColumns.second.empty()) ? false : ((qIdx == std::get<0>(gappyColumns.second.front()) - offset) && (aqIdx < maxqIdx));
        int gapRLen = (!gapR) ? 0 : std::get<1>(gappyColumns.first.front());
        int gapQLen = (!gapQ) ? 0 : std::get<1>(gappyColumns.second.front());
        bool allNR = gapRLen < 0;
        bool allNQ = gapQLen < 0;
        gapRLen = std::abs(gapRLen);
        gapQLen = std::abs(gapQLen);
        if (gapR || gapQ) {
            if (gapR && !gapQ) {
                if (!allNR) {
                    for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                    delete [] std::get<2>(gappyColumns.first.front());
                    gappyColumns.first.pop();
                    preN = allNR;
                }
                else {
                    int delLen = 0, qTerminal = (gappyColumns.second.empty()) ? aln_old.size() : std::get<0>(gappyColumns.second.front()); 
                    for (int k = j; k < aln_old.size(); ++k) {
                        if ((aln_old[k] & 0xFFFF) != 1 || (qIdx+delLen) >= qTerminal) break;
                        ++delLen; 
                    }
                    if (delLen > gapRLen) {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        delete [] std::get<2>(gappyColumns.first.front());
                        gappyColumns.first.pop();
                        preN = false;
                    }
                    else {
                        for (int g = 0; g < delLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        for (int g = delLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                        delete [] std::get<2>(gappyColumns.first.front());
                        gappyColumns.first.pop();
                        preN = true;
                    }
                }   
            }
            else if (!gapR && gapQ) {
                if (!allNQ) {
                    for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                    delete [] std::get<2>(gappyColumns.second.front());
                    gappyColumns.second.pop();
                    preN = allNQ;
                }
                else {
                    int delLen = 0, rTerminal = (gappyColumns.first.empty()) ? aln_old.size() : std::get<0>(gappyColumns.first.front()); 
                    for (int k = j; k < aln_old.size(); ++k) {
                        if ((aln_old[k] & 0xFFFF) != 2 || (rIdx+delLen) >= rTerminal) break;
                        ++delLen; 
                    }
                    if (delLen > gapQLen) {
                        for (int g = 0; g < gapQLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        delete [] std::get<2>(gappyColumns.second.front());
                        gappyColumns.second.pop();
                        preN = false;
                    }
                    else {
                        for (int g = 0; g < delLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        for (int g = delLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                        delete [] std::get<2>(gappyColumns.second.front());
                        gappyColumns.second.pop();
                        preN = true;
                    }
                }   
            }
            else if (gapR && gapQ) {
                if (!allNR && !allNQ) {
                    if (mergeCols) {  
                        std::string consRef = "", consQry = "";
                        std::vector<int8_t> gapAln;
                        if (param.matrixSize == 5) {   
                            consRef = getConsensusDNA(std::get<2>(gappyColumns.first.front()), gapRLen);
                            consQry = getConsensusDNA(std::get<2>(gappyColumns.second.front()), gapQLen);
                            global_alignmentDNA(consRef, consQry, gapAln);
                        }
                        else {
                            consRef = getConsensusProtein(std::get<2>(gappyColumns.first.front()), gapRLen);
                            consQry = getConsensusProtein(std::get<2>(gappyColumns.second.front()), gapQLen);
                            global_alignmentProtein(consRef, consQry, gapAln);
                        }
                        
                        for (auto a: gapAln) {
                            if (a == 0) {
                                ++rIdx; ++qIdx; aln.push_back(0);
                            }
                            else if (a == 1) {
                                ++qIdx;         aln.push_back(1);
                            }
                            else if (a == 2) {
                                ++rIdx;         aln.push_back(2);
                            }
                        }
                        // for (auto a: gapAln) std::cout << (a & 0xFFFF);
                        // std::cout << '\n';
                        // if (gapRLen >= gapQLen) {
                        //     for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                        //     for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                        // }
                        // else {
                        //     for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                        //     for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                        // }
                        // std::cout << '\n';
                        
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                    }
                    delete [] std::get<2>(gappyColumns.first.front());
                    gappyColumns.first.pop();
                    delete [] std::get<2>(gappyColumns.second.front());
                    gappyColumns.second.pop();
                    preN = false;
                }
                else if (allNR && !allNQ) {
                    if (!preN) {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                        delete [] std::get<2>(gappyColumns.second.front());
                        gappyColumns.second.pop();
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        delete [] std::get<2>(gappyColumns.first.front());
                        gappyColumns.first.pop();
                    }   
                }
                else if (!allNR && allNQ) {
                    if (!preN) {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        delete [] std::get<2>(gappyColumns.first.front());
                        gappyColumns.first.pop();
                    }
                    else {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                        delete [] std::get<2>(gappyColumns.second.front());
                        gappyColumns.second.pop();
                    }
                }
                else {
                    if (gapRLen > gapQLen) {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                        std::get<0>(gappyColumns.first.front()) = rIdx;
                        std::get<1>(gappyColumns.first.front()) = gapRLen - gapQLen;
                        delete [] std::get<2>(gappyColumns.second.front());
                        gappyColumns.second.pop();
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                        std::get<0>(gappyColumns.second.front()) = qIdx;
                        std::get<1>(gappyColumns.second.front()) = gapQLen - gapRLen;
                        delete [] std::get<2>(gappyColumns.first.front());
                        gappyColumns.first.pop();
                    }
                    preN = true;
                }
            }
        }
        else {
            switch ((aln_old[j] & 0xFFFF)) {
                case 0: ++rIdx; ++qIdx; ++arIdx; ++aqIdx; aln.push_back(0); break;
                case 1: ++qIdx;         ++aqIdx;          aln.push_back(1); break;
                case 2: ++rIdx;         ++arIdx;          aln.push_back(2); break;
            }
            ++j;
        }
        if (gappyColumns.first.empty() && gappyColumns.second.empty()) {
            for (j = j; j < aln_old.size(); ++j) {
                switch ((aln_old[j] & 0xFFFF)) {
                    case 0: ++rIdx; ++qIdx; ++arIdx; ++aqIdx; aln.push_back(0); break;
                    case 1: ++qIdx;         ++aqIdx;          aln.push_back(1); break;
                    case 2: ++rIdx;         ++arIdx;          aln.push_back(2); break;
                }
            }
        }
    }
    
    if (endAln && (!gappyColumns.first.empty() || !gappyColumns.second.empty())) {
        bool gapR = !gappyColumns.first.empty();
        bool gapQ = !gappyColumns.second.empty();
        int gapRLen = (!gapR) ? 0 : std::get<1>(gappyColumns.first.front());
        int gapQLen = (!gapQ) ? 0 : std::get<1>(gappyColumns.second.front());
        if ((gapR && !gapQ)) {
            for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
            delete [] std::get<2>(gappyColumns.first.front());
            gappyColumns.first.pop();
        }
        else if ((gapQ && !gapR)) {
            for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            delete [] std::get<2>(gappyColumns.second.front());
            gappyColumns.second.pop();
        }
        else {
            if (mergeCols) {
                if (gapRLen >= gapQLen) {
                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                }
                else {
                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                }
            }
            else {
                for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            }
            delete [] std::get<2>(gappyColumns.first.front());
            gappyColumns.first.pop();
            delete [] std::get<2>(gappyColumns.second.front());
            gappyColumns.second.pop();
        }
    }
    debugIdx.first = rIdx; debugIdx.second = qIdx;
    return;
}


void mergeInsertions (Tree* tree, Node* nodeRef, std::vector<Node*>& nodes, msa::utility* util, std::vector<std::vector<int8_t>>& alnBad) {
    assert(nodes.size() == alnBad.size());
    std::vector<std::vector<int8_t>> alnPaths (alnBad.size()); 
    std::map<int, std::vector<std::pair<std::string, int>>> insertionGroup;
    std::vector<std::tuple<int, std::string, int>> insert2backbone;
    int32_t refLen = util->seqsLen[nodeRef->identifier];

    
    tbb::spin_rw_mutex mutex;
    tbb::parallel_for(tbb::blocked_range<int>(0, alnBad.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        int refIdx = 0, alnIdx = 0;
        // 0: match
        // 1: gaps 
        while (alnIdx < alnBad[s].size()) {
            if (alnBad[s][alnIdx] == 0) {
                alnPaths[s].push_back(0);
                alnIdx++;
                refIdx++;
            }
            else if (alnBad[s][alnIdx] == 2) {
                alnPaths[s].push_back(1);
                alnIdx++;
                refIdx++;
            }
            else if (alnBad[s][alnIdx] == 1) {
                int num = 0;
                while (alnBad[s][alnIdx] == 1) {
                    num++;
                    alnIdx++;
                }
                // for (int i = 0; i < num; ++i) alnPaths[s].push_back(3);
                {
                    tbb::spin_rw_mutex::scoped_lock lock(mutex, true);
                    insert2backbone.push_back({refIdx, nodes[s]->identifier, num});
                }
            }
        }
        if (refLen != refIdx) {
            std::cerr << "ERROR: Length not match.\t" << refLen << '-' << refIdx << '\n';
        }
    }
    });

    std::sort(insert2backbone.begin(), insert2backbone.end(), [](auto &a, auto &b){
        if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
        return std::get<1>(a) < std::get<1>(b);
    });
    for (auto& [pos, seqID, gap] : insert2backbone) {
        if (insertionGroup.find(pos) == insertionGroup.end()) insertionGroup.emplace(pos, std::vector<std::pair<std::string,int>>(0));
        insertionGroup[pos].push_back({seqID, gap});
    }

    for (auto it = insertionGroup.begin(); it != insertionGroup.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), [](auto &a, auto &b){
            return a.second > b.second;
        });
    }

    int32_t insLen = 0;
    for (auto it = insertionGroup.begin(); it != insertionGroup.end(); ++it) {
        insLen += it->second[0].second;
    }

    
    // Update sequences
    tbb::parallel_for(tbb::blocked_range<int>(0, alnBad.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        for (int ss = 0; ss < tree->allNodes[nodes[s]->identifier]->msaIdx.size(); ++ss) {
            int sIdx = tree->allNodes[nodes[s]->identifier]->msaIdx[ss];  
            int alnIdx = 0, newIdx = 0, orgIdx = 0;
            auto insertIt = insertionGroup.begin();
            util->memCheck((refLen+insLen), sIdx);
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - storeFrom;

            while (alnIdx < alnPaths[s].size()) {
                bool inserted = false;
                if (!insertionGroup.empty()) inserted = (alnIdx == insertIt->first);
                if (inserted) {
                    int maxGapLen = insertIt->second[0].second, gapLen = -1;
                    for (auto a: insertIt->second) {
                        if (a.first == nodes[s]->identifier) {
                            gapLen = a.second;
                            break;
                        }
                    }
                    if (gapLen == -1) {
                        for (int j = 0; j < maxGapLen; ++j) {
                            util->alnStorage[storeTo][sIdx][newIdx] = '-';
                            newIdx++;
                        }
                    }
                    else {
                        int j;
                        for (j = 0; j < gapLen; ++j) {
                            util->alnStorage[storeTo][sIdx][newIdx] = util->alnStorage[storeFrom][sIdx][orgIdx];
                            orgIdx++; newIdx++;
                        }
                        for (j = gapLen; j < maxGapLen; ++j) {
                            util->alnStorage[storeTo][sIdx][newIdx] = '-';
                            newIdx++;
                        }
                    }
                    ++insertIt;
                }
                else {
                    if (alnPaths[s][alnIdx] == 0) {
                        util->alnStorage[storeTo][sIdx][newIdx] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++; newIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][newIdx] = '-';
                        newIdx++;
                    }
                    alnIdx++;
                }
            }
            util->changeStorage(sIdx);
            if (newIdx != (refLen + insLen)) {
                std::cerr << "ERROR: Length not match after insertion.\t" << (refLen + insLen) << '-' << newIdx << '\n';
            }
            
        }
        
        
    }
    });

    auto insertIt = insertionGroup.begin();
    std::vector<int8_t> refUpdatedPath;
    int r_alnIdx = 0;
    while (r_alnIdx < refLen) {
        bool inserted = false;
        if (!insertionGroup.empty()) inserted = (r_alnIdx == insertIt->first);
        if (inserted) {
            for (int j = 0; j < insertIt->second[0].second; ++j) refUpdatedPath.push_back(1);
            ++insertIt;
        }
        else {
            refUpdatedPath.push_back(0);
            r_alnIdx++;
        }
    }

    // update reference alignment
    tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodeRef->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        int alnIdx = 0, newIdx = 0, orgIdx = 0;
        int sIdx = tree->allNodes[nodeRef->identifier]->msaIdx[s];
        util->memCheck((refLen+insLen), sIdx);
        int storeFrom = util->seqsStorage[sIdx];
        int storeTo = 1 - storeFrom;
        for (auto a: refUpdatedPath) {
            if (a == 0) {
                util->alnStorage[storeTo][sIdx][newIdx] = util->alnStorage[storeFrom][sIdx][orgIdx];
                orgIdx++; newIdx++;
            }
            else {
                util->alnStorage[storeTo][sIdx][newIdx] = '-';
                newIdx++;
            }
        }
        util->changeStorage(sIdx);
        if (newIdx != (refLen + insLen)) {
            std::cerr << "ERROR: Length not match after insertion.\t" << (refLen + insLen) << '-' << newIdx << '\n';
        }
        util->seqsLen[nodeRef->identifier] = newIdx;
    }
    });

    for (int i = 0; i < nodes.size(); ++i) {
        if (!alnBad[i].empty()) {
            for (auto q: tree->allNodes[nodes[i]->identifier]->msaIdx) {
                tree->allNodes[nodeRef->identifier]->msaIdx.push_back(q);
            }
            tree->allNodes[nodes[i]->identifier]->msaIdx.clear();
        }
    }
    tree->allNodes[nodeRef->identifier]->msaFreq.clear();

    return;
}

void createAlnPairs(Tree* tree, msa::utility* util, msa::option* option, std::vector<std::pair<Node*, Node*>>& alnPairs) {
    
    if (option->treeFile == "") {
        int profileSize = (option->type == 'n') ? 6 : 22;
        for (auto seq: util->placedSeqs) {
            alnPairs.push_back({tree->root, tree->allNodes[seq.first]});
        }
        for (auto seq: util->placedSeqs) {
            tree->allNodes[seq.first]->msaFreq = std::vector<std::vector<float>> (seq.second.size(), std::vector<float> (profileSize, 0.0));
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, seq.second.size()), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
                int letterIndex = letterIdx(option->type, toupper(seq.second[j]));
                tree->allNodes[seq.first]->msaFreq[j][letterIndex] += 1.0;
            }
            });
            });
        }
        return;
    }
    
    std::vector<std::string> newSeqs;
    
    for (auto node: util->placedSeqs) newSeqs.push_back(node.first);
    
    std::vector<std::string> parentNodes (newSeqs.size(), "");
    alnPairs = std::vector<std::pair<Node*, Node*>> (newSeqs.size(), {nullptr, nullptr});
    // find the nearest node that has aligned descendants
    tbb::parallel_for(tbb::blocked_range<int>(0, newSeqs.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        Node* node = tree->allNodes[newSeqs[s]];
        Node* current = node->parent;
        while (true) {
            if (checkAlignedDescendant(current, util)) {
                parentNodes[s] = current->identifier;
                alnPairs[s].first = current;
                alnPairs[s].second = node;
                break;
            } 
            if (current->parent == nullptr) break;
            current = current->parent;
        }
    }
    });

    
    // Calculate profile for each parent node
    int profileSize = (option->type == 'n') ? 6 : 22;
    std::unordered_set<std::string> parentNodesSet;
    for (auto p: parentNodes) {
        if (parentNodesSet.find(p) == parentNodesSet.end()) parentNodesSet.insert(p);
    }
    std::vector<std::string> parentNodesList;
    for (auto a: parentNodesSet) parentNodesList.push_back(a);

    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, parentNodesList.size()), [&](tbb::blocked_range<int> rr) {
    for (int s = rr.begin(); s < rr.end(); ++s) {
        // Remove all-gap columns
        Node* node = tree->allNodes[parentNodesList[s]];
        int backboneLen = util->backboneAln.begin()->second.size(); 
        collectAlignedDescendant(node, util, node->msa);
        std::vector<bool> allGaps (backboneLen, true);
        // tbb::this_task_arena::isolate( [&]{
        // tbb::parallel_for(tbb::blocked_range<int>(0, backboneLen), [&](tbb::blocked_range<int> r) {
        // for (int t = r.begin(); t < r.end(); ++t) {
        for (int t = 0; t < backboneLen; ++t) {
            for (auto node_id: node->msa) {
                if (util->backboneAln[node_id][t] != '-') {
                    allGaps[t] = false;
                    break;
                }
            }
        }
        // });
        // });
        std::vector<int> keptCols;
        for (int i = 0; i < allGaps.size(); ++i) if (!allGaps[i]) keptCols.push_back(i);
        
        // Calculate Profile
        node->msaFreq = std::vector<std::vector<float>>(keptCols.size(), std::vector<float>(profileSize,0.0));
        int refNum = node->msa.size();
        float refWeight = 0.0;
        for (auto node_id: node->msa) { 
            float w = tree->allNodes[node_id]->weight;
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, keptCols.size()), [&](tbb::blocked_range<int> r) {
            for (int ss = r.begin(); ss < r.end(); ++ss) {
                int t = keptCols[ss];
                int letterIndex = letterIdx(option->type, toupper(util->backboneAln[node_id][t]));
                tree->allNodes[node->identifier]->msaFreq[ss][letterIndex] += 1.0 * w;
            }
            });
            });
        }
        for (auto i: keptCols) {
            if (util->backboneAln[node->msa[0]][i] == '-') node->msaAln.push_back(1);
            else                                           node->msaAln.push_back(0);
        }
    }
    });
    });
    // Calculate profile for new sequences
    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, newSeqs.size()), [&](tbb::blocked_range<int> rr) {
    for (int s = rr.begin(); s < rr.end(); ++s) {
        Node* node = tree->allNodes[newSeqs[s]];
        int seqLen = util->placedSeqs[newSeqs[s]].size();
        float w = tree->allNodes[newSeqs[s]]->weight;
        tree->allNodes[node->identifier]->msaFreq = std::vector<std::vector<float>>(seqLen, std::vector<float>(profileSize,0.0));
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
        for (int t = r.begin(); t < r.end(); ++t) {
            int letterIndex = letterIdx(option->type, toupper(util->placedSeqs[newSeqs[s]][t]));
            tree->allNodes[node->identifier]->msaFreq[t][letterIndex] += 1.0 * w;
        }
        });
        });
    }
    });
    });
    return;
}

bool checkAlignedDescendant(Node* node, msa::utility* util) {
    std::queue<Node*> nodeQueue;
    nodeQueue.push(node);
    
    while (!nodeQueue.empty()) {
        Node* currentNode = nodeQueue.front();
        nodeQueue.pop();
        if (currentNode->is_leaf()) {
            if (util->backboneAln.find(currentNode->identifier) != util->backboneAln.end()) return true;
        }
        for (Node* child : currentNode->children) {
            nodeQueue.push(child);
        }
    }
    return false;
}

void collectAlignedDescendant(Node* node, msa::utility* util, std::vector<std::string>& descendants) {
    std::queue<Node*> nodeQueue;
    nodeQueue.push(node);
    while (!nodeQueue.empty()) {
        Node* currentNode = nodeQueue.front();
        nodeQueue.pop();
        if (currentNode->is_leaf()) {
            if (util->backboneAln.find(currentNode->identifier) != util->backboneAln.end()) descendants.push_back(currentNode->identifier);
        }
        for (Node* child : currentNode->children) {
            nodeQueue.push(child);
        }
    }
    return;
}

void mergedAlignedSeqs(Tree* tree, msa::utility* util, msa::option* option, const std::vector<std::pair<Node*, Node*>>& alnPairs) {
    
    std::vector<std::vector<int8_t>> alnPaths (alnPairs.size()); 
    std::map<std::pair<int, int>, std::unordered_set<std::string>> insertionGroup;
    std::vector<std::tuple<int, std::string, int>> insert2backbone;
    if (option->treeFile != "") {
        tbb::spin_rw_mutex mutex;
        tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs.size()), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            Node* refNode = alnPairs[s].first;
            Node* newSeq = alnPairs[s].second;
            std::string backboneSeq = util->backboneAln[refNode->msa[0]];

            int backboneIdx = 0, refIdx = 0, alnIdx = 0;
            // 0: match
            // 1: gaps 
            while (alnIdx < newSeq->msaAln.size()) {
                if (newSeq->msaAln[alnIdx] == 0 || newSeq->msaAln[alnIdx] == 2) {
                    if ((backboneSeq[backboneIdx] != '-' && refNode->msaAln[refIdx] == 0) ||
                        (backboneSeq[backboneIdx] == '-' && refNode->msaAln[refIdx] == 1)) {
                        if (newSeq->msaAln[alnIdx] == 0) alnPaths[s].push_back(0);
                        else                             alnPaths[s].push_back(1);
                        backboneIdx++; refIdx++; alnIdx++;
                    }
                    else if ((backboneSeq[backboneIdx] == '-' && refNode->msaAln[refIdx] == 0)) {
                        alnPaths[s].push_back(1);
                        backboneIdx++;
                    }
                    else {
                        std::cerr << "Something might wrong!\n";
                    }
                }
                else if (newSeq->msaAln[alnIdx] == 1) {
                    int num = 0;
                    while (newSeq->msaAln[alnIdx] == 1) {
                        num++;
                        alnIdx++;
                    }
                    // for (int i = 0; i < num; ++i) alnPaths[s].push_back(3);
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(mutex, true);
                        insert2backbone.push_back({backboneIdx, newSeq->identifier, num});
                    }
                }
            }
            while (alnPaths[s].size() < backboneSeq.size()) alnPaths[s].push_back(1);
            if (backboneSeq.size() != alnPaths[s].size()) {
                std::cerr << "ERROR: Length not match.\t" << backboneSeq.size() << '-' << alnPaths[s].size() << '\n';
            }
        }
        });

        std::sort(insert2backbone.begin(), insert2backbone.end(), [](auto &a, auto &b){
            if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
            return std::get<1>(a) < std::get<1>(b);
        });
        for (auto& [pos, seqID, gap] : insert2backbone) {
            insertionGroup[{pos, gap}].insert(seqID);
        }
    }
    else {
        tbb::spin_rw_mutex mutex;
        tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs.size()), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            Node* newSeq = alnPairs[s].second;
            int backboneIdx = 0, alnIdx = 0;
            // 0: match
            // 1: gaps 
            while (alnIdx < newSeq->msaAln.size()) {
                if (newSeq->msaAln[alnIdx] == 0 || newSeq->msaAln[alnIdx] == 2) {
                    alnPaths[s].push_back(newSeq->msaAln[alnIdx]);
                    alnIdx++;
                    backboneIdx++;
                }
                else if (newSeq->msaAln[alnIdx] == 1) {
                    int num = 0;
                    while (newSeq->msaAln[alnIdx] == 1) {
                        num++;
                        alnIdx++;
                    }
                    // for (int i = 0; i < num; ++i) alnPaths[s].push_back(3);
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(mutex, true);
                        insert2backbone.push_back({backboneIdx, newSeq->identifier, num});
                    }
                }
            }
            if (tree->root->msaFreq.size() != alnPaths[s].size()) {
                std::cerr << "ERROR: Length not match.\t" << tree->root->msaFreq.size() << " != " << alnPaths[s].size() << '\n';
            }
        }
        });
    }
    std::sort(insert2backbone.begin(), insert2backbone.end(), [](auto &a, auto &b){
        if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
        return std::get<1>(a) < std::get<1>(b);
    });
    for (auto& [pos, seqID, gap] : insert2backbone) {
        insertionGroup[{pos, gap}].insert(seqID);
    }
    
    // Update sequences
    tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        int alnIdx = 0, seqIdx = 0;
        auto insertIt = insertionGroup.begin();
        std::string updated_seq = "";
        while (alnIdx < alnPaths[s].size()) {
            if (alnIdx == insertIt->first.first) {
                if (insertIt->second.find(alnPairs[s].second->identifier) != insertIt->second.end()) {
                    for (int j = 0; j < insertIt->first.second; ++j) {
                        updated_seq.push_back(util->placedSeqs[alnPairs[s].second->identifier][seqIdx]);
                        seqIdx++;
                    }
                }
                else {
                    for (int j = 0; j < insertIt->first.second; ++j) updated_seq.push_back('-');
                }
                ++insertIt;
            }
            else {
                if (alnPaths[s][alnIdx] == 0) {
                    updated_seq.push_back(util->placedSeqs[alnPairs[s].second->identifier][seqIdx]);
                    seqIdx++;
                }
                else {
                    updated_seq.push_back('-');
                }
                alnIdx++;
            }
        }
        util->placedSeqs[alnPairs[s].second->identifier] = updated_seq;
    }
    });

    auto insertIt = insertionGroup.begin();
    std::vector<int8_t> backbonePath;
    int backboneLen = (option->treeFile == "") ? tree->root->msaFreq.size() : util->backboneAln.begin()->second.size();
    int b_alnIdx = 0;
    while (b_alnIdx < backboneLen) {
        if (b_alnIdx == insertIt->first.first) {
            for (int j = 0; j < insertIt->first.second; ++j) backbonePath.push_back(1);
            ++insertIt;
        }
        else {
            backbonePath.push_back(0);
            b_alnIdx++;
        }
    }

    // std::cout << backboneLen << " -> " << backbonePath.size() << '\n';

    if (option->treeFile == "") {
        tree->root->msaAln = backbonePath;
        return;
    }

    std::vector<std::string> backboneName;
    for (auto a: util->backboneAln) backboneName.push_back(a.first);
    tbb::parallel_for(tbb::blocked_range<int>(0, util->backboneAln.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        std::string updated_msa = "";
        std::string org_aln = util->backboneAln[backboneName[s]];
        int oIdx = 0;
        for (auto p: backbonePath) {
            if (p == 0) {
                updated_msa.push_back(org_aln[oIdx]);
                oIdx++;
            }
            else {
                updated_msa.push_back('-');
            }
        }
        util->backboneAln[backboneName[s]] = updated_msa;
    }
    });

    return;
}


void global_alignmentDNA(const std::string &seq1, const std::string &seq2, std::vector<int8_t>& alnPath) {
    int m = seq1.size();
    int n = seq2.size();
    // Score matrix
    std::vector<std::vector<float>> score_matrix(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<int8_t>> tb_matrix(m + 1, std::vector<int8_t>(n + 1, 0));
    // BLAST
    float mat = 1.0, mis = -2.0, gap = -2.0;
    // Initialize first row and column
    for (int i = 1; i <= m; i++) {score_matrix[i][0] = gap * i; tb_matrix[i][0] = 1;}
    for (int j = 1; j <= n; j++) {score_matrix[0][j] = gap * j; tb_matrix[0][j] = 2;}
    // Fill the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match_score = score_matrix[i - 1][j - 1] + ((seq1[i-1] == 'N' || seq2[j-1] == 'N') ? 0 : (seq1[i - 1] == seq2[j - 1] ? mat : mis));
            int delete_score = score_matrix[i - 1][j] + gap;
            int insert_score = score_matrix[i][j - 1] + gap;
            if (match_score >= delete_score && match_score >= insert_score) {
                score_matrix[i][j] = match_score;
                tb_matrix[i][j] = 0;
            }
            else if (delete_score >= insert_score) {
                score_matrix[i][j] = delete_score;
                tb_matrix[i][j] = 1;
            }
            else {
                score_matrix[i][j] = insert_score;
                tb_matrix[i][j] = 2;
            }
            
        }
    }
    // Traceback
    std::string aligned_seq1, aligned_seq2;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        switch (tb_matrix[i][j])
        {
        case 0:
            alnPath.push_back(0);
            i--;j--;
            break;
        case 1:
            alnPath.push_back(2);
            i--;
            break;
        case 2:
            alnPath.push_back(1);
            j--;
            break;
        }
    }
    std::reverse(alnPath.begin(), alnPath.end());
    return;
}

void global_alignmentProtein(const std::string &seq1, const std::string &seq2, std::vector<int8_t>& alnPath) {
    int m = seq1.size();
    int n = seq2.size();
    // Score matrix
    std::vector<std::vector<float>> score_matrix(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<int8_t>> tb_matrix(m + 1, std::vector<int8_t>(n + 1, 0));
    // BLAST
    float mat = 1.0, mis = -2.0, gap = -8.0;
    // Initialize first row and column
    for (int i = 1; i <= m; i++) {score_matrix[i][0] = gap * i; tb_matrix[i][0] = 1;}
    for (int j = 1; j <= n; j++) {score_matrix[0][j] = gap * j; tb_matrix[0][j] = 2;}
    // Fill the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int a = letterIdx('p', seq1[i-1]), b = letterIdx('p', seq2[j-1]);
            int match_score = score_matrix[i - 1][j - 1] + ((seq1[i-1] == 'X' || seq2[j-1] == 'X') ? 0 : BLOSUM62[a][b]);
            int delete_score = score_matrix[i - 1][j] + gap;
            int insert_score = score_matrix[i][j - 1] + gap;
            if (match_score >= delete_score && match_score >= insert_score) {
                score_matrix[i][j] = match_score;
                tb_matrix[i][j] = 0;
            }
            else if (delete_score >= insert_score) {
                score_matrix[i][j] = delete_score;
                tb_matrix[i][j] = 1;
            }
            else {
                score_matrix[i][j] = insert_score;
                tb_matrix[i][j] = 2;
            }
            
        }
    }
    // Traceback
    std::string aligned_seq1, aligned_seq2;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        switch (tb_matrix[i][j])
        {
        case 0:
            alnPath.push_back(0);
            i--;j--;
            break;
        case 1:
            alnPath.push_back(2);
            i--;
            break;
        case 2:
            alnPath.push_back(1);
            j--;
            break;
        }
    }
    std::reverse(alnPath.begin(), alnPath.end());
    return;
}


std::string getConsensusDNA(float* profile, int len) {
    const char bases[5] = {'A', 'C', 'G', 'T', 'N'};
    std::string consensus = "";
    for (int i = 0; i < len; ++i) {
        int max_idx = 0;
        float max_count = 0;
        for (int j = 0; j < 5; ++j) {
            if (profile[6*i+j] > max_count) {
                max_count = profile[6*i+j];
                max_idx = j;
            }
        }
        consensus.push_back(bases[max_idx]);
    }
    return consensus;
}

std::string getConsensusProtein(float* profile, int len) {
    const char acids[21] = 
    {'A','C','D','E','F',
     'G','H','I','K','L',
     'M','N','P','Q','R',
     'S','T','V','W','Y','X'};

    std::string consensus = "";
    for (int i = 0; i < len; ++i) {
        int max_idx = 0;
        float max_count = 0;
        for (int j = 0; j < 21; ++j) {
            if (profile[22*i+j] > max_count) {
                max_count = profile[22*i+j];
                max_idx = j;
            }
        }
        consensus.push_back(acids[max_idx]);
    }
    return consensus;
}
