#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include <chrono>
#include <functional>


void msa::progressive::getProgressivePairs(std::vector<std::pair<NodePair,int>>& alnOrder, std::stack<Node*> postStack, int grpID, int currentTask) {
    std::map<std::string, int> NodeAlnOrder;
    if (currentTask == 0) {
        while(!postStack.empty()) {
            Node* node = postStack.top();
            // Leaf or not belongs to the subtree
            if (!(node->grpID==-1 || node->grpID==grpID) || node->is_leaf()) {
                postStack.pop();
                continue;
            }
            // Collect children that belong to the subtree
            std::vector<Node*> children;
            for (int i = 0; i < node->children.size(); ++i) {
                if (node->children[i]->grpID == grpID) children.push_back(node->children[i]);
            }
            // Useless node, remove from the subtree
            if (children.empty() && node->seqsIncluded.empty()) {
                node->grpID = -2;
                postStack.pop();
                for (auto it = node->parent->children.begin(); it != node->parent->children.end(); ++it) {
                    if ((*it)->identifier == node->identifier) {
                        node->parent->children.erase(it);
                        break;
                    }
                }
                continue;
            }
            // Only one child, merge the child and the parent 
            if (children.size() == 1 && node->parent != nullptr && node->seqsIncluded.empty()) {
                if (node->parent->grpID == grpID) {
                    for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                        if (node->parent->children[chIdx]->identifier == node->identifier) {
                            node->parent->children[chIdx] = children[0];
                            children[0]->branchLength += node->branchLength;
                            children[0]->parent = node->parent;
                            break;
                        }
                    }
                    postStack.pop();
                    continue;
                }
            }
            // Pair children
            if (children.size() > 1) {
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
            }
            if (children.size() == 1 && !node->seqsIncluded.empty()) {
                int firstIdx  = (NodeAlnOrder.find(node->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[node->identifier]+1 : 0;
                int secondIdx = (NodeAlnOrder.find(node->children[0]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[node->children[0]->identifier]+1 : 0;
                int maxIdx = std::max(firstIdx, secondIdx);
                NodeAlnOrder[node->identifier] = maxIdx;
                NodeAlnOrder[node->children[0]->identifier] = maxIdx;
                alnOrder.push_back(std::make_pair(std::make_pair(node, node->children[0]), maxIdx));
            }
            NodeAlnOrder[node->identifier] = NodeAlnOrder[children[0]->identifier];
            postStack.pop();
        }
    }
    else {
        while(!postStack.empty()) {
            Node* node = postStack.top();
            // Pair children with their parent
            if (node->parent != nullptr) {
                int firstIdx  = (NodeAlnOrder.find(node->identifier)         != NodeAlnOrder.end()) ? NodeAlnOrder[node->identifier]+1 : 0;
                int secondIdx = (NodeAlnOrder.find(node->parent->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[node->parent->identifier]+1 : 0;
                int maxIdx = std::max(firstIdx, secondIdx);
                NodeAlnOrder[node->identifier] = maxIdx;
                NodeAlnOrder[node->parent->identifier] = maxIdx;
                alnOrder.push_back(std::make_pair(std::make_pair(node->parent, node), maxIdx));
            }
            postStack.pop();
        }
    }
    return;
}

void msa::progressive::scheduling(Node* root, std::vector<NodePairVec>& levels, int currentTask) {
    levels.clear(); // Ensure levels is empty
    std::stack<Node*> msaStack;
    root->collectPostOrder(msaStack);
    std::vector<std::pair<NodePair, int>> prgressivePairs;
    int grpID = root->grpID;
    getProgressivePairs(prgressivePairs, msaStack, grpID, currentTask);
    for (auto h : prgressivePairs) {
        while (levels.size() < h.second + 1) {
            std::vector<std::pair<Node *, Node *>> temp;
            levels.push_back(temp);
        }
        levels[h.second].push_back(h.first);
    }
    return;
}

void msa::progressive::updateNode(Tree* tree, NodePairVec& nodes, SequenceDB* database) {
    for (auto n: nodes) {
        if (n.first->is_leaf() && n.first->seqsIncluded.empty()) {
            n.first->seqsIncluded.push_back(database->name_map[n.first->identifier]->id);
            n.first->alnLen = database->name_map[n.first->identifier]->len;
            n.first->alnNum = 1;
            n.first->alnWeight = database->name_map[n.first->identifier]->weight;
        }
        else if (n.first->seqsIncluded.size() == 0) {
            Node* node = n.first;
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    n.first->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    n.first->seqsIncluded = node->children[childIndex]->seqsIncluded;
                    n.first->alnLen = node->children[childIndex]->alnLen;
                    n.first->alnNum = node->children[childIndex]->alnNum;
                    n.first->alnWeight = node->children[childIndex]->alnWeight;
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->seqsIncluded.empty()) {
            n.second->seqsIncluded.push_back(database->name_map[n.second->identifier]->id);
            n.second->alnLen = database->name_map[n.second->identifier]->len;
            n.second->alnNum = 1;
            n.second->alnWeight = database->name_map[n.second->identifier]->weight;
        }
        else if (n.second->seqsIncluded.size() == 0) {
            Node* node = n.second;
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    n.second->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    n.second->seqsIncluded = node->children[childIndex]->seqsIncluded;
                    n.second->alnLen = node->children[childIndex]->alnLen;
                    n.second->alnNum = node->children[childIndex]->alnNum;
                    n.second->alnWeight = node->children[childIndex]->alnWeight;
                    break;
                }
            }
        }
    }
    return;
}

void msa::progressive::progressiveAlignment(Tree *T, SequenceDB *database, Option *option, std::vector<NodePairVec>& alnPairsPerLevel, Params &param, alnFunction alignmentKernel) {
    int level = 0;
    if (option->printDetail) std::cerr << "Total " << alnPairsPerLevel.size() << " levels.\n";
    for (auto m : alnPairsPerLevel) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        updateNode(T, m, database);
        alignmentKernel(T, m, database, option, param);
        // database->debug();
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (option->printDetail) {
            if (m.size() > 1)
                std::cerr << "Level " << level+1 << ", aligned " << m.size() << " pairs in " << alnTime.count() / 1000000 << " ms\n";
            else
                std::cerr << "Level " << level+1 << ", aligned " << m.size() << " pair in " << alnTime.count() / 1000000 << " ms\n";
        }
        ++level;
    }
}

        
void msa::progressive::msaOnSubtree(Tree *T, SequenceDB *database, Option *option, Params &param, alnFunction alignmentKernel, int subtree) {
    auto progressiveStart = std::chrono::high_resolution_clock::now();
    // Scheduling
    std::vector<msa::NodePairVec> alnPairsPerLevel; 
    scheduling(T->root, alnPairsPerLevel, database->currentTask);
    auto scheduleEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scheduleTime = scheduleEnd - progressiveStart;
    if (option->printDetail) std::cerr << "Scheduling in " << scheduleTime.count() / 1000 << " us\n";

    // Start Progressive Alignment
    msa::progressive::progressiveAlignment(T, database, option, alnPairsPerLevel, param, alignmentKernel);
    // Push msa results to roots of the tree
    if (database->currentTask == 0) {
        Node* lastAligned = alnPairsPerLevel.back()[0].first;
        T->root->seqsIncluded = lastAligned->seqsIncluded;
        if (!lastAligned->msaFreq.empty()) T->root->msaFreq = lastAligned->msaFreq;
        T->root->alnLen = lastAligned->alnLen;
        T->root->alnNum = lastAligned->alnNum;
        T->root->alnWeight = lastAligned->alnWeight;
        lastAligned->seqsIncluded.clear();
        lastAligned->msaFreq.clear();
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progressiveTime = progressiveEnd - progressiveStart;
    if (database->currentTask != 2) std::cerr << "Progressive alignment (length: " << T->root->alnLen << ") in " << progressiveTime.count() / 1000000000 << " s\n";
    else std::cerr<< "Progressive alignment on " << T->allNodes.size() << " subalignments (length: " << T->root->getAlnLen(database->currentTask) << ") in " << progressiveTime.count() / 1000000 << " ms\n";

    if (database->fallback_nodes.empty())
        return;

    // Adding bad sequences back
    auto badStart = std::chrono::high_resolution_clock::now();
    database->currentTask = 1; 
    int badSeqBefore = 0;
    int badProfileBefore = 0;
    alnPairsPerLevel.clear();
    for (auto bad: database->fallback_nodes) badSeqBefore += bad->seqsIncluded.size();
    std::sort(database->fallback_nodes.begin(), database->fallback_nodes.end(), [&](Node* &a, Node* &b) {
        if (a->seqsIncluded.size() == b->seqsIncluded.size()) return a->getAlnLen(database->currentTask) > b->getAlnLen(database->currentTask);
        return a->seqsIncluded.size() > b->seqsIncluded.size();  // descending
    });
    for (auto bad: database->fallback_nodes) {
        alnPairsPerLevel.push_back(std::vector<NodePair>(1, {T->root, bad}));
    }
    std::cerr << "Realign profiles that have been deferred. Total profiles/sequences: " << database->fallback_nodes.size() << " / " << badSeqBefore << '\n';
    database->fallback_nodes.clear();
    progressiveAlignment(T, database, option, alnPairsPerLevel,param, cpu::alignmentKernel_CPU);
    // Reset currentTask
    database->currentTask = 0;
    auto badEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds badTime = badEnd - badStart;
    std::cerr << "Realigned profiles/sequences in " << badTime.count() / 1000000000 << " s\n";
    return;
}
