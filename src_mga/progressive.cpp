#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <chrono>
#include <functional>
#include <stack>
#include <vector>
#include <map>
#include <algorithm>

namespace mga {
namespace progressive {

void getProgressivePairs(std::vector<std::pair<NodePair, int>>& alnOrder, std::stack<Node*>& postStack, int grpID, int mode) {
    std::map<std::string, int> nodeAlnOrder;
    if (mode == 0) {
        while (!postStack.empty()) {
            Node* node = postStack.top();
            postStack.pop();

            if (!(node->grpID == -1 || node->grpID == grpID) || node->is_leaf()) {
                continue;
            }

            std::vector<Node*> children;
            for (const auto& child_ptr : node->children) {
                if (child_ptr->grpID == grpID) {
                    children.push_back(child_ptr.get());
                }
            }

            if (children.empty() && node->blockId == 0) {
                node->grpID = -2;
                if (node->parent) {
                    auto& parent_children = node->parent->children;
                    parent_children.erase(std::remove_if(parent_children.begin(), parent_children.end(),
                                                       [node](const std::unique_ptr<Node>& p) { return p.get() == node; }),
                                          parent_children.end());
                }
                continue;
            }
            
            if (children.size() == 1 && node->parent != nullptr && node->blockId == 0) {
                if (node->parent->grpID == grpID) {
                    for (auto& child_ptr : node->parent->children) {
                        if (child_ptr.get() == node) {
                            children[0]->branchLength += node->branchLength;
                            children[0]->parent = node->parent;
                            child_ptr = std::move(node->children[0]); // Transfer ownership
                            break;
                        }
                    }
                    continue;
                }
            }

            if (children.size() > 1) {
                while (children.size() > 1) {
                    std::vector<Node*> nodeLeft;
                    for (size_t i = 0; i < children.size() - 1; i += 2) {
                        int firstIdx = (nodeAlnOrder.count(children[i]->identifier)) ? nodeAlnOrder[children[i]->identifier] + 1 : 0;
                        int secondIdx = (nodeAlnOrder.count(children[i+1]->identifier)) ? nodeAlnOrder[children[i+1]->identifier] + 1 : 0;
                        int maxIdx = std::max(firstIdx, secondIdx);
                        nodeAlnOrder[children[i]->identifier] = maxIdx;
                        nodeAlnOrder[children[i+1]->identifier] = maxIdx;
                        alnOrder.push_back({{children[i], children[i+1]}, maxIdx});
                        nodeLeft.push_back(children[i]);
                    }
                    if (children.size() % 2 == 1) {
                        nodeLeft.push_back(children.back());
                    }
                    children = nodeLeft;
                }
            }

            if (children.size() == 1 && node->blockId != 0) {
                 int firstIdx = (nodeAlnOrder.count(node->identifier)) ? nodeAlnOrder[node->identifier] + 1 : 0;
                 int secondIdx = (nodeAlnOrder.count(children[0]->identifier)) ? nodeAlnOrder[children[0]->identifier] + 1 : 0;
                 int maxIdx = std::max(firstIdx, secondIdx);
                 nodeAlnOrder[node->identifier] = maxIdx;
                 nodeAlnOrder[children[0]->identifier] = maxIdx;
                 alnOrder.push_back({{node, children[0]}, maxIdx});
            }

            if (!children.empty()) {
                nodeAlnOrder[node->identifier] = nodeAlnOrder[children[0]->identifier];
            }
        }
    } else if (mode == 1) {
        while(!postStack.empty()) {
            Node* node = postStack.top();
            postStack.pop();
            if (node->parent != nullptr) {
                int firstIdx  = (nodeAlnOrder.count(node->identifier)) ? nodeAlnOrder[node->identifier]+1 : 0;
                int secondIdx = (nodeAlnOrder.count(node->parent->identifier)) ? nodeAlnOrder[node->parent->identifier]+1 : 0;
                int maxIdx = std::max(firstIdx, secondIdx);
                nodeAlnOrder[node->identifier] = maxIdx;
                nodeAlnOrder[node->parent->identifier] = maxIdx;
                alnOrder.push_back({ {node->parent, node}, maxIdx });
            }
        }
    } else {
        while(!postStack.empty()) {
            Node* node = postStack.top();
            postStack.pop();
            if (node->parent != nullptr) {
                alnOrder.push_back({ {node->parent, node}, 0 });
            }
        }
    }
}


void scheduling(Node* root, std::vector<NodePairVec>& levels, int mode) {
    levels.clear();
    std::stack<Node*> msaStack;
    if (root) {
        root->collectPostOrder(msaStack);
    }
    std::vector<std::pair<NodePair, int>> progressivePairs;
    int grpID = root ? root->grpID : -1;
    getProgressivePairs(progressivePairs, msaStack, grpID, mode);
    for (const auto& h : progressivePairs) {
        if (levels.size() < (size_t)h.second + 1) {
            levels.resize(h.second + 1);
        }
        levels[h.second].push_back(h.first);
    }
}

void updateNode(NodePairVec& nodes, BlockManager& blockManager) {
    for (auto& n : nodes) {
        for (int i = 0; i < 2; ++i) {
            Node* currentNode = (i == 0) ? n.first : n.second;
            if (currentNode->is_leaf() && currentNode->blockId == 0) {
                // Assuming the node identifier is the block ID (or can be mapped to it)
                try {
                    uint64_t blockId_val = std::stoull(currentNode->identifier);
                    auto block = blockManager.getBlock(blockId_val);
                    if (block) {
                        currentNode->blockId = blockId_val;
                    }
                } catch (const std::invalid_argument& e) {
                    // Handle cases where identifier is not a numeric ID
                }
            } else if (currentNode->blockId == 0) {
                int grpID = currentNode->grpID;
                 for (const auto& child_ptr : currentNode->children) {
                    Node* child = child_ptr.get();
                     if ((child->grpID == -1 || child->grpID == grpID) && (child != ((i == 0) ? n.second : n.first))) {
                        currentNode->blockId = child->blockId;
                        break;
                    }
                }
            }
        }
    }
}


void progressiveAlignment(Tree& T, Option& option, std::vector<NodePairVec>& alnPairsPerLevel, BlockManager& blockManager) {
    int level = 0;
    if (option.verbose) {
        std::cerr << "Total " << alnPairsPerLevel.size() << " levels.\n";
    }
    for (auto& m : alnPairsPerLevel) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        
        updateNode(m, blockManager);

        alignmentKernel(m, blockManager, option);
        
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> alnTime = alnEnd - alnStart;
        
        if (option.verbose) {
            std::cerr << "Level " << ++level << ", aligned " << m.size() 
                      << " pair(s) in " << alnTime.count() << " ms\n";
        }
    }
}


void msaOnSubtree(Tree& T, Option& option, BlockManager& blockManager, int subtree) {

    auto progressiveStart = std::chrono::high_resolution_clock::now();
    std::cerr << "============================\n";

    std::vector<NodePairVec> alnPairsPerLevel;
    // Mode based on alignment strategy
    int mode = 0; // Simplified logic, assuming 2 is PLACE_WO_TREE
    scheduling(T.root.get(), alnPairsPerLevel, mode);
    
    auto scheduleEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> scheduleTime = scheduleEnd - progressiveStart;
    if (option.verbose) {
        std::cerr << "Scheduling in " << scheduleTime.count() << " us\n";
    }

    progressiveAlignment(T, option, alnPairsPerLevel, blockManager);
    
    // Push MSA results to the root of the tree
    if (!alnPairsPerLevel.empty() && !alnPairsPerLevel.back().empty()) {
        Node* lastAligned = alnPairsPerLevel.back()[0].first;
        T.root->blockId = lastAligned->blockId;
        lastAligned->blockId = 0;
        // lastAligned->msaFreq.clear();
    }
    
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> progressiveTime = progressiveEnd - progressiveStart;

    std::cerr << "Alignment completed in " << progressiveTime.count() << " s\n";
}

} // namespace progressive
} // namespace mga
