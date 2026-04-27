#ifndef PYYLO_HPP
#include "phylogeny.hpp"
#endif

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace {
    void assignGroup(phylogeny::Node* node, int grp_id, std::unordered_map<phylogeny::Node*, size_t>& unassigned_map) {
        if (unassigned_map[node] == 0) return;
        node->grpID = grp_id;
        for (auto ch : node->children) {
            assignGroup(ch, grp_id, unassigned_map);
        }
        unassigned_map[node] = 0;
    }

    size_t greedyPartition(phylogeny::Node* node, int& next_grp_id, size_t max_size, 
                           std::unordered_map<std::string, std::pair<phylogeny::Node*, size_t>>& partitionsRoot, 
                           std::unordered_map<phylogeny::Node*, size_t>& unassigned_map) {
        if (node->is_leaf()) {
            unassigned_map[node] = 1;
            return 1;
        }
        
        size_t total_leaves = 0;
        std::vector<phylogeny::Node*> children_with_leaves;
        
        for (auto child : node->children) {
            size_t leaves = greedyPartition(child, next_grp_id, max_size, partitionsRoot, unassigned_map);
            if (leaves > 0) {
                children_with_leaves.push_back(child);
                total_leaves += leaves;
            }
        }
        
        if (total_leaves > max_size) {
            // 從包含最多 unassigned leaves 的 child 開始切，盡量減少 subtree 數量
            std::sort(children_with_leaves.begin(), children_with_leaves.end(), 
                      [&](phylogeny::Node* a, phylogeny::Node* b) {
                          return unassigned_map[a] > unassigned_map[b];
                      });
            
            for (auto child : children_with_leaves) {
                if (total_leaves <= max_size) break;
                
                size_t cut_leaves = unassigned_map[child];
                
                // 尋找 LCA (能包含所有 unassigned tips 的最下面 node)
                phylogeny::Node* lca = child;
                while (!lca->is_leaf()) {
                    phylogeny::Node* next_lca = nullptr;
                    for (auto c : lca->children) {
                        if (unassigned_map[c] == cut_leaves) {
                            next_lca = c;
                            break;
                        }
                    }
                    if (next_lca) {
                        unassigned_map[lca] = 0; // 將 bypassed node 歸零
                        lca = next_lca;
                    }
                    else break;
                }
                
                int grp_id = next_grp_id++;
                assignGroup(lca, grp_id, unassigned_map);
                partitionsRoot[lca->identifier] = std::make_pair(lca, cut_leaves);
                
                total_leaves -= cut_leaves;
            }
        }
        
        unassigned_map[node] = total_leaves;
        return total_leaves;
    }
}

void phylogeny::PartitionInfo::partitionTree(Node* root) {
    std::unordered_map<Node*, size_t> unassigned_map;
    int next_grp_id = 1;
    
    size_t remaining_leaves = greedyPartition(root, next_grp_id, this->maxPartitionSize, this->partitionsRoot, unassigned_map);
    
    if (remaining_leaves > 0) {
        Node* lca = root;
        while (!lca->is_leaf()) {
            Node* next_lca = nullptr;
            for (auto c : lca->children) {
                if (unassigned_map[c] == remaining_leaves) {
                    next_lca = c;
                    break;
                }
            }
            if (next_lca) {
                unassigned_map[lca] = 0;
                lca = next_lca;
            }
            else break;
        }
        assignGroup(lca, 0, unassigned_map);
        this->partitionsRoot[lca->identifier] = std::make_pair(lca, remaining_leaves);
    }
    
    this->numPartitions = this->partitionsRoot.size();
}

void phylogeny::PartitionInfo::bipartition(Node* root, Node* edge, Node*& tree1Root, Node*& tree2Root) {
}

void phylogeny::PartitionInfo::bipartition(phylogeny::Node* rootOfPartition, phylogeny::Node* edgeToCut, int newGrpID) {
}

phylogeny::PartitionInfo::~PartitionInfo() {
    this->partitionsRoot.clear();
}