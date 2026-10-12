#ifndef PYYLO_HPP
#include "phylogeny.hpp"
#endif

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace {
    void assignGroup(phylogeny::Node* rootNode, int grp_id, std::unordered_map<phylogeny::Node*, size_t>& unassigned_map) {
        std::vector<phylogeny::Node*> stack;
        if (rootNode) stack.push_back(rootNode);
        while (!stack.empty()) {
            phylogeny::Node* node = stack.back();
            stack.pop_back();
            if (unassigned_map[node] == 0) continue;
            node->grpID = grp_id;
            for (auto ch : node->children) {
                stack.push_back(ch);
            }
            unassigned_map[node] = 0;
        }
    }
}

void phylogeny::PartitionInfo::partitionTree(Node* root) {
    if (!root) return;

    std::unordered_map<Node*, size_t> unassigned_map;
    int next_grp_id = 1;

    // Iterative post-order traversal to avoid call stack overflow on deep trees
    std::vector<phylogeny::Node*> postOrder;
    std::vector<phylogeny::Node*> stack;
    stack.push_back(root);
    while (!stack.empty()) {
        phylogeny::Node* curr = stack.back();
        stack.pop_back();
        postOrder.push_back(curr);
        for (auto child : curr->children) {
            stack.push_back(child);
        }
    }
    std::reverse(postOrder.begin(), postOrder.end());

    for (phylogeny::Node* node : postOrder) {
        if (node->is_leaf()) {
            unassigned_map[node] = 1;
            continue;
        }

        size_t total_leaves = 0;
        std::vector<phylogeny::Node*> children_with_leaves;

        for (auto child : node->children) {
            size_t leaves = unassigned_map[child];
            if (leaves > 0) {
                children_with_leaves.push_back(child);
                total_leaves += leaves;
            }
        }

        if (total_leaves > this->maxPartitionSize) {
            // Cut from the child with the most unassigned leaves to minimize subtree count
            std::sort(children_with_leaves.begin(), children_with_leaves.end(), 
                      [&](phylogeny::Node* a, phylogeny::Node* b) {
                          return unassigned_map[a] > unassigned_map[b];
                      });

            for (auto child : children_with_leaves) {
                if (total_leaves <= this->maxPartitionSize) break;

                size_t cut_leaves = unassigned_map[child];

                // Find LCA (deepest node that covers all unassigned leaves in this branch)
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
                        unassigned_map[lca] = 0; // Clear bypassed intermediate node
                        lca = next_lca;
                    }
                    else break;
                }

                int grp_id = next_grp_id++;
                assignGroup(lca, grp_id, unassigned_map);
                this->partitionsRoot[lca->identifier] = std::make_pair(lca, cut_leaves);

                total_leaves -= cut_leaves;
            }
        }

        unassigned_map[node] = total_leaves;
    }

    size_t remaining_leaves = unassigned_map[root];
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