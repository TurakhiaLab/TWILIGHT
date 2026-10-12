#include "msa.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <stack>

namespace msa {
namespace accurate {

static void collectSubtreeRootsDFS(Node* node,
                                   const std::unordered_map<std::string, std::pair<Node*, size_t>>& partitionsRoot,
                                   std::vector<std::string>& orderedRootIDs,
                                   std::unordered_set<std::string>& visited) {
    if (!node) return;
    if (partitionsRoot.find(node->identifier) != partitionsRoot.end()) {
        if (visited.insert(node->identifier).second) {
            orderedRootIDs.push_back(node->identifier);
        }
    }
    for (auto child : node->children) {
        collectSubtreeRootsDFS(child, partitionsRoot, orderedRootIDs, visited);
    }
}

HierarchyPlan computeHierarchyPlan(Tree* T, PartitionInfo* P, Tree* subRoot_T, Option* option) {
    HierarchyPlan plan;
    plan.sampleRate = option->accSampleRate;
    plan.maxGroupReps = option->accMaxGroup;
    plan.totalSequences = 0;
    plan.totalLayers = 0;

    bool freeLocalP = false;
    PartitionInfo* effectiveP = P;
    Tree* effectiveSubRoot = subRoot_T;

    if (!effectiveP || effectiveP->partitionsRoot.empty()) {
        if (T && T->root && T->root->numLeaves > (size_t)option->accMaxGroup) {
            effectiveP = new PartitionInfo(option->accMaxGroup, 0, 0);
            effectiveP->partitionTree(T->root);
            effectiveSubRoot = constructTreeFromPartitions_new(T->root, effectiveP);
            freeLocalP = true;
        } else if (T && T->root) {
            HierarchyGroup g;
            g.layer = 1;
            g.groupID = 0;
            g.rootIdentifier = T->root->identifier;
            g.subtreeIDs = { 0 };
            g.totalSequences = T->root->numLeaves;
            g.repCount = T->root->numLeaves;
            g.nextRepCount = std::max((size_t)1, (size_t)std::ceil(g.repCount * plan.sampleRate));
            plan.layers.push_back({ g });
            plan.totalSequences = g.totalSequences;
            plan.totalLayers = 1;
            return plan;
        } else {
            return plan;
        }
    }

    // Step 1: Collect partition roots in tree topological proximity order
    std::vector<std::string> orderedRootIDs;
    std::unordered_set<std::string> visited;

    if (effectiveSubRoot && effectiveSubRoot->root) {
        collectSubtreeRootsDFS(effectiveSubRoot->root, effectiveP->partitionsRoot, orderedRootIDs, visited);
    }

    // Safety fallback: append any partition roots not visited in DFS
    for (const auto& kv : effectiveP->partitionsRoot) {
        if (visited.insert(kv.first).second) {
            orderedRootIDs.push_back(kv.first);
        }
    }

    // Step 2: Construct Layer 1 (Base Subtree Groups)
    std::vector<HierarchyGroup> layer1;
    layer1.reserve(orderedRootIDs.size());

    for (size_t i = 0; i < orderedRootIDs.size(); ++i) {
        const std::string& rootID = orderedRootIDs[i];
        size_t leafCount = effectiveP->partitionsRoot.at(rootID).second;

        int subtreeID = -1;
        if (effectiveP->partitionsRoot.size() > 1 && T && T->allNodes.find(rootID) != T->allNodes.end()) {
            subtreeID = T->allNodes[rootID]->grpID;
        } else {
            subtreeID = (int)i;
        }

        HierarchyGroup g;
        g.layer = 1;
        g.groupID = (int)i;
        g.rootIdentifier = rootID;
        g.subtreeIDs = { subtreeID };
        g.childGroupIDs = {};
        g.totalSequences = leafCount;
        g.repCount = leafCount; // In Layer 1, all sequences undergo intra-subtree alignment
        g.nextRepCount = std::max((size_t)1, (size_t)std::ceil(leafCount * plan.sampleRate));

        layer1.push_back(std::move(g));
        plan.totalSequences += leafCount;
    }

    plan.layers.push_back(std::move(layer1));

    // Step 3: Hierarchically group adjacent subtrees/groups until root is reached
    int currentLayerIndex = 1; // 1-based index (Layer 1 is already in plan.layers[0])

    // Helper lambdas for clade-based hierarchical grouping on effectiveSubRoot
    auto findCladeGroups = [](auto& self, phylogeny::Node* node, size_t maxGroupReps,
                              const std::unordered_map<phylogeny::Node*, size_t>& cladeReps,
                              std::vector<phylogeny::Node*>& groupRoots) -> void {
        auto it = cladeReps.find(node);
        if (it == cladeReps.end() || it->second == 0) return;
        if (it->second <= maxGroupReps) {
            groupRoots.push_back(node);
            return;
        }
        bool anyChildValid = false;
        for (auto ch : node->children) {
            auto itCh = cladeReps.find(ch);
            if (itCh != cladeReps.end() && itCh->second > 0) {
                anyChildValid = true;
                self(self, ch, maxGroupReps, cladeReps, groupRoots);
            }
        }
        if (!anyChildValid) {
            groupRoots.push_back(node);
        }
    };

    auto collectCladeTips = [](auto& self, phylogeny::Node* node,
                               const std::unordered_map<std::string, int>& rootIDToGroupIdx,
                               std::vector<int>& childGroupIDs) -> void {
        auto it = rootIDToGroupIdx.find(node->identifier);
        if (it != rootIDToGroupIdx.end()) {
            childGroupIDs.push_back(it->second);
            return;
        }
        for (auto ch : node->children) {
            self(self, ch, rootIDToGroupIdx, childGroupIDs);
        }
    };

    while (plan.layers.back().size() > 1) {
        const auto& prevLayer = plan.layers.back();
        std::unordered_map<std::string, int> rootIDToGroupIdx;
        for (size_t i = 0; i < prevLayer.size(); ++i) {
            rootIDToGroupIdx[prevLayer[i].rootIdentifier] = static_cast<int>(i);
        }

        // Post-order calculation of representative count under each node in effectiveSubRoot
        std::unordered_map<phylogeny::Node*, size_t> cladeReps;
        std::stack<phylogeny::Node*> postStack;
        if (effectiveSubRoot && effectiveSubRoot->root) {
            effectiveSubRoot->root->collectPostOrder(postStack);
        }

        while (!postStack.empty()) {
            phylogeny::Node* n = postStack.top();
            postStack.pop();
            size_t total = 0;
            auto it = rootIDToGroupIdx.find(n->identifier);
            if (it != rootIDToGroupIdx.end()) {
                total = prevLayer[it->second].nextRepCount;
            } else {
                for (auto ch : n->children) {
                    auto itCh = cladeReps.find(ch);
                    if (itCh != cladeReps.end()) {
                        total += itCh->second;
                    }
                }
            }
            cladeReps[n] = total;
        }

        size_t rootReps = (effectiveSubRoot && effectiveSubRoot->root && cladeReps.find(effectiveSubRoot->root) != cladeReps.end())
                        ? cladeReps[effectiveSubRoot->root] : 0;

        if (rootReps <= static_cast<size_t>(plan.maxGroupReps)) {
            // Entire remaining tree fits in a single final merge group
            HierarchyGroup curGroup;
            curGroup.layer = currentLayerIndex + 1;
            curGroup.groupID = 0;
            curGroup.rootIdentifier = effectiveSubRoot->root->identifier;
            curGroup.totalSequences = 0;
            curGroup.repCount = 0;

            for (const auto& child : prevLayer) {
                curGroup.childGroupIDs.push_back(child.groupID);
                curGroup.subtreeIDs.insert(curGroup.subtreeIDs.end(), child.subtreeIDs.begin(), child.subtreeIDs.end());
                curGroup.totalSequences += child.totalSequences;
                curGroup.repCount += child.nextRepCount;
            }
            curGroup.nextRepCount = std::max((size_t)1, (size_t)std::ceil(curGroup.repCount * plan.sampleRate));
            plan.layers.push_back({ std::move(curGroup) });
            break;
        } else {
            // Split into maximal clades
            std::vector<phylogeny::Node*> groupRoots;
            findCladeGroups(findCladeGroups, effectiveSubRoot->root, plan.maxGroupReps, cladeReps, groupRoots);

            // Fallback safety if findCladeGroups returned 1 or 0 roots
            if (groupRoots.size() <= 1) {
                groupRoots.clear();
                for (auto ch : effectiveSubRoot->root->children) {
                    if (cladeReps.find(ch) != cladeReps.end() && cladeReps[ch] > 0) {
                        groupRoots.push_back(ch);
                    }
                }
            }

            std::vector<HierarchyGroup> nextLayer;
            for (size_t gIdx = 0; gIdx < groupRoots.size(); ++gIdx) {
                phylogeny::Node* gRoot = groupRoots[gIdx];
                HierarchyGroup curGroup;
                curGroup.layer = currentLayerIndex + 1;
                curGroup.groupID = static_cast<int>(gIdx);
                curGroup.rootIdentifier = gRoot->identifier;
                curGroup.totalSequences = 0;
                curGroup.repCount = 0;

                std::vector<int> childIDs;
                collectCladeTips(collectCladeTips, gRoot, rootIDToGroupIdx, childIDs);
                for (int cID : childIDs) {
                    const auto& child = prevLayer[cID];
                    curGroup.childGroupIDs.push_back(child.groupID);
                    curGroup.subtreeIDs.insert(curGroup.subtreeIDs.end(), child.subtreeIDs.begin(), child.subtreeIDs.end());
                    curGroup.totalSequences += child.totalSequences;
                    curGroup.repCount += child.nextRepCount;
                }
                curGroup.nextRepCount = std::max((size_t)1, (size_t)std::ceil(curGroup.repCount * plan.sampleRate));
                nextLayer.push_back(std::move(curGroup));
            }
            plan.layers.push_back(std::move(nextLayer));
            currentLayerIndex++;
        }
    }

    plan.totalLayers = (int)plan.layers.size();
    if (freeLocalP) {
        if (effectiveSubRoot) delete effectiveSubRoot;
        if (effectiveP) delete effectiveP;
    }
    return plan;
}

void HierarchyPlan::printSummary(bool verbose) const {
    std::cerr << "\n================================================================================\n";
    std::cerr << "[Accurate Mode] Hierarchical Consistency Plan \n";
    std::cerr << "================================================================================\n";
    std::cerr << "Total Sequences         : " << totalSequences << "\n";
    std::cerr << "Max Reps Per Group      : " << maxGroupReps << "\n";
    std::cerr << "Sample Rate             : " << sampleRate << " (" << std::fixed << std::setprecision(1) << (sampleRate * 100.0f) << "%)\n";
    std::cerr << "Total Hierarchy Layers  : " << totalLayers << "\n";
    std::cerr << "--------------------------------------------------------------------------------\n";

    for (size_t l = 0; l < layers.size(); ++l) {
        int layerNum = (int)l + 1;
        const auto& layer = layers[l];
        std::cerr << "Layer " << layerNum << ": " << layer.size() << " consistency group" 
                  << (layer.size() > 1 ? "s" : "") << " (";
        
        if (layerNum == 1) {
            std::cerr << "Base subtrees, intra-subtree dense alignment <= " << maxGroupReps << " seqs/group)\n";
        } else if (l == layers.size() - 1) {
            std::cerr << "Final Root Merge Level, " << layer[0].repCount << " reps total)\n";
        } else {
            std::cerr << "Inter-Subtree Merge Level " << (layerNum - 1) << ", <= " << maxGroupReps << " reps/group)\n";
        }

        if (verbose || layer.size() <= 10) {
            for (const auto& g : layer) {
                std::cerr << "  Group " << std::setw(2) << g.groupID << ": covers " 
                          << std::setw(3) << g.subtreeIDs.size() << " subtree(s) [";
                if (g.subtreeIDs.size() <= 4) {
                    for (size_t s = 0; s < g.subtreeIDs.size(); ++s) {
                        std::cerr << (s > 0 ? ", " : "") << g.subtreeIDs[s];
                    }
                } else {
                    std::cerr << g.subtreeIDs.front() << " .. " << g.subtreeIDs.back();
                }
                std::cerr << "], total " << g.totalSequences << " seqs, " 
                          << g.repCount << " reps (next level: " << g.nextRepCount << " reps)\n";
            }
        } else {
            std::cerr << "  [Detailed breakdown of " << layer.size() << " groups omitted, use --verbose to show all]\n";
        }
    }
    std::cerr << "================================================================================\n\n";
}

} // namespace accurate
} // namespace msa
