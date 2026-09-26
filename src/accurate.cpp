#include "msa.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <atomic>
#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

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
    std::cerr << "[Accurate Mode] Hierarchical Consistency Plan (Plan B)\n";
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

// Collect all leaves under a node in topological DFS order
static void collectLeavesInOrder(Node* node, std::vector<std::string>& leaves) {
    if (!node) return;
    if (node->is_leaf() || node->children.empty()) {
        leaves.push_back(node->identifier);
        return;
    }
    for (auto* child : node->children) {
        collectLeavesInOrder(child, leaves);
    }
}

// Systematically sample targetReps representative sequences evenly across tree clades
std::vector<int> selectRepresentativeIndices(Node* root, size_t targetReps, SequenceDB* database) {
    std::vector<std::string> leaves;
    collectLeavesInOrder(root, leaves);

    std::vector<int> selectedIndices;
    if (leaves.empty() || targetReps == 0) return selectedIndices;

    if (leaves.size() <= targetReps) {
        for (const auto& name : leaves) {
            auto it = database->name_map.find(name);
            if (it != database->name_map.end()) {
                selectedIndices.push_back(it->second->id);
            }
        }
        return selectedIndices;
    }

    // Uniform stride sampling across phylogenetic clades to avoid clade clustering
    double stride = static_cast<double>(leaves.size()) / static_cast<double>(targetReps);
    for (size_t r = 0; r < targetReps; ++r) {
        size_t idx = static_cast<size_t>(r * stride);
        if (idx >= leaves.size()) idx = leaves.size() - 1;
        const auto& name = leaves[idx];
        auto it = database->name_map.find(name);
        if (it != database->name_map.end()) {
            selectedIndices.push_back(it->second->id);
        }
    }
    return selectedIndices;
}

// Extract representative sequences using centrality scores across phylogenetic clades
std::vector<Representative> extractRepresentativesByCentrality(
    Node* root,
    size_t targetReps,
    SequenceDB* database,
    int subtreeID,
    const std::vector<float>& centralityScores,
    const LocalHomTable& localHomTable,
    int alnLen
) {
    std::vector<std::string> leaves;
    collectLeavesInOrder(root, leaves);

    std::vector<Representative> reps;
    if (leaves.empty() || targetReps == 0) return reps;

    std::vector<std::string> selectedNames;
    if (leaves.size() <= targetReps) {
        selectedNames = leaves;
    } else {
        // Uniform stride partitioning across phylogenetic clades to ensure diversity
        double stride = static_cast<double>(leaves.size()) / static_cast<double>(targetReps);
        for (size_t r = 0; r < targetReps; ++r) {
            size_t startIdx = static_cast<size_t>(r * stride);
            size_t endIdx = static_cast<size_t>((r + 1) * stride);
            if (endIdx > leaves.size() || r == targetReps - 1) endIdx = leaves.size();

            // Select candidate with highest centrality score in [startIdx, endIdx)
            std::string bestName = leaves[startIdx];
            float bestScore = -1.0f;
            for (size_t idx = startIdx; idx < endIdx; ++idx) {
                const auto& name = leaves[idx];
                auto it = database->name_map.find(name);
                if (it != database->name_map.end()) {
                    int globalId = it->second->id;
                    int locIdx = localHomTable.getLocalIndex(globalId);
                    float score = (locIdx >= 0 && locIdx < static_cast<int>(centralityScores.size()))
                                  ? centralityScores[locIdx] : 0.0f;
                    if (score > bestScore) {
                        bestScore = score;
                        bestName = name;
                    }
                }
            }
            selectedNames.push_back(bestName);
        }
    }

    reps.reserve(selectedNames.size());
    for (size_t rIdx = 0; rIdx < selectedNames.size(); ++rIdx) {
        const auto& name = selectedNames[rIdx];
        auto it = database->name_map.find(name);
        if (it != database->name_map.end()) {
            const auto* sInfo = it->second;
            Representative rep;
            rep.globalSeqID = sInfo->id;
            rep.subtreeID = subtreeID;
            rep.name = name;
            rep.originGroupID = -1; // To be set by caller
            rep.localIdxInChild = static_cast<int>(rIdx);

            int locIdx = localHomTable.getLocalIndex(sInfo->id);
            rep.centralityScore = (locIdx >= 0 && locIdx < static_cast<int>(centralityScores.size()))
                                  ? centralityScores[locIdx] : 0.0f;

            // Extract unaligned sequence string (strip '-' and '.') and record column mapping
            const char* rawChars = sInfo->alnStorage[sInfo->storage];
            rep.unalignedSeq.reserve(sInfo->len);
            rep.colInSubAln.clear();
            for (int k = 0; k < sInfo->len && k < alnLen; ++k) {
                if (rawChars[k] != '-' && rawChars[k] != '.') {
                    rep.unalignedSeq.push_back(rawChars[k]);
                    rep.colInSubAln.push_back(k);
                }
            }
            reps.push_back(std::move(rep));
        }
    }
    return reps;
}

// Extract full representative sequence objects (name, unaligned string, subtree ID) for inter-layer consistency
std::vector<Representative> extractRepresentatives(Node* root, size_t targetReps, SequenceDB* database, int subtreeID) {
    std::vector<std::string> leaves;
    collectLeavesInOrder(root, leaves);

    std::vector<Representative> reps;
    if (leaves.empty() || targetReps == 0) return reps;

    std::vector<std::string> selectedNames;
    if (leaves.size() <= targetReps) {
        selectedNames = leaves;
    } else {
        double stride = static_cast<double>(leaves.size()) / static_cast<double>(targetReps);
        for (size_t r = 0; r < targetReps; ++r) {
            size_t idx = static_cast<size_t>(r * stride);
            if (idx >= leaves.size()) idx = leaves.size() - 1;
            selectedNames.push_back(leaves[idx]);
        }
    }

    reps.reserve(selectedNames.size());
    for (const auto& name : selectedNames) {
        auto it = database->name_map.find(name);
        if (it != database->name_map.end()) {
            const auto* sInfo = it->second;
            Representative rep;
            rep.globalSeqID = sInfo->id;
            rep.subtreeID = subtreeID;
            rep.name = name;

            // Extract unaligned sequence string (strip '-' and '.')
            const char* rawChars = sInfo->alnStorage[sInfo->storage];
            rep.unalignedSeq.reserve(sInfo->len);
            for (int k = 0; k < sInfo->len; ++k) {
                if (rawChars[k] != '-' && rawChars[k] != '.') {
                    rep.unalignedSeq.push_back(rawChars[k]);
                }
            }
            reps.push_back(std::move(rep));
        }
    }
    return reps;
}

void msaOnSubtree_accurate(Tree *T, SequenceDB *database, Option *option, Params &param, alnFunction alignmentKernel, int subtree) {
    phylogeny::PartitionInfo* P = nullptr;
    phylogeny::Tree* subRoot_T = nullptr;
    T->convert2binaryTree();
    T->calLeafNum();
    int groupSize = std::min(option->maxSubtree, option->accMaxGroup);
    if (T->root && T->root->numLeaves > (size_t)groupSize) {
        P = new phylogeny::PartitionInfo(groupSize, 0, 0); 
        P->partitionTree(T->root);
        subRoot_T = phylogeny::constructTreeFromPartitions_new(T->root, P);
        if (P->partitionsRoot.size() > 1) {
            if (option->maxSubtree < INT32_MAX) {
                msa::io::writeSubtrees(T, P, option);
            }
        }
    }

    // Calculate and display hierarchical consistency alignment plan (Plan B)
    auto hierarchyPlan = computeHierarchyPlan(T, P, subRoot_T, option);
    hierarchyPlan.printSummary(option->printDetail);

    // Execute hierarchical consistency alignment plan
    executeHierarchyPlan(hierarchyPlan, T, P, subRoot_T, database, option, param, alignmentKernel);

    if (subRoot_T) delete subRoot_T;
    if (P) delete P;
}

// -----------------------------------------------------------------------------
// 1. All-to-all pairwise alignment for consistency score
// -----------------------------------------------------------------------------
std::shared_ptr<LocalHomTable> alignGroupPairwiseAllToAll(
    HierarchyGroup& group,
    Tree* T,
    SequenceDB* database,
    Option* option,
    Params& param,
    const std::vector<HierarchyGroup>* prevLayerGroups)
{
    auto alignStart = std::chrono::high_resolution_clock::now();

    // 1. Gather active sequence strings and indices for this group
    std::vector<std::string> seqStrings;
    std::vector<int> activeSeqIndices;

    if (!group.representatives.empty()) {
        // Inter-subtree layer (Layer 2+): Align representative sequences
        seqStrings.reserve(group.representatives.size());
        activeSeqIndices.reserve(group.representatives.size());
        for (size_t i = 0; i < group.representatives.size(); ++i) {
            seqStrings.push_back(group.representatives[i].unalignedSeq);
            activeSeqIndices.push_back(group.representatives[i].globalSeqID >= 0 ? group.representatives[i].globalSeqID : (int)i);
        }
    } else {
        // Base layer (Layer 1): Intra-subtree alignment using sequences currently loaded in SequenceDB
        std::vector<std::string> leafNames;
        if (T && T->root) {
            collectLeavesInOrder(T->root, leafNames);
        }
        if (!leafNames.empty()) {
            activeSeqIndices.reserve(leafNames.size());
            for (const auto& name : leafNames) {
                auto it = database->name_map.find(name);
                if (it != database->name_map.end()) {
                    int s = it->second->id;
                    if (!database->sequences[s]->lowQuality || option->noFilter) {
                        activeSeqIndices.push_back(s);
                    }
                }
            }
        } else {
            activeSeqIndices.reserve(database->sequences.size());
            for (size_t s = 0; s < database->sequences.size(); ++s) {
                if (!database->sequences[s]->lowQuality || option->noFilter) {
                    activeSeqIndices.push_back((int)s);
                }
            }
        }
        seqStrings.resize(activeSeqIndices.size());
        for (size_t i = 0; i < activeSeqIndices.size(); ++i) {
            const auto* sInfo = database->sequences[activeSeqIndices[i]];
            const char* rawChars = sInfo->alnStorage[sInfo->storage];
            std::string s;
            s.reserve(sInfo->len);
            for (int k = 0; k < sInfo->len; ++k) {
                if (rawChars[k] != '-' && rawChars[k] != '.') {
                    s.push_back(rawChars[k]);
                }
            }
            seqStrings[i] = std::move(s);
        }
    }

    int N = static_cast<int>(seqStrings.size());
    auto localHomTable = std::make_shared<LocalHomTable>(N, activeSeqIndices);

    if (N < 2) {
        return localHomTable;
    }

    // 2. Prepare all-to-all pair jobs (upper triangle: i < j)
    // Check for cached intra-group pairwise results from child groups (Scheme B)
    size_t totalPairs = (size_t)N * (N - 1) / 2;
    std::vector<std::pair<int, int>> pairJobs;
    pairJobs.reserve(totalPairs);
    size_t reusedCount = 0;
    size_t reusedSegments = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            bool reused = false;
            if (!group.representatives.empty() && prevLayerGroups != nullptr) {
                int orig_i = group.representatives[i].originGroupID;
                int orig_j = group.representatives[j].originGroupID;
                if (orig_i >= 0 && orig_i == orig_j && orig_i < static_cast<int>(prevLayerGroups->size())) {
                    uint64_t key = makeGlobalPairKey(group.representatives[i].globalSeqID, group.representatives[j].globalSeqID);
                    const auto& cache = (*prevLayerGroups)[orig_i].repPairwiseCache;
                    auto it = cache.find(key);
                    if (it != cache.end()) {
                        // Directly reuse cached pairwise result from child group (0 SW computation)
                        localHomTable->setPairResult(i, j, LocalHomTable::PairResult(it->second));
                        reusedSegments += it->second.segments.size();
                        ++reusedCount;
                        reused = true;
                    }
                }
            }
            if (!reused) {
                pairJobs.push_back({i, j});
            }
        }
    }

    std::atomic<size_t> completedCount{0};
    std::atomic<size_t> totalSegments{reusedSegments};
    size_t numJobs = pairJobs.size();

    // 3. Parallel execution via TBB for remaining inter-group jobs
    if (numJobs > 0) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, numJobs),
            [&](const tbb::blocked_range<size_t>& range) {
                Aligner aligner; // Thread-local aligner
                size_t localSegments = 0;

                for (size_t p = range.begin(); p < range.end(); ++p) {
                    int i = pairJobs[p].first;
                    int j = pairJobs[p].second;

                    int len1 = static_cast<int>(seqStrings[i].size());
                    int len2 = static_cast<int>(seqStrings[j].size());
                    int minL = std::min(len1, len2);
                    int maxL = std::max(len1, len2);
                    int diff = maxL - minL;

                    LocalHomTable::PairResult pairRes;

                    // Criterion 3: Banded optimization for long sequences (minL >= 500) with similar lengths (diff <= 15% of minL)
                    if (minL >= 500 && (static_cast<float>(diff) / static_cast<float>(minL) <= 0.15f)) {
                        int bandWidth = diff + std::max(128, static_cast<int>(0.10f * minL));
                        if (bandWidth < minL) {
                            pairRes = aligner.align_affine_local_segments_banded(
                                seqStrings[i], seqStrings[j], option->type, param, bandWidth
                            );
                        } else {
                            pairRes = aligner.align_affine_local_segments(
                                seqStrings[i], seqStrings[j], option->type, param
                            );
                        }
                    } else {
                        pairRes = aligner.align_affine_local_segments(
                            seqStrings[i], seqStrings[j], option->type, param
                        );
                    }

                    localSegments += pairRes.segments.size();

                    // Assign result into upper-triangular matrix via move (zero-copy)
                    localHomTable->setPairResult(i, j, std::move(pairRes));

                    size_t done = ++completedCount;
                    if (option->printDetail && (done % 500 == 0 || done == numJobs)) {
                        double pct = 100.0 * done / numJobs;
                        std::cerr << "\r  [All-to-All Pairwise SW] " << done << "/" << numJobs 
                                  << " (" << std::fixed << std::setprecision(1) << pct << "%) completed..." << std::flush;
                    }
                }
                totalSegments += localSegments;
            });

        if (option->printDetail && numJobs > 0) {
            std::cerr << "\n";
        }
    }

    auto alignEnd = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(alignEnd - alignStart);

    if (option->printDetail) {
        std::cerr << "  Group " << group.groupID << ": Completed " << totalPairs 
                  << " pairs (" << reusedCount << " reused from child groups, " 
                  << numJobs << " computed with SW, total " << totalSegments.load() << " segments) in " 
                  << elapsed.count() << " ms.\n";
    }

    return localHomTable;
}

// -----------------------------------------------------------------------------
// 2. Build consistency score library & calculate centrality scores
// -----------------------------------------------------------------------------
std::vector<float> buildGroupConsistencyLibrary(const HierarchyGroup& group, std::shared_ptr<LocalHomTable> localHomTable, Tree* T, SequenceDB* database, Option* option, Params& param) {
    auto libStart = std::chrono::high_resolution_clock::now();
    int N = localHomTable->getNumSeqs();
    if (N < 2) return std::vector<float>(N, 0.0f);

    if (option->printDetail) {
        std::cerr << "  [Step 2/3] Building consistency library for Group " << group.groupID 
                  << " (N = " << N << " sequences)...\n";
    }

    // 1. Collect sequence lengths and sequence weights
    std::vector<int> seqLengths(N, 0);
    std::vector<float> seqWeights(N, 1.0f);

    if (!group.representatives.empty()) {
        for (int i = 0; i < N; ++i) {
            seqLengths[i] = static_cast<int>(group.representatives[i].unalignedSeq.size());
            seqWeights[i] = 1.0f; // Equal weight for reps
        }
    } else {
        for (int i = 0; i < N; ++i) {
            int gIdx = localHomTable->getGlobalIndex(i);
            if (gIdx >= 0 && gIdx < static_cast<int>(database->sequences.size())) {
                const auto* sInfo = database->sequences[gIdx];
                const char* rawChars = sInfo->alnStorage[sInfo->storage];
                int unalignedLen = 0;
                for (int k = 0; k < sInfo->len; ++k) {
                    if (rawChars[k] != '-' && rawChars[k] != '.') ++unalignedLen;
                }
                seqLengths[i] = unalignedLen;
                seqWeights[i] = (sInfo->weight > 0.0f) ? sInfo->weight : 1.0f;
            }
        }
    }

    // Normalize weights so sum(weights) == 1.0
    float totalWeight = 0.0f;
    for (float w : seqWeights) totalWeight += w;
    if (totalWeight > 0.0f) {
        for (float& w : seqWeights) w /= totalWeight;
    } else {
        for (float& w : seqWeights) w = 1.0f / N;
    }

    // 2. Allocate position-specific importance table (Scratchpad: N arrays of length L_i)
    // Memory footprint: sum(L_i) * sizeof(float) <= N * L_max * 4 bytes (~a few MBs)
    std::vector<std::vector<float>> positionImportance(N);
    size_t totalResidues = 0;
    for (int i = 0; i < N; ++i) {
        positionImportance[i].assign(seqLengths[i], 0.0f);
        totalResidues += seqLengths[i];
    }

    // 3. Parallel position-specific coverage weight accumulation (Lock-free per sequence i)
    tbb::parallel_for(0, N, [&](int i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            const auto& pairRes = localHomTable->getPairResult(i, j);
            float w_j = seqWeights[j];
            if (i < j) {
                for (const auto& seg : pairRes.segments) {
                    int s1 = std::max(0, seg.start1);
                    int e1 = std::min(seqLengths[i] - 1, seg.end1);
                    for (int pos = s1; pos <= e1; ++pos) {
                        positionImportance[i][pos] += w_j;
                    }
                }
            } else {
                for (const auto& seg : pairRes.segments) {
                    int s2 = std::max(0, seg.start2);
                    int e2 = std::min(seqLengths[i] - 1, seg.end2);
                    for (int pos = s2; pos <= e2; ++pos) {
                        positionImportance[i][pos] += w_j;
                    }
                }
            }
        }
    });

    // 4. Compute weighted importance for each pairwise local alignment segment
    size_t totalPairs = (size_t)N * (N - 1) / 2;
    std::vector<std::pair<int, int>> pairJobs;
    pairJobs.reserve(totalPairs);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            pairJobs.push_back({i, j});
        }
    }

    std::atomic<size_t> totalSegmentsWeighted{0};
    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPairs),
        [&](const tbb::blocked_range<size_t>& range) {
            size_t localCount = 0;
            for (size_t p = range.begin(); p < range.end(); ++p) {
                int i = pairJobs[p].first;
                int j = pairJobs[p].second;
                auto& pairRes = localHomTable->getPairResult(i, j);
                for (auto& seg : pairRes.segments) {
                    // Evaluate Sequence 1 positional coverage across [start1, end1]
                    float sum_imp1 = 0.0f;
                    int s1 = std::max(0, seg.start1);
                    int e1 = std::min(seqLengths[i] - 1, seg.end1);
                    int len1 = e1 - s1 + 1;
                    for (int pos = s1; pos <= e1; ++pos) {
                        sum_imp1 += positionImportance[i][pos];
                    }
                    float mean_imp1 = (len1 > 0) ? (sum_imp1 / static_cast<float>(len1)) : 0.0f;

                    // Evaluate Sequence 2 positional coverage across [start2, end2]
                    float sum_imp2 = 0.0f;
                    int s2 = std::max(0, seg.start2);
                    int e2 = std::min(seqLengths[j] - 1, seg.end2);
                    int len2 = e2 - s2 + 1;
                    for (int pos = s2; pos <= e2; ++pos) {
                        sum_imp2 += positionImportance[j][pos];
                    }
                    float mean_imp2 = (len2 > 0) ? (sum_imp2 / static_cast<float>(len2)) : 0.0f;

                    // Forward and reverse directional importance
                    seg.importance = mean_imp1 * seg.optScore;
                    seg.rimportance = mean_imp2 * seg.optScore;

                    // MAFFT symmetric importance: average both sequences' positional reliabilities
                    float sym_imp = 0.5f * (seg.importance + seg.rimportance);
                    seg.importance = sym_imp;
                    seg.rimportance = sym_imp;
                }
                localCount += pairRes.segments.size();
            }
            totalSegmentsWeighted += localCount;
        });

    // 5. Compute summary statistics for verification
    float minImp = 1e9f, maxImp = -1e9f;
    double sumImp = 0.0;
    size_t sampledPositions = 0;
    for (int i = 0; i < N; ++i) {
        for (float val : positionImportance[i]) {
            if (val < minImp) minImp = val;
            if (val > maxImp) maxImp = val;
            sumImp += val;
            ++sampledPositions;
        }
    }
    float avgImp = (sampledPositions > 0) ? static_cast<float>(sumImp / sampledPositions) : 0.0f;

    auto libEnd = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(libEnd - libStart);

    if (option->printDetail) {
        std::cerr << "  Consistency Library built in " << elapsed.count() << " ms:\n";
        std::cerr << "    Total Residues Evaluated: " << totalResidues 
                  << " (Table Memory: " << std::fixed << std::setprecision(2) 
                  << (totalResidues * sizeof(float) / (1024.0 * 1024.0)) << " MB)\n";
        std::cerr << "    Position Importance Range: [" << minImp << ", " << maxImp << "], Mean: " << avgImp << "\n";
        std::cerr << "    Total Weighted Segments : " << totalSegmentsWeighted.load() 
                  << " (Symmetric Mean Importance Applied)\n";
    }

    // 6. Compute sequence centrality scores: mean(positionImportance[i]) * min(1.0, len_i / medianLen)
    std::vector<int> sortedLens = seqLengths;
    std::sort(sortedLens.begin(), sortedLens.end());
    float medianLen = (N > 0) ? static_cast<float>(sortedLens[N / 2]) : 1.0f;
    if (medianLen <= 0.0f) medianLen = 1.0f;

    std::vector<float> centralityScores(N, 0.0f);
    for (int i = 0; i < N; ++i) {
        if (seqLengths[i] == 0) continue;
        double sumImpI = 0.0;
        for (float v : positionImportance[i]) {
            sumImpI += v;
        }
        float meanImpI = static_cast<float>(sumImpI / seqLengths[i]);
        float lengthFactor = std::min(1.0f, static_cast<float>(seqLengths[i]) / medianLen);
        centralityScores[i] = meanImpI * lengthFactor;
    }

    // 7. Store positional importance inside LocalHomTable for downstream profile DP lookups
    localHomTable->setPositionImportance(std::move(positionImportance));

    return centralityScores;
}

// -----------------------------------------------------------------------------
// Consistency Table Builder: Sparse Segment Scatter-Add (MAFFT fillimp Phase 3)
// -----------------------------------------------------------------------------
std::vector<std::vector<float>> buildConsistencyTableFromLocalHom(
    const ColumnProvenance& refProvenance,
    const ColumnProvenance& qryProvenance,
    const LocalHomTable& localHomTable)
{
    int refLen = static_cast<int>(refProvenance.size());
    int qryLen = static_cast<int>(qryProvenance.size());
    if (refLen == 0 || qryLen == 0) return {};

    std::vector<std::vector<float>> table(refLen, std::vector<float>(qryLen, 0.0f));

    // Fast mapping from (localSeqIdx, residueIndex) to profile column
    std::unordered_map<int, std::vector<int>> refPosToCol;
    for (int col = 0; col < refLen; ++col) {
        for (const auto& r : refProvenance[col]) {
            int localA = localHomTable.getLocalIndex(r.seqId);
            if (localA < 0) continue;
            auto& vec = refPosToCol[localA];
            if (static_cast<int>(vec.size()) <= r.residueIndex) {
                vec.resize(r.residueIndex + 1, -1);
            }
            vec[r.residueIndex] = col;
        }
    }

    std::unordered_map<int, std::vector<int>> qryPosToCol;
    for (int col = 0; col < qryLen; ++col) {
        for (const auto& q : qryProvenance[col]) {
            int localB = localHomTable.getLocalIndex(q.seqId);
            if (localB < 0) continue;
            auto& vec = qryPosToCol[localB];
            if (static_cast<int>(vec.size()) <= q.residueIndex) {
                vec.resize(q.residueIndex + 1, -1);
            }
            vec[q.residueIndex] = col;
        }
    }

    if (refPosToCol.empty() || qryPosToCol.empty()) return table;

    // Sparse Segment Scatter-Add (MAFFT fillimp Phase 3)
    for (const auto& refEntry : refPosToCol) {
        int localA = refEntry.first;
        const auto& refMap = refEntry.second;

        for (const auto& qryEntry : qryPosToCol) {
            int localB = qryEntry.first;
            if (localA == localB) continue;

            const auto& qryMap = qryEntry.second;
            const auto& pairRes = localHomTable.getPairResult(localA, localB);
            if (pairRes.segments.empty()) continue;

            if (localA < localB) {
                for (const auto& seg : pairRes.segments) {
                    float imp = seg.importance;
                    if (imp <= 0.0f) continue;
                    int segLen = seg.end1 - seg.start1 + 1;
                    for (int k = 0; k < segLen; ++k) {
                        int posA = seg.start1 + k;
                        int posB = seg.start2 + k;
                        if (posA < static_cast<int>(refMap.size()) && posB < static_cast<int>(qryMap.size())) {
                            int cA = refMap[posA];
                            int cB = qryMap[posB];
                            if (cA >= 0 && cB >= 0) {
                                table[cA][cB] += imp;
                            }
                        }
                    }
                }
            } else {
                for (const auto& seg : pairRes.segments) {
                    float imp = seg.importance;
                    if (imp <= 0.0f) continue;
                    int segLen = seg.end1 - seg.start1 + 1;
                    for (int k = 0; k < segLen; ++k) {
                        int posB = seg.start1 + k; // Sequence with smaller index is localB
                        int posA = seg.start2 + k; // Sequence with larger index is localA
                        if (posA < static_cast<int>(refMap.size()) && posB < static_cast<int>(qryMap.size())) {
                            int cA = refMap[posA];
                            int cB = qryMap[posB];
                            if (cA >= 0 && cB >= 0) {
                                table[cA][cB] += imp;
                            }
                        }
                    }
                }
            }
        }
    }

    // Normalize consistency bonus by the number of sequence pairs crossing the two sub-profiles
    float norm = 1.0f / (static_cast<float>(refPosToCol.size()) * static_cast<float>(qryPosToCol.size()));
    for (int i = 0; i < refLen; ++i) {
        for (int j = 0; j < qryLen; ++j) {
            table[i][j] *= norm;
        }
    }

    return table;
}

// -----------------------------------------------------------------------------
// Accurate Alignment Kernel (CPU): Profile DP with Consistency Guidance
// -----------------------------------------------------------------------------
void alignmentKernel_Accurate_CPU(
    Tree* tree,
    NodePairVec& nodes,
    SequenceDB* database,
    Option* option,
    Params& param,
    std::shared_ptr<LocalHomTable> localHomTable,
    std::unordered_map<std::string, ColumnProvenance>* subrootProvenance)
{
    int profileSize = param.matrixSize + 1;

    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range) {
        for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
            int32_t refLen = nodes[nIdx].first->getAlnLen(database->currentTask);
            int32_t qryLen = nodes[nIdx].second->getAlnLen(database->currentTask);
            int32_t refNum = nodes[nIdx].first->getAlnNum(database->currentTask);
            int32_t qryNum = nodes[nIdx].second->getAlnNum(database->currentTask);
            int32_t memLen = std::max(refLen, qryLen);

            float* hostFreq = nullptr;
            float* hostGapOp = nullptr;
            float* hostGapEx = nullptr;
            msa::progressive::cpu::allocateMemory_and_Initialize(hostFreq, hostGapOp, hostGapEx, memLen, profileSize);

            std::pair<IntPairVec, IntPairVec> gappyColumns;
            stringPair consensus({"", ""});
            IntPair lens = {refLen, qryLen};

            accurate::ColumnProvenance refProvenance, qryProvenance;
            std::vector<std::vector<float>> consistencyTable;

            alignment_helper::calculateProfile(hostFreq, nodes[nIdx], database, option, memLen);

            // Extract column provenance for reference and query sub-MSAs
            if (subrootProvenance != nullptr) {
                auto itRef = subrootProvenance->find(nodes[nIdx].first->identifier);
                if (itRef != subrootProvenance->end()) {
                    refProvenance = itRef->second;
                }
                auto itQry = subrootProvenance->find(nodes[nIdx].second->identifier);
                if (itQry != subrootProvenance->end()) {
                    qryProvenance = itQry->second;
                }
            } else if (database->currentTask == 0) {
                alignment_helper::extractColumnProvenance(nodes[nIdx].first, database, refProvenance);
                alignment_helper::extractColumnProvenance(nodes[nIdx].second, database, qryProvenance);
            }

            alignment_helper::getConsensus(option, hostFreq, consensus.first, refLen);
            alignment_helper::getConsensus(option, hostFreq + profileSize * memLen, consensus.second, qryLen);
            alignment_helper::removeGappyColumns(hostFreq, nodes[nIdx], option, gappyColumns, memLen, lens, database->currentTask);

            // Remove gappy columns to align coordinate systems with the DP matrix
            accurate::removeColumns(refProvenance, gappyColumns.first);
            accurate::removeColumns(qryProvenance, gappyColumns.second);

            // Build consistency lookup table from localHomTable
            if (localHomTable && option->accurate) {
                consistencyTable = buildConsistencyTableFromLocalHom(refProvenance, qryProvenance, *localHomTable);
            }

            // alignment_helper::calculatePSGP_MAFFT_new(hostGapOp, hostGapEx, nodes[nIdx], database, memLen, lens, param);
            msa::FAMSAProfile profRef, profQry;
            alignment_helper::calculatePSGP_FAMSA(profRef, profQry, nodes[nIdx], database, option, param, lens, gappyColumns);

            // Prepare profiles for DP alignment
            std::vector<int8_t> aln_wo_gc;
            Profile freqRef(lens.first, std::vector<float>(profileSize, 0.0f));
            Profile freqQry(lens.second, std::vector<float>(profileSize, 0.0f));
            Profile gapOp(2), gapEx(2);
            for (int s = 0; s < lens.first; ++s) {
                for (int t = 0; t < profileSize; ++t) {
                    freqRef[s][t] = hostFreq[profileSize * s + t];
                }
            }
            for (int s = 0; s < lens.second; ++s) {
                for (int t = 0; t < profileSize; ++t) {
                    freqQry[s][t] = hostFreq[profileSize * (memLen + s) + t];
                }
            }
            for (int r = 0; r < lens.first; ++r) {
                gapOp[0].push_back(hostGapOp[r]);
                gapEx[0].push_back(hostGapEx[r]);
            }
            for (int q = 0; q < lens.second; ++q) {
                gapOp[1].push_back(hostGapOp[memLen + q]);
                gapEx[1].push_back(hostGapEx[memLen + q]);
            }
            msa::progressive::cpu::freeMemory(hostFreq, hostGapOp, hostGapEx);

            std::pair<float, float> num = std::make_pair(static_cast<float>(refNum), static_cast<float>(qryNum));

            if (lens.first == 0) {
                for (int j = 0; j < lens.second; ++j) aln_wo_gc.push_back(1);
            } else if (lens.second == 0) {
                for (int j = 0; j < lens.first; ++j) aln_wo_gc.push_back(2);
            } else {
                int minL = std::min(lens.first, lens.second);
                int maxL = std::max(lens.first, lens.second);
                int diff = maxL - minL;

                // Criterion 3: Banded optimization for long profiles (minL >= 500) with similar lengths (diff <= 15% of minL)
                /*
                if (0 > 1) {
                // if (minL >= 500 && (static_cast<float>(diff) / static_cast<float>(minL) <= 0.15f)) {
                    int bandWidth = diff + std::max(128, static_cast<int>(0.10f * minL));
                    if (bandWidth < minL) {
                        aln_wo_gc = alignProfile_global_banded(
                            freqRef,
                            freqQry,
                            gapOp,
                            gapEx,
                            num,
                            param,
                            bandWidth,
                            (!consistencyTable.empty()) ? &consistencyTable : nullptr,
                            option->consistencyWeight
                        );
                    } else {
                        aln_wo_gc = alignProfile_global(
                            freqRef,
                            freqQry,
                            gapOp,
                            gapEx,
                            num,
                            param,
                            (!consistencyTable.empty()) ? &consistencyTable : nullptr,
                            option->consistencyWeight
                        );
                    }
                } else {
                    aln_wo_gc = alignProfile_global(
                        freqRef,
                        freqQry,
                        gapOp,
                        gapEx,
                        num,
                        param,
                        (!consistencyTable.empty()) ? &consistencyTable : nullptr,
                        option->consistencyWeight
                    );
                }
                */
                // FAMSA profile alignment DP kernel with consistency score
                const auto* cTable = (!consistencyTable.empty()) ? &consistencyTable : nullptr;
                float cWeight = option->consistencyWeight;
                aln_wo_gc = alignProfile_FAMSA(profRef, profQry, param, cTable, cWeight);
            }

            if (!aln_wo_gc.empty()) {
                alnPath aln_w_gc;
                int alnRef = 0, alnQry = 0;
                for (auto a : aln_wo_gc) {
                    if (a == 0) { alnRef += 1; alnQry += 1; }
                    if (a == 1) { alnQry += 1; }
                    if (a == 2) { alnRef += 1; }
                }
                alignment_helper::addGappyColumnsBack(aln_wo_gc, aln_w_gc, gappyColumns, param, {alnRef, alnQry}, consensus);

                float refWeight = nodes[nIdx].first->alnWeight;
                float qryWeight = nodes[nIdx].second->alnWeight;

                if (option->alnMode != PLACE_WO_TREE) {
                    alignment_helper::updateFrequency(nodes[nIdx], database, aln_w_gc, {refWeight, qryWeight});
                    alignment_helper::updateAlignment(nodes[nIdx], database, option, aln_w_gc);
                } else {
                    database->subtreeAln[nodes[nIdx].second->seqsIncluded[0]] = aln_w_gc;
                }

                // If running inter-subtree merge, update column provenance for the parent node
                if (subrootProvenance != nullptr) {
                    ColumnProvenance mergedProv(aln_w_gc.size());
                    int rCol = 0, qCol = 0;
                    for (size_t c = 0; c < aln_w_gc.size(); ++c) {
                        int8_t code = aln_w_gc[c];
                        if (code == 0) {
                            if (rCol < (int)refProvenance.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), refProvenance[rCol].begin(), refProvenance[rCol].end());
                            }
                            if (qCol < (int)qryProvenance.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), qryProvenance[qCol].begin(), qryProvenance[qCol].end());
                            }
                            rCol++; qCol++;
                        } else if (code == 1) {
                            if (qCol < (int)qryProvenance.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), qryProvenance[qCol].begin(), qryProvenance[qCol].end());
                            }
                            qCol++;
                        } else if (code == 2) {
                            if (rCol < (int)refProvenance.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), refProvenance[rCol].begin(), refProvenance[rCol].end());
                            }
                            rCol++;
                        }
                    }
                    Node* parentNode = nodes[nIdx].first->parent ? nodes[nIdx].first->parent : nodes[nIdx].second->parent;
                    if (parentNode != nullptr) {
                        static std::mutex provMutex;
                        std::lock_guard<std::mutex> lock(provMutex);
                        (*subrootProvenance)[parentNode->identifier] = std::move(mergedProv);
                    }
                }
            }
        }
    });
}

// -----------------------------------------------------------------------------
// 3. Progressive alignment with consistency score guidance
// -----------------------------------------------------------------------------
void progressiveAlignGroupWithConsistency(
    HierarchyGroup& group,
    std::shared_ptr<LocalHomTable> localHomTable,
    Tree* T,
    SequenceDB* database,
    Option* option,
    Params& param,
    alnFunction alignmentKernel,
    std::unordered_map<std::string, ColumnProvenance>* subrootProvenance)
{
    if (option->printDetail) {
        std::cerr << "  [Step 3/3] Progressive alignment with consistency bonus for Group " << group.groupID << "...\n";
    }

    // Wrap the CPU-based accurate alignment kernel with captured localHomTable and subrootProvenance
    alnFunction accurateAlnKernel = [localHomTable, subrootProvenance](Tree* tree, NodePairVec& nodes, SequenceDB* database, Option* option, Params& param) {
        alignmentKernel_Accurate_CPU(tree, nodes, database, option, param, localHomTable, subrootProvenance);
    };

    // Reuse msa::progressive::msaOnSubtree with accurateAlnKernel
    int subtree = (group.layer == 1 && !group.subtreeIDs.empty()) ? group.subtreeIDs[0] : -1;
    msa::progressive::msaOnSubtree(T, database, option, param, accurateAlnKernel, subtree);
}

// -----------------------------------------------------------------------------
// Orchestrator: Main progressive alignment pipeline for a consistency group
// -----------------------------------------------------------------------------
void msaOnGroup_accurate(
    HierarchyGroup& group,
    Tree* T,
    SequenceDB* database,
    Option* option,
    Params& param,
    alnFunction alignmentKernel,
    const std::vector<HierarchyGroup>* prevLayerGroups,
    std::unordered_map<std::string, ColumnProvenance>* subrootProvenance)
{
    auto groupStart = std::chrono::high_resolution_clock::now();
    if (option->printDetail) {
        std::cerr << "\n>>> Starting Progressive Alignment on Group " << group.groupID 
                  << " (Layer " << group.layer << ", Subtrees: " << group.subtreeIDs.size()
                  << ", Seqs: " << group.totalSequences << ", Reps: " << group.repCount << ") <<<\n";
    }

    // 1. All-to-all alignment for consistency score (reusing intra-group pairwise results from prevLayerGroups if available)
    auto localHomTable = alignGroupPairwiseAllToAll(group, T, database, option, param, prevLayerGroups);

    // 2. Build consistency score library & calculate centrality scores
    std::vector<float> centralityScores = buildGroupConsistencyLibrary(group, localHomTable, T, database, option, param);

    // 3. Progressive alignment with consistency score
    progressiveAlignGroupWithConsistency(group, localHomTable, T, database, option, param, alignmentKernel, subrootProvenance);

    // 4. In Layer 1, extract representative sequences using centrality scores
    if (group.layer == 1 && T && T->root) {
        int subtree = (!group.subtreeIDs.empty()) ? group.subtreeIDs[0] : -1;
        int alnLen = T->root->getAlnLen(database->currentTask);
        group.representatives = extractRepresentativesByCentrality(
            T->root, group.nextRepCount, database, subtree, centralityScores, *localHomTable, alnLen
        );
        for (auto& rep : group.representatives) {
            rep.originGroupID = group.groupID;
        }
        group.repSeqNames.clear();
        group.repSeqIndices.clear();
        for (const auto& rep : group.representatives) {
            group.repSeqNames.push_back(rep.name);
            group.repSeqIndices.push_back(rep.globalSeqID);
        }
    }

    // 5. Cache intra-representative pairwise results for higher layers (Pass-Forward Cache)
    if (!group.representatives.empty()) {
        group.repPairwiseCache.clear();
        for (size_t i = 0; i < group.representatives.size(); ++i) {
            int g1 = group.representatives[i].globalSeqID;
            int loc1 = localHomTable->getLocalIndex(g1);
            if (loc1 < 0) continue;
            for (size_t j = i + 1; j < group.representatives.size(); ++j) {
                int g2 = group.representatives[j].globalSeqID;
                if (g1 == g2) continue;
                int loc2 = localHomTable->getLocalIndex(g2);
                if (loc2 < 0 || loc1 == loc2) continue;
                uint64_t key = makeGlobalPairKey(g1, g2);
                group.repPairwiseCache[key] = localHomTable->getPairResult(loc1, loc2);
            }
        }
        if (option->printDetail) {
            std::cerr << "  Cached " << group.repPairwiseCache.size() 
                      << " intra-representative pairs for Group " << group.groupID << ".\n";
        }
    }

    auto groupEnd = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds groupElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(groupEnd - groupStart);
    if (option->printDetail) {
        std::cerr << ">>> Finished Group " << group.groupID << " in " << groupElapsed.count() << " ms <<<\n";
    }
}

// -----------------------------------------------------------------------------
// Top-Level Hierarchy Executor: Runs progressive consistency alignment across layers
// -----------------------------------------------------------------------------
void executeHierarchyPlan(HierarchyPlan& plan, Tree* T, PartitionInfo* P, Tree* subRoot_T, 
                          SequenceDB* database, Option* option, Params& param, alnFunction alignmentKernel) {
    if (option->printDetail) {
        std::cerr << "\n================================================================================\n";
        std::cerr << "[Accurate Mode] Executing Hierarchical Alignment Plan (Total Layers: " << plan.totalLayers << ")\n";
        std::cerr << "================================================================================\n";
    }

    if (plan.layers.empty()) return;

    // -------------------------------------------------------------------------
    // Layer 1: Base Subtree Progressive Alignment with Intra-Subtree Consistency
    // -------------------------------------------------------------------------
    auto& layer1 = plan.layers[0];
    int proceeded = 0;
    auto alnSubtreeStart = std::chrono::high_resolution_clock::now();

    bool inMemoryMode = (option->maxSubtree == INT32_MAX);
    if (inMemoryMode && database->sequences.empty()) {
        msa::io::readSequences(option->seqFile, database, option, T, -1);
    }

    struct BatchRange {
        size_t start;
        size_t end;
    };

    auto computeBalancedBatches = [](size_t totalGroups, size_t maxBatch) -> std::vector<BatchRange> {
        std::vector<BatchRange> batches;
        if (totalGroups == 0) return batches;
        if (totalGroups <= maxBatch) {
            batches.push_back({0, totalGroups});
            return batches;
        }
        size_t K = (totalGroups + maxBatch - 1) / maxBatch;
        size_t baseSize = totalGroups / K;
        size_t remainder = totalGroups % K;

        size_t currentStart = 0;
        for (size_t b = 0; b < K; ++b) {
            size_t bSize = baseSize + (b < remainder ? 1 : 0);
            batches.push_back({currentStart, currentStart + bSize});
            currentStart += bSize;
        }
        return batches;
    };

    struct GroupBatchContext {
        HierarchyGroup* group;
        phylogeny::Node* subRootNode;
        phylogeny::Tree* subT;
        int subtree;
        std::shared_ptr<LocalHomTable> localHomTable;
        std::vector<float> centralityScores;
    };

    if (!inMemoryMode) {
        // Disk mode: Process subtrees sequentially due to shared disk/database clearing
        for (auto& group : layer1) {
            auto subtreeStart = std::chrono::high_resolution_clock::now();
            ++proceeded;
            const std::string& rootID = group.rootIdentifier;

            int subtree = -1;
            phylogeny::Node* subRootNode = nullptr;
            if (P && P->partitionsRoot.find(rootID) != P->partitionsRoot.end()) {
                subRootNode = P->partitionsRoot.at(rootID).first;
                subtree = (P->partitionsRoot.size() > 1 && T && T->allNodes.find(rootID) != T->allNodes.end()) 
                          ? T->allNodes[rootID]->grpID : -1;
            } else if (T && T->allNodes.find(rootID) != T->allNodes.end()) {
                subRootNode = T->allNodes[rootID];
                subtree = (T->allNodes.size() > 1) ? T->allNodes[rootID]->grpID : -1;
            } else if (T && T->root) {
                subRootNode = T->root;
                subtree = 0;
            }

            if (!subRootNode) continue;

            if (P && P->partitionsRoot.size() > 1) {
                std::cerr << "Start processing subalignment No. " << subtree << ". (" 
                          << proceeded << '/' << layer1.size() << ")\n";
            }

            phylogeny::Tree* subT = new phylogeny::Tree(subRootNode, option->reroot);
            msa::io::readSequences(option->seqFile, database, option, subT, subtree);

            msaOnGroup_accurate(group, subT, database, option, param, alignmentKernel);

            if (option->debug) database->debug();

            if (P && P->partitionsRoot.size() > 1) {
                auto storeStart = std::chrono::high_resolution_clock::now();
                database->storeSubtreeProfile(subT, option->type, subtree);
                msa::io::writeSubAlignments(database, option, subtree, subT->root->getAlnLen(database->currentTask));
                std::string subNodeID = (subRoot_T && subRoot_T->allNodes.find(rootID) != subRoot_T->allNodes.end()) ? 
                                        rootID : (subRootNode && subRoot_T->allNodes.find(subRootNode->identifier) != subRoot_T->allNodes.end() ? 
                                                  subRootNode->identifier : subT->root->identifier);
                if (subRoot_T && subRoot_T->allNodes.find(subNodeID) != subRoot_T->allNodes.end()) {
                    phylogeny::updateSubrootInfo(subRoot_T->allNodes[subNodeID], subT, subtree);
                }
                group.rootIdentifier = subNodeID;
                database->cleanSubtreeDB();
                auto storeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                if (option->printDetail) {
                    std::cerr << "Stored the subalignments in " << storeTime.count() / 1000000 << " ms.\n";
                }
            } else {
                auto outStart = std::chrono::high_resolution_clock::now();
                msa::io::writeFinalMSA(database, option, subT->root->getAlnLen(database->currentTask));
                auto outEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds outTime = outEnd - outStart;
                std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
                std::cerr << "Wrote alignment to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
            }

            delete subT;

            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            if (P && P->partitionsRoot.size() > 1) {
                std::cerr << "Finished subalignment No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
            } else {
                std::cerr << "Finished the alignment in " << subtreeTime.count() / 1000000000 << " s\n";
            }
        }
    } else {
        // In-memory mode: Balanced batched execution with group-level parallel progressive alignment
        const size_t TARGET_BATCH_GROUPS = 1000;
        std::vector<BatchRange> batches = computeBalancedBatches(layer1.size(), TARGET_BATCH_GROUPS);

        for (size_t b = 0; b < batches.size(); ++b) {
            auto batchStart = std::chrono::high_resolution_clock::now();
            size_t bStart = batches[b].start;
            size_t bEnd = batches[b].end;
            size_t bCount = bEnd - bStart;

            if (option->printDetail || batches.size() > 1) {
                std::cerr << "\n[Accurate Mode] Processing Layer 1 Batch " << (b + 1) << "/" << batches.size() 
                          << " (Groups " << bStart << ".." << (bEnd - 1) << ", Total " << bCount << " groups)...\n";
            }

            // 1. Prepare contexts for groups in this batch
            std::vector<GroupBatchContext> batchContexts;
            batchContexts.reserve(bCount);

            for (size_t g = bStart; g < bEnd; ++g) {
                auto& group = layer1[g];
                const std::string& rootID = group.rootIdentifier;

                int subtree = -1;
                phylogeny::Node* subRootNode = nullptr;
                if (P && P->partitionsRoot.find(rootID) != P->partitionsRoot.end()) {
                    subRootNode = P->partitionsRoot.at(rootID).first;
                    subtree = (P->partitionsRoot.size() > 1 && T && T->allNodes.find(rootID) != T->allNodes.end()) 
                              ? T->allNodes[rootID]->grpID : -1;
                } else if (T && T->allNodes.find(rootID) != T->allNodes.end()) {
                    subRootNode = T->allNodes[rootID];
                    subtree = (T->allNodes.size() > 1) ? T->allNodes[rootID]->grpID : -1;
                } else if (T && T->root) {
                    subRootNode = T->root;
                    subtree = 0;
                }

                if (!subRootNode) continue;

                phylogeny::Tree* subT = new phylogeny::Tree(subRootNode, option->reroot);
                batchContexts.push_back({&group, subRootNode, subT, subtree, nullptr, {}});
            }

            // Phase 1: All-to-all alignment and consistency library generation for all groups in this batch
            for (auto& ctx : batchContexts) {
                ctx.localHomTable = alignGroupPairwiseAllToAll(*ctx.group, ctx.subT, database, option, param, nullptr);
                ctx.centralityScores = buildGroupConsistencyLibrary(*ctx.group, ctx.localHomTable, ctx.subT, database, option, param);
            }

            // Phase 2: Group-Level Parallel Progressive Alignment
            // Saturate all CPU threads by aligning multiple groups concurrently
            tbb::this_task_arena::isolate([&] {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, batchContexts.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t i = r.begin(); i < r.end(); ++i) {
                        auto& ctx = batchContexts[i];
                        progressiveAlignGroupWithConsistency(
                            *ctx.group, ctx.localHomTable, ctx.subT, database, option, param, alignmentKernel, nullptr
                        );
                    }
                });
            });

            // Phase 3: Post-processing (representatives, subroot updates, and memory release)
            for (auto& ctx : batchContexts) {
                ++proceeded;

                // Representative extraction using centrality scores
                if (ctx.group->layer == 1 && ctx.subT && ctx.subT->root) {
                    int alnLen = ctx.subT->root->getAlnLen(database->currentTask);
                    ctx.group->representatives = extractRepresentativesByCentrality(
                        ctx.subT->root, ctx.group->nextRepCount, database, ctx.subtree, ctx.centralityScores, *ctx.localHomTable, alnLen
                    );
                    for (auto& rep : ctx.group->representatives) {
                        rep.originGroupID = ctx.group->groupID;
                    }
                    ctx.group->repSeqNames.clear();
                    ctx.group->repSeqIndices.clear();
                    for (const auto& rep : ctx.group->representatives) {
                        ctx.group->repSeqNames.push_back(rep.name);
                        ctx.group->repSeqIndices.push_back(rep.globalSeqID);
                    }
                }

                // Cache intra-representative pairwise results for higher layers
                if (!ctx.group->representatives.empty()) {
                    ctx.group->repPairwiseCache.clear();
                    for (size_t i = 0; i < ctx.group->representatives.size(); ++i) {
                        int g1 = ctx.group->representatives[i].globalSeqID;
                        int loc1 = ctx.localHomTable->getLocalIndex(g1);
                        if (loc1 < 0) continue;
                        for (size_t j = i + 1; j < ctx.group->representatives.size(); ++j) {
                            int g2 = ctx.group->representatives[j].globalSeqID;
                            if (g1 == g2) continue;
                            int loc2 = ctx.localHomTable->getLocalIndex(g2);
                            if (loc2 < 0 || loc1 == loc2) continue;
                            uint64_t key = makeGlobalPairKey(g1, g2);
                            ctx.group->repPairwiseCache[key] = ctx.localHomTable->getPairResult(loc1, loc2);
                        }
                    }
                }

                // Update subRoot_T's leaf node with the group's sequences and profile
                if (P && P->partitionsRoot.size() > 1) {
                    std::vector<std::string> leafNames;
                    collectLeavesInOrder(ctx.subT->root, leafNames);
                    std::vector<int> groupSeqIDs;
                    groupSeqIDs.reserve(leafNames.size());
                    for (const auto& name : leafNames) {
                        auto it = database->name_map.find(name);
                        if (it != database->name_map.end()) {
                            groupSeqIDs.push_back(it->second->id);
                        }
                    }

                    const std::string& rootID = ctx.group->rootIdentifier;
                    std::string subNodeID = (subRoot_T && subRoot_T->allNodes.find(rootID) != subRoot_T->allNodes.end()) ? 
                                            rootID : (ctx.subRootNode && subRoot_T->allNodes.find(ctx.subRootNode->identifier) != subRoot_T->allNodes.end() ? 
                                                      ctx.subRootNode->identifier : ctx.subT->root->identifier);
                    if (subRoot_T && subRoot_T->allNodes.find(subNodeID) != subRoot_T->allNodes.end()) {
                        Node* subLeaf = subRoot_T->allNodes[subNodeID];
                        subLeaf->seqsIncluded = std::move(groupSeqIDs);
                        subLeaf->alnLen = ctx.subT->root->alnLen;
                        subLeaf->alnNum = ctx.subT->root->alnNum;
                        subLeaf->alnWeight = ctx.subT->root->alnWeight;
                        if (ctx.subT->root->msaFreq.empty()) {
                            int refLen = ctx.subT->root->alnLen;
                            int profileSize = (option->type == 'n') ? 6 : 22;
                            ctx.subT->root->msaFreq.assign(refLen, std::vector<float>(profileSize, 0.0f));
                            for (int sIdx : subLeaf->seqsIncluded) {
                                float w = database->sequences[sIdx]->weight;
                                int storage = database->sequences[sIdx]->storage;
                                for (int t = 0; t < refLen; ++t) {
                                    int letterIndex = letterIdx(option->type, toupper(database->sequences[sIdx]->alnStorage[storage][t]));
                                    ctx.subT->root->msaFreq[t][letterIndex] += w;
                                }
                            }
                        }
                        subLeaf->msaFreq = ctx.subT->root->msaFreq;
                    }
                    ctx.group->rootIdentifier = subNodeID;
                } else {
                    // Single tree alignment finished
                    auto outStart = std::chrono::high_resolution_clock::now();
                    msa::io::writeFinalMSA(database, option, ctx.subT->root->getAlnLen(database->currentTask));
                    auto outEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds outTime = outEnd - outStart;
                    std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
                    std::cerr << "Wrote alignment to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
                }

                // Immediately release consistency table memory to bound RAM usage
                ctx.localHomTable.reset();
                ctx.centralityScores.clear();
                delete ctx.subT;
            }

            if (option->debug) database->debug();

            auto batchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::milliseconds batchElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(batchEnd - batchStart);
            if (P && P->partitionsRoot.size() > 1) {
                std::cerr << "[Batch " << (b + 1) << "/" << batches.size() 
                          << "] Finished " << bCount << " subalignments (" 
                          << proceeded << "/" << layer1.size() << ") in " 
                          << batchElapsed.count() / 1000.0f << " s.\n";
            } else {
                std::cerr << "Finished the alignment in " << batchElapsed.count() / 1000.0f << " s.\n";
            }
        }
    }

    if (!P || P->partitionsRoot.size() <= 1 || plan.layers.size() <= 1) {
        return;
    }

    auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
    std::cerr << "Finished all subalignments in " << alnSubtreeTime.count() / 1000000000 << " s.\n";

    // -------------------------------------------------------------------------
    // Layer 2+: Inter-Subtree Hierarchical Merges
    // -------------------------------------------------------------------------
    for (size_t l = 1; l < plan.layers.size(); ++l) {
        auto& currentLayer = plan.layers[l];
        const auto& prevLayer = plan.layers[l - 1];

        if (option->printDetail) {
            std::cerr << "\n--------------------------------------------------------------------------------\n";
            std::cerr << "[Accurate Mode] Starting Layer " << (l + 1) << " (" << currentLayer.size() << " merge group(s))\n";
            std::cerr << "--------------------------------------------------------------------------------\n";
        }

        for (auto& group : currentLayer) {
            // Aggregate representative sequences from child groups
            group.representatives.clear();
            group.repSeqNames.clear();
            group.repSeqIndices.clear();

            for (int childID : group.childGroupIDs) {
                if (childID >= 0 && childID < static_cast<int>(prevLayer.size())) {
                    const auto& child = prevLayer[childID];
                    for (auto rep : child.representatives) {
                        rep.originGroupID = childID;
                        group.representatives.push_back(std::move(rep));
                    }
                }
            }

            group.repCount = group.representatives.size();
            for (const auto& rep : group.representatives) {
                group.repSeqNames.push_back(rep.name);
                group.repSeqIndices.push_back(rep.globalSeqID);
            }

            if (group.childGroupIDs.size() == 1) {
                // Trivial pass-through: single child group
                int childID = group.childGroupIDs[0];
                if (childID >= 0 && childID < static_cast<int>(prevLayer.size())) {
                    group.rootIdentifier = prevLayer[childID].rootIdentifier;
                }
            } else {
                // Extract induced subtree for this group's clade
                std::unordered_map<std::string, std::pair<phylogeny::Node*, size_t>> targetTips;
                for (int childID : group.childGroupIDs) {
                    if (childID >= 0 && childID < static_cast<int>(prevLayer.size())) {
                        const std::string& tipID = prevLayer[childID].rootIdentifier;
                        if (subRoot_T->allNodes.find(tipID) != subRoot_T->allNodes.end()) {
                            targetTips[tipID] = std::make_pair(subRoot_T->allNodes[tipID], 1);
                        }
                    }
                }

                phylogeny::Node* cladeRootInSubRoot = nullptr;
                auto itRoot = subRoot_T->allNodes.find(group.rootIdentifier);
                if (itRoot != subRoot_T->allNodes.end()) {
                    cladeRootInSubRoot = itRoot->second;
                } else {
                    cladeRootInSubRoot = subRoot_T->root;
                }

                phylogeny::Tree* group_T = new phylogeny::Tree();
                group_T->root = phylogeny::buildInducedTree(cladeRootInSubRoot, targetTips, group_T);
                if (group_T->root) {
                    group_T->root->parent = nullptr;
                }

                // Copy alignment profile and sequence membership to the leaf tips of group_T
                for (auto& kv : group_T->allNodes) {
                    phylogeny::Node* n = kv.second;
                    if (targetTips.find(n->identifier) != targetTips.end()) {
                        phylogeny::Node* orig = subRoot_T->allNodes[n->identifier];
                        n->seqsIncluded = orig->seqsIncluded;
                        n->alnLen = orig->alnLen;
                        n->alnNum = orig->alnNum;
                        n->alnWeight = orig->alnWeight;
                        n->msaFreq = orig->msaFreq;
                    }
                }

                // Perform consistency-guided progressive alignment on this group's induced subtree
                database->currentTask = inMemoryMode ? 0 : 2;
                msaOnGroup_accurate(group, group_T, database, option, param, alignmentKernel, &prevLayer, nullptr);

                // Copy merged alignment attributes back to clade root in subRoot_T
                if (group_T->root) {
                    phylogeny::Node* actualRootInSubRoot = subRoot_T->allNodes[group_T->root->identifier];
                    actualRootInSubRoot->seqsIncluded = std::move(group_T->root->seqsIncluded);
                    actualRootInSubRoot->msaFreq = std::move(group_T->root->msaFreq);
                    actualRootInSubRoot->alnLen = group_T->root->alnLen;
                    actualRootInSubRoot->alnNum = group_T->root->alnNum;
                    actualRootInSubRoot->alnWeight = group_T->root->alnWeight;
                    group.rootIdentifier = actualRootInSubRoot->identifier;
                }

                delete group_T;
            }

            // If not top layer, sample representatives for the next level
            if (l + 1 < plan.layers.size() && group.nextRepCount < group.representatives.size()) {
                std::vector<Representative> nextReps;
                double stride = static_cast<double>(group.representatives.size()) / static_cast<double>(group.nextRepCount);
                for (size_t r = 0; r < group.nextRepCount; ++r) {
                    size_t idx = static_cast<size_t>(r * stride);
                    if (idx < group.representatives.size()) {
                        nextReps.push_back(std::move(group.representatives[idx]));
                    }
                }
                group.representatives = std::move(nextReps);
            }
        }
    }

    // Determine final alignment length from top-layer root
    int finalAlnLen = 0;
    if (!plan.layers.empty() && !plan.layers.back().empty()) {
        const std::string& topRootID = plan.layers.back()[0].rootIdentifier;
        auto it = subRoot_T->allNodes.find(topRootID);
        if (it != subRoot_T->allNodes.end()) {
            finalAlnLen = it->second->getAlnLen(inMemoryMode ? 0 : database->currentTask);
        }
    }
    if (finalAlnLen == 0 && subRoot_T && subRoot_T->root) {
        finalAlnLen = subRoot_T->root->getAlnLen(inMemoryMode ? 0 : database->currentTask);
    }

    // Write final merged MSA
    if (inMemoryMode) {
        database->currentTask = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        msa::io::writeFinalMSA(database, option, finalAlnLen);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
        std::cerr << "Wrote final alignment (total " << database->sequences.size() 
                  << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
    } else {
        database->currentTask = 2;
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        msa::io::update_and_writeAlignments(database, option, totalSeqs);
        msa::io::writeFinalMSA(database, option, finalAlnLen);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
        std::cerr << "Wrote " << subRoot_T->allNodes.size() << " subalignments (total " << totalSeqs 
                  << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
    }
}

} // namespace accurate
} // namespace msa
