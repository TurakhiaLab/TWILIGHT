#include "msa.hpp"
#include "align-method.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>

namespace msa {
namespace accurate {

void msaOnSubtree_accurate(Tree *T, SequenceDB *database, Option *option, Params &param, alnFunction alignmentKernel, int subtree) {
    if (!alignmentKernel) {
        alignmentKernel = [](Tree* t, NodePairVec& n, SequenceDB* d, Option* o, Params& p) {
            msa::accurate::alignmentKernel_Accurate_CPU(t, n, d, o, p);
        };
    }

    int seqCount = (T && T->root) ? T->root->numLeaves : (int)database->sequences.size();
    option->resolveAccurateParams(seqCount);

    // If leaf count is within accMaxGroup, align as a single group with consistency table
    if (T->root && T->root->numLeaves <= (size_t)option->accMaxGroup) {
        HierarchyGroup group;
        group.layer = 1;
        group.groupID = (subtree >= 0) ? subtree : 0;
        group.rootIdentifier = T->root->identifier;
        if (subtree >= 0) group.subtreeIDs = { subtree };
        group.totalSequences = T->root->numLeaves;
        group.repCount = T->root->numLeaves;
        group.nextRepCount = std::max((size_t)1, (size_t)std::ceil(group.repCount * option->accSampleRate));

        auto singleGroupStart = std::chrono::high_resolution_clock::now();
        msaOnGroup_accurate(group, T, database, option, param, alignmentKernel);
        auto singleGroupEnd = std::chrono::high_resolution_clock::now();
        auto singleGroupMs = std::chrono::duration_cast<std::chrono::milliseconds>(singleGroupEnd - singleGroupStart).count();
        if (option->printDetail) {
            std::cerr << "Progressive alignment on 1 groups in " << singleGroupMs << " ms\n";
        }
        if (subtree >= 0 && !group.representatives.empty()) {
            database->subtreeRepresentatives[T->root->identifier] = group.representatives;
            database->subtreeRepresentatives[std::to_string(subtree)] = group.representatives;
        }
    } else {
        // Divide-and-conquer into inner groups to bound consistency score runtime within this subtree
        phylogeny::PartitionInfo* innerP = new phylogeny::PartitionInfo(option->accMaxGroup, 0, 0);
        innerP->partitionTree(T->root);
        phylogeny::Tree* innerSubRoot_T = phylogeny::constructTreeFromPartitions_new(T->root, innerP);

        auto hierarchyPlan = computeHierarchyPlan(T, innerP, innerSubRoot_T, option);
        if (option->printDetail) {
            hierarchyPlan.printSummary(option->printDetail);
        }

        executeHierarchyPlan(hierarchyPlan, T, innerP, innerSubRoot_T, database, option, param, alignmentKernel, false);

        if (subtree >= 0 && !hierarchyPlan.layers.empty() && !hierarchyPlan.layers.back().empty()) {
            auto topReps = hierarchyPlan.layers.back()[0].representatives;
            if (!topReps.empty()) {
                int finalAlnLen = T->root ? T->root->alnLen : 0;
                for (auto& rep : topReps) {
                    auto it = database->name_map.find(rep.name);
                    if (it != database->name_map.end()) {
                        const auto* sInfo = it->second;
                        const char* rawChars = sInfo->alnStorage[sInfo->storage];
                        rep.colInSubAln.clear();
                        for (int k = 0; k < sInfo->len && (finalAlnLen == 0 || k < finalAlnLen); ++k) {
                            if (rawChars[k] != '-' && rawChars[k] != '.' && rawChars[k] != '\0') {
                                rep.colInSubAln.push_back(k);
                            }
                        }
                    }
                }
                database->subtreeRepresentatives[T->root->identifier] = topReps;
                database->subtreeRepresentatives[std::to_string(subtree)] = topReps;
            }
        }

        delete innerSubRoot_T;
        delete innerP;
    }
}


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
        group.representatives = extractRepresentatives(
            T->root, group.nextRepCount, database, group.groupID, centralityScores, *localHomTable, alnLen, option, nullptr, subtree
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
    }
}



// Top-Level Hierarchy Executor: Runs progressive consistency alignment across layers
// -----------------------------------------------------------------------------
void executeHierarchyPlan(HierarchyPlan& plan, Tree* T, PartitionInfo* P, Tree* subRoot_T, 
                          SequenceDB* database, Option* option, Params& param, alnFunction alignmentKernel,
                          bool writeOutput) {
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

    if (option->printDetail) {
        std::cerr << "--------------------------------------------------------------------------------\n";
        std::cerr << "[Accurate Mode] Starting Layer 1 (" << layer1.size() << " base group(s))\n";
        std::cerr << "--------------------------------------------------------------------------------\n";
    }

    bool inMemoryMode = !database->sequences.empty() || (option->maxSubtree == INT32_MAX);
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

            // Disable verbose tree depth printing for accurate mode per-group subtrees
            phylogeny::Tree* subT = new phylogeny::Tree(subRootNode, option->reroot, false);
            msa::io::readSequences(option->seqFile, database, option, subT, subtree);

            msaOnGroup_accurate(group, subT, database, option, param, alignmentKernel);

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
        auto layer1Start = std::chrono::high_resolution_clock::now();
        const size_t TARGET_BATCH_GROUPS = 1000;
        std::vector<BatchRange> batches = computeBalancedBatches(layer1.size(), TARGET_BATCH_GROUPS);

        for (size_t b = 0; b < batches.size(); ++b) {
            auto batchStart = std::chrono::high_resolution_clock::now();
            size_t bStart = batches[b].start;
            size_t bEnd = batches[b].end;
            size_t bCount = bEnd - bStart;

            if (batches.size() > 1) {
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

                // Disable verbose tree depth printing for accurate mode per-group subtrees
                phylogeny::Tree* subT = new phylogeny::Tree(subRootNode, option->reroot, false);
                if (subT->root && (subT->root->is_leaf() || subT->root->numLeaves <= 1)) {
                    auto itSeq = database->name_map.find(subT->root->identifier);
                    if (itSeq != database->name_map.end()) {
                        subT->root->seqsIncluded = { itSeq->second->id };
                        subT->root->alnLen = itSeq->second->len;
                        subT->root->alnNum = 1;
                        subT->root->alnWeight = itSeq->second->weight;
                    }
                }
                batchContexts.push_back({&group, subRootNode, subT, subtree, nullptr, {}});
            }

            // Phase 1: All-to-all alignment and consistency library generation for all groups in this batch
            auto phase1Start = std::chrono::high_resolution_clock::now();
            for (auto& ctx : batchContexts) {
                ctx.localHomTable = alignGroupPairwiseAllToAll(*ctx.group, ctx.subT, database, option, param, nullptr);
                ctx.centralityScores = buildGroupConsistencyLibrary(*ctx.group, ctx.localHomTable, ctx.subT, database, option, param);
            }
            auto phase1End = std::chrono::high_resolution_clock::now();
            auto phase1Ms = std::chrono::duration_cast<std::chrono::milliseconds>(phase1End - phase1Start).count();
            if (option->printDetail) {
                std::cerr << "Consistency library built on " << batchContexts.size() << " groups in " << phase1Ms << " ms\n";
            }

            // Phase 2: Group-Level Parallel Progressive Alignment
            // Saturate all CPU threads by aligning multiple groups concurrently
            auto alignGroupsStart = std::chrono::high_resolution_clock::now();
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
            auto alignGroupsEnd = std::chrono::high_resolution_clock::now();
            auto alignGroupsMs = std::chrono::duration_cast<std::chrono::milliseconds>(alignGroupsEnd - alignGroupsStart).count();
            if (option->printDetail) {
                std::cerr << "Progressive alignment on " << batchContexts.size() << " groups in " << alignGroupsMs << " ms\n";
            }

            // Phase 3: Post-processing (representatives, subroot updates, and memory release)
            auto samplingStart = std::chrono::high_resolution_clock::now();
            std::vector<std::string> samplingLogs(batchContexts.size());

            // Parallel sampling and intra-representative cache computation across groups in this batch
            tbb::this_task_arena::isolate([&] {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, batchContexts.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t i = r.begin(); i < r.end(); ++i) {
                        auto& ctx = batchContexts[i];
                        if (ctx.group->layer == 1 && ctx.subT && ctx.subT->root) {
                            if (ctx.subT->root->alnNum == 0) {
                                if (!ctx.subT->root->seqsIncluded.empty()) {
                                    int s0 = ctx.subT->root->seqsIncluded[0];
                                    ctx.subT->root->alnLen = database->sequences[s0]->len;
                                    ctx.subT->root->alnNum = ctx.subT->root->seqsIncluded.size();
                                    ctx.subT->root->alnWeight = database->sequences[s0]->weight;
                                } else {
                                    auto itSeq = database->name_map.find(ctx.subT->root->identifier);
                                    if (itSeq != database->name_map.end()) {
                                        ctx.subT->root->seqsIncluded = { itSeq->second->id };
                                        ctx.subT->root->alnLen = itSeq->second->len;
                                        ctx.subT->root->alnNum = 1;
                                        ctx.subT->root->alnWeight = itSeq->second->weight;
                                    }
                                }
                            }
                            int alnLen = ctx.subT->root->getAlnLen(database->currentTask);
                            std::string logMsg;
                            ctx.group->representatives = extractRepresentatives(
                                ctx.subT->root, ctx.group->nextRepCount, database, ctx.group->groupID, ctx.centralityScores, *ctx.localHomTable, alnLen, option, &logMsg, ctx.subtree
                            );
                            samplingLogs[i] = std::move(logMsg);
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
                            for (size_t a = 0; a < ctx.group->representatives.size(); ++a) {
                                int g1 = ctx.group->representatives[a].globalSeqID;
                                int loc1 = ctx.localHomTable->getLocalIndex(g1);
                                if (loc1 < 0) continue;
                                for (size_t b = a + 1; b < ctx.group->representatives.size(); ++b) {
                                    int g2 = ctx.group->representatives[b].globalSeqID;
                                    if (g1 == g2) continue;
                                    int loc2 = ctx.localHomTable->getLocalIndex(g2);
                                    if (loc2 < 0 || loc1 == loc2) continue;
                                    uint64_t key = makeGlobalPairKey(g1, g2);
                                    ctx.group->repPairwiseCache[key] = ctx.localHomTable->getPairResult(loc1, loc2);
                                }
                            }
                        }
                    }
                });
            });
            auto samplingEnd = std::chrono::high_resolution_clock::now();
            auto samplingMs = std::chrono::duration_cast<std::chrono::milliseconds>(samplingEnd - samplingStart).count();

            // Print sampling logs strictly in groupID order
            if (option->printDetail) {
                for (const auto& log : samplingLogs) {
                    if (!log.empty()) std::cerr << log;
                }
                std::cerr << "Sampling on " << batchContexts.size() << " groups in " << samplingMs << " ms\n";
            }

            for (auto& ctx : batchContexts) {
                ++proceeded;

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
                        if (ctx.subT->root->alnNum == 0 && !subLeaf->seqsIncluded.empty()) {
                            int s0 = subLeaf->seqsIncluded[0];
                            ctx.subT->root->alnLen = database->sequences[s0]->len;
                            ctx.subT->root->alnNum = subLeaf->seqsIncluded.size();
                            ctx.subT->root->alnWeight = database->sequences[s0]->weight;
                        }
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

            auto batchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::milliseconds batchElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(batchEnd - batchStart);
            if (batches.size() > 1) {
                std::cerr << "Layer 1 Batch " << (b + 1) << "/" << batches.size() 
                          << ": Finished " << bCount << " groups in " 
                          << std::fixed << std::setprecision(2) << (batchElapsed.count() / 1000.0f) << " s\n";
            }
        }
        auto layer1End = std::chrono::high_resolution_clock::now();
        std::chrono::milliseconds layer1Elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(layer1End - layer1Start);
        std::cerr << "Layer 1: Finished " << layer1.size() << " groups in " 
                  << std::fixed << std::setprecision(2) << (layer1Elapsed.count() / 1000.0f) << " s\n";
    }

    if (!P || P->partitionsRoot.size() <= 1 || plan.layers.size() <= 1) {
        return;
    }



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

        auto layerStart = std::chrono::high_resolution_clock::now();

        struct LayerGroupContext {
            HierarchyGroup* group;
            phylogeny::Tree* group_T = nullptr;
            std::shared_ptr<LocalHomTable> localHomTable;
            std::vector<float> centralityScores;
        };

        // Step 1: Pre-aggregate representative sequences from child groups
        for (auto& group : currentLayer) {
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
        }

        // Step 2: Batched execution for groups in this layer to bound peak memory usage
        const size_t TARGET_BATCH_GROUPS = 1000;
        std::vector<BatchRange> batches = computeBalancedBatches(currentLayer.size(), TARGET_BATCH_GROUPS);

        for (size_t b = 0; b < batches.size(); ++b) {
            auto batchStart = std::chrono::high_resolution_clock::now();
            size_t bStart = batches[b].start;
            size_t bEnd = batches[b].end;
            size_t bCount = bEnd - bStart;

            if (batches.size() > 1) {
                std::cerr << "\n[Accurate Mode] Processing Layer " << (l + 1) << " Batch " << (b + 1) << "/" << batches.size() 
                          << " (Groups " << bStart << ".." << (bEnd - 1) << ", Total " << bCount << " groups)...\n";
            }

            std::vector<LayerGroupContext> batchContexts;
            batchContexts.reserve(bCount);

            for (size_t g = bStart; g < bEnd; ++g) {
                auto& group = currentLayer[g];
                if (group.childGroupIDs.size() == 1) {
                    // Trivial pass-through: single child group
                    int childID = group.childGroupIDs[0];
                    if (childID >= 0 && childID < static_cast<int>(prevLayer.size())) {
                        group.rootIdentifier = prevLayer[childID].rootIdentifier;
                    }
                    batchContexts.push_back({&group, nullptr, nullptr, {}});
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
                            if (orig->alnNum == 0 && !orig->seqsIncluded.empty()) {
                                int s0 = orig->seqsIncluded[0];
                                orig->alnNum = orig->seqsIncluded.size();
                                orig->alnLen = database->sequences[s0]->len;
                                orig->alnWeight = database->sequences[s0]->weight;
                            }
                            n->seqsIncluded = orig->seqsIncluded;
                            n->alnLen = orig->alnLen;
                            n->alnNum = orig->alnNum;
                            n->alnWeight = orig->alnWeight;
                            n->msaFreq = orig->msaFreq;
                        }
                    }

                    batchContexts.push_back({&group, group_T, nullptr, {}});
                }
            }

            // Phase 1: All-to-all alignment and consistency library generation
            auto phase1Start = std::chrono::high_resolution_clock::now();
            size_t nonTrivialCount = 0;
            for (auto& ctx : batchContexts) {
                if (!ctx.group_T) continue;
                ++nonTrivialCount;
                ctx.localHomTable = alignGroupPairwiseAllToAll(*ctx.group, ctx.group_T, database, option, param, &prevLayer);
                ctx.centralityScores = buildGroupConsistencyLibrary(*ctx.group, ctx.localHomTable, ctx.group_T, database, option, param);
            }
            auto phase1End = std::chrono::high_resolution_clock::now();
            auto phase1Ms = std::chrono::duration_cast<std::chrono::milliseconds>(phase1End - phase1Start).count();
            if (option->printDetail) {
                std::cerr << "Consistency library built on " << nonTrivialCount << " groups in " << phase1Ms << " ms\n";
            }

            // Phase 2: Group-level parallel progressive alignment
            database->currentTask = inMemoryMode ? 0 : 2;
            auto phase2Start = std::chrono::high_resolution_clock::now();
            tbb::this_task_arena::isolate([&] {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, batchContexts.size()), [&](const tbb::blocked_range<size_t>& r) {
                    for (size_t i = r.begin(); i < r.end(); ++i) {
                        auto& ctx = batchContexts[i];
                        if (!ctx.group_T) continue;
                        progressiveAlignGroupWithConsistency(
                            *ctx.group, ctx.localHomTable, ctx.group_T, database, option, param, alignmentKernel, nullptr
                        );
                    }
                });
            });
            auto phase2End = std::chrono::high_resolution_clock::now();
            auto phase2Ms = std::chrono::duration_cast<std::chrono::milliseconds>(phase2End - phase2Start).count();
            if (option->printDetail) {
                std::cerr << "Progressive alignment on " << nonTrivialCount << " groups in " << phase2Ms << " ms\n";
            }

            // Copy merged alignment attributes back to clade root in subRoot_T and clean up group_T
            for (auto& ctx : batchContexts) {
                if (!ctx.group_T) continue;
                if (ctx.group_T->root) {
                    phylogeny::Node* actualRootInSubRoot = subRoot_T->allNodes[ctx.group_T->root->identifier];
                    actualRootInSubRoot->seqsIncluded = std::move(ctx.group_T->root->seqsIncluded);
                    actualRootInSubRoot->msaFreq = std::move(ctx.group_T->root->msaFreq);
                    actualRootInSubRoot->alnLen = ctx.group_T->root->alnLen;
                    actualRootInSubRoot->alnNum = ctx.group_T->root->alnNum;
                    actualRootInSubRoot->alnWeight = ctx.group_T->root->alnWeight;
                    ctx.group->rootIdentifier = actualRootInSubRoot->identifier;
                }
                delete ctx.group_T;
                ctx.localHomTable.reset();
            }

            auto batchEnd = std::chrono::high_resolution_clock::now();
            std::chrono::milliseconds batchElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(batchEnd - batchStart);
            if (batches.size() > 1) {
                std::cerr << "Layer " << (l + 1) << " Batch " << (b + 1) << "/" << batches.size() 
                          << ": Finished " << bCount << " groups in " 
                          << std::fixed << std::setprecision(2) << (batchElapsed.count() / 1000.0f) << " s\n";
            }
        }

        // Phase 3: Sample representatives for next level (if not top layer)
        if (l + 1 < plan.layers.size()) {
            auto samplingStart = std::chrono::high_resolution_clock::now();
            size_t sampledCount = 0;
            for (auto& group : currentLayer) {
                if (group.nextRepCount < group.representatives.size()) {
                    std::vector<Representative> nextReps;
                    double stride = static_cast<double>(group.representatives.size()) / static_cast<double>(group.nextRepCount);
                    for (size_t r = 0; r < group.nextRepCount; ++r) {
                        size_t idx = static_cast<size_t>(r * stride);
                        if (idx < group.representatives.size()) {
                            nextReps.push_back(std::move(group.representatives[idx]));
                        }
                    }
                    group.representatives = std::move(nextReps);
                    ++sampledCount;
                }
            }
            auto samplingEnd = std::chrono::high_resolution_clock::now();
            auto samplingMs = std::chrono::duration_cast<std::chrono::milliseconds>(samplingEnd - samplingStart).count();
            if (option->printDetail && sampledCount > 0) {
                std::cerr << "Sampling on " << sampledCount << " groups in " << samplingMs << " ms\n";
            }
        }

        auto layerEnd = std::chrono::high_resolution_clock::now();
        std::chrono::milliseconds layerElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(layerEnd - layerStart);
        std::cerr << "Layer " << (l + 1) << ": Finished " << currentLayer.size() << " groups in " 
                  << std::fixed << std::setprecision(2) << (layerElapsed.count() / 1000.0f) << " s\n";
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
    if (T && T->root) {
        T->root->alnLen = finalAlnLen;
        T->root->alnNum = T->root->numLeaves;
        if (subRoot_T && subRoot_T->root) {
            T->root->alnWeight = subRoot_T->root->alnWeight;
            if (!subRoot_T->root->msaFreq.empty()) {
                T->root->msaFreq = subRoot_T->root->msaFreq;
            }
        }
        if (T->root->seqsIncluded.empty()) {
            for (const auto& s : database->sequences) {
                if (s) T->root->seqsIncluded.push_back(s->id);
            }
        }
    }

    if (inMemoryMode && T && T->root) {
        msa::progressive::updateAlignment(T->root, database);
    }

    // Write final merged MSA
    if (writeOutput) {
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
}


void mergeSubtreesAccurate(
    Tree* subRoot_T,
    SequenceDB* database,
    Option* option,
    Params& param)
{
    if (!subRoot_T || !subRoot_T->root) return;

    // Collect leaves of subRoot_T in phylogenetic order
    std::vector<std::string> tipNames;
    collectLeavesInOrder(subRoot_T->root, tipNames);

    // 1. Gather all representatives from the subtrees in subRoot_T
    HierarchyGroup group;
    group.layer = 2;
    group.groupID = 0;
    group.rootIdentifier = subRoot_T->root->identifier;

    size_t numTips = tipNames.size();
    size_t maxTotalReps = 500;
    size_t maxPerSubtree = (numTips > 0) ? std::max((size_t)20, maxTotalReps / numTips) : 100;

    for (const auto& tipName : tipNames) {
        auto it = database->subtreeRepresentatives.find(tipName);
        if (it != database->subtreeRepresentatives.end()) {
            size_t count = 0;
            for (const auto& rep : it->second) {
                if (count++ >= maxPerSubtree) break;
                group.representatives.push_back(rep);
            }
        } else {
            // Check by seqsIncluded
            auto itNode = subRoot_T->allNodes.find(tipName);
            if (itNode != subRoot_T->allNodes.end() && !itNode->second->seqsIncluded.empty()) {
                auto it2 = database->subtreeRepresentatives.find(std::to_string(itNode->second->seqsIncluded[0]));
                if (it2 != database->subtreeRepresentatives.end()) {
                    size_t count = 0;
                    for (const auto& rep : it2->second) {
                        if (count++ >= maxPerSubtree) break;
                        group.representatives.push_back(rep);
                    }
                }
            }
        }
    }

    if (group.representatives.empty()) {
        // Fallback to vanilla progressive alignment if no representatives found
        msa::progressive::msaOnSubtree(subRoot_T, database, option, param, msa::progressive::cpu::alignmentKernel_CPU);
        return;
    }

    group.totalSequences = group.representatives.size();
    group.repCount = group.representatives.size();
    for (const auto& rep : group.representatives) {
        group.repSeqNames.push_back(rep.name);
        group.repSeqIndices.push_back(rep.globalSeqID);
    }

    std::cerr << "\n>>> Starting Accurate Subtree Merge on " << tipNames.size()
              << " Subtrees (Total Reps: " << group.repCount << ") <<<\n";

    // 2. All-to-all alignment for consistency among representatives
    auto localHomTable = alignGroupPairwiseAllToAll(group, subRoot_T, database, option, param);

    // 3. Build consistency library & calculate centrality
    auto centralityScores = buildGroupConsistencyLibrary(group, localHomTable, subRoot_T, database, option, param);

    // 4. Progressive alignment on subRoot_T with consistency score guidance
    std::unordered_map<std::string, ColumnProvenance> subrootProvenance;
    progressiveAlignGroupWithConsistency(group, localHomTable, subRoot_T, database, option, param, nullptr, &subrootProvenance);

    std::cerr << ">>> Finished Accurate Subtree Merge <<<\n";
}

// -----------------------------------------------------------------------------
// Accurate Placement without Guide Tree: Streaming On-The-Fly Consistency
// -----------------------------------------------------------------------------
void placeWithoutTree_Accurate(Tree* T, SequenceDB* database, Option* option, Params& param) {
    auto placeStart = std::chrono::high_resolution_clock::now();

    // 1. Extract hybrid representatives from backbone alignment
    int targetReps = (option->accMaxGroup > 1) ? option->accMaxGroup : 100;
    std::vector<BackboneRep> backboneReps = extractBackboneRepresentatives(option->backboneAlnFile, targetReps, option);
    if (backboneReps.empty()) {
        std::cerr << "WARNING: No backbone representatives selected. Falling back to standard placement.\n";
        msa::progressive::msaOnSubtree(T, database, option, param, msa::progressive::cpu::alignmentKernel_CPU);
        return;
    }

    // 2. Collect query nodes (all children of T->root)
    std::vector<Node*> queryNodes;
    queryNodes.reserve(T->root->children.size());
    for (auto* child : T->root->children) {
        if (child) {
            if (child->seqsIncluded.empty() && database->name_map.find(child->identifier) != database->name_map.end()) {
                auto* s = database->name_map[child->identifier];
                child->seqsIncluded.push_back(s->id);
                child->alnLen = s->len;
                child->alnNum = 1;
                child->alnWeight = s->weight;
            }
            if (!child->seqsIncluded.empty()) {
                queryNodes.push_back(child);
            }
        }
    }

    const int totalQueries = static_cast<int>(queryNodes.size());
    const int bbLen = T->root->getAlnLen(database->currentTask);
    const int profileSize = (option->type == 'n') ? 6 : 22;

    std::atomic<int> completedQueries{0};

    // 3. Parallel Streamed Alignment: On-the-fly consistency table per query sequence
    tbb::parallel_for(tbb::blocked_range<int>(0, totalQueries), [&](const tbb::blocked_range<int>& range) {
        Aligner aligner; // Thread-local scratchpad for pairwise SW
        for (int qIdx = range.begin(); qIdx < range.end(); ++qIdx) {
            Node* qNode = queryNodes[qIdx];
            int sIdx = qNode->seqsIncluded[0];
            if (sIdx < 0 || sIdx >= static_cast<int>(database->sequences.size())) continue;
            auto* qSeqInfo = database->sequences[sIdx];
            if (!qSeqInfo || qSeqInfo->lowQuality) continue;

            const char* qRaw = qSeqInfo->alnStorage[qSeqInfo->storage];
            int qLen = qSeqInfo->len;
            std::string qSeq(qRaw, qLen);

            // A. Compute Pairwise SW matches with each backbone representative
            // and scatter into consistency table (bbLen x qLen)
            std::vector<std::vector<float>> consTable(bbLen, std::vector<float>(qLen, 0.0f));

            for (const auto& rep : backboneReps) {
                LocalHomTable::PairResult pairRes = aligner.align_affine_local_segments(
                    rep.unalignedSeq, qSeq, option->type, param
                );
                for (const auto& seg : pairRes.segments) {
                    int segLen = seg.length();
                    if (segLen <= 0) continue;
                    float resScore = seg.optScore / static_cast<float>(segLen);
                    for (int k = 0; k < segLen; ++k) {
                        int rPos = seg.start1 + k;
                        int qPos = seg.start2 + k;
                        if (rPos < static_cast<int>(rep.colMap.size()) && qPos < qLen) {
                            int bbCol = rep.colMap[rPos];
                            if (bbCol >= 0 && bbCol < bbLen) {
                                consTable[bbCol][qPos] += resScore;
                            }
                        }
                    }
                }
            }

            // Normalize consistency bonus: divide by (number of reps * backbone cardinality)
            // to align with profile-profile substitution matrix score scale
            float bbCard = std::max(1.0f, static_cast<float>(T->root->alnNum));
            float norm = 1.0f / (static_cast<float>(backboneReps.size()) * bbCard);
            for (int i = 0; i < bbLen; ++i) {
                for (int j = 0; j < qLen; ++j) {
                    consTable[i][j] *= norm;
                }
            }

            // B. Profile preparation and Gappy Columns detection
            NodePair nodePair = {T->root, qNode};
            int32_t memLen = std::max(bbLen, qLen);
            IntPair lens = {bbLen, qLen};
            std::pair<IntPairVec, IntPairVec> gappyColumns;
            stringPair consensus({"", ""});

            float* hostFreq = nullptr;
            float* hostGapOp = nullptr;
            float* hostGapEx = nullptr;
            msa::progressive::cpu::allocateMemory_and_Initialize(hostFreq, hostGapOp, hostGapEx, memLen, profileSize);

            alignment_helper::calculateProfile(hostFreq, nodePair, database, option, memLen);
            alignment_helper::getConsensus(option, hostFreq, consensus.first, bbLen);
            alignment_helper::getConsensus(option, hostFreq + profileSize * memLen, consensus.second, qLen);
            alignment_helper::removeGappyColumns(hostFreq, nodePair, option, gappyColumns, memLen, lens, database->currentTask);

            // C. Sync consistency table coordinates if gappy columns were removed from backbone
            std::vector<std::vector<float>> trimmedConsTable;
            const std::vector<std::vector<float>>* activeConsTable = &consTable;
            if (!gappyColumns.first.empty()) {
                trimmedConsTable.reserve(lens.first);
                int rawIdx = 0;
                int gcIdx = 0;
                const int numGc = static_cast<int>(gappyColumns.first.size());
                for (int i = 0; i < lens.first; ++i) {
                    while (gcIdx < numGc && rawIdx == gappyColumns.first[gcIdx].first) {
                        rawIdx += gappyColumns.first[gcIdx].second;
                        gcIdx++;
                    }
                    if (rawIdx < bbLen) {
                        trimmedConsTable.push_back(std::move(consTable[rawIdx]));
                    } else {
                        trimmedConsTable.push_back(std::vector<float>(qLen, 0.0f));
                    }
                    rawIdx++;
                }
                activeConsTable = &trimmedConsTable;
            }

            // D. Profile-Profile Alignment with consistency bonus
            msa::AlnProfile profRef, profQry;
            alignment_helper::calculatePSGP(profRef, profQry, nodePair, database, option, param, lens, gappyColumns);
            msa::progressive::cpu::freeMemory(hostFreq, hostGapOp, hostGapEx);

            std::vector<int8_t> aln_wo_gc = alignProfile(profRef, profQry, param, activeConsTable, option->consistencyWeight);

            // E. Add gappy columns back
            alnPath aln_w_gc;
            int alnRef = 0, alnQry = 0;
            for (auto a : aln_wo_gc) {
                if (a == 0) { alnRef += 1; alnQry += 1; }
                if (a == 1) { alnQry += 1; }
                if (a == 2) { alnRef += 1; }
            }
            alignment_helper::addGappyColumnsBack(aln_wo_gc, aln_w_gc, gappyColumns, param, {alnRef, alnQry}, consensus);

            // F. Store individual placement path into subtreeAln
            database->subtreeAln[sIdx] = std::move(aln_w_gc);

            completedQueries.fetch_add(1);
        }
    });

    // 4. Merge all insertions into final aligned storage
    msa::alignment_helper::mergeInsertions(database, T->root);

    int placedNum = 0;
    for (auto* s : database->sequences) {
        if (!s->lowQuality) placedNum++;
    }
    T->root->alnNum += placedNum;

    auto placeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds placeElapsed = placeEnd - placeStart;
    std::cerr << "Placed " << placedNum << " sequences with consistency bonus in " 
              << placeElapsed.count() / 1000000000 << " s\n";
}

} // namespace accurate
} // namespace msa
