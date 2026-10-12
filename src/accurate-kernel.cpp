#include "msa.hpp"
#include "align-method.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <atomic>
#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace msa {
namespace accurate {

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

                }
                totalSegments += localSegments;
            });
    }

    auto alignEnd = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(alignEnd - alignStart);
    localHomTable->pairwiseTimeMs = elapsed.count();

    return localHomTable;
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

            // Extract column provenance for reference and query sub-MSAs (hierarchical lookup with fallback)
            auto findProvenance = [&](phylogeny::Node* node) -> accurate::ColumnProvenance {
                if (!subrootProvenance) return {};
                static std::mutex provMutex;
                std::lock_guard<std::mutex> lock(provMutex);

                // 1. Direct match by node identifier
                auto it = subrootProvenance->find(node->identifier);
                if (it != subrootProvenance->end()) {
                    accurate::ColumnProvenance res = std::move(it->second);
                    subrootProvenance->erase(it);
                    return res;
                }
                // 2. Child match (if this node is an internal node whose child finished alignment)
                for (auto* ch : node->children) {
                    it = subrootProvenance->find(ch->identifier);
                    if (it != subrootProvenance->end()) {
                        accurate::ColumnProvenance res = std::move(it->second);
                        subrootProvenance->erase(it);
                        return res;
                    }
                }
                // 3. Parent match (if this node is a child continuing alignment)
                if (node->parent) {
                    it = subrootProvenance->find(node->parent->identifier);
                    if (it != subrootProvenance->end()) {
                        accurate::ColumnProvenance res = std::move(it->second);
                        subrootProvenance->erase(it);
                        return res;
                    }
                }
                return {};
            };

            if (option->accurate) {
                refProvenance = findProvenance(nodes[nIdx].first);
                qryProvenance = findProvenance(nodes[nIdx].second);

                auto getTipReps = [&](phylogeny::Node* n) -> const std::vector<accurate::Representative>* {
                    auto it = database->subtreeRepresentatives.find(n->identifier);
                    if (it != database->subtreeRepresentatives.end()) return &(it->second);
                    if (!n->seqsIncluded.empty()) {
                        auto it2 = database->subtreeRepresentatives.find(std::to_string(n->seqsIncluded[0]));
                        if (it2 != database->subtreeRepresentatives.end()) return &(it2->second);
                    }
                    return nullptr;
                };

                if (refProvenance.empty()) {
                    if (database->currentTask == 0) {
                        alignment_helper::extractColumnProvenance(nodes[nIdx].first, database, refProvenance);
                    } else {
                        const auto* reps = getTipReps(nodes[nIdx].first);
                        if (reps) {
                            extractProvenanceFromRepresentatives(*reps, refLen, refProvenance);
                        }
                    }
                }
                if (qryProvenance.empty()) {
                    if (database->currentTask == 0) {
                        alignment_helper::extractColumnProvenance(nodes[nIdx].second, database, qryProvenance);
                    } else {
                        const auto* reps = getTipReps(nodes[nIdx].second);
                        if (reps) {
                            extractProvenanceFromRepresentatives(*reps, qryLen, qryProvenance);
                        }
                    }
                }
            }

            alignment_helper::getConsensus(option, hostFreq, consensus.first, refLen);
            alignment_helper::getConsensus(option, hostFreq + profileSize * memLen, consensus.second, qryLen);
            alignment_helper::removeGappyColumns(hostFreq, nodes[nIdx], option, gappyColumns, memLen, lens, database->currentTask);

            // Keep un-trimmed provenance copies for traceback propagation
            accurate::ColumnProvenance unTrimmedRefProv = refProvenance;
            accurate::ColumnProvenance unTrimmedQryProv = qryProvenance;

            if (option->accurate) {
                // Remove gappy columns to align coordinate systems with the DP matrix
                accurate::removeColumns(refProvenance, gappyColumns.first);
                accurate::removeColumns(qryProvenance, gappyColumns.second);
            }

            // Build consistency lookup table from localHomTable
            if (localHomTable && option->accurate) {
                if (option->accurateMethod == "tcs") {
                    consistencyTable = buildConsistencyTableFromTCS(refProvenance, qryProvenance, *localHomTable);
                } else {
                    consistencyTable = buildConsistencyTableFromLocalHom(refProvenance, qryProvenance, *localHomTable);
                }
            }

            msa::AlnProfile profRef, profQry;
            alignment_helper::calculatePSGP(profRef, profQry, nodes[nIdx], database, option, param, lens, gappyColumns);

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
            (void)num;

            if (lens.first == 0) {
                for (int j = 0; j < lens.second; ++j) aln_wo_gc.push_back(1);
            } else if (lens.second == 0) {
                for (int j = 0; j < lens.first; ++j) aln_wo_gc.push_back(2);
            } else {
                int minL = std::min(lens.first, lens.second);
                int maxL = std::max(lens.first, lens.second);
                int diff = maxL - minL;

                // Profile alignment DP kernel with consistency score guidance
                const auto* cTable = (!consistencyTable.empty()) ? &consistencyTable : nullptr;
                float cWeight = option->consistencyWeight;
                aln_wo_gc = alignProfile(profRef, profQry, param, cTable, cWeight);
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

                // If tracking column provenance, update column provenance for the parent node
                if (subrootProvenance != nullptr) {
                    ColumnProvenance mergedProv(aln_w_gc.size());
                    int rCol = 0, qCol = 0;
                    for (size_t c = 0; c < aln_w_gc.size(); ++c) {
                        int8_t code = aln_w_gc[c];
                        if (code == 0) {
                            if (rCol < (int)unTrimmedRefProv.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), unTrimmedRefProv[rCol].begin(), unTrimmedRefProv[rCol].end());
                            }
                            if (qCol < (int)unTrimmedQryProv.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), unTrimmedQryProv[qCol].begin(), unTrimmedQryProv[qCol].end());
                            }
                            rCol++; qCol++;
                        } else if (code == 1) {
                            if (qCol < (int)unTrimmedQryProv.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), unTrimmedQryProv[qCol].begin(), unTrimmedQryProv[qCol].end());
                            }
                            qCol++;
                        } else if (code == 2) {
                            if (rCol < (int)unTrimmedRefProv.size()) {
                                mergedProv[c].insert(mergedProv[c].end(), unTrimmedRefProv[rCol].begin(), unTrimmedRefProv[rCol].end());
                            }
                            rCol++;
                        }
                    }
                    Node* parentNode = nodes[nIdx].first->parent ? nodes[nIdx].first->parent : nodes[nIdx].second->parent;
                    static std::mutex provMutex;
                    std::lock_guard<std::mutex> lock(provMutex);
                    std::string saveID = parentNode ? parentNode->identifier : nodes[nIdx].first->identifier;
                    (*subrootProvenance)[saveID] = std::move(mergedProv);
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
    // Maintain a local provenance map if no external subrootProvenance is provided
    std::unordered_map<std::string, ColumnProvenance> localProvenanceMap;
    std::unordered_map<std::string, ColumnProvenance>* activeProvMap = subrootProvenance ? subrootProvenance : &localProvenanceMap;

    // Wrap the CPU-based accurate alignment kernel with captured localHomTable and activeProvMap
    alnFunction accurateAlnKernel = [localHomTable, activeProvMap](Tree* tree, NodePairVec& nodes, SequenceDB* database, Option* option, Params& param) {
        alignmentKernel_Accurate_CPU(tree, nodes, database, option, param, localHomTable, activeProvMap);
    };

    // Reuse msa::progressive::msaOnSubtree with accurateAlnKernel
    int subtree = (group.layer == 1 && !group.subtreeIDs.empty()) ? group.subtreeIDs[0] : -1;
    msa::progressive::msaOnSubtree(T, database, option, param, accurateAlnKernel, subtree);
}

} // namespace accurate
} // namespace msa
