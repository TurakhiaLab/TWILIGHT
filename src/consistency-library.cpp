/*
 * =============================================================================
 * Consistency Library and Reliability Scoring for Accurate Alignment
 * -----------------------------------------------------------------------------
 * References:
 * 1. MAFFT Consistency Library (Default --accurate):
 *    Kazutaka Katoh, Daron M. Standley,
 *    MAFFT Multiple Sequence Alignment Software Version 7: Improvements in
 *    Performance and Usability,
 *    Molecular Biology and Evolution, Volume 30, Issue 4, April 2013, Pages 772-780,
 *    https://doi.org/10.1093/molbev/mst010
 *
 * 2. Transitive Consistency Score (TCS) Reliability Measure (--accurate tcs):
 *    Jia-Ming Chang, Paolo Di Tommaso, Cedric Notredame,
 *    TCS: A New Multiple Sequence Alignment Reliability Measure to Estimate
 *    Alignment Accuracy and Improve Phylogenetic Tree Reconstruction,
 *    Molecular Biology and Evolution, Volume 31, Issue 6, June 2014, Pages 1625-1637,
 *    https://doi.org/10.1093/molbev/msu117
 * =============================================================================
 */

#include "msa.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace msa {
namespace accurate {

// -----------------------------------------------------------------------------
// Transitive Consistency Score (TCS) Library Builder
// -----------------------------------------------------------------------------
std::vector<float> buildGroupTCSLibrary(
    const HierarchyGroup& group,
    std::shared_ptr<LocalHomTable> localHomTable,
    Tree* T,
    SequenceDB* database,
    Option* option,
    Params& param)
{
    auto libStart = std::chrono::high_resolution_clock::now();
    int N = localHomTable->getNumSeqs();
    if (N < 2) return std::vector<float>(N, 0.0f);

    // 1. Collect sequence lengths
    std::vector<int> seqLengths(N, 0);
    if (!group.representatives.empty()) {
        for (int i = 0; i < N; ++i) {
            seqLengths[i] = static_cast<int>(group.representatives[i].unalignedSeq.size());
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
            }
        }
    }

    // Precompute average pairwise segment score across the group as reference scale
    float sumGroupOptScore = 0.0f;
    size_t groupSegCount = 0;
    for (int u = 0; u < N; ++u) {
        for (int v = u + 1; v < N; ++v) {
            const auto& pairRes = localHomTable->getPairResult(u, v);
            for (const auto& seg : pairRes.segments) {
                if (seg.optScore > 0.0f) {
                    sumGroupOptScore += seg.optScore;
                    ++groupSegCount;
                }
            }
        }
    }
    float defaultOptScore = (groupSegCount > 0) ? (sumGroupOptScore / groupSegCount) : 50.0f;

    // 2. Build fast pairwise residue-to-residue mapping tables:
    // alignMap[z][x][k] = residue index in sequence x that residue k of sequence z aligns to (-1 if unaligned)
    // weightMap[z][x][k] = pairwise alignment confidence weight
    std::vector<std::vector<std::vector<int>>> alignMap(N, std::vector<std::vector<int>>(N));
    std::vector<std::vector<std::vector<float>>> weightMap(N, std::vector<std::vector<float>>(N));

    for (int z = 0; z < N; ++z) {
        int L_z = seqLengths[z];
        for (int x = 0; x < N; ++x) {
            alignMap[z][x].assign(L_z, -1);
            weightMap[z][x].assign(L_z, 0.0f);
        }
    }

    for (int u = 0; u < N; ++u) {
        for (int v = u + 1; v < N; ++v) {
            const auto& pairRes = localHomTable->getPairResult(u, v);
            for (const auto& seg : pairRes.segments) {
                float w = (seg.optScore > 0.0f) ? seg.optScore : defaultOptScore;
                int len = seg.end1 - seg.start1 + 1;
                for (int k = 0; k < len; ++k) {
                    int p1 = seg.start1 + k;
                    int p2 = seg.start2 + k;
                    if (p1 < seqLengths[u] && p2 < seqLengths[v]) {
                        alignMap[u][v][p1] = p2;
                        weightMap[u][v][p1] = w;
                        alignMap[v][u][p2] = p1;
                        weightMap[v][u][p2] = w;
                    }
                }
            }
        }
    }

    // 3. Allocate TCS storage in localHomTable
    localHomTable->initTCS();

    // 4. Parallel Transitivity Computation across all pairs (x, y) with x < y
    size_t totalPairs = (size_t)N * (N - 1) / 2;
    std::vector<std::pair<int, int>> pairJobs;
    pairJobs.reserve(totalPairs);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            pairJobs.push_back({i, j});
        }
    }

    std::atomic<size_t> totalTcsEntries{0};

    tbb::parallel_for(tbb::blocked_range<size_t>(0, totalPairs),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t p = range.begin(); p < range.end(); ++p) {
                int x = pairJobs[p].first;
                int y = pairJobs[p].second;
                int Lx = seqLengths[x];
                int Ly = seqLengths[y];
                if (Lx == 0 || Ly == 0) continue;

                // Determine pairOptScore for DP scaling to match scoring matrix
                float pairOptScore = defaultOptScore;
                const auto& directPair = localHomTable->getPairResult(x, y);
                if (!directPair.segments.empty()) {
                    float sumOpt = 0.0f;
                    for (const auto& s : directPair.segments) sumOpt += s.optScore;
                    pairOptScore = sumOpt / directPair.segments.size();
                }

                // Temporary accumulators for pair (x, y)
                std::unordered_map<uint64_t, float> numMap;
                std::vector<float> denomX(Lx, 0.0f);
                std::vector<float> denomY(Ly, 0.0f);

                // Direct pairwise alignment contribution
                for (int i = 0; i < Lx; ++i) {
                    int j = alignMap[x][y][i];
                    if (j >= 0) {
                        float w = weightMap[x][y][i];
                        uint64_t key = (static_cast<uint64_t>(i) << 32) | static_cast<uint32_t>(j);
                        numMap[key] += w;
                        denomX[i] += w;
                        denomY[j] += w;
                    }
                }

                // Transitive alignments through all intermediate sequences z != x, y
                for (int z = 0; z < N; ++z) {
                    if (z == x || z == y) continue;
                    int Lz = seqLengths[z];
                    if (Lz == 0) continue;
                    const auto& mapZX = alignMap[z][x];
                    const auto& mapZY = alignMap[z][y];
                    const auto& wtZX = weightMap[z][x];
                    const auto& wtZY = weightMap[z][y];

                    for (int k = 0; k < Lz; ++k) {
                        int i = mapZX[k];
                        int j = mapZY[k];
                        if (i >= 0 && j >= 0) {
                            float w = std::min(wtZX[k], wtZY[k]);
                            uint64_t key = (static_cast<uint64_t>(i) << 32) | static_cast<uint32_t>(j);
                            numMap[key] += w;
                            denomX[i] += w;
                            denomY[j] += w;
                        } else if (i >= 0) {
                            denomX[i] += wtZX[k];
                        } else if (j >= 0) {
                            denomY[j] += wtZY[k];
                        }
                    }
                }

                // Compute PairTCS using Formula (1):
                // PairTCS(R^x_i, R^y_j) = 2 * sum_z min(W_xz, W_zy) / (sum_z min(W_xz, W_z*) + sum_z min(W_*z, W_zy))
                LocalHomTable::TCSPairResult tcsRes;
                tcsRes.entries.reserve(numMap.size());
                for (const auto& kv : numMap) {
                    int32_t i = static_cast<int32_t>(kv.first >> 32);
                    int32_t j = static_cast<int32_t>(kv.first & 0xFFFFFFFF);
                    float numVal = kv.second;
                    float denomVal = denomX[i] + denomY[j];
                    if (denomVal > 1e-6f) {
                        float pairTCS = (2.0f * numVal) / denomVal;
                        pairTCS = std::min(1.0f, std::max(0.0f, pairTCS));
                        // Filter weak alignment noise below 0.10
                        if (pairTCS >= 0.10f) {
                            float scaledScore = pairTCS * pairOptScore;
                            tcsRes.entries.push_back({i, j, scaledScore, pairTCS});
                        }
                    }
                }

                totalTcsEntries += tcsRes.entries.size();
                localHomTable->setTCSPairResult(x, y, std::move(tcsRes));
            }
        });

    // 5. Compute ResidueTCS and Sequence centrality scores (Formulas 3 & 4)
    // ResidueTCS(x_i) = \frac{1}{|S| - 1} \sum_{y \neq x} PairTCS(x_i, y_j)
    std::vector<std::vector<float>> residueTCS(N);
    for (int i = 0; i < N; ++i) {
        residueTCS[i].assign(seqLengths[i], 0.0f);
    }

    float normDenom = (N > 1) ? (1.0f / (N - 1)) : 1.0f;
    for (int x = 0; x < N; ++x) {
        for (int y = x + 1; y < N; ++y) {
            const auto& tcsRes = localHomTable->getTCSPairResult(x, y);
            for (const auto& entry : tcsRes.entries) {
                if (entry.pos1 < seqLengths[x]) {
                    residueTCS[x][entry.pos1] += entry.tcs * normDenom;
                }
                if (entry.pos2 < seqLengths[y]) {
                    residueTCS[y][entry.pos2] += entry.tcs * normDenom;
                }
            }
        }
    }

    std::vector<int> sortedLens = seqLengths;
    std::sort(sortedLens.begin(), sortedLens.end());
    float medianLen = (N > 0) ? static_cast<float>(sortedLens[N / 2]) : 1.0f;
    if (medianLen <= 0.0f) medianLen = 1.0f;

    std::vector<float> centralityScores(N, 0.0f);
    double totalResTCS = 0.0;
    size_t totalPositions = 0;
    for (int i = 0; i < N; ++i) {
        if (seqLengths[i] == 0) continue;
        double sumResTCS = 0.0;
        for (float v : residueTCS[i]) {
            sumResTCS += v;
        }
        totalResTCS += sumResTCS;
        totalPositions += seqLengths[i];
        float meanResTCS = static_cast<float>(sumResTCS / seqLengths[i]);
        float lengthFactor = std::min(1.0f, static_cast<float>(seqLengths[i]) / medianLen);
        centralityScores[i] = meanResTCS * lengthFactor;
    }
    float avgTCS = (totalPositions > 0) ? static_cast<float>(totalResTCS / totalPositions) : 0.0f;

    // Store residue reliability scores into LocalHomTable
    localHomTable->setPositionImportance(std::move(residueTCS));

    auto libEnd = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(libEnd - libStart);
    localHomTable->libTimeMs = elapsed.count();
    localHomTable->meanImportance = avgTCS;
    long long totalTimeMs = localHomTable->pairwiseTimeMs + localHomTable->libTimeMs;

    if (option->printDetail) {
        std::cerr << "[Consistency Group " << group.groupID << "] Seqs=" << N 
                  << " | Time=" << totalTimeMs << " ms"
                  << " | Mean Positional Importance=" << std::fixed << std::setprecision(2) << avgTCS 
                  << "\n";
    }

    return centralityScores;
}

// -----------------------------------------------------------------------------
// 2. Build consistency score library & calculate centrality scores
// -----------------------------------------------------------------------------
std::vector<float> buildGroupConsistencyLibrary(const HierarchyGroup& group, std::shared_ptr<LocalHomTable> localHomTable, Tree* T, SequenceDB* database, Option* option, Params& param) {
    if (option->accurateMethod == "tcs") {
        return buildGroupTCSLibrary(group, localHomTable, T, database, option, param);
    }

    auto libStart = std::chrono::high_resolution_clock::now();
    int N = localHomTable->getNumSeqs();
    if (N < 2) return std::vector<float>(N, 0.0f);

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
    localHomTable->libTimeMs = elapsed.count();
    localHomTable->meanImportance = avgImp;
    long long totalTimeMs = localHomTable->pairwiseTimeMs + localHomTable->libTimeMs;

    if (option->printDetail) {
        std::cerr << "[Consistency Group " << group.groupID << "] Seqs=" << N 
                  << " | Time=" << totalTimeMs << " ms"
                  << " | Mean Positional Importance=" << std::fixed << std::setprecision(2) << avgImp 
                  << "\n";
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
// Consistency Table Builder: Transitive Consistency Score (TCS)
// Maps precomputed PairTCS cross-profile connections into Profile DP matrix
// -----------------------------------------------------------------------------
std::vector<std::vector<float>> buildConsistencyTableFromTCS(
    const ColumnProvenance& refProvenance,
    const ColumnProvenance& qryProvenance,
    const LocalHomTable& localHomTable,
    float scaleFactor)
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

    for (const auto& refEntry : refPosToCol) {
        int localA = refEntry.first;
        const auto& refMap = refEntry.second;

        for (const auto& qryEntry : qryPosToCol) {
            int localB = qryEntry.first;
            if (localA == localB) continue;

            const auto& qryMap = qryEntry.second;
            const auto& tcsRes = localHomTable.getTCSPairResult(localA, localB);
            if (tcsRes.entries.empty()) continue;

            if (localA < localB) {
                for (const auto& entry : tcsRes.entries) {
                    int posA = entry.pos1;
                    int posB = entry.pos2;
                    if (posA < static_cast<int>(refMap.size()) && posB < static_cast<int>(qryMap.size())) {
                        int cA = refMap[posA];
                        int cB = qryMap[posB];
                        if (cA >= 0 && cB >= 0) {
                            table[cA][cB] += entry.score;
                        }
                    }
                }
            } else {
                for (const auto& entry : tcsRes.entries) {
                    int posB = entry.pos1;
                    int posA = entry.pos2;
                    if (posA < static_cast<int>(refMap.size()) && posB < static_cast<int>(qryMap.size())) {
                        int cA = refMap[posA];
                        int cB = qryMap[posB];
                        if (cA >= 0 && cB >= 0) {
                            table[cA][cB] += entry.score;
                        }
                    }
                }
            }
        }
    }

    // Normalize cross-cluster consistency and apply scale factor
    // Resulting score represents the average PairTCS of cross-profile pairs, scaled for DP
    float norm = scaleFactor / (static_cast<float>(refPosToCol.size()) * static_cast<float>(qryPosToCol.size()));
    for (int i = 0; i < refLen; ++i) {
        for (int j = 0; j < qryLen; ++j) {
            table[i][j] *= norm;
        }
    }

    return table;
}

// -----------------------------------------------------------------------------
// ColumnTCS: Formula (2) from TCS Paper (Chang et al., MBE 2014)
// Computes reliability score for an alignment column across all aligned residues
// -----------------------------------------------------------------------------
float computeColumnTCS(
    const std::vector<std::pair<int, int>>& colResidues,
    const LocalHomTable& localHomTable)
{
    size_t k = colResidues.size();
    if (k < 2) return 0.0f;

    double sumPairTCS = 0.0;
    size_t pairCount = 0;

    for (size_t x = 0; x < k; ++x) {
        int seqX = colResidues[x].first;
        int posX = colResidues[x].second;
        for (size_t y = x + 1; y < k; ++y) {
            int seqY = colResidues[y].first;
            int posY = colResidues[y].second;
            if (seqX == seqY) continue;

            const auto& tcsRes = localHomTable.getTCSPairResult(seqX, seqY);
            float score = 0.0f;
            if (seqX < seqY) {
                for (const auto& entry : tcsRes.entries) {
                    if (entry.pos1 == posX && entry.pos2 == posY) {
                        score = entry.tcs;
                        break;
                    }
                }
            } else {
                for (const auto& entry : tcsRes.entries) {
                    if (entry.pos1 == posY && entry.pos2 == posX) {
                        score = entry.tcs;
                        break;
                    }
                }
            }
            sumPairTCS += score;
            ++pairCount;
        }
    }

    return (pairCount > 0) ? static_cast<float>(sumPairTCS / pairCount) : 0.0f;
}

} // namespace accurate
} // namespace msa
