#include "msa.hpp"
#include "kseq.h"

#include <zlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>

KSEQ_INIT2(static, gzFile, gzread);

namespace msa {
namespace accurate {

// Collect all leaves under a node in topological DFS order
void collectLeavesInOrder(Node* node, std::vector<std::string>& leaves) {
    if (!node) return;
    if (node->is_leaf() || node->children.empty()) {
        leaves.push_back(node->identifier);
        return;
    }
    for (auto* child : node->children) {
        collectLeavesInOrder(child, leaves);
    }
}

// Extract representative sequence objects using centrality scores across phylogenetic clades
std::vector<Representative> extractRepresentatives(
    Node* root,
    size_t targetReps,
    SequenceDB* database,
    int groupID,
    const std::vector<float>& centralityScores,
    const LocalHomTable& localHomTable,
    int alnLen,
    Option* option,
    std::string* logMsg,
    int subtreeID
) {
    std::vector<std::string> leaves;
    collectLeavesInOrder(root, leaves);

    std::vector<Representative> reps;
    if (leaves.empty() || targetReps == 0) return reps;

    float sampleTop = (option != nullptr) ? option->accSampleTop : 0.3f;

    struct Candidate {
        size_t leafIdx;
        std::string name;
        int globalId;
        float centrality;
        std::vector<int> coveredCols;
    };

    std::vector<Candidate> candidates;
    candidates.reserve(leaves.size());
    std::vector<int> colCoverage(alnLen > 0 ? alnLen : 0, 0);

    for (size_t idx = 0; idx < leaves.size(); ++idx) {
        const auto& name = leaves[idx];
        auto it = database->name_map.find(name);
        if (it == database->name_map.end()) continue;

        const auto* sInfo = it->second;
        int globalId = sInfo->id;
        int locIdx = localHomTable.getLocalIndex(globalId);
        float score = (locIdx >= 0 && locIdx < static_cast<int>(centralityScores.size()))
                      ? centralityScores[locIdx] : 0.0f;

        Candidate cand;
        cand.leafIdx = idx;
        cand.name = name;
        cand.globalId = globalId;
        cand.centrality = score;

        const char* rawChars = sInfo->alnStorage[sInfo->storage];
        for (int k = 0; k < sInfo->len && k < alnLen; ++k) {
            if (rawChars[k] != '-' && rawChars[k] != '.' && rawChars[k] != '\0') {
                cand.coveredCols.push_back(k);
                if (k < static_cast<int>(colCoverage.size())) {
                    colCoverage[k]++;
                }
            }
        }
        candidates.push_back(std::move(cand));
    }

    int totalCoverableCols = 0;
    for (int k = 0; k < alnLen; ++k) {
        if (colCoverage[k] > 0) ++totalCoverableCols;
    }

    std::vector<std::string> selectedNames;
    if (candidates.size() <= targetReps) {
        for (const auto& c : candidates) {
            selectedNames.push_back(c.name);
        }
    } else {
        // Sort candidate indices by centrality score descending
        std::vector<size_t> candOrder(candidates.size());
        std::iota(candOrder.begin(), candOrder.end(), 0);
        std::sort(candOrder.begin(), candOrder.end(), [&](size_t a, size_t b) {
            return candidates[a].centrality > candidates[b].centrality;
        });

        std::vector<bool> isCandidateSelected(candidates.size(), false);
        std::vector<bool> isColCovered(alnLen > 0 ? alnLen : 0, false);
        std::vector<size_t> selectedCandIndices;
        selectedCandIndices.reserve(targetReps);

        // Step 1: Select top T% highest centrality sequences as core anchors
        size_t coreCount = std::max((size_t)1, static_cast<size_t>(std::round(sampleTop * targetReps)));
        if (coreCount > targetReps) coreCount = targetReps;
        if (coreCount > candidates.size()) coreCount = candidates.size();

        for (size_t i = 0; i < coreCount; ++i) {
            size_t candIdx = candOrder[i];
            isCandidateSelected[candIdx] = true;
            selectedCandIndices.push_back(candIdx);
            for (int col : candidates[candIdx].coveredCols) {
                if (col >= 0 && col < static_cast<int>(isColCovered.size())) {
                    isColCovered[col] = true;
                }
            }
        }

        // Step 2: Feature-driven expansion: sort unmarked columns by coverage descending
        if (selectedCandIndices.size() < targetReps && alnLen > 0) {
            std::vector<int> unmarkedCols;
            unmarkedCols.reserve(alnLen);
            for (int col = 0; col < alnLen; ++col) {
                if (!isColCovered[col] && colCoverage[col] > 0) {
                    unmarkedCols.push_back(col);
                }
            }

            std::sort(unmarkedCols.begin(), unmarkedCols.end(), [&](int c1, int c2) {
                if (colCoverage[c1] != colCoverage[c2]) {
                    return colCoverage[c1] > colCoverage[c2];
                }
                return c1 < c2;
            });

            // Inverted index: column -> candidate indices that cover this column (ordered by centrality desc)
            std::vector<std::vector<size_t>> colToCandidates(alnLen);
            for (size_t candIdx : candOrder) {
                for (int col : candidates[candIdx].coveredCols) {
                    if (col >= 0 && col < alnLen) {
                        colToCandidates[col].push_back(candIdx);
                    }
                }
            }

            for (int col : unmarkedCols) {
                if (selectedCandIndices.size() >= targetReps) break;
                if (isColCovered[col]) continue; // Covered by an earlier selected candidate

                // Find first unselected candidate covering this column (highest centrality)
                for (size_t candIdx : colToCandidates[col]) {
                    if (!isCandidateSelected[candIdx]) {
                        isCandidateSelected[candIdx] = true;
                        selectedCandIndices.push_back(candIdx);
                        for (int c : candidates[candIdx].coveredCols) {
                            if (c >= 0 && c < static_cast<int>(isColCovered.size())) {
                                isColCovered[c] = true;
                            }
                        }
                        break;
                    }
                }
            }
        }

        // Step 3: Backfill remaining quota with highest centrality candidates
        if (selectedCandIndices.size() < targetReps) {
            for (size_t candIdx : candOrder) {
                if (!isCandidateSelected[candIdx]) {
                    isCandidateSelected[candIdx] = true;
                    selectedCandIndices.push_back(candIdx);
                    if (selectedCandIndices.size() == targetReps) break;
                }
            }
        }

        // Sort selected candidates by original tree leaf index to preserve phylogenetic ordering
        std::sort(selectedCandIndices.begin(), selectedCandIndices.end(), [&](size_t a, size_t b) {
            return candidates[a].leafIdx < candidates[b].leafIdx;
        });

        for (size_t cIdx : selectedCandIndices) {
            selectedNames.push_back(candidates[cIdx].name);
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
            rep.subtreeID = (subtreeID >= 0) ? subtreeID : groupID;
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

    // Compute post-sampling column coverage and mean centrality logging
    if (alnLen > 0 && option && option->printDetail) {
        std::vector<bool> repColCovered(alnLen, false);
        double sumRepCentrality = 0.0;
        for (const auto& rep : reps) {
            sumRepCentrality += rep.centralityScore;
            for (int col : rep.colInSubAln) {
                if (col >= 0 && col < alnLen) {
                    repColCovered[col] = true;
                }
            }
        }
        int sampledCoveredCols = 0;
        for (int c = 0; c < alnLen; ++c) {
            if (repColCovered[c]) ++sampledCoveredCols;
        }
        float avgRepCentrality = reps.empty() ? 0.0f : static_cast<float>(sumRepCentrality / reps.size());
        float covPct = (totalCoverableCols > 0) ? (100.0f * sampledCoveredCols / totalCoverableCols) : 100.0f;

        std::ostringstream oss;
        oss << "[Sampling Group " << groupID << "] Selected=" << reps.size()
            << "/" << leaves.size() << " seqs | Column Coverage=" << std::fixed << std::setprecision(1) << covPct
            << "% | Mean Rep Centrality=" << std::fixed << std::setprecision(2) << avgRepCentrality << "\n";

        if (logMsg) {
            *logMsg = oss.str();
        } else {
            std::cerr << oss.str();
        }
    }

    return reps;
}

void extractProvenanceFromRepresentatives(
    const std::vector<accurate::Representative>& reps,
    int alnLen,
    accurate::ColumnProvenance& provenance)
{
    provenance.assign(alnLen, accurate::ColumnResidues{});
    for (const auto& rep : reps) {
        for (size_t k = 0; k < rep.unalignedSeq.size() && k < rep.colInSubAln.size(); ++k) {
            int col = rep.colInSubAln[k];
            if (col >= 0 && col < alnLen) {
                provenance[col].push_back({rep.globalSeqID, static_cast<int>(k), 1.0f});
            }
        }
    }
}

void removeColumns(ColumnProvenance& provenance, const IntPairVec& removedColumns)
{
    if (removedColumns.empty()) return;
    ColumnProvenance filtered;
    filtered.reserve(provenance.size());
    int removeIdx = 0;
    int nextRemoveStart = removedColumns[removeIdx].first;
    int nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
    for (int col = 0; col < static_cast<int>(provenance.size()); ++col) {
        while (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveEnd) {
            ++removeIdx;
            if (removeIdx < static_cast<int>(removedColumns.size())) {
                nextRemoveStart = removedColumns[removeIdx].first;
                nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
            }
        }
        const bool removed = (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveStart && col < nextRemoveEnd);
        if (!removed) filtered.push_back(std::move(provenance[col]));
    }
    provenance = std::move(filtered);
}

// -----------------------------------------------------------------------------
// Hybrid Representative Selection for Backbone Alignment (Placement mode)
// Strategy: Select N/2 longest sequences (least gaps), then N/2 greedy column coverage
// -----------------------------------------------------------------------------
std::vector<BackboneRep> extractBackboneRepresentatives(
    const std::string& backboneFile,
    int targetReps,
    Option* option)
{
    gzFile f_rd = gzopen(backboneFile.c_str(), "r");
    if (!f_rd) {
        std::cerr << "ERROR: Cannot open backbone file: " << backboneFile << "\n";
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);

    struct RawSeqEntry {
        std::string name;
        std::string alignedSeq;
        int validLength = 0;
    };

    std::vector<RawSeqEntry> rawSeqs;
    int msaLen = 0;
    while (kseq_read(kseq_rd) >= 0) {
        RawSeqEntry entry;
        entry.name = kseq_rd->name.s;
        entry.alignedSeq = std::string(kseq_rd->seq.s, kseq_rd->seq.l);
        if (rawSeqs.empty()) {
            msaLen = kseq_rd->seq.l;
        }
        int vLen = 0;
        for (char ch : entry.alignedSeq) {
            if (ch != '-' && ch != '.') ++vLen;
        }
        entry.validLength = vLen;
        rawSeqs.push_back(std::move(entry));
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    int totalSeqs = static_cast<int>(rawSeqs.size());
    if (totalSeqs == 0 || msaLen == 0) return {};

    targetReps = std::min(targetReps, totalSeqs);
    int topQuota = std::max(1, targetReps / 2);
    int colQuota = targetReps - topQuota;

    std::vector<int> sortedIndices(totalSeqs);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&](int a, int b) {
        return rawSeqs[a].validLength > rawSeqs[b].validLength; // Descending
    });

    std::vector<int> chosenIndices;
    chosenIndices.reserve(targetReps);
    std::vector<bool> isChosen(totalSeqs, false);
    std::vector<bool> coveredCols(msaLen, false);

    // 1. Pick top N/2 longest sequences
    for (int i = 0; i < topQuota; ++i) {
        int idx = sortedIndices[i];
        chosenIndices.push_back(idx);
        isChosen[idx] = true;
        for (int c = 0; c < msaLen; ++c) {
            char ch = rawSeqs[idx].alignedSeq[c];
            if (ch != '-' && ch != '.') coveredCols[c] = true;
        }
    }

    // 2. Pick remaining N/2 sequences by maximizing greedy column coverage
    for (int step = 0; step < colQuota; ++step) {
        int bestIdx = -1;
        int bestGain = -1;
        for (int i = 0; i < totalSeqs; ++i) {
            if (isChosen[i]) continue;
            int gain = 0;
            const auto& seqStr = rawSeqs[i].alignedSeq;
            for (int c = 0; c < msaLen; ++c) {
                if (!coveredCols[c] && seqStr[c] != '-' && seqStr[c] != '.') {
                    gain++;
                }
            }
            if (gain > bestGain) {
                bestGain = gain;
                bestIdx = i;
            }
        }
        if (bestIdx == -1 || bestGain == 0) {
            // No more columns to cover, fill remaining with next longest
            for (int i = 0; i < totalSeqs && static_cast<int>(chosenIndices.size()) < targetReps; ++i) {
                int cand = sortedIndices[i];
                if (!isChosen[cand]) {
                    chosenIndices.push_back(cand);
                    isChosen[cand] = true;
                }
            }
            break;
        }
        chosenIndices.push_back(bestIdx);
        isChosen[bestIdx] = true;
        const auto& seqStr = rawSeqs[bestIdx].alignedSeq;
        for (int c = 0; c < msaLen; ++c) {
            if (seqStr[c] != '-' && seqStr[c] != '.') coveredCols[c] = true;
        }
    }

    // 3. Construct BackboneRep objects with unalignedSeq and colMap
    std::vector<BackboneRep> reps;
    reps.reserve(chosenIndices.size());
    for (size_t r = 0; r < chosenIndices.size(); ++r) {
        int idx = chosenIndices[r];
        BackboneRep rep;
        rep.id = static_cast<int>(r);
        rep.name = rawSeqs[idx].name;
        rep.unalignedSeq.reserve(rawSeqs[idx].validLength);
        rep.colMap.reserve(rawSeqs[idx].validLength);

        for (int c = 0; c < msaLen; ++c) {
            char ch = rawSeqs[idx].alignedSeq[c];
            if (ch != '-' && ch != '.') {
                rep.unalignedSeq.push_back(ch);
                rep.colMap.push_back(c);
            }
        }
        reps.push_back(std::move(rep));
    }

    int totalCovered = 0;
    for (bool b : coveredCols) if (b) ++totalCovered;

    std::cerr << "  [Accurate Placement] Selected " << reps.size() << "/" << totalSeqs 
              << " backbone representatives (" << topQuota << " longest, " 
              << (reps.size() - topQuota) << " max coverage): Column Coverage = " 
              << totalCovered << "/" << msaLen << " (" 
              << std::fixed << std::setprecision(1) << (100.0 * totalCovered / msaLen) << "%)\n";

    return reps;
}

} // namespace accurate
} // namespace msa
