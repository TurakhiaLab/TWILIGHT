#include "merge_debug.hpp"
#include "message_manager.hpp"
#include "block.hpp"
#include "block_set.hpp"
#include "block_manager.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>

namespace mga {

MergeDebugger::MergeDebugger(bool enabled) {
    if (enabled) {
        setPreset(MergeDebugPreset::STANDARD);
    } else {
        setPreset(MergeDebugPreset::OFF);
    }
}

MergeDebugger::MergeDebugger(MergeDebugPreset preset) {
    setPreset(preset);
}

MergeDebugger::MergeDebugger(const MergeDebugConfig &config)
    : config_(config) {}

void MergeDebugger::setPreset(MergeDebugPreset preset) {
    switch (preset) {
        case MergeDebugPreset::OFF:
            config_.enabled = false;
            break;
        case MergeDebugPreset::SUMMARY:
            config_.enabled = true;
            config_.phases = true;
            config_.aln_details = false;
            config_.skips = false;
            config_.energy_check = false;
            config_.core_extract = false;
            config_.crossing_check = false;
            config_.scenarios = false;
            config_.merge_results = false;
            config_.case_d_dump = false;
            config_.validations = false;
            break;
        case MergeDebugPreset::STANDARD:
            config_.enabled = true;
            config_.phases = true;
            config_.aln_details = false;
            config_.skips = true;
            config_.energy_check = true;
            config_.core_extract = false;
            config_.crossing_check = false;
            config_.scenarios = true;
            config_.merge_results = true;
            config_.case_d_dump = false;
            config_.validations = false;
            break;
        case MergeDebugPreset::DETAILED:
            config_.enabled = true;
            config_.phases = true;
            config_.aln_details = true;
            config_.skips = true;
            config_.energy_check = true;
            config_.core_extract = true;
            config_.crossing_check = true;
            config_.scenarios = true;
            config_.merge_results = true;
            config_.case_d_dump = false;
            config_.validations = false;
            break;
        case MergeDebugPreset::FULL:
            config_.enabled = true;
            config_.phases = true;
            config_.aln_details = true;
            config_.skips = true;
            config_.energy_check = true;
            config_.core_extract = true;
            config_.crossing_check = true;
            config_.scenarios = true;
            config_.merge_results = true;
            config_.case_d_dump = true;
            config_.validations = true;
            break;
    }
}

void MergeDebugger::logMergeHeader(const BlockSetID &refId, const BlockSetID &qryId) {
    if (!config_.enabled || !config_.phases) return;
    std::ostringstream oss;
    oss << "\n========================================================\n"
        << "=== TWILIGHT-MGA DYNAMIC MERGE: " << refId << " + " << qryId << " ===\n"
        << "========================================================";
    mga::msg::debug("%s", oss.str().c_str());
}

void MergeDebugger::logPhase(int phaseNumber, const std::string &description) {
    if (!config_.enabled || !config_.phases) return;
    mga::msg::debug("[Phase %d] %s", phaseNumber, description.c_str());
}

void MergeDebugger::logProcessedAlignments(int count) {
    if (!config_.enabled || !config_.phases) return;
    mga::msg::debug("\n  -> Processed %d valid alignments.", count);
}

void MergeDebugger::logMergeSummary(std::chrono::high_resolution_clock::time_point startTime) {
    if (!config_.enabled || !config_.phases) return;
    auto timeEnd = std::chrono::high_resolution_clock::now();
    long long durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - startTime).count();
    std::ostringstream oss;
    oss << "\n========================================================\n"
        << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
        << "Total Execution Time: " << durationMs << " ms\n"
        << "========================================================\n";
    mga::msg::debug("%s", oss.str().c_str());
}

void MergeDebugger::logAlignmentDetails(const CigarString &cigar, bool isInverse,
                                       int r_start, int r_end, int q_start, int q_end,
                                       bool r_merged, bool q_merged,
                                       const std::vector<BlockID> &r_overlaps,
                                       const std::vector<BlockID> &q_overlaps) {
    if (!config_.enabled || !config_.aln_details) return;

    auto formatOverlaps = [](const std::vector<BlockID> &list) -> std::string {
        std::ostringstream oss;
        oss << "[ ";
        for (size_t i = 0; i < list.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << list[i];
        }
        oss << " ]";
        return oss.str();
    };

    std::ostringstream oss;
    oss << "\n========================================\n"
        << "🔍 [ALIGNMENT DETAILS]\n"
        << "  ├─ CIGAR : " << cigarToStr(cigar) << (isInverse ? "(-)" : "(+)") << "\n"
        << "  ├─ Ref   : [" << r_start << ", " << r_end << ")\n"
        << "  ├─ Qry   : [" << q_start << ", " << q_end << ")";

    if (r_merged) {
        oss << "\n  ├─ Ref Overlap : " << formatOverlaps(r_overlaps);
    }
    if (q_merged) {
        oss << "\n  └─ Qry Overlap : " << formatOverlaps(q_overlaps);
    }
    oss << "\n========================================";

    mga::msg::debug("%s", oss.str().c_str());
}

void MergeDebugger::logSkipIdenticalBlock(BlockID blkId) {
    if (!config_.enabled || !config_.skips) return;
    mga::msg::debug("  -> SKIP. Both paths resolve to the exact same block %d.", blkId);
}

void MergeDebugger::logNullAlignment() {
    if (!config_.enabled || !config_.skips) return;
    mga::msg::debug("  [Null Alignment] Skip.");
}

void MergeDebugger::logPureIndelBypass() {
    if (!config_.enabled || !config_.skips) return;
    mga::msg::debug("  -> ⚡ [PURE INDEL] Bypassing logic. Block is already safely in mergedSet.");
}

void MergeDebugger::logSecondaryScoreSkip(int score, int threshold) {
    if (!config_.enabled || !config_.skips) return;
    mga::msg::debug("  -> SKIP. Secondary alignment score below threshold (%d < %d).", score, threshold);
}

void MergeDebugger::logOccurrenceCrossingReject(bool isExtractCore) {
    if (!config_.enabled || !config_.energy_check) return;
    if (isExtractCore) {
        mga::msg::debug("  🚨 [REJECT] extractCoreBlock: Alignment crosses occurrence boundaries of the same block. Rejecting merge.");
    } else {
        mga::msg::debug("      ❌ [ALIGNMENT REJECTED] Crosses occurrence boundary of the same block.");
    }
}

void MergeDebugger::logEnergyCalculation(bool isRef, int length, int Nc, int energy) {
    if (!config_.enabled || !config_.energy_check) return;
    mga::msg::debug("      📊 [ENERGY-CHECK] %s | L: %d bp, Nc: %d -> Energy: %d",
                    isRef ? "Ref" : "Qry", length, Nc, energy);
}

void MergeDebugger::logEnergyReject() {
    if (!config_.enabled || !config_.energy_check) return;
    mga::msg::debug("      ❌ [ALIGNMENT REJECTED] High topology energy. Fragmentation prevented.");
}

void MergeDebugger::logMultiBlockSpanningReject(size_t r_count, size_t q_count) {
    if (!config_.enabled || !config_.energy_check) return;
    if (q_count == 0) {
        mga::msg::debug("  🚨 [REJECT] extractCoreBlock: Alignment spans %zu blocks. Rejecting merge.", r_count);
    } else {
        mga::msg::debug("  🚨 [ALIGNMENT REJECTED] Alignment spans multiple blocks (Ref: %zu, Qry: %zu). Bypassing merge.",
                        r_count, q_count);
    }
}

void MergeDebugger::logTopologySkip() {
    if (!config_.enabled || !config_.energy_check) return;
    mga::msg::debug("  🚨 [ALIGNMENT REJECTED] Skipping this alignment to maintain graph topology.");
}

void MergeDebugger::logLowIdentityReject(double identity) {
    if (!config_.enabled || !config_.energy_check) return;
    mga::msg::debug("  🚨 [REJECT] Alignment Identity too low (%.2f%%). Graph untouched. Skipping.", identity * 100.0);
}

void MergeDebugger::logInvalidCoreBlocks(BlockID rCore, BlockID qCore) {
    if (!config_.enabled || !config_.energy_check) return;
    mga::msg::debug("  🚨 [SKIP] Invalid or identical core blocks (rCore: %d, qCore: %d). Bypassing merge.", rCore, qCore);
}

void MergeDebugger::logCoreExtraction(const std::string &trackName, BlockID coreId,
                                     const CoordinateManager::GlobalCoords &coords) {
    if (!config_.enabled || !config_.core_extract) return;
    std::ostringstream oss;
    oss << "  [" << trackName << "] Extracting core block " << coreId << "...\n";

    auto printIntervals = [&](const std::string &label, const auto &list) {
        oss << "         ├─ " << label << ": ";
        if (list.empty()) {
            oss << "None";
        } else {
            for (const auto &iv : list) {
                oss << "[" << iv.start << ", " << iv.end << ") " << (iv.isInverse ? "(-)" : "(+)") << " ";
            }
        }
        oss << "\n";
    };

    printIntervals("Ref", coords.refIntervals);
    printIntervals("Qry", coords.qryIntervals);

    mga::msg::debug("%s", oss.str().c_str());
}

void MergeDebugger::logAnchorNeighbors(const std::pair<int, int> &ref_before,
                                      const std::pair<int, int> &ref_after,
                                      const std::pair<int, int> &qry_before,
                                      const std::pair<int, int> &qry_after) {
    if (!config_.enabled || !config_.crossing_check) return;
    mga::msg::debug("Ref: (%d,%d)\t(%d,%d)\nQry: (%d,%d)\t(%d,%d)",
                    ref_before.first, ref_before.second, ref_after.first, ref_after.second,
                    qry_before.first, qry_before.second, qry_after.first, qry_after.second);
}

void MergeDebugger::logScenario(const std::string &scenarioName, const std::string &details) {
    if (!config_.enabled || !config_.scenarios) return;
    if (details.empty()) {
        mga::msg::debug("      [%s]", scenarioName.c_str());
    } else {
        mga::msg::debug("  [%s] %s", scenarioName.c_str(), details.c_str());
    }
}

void MergeDebugger::logCaseAPatchedCigar(const CigarString &cigar) {
    if (!config_.enabled || !config_.scenarios) return;
    std::string cigarStr = cigarToStr(cigar);
    mga::msg::debug("      -> Patched CIGAR for Case A: %s", cigarStr.c_str());
}

void MergeDebugger::logRealigning(size_t rLen, size_t qLen) {
    if (!config_.enabled || !config_.scenarios) return;
    mga::msg::debug("      -> Re-aligning Consensus: Ref(%zubp) vs Qry(%zubp)...", rLen, qLen);
}

void MergeDebugger::logRealignedCigar(const CigarString &cigar, double identity) {
    if (!config_.enabled || !config_.scenarios) return;
    std::string cigarStr = cigarToStr(cigar);
    mga::msg::debug("      -> Re-aligned CIGAR: %s\n      -> Alignment Identity: %.2f%%",
                    cigarStr.c_str(), identity * 100.0);
}

void MergeDebugger::logCaseDConsensus(const std::string &rCons, const std::string &qCons) {
    if (!config_.enabled || !config_.scenarios) return;
    mga::msg::debug("%s\n%s", rCons.c_str(), qCons.c_str());
}

void MergeDebugger::logMergedBlock(BlockID mId, size_t length) {
    if (!config_.enabled || !config_.merge_results) return;
    mga::msg::debug("      -> Block %d created/merged (Length: %zu)", mId, length);
}

void MergeDebugger::logCaseDDetails(Block *mBlk) {
    if (!config_.enabled || !config_.case_d_dump || !mBlk) return;

    std::ostringstream oss;
    oss << "================================================================================\n"
        << "  [Case D] Merged Block ID: " << mBlk->getId()
        << " (Consensus Length: " << mBlk->getConsensus().length() << " bp)\n"
        << "--------------------------------------------------------------------------------\n"
        << "  Sequences & Segments Detail:\n";

    auto &seqs = mBlk->getSequences();
    std::vector<std::string> seq_ids;
    for (const auto &p : seqs) {
        seq_ids.push_back(p.first);
    }
    std::sort(seq_ids.begin(), seq_ids.end());

    for (size_t i = 0; i < seq_ids.size(); ++i) {
        const std::string &seq_id = seq_ids[i];
        bool is_last_seq = (i == seq_ids.size() - 1);
        std::string seq_branch = is_last_seq ? "  └─ " : "  ├─ ";
        std::string seq_indent = is_last_seq ? "     " : "  │  ";

        oss << seq_branch << "Sequence: " << seq_id << "\n";

        auto &segments = seqs.at(seq_id).getSegments();
        size_t seg_idx = 0;
        for (auto &seg_pair : segments) {
            bool is_last_seg = (++seg_idx == segments.size());
            std::string seg_branch = is_last_seg ? "└─ " : "├─ ";
            std::string seg_indent = is_last_seg ? "   " : "│  ";

            auto &seg = seg_pair.second;
            oss << seq_indent << seg_branch << "Segment: [" << seg.getStart() << " -> " << seg.getEnd() << "]"
                << " | Strand: " << (seg.isReverse() ? "-" : "+")
                << " | Copy: " << seg.getCopyCount() << "\n";

            auto &variants = seg.getVariants();
            if (variants.empty()) {
                oss << seq_indent << seg_indent << "└─ Variants: (none)\n";
            } else {
                oss << seq_indent << seg_indent << "├─ Variants (" << variants.size() << "):\n";
                for (size_t v_idx = 0; v_idx < variants.size(); ++v_idx) {
                    auto &var = variants[v_idx];
                    bool is_last_var = (v_idx == variants.size() - 1);
                    std::string var_branch = is_last_var ? "└─ " : "├─ ";

                    oss << seq_indent << seg_indent << "│  " << var_branch;
                    if (var.getType() == VariantType::SNV) {
                        oss << "SNV pos: " << var.getStart() << ", alt: '" << var.getAlt() << "'\n";
                    } else if (var.getType() == VariantType::GAP) {
                        oss << "GAP range: [" << var.getStart() << " -> " << var.getEnd() << "]\n";
                    } else {
                        oss << "Variant range: [" << var.getStart() << " -> " << var.getEnd() << "]\n";
                    }
                }
            }
        }
    }
    oss << "================================================================================";

    mga::msg::debug("%s", oss.str().c_str());
}

void MergeDebugger::validateAfterMerge(BlockSet *mergedSet, BlockManager *mgr) {
    if (!config_.enabled || !config_.validations || !mergedSet) return;
    mergedSet->debugValidateCopies(true);
    if (mgr) {
        mergedSet->debugValidateSequences(mgr);
    }
}

void MergeDebugger::validateCheckpoint(const std::string &label, BlockSet *mergedSet, BlockManager *mgr) {
    if (!config_.enabled || !config_.validations || !mergedSet) return;
    mga::msg::debug("\n  %s", label.c_str());
    if (mgr) {
        mergedSet->debugValidateSequences(mgr);
    }
}

} // namespace mga
