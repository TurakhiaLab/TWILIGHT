#include "alignment_debug.hpp"
#include "message_manager.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>

namespace mga {

AlignmentDebugger::AlignmentDebugger(bool enabled) {
    if (enabled) {
        setPreset(AlignmentDebugPreset::STANDARD);
    } else {
        setPreset(AlignmentDebugPreset::OFF);
    }
}

AlignmentDebugger::AlignmentDebugger(AlignmentDebugPreset preset) {
    setPreset(preset);
}

AlignmentDebugger::AlignmentDebugger(const AlignmentDebugConfig &config)
    : config_(config) {}

void AlignmentDebugger::setPreset(AlignmentDebugPreset preset) {
    switch (preset) {
        case AlignmentDebugPreset::OFF:
            config_.enabled = false;
            break;
        case AlignmentDebugPreset::SUMMARY:
            config_.enabled = true;
            config_.dequeued = true;
            config_.stage0_snaps = false;
            config_.strict_borders = false;
            config_.cut_points = false;
            config_.clusters = false;
            config_.sandboxes = false;
            config_.candidate_scores = false;
            config_.winners = true;
            config_.final_fragments = false;
            break;
        case AlignmentDebugPreset::STANDARD:
            config_.enabled = true;
            config_.dequeued = true;
            config_.stage0_snaps = false;
            config_.strict_borders = false;
            config_.cut_points = false;
            config_.clusters = true;
            config_.sandboxes = false;
            config_.candidate_scores = false;
            config_.winners = true;
            config_.final_fragments = false;
            break;
        case AlignmentDebugPreset::DETAILED:
            config_.enabled = true;
            config_.dequeued = true;
            config_.stage0_snaps = true;
            config_.strict_borders = true;
            config_.cut_points = true;
            config_.clusters = true;
            config_.sandboxes = true;
            config_.candidate_scores = true;
            config_.winners = true;
            config_.final_fragments = true;
            break;
        case AlignmentDebugPreset::FULL:
            config_.enabled = true;
            config_.dequeued = true;
            config_.stage0_snaps = true;
            config_.strict_borders = true;
            config_.cut_points = true;
            config_.clusters = true;
            config_.sandboxes = true;
            config_.candidate_scores = true;
            config_.winners = true;
            config_.final_fragments = true;
            break;
    }
}

void AlignmentDebugger::logDequeuedAlignment(double energy, int id, bool primary,
                                             int r_start, int r_end, int q_start, int q_end,
                                             bool inverse, int alnLength, const CigarString &cigar) {
    if (!config_.enabled || !config_.dequeued) return;
    int r_len = r_end - r_start;
    int q_len = q_end - q_start;

    std::ostringstream oss;
    oss << "\n=================================================================\n"
        << "🎯 [NEW ALIGNMENT DEQUEUED] Score/Energy: " << energy
        << " | ID: " << id << " | TYPE: " << (primary ? "PRIMARY" : "SECONDARY") << "\n"
        << "   ├─ Ref Range: [" << r_start << ", " << r_end << ") (Len: " << r_len << ")\n"
        << "   ├─ Qry Range: [" << q_start << ", " << q_end << ") (Len: " << q_len << ")\n"
        << "   ├─ Strand   : " << (inverse ? "Reverse (-)" : "Forward (+)") << "\n"
        << "   ├─ Aln Length: " << alnLength << "\n"
        << "   └─ Raw CIGAR: " << cigarToStr(cigar);
    mga::msg::debug("%s", oss.str().c_str());
}

void AlignmentDebugger::logSealingWall(bool isLeft, int rDist, int qDist, bool useSemi, const CigarString &patchCigar) {
    if (!config_.enabled || !config_.stage0_snaps) return;
    std::string side = isLeft ? "Left" : "Right";
    std::string method = useSemi ? " [Using Semi-Global]" : " [Using Global]";
    std::string cigarStr = cigarToStr(patchCigar);
    mga::msg::debug("  -> [STAGE 0] 🔗 Micro-overlap resolved. Sealing %s Wall (%d,%d)...%s -> CIGAR: %s",
                    side.c_str(), rDist, qDist, method.c_str(), cigarStr.c_str());
}

void AlignmentDebugger::logEdgeGapPadding(bool isLeft, bool isInsertion) {
    if (!config_.enabled || !config_.stage0_snaps) return;
    std::string side = isLeft ? "(Left)" : "(Right)";
    std::string type = isInsertion ? "Insertion Gap" : "Deletion Gap";
    mga::msg::debug("  -> [STAGE 0] 🩹 Appending %s %s to edge...", type.c_str(), side.c_str());
}

void AlignmentDebugger::logFinalFragmentSnap(int r_pad_L, int r_pad_R, int q_pad_L, int q_pad_R) {
    if (!config_.enabled || !config_.final_fragments) return;
    mga::msg::debug("    -> [DEBUG] 🧲 SNAPPING Final Fragment! Ref: (+%d, +%d), Qry: (+%d, +%d)",
                    r_pad_L, r_pad_R, q_pad_L, q_pad_R);
}

void AlignmentDebugger::logStrictCoverageDetected() {
    if (!config_.enabled || !config_.strict_borders) return;
    mga::msg::debug("  [DEBUG-STAGE 1-1] 💎 STRICT Coverage Edge Detected!");
}

void AlignmentDebugger::logStrictBorderHit(bool isRef, int pos, int mappedPos) {
    if (!config_.enabled || !config_.strict_borders) return;
    std::string mappedStr = (mappedPos != -1 ? std::to_string(mappedPos) : "Gap/Out");
    if (isRef) {
        mga::msg::debug("     [STRICT-REF] Hit Coverage Border @ Ref: %d | Mapped Qry: %s (Likely due to trailing/heading Overlap)",
                        pos, mappedStr.c_str());
    } else {
        mga::msg::debug("     [STRICT-QRY] Hit Coverage Border @ Qry: %d | Mapped Ref: %s (Likely due to trailing/heading Overlap)",
                        pos, mappedStr.c_str());
    }
}

void AlignmentDebugger::logCollectingCutPoints() {
    if (!config_.enabled || !config_.cut_points) return;
    mga::msg::debug("  [DEBUG-STAGE 1&2] Collecting Unified Cut Points...");
}

void AlignmentDebugger::logCollectingDynamicBoundaries() {
    if (!config_.enabled || !config_.cut_points) return;
    mga::msg::debug("  [DEBUG-STAGE 1-3] Collecting Dynamic Graph Boundaries from CoordinateManager...");
}

void AlignmentDebugger::logDynamicBoundaryRecovered(bool isRef, int bound, int mappedPos) {
    if (!config_.enabled || !config_.strict_borders) return;
    std::string mappedStr = (mappedPos != -1 ? std::to_string(mappedPos) : "Gap/Out");
    if (isRef) {
        mga::msg::debug("    ⚠️ [STRICT-DYNAMIC-REF] Recovered Topological Boundary @ Ref: %d | Mapped Qry: %s",
                        bound, mappedStr.c_str());
    } else {
        mga::msg::debug("    ⚠️ [STRICT-DYNAMIC-QRY] Recovered Topological Boundary @ Qry: %d | Mapped Ref: %s",
                        bound, mappedStr.c_str());
    }
}

void AlignmentDebugger::logClusterDetected(int cluster_size, int L_min,
                                          const std::vector<std::pair<std::pair<int, int>, CutPoint>> &cuts,
                                          int start_idx, int end_idx) {
    if (!config_.enabled || !config_.clusters || cluster_size <= 1) return;
    std::ostringstream oss;
    oss << "    -> 🎯 [CLUSTER DETECTED] B-D Range contains " << cluster_size << " candidate(s) [L_min: " << L_min << "]:\n";
    for (int k = start_idx; k <= end_idx; ++k) {
        const auto &cut = cuts[k];
        bool is_last = (k == end_idx);
        oss << "       " << (is_last ? "└─ " : "├─ ")
            << "[" << (k - start_idx + 1) << "/" << cluster_size << "] "
            << "Ref: " << cut.first.first << ", Qry: " << cut.first.second
            << " | Type: " << (cut.second.is_strict ? "STRICT" : "FLEXIBLE")
            << " | Reason: " << cut.second.reasons;
        if (!is_last) oss << "\n";
    }
    mga::msg::debug("%s", oss.str().c_str());
}

void AlignmentDebugger::logBoundaryCutStrict(int r, int q) {
    if (!config_.enabled || !config_.clusters) return;
    mga::msg::debug("    -> ⚠️ [BOUNDARY-STRICT] Cut @ R:%d, Q:%d is near boundary but STRICT. Fast-Path granted.", r, q);
}

void AlignmentDebugger::logBoundaryCutTrash(int r, int q) {
    if (!config_.enabled || !config_.clusters) return;
    mga::msg::debug("    -> 🗑️ [BOUNDARY-TRASH] Cut @ R:%d, Q:%d is fuzzy edge noise (< 50bp). Discarding completely!", r, q);
}

void AlignmentDebugger::logIsolatedFastPath(int r, int q) {
    if (!config_.enabled || !config_.clusters) return;
    mga::msg::debug("    -> ⚡ [FAST-PATH] Isolated WGA Cut @ R:%d, Q:%d", r, q);
}

void AlignmentDebugger::logSandboxDetails(int A_R, int E_R, int qry_box_start, int qry_box_end, int merged_consensus_len) {
    if (!config_.enabled || !config_.sandboxes) return;
    std::ostringstream oss;
    oss << "    -> 📦 Extracting and Merging A-E Sandbox:\n"
        << "       ├─ Ref Bounds: [" << A_R << " -> " << E_R << "] (Len: " << (E_R - A_R) << ")\n"
        << "       ├─ Qry Bounds: [" << qry_box_start << " -> " << qry_box_end << "] (Len: " << (qry_box_end - qry_box_start) << ")\n"
        << "       └─ Merged Consensus Length: " << merged_consensus_len;
    mga::msg::debug("%s", oss.str().c_str());
}

void AlignmentDebugger::logCandidateScore(int r, int q, double final_score) {
    if (!config_.enabled || !config_.candidate_scores) return;
    if (final_score < -5000) {
        mga::msg::debug("      [CANDIDATE] R:%d Q:%d -> Final: -INF", r, q);
    } else {
        mga::msg::debug("      [CANDIDATE] R:%d Q:%d -> Final: %.2f", r, q, final_score);
    }
}

void AlignmentDebugger::logStrictOverride(const std::vector<std::pair<int, int>> &strict_points) {
    if (!config_.enabled || !config_.winners) return;
    mga::msg::debug("      ---> ⚠️ [STRICT OVERRIDE] %zu STRICT point(s) found in cluster. Forcing cuts to prevent topology collapse.",
                    strict_points.size());
}

void AlignmentDebugger::logClusterWinner(int best_cut, int best_q_cut, double max_score) {
    if (!config_.enabled || !config_.winners) return;
    mga::msg::debug("      ---> 🎉 [WINNER] Snapped to R:%d Q:%d with Max Score: %.2f (Best in Cluster)",
                    best_cut, best_q_cut, max_score);
}

void AlignmentDebugger::logProcessingCutArray() {
    if (!config_.enabled || !config_.cut_points) return;
    mga::msg::debug("  [DEBUG-STAGE 3] Processing 1D Flattened Cut Array...");
}

} // namespace mga
