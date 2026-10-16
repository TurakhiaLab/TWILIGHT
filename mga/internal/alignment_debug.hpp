#pragma once

#include "type.hpp"
#include "cigar_util.hpp"
#include <string>
#include <vector>
#include <utility>

struct CutPoint {
    bool is_strict = false;
    std::string reasons = "";
};

namespace mga {

/**
 * @brief Preset debugging profiles for alignment segmentation and cut cluster evaluation.
 */
enum class AlignmentDebugPreset {
    OFF,        // Fully disabled (zero overhead)
    SUMMARY,    // High-level steps & cluster winners only
    STANDARD,   // Standard debugging (dequeued alignments, clusters, overrides, winners)
    DETAILED,   // Standard + Stage 0 snaps, candidate scores, sandboxes
    FULL        // All logging including coverage border hits and dynamic boundary recoveries
};

/**
 * @brief Granular category flags for fine-grained alignment debug output.
 */
struct AlignmentDebugConfig {
    bool enabled = false;           // Master switch
    bool dequeued = false;          // Dequeued alignment headers & ranges
    bool stage0_snaps = false;      // Micro-overlap resolving, flank alignment, edge gap padding
    bool strict_borders = false;    // STRICT coverage edge hits and dynamic graph boundary detections
    bool cut_points = false;        // Unified cut points collection stages
    bool clusters = false;          // Cluster grouping, branch displays, and fast-paths
    bool sandboxes = false;         // A-E sandbox ranges and merged consensus extraction
    bool candidate_scores = false;  // Per-candidate split scores
    bool winners = false;           // Cluster winners and strict overrides
    bool final_fragments = false;   // Final fragment snapping
};

/**
 * @brief Logger and diagnostic manager for AlignmentCollection and alignment segmentation.
 */
class AlignmentDebugger {
public:
    AlignmentDebugger(bool enabled = false);
    AlignmentDebugger(AlignmentDebugPreset preset);
    AlignmentDebugger(const AlignmentDebugConfig &config);

    void setConfig(const AlignmentDebugConfig &config) { config_ = config; }
    void setPreset(AlignmentDebugPreset preset);
    const AlignmentDebugConfig& getConfig() const { return config_; }
    bool isEnabled() const { return config_.enabled; }

    // Dequeued alignment logs
    void logDequeuedAlignment(double energy, int id, bool primary,
                              int r_start, int r_end, int q_start, int q_end,
                              bool inverse, int alnLength, const CigarString &cigar);

    // Stage 0 edge sealing & gap padding
    void logSealingWall(bool isLeft, int rL, int qL, bool useSemi, const CigarString &patchCigar);
    void logEdgeGapPadding(bool isLeft, bool isInsertion);
    void logFinalFragmentSnap(int r_pad_L, int r_pad_R, int q_pad_L, int q_pad_R);

    // Strict coverage and dynamic boundaries
    void logStrictCoverageDetected();
    void logStrictBorderHit(bool isRef, int pos, int mappedPos);
    void logCollectingCutPoints();
    void logCollectingDynamicBoundaries();
    void logDynamicBoundaryRecovered(bool isRef, int bound, int mappedPos);

    // Cut clusters & candidates
    void logClusterDetected(int cluster_size, int L_min,
                            const std::vector<std::pair<std::pair<int, int>, CutPoint>> &cuts,
                            int start_idx, int end_idx);
    void logBoundaryCutStrict(int r, int q);
    void logBoundaryCutTrash(int r, int q);
    void logIsolatedFastPath(int r, int q);

    // Sandbox & score selection
    void logSandboxDetails(int A_R, int E_R, int qry_box_start, int qry_box_end, int merged_consensus_len);
    void logCandidateScore(int r, int q, double final_score);
    void logStrictOverride(const std::vector<std::pair<int, int>> &strict_points);
    void logClusterWinner(int best_cut, int best_q_cut, double max_score);

    // Processing steps
    void logProcessingCutArray();

private:
    AlignmentDebugConfig config_;
};

} // namespace mga
