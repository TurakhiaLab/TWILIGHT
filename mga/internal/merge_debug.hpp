#pragma once

#include "type.hpp"
#include "cigar_util.hpp"
#include "coordinate_manager.hpp"
#include <string>
#include <vector>
#include <chrono>
#include <utility>

class Block;
class BlockSet;
class BlockManager;

namespace mga {

/**
 * @brief Preset debugging profiles for block merging.
 */
enum class MergeDebugPreset {
    OFF,        // Fully disabled (zero overhead)
    SUMMARY,    // High-level phases and execution summary only
    STANDARD,   // Standard debugging (phases, rejections, scenario classification, merge results)
    DETAILED,   // Standard + per-alignment details, energy values, interval dumps
    FULL        // All logging including deep Case D dump and internal graph validations
};

/**
 * @brief Granular category flags for fine-grained debug message control.
 */
struct MergeDebugConfig {
    bool enabled = false;          // Global master switch
    bool phases = true;           // Phase progress headers & timing
    bool aln_details = false;     // Per-alignment coordinates, CIGAR, overlaps
    bool skips = false;           // Null alignments, pure indels, low secondary scores
    bool energy_check = false;    // Energy calculations and graph fragmentation rejects
    bool core_extract = false;    // Core block extraction and coordinate intervals
    bool crossing_check = false;  // Anchor crossing checks and shared sequence detection
    bool scenarios = false;       // Scenario A1-A3, B, C, D classification & CIGAR patching
    bool merge_results = false;   // Created/merged block IDs and lengths
    bool case_d_dump = false;     // Comprehensive segment & variant dump for Case D
    bool validations = false;     // Internal consistency validations (copies & sequences)
};

/**
 * @brief Logger and diagnostic manager for BlockManager::merge.
 * Encapsulates formatted outputs using mga::msg::debug and eliminates scattered if (DEBUG_MODE) checks.
 */
class MergeDebugger {
public:
    MergeDebugger(bool enabled = false);
    MergeDebugger(MergeDebugPreset preset);
    MergeDebugger(const MergeDebugConfig &config);

    void setConfig(const MergeDebugConfig &config) { config_ = config; }
    void setPreset(MergeDebugPreset preset);
    const MergeDebugConfig& getConfig() const { return config_; }
    bool isEnabled() const { return config_.enabled; }

    // Phase & lifecycle logging
    void logMergeHeader(const BlockSetID &refId, const BlockSetID &qryId);
    void logPhase(int phaseNumber, const std::string &description);
    void logProcessedAlignments(int count);
    void logMergeSummary(std::chrono::high_resolution_clock::time_point startTime);

    // Alignment inspection
    void logAlignmentDetails(const CigarString &cigar, bool isInverse,
                             int r_start, int r_end, int q_start, int q_end,
                             bool r_merged, bool q_merged,
                             const std::vector<BlockID> &r_overlaps,
                             const std::vector<BlockID> &q_overlaps);

    // Skip & bypass logging
    void logSkipIdenticalBlock(BlockID blkId);
    void logNullAlignment();
    void logPureIndelBypass();
    void logSecondaryScoreSkip(int score, int threshold);

    // Energy check & rejection logging
    void logOccurrenceCrossingReject(bool isExtractCore = false);
    void logEnergyCalculation(bool isRef, int length, int Nc, int energy);
    void logEnergyReject();
    void logMultiBlockSpanningReject(size_t r_count, size_t q_count);
    void logTopologySkip();
    void logLowIdentityReject(double identity);
    void logInvalidCoreBlocks(BlockID rCore, BlockID qCore);

    // Core block extraction & intervals
    void logCoreExtraction(const std::string &trackName, BlockID coreId,
                           const CoordinateManager::GlobalCoords &coords);

    // Topology & crossing check
    void logAnchorNeighbors(const std::pair<int, int> &ref_before,
                            const std::pair<int, int> &ref_after,
                            const std::pair<int, int> &qry_before,
                            const std::pair<int, int> &qry_after);

    // Scenario classification & CIGAR patching/realigning
    void logScenario(const std::string &scenarioName, const std::string &details = "");
    void logCaseAPatchedCigar(const CigarString &cigar);
    void logRealigning(size_t rLen, size_t qLen);
    void logRealignedCigar(const CigarString &cigar, double identity);
    void logCaseDConsensus(const std::string &rCons, const std::string &qCons);

    // Merged block results & deep dumps
    void logMergedBlock(BlockID mId, size_t length);
    void logCaseDDetails(Block *mBlk);

    // Diagnostics & validation
    void validateAfterMerge(BlockSet *mergedSet, BlockManager *mgr);
    void validateCheckpoint(const std::string &label, BlockSet *mergedSet, BlockManager *mgr);

private:
    MergeDebugConfig config_;
};

} // namespace mga
