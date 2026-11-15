#ifndef SCORE_HPP
#define SCORE_HPP

#include "msa.hpp"
#include <string>

namespace msa {
namespace score {

struct AlignmentScoreResult {
    int mode = 1;                    // 1: Consistency objective score, 2: WSP
    double score = 0.0;              // Primary score (for maximization and convergence)
    double rawScore = 0.0;           // Unnormalized raw score sum
    double normalizedScore = 0.0;    // Normalized score
    double consistencyScore = 0.0;   // Consistency objective score [0.0, 1.0]
    double wspScore = 0.0;           // Normalized WSP score
    double rawWspScore = 0.0;        // Raw total WSP score
    int alnLen = 0;                  // Alignment column length
    int numSeqs = 0;                 // Total sequence count
    double gapRatio = 0.0;           // Overall gap percentage [0.0, 1.0]
    double conservedColsRatio = 0.0; // Fraction of columns with high conservation
    double meanSeqLen = 0.0;         // Average unaligned sequence length
};

/**
 * Compute alignment score for the final MSA in SequenceDB.
 * mode: 1 = Consistency Objective Score, 2 = Profile WSP Score
 */
AlignmentScoreResult computeAlignmentScore(
    SequenceDB* database,
    Option* option,
    const Params* param,
    int alnLen,
    int scoreMode
);

/**
 * Write the score result to scoreFilePath.
 * Line 1 contains the primary score numeric value for instant single-line parsing in scripts.
 * Remaining lines contain structured metadata.
 */
void writeScoreToFile(
    const std::string& scoreFilePath,
    const AlignmentScoreResult& result
);

/**
 * High-level evaluator called at MSA finalization.
 */
void evaluateAndOutputScore(
    SequenceDB* database,
    Option* option,
    const Params* param,
    int alnLen
);

} // namespace score
} // namespace msa

#endif // SCORE_HPP
