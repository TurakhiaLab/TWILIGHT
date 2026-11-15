#include "score.hpp"
#include "msa.hpp"
#include "blosum.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <zlib.h>

namespace msa {
namespace score {

namespace {

inline bool isGapChar(char c) {
    return c == '-' || c == '.' || c == ' ' || c == '\0' || c == '\n' || c == '\r';
}

struct ColumnStats {
    double identicalPairs = 0.0;
    double wspSubScore = 0.0;
    double wspGapScore = 0.0;
    int conservedCols = 0;
    int totalGaps = 0;
};

// Fallback FASTA reader in case SequenceDB has no memory buffers
std::vector<std::string> loadAlignmentFromFasta(const std::string& filePath) {
    std::vector<std::string> sequences;
    bool isGz = (filePath.size() >= 3 && filePath.substr(filePath.size() - 3) == ".gz");

    if (isGz) {
        gzFile gz = gzopen(filePath.c_str(), "rb");
        if (!gz) return sequences;
        std::string currentSeq = "";
        char buf[8192];
        std::string line = "";
        while (gzgets(gz, buf, sizeof(buf)) != NULL) {
            std::string s(buf);
            for (char c : s) {
                if (c == '\n') {
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    if (!line.empty()) {
                        if (line[0] == '>') {
                            if (!currentSeq.empty()) {
                                sequences.push_back(currentSeq);
                                currentSeq.clear();
                            }
                        } else {
                            currentSeq += line;
                        }
                    }
                    line.clear();
                } else {
                    line += c;
                }
            }
        }
        if (!line.empty() && line[0] != '>') currentSeq += line;
        if (!currentSeq.empty()) sequences.push_back(currentSeq);
        gzclose(gz);
    } else {
        std::ifstream file(filePath);
        if (!file.is_open()) return sequences;
        std::string line;
        std::string currentSeq = "";
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (!currentSeq.empty()) {
                    sequences.push_back(currentSeq);
                    currentSeq.clear();
                }
            } else {
                currentSeq += line;
            }
        }
        if (!currentSeq.empty()) sequences.push_back(currentSeq);
        file.close();
    }
    return sequences;
}

} // anonymous namespace

AlignmentScoreResult computeAlignmentScore(
    SequenceDB* database,
    Option* option,
    const Params* param,
    int alnLen,
    int scoreMode
) {
    AlignmentScoreResult result;
    result.mode = scoreMode;

    // Collect aligned sequence buffers
    std::vector<const char*> seqBuffers;
    std::vector<std::string> fallbackSeqs;

    if (database && !database->sequences.empty()) {
        seqBuffers.reserve(database->sequences.size());
        for (auto seq : database->sequences) {
            if (seq && !seq->lowQuality && seq->alnStorage[seq->storage] != nullptr) {
                seqBuffers.push_back(&seq->alnStorage[seq->storage][0]);
            }
        }
    }

    // Fallback: read from output file if memory buffer is unavailable
    if (seqBuffers.empty() && option && !option->outFile.empty()) {
        std::string path = option->outFile;
        if (option->compressed && path.substr(path.size() - 3) != ".gz") {
            path += ".gz";
        }
        fallbackSeqs = loadAlignmentFromFasta(path);
        seqBuffers.reserve(fallbackSeqs.size());
        for (const auto& s : fallbackSeqs) {
            seqBuffers.push_back(s.data());
            if (alnLen <= 0) alnLen = static_cast<int>(s.length());
        }
    }

    const int numSeqs = static_cast<int>(seqBuffers.size());
    result.numSeqs = numSeqs;
    result.alnLen = alnLen;

    if (numSeqs < 2 || alnLen <= 0) {
        return result;
    }

    const char type = (option) ? option->type : 'n';
    const int alphabetSize = (type == 'p') ? 20 : 4;
    const int maxMatrixDim = (type == 'p') ? 21 : 5;

    // Calculate unaligned sequence lengths & total non-gap characters
    double totalNonGaps = 0.0;
    for (int s = 0; s < numSeqs; ++s) {
        int nonGapCount = 0;
        const char* buf = seqBuffers[s];
        for (int c = 0; c < alnLen; ++c) {
            if (!isGapChar(buf[c])) nonGapCount++;
        }
        totalNonGaps += nonGapCount;
    }

    const double meanSeqLen = (numSeqs > 0) ? (totalNonGaps / numSeqs) : 0.0;
    result.meanSeqLen = meanSeqLen;
    const double numPairs = (static_cast<double>(numSeqs) * (numSeqs - 1)) / 2.0;
    const double theoreticalMaxPairs = numPairs * meanSeqLen;

    // Setup substitution matrix and gap parameters for Mode 2 (WSP)
    float localMatrix[21][21];
    float gapOpen = (type == 'n') ? -50.0f : -1080.0f;
    float gapExtend = (type == 'n') ? -5.0f : -108.0f;
    int matrixDim = maxMatrixDim;

    if (param && param->scoringMatrix) {
        gapOpen = param->gapOpen;
        gapExtend = param->gapExtend;
        matrixDim = std::min(param->matrixSize, 21);
        for (int i = 0; i < matrixDim; ++i) {
            for (int j = 0; j < matrixDim; ++j) {
                localMatrix[i][j] = param->scoringMatrix[i][j];
            }
        }
    } else {
        // Fallback scoring matrix
        if (type == 'n') {
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    if (i == 4 || j == 4) localMatrix[i][j] = 0.0f;
                    else if (i == j)      localMatrix[i][j] = 18.0f;
                    else if (std::abs(i - j) == 2) localMatrix[i][j] = -4.0f;
                    else                  localMatrix[i][j] = -8.0f;
                }
            }
        } else {
            for (int i = 0; i < 20; ++i) {
                for (int j = 0; j < 20; ++j) {
                    localMatrix[i][j] = BLOSUM62[i][j];
                }
                localMatrix[i][20] = 0.0f;
                localMatrix[20][i] = 0.0f;
            }
            localMatrix[20][20] = 0.0f;
        }
    }

    // Parallel column evaluation across alignment
    ColumnStats stats = tbb::parallel_reduce(
        tbb::blocked_range<int>(0, alnLen, 256),
        ColumnStats(),
        [&](const tbb::blocked_range<int>& range, ColumnStats localStats) -> ColumnStats {
            int counts[21];
            for (int c = range.begin(); c != range.end(); ++c) {
                std::fill(counts, counts + maxMatrixDim, 0);
                int numGaps = 0;
                int numRes = 0;
                int numGapOpen = 0;
                int numGapExt = 0;

                for (int s = 0; s < numSeqs; ++s) {
                    char ch = seqBuffers[s][c];
                    if (isGapChar(ch)) {
                        numGaps++;
                        char prev = (c > 0) ? seqBuffers[s][c - 1] : 'A';
                        if (c == 0 || !isGapChar(prev)) {
                            numGapOpen++;
                        } else {
                            numGapExt++;
                        }
                    } else {
                        numRes++;
                        int idx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                        if (idx >= 0 && idx < maxMatrixDim) {
                            counts[idx]++;
                        }
                    }
                }

                localStats.totalGaps += numGaps;

                // Conservation check
                int maxResCount = 0;
                for (int a = 0; a < alphabetSize; ++a) {
                    if (counts[a] > maxResCount) maxResCount = counts[a];
                }
                if (numSeqs > 0 && (static_cast<double>(maxResCount) / numSeqs) >= 0.70) {
                    localStats.conservedCols++;
                }

                // Mode 1: Consistency column calculation
                double idPairs = 0.0;
                for (int a = 0; a < alphabetSize; ++a) {
                    if (counts[a] > 1) {
                        idPairs += (static_cast<double>(counts[a]) * (counts[a] - 1)) / 2.0;
                    }
                }
                localStats.identicalPairs += idPairs;

                // Mode 2: WSP Profile calculation (always computed for comprehensive scoring)
                double colSub = 0.0;
                for (int a = 0; a < matrixDim; ++a) {
                    if (counts[a] == 0) continue;
                    // a == b self-pairs
                    if (counts[a] > 1) {
                        colSub += (static_cast<double>(counts[a]) * (counts[a] - 1) / 2.0) * localMatrix[a][a];
                    }
                    // a < b cross-pairs
                    for (int b = a + 1; b < matrixDim; ++b) {
                        if (counts[b] > 0) {
                            colSub += static_cast<double>(counts[a]) * counts[b] * localMatrix[a][b];
                        }
                    }
                }
                localStats.wspSubScore += colSub;

                // Affine gap penalties against non-gap residues
                if (numRes > 0) {
                    double gapPen = static_cast<double>(numRes) * (numGapOpen * gapOpen + numGapExt * gapExtend);
                    localStats.wspGapScore += gapPen;
                }
            }
            return localStats;
        },
        [](ColumnStats a, const ColumnStats& b) -> ColumnStats {
            a.identicalPairs += b.identicalPairs;
            a.wspSubScore += b.wspSubScore;
            a.wspGapScore += b.wspGapScore;
            a.conservedCols += b.conservedCols;
            a.totalGaps += b.totalGaps;
            return a;
        }
    );

    const double totalCells = static_cast<double>(numSeqs) * alnLen;
    result.gapRatio = (totalCells > 0.0) ? (stats.totalGaps / totalCells) : 0.0;
    result.conservedColsRatio = (alnLen > 0) ? (static_cast<double>(stats.conservedCols) / alnLen) : 0.0;

    // Consistency Objective Score
    result.consistencyScore = (theoreticalMaxPairs > 0.0) ? (stats.identicalPairs / theoreticalMaxPairs) : 0.0;

    // WSP (Weighted Profile Sum-of-Pairs)
    double totalWSP = stats.wspSubScore + stats.wspGapScore;
    result.rawWspScore = totalWSP;
    result.wspScore = (theoreticalMaxPairs > 0.0) ? (totalWSP / theoreticalMaxPairs) : 0.0;

    if (scoreMode == 1) {
        result.rawScore = stats.identicalPairs;
        result.normalizedScore = result.consistencyScore;
        result.score = result.consistencyScore;
    } else {
        result.rawScore = totalWSP;
        result.normalizedScore = result.wspScore;
        result.score = result.wspScore;
    }

    return result;
}

void writeScoreToFile(
    const std::string& scoreFilePath,
    const AlignmentScoreResult& result
) {
    if (scoreFilePath.empty()) return;

    std::ofstream out(scoreFilePath);
    if (!out.is_open()) {
        std::cerr << "WARNING: Could not open score file for writing: " << scoreFilePath << "\n";
        return;
    }

    // Line 1: Primary score numeric value (for instant 1-line script reading: float(open("out.score").readline().strip()))
    out << std::fixed << std::setprecision(6) << result.score << "\n";

    // Structured metadata
    out << "SCORE=" << std::fixed << std::setprecision(6) << result.score << "\n";
    out << "RAW_SCORE=" << std::fixed << std::setprecision(6) << result.rawScore << "\n";
    out << "NORMALIZED_SCORE=" << std::fixed << std::setprecision(6) << result.normalizedScore << "\n";
    out << "CONSISTENCY_SCORE=" << std::fixed << std::setprecision(6) << result.consistencyScore << "\n";
    out << "WSP_SCORE=" << std::fixed << std::setprecision(6) << result.wspScore << "\n";
    out << "WSP_RAW_SCORE=" << std::fixed << std::setprecision(6) << result.rawWspScore << "\n";
    out << "METHOD=" << (result.mode == 1 ? "CONSISTENCY" : "WSP") << "\n";
    out << "ALIGNMENT_LENGTH=" << result.alnLen << "\n";
    out << "NUM_SEQUENCES=" << result.numSeqs << "\n";
    out << "GAP_RATIO=" << std::fixed << std::setprecision(6) << result.gapRatio << "\n";
    out << "CONSERVED_COLS_RATIO=" << std::fixed << std::setprecision(6) << result.conservedColsRatio << "\n";
    out << "MEAN_SEQUENCE_LENGTH=" << std::fixed << std::setprecision(2) << result.meanSeqLen << "\n";

    out.close();
}

void evaluateAndOutputScore(
    SequenceDB* database,
    Option* option,
    const Params* param,
    int alnLen
) {
    if (!option || option->scoreMode <= 0) return;

    auto start = std::chrono::high_resolution_clock::now();
    AlignmentScoreResult result = computeAlignmentScore(database, option, param, alnLen, option->scoreMode);
    writeScoreToFile(option->scoreFile, result);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cerr << "Alignment Score (" << (result.mode == 1 ? "Consistency" : "WSP") 
              << "): " << std::fixed << std::setprecision(6) << result.score 
              << " (written to " << option->scoreFile << " in " << ms.count() << " ms)\n";
}

} // namespace score
} // namespace msa

