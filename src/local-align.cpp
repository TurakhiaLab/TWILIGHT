#include "msa.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>

namespace msa {
namespace accurate {

static int residueScore(char referenceBase, char queryBase, char type, Params& params)
{
    const char ref = static_cast<char>(std::toupper(static_cast<unsigned char>(referenceBase)));
    const char qry = static_cast<char>(std::toupper(static_cast<unsigned char>(queryBase)));
    int refIndex = letterIdx(type, toupper(referenceBase));
    int qryIndex = letterIdx(type, toupper(queryBase));
    return params.scoringMatrix[refIndex][qryIndex];
}

AlignmentResult Aligner::align_affine (const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    const int gOpen = static_cast<int>(params.gapOpen);
    const int gExt  = static_cast<int>(params.gapExtend);

    const int MIN_INF = -100000000; 

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    // 1. Tiling Parameters
    const int T = 1000;
    const int O = 100;

    int ref_idx = 0;
    int qry_idx = 0;

    int max_ref_tile = std::min(T, totalRef);
    int max_qry_tile = std::min(T, totalQry);

    std::vector<std::vector<int>> score(max_ref_tile + 1, std::vector<int>(max_qry_tile + 1, 0));
    std::vector<std::vector<uint8_t>> traceback(max_ref_tile + 1, std::vector<uint8_t>(max_qry_tile + 1, 0));
    std::vector<int> E(max_qry_tile + 1, MIN_INF);

    std::vector<std::pair<int, int>> global_pairs;
    
    int last_max_score = 0;
    int last_max_state = 0;
    int final_alignment_score = 0;

    while (ref_idx < totalRef && qry_idx < totalQry) {
        int refLen = std::min(T, totalRef - ref_idx);
        int qryLen = std::min(T, totalQry - qry_idx);

        std::fill(E.begin(), E.end(), MIN_INF);

        if (ref_idx == 0 && qry_idx == 0) {
            // Semi-global: No gap penalty at the begining for the first tile
            for (int i = 0; i <= refLen; ++i) { score[i][0] = 0; traceback[i][0] = 0; }
            for (int j = 0; j <= qryLen; ++j) { score[0][j] = 0; traceback[0][j] = 0; }
        } else {
            score[0][0] = last_max_score;
            for (int i = 1; i <= refLen; ++i) {
                int penalty = (last_max_state == 1 && i == 1) ? gExt : (i == 1 ? gOpen : gExt);
                score[i][0] = score[i-1][0] + penalty;
                uint8_t tb = 2; // UP
                if (i > 1 || last_max_state == 1) tb |= 0x04; // affine gap extension
                traceback[i][0] = tb;
            }
            for (int j = 1; j <= qryLen; ++j) {
                int penalty = (last_max_state == 2 && j == 1) ? gExt : (j == 1 ? gOpen : gExt);
                score[0][j] = score[0][j-1] + penalty;
                uint8_t tb = 3; // LEFT
                if (j > 1 || last_max_state == 2) tb |= 0x08; // affine gap extension
                traceback[0][j] = tb;
            }
        }

        for (int i = 1; i <= refLen; ++i) {
            int F_val = MIN_INF;
            for (int j = 1; j <= qryLen; ++j) {
                uint8_t tb = 0;

                int e_open = score[i - 1][j] + gOpen; 
                int e_ext  = E[j] + gExt;
                if (e_open >= e_ext) {
                    E[j] = e_open;
                } else {
                    E[j] = e_ext;
                    tb |= 0x04;
                }

                int f_open = score[i][j - 1] + gOpen;
                int f_ext  = F_val + gExt;
                if (f_open >= f_ext) {
                    F_val = f_open;
                } else {
                    F_val = f_ext;
                    tb |= 0x08;
                }

                int diag = score[i - 1][j - 1] + residueScore(reference[ref_idx + i - 1], query[qry_idx + j - 1], type, params);

                int best = MIN_INF;
                uint8_t h_src = 0;

                if (diag > best) { best = diag; h_src = 1; }
                if (E[j] > best) { best = E[j]; h_src = 2; }
                if (F_val > best) { best = F_val; h_src = 3; }

                score[i][j] = best;
                tb |= h_src;
                traceback[i][j] = tb;
            }
        }

        bool is_last_tile_ref = (ref_idx + refLen == totalRef);
        bool is_last_tile_qry = (qry_idx + qryLen == totalQry);
        bool is_last_tile = is_last_tile_ref || is_last_tile_qry;

        int best_i = refLen;
        int best_j = qryLen;
        int max_score = MIN_INF;
        uint8_t best_tb = 0;

        if (is_last_tile) {
            // Semi-global: Find maximum at the border line (Free end gap)
            if (is_last_tile_ref) {
                for (int j = 0; j <= qryLen; ++j) {
                    if (score[refLen][j] > max_score) {
                        max_score = score[refLen][j];
                        best_i = refLen; best_j = j;
                        best_tb = traceback[refLen][j];
                    }
                }
            }
            if (is_last_tile_qry) {
                for (int i = 0; i <= refLen; ++i) {
                    if (score[i][qryLen] > max_score) {
                        max_score = score[i][qryLen];
                        best_i = i; best_j = qryLen;
                        best_tb = traceback[i][qryLen];
                    }
                }
            }
            final_alignment_score = max_score;
        } else {
            int start_i = std::max(1, refLen - O);
            int start_j = std::max(1, qryLen - O);
            
            // 1. Search for the L-shaped "bottom horizontal strip region" (Bottom Region)
            for (int i = start_i; i <= refLen; ++i) {
                for (int j = 1; j <= qryLen; ++j) {
                    if (score[i][j] > max_score) {
                        max_score = score[i][j];
                        best_i = i; best_j = j;
                        best_tb = traceback[i][j];
                    }
                }
            }
            
            // 2. Search for the L-shaped "right vertical strip region" (Right Region)
            for (int i = 1; i < start_i; ++i) {
                for (int j = start_j; j <= qryLen; ++j) {
                    if (score[i][j] > max_score) {
                        max_score = score[i][j];
                        best_i = i; best_j = j;
                        best_tb = traceback[i][j];
                    }
                }
            }
        }

        int i = best_i;
        int j = best_j;
        int currentState = 0;
        std::vector<std::pair<int, int>> local_pairs;

        while (i > 0 || j > 0) {
            if (ref_idx == 0 && qry_idx == 0 && score[i][j] == 0 && traceback[i][j] == 0) break;
            if (i == 0 && j == 0) break;
            uint8_t tb = traceback[i][j];
            if (currentState == 0) { 
                uint8_t h_src = tb & 0x03;
                if (h_src == 1) {
                    local_pairs.push_back({ref_idx + i - 1, qry_idx + j - 1});
                    --i; --j;
                } else if (h_src == 2) { 
                    currentState = 1; 
                } else if (h_src == 3) { 
                    currentState = 2;
                } else { 
                    if (i > 0 && j == 0) currentState = 1;
                    else if (j > 0 && i == 0) currentState = 2;
                    else break;
                }
            } 
            else if (currentState == 1) { 
                bool e_from_e = (tb & 0x04) != 0; 
                --i; 
                if (!e_from_e) currentState = 0; 
            } 
            else if (currentState == 2) { 
                bool f_from_f = (tb & 0x08) != 0; 
                --j; 
                if (!f_from_f) currentState = 0; 
            }
        }

        std::reverse(local_pairs.begin(), local_pairs.end());
        global_pairs.insert(global_pairs.end(), local_pairs.begin(), local_pairs.end());

        if (is_last_tile) break;

        ref_idx += best_i;
        qry_idx += best_j;
        last_max_score = max_score;

        uint8_t h_src = best_tb & 0x03;
        if (h_src == 1) last_max_state = 0;
        else if (h_src == 2) last_max_state = 1;
        else if (h_src == 3) last_max_state = 2;
        else last_max_state = 0;
    }

    AlignmentResult result;
    result.score = final_alignment_score;

    int identicalPairs = 0;
    for (const auto& alignedPair : global_pairs) {
        result.alignedPairs.push_back({alignedPair.first, alignedPair.second});
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    
    if (!result.alignedPairs.empty()) {
        result.identity = static_cast<float>(identicalPairs) / static_cast<float>(result.alignedPairs.size());
    }
    
    return result;
}

LocalHomTable::PairResult Aligner::align_affine_local_segments_banded(
    const std::string& reference,
    const std::string& query,
    char type,
    msa::Params& params,
    int bandWidth)
{
    const int gOpen = static_cast<int>(params.localGapOpen);
    const int gExt  = static_cast<int>(params.localGapExtend);
    const int MIN_INF = -100000000;

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    LocalHomTable::PairResult result;
    if (totalRef == 0 || totalQry == 0) {
        return result;
    }

    const size_t qrySize = static_cast<size_t>(totalQry + 1);
    if (score_0.size() < qrySize) score_0.resize(qrySize);
    if (score_1.size() < qrySize) score_1.resize(qrySize);
    if (score_2.size() < qrySize) score_2.resize(qrySize);
    if (E_buf.size() < qrySize) E_buf.resize(qrySize);
    if (qry_idx_buf.size() < static_cast<size_t>(totalQry)) qry_idx_buf.resize(totalQry);

    const size_t pitch = qrySize;
    const size_t tbSize = static_cast<size_t>(totalRef + 1) * pitch;
    if (tb_buf.size() < tbSize) tb_buf.resize(tbSize);

    // Initialize DP state arrays
    std::fill_n(score_0.data(), qrySize, 0);
    std::fill_n(score_1.data(), qrySize, 0);
    std::fill_n(score_2.data(), qrySize, 0);
    std::fill_n(E_buf.data(), qrySize, MIN_INF);

    int* scoreRows[3] = { score_0.data(), score_1.data(), score_2.data() };
    int* E = E_buf.data();
    uint8_t* traceback = tb_buf.data();

    // Precompute query residue indices
    for (int j = 0; j < totalQry; ++j) {
        qry_idx_buf[j] = letterIdx(type, query[j]);
    }
    const int* qryIndices = qry_idx_buf.data();

    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    const int halfBand = std::max(1, bandWidth / 2);
    int prev_j_min = 1;
    int prev_j_max = std::min(totalQry, halfBand);

    for (int i = 1; i <= totalRef; ++i) {
        int curr = i % 3;
        int prev = (i - 1 + 3) % 3;
        int* score_curr = scoreRows[curr];
        int* score_prev = scoreRows[prev];

        // Slanted center tracking along the main diagonal
        int j_center = static_cast<int>((static_cast<int64_t>(i) * totalQry) / totalRef);
        int j_min = std::max(1, j_center - halfBand);
        int j_max = std::min(totalQry, j_center + halfBand);

        // Sanitize score_curr: zero out boundary cells
        score_curr[0] = 0;
        if (j_min > 0) score_curr[j_min - 1] = 0;

        // Clean up E and score_prev for indices that fell out of the active band
        if (j_min > prev_j_min) {
            for (int k = prev_j_min; k < j_min; ++k) {
                E[k] = MIN_INF;
                score_prev[k] = 0;
            }
        }
        if (prev_j_max + 1 <= totalQry) {
            score_prev[prev_j_max + 1] = 0;
        }

        int refIdx = letterIdx(type, reference[i - 1]);
        const float* scoreRow = params.scoringMatrix[refIdx];
        size_t rowOffset = static_cast<size_t>(i) * pitch;

        int F_val = MIN_INF;
        for (int j = j_min; j <= j_max; ++j) {
            uint8_t tb = 0;

            int e_open = score_prev[j] + gOpen; 
            int e_ext  = E[j] + gExt;
            if (e_open >= e_ext) {
                E[j] = e_open;
            } else {
                E[j] = e_ext;
                tb |= 0x04;
            }

            int f_open = score_curr[j - 1] + gOpen;
            int f_ext  = F_val + gExt;
            if (f_open >= f_ext) {
                F_val = f_open;
            } else {
                F_val = f_ext;
                tb |= 0x08;
            }

            int match = static_cast<int>(scoreRow[qryIndices[j - 1]]);
            int diag = score_prev[j - 1] + match;

            int best = 0; 
            uint8_t h_src = 0;

            if (diag > best) { best = diag; h_src = 1; }
            if (E[j] > best) { best = E[j]; h_src = 2; }
            if (F_val > best) { best = F_val; h_src = 3; }

            score_curr[j] = best;
            tb |= h_src;
            traceback[rowOffset + j] = tb;

            if (best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }

        if (j_max + 1 <= totalQry) {
            score_curr[j_max + 1] = 0;
            E[j_max + 1] = MIN_INF;
        }

        prev_j_min = j_min;
        prev_j_max = j_max;
    }

    result.score = max_score;
    if (max_score <= 0) {
        return result;
    }

    // Traceback within band limits
    int i = max_i;
    int j = max_j;
    int currentState = 0;
    pair_idx_buf.clear();

    while (i > 0 && j > 0) {
        int j_center = static_cast<int>((static_cast<int64_t>(i) * totalQry) / totalRef);
        int j_min = std::max(1, j_center - halfBand);
        int j_max = std::min(totalQry, j_center + halfBand);
        if (j < j_min || j > j_max) {
            break; // Reached boundary of band, terminate local alignment gracefully
        }

        uint8_t tb = traceback[static_cast<size_t>(i) * pitch + j];
        
        if (currentState == 0) { 
            uint8_t h_src = tb & 0x03;
            if (h_src == 1) { // Diagonal
                pair_idx_buf.push_back({i - 1, j - 1});
                --i; --j;
            } else if (h_src == 2) { // Up (E)
                currentState = 1; 
            } else if (h_src == 3) { // Left (F)
                currentState = 2;
            } else { 
                break; // Local alignment ends when h_src is 0
            }
        } 
        else if (currentState == 1) { // Up (E)
            bool e_from_e = (tb & 0x04) != 0; 
            --i; 
            if (!e_from_e) currentState = 0;
        } 
        else if (currentState == 2) { // Left (F)
            bool f_from_f = (tb & 0x08) != 0; 
            --j; 
            if (!f_from_f) currentState = 0;
        }
    }

    if (pair_idx_buf.empty()) {
        return result;
    }

    std::reverse(pair_idx_buf.begin(), pair_idx_buf.end());

    // Compute percent identity
    int identicalPairs = 0;
    for (const auto& alignedPair : pair_idx_buf) {
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    // Compute ungapped substitution score sum across aligned pairs (MAFFT-compatible)
    int sum_subst_score = 0;
    for (const auto& p : pair_idx_buf) {
        int refIdx = letterIdx(type, reference[p.first]);
        int qryIdx = letterIdx(type, query[p.second]);
        sum_subst_score += static_cast<int>(params.scoringMatrix[refIdx][qryIdx]);
    }

    // Normalize pairwise substitution score by overlap length to obtain average per-residue score
    const float normOptScore = (!pair_idx_buf.empty()) 
        ? (static_cast<float>(sum_subst_score) / static_cast<float>(pair_idx_buf.size())) 
        : 0.0f;

    // Decompose contiguous match pairs into ungapped segments (putlocalhom2 equivalent)
    int segStart1 = pair_idx_buf[0].first;
    int segStart2 = pair_idx_buf[0].second;
    int prev1 = segStart1;
    int prev2 = segStart2;

    for (size_t k = 1; k < pair_idx_buf.size(); ++k) {
        int cur1 = pair_idx_buf[k].first;
        int cur2 = pair_idx_buf[k].second;

        if (cur1 == prev1 + 1 && cur2 == prev2 + 1) {
            prev1 = cur1;
            prev2 = cur2;
        } else {
            LocalHomTable::Segment seg;
            seg.start1 = segStart1;
            seg.end1 = prev1;
            seg.start2 = segStart2;
            seg.end2 = prev2;
            seg.optScore = normOptScore;
            seg.importance = 0.0f;
            result.segments.push_back(seg);

            segStart1 = cur1;
            segStart2 = cur2;
            prev1 = cur1;
            prev2 = cur2;
        }
    }

    LocalHomTable::Segment lastSeg;
    lastSeg.start1 = segStart1;
    lastSeg.end1 = prev1;
    lastSeg.start2 = segStart2;
    lastSeg.end2 = prev2;
    lastSeg.optScore = normOptScore;
    lastSeg.importance = 0.0f;
    result.segments.push_back(lastSeg);

    return result;
}

AlignmentResult Aligner::align_affine_local (const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    const int gOpen = static_cast<int>(params.gapOpen);
    const int gExt  = static_cast<int>(params.gapExtend);
    const int MIN_INF = -100000000; 

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    AlignmentResult result;
    if (totalRef == 0 || totalQry == 0) {
        return result;
    }

    const size_t qrySize = static_cast<size_t>(totalQry + 1);
    if (score_0.size() < qrySize) score_0.resize(qrySize);
    if (score_1.size() < qrySize) score_1.resize(qrySize);
    if (score_2.size() < qrySize) score_2.resize(qrySize);
    if (E_buf.size() < qrySize) E_buf.resize(qrySize);
    if (qry_idx_buf.size() < static_cast<size_t>(totalQry)) qry_idx_buf.resize(totalQry);

    const size_t pitch = qrySize;
    const size_t tbSize = static_cast<size_t>(totalRef + 1) * pitch;
    if (tb_buf.size() < tbSize) tb_buf.resize(tbSize);

    // Fast O(totalQry) boundary initialization (no memset on tb_buf!)
    std::fill_n(score_0.data(), qrySize, 0);
    std::fill_n(E_buf.data(), qrySize, MIN_INF);

    int* scoreRows[3] = { score_0.data(), score_1.data(), score_2.data() };
    int* E = E_buf.data();
    uint8_t* traceback = tb_buf.data();

    // Precompute query residue indices
    for (int j = 0; j < totalQry; ++j) {
        qry_idx_buf[j] = letterIdx(type, query[j]);
    }
    const int* qryIndices = qry_idx_buf.data();

    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    // Fill DP table
    for (int i = 1; i <= totalRef; ++i) {
        int curr = i % 3;
        int prev = (i - 1 + 3) % 3;
        int* score_curr = scoreRows[curr];
        const int* score_prev = scoreRows[prev];
        score_curr[0] = 0;

        int refIdx = letterIdx(type, reference[i - 1]);
        const float* scoreRow = params.scoringMatrix[refIdx];
        size_t rowOffset = static_cast<size_t>(i) * pitch;

        int F_val = MIN_INF;
        for (int j = 1; j <= totalQry; ++j) {
            uint8_t tb = 0;

            // Up (E state - delete in query)
            int e_open = score_prev[j] + gOpen; 
            int e_ext  = E[j] + gExt;
            if (e_open >= e_ext) {
                E[j] = e_open;
            } else {
                E[j] = e_ext;
                tb |= 0x04;
            }

            // Left (F state - insert in query)
            int f_open = score_curr[j - 1] + gOpen;
            int f_ext  = F_val + gExt;
            if (f_open >= f_ext) {
                F_val = f_open;
            } else {
                F_val = f_ext;
                tb |= 0x08;
            }

            // Diagonal (Match / Mismatch)
            int match = static_cast<int>(scoreRow[qryIndices[j - 1]]);
            int diag = score_prev[j - 1] + match;

            // Local alignment score floor is 0
            int best = 0; 
            uint8_t h_src = 0;

            if (diag > best) { best = diag; h_src = 1; }
            if (E[j] > best) { best = E[j]; h_src = 2; }
            if (F_val > best) { best = F_val; h_src = 3; }

            score_curr[j] = best;
            tb |= h_src;
            traceback[rowOffset + j] = tb;

            if (best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    result.score = max_score;
    if (max_score <= 0) {
        return result;
    }

    // Traceback phase
    int i = max_i;
    int j = max_j;
    int currentState = 0;
    pair_idx_buf.clear();

    while (i > 0 && j > 0) {
        uint8_t tb = traceback[static_cast<size_t>(i) * pitch + j];
        
        if (currentState == 0) { 
            uint8_t h_src = tb & 0x03;
            if (h_src == 1) { // Diagonal
                pair_idx_buf.push_back({i - 1, j - 1});
                --i; --j;
            } else if (h_src == 2) { // Up (E)
                currentState = 1; 
            } else if (h_src == 3) { // Left (F)
                currentState = 2;
            } else { 
                break; // Local alignment ends when h_src is 0
            }
        } 
        else if (currentState == 1) { // Up (E state)
            bool e_from_e = (tb & 0x04) != 0; 
            --i; 
            if (!e_from_e) currentState = 0;
        } 
        else if (currentState == 2) { // Left (F state)
            bool f_from_f = (tb & 0x08) != 0; 
            --j; 
            if (!f_from_f) currentState = 0;
        }
    }

    std::reverse(pair_idx_buf.begin(), pair_idx_buf.end());

    int identicalPairs = 0;
    result.alignedPairs.clear();
    result.alignedPairs.reserve(pair_idx_buf.size());
    for (const auto& alignedPair : pair_idx_buf) {
        result.alignedPairs.push_back({alignedPair.first, alignedPair.second});
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    
    if (!result.alignedPairs.empty()) {
        result.identity = static_cast<float>(identicalPairs) / static_cast<float>(result.alignedPairs.size());
    }
    
    return result;
}

LocalHomTable::PairResult Aligner::align_affine_local_segments(const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    const int gOpen = static_cast<int>(params.localGapOpen);
    const int gExt  = static_cast<int>(params.localGapExtend);
    const int MIN_INF = -100000000; 

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    LocalHomTable::PairResult result;
    if (totalRef == 0 || totalQry == 0) {
        return result;
    }

    const size_t qrySize = static_cast<size_t>(totalQry + 1);
    if (score_0.size() < qrySize) score_0.resize(qrySize);
    if (score_1.size() < qrySize) score_1.resize(qrySize);
    if (score_2.size() < qrySize) score_2.resize(qrySize);
    if (E_buf.size() < qrySize) E_buf.resize(qrySize);
    if (qry_idx_buf.size() < static_cast<size_t>(totalQry)) qry_idx_buf.resize(totalQry);

    const size_t pitch = qrySize;
    const size_t tbSize = static_cast<size_t>(totalRef + 1) * pitch;
    if (tb_buf.size() < tbSize) tb_buf.resize(tbSize);

    // Fast O(totalQry) boundary initialization (no memset on tb_buf!)
    std::fill_n(score_0.data(), qrySize, 0);
    std::fill_n(E_buf.data(), qrySize, MIN_INF);

    int* scoreRows[3] = { score_0.data(), score_1.data(), score_2.data() };
    int* E = E_buf.data();
    uint8_t* traceback = tb_buf.data();

    // Precompute query residue indices
    for (int j = 0; j < totalQry; ++j) {
        qry_idx_buf[j] = letterIdx(type, query[j]);
    }
    const int* qryIndices = qry_idx_buf.data();

    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    for (int i = 1; i <= totalRef; ++i) {
        int curr = i % 3;
        int prev = (i - 1 + 3) % 3;
        int* score_curr = scoreRows[curr];
        const int* score_prev = scoreRows[prev];
        score_curr[0] = 0;

        int refIdx = letterIdx(type, reference[i - 1]);
        const float* scoreRow = params.scoringMatrix[refIdx];
        size_t rowOffset = static_cast<size_t>(i) * pitch;

        int F_val = MIN_INF;
        for (int j = 1; j <= totalQry; ++j) {
            uint8_t tb = 0;

            int e_open = score_prev[j] + gOpen; 
            int e_ext  = E[j] + gExt;
            if (e_open >= e_ext) {
                E[j] = e_open;
            } else {
                E[j] = e_ext;
                tb |= 0x04;
            }

            int f_open = score_curr[j - 1] + gOpen;
            int f_ext  = F_val + gExt;
            if (f_open >= f_ext) {
                F_val = f_open;
            } else {
                F_val = f_ext;
                tb |= 0x08;
            }

            int match = static_cast<int>(scoreRow[qryIndices[j - 1]]);
            int diag = score_prev[j - 1] + match;

            int best = 0; 
            uint8_t h_src = 0;

            if (diag > best) { best = diag; h_src = 1; }
            if (E[j] > best) { best = E[j]; h_src = 2; }
            if (F_val > best) { best = F_val; h_src = 3; }

            score_curr[j] = best;
            tb |= h_src;
            traceback[rowOffset + j] = tb;

            if (best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    result.score = max_score;
    if (max_score <= 0) {
        return result;
    }

    // 1. Traceback Primary Local Alignment
    int i = max_i;
    int j = max_j;
    int currentState = 0;
    pair_idx_buf.clear();

    while (i > 0 && j > 0) {
        uint8_t tb = traceback[static_cast<size_t>(i) * pitch + j];
        
        if (currentState == 0) { 
            uint8_t h_src = tb & 0x03;
            if (h_src == 1) { // Diagonal
                pair_idx_buf.push_back({i - 1, j - 1});
                --i; --j;
            } else if (h_src == 2) { // Up (E)
                currentState = 1; 
            } else if (h_src == 3) { // Left (F)
                currentState = 2;
            } else { 
                break; // Local alignment ends when h_src is 0
            }
        } 
        else if (currentState == 1) { // Up (E)
            bool e_from_e = (tb & 0x04) != 0; 
            --i; 
            if (!e_from_e) currentState = 0;
        } 
        else if (currentState == 2) { // Left (F)
            bool f_from_f = (tb & 0x08) != 0; 
            --j; 
            if (!f_from_f) currentState = 0;
        }
    }

    if (pair_idx_buf.empty()) {
        return result;
    }

    std::reverse(pair_idx_buf.begin(), pair_idx_buf.end());

    // Compute percent identity for primary alignment
    int identicalPairs = 0;
    for (const auto& alignedPair : pair_idx_buf) {
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    result.identity = static_cast<float>(identicalPairs) / static_cast<float>(pair_idx_buf.size());

    // Helper lambda to append segments from aligned coordinate pairs
    auto appendSegmentsFromPairs = [&](const std::vector<std::pair<int, int>>& pairs, float optScore) {
        if (pairs.empty()) return;
        int segStart1 = pairs[0].first;
        int segStart2 = pairs[0].second;
        int prev1 = segStart1;
        int prev2 = segStart2;

        for (size_t k = 1; k < pairs.size(); ++k) {
            int cur1 = pairs[k].first;
            int cur2 = pairs[k].second;

            if (cur1 == prev1 + 1 && cur2 == prev2 + 1) {
                prev1 = cur1;
                prev2 = cur2;
            } else {
                LocalHomTable::Segment seg;
                seg.start1 = segStart1;
                seg.end1 = prev1;
                seg.start2 = segStart2;
                seg.end2 = prev2;
                seg.optScore = optScore;
                seg.importance = 0.0f;
                result.segments.push_back(seg);

                segStart1 = cur1;
                segStart2 = cur2;
                prev1 = cur1;
                prev2 = cur2;
            }
        }

        LocalHomTable::Segment lastSeg;
        lastSeg.start1 = segStart1;
        lastSeg.end1 = prev1;
        lastSeg.start2 = segStart2;
        lastSeg.end2 = prev2;
        lastSeg.optScore = optScore;
        lastSeg.importance = 0.0f;
        result.segments.push_back(lastSeg);
    };

    // Calculate sum of substitution scores for aligned ungapped residue pairs (MAFFT-compatible)
    int sum_subst_score = 0;
    for (const auto& p : pair_idx_buf) {
        int refIdx = letterIdx(type, reference[p.first]);
        int qryIdx = letterIdx(type, query[p.second]);
        sum_subst_score += static_cast<int>(params.scoringMatrix[refIdx][qryIdx]);
    }

    // Normalize primary pairwise score by overlap length (ungapped average match score)
    const float primaryNormOpt = (!pair_idx_buf.empty())
        ? (static_cast<float>(sum_subst_score) / static_cast<float>(pair_idx_buf.size()))
        : 0.0f;
    appendSegmentsFromPairs(pair_idx_buf, primaryNormOpt);

    return result;
}

AlignmentResult Aligner::align_linear_local (const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    // Linear gap only uses gapExtend
    const int gExt = static_cast<int>(params.gapExtend);

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    // Allocate full matrix
    std::vector<std::vector<int>> score(totalRef + 1, std::vector<int>(totalQry + 1, 0));
    std::vector<std::vector<uint8_t>> traceback(totalRef + 1, std::vector<uint8_t>(totalQry + 1, 0));
    
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    // Fill DP table
    for (int i = 1; i <= totalRef; ++i) {
        for (int j = 1; j <= totalQry; ++j) {
            
            // Calculate scores for three directions (Linear Gap)
            int diag = score[i - 1][j - 1] + residueScore(reference[i - 1], query[j - 1], type, params);
            int up   = score[i - 1][j] + gExt; // Delete in query
            int left = score[i][j - 1] + gExt; // Insert in query

            // Local alignment score floor is 0
            int best = 0; 
            uint8_t h_src = 0; // 0 indicates termination

            if (diag > best) { best = diag; h_src = 1; }
            if (up > best)   { best = up; h_src = 2; }
            if (left > best) { best = left; h_src = 3; }

            score[i][j] = best;
            traceback[i][j] = h_src;

            // Track global max score as traceback starting point
            if (best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback stage
    int i = max_i;
    int j = max_j;
    std::vector<std::pair<int, int>> aligned_pairs;

    // Linear gap traceback without state machine
    while (i > 0 && j > 0) {
        // Stop when local alignment score reaches 0
        if (score[i][j] <= 0) break; 
        
        uint8_t tb = traceback[i][j];
        
        if (tb == 1) { // From Diagonal (Match/Mismatch)
            aligned_pairs.push_back({i - 1, j - 1});
            --i; --j;
        } else if (tb == 2) { // From Up (Gap in query)
            --i;
        } else if (tb == 3) { // From Left (Gap in reference)
            --j;
        } else { 
            break; // Stop at 0 (safety guard)
        }
    }

    // Traceback collected pairs backwards, reverse to forward order
    std::reverse(aligned_pairs.begin(), aligned_pairs.end());

    // Build alignment result
    AlignmentResult result;
    result.score = max_score;

    int identicalPairs = 0;
    for (const auto& alignedPair : aligned_pairs) {
        result.alignedPairs.push_back({alignedPair.first, alignedPair.second});
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    
    if (!result.alignedPairs.empty()) {
        result.identity = static_cast<float>(identicalPairs) / static_cast<float>(result.alignedPairs.size());
    }
    
    return result;
}

/*
// Tiling Local (Not good)
AlignmentResult Aligner::align_affine_local (const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    const int gOpen = static_cast<int>(params.gapOpen);
    const int gExt  = static_cast<int>(params.gapExtend);

    const int MIN_INF = -100000000; 

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    // 1. Tiling Parameters
    const int T = 1000;
    const int O = 100;

    int ref_idx = 0;
    int qry_idx = 0;

    int max_ref_tile = std::min(T, totalRef);
    int max_qry_tile = std::min(T, totalQry);

    std::vector<std::vector<int>> score(max_ref_tile + 1, std::vector<int>(max_qry_tile + 1, 0));
    std::vector<std::vector<uint8_t>> traceback(max_ref_tile + 1, std::vector<uint8_t>(max_qry_tile + 1, 0));
    std::vector<int> E(max_qry_tile + 1, MIN_INF);

    std::vector<std::pair<int, int>> global_pairs;
    
    int last_max_score = 0;
    int last_max_state = 0;
    int final_alignment_score = 0;

    while (ref_idx < totalRef && qry_idx < totalQry) {
        int refLen = std::min(T, totalRef - ref_idx);
        int qryLen = std::min(T, totalQry - qry_idx);

        std::fill(E.begin(), E.end(), MIN_INF);

        if (ref_idx == 0 && qry_idx == 0) {
            for (int i = 0; i <= refLen; ++i) { score[i][0] = 0; traceback[i][0] = 0; }
            for (int j = 0; j <= qryLen; ++j) { score[0][j] = 0; traceback[0][j] = 0; }
        } else {
            score[0][0] = last_max_score;
            for (int i = 1; i <= refLen; ++i) {
                int penalty = (last_max_state == 1 && i == 1) ? gExt : (i == 1 ? gOpen : gExt);
                score[i][0] = score[i-1][0] + penalty;
                uint8_t tb = 2; // UP
                if (i > 1 || last_max_state == 1) tb |= 0x04; // affine gap extension
                traceback[i][0] = tb;
            }
            for (int j = 1; j <= qryLen; ++j) {
                int penalty = (last_max_state == 2 && j == 1) ? gExt : (j == 1 ? gOpen : gExt);
                score[0][j] = score[0][j-1] + penalty;
                uint8_t tb = 3; // LEFT
                if (j > 1 || last_max_state == 2) tb |= 0x08; // affine gap extension
                traceback[0][j] = tb;
            }
        }

        for (int i = 1; i <= refLen; ++i) {
            int F_val = MIN_INF;
            for (int j = 1; j <= qryLen; ++j) {
                uint8_t tb = 0;

                int e_open = score[i - 1][j] + gOpen; 
                int e_ext  = E[j] + gExt;
                if (e_open >= e_ext) {
                    E[j] = e_open;
                } else {
                    E[j] = e_ext;
                    tb |= 0x04;
                }

                int f_open = score[i][j - 1] + gOpen;
                int f_ext  = F_val + gExt;
                if (f_open >= f_ext) {
                    F_val = f_open;
                } else {
                    F_val = f_ext;
                    tb |= 0x08;
                }

                int diag = score[i - 1][j - 1] + residueScore(reference[ref_idx + i - 1], query[qry_idx + j - 1], type, params);

                int best = 0; 
                uint8_t h_src = 0;

                if (diag > best) { best = diag; h_src = 1; }
                if (E[j] > best) { best = E[j]; h_src = 2; }
                if (F_val > best) { best = F_val; h_src = 3; }

                score[i][j] = best;
                tb |= h_src;
                traceback[i][j] = tb;
            }
        }

        bool is_last_tile_ref = (ref_idx + refLen == totalRef);
        bool is_last_tile_qry = (qry_idx + qryLen == totalQry);
        bool is_last_tile = is_last_tile_ref || is_last_tile_qry;

        int best_i = refLen;
        int best_j = qryLen;
        int max_score = MIN_INF;
        uint8_t best_tb = 0;

        if (is_last_tile) {
            for (int i = 0; i <= refLen; ++i) {
                for (int j = 0; j <= qryLen; ++j) {
                    if (score[i][j] > max_score) {
                        max_score = score[i][j];
                        best_i = i; best_j = j;
                        best_tb = traceback[i][j];
                    }
                }
            }
            final_alignment_score = max_score;
        } else {
            int start_i = std::max(1, refLen - O);
            int start_j = std::max(1, qryLen - O);
            
            // 1. Search for the L-shaped "bottom horizontal strip region" (Bottom Region)
            for (int i = start_i; i <= refLen; ++i) {
                for (int j = 1; j <= qryLen; ++j) {
                    if (score[i][j] > max_score) {
                        max_score = score[i][j];
                        best_i = i; best_j = j;
                        best_tb = traceback[i][j];
                    }
                }
            }
            
            // 2. Search for the L-shaped "right vertical strip region" (Right Region)
            for (int i = 1; i < start_i; ++i) {
                for (int j = start_j; j <= qryLen; ++j) {
                    if (score[i][j] > max_score) {
                        max_score = score[i][j];
                        best_i = i; best_j = j;
                        best_tb = traceback[i][j];
                    }
                }
            }
        }

        int i = best_i;
        int j = best_j;
        int currentState = 0;
        std::vector<std::pair<int, int>> local_pairs;

        while (i > 0 || j > 0) {
            if (score[i][j] <= 0) break; 
            if (ref_idx == 0 && qry_idx == 0 && score[i][j] == 0 && traceback[i][j] == 0) break;
            
            if (i == 0 && j == 0) break;
            uint8_t tb = traceback[i][j];
            if (currentState == 0) { 
                uint8_t h_src = tb & 0x03;
                if (h_src == 1) {
                    local_pairs.push_back({ref_idx + i - 1, qry_idx + j - 1});
                    --i; --j;
                } else if (h_src == 2) { 
                    currentState = 1; 
                } else if (h_src == 3) { 
                    currentState = 2;
                } else { 
                    if (i > 0 && j == 0) currentState = 1;
                    else if (j > 0 && i == 0) currentState = 2;
                    else break;
                }
            } 
            else if (currentState == 1) { 
                bool e_from_e = (tb & 0x04) != 0; 
                --i; 
                if (!e_from_e) currentState = 0; 
            } 
            else if (currentState == 2) { 
                bool f_from_f = (tb & 0x08) != 0; 
                --j; 
                if (!f_from_f) currentState = 0; 
            }
        }

        std::reverse(local_pairs.begin(), local_pairs.end());
        global_pairs.insert(global_pairs.end(), local_pairs.begin(), local_pairs.end());

        if (is_last_tile) break;

        ref_idx += best_i;
        qry_idx += best_j;
        last_max_score = max_score;

        uint8_t h_src = best_tb & 0x03;
        if (h_src == 1) last_max_state = 0;
        else if (h_src == 2) last_max_state = 1;
        else if (h_src == 3) last_max_state = 2;
        else last_max_state = 0;
    }

    AlignmentResult result;
    result.score = final_alignment_score;

    int identicalPairs = 0;
    for (const auto& alignedPair : global_pairs) {
        result.alignedPairs.push_back({alignedPair.first, alignedPair.second});
        const char ref_c = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.first])));
        const char qry_c = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.second])));
        if (ref_c == qry_c) ++identicalPairs;
    }
    
    if (!result.alignedPairs.empty()) {
        result.identity = static_cast<float>(identicalPairs) / static_cast<float>(result.alignedPairs.size());
    }
    
    return result;
}
*/

} // namespace accurate
} // namespace msa



// Helper function to compute the similarity score between two profile columns 
// mimicking the numerator/denominator logic from the provided code.
inline float scoreProfileOptimized(
    const std::vector<float>& refCol,
    const std::vector<float>& transQryCol,
    int alphabetSize)
{
    float score = 0.0f;
    for (int l = 0; l < alphabetSize; ++l) {
        score += refCol[l] * transQryCol[l];
    }
    return score;
}

std::vector<int8_t> msa::alignProfile_semi_global(
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp, 
    const std::vector<std::vector<float>>& gapEx, 
    const std::pair<float, float>& num,           
    msa::Params& param,
    const std::vector<std::vector<float>>* consistencyTable,
    float consistencyWeight)
{
    const int N = static_cast<int>(reference.size());
    const int M = static_cast<int>(query.size());
    const float MIN_INF = -1e9f;

    bool isProtein = (reference[0].size() != 6);
    int alphabetSize = isProtein ? 22 : 6;
    int gapIdx = alphabetSize - 1; 
    float denominator = num.first * num.second;

    // =================================================================
    // Optimization core: Pre-calculate Transformed Query Profile (O(M * K^2) instead of O(N * M * K^2))
    // Perform denominator division here as well to save division overhead within DP
    // =================================================================
    std::vector<std::vector<float>> transQry(M, std::vector<float>(alphabetSize, 0.0f));
    for (int j = 0; j < M; ++j) {
        for (int l = 0; l < alphabetSize; ++l) {
            float sum = 0.0f;
            for (int m = 0; m < alphabetSize; ++m) {
                if (m != gapIdx && l != gapIdx) {
                    sum += query[j][m] * param.scoringMatrix[m][l];
                }
            }
            transQry[j][l] = sum / denominator;
        }
    }
    // =================================================================
    // Memory optimization: Use 1D vector to flatten 2D matrix, ensuring contiguous memory and eliminating allocation overhead
    // For floating-point scores, we only need to keep the previous (prev) and current (curr) rows
    // =================================================================
    std::vector<float> M_prev(M + 1, 0.0f), M_curr(M + 1, MIN_INF);
    std::vector<float> I_prev(M + 1, 0.0f), I_curr(M + 1, MIN_INF);
    std::vector<float> D_prev(M + 1, 0.0f), D_curr(M + 1, MIN_INF);

    // Traceback 
    const uint8_t STATE_M = 0, STATE_I = 1, STATE_D = 2;
    int totalCells = (N + 1) * (M + 1); // Traceback needs to record the full path, using a 1D array to calculate Index (i * (M+1) + j)
    std::vector<uint8_t> tb_M(totalCells, 0);
    std::vector<uint8_t> tb_I(totalCells, 0);
    std::vector<uint8_t> tb_D(totalCells, 0);


    M_prev[0] = 0.0f;
    for (int j = 1; j <= M; ++j) {
        M_prev[j] = 0.0f;
        I_prev[j] = 0.0f;
        D_prev[j] = MIN_INF; 
    }

    for (int i = 1; i <= N; ++i) {
        M_curr[0] = 0.0f;
        D_curr[0] = 0.0f;
        I_curr[0] = MIN_INF;

        float pos_gapOpen_ref   = gapOp[0][i - 1];
        float pos_gapExtend_ref = gapEx[0][i - 1];

        for (int j = 1; j <= M; ++j) {
            float pos_gapOpen_qry   = gapOp[1][j - 1];
            float pos_gapExtend_qry = gapEx[1][j - 1];
            
            int idx = i * (M + 1) + j;

            // -- Calculate I_mat --
            float i_open = M_curr[j - 1] + pos_gapOpen_qry;
            float i_ext  = I_curr[j - 1] + pos_gapExtend_qry;
            if (i_open >= i_ext) { I_curr[j] = i_open; tb_I[idx] = STATE_M; } 
            else                 { I_curr[j] = i_ext;  tb_I[idx] = STATE_I; }

            // -- Calculate D_mat --
            float d_open = M_prev[j] + pos_gapOpen_ref;
            float d_ext  = D_prev[j] + pos_gapExtend_ref;
            if (d_open >= d_ext) { D_curr[j] = d_open; tb_D[idx] = STATE_M; } 
            else                 { D_curr[j] = d_ext;  tb_D[idx] = STATE_D; }

            // -- Calculate M_mat --
            // Changed here to call optimized O(K) dot product
            // float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize);
            // Calculate Consistency Bonus
            float consistencyBonus = 0.0f;
            if (consistencyTable != nullptr &&
                i - 1 < static_cast<int32_t>(consistencyTable->size()) &&
                j - 1 < static_cast<int32_t>((*consistencyTable)[i - 1].size())) {
                consistencyBonus = consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            // M_mat = Dot Product + Consistency Bonus
            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;

            // if (consistencyBonus > 0 && (i == j)) std::cout << match_score  << " -> " << consistencyBonus << std::endl;
            
            
            float m_from_m = M_prev[j - 1];
            float m_from_i = I_prev[j - 1];
            float m_from_d = D_prev[j - 1];

            float max_m = m_from_m; 
            tb_M[idx] = STATE_M;
            
            if (m_from_i > max_m) { max_m = m_from_i; tb_M[idx] = STATE_I; }
            if (m_from_d > max_m) { max_m = m_from_d; tb_M[idx] = STATE_D; }
            
            M_curr[j] = match_score + max_m;
        }
        
        M_prev = M_curr;
        I_prev = I_curr;
        D_prev = D_curr;
    }

    // =================================================================
    // Traceback logic (finding the maximum value in the last row or column)
    // =================================================================
    float best_score = MIN_INF;
    int best_i = N, best_j = M;
    uint8_t best_state = STATE_M;

    
    std::vector<int8_t> path;
    int curr_i = best_i;
    int curr_j = best_j;

    while (curr_i < N) { path.push_back(2); curr_i++; } 
    while (curr_j < M) { path.push_back(1); curr_j++; } 
    
    curr_i = N; // Theoretically, best_i is already N
    curr_j = best_j;
    uint8_t state = best_state;

    while (curr_i > 0 || curr_j > 0) {
        if (curr_i == 0) {
            path.push_back(1);
            curr_j--;
        } 
        else if (curr_j == 0) {
            path.push_back(2);
            curr_i--;
        } 
        else {
            int idx = curr_i * (M + 1) + curr_j;
            if (state == STATE_M) {
                path.push_back(0);
                state = tb_M[idx];
                curr_i--;
                curr_j--;
            } 
            else if (state == STATE_I) {
                path.push_back(1);
                state = tb_I[idx];
                curr_j--;
            } 
            else if (state == STATE_D) {
                path.push_back(2);
                state = tb_D[idx];
                curr_i--;
            }
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}


std::vector<int8_t> msa::alignProfile_global(
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp, 
    const std::vector<std::vector<float>>& gapEx, 
    const std::pair<float, float>& num,           
    msa::Params& param,
    const std::vector<std::vector<float>>* consistencyTable, 
    float consistencyWeight)
{
    const int N = static_cast<int>(reference.size());
    const int M = static_cast<int>(query.size());
    const float MIN_INF = -1e9f;

    bool isProtein = (reference[0].size() != 6);
    int alphabetSize = isProtein ? 22 : 6;
    int gapIdx = alphabetSize - 1; 
    float denominator = num.first * num.second;

    // =================================================================
    // Pre-calculate Transformed Query Profile
    // =================================================================
    std::vector<std::vector<float>> transQry(M, std::vector<float>(alphabetSize, 0.0f));
    for (int j = 0; j < M; ++j) {
        for (int l = 0; l < alphabetSize; ++l) {
            float sum = 0.0f;
            for (int m = 0; m < alphabetSize; ++m) {
                if (m != gapIdx && l != gapIdx) {
                    sum += query[j][m] * param.scoringMatrix[m][l];
                }
            }
            transQry[j][l] = sum / denominator;
        }
    }

    // =================================================================
    // DP Matrices (1D vectors for memory optimization)
    // =================================================================
    std::vector<float> M_prev(M + 1, MIN_INF), M_curr(M + 1, MIN_INF);
    std::vector<float> I_prev(M + 1, MIN_INF), I_curr(M + 1, MIN_INF);
    std::vector<float> D_prev(M + 1, MIN_INF), D_curr(M + 1, MIN_INF);

    const uint8_t STATE_M = 0, STATE_I = 1, STATE_D = 2;
    size_t totalCells = static_cast<size_t>(N + 1) * (M + 1); 
    std::vector<uint8_t> tb(totalCells, 0);

    // =================================================================
    // 1. Global Alignment Initialization (Top Row)
    // Applies linear penalty for leading query gaps (Insertion state)
    // =================================================================
    M_prev[0] = 0.0f;
    for (int j = 1; j <= M; ++j) {
        M_prev[j] = MIN_INF;
        I_prev[j] = j * param.gapTerminal; // Linear leading gap
        D_prev[j] = MIN_INF;
    }

    // =================================================================
    // DP Loop
    // =================================================================
    for (int i = 1; i <= N; ++i) {
        // 1. Global Alignment Initialization (Left Column)
        // Applies linear penalty for leading reference gaps (Deletion state)
        M_curr[0] = MIN_INF;
        D_curr[0] = i * param.gapTerminal; // Linear leading gap
        I_curr[0] = MIN_INF;

        for (int j = 1; j <= M; ++j) {
            size_t idx = static_cast<size_t>(i) * (M + 1) + j;

            float nongap_ref = 1.0f - (reference[i-1][gapIdx] / num.first);
            float nongap_qry = 1.0f - (query[j-1][gapIdx] / num.second);


            // 2. Dynamic Border Penalties for Trailing Gaps
            // If at the last row of Reference (i==N), horizontal movement is a trailing insertion in Reference
            // If at the last col of Query (j==M), vertical movement is a trailing deletion in Query
            float pos_gapOpen_qry   = (i == N) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[1][j - 1]));
            float pos_gapExtend_qry = (i == N) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[1][j - 1]));

            float pos_gapOpen_ref   = (j == M) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[0][i - 1]));
            float pos_gapExtend_ref = (j == M) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[0][i - 1]));
            
            // -- Calculate I_mat (Insertion in Reference: Query residue j, Reference gap i) --
            float i_open = M_curr[j - 1] + (pos_gapOpen_qry * nongap_qry);
            float i_ext  = I_curr[j - 1] + (pos_gapExtend_qry * nongap_qry);
            uint8_t tb_i_val = 0;
            if (i_open >= i_ext) { I_curr[j] = i_open; } 
            else                 { I_curr[j] = i_ext;  tb_i_val = 0x04; }

            // -- Calculate D_mat (Deletion in Query: Reference residue i, Query gap j) --
            float d_open = M_prev[j] + (pos_gapOpen_ref * nongap_ref);
            float d_ext  = D_prev[j] + (pos_gapExtend_ref * nongap_ref);
            uint8_t tb_d_val = 0;
            if (d_open >= d_ext) { D_curr[j] = d_open; } 
            else                 { D_curr[j] = d_ext;  tb_d_val = 0x08; }

            // -- Calculate M_mat --
            float consistencyBonus = 0.0f;
            if (consistencyTable != nullptr &&
                i - 1 < static_cast<int32_t>(consistencyTable->size()) &&
                j - 1 < static_cast<int32_t>((*consistencyTable)[i - 1].size())) {
                consistencyBonus = consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;
            
            float m_from_m = M_prev[j - 1];
            float m_from_i = I_prev[j - 1];
            float m_from_d = D_prev[j - 1];

            float max_m = m_from_m; 
            uint8_t tb_m_val = STATE_M;
            
            if (m_from_i > max_m) { max_m = m_from_i; tb_m_val = STATE_I; }
            if (m_from_d > max_m) { max_m = m_from_d; tb_m_val = STATE_D; }
            
            M_curr[j] = match_score + max_m;
            tb[idx] = tb_m_val | tb_i_val | tb_d_val;
        }
        
        M_prev = M_curr;
        I_prev = I_curr;
        D_prev = D_curr;
    }

    // =================================================================
    // 3. Global Traceback strictly starts from the bottom-right corner (N, M)
    // =================================================================
    // Note: After the loop finishes, M_prev holds the values of the N-th row
    float best_score = M_prev[M]; 
    uint8_t state = STATE_M;

    if (I_prev[M] > best_score) { best_score = I_prev[M]; state = STATE_I; }
    if (D_prev[M] > best_score) { best_score = D_prev[M]; state = STATE_D; }

    std::vector<int8_t> path;
    int curr_i = N;
    int curr_j = M;

    // Removed the "while (curr_i < N)" padding logic since Global Alignment naturally covers the whole sequence
    while (curr_i > 0 || curr_j > 0) {
        if (curr_i == 0) {
            // Hit the top border, forced to track left (Insertion)
            path.push_back(1);
            curr_j--;
        } 
        else if (curr_j == 0) {
            // Hit the left border, forced to track up (Deletion)
            path.push_back(2);
            curr_i--;
        } 
        else {
            size_t idx = static_cast<size_t>(curr_i) * (M + 1) + curr_j;
            uint8_t cell = tb[idx];
            if (state == STATE_M) {
                path.push_back(0);
                state = cell & 0x03;
                curr_i--;
                curr_j--;
            } 
            else if (state == STATE_I) {
                path.push_back(1);
                state = (cell & 0x04) ? STATE_I : STATE_M;
                curr_j--;
            } 
            else if (state == STATE_D) {
                path.push_back(2);
                state = (cell & 0x08) ? STATE_D : STATE_M;
                curr_i--;
            }
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<int8_t> msa::alignProfile_global_banded(
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp, 
    const std::vector<std::vector<float>>& gapEx, 
    const std::pair<float, float>& num,           
    msa::Params& param,
    int bandWidth,
    const std::vector<std::vector<float>>* consistencyTable, 
    float consistencyWeight)
{
    const int N = static_cast<int>(reference.size());
    const int M = static_cast<int>(query.size());
    const float MIN_INF = -1e9f;

    if (N == 0) {
        return std::vector<int8_t>(M, 1);
    }
    if (M == 0) {
        return std::vector<int8_t>(N, 2);
    }

    // Fall back to full global alignment if band covers the entire matrix
    if (bandWidth >= std::min(N, M)) {
        return alignProfile_global(reference, query, gapOp, gapEx, num, param, consistencyTable, consistencyWeight);
    }

    bool isProtein = (reference[0].size() != 6);
    int alphabetSize = isProtein ? 22 : 6;
    int gapIdx = alphabetSize - 1; 
    float denominator = num.first * num.second;

    // =================================================================
    // Pre-calculate Transformed Query Profile
    // =================================================================
    std::vector<std::vector<float>> transQry(M, std::vector<float>(alphabetSize, 0.0f));
    for (int j = 0; j < M; ++j) {
        for (int l = 0; l < alphabetSize; ++l) {
            float sum = 0.0f;
            for (int m = 0; m < alphabetSize; ++m) {
                if (m != gapIdx && l != gapIdx) {
                    sum += query[j][m] * param.scoringMatrix[m][l];
                }
            }
            transQry[j][l] = sum / denominator;
        }
    }

    // =================================================================
    // DP Matrices (1D vectors for memory optimization)
    // =================================================================
    std::vector<float> M_prev(M + 1, MIN_INF), M_curr(M + 1, MIN_INF);
    std::vector<float> I_prev(M + 1, MIN_INF), I_curr(M + 1, MIN_INF);
    std::vector<float> D_prev(M + 1, MIN_INF), D_curr(M + 1, MIN_INF);

    const uint8_t STATE_M = 0, STATE_I = 1, STATE_D = 2;
    const int halfBand = std::max(1, bandWidth / 2);
    const int bandSpan = 2 * halfBand + 1;
    size_t totalCells = static_cast<size_t>(N + 1) * bandSpan; 
    std::vector<uint8_t> tb_band(totalCells, 0);

    // =================================================================
    // 1. Global Alignment Initialization (Top Row)
    // =================================================================
    M_prev[0] = 0.0f;
    int top_max = std::min(M, halfBand);
    for (int j = 1; j <= top_max; ++j) {
        I_prev[j] = j * param.gapTerminal; // Linear leading gap within the band
    }

    int prev_j_min = 1;
    int prev_j_max = top_max;

    // =================================================================
    // DP Loop
    // =================================================================
    for (int i = 1; i <= N; ++i) {
        // Slanted center tracking along the main diagonal
        int j_center = static_cast<int>((static_cast<int64_t>(i) * M) / N);
        int j_min = std::max(1, j_center - halfBand);
        int j_max = std::min(M, j_center + halfBand);
        int j_base = j_center - halfBand;

        // Sanitize M_curr, I_curr, D_curr: reset previously active cells to MIN_INF
        for (int k = prev_j_min; k <= prev_j_max; ++k) {
            M_curr[k] = MIN_INF;
            I_curr[k] = MIN_INF;
            D_curr[k] = MIN_INF;
        }
        M_curr[0] = MIN_INF;
        I_curr[0] = MIN_INF;

        // Global Alignment Initialization (Left Column)
        // Applies linear penalty for leading reference gaps if column 0 is in the band
        if (j_base <= 0) {
            D_curr[0] = i * param.gapTerminal;
        } else {
            D_curr[0] = MIN_INF;
        }

        float nongap_ref = 1.0f - (reference[i - 1][gapIdx] / num.first);

        for (int j = j_min; j <= j_max; ++j) {
            size_t idx = static_cast<size_t>(i) * bandSpan + (j - j_base);

            float nongap_qry = 1.0f - (query[j - 1][gapIdx] / num.second);

            // Dynamic Border Penalties for Trailing Gaps
            float pos_gapOpen_qry   = (i == N) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[1][j - 1]));
            float pos_gapExtend_qry = (i == N) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[1][j - 1]));

            float pos_gapOpen_ref   = (j == M) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[0][i - 1]));
            float pos_gapExtend_ref = (j == M) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[0][i - 1]));

            // -- Calculate I_mat (Insertion in Reference: Query residue j, Reference gap i) --
            float i_open = M_curr[j - 1] + (pos_gapOpen_qry * nongap_qry);
            float i_ext  = I_curr[j - 1] + (pos_gapExtend_qry * nongap_qry);
            uint8_t tb_i_val = 0;
            if (i_open >= i_ext) { I_curr[j] = i_open; } 
            else                 { I_curr[j] = i_ext;  tb_i_val = 0x04; }

            // -- Calculate D_mat (Deletion in Query: Reference residue i, Query gap j) --
            float d_open = M_prev[j] + (pos_gapOpen_ref * nongap_ref);
            float d_ext  = D_prev[j] + (pos_gapExtend_ref * nongap_ref);
            uint8_t tb_d_val = 0;
            if (d_open >= d_ext) { D_curr[j] = d_open; } 
            else                 { D_curr[j] = d_ext;  tb_d_val = 0x08; }

            // -- Calculate M_mat --
            float consistencyBonus = 0.0f;
            if (consistencyTable != nullptr &&
                i - 1 < static_cast<int32_t>(consistencyTable->size()) &&
                j - 1 < static_cast<int32_t>((*consistencyTable)[i - 1].size())) {
                consistencyBonus = consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;
            float m_from_m = M_prev[j - 1];
            float m_from_i = I_prev[j - 1];
            float m_from_d = D_prev[j - 1];

            float max_m = m_from_m; 
            uint8_t tb_m_val = STATE_M;

            if (m_from_i > max_m) { max_m = m_from_i; tb_m_val = STATE_I; }
            if (m_from_d > max_m) { max_m = m_from_d; tb_m_val = STATE_D; }

            M_curr[j] = match_score + max_m;
            tb_band[idx] = tb_m_val | tb_i_val | tb_d_val;
        }

        M_prev = M_curr;
        I_prev = I_curr;
        D_prev = D_curr;
        prev_j_min = j_min;
        prev_j_max = j_max;
    }

    // =================================================================
    // 2. Global Traceback strictly starts from the bottom-right corner (N, M)
    // =================================================================
    float best_score = M_prev[M]; 
    uint8_t state = STATE_M;

    if (I_prev[M] > best_score) { best_score = I_prev[M]; state = STATE_I; }
    if (D_prev[M] > best_score) { best_score = D_prev[M]; state = STATE_D; }

    std::vector<int8_t> path;
    int curr_i = N;
    int curr_j = M;

    while (curr_i > 0 || curr_j > 0) {
        if (curr_i == 0) {
            // Hit the top border, forced to track left (Insertion)
            path.push_back(1);
            curr_j--;
        } 
        else if (curr_j == 0) {
            // Hit the left border, forced to track up (Deletion)
            path.push_back(2);
            curr_i--;
        } 
        else {
            int j_center = static_cast<int>((static_cast<int64_t>(curr_i) * M) / N);
            int j_base = j_center - halfBand;
            int offset = curr_j - j_base;

            if (offset < 0) {
                // To the left of the band, forced to step Up (Deletion in query)
                path.push_back(2);
                state = STATE_D;
                curr_i--;
            } 
            else if (offset >= bandSpan) {
                // To the right of the band, forced to step Left (Insertion in reference)
                path.push_back(1);
                state = STATE_I;
                curr_j--;
            } 
            else {
                size_t idx = static_cast<size_t>(curr_i) * bandSpan + offset;
                uint8_t cell = tb_band[idx];
                if (state == STATE_M) {
                    path.push_back(0);
                    state = cell & 0x03;
                    curr_i--;
                    curr_j--;
                } 
                else if (state == STATE_I) {
                    path.push_back(1);
                    state = (cell & 0x04) ? STATE_I : STATE_M;
                    curr_j--;
                } 
                else if (state == STATE_D) {
                    path.push_back(2);
                    state = (cell & 0x08) ? STATE_D : STATE_M;
                    curr_i--;
                }
            }
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

std::vector<int8_t> msa::alignProfile_global_tiling(
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp, 
    const std::vector<std::vector<float>>& gapEx, 
    const std::pair<float, float>& num,           
    msa::Params& param,
    const std::vector<std::vector<float>>* consistencyTable, 
    float consistencyWeight)
{
    const int N = static_cast<int>(reference.size());
    const int M = static_cast<int>(query.size());
    const float MIN_INF = -1e9f;
    const int TILE_SIZE = 500;
    const int O = 100;

    // If both sequence lengths are <= TILE_SIZE, run standard full-matrix alignment directly
    if (N <= TILE_SIZE && M <= TILE_SIZE) {
        return alignProfile_global(reference, query, gapOp, gapEx, num, param, consistencyTable, consistencyWeight);
    }

    bool isProtein = (reference[0].size() != 6);
    int alphabetSize = isProtein ? 22 : 6;
    int gapIdx = alphabetSize - 1; 
    float denominator = num.first * num.second;

    // =================================================================
    // Pre-calculate Transformed Query Profile (One-time global computation)
    // =================================================================
    std::vector<std::vector<float>> transQry(M, std::vector<float>(alphabetSize, 0.0f));
    for (int j = 0; j < M; ++j) {
        for (int l = 0; l < alphabetSize; ++l) {
            float sum = 0.0f;
            for (int m = 0; m < alphabetSize; ++m) {
                if (m != gapIdx && l != gapIdx) {
                    sum += query[j][m] * param.scoringMatrix[m][l];
                }
            }
            transQry[j][l] = sum / denominator;
        }
    }

    const uint8_t STATE_M = 0, STATE_I = 1, STATE_D = 2;
    std::vector<int8_t> full_path;

    int ref_idx = 0;
    int qry_idx = 0;

    float last_score = 0.0f;
    uint8_t last_state = STATE_M;

    // Preallocate reusable scratchpad vectors within each tile
    int maxTileCells = (TILE_SIZE + 1) * (TILE_SIZE + 1);
    std::vector<float> M_prev(TILE_SIZE + 1, MIN_INF), M_curr(TILE_SIZE + 1, MIN_INF);
    std::vector<float> I_prev(TILE_SIZE + 1, MIN_INF), I_curr(TILE_SIZE + 1, MIN_INF);
    std::vector<float> D_prev(TILE_SIZE + 1, MIN_INF), D_curr(TILE_SIZE + 1, MIN_INF);
    std::vector<uint8_t> tb_M(maxTileCells, 0);
    std::vector<uint8_t> tb_I(maxTileCells, 0);
    std::vector<uint8_t> tb_D(maxTileCells, 0);

    while (ref_idx < N || qry_idx < M) {
        // If either sequence is exhausted, fill remaining residues as gaps
        if (ref_idx == N) {
            while (qry_idx < M) {
                full_path.push_back(1); // Insertion in ref (gap in ref)
                qry_idx++;
            }
            break;
        }
        if (qry_idx == M) {
            while (ref_idx < N) {
                full_path.push_back(2); // Deletion in qry (gap in qry)
                ref_idx++;
            }
            break;
        }

        int rem_ref = N - ref_idx;
        int rem_qry = M - qry_idx;

        int refLen, qryLen;
        if (rem_ref >= rem_qry) {
            refLen = std::min(TILE_SIZE, rem_ref);
            qryLen = std::max(1, static_cast<int>(std::round(static_cast<double>(refLen) * rem_qry / rem_ref)));
            qryLen = std::min(qryLen, rem_qry);
        } else {
            qryLen = std::min(TILE_SIZE, rem_qry);
            refLen = std::max(1, static_cast<int>(std::round(static_cast<double>(qryLen) * rem_ref / rem_qry)));
            refLen = std::min(refLen, rem_ref);
        }
        bool is_last_tile = (ref_idx + refLen == N && qry_idx + qryLen == M);

        // Reset DP vectors for current tile
        std::fill(M_prev.begin(), M_prev.begin() + qryLen + 1, MIN_INF);
        std::fill(I_prev.begin(), I_prev.begin() + qryLen + 1, MIN_INF);
        std::fill(D_prev.begin(), D_prev.begin() + qryLen + 1, MIN_INF);

        // --- 1. Initialize tile (0, 0) and Top Row (i = 0) ---
        if (ref_idx == 0 && qry_idx == 0) {
            M_prev[0] = 0.0f;
            for (int j = 1; j <= qryLen; ++j) {
                I_prev[j] = j * param.gapTerminal;
                tb_I[0 * (qryLen + 1) + j] = STATE_I;
            }
        } else {
            // Inherit boundary score and state from previous tile
            if (last_state == STATE_M)      M_prev[0] = last_score;
            else if (last_state == STATE_I) I_prev[0] = last_score;
            else if (last_state == STATE_D) D_prev[0] = last_score;
            else M_prev[0] = last_score;

            for (int j = 1; j <= qryLen; ++j) {
                int g_j = qry_idx + j;
                float nongap_qry = 1.0f - (query[g_j - 1][gapIdx] / num.second);
                float pos_gapOpen_ref   = (g_j == M) ? param.gapTerminal : param.gapOpen * (1 - gapOp[0][ref_idx]);
                float pos_gapExtend_ref = (g_j == M) ? param.gapTerminal : param.gapExtend * (1 - gapEx[0][ref_idx]);

                if (j == 1) {
                    if (last_state == STATE_I) {
                        I_prev[1] = I_prev[0] + (pos_gapExtend_ref * nongap_qry);
                        tb_I[0 * (qryLen + 1) + 1] = STATE_I;
                    } else {
                        I_prev[1] = last_score + (pos_gapOpen_ref * nongap_qry);
                        tb_I[0 * (qryLen + 1) + 1] = last_state;
                    }
                } else {
                    I_prev[j] = I_prev[j - 1] + (pos_gapExtend_ref * nongap_qry);
                    tb_I[0 * (qryLen + 1) + j] = STATE_I;
                }
            }
        }

        // --- 2. Comparator: check if Candidate (score_a, i_a, j_a) is better than (score_b, i_b, j_b) ---
        auto is_better = [&](float score_a, int i_a, int j_a, float score_b, int i_b, int j_b) -> bool {
            if (std::abs(score_a - score_b) > 1e-6f) {
                return score_a > score_b;
            }
            // (1) Tie-breaker 1: choose larger i + j (further towards bottom-right)
            int sum_a = i_a + j_a;
            int sum_b = i_b + j_b;
            if (sum_a != sum_b) return sum_a > sum_b;

            // (2) Tie-breaker 2: choose smaller |i - j| (closer to main diagonal)
            int diff_a = std::abs(i_a - j_a);
            int diff_b = std::abs(i_b - j_b);
            if (diff_a != diff_b) return diff_a < diff_b;

            // (3) Tie-breaker 3: consume whichever sequence has more remaining residues
            int rem_ref = N - ref_idx;
            int rem_qry = M - qry_idx;
            if (rem_ref > rem_qry) {
                return i_a > i_b; // Consume more reference
            } else if (rem_qry > rem_ref) {
                return j_a > j_b; // Consume more query
            }
            return i_a > i_b;
        };

        float best_score = MIN_INF;
        int best_i = refLen;
        int best_j = qryLen;
        uint8_t best_state = STATE_M;

        int o_i = std::min(O, refLen / 2);
        int o_j = std::min(O, qryLen / 2);
        int start_i = std::max(1, refLen - o_i);
        int start_j = std::max(1, qryLen - o_j);

        // --- 3. DP Table Fill Loop ---
        for (int i = 1; i <= refLen; ++i) {
            int g_i = ref_idx + i;
            float nongap_ref = 1.0f - (reference[g_i - 1][gapIdx] / num.first);
            float pos_gapOpen_qry_c0   = (g_i == N) ? param.gapTerminal : param.gapOpen * (1 - gapOp[1][qry_idx]);
            float pos_gapExtend_qry_c0 = (g_i == N) ? param.gapTerminal : param.gapExtend * (1 - gapEx[1][qry_idx]);

            // Initialize Left Column (j = 0)
            M_curr[0] = MIN_INF;
            I_curr[0] = MIN_INF;
            if (ref_idx == 0 && qry_idx == 0) {
                D_curr[0] = i * param.gapTerminal;
                tb_D[i * (qryLen + 1) + 0] = STATE_D;
            } else {
                if (i == 1) {
                    if (last_state == STATE_D) {
                        D_curr[0] = D_prev[0] + (pos_gapExtend_qry_c0 * nongap_ref);
                        tb_D[1 * (qryLen + 1) + 0] = STATE_D;
                    } else {
                        D_curr[0] = last_score + (pos_gapOpen_qry_c0 * nongap_ref);
                        tb_D[1 * (qryLen + 1) + 0] = last_state;
                    }
                } else {
                    D_curr[0] = D_prev[0] + (pos_gapExtend_qry_c0 * nongap_ref);
                    tb_D[i * (qryLen + 1) + 0] = STATE_D;
                }
            }

            for (int j = 1; j <= qryLen; ++j) {
                int g_j = qry_idx + j;
                int idx = i * (qryLen + 1) + j;

                float nongap_qry = 1.0f - (query[g_j - 1][gapIdx] / num.second);

                float pos_gapOpen_qry   = (g_i == N) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[1][g_j - 1]));
                float pos_gapExtend_qry = (g_i == N) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[1][g_j - 1]));

                float pos_gapOpen_ref   = (g_j == M) ? param.gapTerminal : 0.5f * param.gapOpen * std::max(0.0f, (1.0f - gapOp[0][g_i - 1]));
                float pos_gapExtend_ref = (g_j == M) ? param.gapTerminal : param.gapExtend * std::max(0.1f, (1.0f - gapEx[0][g_i - 1]));

                // -- Calculate I_mat (Insertion in Reference: Query residue j, Reference gap i) --
                float i_open = M_curr[j - 1] + (pos_gapOpen_qry * nongap_qry);
                float i_ext  = I_curr[j - 1] + (pos_gapExtend_qry * nongap_qry);
                if (i_open >= i_ext) { I_curr[j] = i_open; tb_I[idx] = STATE_M; } 
                else                 { I_curr[j] = i_ext;  tb_I[idx] = STATE_I; }

                // -- Calculate D_mat (Deletion in Query: Reference residue i, Query gap j) --
                float d_open = M_prev[j] + (pos_gapOpen_ref * nongap_ref);
                float d_ext  = D_prev[j] + (pos_gapExtend_ref * nongap_ref);
                if (d_open >= d_ext) { D_curr[j] = d_open; tb_D[idx] = STATE_M; } 
                else                 { D_curr[j] = d_ext;  tb_D[idx] = STATE_D; }

                // -- Calculate M_mat --
                float consistencyBonus = 0.0f;
                if (consistencyTable != nullptr &&
                    g_i - 1 < static_cast<int32_t>(consistencyTable->size()) &&
                    g_j - 1 < static_cast<int32_t>((*consistencyTable)[g_i - 1].size())) {
                    consistencyBonus = consistencyWeight * (*consistencyTable)[g_i - 1][g_j - 1];
                }

                float match_score = scoreProfileOptimized(reference[g_i - 1], transQry[g_j - 1], alphabetSize) + consistencyBonus;

                float m_from_m = M_prev[j - 1];
                float m_from_i = I_prev[j - 1];
                float m_from_d = D_prev[j - 1];

                float max_m = m_from_m; 
                tb_M[idx] = STATE_M;
                
                if (m_from_i > max_m) { max_m = m_from_i; tb_M[idx] = STATE_I; }
                if (m_from_d > max_m) { max_m = m_from_d; tb_M[idx] = STATE_D; }
                
                M_curr[j] = match_score + max_m;

                // --- 4. For non-last tile, record maximum within overlap window O at bottom-right corner ---
                if (!is_last_tile && i >= start_i && j >= start_j) {
                    float cell_score = M_curr[j];
                    uint8_t cell_state = STATE_M;
                    if (I_curr[j] > cell_score) { cell_score = I_curr[j]; cell_state = STATE_I; }
                    if (D_curr[j] > cell_score) { cell_score = D_curr[j]; cell_state = STATE_D; }

                    if (best_score == MIN_INF || is_better(cell_score, i, j, best_score, best_i, best_j)) {
                        best_score = cell_score;
                        best_i = i;
                        best_j = j;
                        best_state = cell_state;
                    }
                }
            }

            M_prev = M_curr;
            I_prev = I_curr;
            D_prev = D_curr;
        }

        // --- 5. For the last tile, strictly anchor to endpoint (refLen, qryLen) ---
        if (is_last_tile) {
            best_i = refLen;
            best_j = qryLen;
            best_score = M_prev[qryLen];
            best_state = STATE_M;
            if (I_prev[qryLen] > best_score) { best_score = I_prev[qryLen]; best_state = STATE_I; }
            if (D_prev[qryLen] > best_score) { best_score = D_prev[qryLen]; best_state = STATE_D; }
        }

        // --- 6. Tile Traceback: trace backward from (best_i, best_j) to (0, 0) ---
        std::vector<int8_t> tile_path;
        int curr_i = best_i;
        int curr_j = best_j;
        uint8_t state = best_state;

        while (curr_i > 0 || curr_j > 0) {
            if (curr_i == 0) {
                tile_path.push_back(1);
                curr_j--;
            } else if (curr_j == 0) {
                tile_path.push_back(2);
                curr_i--;
            } else {
                int idx = curr_i * (qryLen + 1) + curr_j;
                if (state == STATE_M) {
                    tile_path.push_back(0);
                    state = tb_M[idx];
                    curr_i--;
                    curr_j--;
                } else if (state == STATE_I) {
                    tile_path.push_back(1);
                    state = tb_I[idx];
                    curr_j--;
                } else if (state == STATE_D) {
                    tile_path.push_back(2);
                    state = tb_D[idx];
                    curr_i--;
                }
            }
        }

        std::reverse(tile_path.begin(), tile_path.end());
        full_path.insert(full_path.end(), tile_path.begin(), tile_path.end());

        if (is_last_tile) break;

        // Advance to next tile
        ref_idx += best_i;
        qry_idx += best_j;
        last_score = best_score;
        last_state = best_state;
    }

    return full_path;
}

std::vector<int8_t> msa::alignProfile_MAFFT(

    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp, 
    const std::vector<std::vector<float>>& gapEx, 
    const std::pair<float, float>& num,           
    msa::Params& param,
    const std::vector<std::vector<float>>* consistencyTable, 
    float consistencyWeight)
{
    const int N = static_cast<int>(reference.size());
    const int M = static_cast<int>(query.size());
    const float MIN_INF = -1e9f;

    bool isProtein = (reference[0].size() != 6);
    int alphabetSize = isProtein ? 22 : 6;
    int gapIdx = alphabetSize - 1; 
    float denominator = num.first * num.second;

    // =================================================================
    // Pre-calculate Transformed Query Profile
    // =================================================================
    std::vector<std::vector<float>> transQry(M, std::vector<float>(alphabetSize, 0.0f));
    for (int j = 0; j < M; ++j) {
        for (int l = 0; l < alphabetSize; ++l) {
            float sum = 0.0f;
            for (int m = 0; m < alphabetSize; ++m) {
                if (m != gapIdx && l != gapIdx) {
                    sum += query[j][m] * param.scoringMatrix[m][l];
                }
            }
            transQry[j][l] = sum / denominator;
        }
    }

    // =================================================================
    // DP Matrices (1D vectors for memory optimization)
    // =================================================================
    std::vector<float> M_prev(M + 1, MIN_INF), M_curr(M + 1, MIN_INF);
    std::vector<float> I_prev(M + 1, MIN_INF), I_curr(M + 1, MIN_INF);
    std::vector<float> D_prev(M + 1, MIN_INF), D_curr(M + 1, MIN_INF);

    const uint8_t STATE_M = 0, STATE_I = 1, STATE_D = 2;
    int totalCells = (N + 1) * (M + 1); 
    std::vector<uint8_t> tb_M(totalCells, 0);
    std::vector<uint8_t> tb_I(totalCells, 0);
    std::vector<uint8_t> tb_D(totalCells, 0);

    std::vector<float> gapOpen_ref (N, 0.0f);
    std::vector<float> gapEnd_ref (N, 0.0f);
    std::vector<float> gapOpen_qry (M, 0.0f);
    std::vector<float> gapEnd_qry (M, 0.0f);

    for(int i = 0; i < N; i++) {
        float nongap_freq = 1.0f - (reference[i][gapIdx] / num.first);
        gapOpen_ref[i] = 0.5f * (1.0f - gapOp[0][i]) * param.gapOpen * nongap_freq;
        gapEnd_ref[i]  = 0.5f * (1.0f - gapEx[0][i]) * param.gapOpen * nongap_freq;
    }

    for(int j = 0; j < M; j++) {
        float nongap_freq = 1.0f - (query[j][gapIdx] / num.second);
        gapOpen_qry[j] = 0.5f * (1.0f - gapOp[1][j]) * param.gapOpen * nongap_freq;
        gapEnd_qry[j]  = 0.5f * (1.0f - gapEx[1][j]) * param.gapOpen * nongap_freq;
    }

    auto get_gapFreq_ref = [&](int idx) -> float {
        if (idx < 0 || idx >= N) return 1.0f; 
        return static_cast<float>(reference[idx][gapIdx]) / num.first;
    };
    auto get_gapFreq_qry = [&](int idx) -> float {
        if (idx < 0 || idx >= M) return 1.0f;
        return static_cast<float>(query[idx][gapIdx]) / num.second;
    };

    // Initial gap probabilities for boundary initialization
    float freq_ref_0 = get_gapFreq_ref(0);
    float freq_qry_0 = get_gapFreq_qry(0);

    // =================================================================
    // 1. Global Alignment Initialization (Top Row)
    // =================================================================
    M_prev[0] = 0.0f;
    I_prev[0] = MIN_INF;
    D_prev[0] = MIN_INF;
    for (int j = 1; j <= M; ++j) {
        M_prev[j] = MIN_INF;
        D_prev[j] = MIN_INF;
        I_prev[j] = j * param.gapTerminal;
    }

    // =================================================================
    // DP Loop
    // =================================================================
    for (int i = 1; i <= N; ++i) {
        // 1. Global Alignment Initialization (Left Column)
        M_curr[0] = MIN_INF;
        I_curr[0] = MIN_INF;
        // Boundary open and end penalties
        float open_pen = 0.5f * param.gapOpen;
        float end_pen  = gapEnd_ref[i - 1] * (1.0f - freq_qry_0);
        D_curr[0] = i * param.gapTerminal;

        for (int j = 1; j <= M; ++j) {
            int idx = i * (M + 1) + j;

            float freq_ref_i = get_gapFreq_ref(i - 1); 
            float freq_qry_j = get_gapFreq_qry(j - 1); 

            // Gap open at current position (j-1 / i-1), end at next position (j / i)
            // -- Penalty for D_mat (Gap in query, residue in reference) --
            float pos_gapOpen_ref = (j == M) ? param.gapTerminal : (gapOpen_ref[i - 1] * (1.0f - freq_qry_j));
            float pos_gapExtend_ref = (j == M) ? param.gapTerminal : (param.gapExtend * (1.0f - freq_ref_i));
            float pos_gapEnd_ref = (i == 1) ? 0.0f : (gapEnd_ref[i - 2] * (1.0f - get_gapFreq_qry(j)));
            
            // -- Penalty for I_mat (Gap in reference, residue in query) --
            float pos_gapOpen_qry = (i == N) ? param.gapTerminal : (gapOpen_qry[j - 1] * (1.0f - freq_ref_i));
            float pos_gapExtend_qry = (i == N) ? param.gapTerminal : (param.gapExtend * (1.0f - freq_qry_j));
            float pos_gapEnd_qry = (j == 1) ? 0.0f : (gapEnd_qry[j - 2] * (1.0f - get_gapFreq_ref(i)));

            // -- Calculate I_mat (Insertion in Reference) --
            float i_open = M_curr[j - 1] + pos_gapOpen_qry;
            float i_ext  = I_curr[j - 1] + pos_gapExtend_qry;
            if (i_open >= i_ext) { I_curr[j] = i_open; tb_I[idx] = STATE_M; } 
            else                 { I_curr[j] = i_ext;  tb_I[idx] = STATE_I; }

            // -- Calculate D_mat (Deletion in Query) --
            float d_open = M_prev[j] + pos_gapOpen_ref;
            float d_ext  = D_prev[j] + pos_gapExtend_ref;
            if (d_open >= d_ext) { D_curr[j] = d_open; tb_D[idx] = STATE_M; } 
            else                 { D_curr[j] = d_ext;  tb_D[idx] = STATE_D; }

            // -- Calculate M_mat --
            float consistencyBonus = 0.0f;
            if (consistencyTable != nullptr &&
                i - 1 < static_cast<int32_t>(consistencyTable->size()) &&
                j - 1 < static_cast<int32_t>((*consistencyTable)[i - 1].size())) {
                consistencyBonus = consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;

            float m_from_m = M_prev[j - 1];
            float m_from_i = I_prev[j - 1] + pos_gapEnd_qry; 
            float m_from_d = D_prev[j - 1] + pos_gapEnd_ref; 

            float max_m = m_from_m; 
            tb_M[idx] = STATE_M;

            if (m_from_i > max_m) { max_m = m_from_i; tb_M[idx] = STATE_I; }
            if (m_from_d > max_m) { max_m = m_from_d; tb_M[idx] = STATE_D; }

            M_curr[j] = match_score + max_m;
        }

        M_prev = M_curr;
        I_prev = I_curr;
        D_prev = D_curr;
    }

    // =================================================================
    // 3. Global Traceback strictly starts from the bottom-right corner (N, M)
    // =================================================================
    float best_score = M_prev[M]; 
    uint8_t state = STATE_M;

    if (I_prev[M] > best_score) { best_score = I_prev[M]; state = STATE_I; }
    if (D_prev[M] > best_score) { best_score = D_prev[M]; state = STATE_D; }

    std::vector<int8_t> path;
    int curr_i = N;
    int curr_j = M;

    while (curr_i > 0 || curr_j > 0) {
        if (curr_i == 0) {
            path.push_back(1);
            curr_j--;
        } 
        else if (curr_j == 0) {
            path.push_back(2);
            curr_i--;
        } 
        else {
            int idx = curr_i * (M + 1) + curr_j;
            if (state == STATE_M) {
                path.push_back(0);
                state = tb_M[idx];
                curr_i--;
                curr_j--;
            } 
            else if (state == STATE_I) {
                path.push_back(1);
                state = tb_I[idx];
                curr_j--;
            } 
            else if (state == STATE_D) {
                path.push_back(2);
                state = tb_D[idx];
                curr_i--;
            }
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

namespace {
    static constexpr int FAMSA_GAP_OPEN       = 25;
    static constexpr int FAMSA_GAP_EXT        = 26;
    static constexpr int FAMSA_GAP_TERM_EXT   = 27;
    static constexpr int FAMSA_GAP_TERM_OPEN  = 28;
    static constexpr int FAMSA_NO_SYMBOLS     = 32;
    static constexpr int FAMSA_NO_AMINOACIDS  = 21;

    struct FAMSA_DPMatrix {
        size_t n_rows;
        size_t n_cols;
        std::vector<uint8_t> raw_data;

        FAMSA_DPMatrix(size_t rows, size_t cols) : n_rows(rows), n_cols(cols), raw_data(rows * cols, 0) {}

        inline int get_dir_D(size_t r, size_t c) const {
            return raw_data[r * n_cols + c] & 0x03;
        }
        inline int get_dir_H(size_t r, size_t c) const {
            return (raw_data[r * n_cols + c] >> 2) & 0x03;
        }
        inline int get_dir_V(size_t r, size_t c) const {
            return (raw_data[r * n_cols + c] >> 4) & 0x03;
        }

        inline void set_dir_D(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xFC) | (dir & 0x03);
        }
        inline void set_dir_H(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xF3) | ((dir & 0x03) << 2);
        }
        inline void set_dir_V(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xCF) | ((dir & 0x03) << 4);
        }
        inline void set_dir_all(size_t r, size_t c, int dir) {
            uint8_t x = dir & 0x03;
            raw_data[r * n_cols + c] = x | (x << 2) | (x << 4);
        }
    };

    inline void famsa_solve_gaps_start(
        size_t source_col_id, size_t prof_width, float prof_size,
        const msa::FAMSAProfile& profile,
        float& n_gap_open, float& n_gap_ext, float& n_gap_term_open, float& n_gap_term_ext)
    {
        n_gap_open = n_gap_ext = n_gap_term_open = n_gap_term_ext = 0.0f;
        if (source_col_id >= prof_width) {
            const float* source_col = profile.get_counters(source_col_id);
            float cnt = source_col[FAMSA_GAP_TERM_OPEN] + source_col[FAMSA_GAP_TERM_EXT];
            n_gap_term_ext = cnt;
            n_gap_term_open += (prof_size - cnt);
        } else {
            const float* next_col = profile.get_counters(source_col_id + 1);
            n_gap_term_open += next_col[FAMSA_GAP_TERM_OPEN];

            const float* source_col = profile.get_counters(source_col_id);
            n_gap_term_ext += source_col[FAMSA_GAP_TERM_OPEN];
            n_gap_term_ext += source_col[FAMSA_GAP_TERM_EXT];

            n_gap_ext = source_col[FAMSA_GAP_OPEN];
            n_gap_ext += source_col[FAMSA_GAP_EXT];

            n_gap_open = prof_size - n_gap_ext - n_gap_term_open - n_gap_term_ext;
        }
    }

    inline void famsa_solve_gaps_cont(
        size_t source_col_id, size_t prof_width, float prof_size,
        const msa::FAMSAProfile& profile,
        float& n_gap_ext, float& n_gap_term_ext)
    {
        n_gap_ext = n_gap_term_ext = 0.0f;
        if (source_col_id == prof_width) {
            n_gap_term_ext = prof_size;
            n_gap_ext = 0.0f;
        } else {
            const float* next_col = profile.get_counters(source_col_id + 1);
            n_gap_term_ext = next_col[FAMSA_GAP_TERM_OPEN];

            const float* source_col = profile.get_counters(source_col_id);
            n_gap_term_ext += source_col[FAMSA_GAP_TERM_OPEN];
            n_gap_term_ext += source_col[FAMSA_GAP_TERM_EXT];

            n_gap_ext = prof_size - n_gap_term_ext;
        }
    }
} // anonymous namespace

std::vector<int8_t> msa::alignProfile_FAMSA(
    const FAMSAProfile& profile1,
    const FAMSAProfile& profile2,
    const Params& param,
    const std::vector<std::vector<float>>* consistencyTable,
    float consistencyWeight)
{
    size_t prof1_width = profile1.width;
    size_t prof2_width = profile2.width;

    if (prof1_width == 0 && prof2_width == 0) return {};
    if (prof1_width == 0) return std::vector<int8_t>(prof2_width, 1);
    if (prof2_width == 0) return std::vector<int8_t>(prof1_width, 2);

    float prof1_card = (profile1.totalWeight > 0.0f) ? profile1.totalWeight : (profile1.numSeqs > 0 ? static_cast<float>(profile1.numSeqs) : 1.0f);
    float prof2_card = (profile2.totalWeight > 0.0f) ? profile2.totalWeight : (profile2.numSeqs > 0 ? static_cast<float>(profile2.numSeqs) : 1.0f);

    float famsa_gap_open = param.gapOpen;
    float famsa_gap_ext = param.gapExtend;
    float famsa_gap_term_open = (param.gapTerminal != 0.0f) ? param.gapTerminal : (famsa_gap_open * (0.66f / 14.85f));
    float famsa_gap_term_ext = (param.gapTerminal != 0.0f) ? param.gapTerminal : (famsa_gap_ext * (0.66f / 1.25f));

    FAMSA_DPMatrix matrix(prof1_width + 1, prof2_width + 1);

    const float infty = 1e20f;
    struct dp_row_elem_t {
        float D, H, V;
    };
    std::vector<dp_row_elem_t> prev_row(prof2_width + 1);
    std::vector<dp_row_elem_t> curr_row(prof2_width + 1);

    struct dp_gap_costs {
        float open, ext, term_open, term_ext;
    };
    struct dp_gap_corrections {
        float n_gap_start_open, n_gap_start_ext, n_gap_start_term_open, n_gap_start_term_ext;
        float n_gap_cont_ext, n_gap_cont_term_ext;
    };

    std::vector<dp_gap_costs> prof2_gaps(prof2_width + 1);
    for (size_t j = 0; j <= prof2_width; ++j) {
        const float* s2 = profile2.get_scores(j);
        prof2_gaps[j].open      = s2[FAMSA_GAP_OPEN];
        prof2_gaps[j].ext       = s2[FAMSA_GAP_EXT];
        prof2_gaps[j].term_open = s2[FAMSA_GAP_TERM_OPEN];
        prof2_gaps[j].term_ext  = s2[FAMSA_GAP_TERM_EXT];
    }

    std::vector<dp_gap_corrections> gap_corrections(prof2_width + 1);
    std::vector<float> n_gaps_prof2_to_change(prof2_width + 1, 0.0f);
    std::vector<float> n_gaps_prof2_term_to_change(prof2_width + 1, 0.0f);
    std::vector<float> gaps_prof2_change(prof2_width + 1, 0.0f);

    for (size_t j = 1; j <= prof2_width; ++j) {
        famsa_solve_gaps_start(j, prof2_width, prof2_card, profile2,
            gap_corrections[j].n_gap_start_open, gap_corrections[j].n_gap_start_ext,
            gap_corrections[j].n_gap_start_term_open, gap_corrections[j].n_gap_start_term_ext);

        famsa_solve_gaps_cont(j, prof2_width, prof2_card, profile2,
            gap_corrections[j].n_gap_cont_ext, gap_corrections[j].n_gap_cont_term_ext);

        const float* c2 = profile2.get_counters(j);
        n_gaps_prof2_to_change[j]      = c2[FAMSA_GAP_OPEN];
        n_gaps_prof2_term_to_change[j] = c2[FAMSA_GAP_TERM_OPEN];

        gaps_prof2_change[j] = n_gaps_prof2_to_change[j] * (famsa_gap_ext - famsa_gap_open) +
                               n_gaps_prof2_term_to_change[j] * (famsa_gap_term_ext - famsa_gap_term_open);
    }

    // Boundary conditions for row 0
    prev_row[0].D = 0.0f;
    prev_row[0].H = -infty;
    prev_row[0].V = -infty;

    if (prof2_width >= 1) {
        prev_row[1].D = -infty;
        prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open * prof1_card;
        prev_row[1].V = -infty;
        matrix.set_dir_all(0, 1, 1);
    }

    for (size_t j = 2; j <= prof2_width; ++j) {
        prev_row[j].D = -infty;
        prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext * prof1_card;
        prev_row[j].V = -infty;
        matrix.set_dir_all(0, j, 1);
    }
    prev_row[prof2_width].H = -infty;

    // Calculate matrix interior
    for (size_t i = 1; i <= prof1_width; ++i) {
        const float* s1_i = profile1.get_scores(i);
        float prof1_gap_open_curr      = s1_i[FAMSA_GAP_OPEN];
        float prof1_gap_term_open_curr = s1_i[FAMSA_GAP_TERM_OPEN];
        float prof1_gap_ext_curr       = s1_i[FAMSA_GAP_EXT];
        float prof1_gap_term_ext_curr  = s1_i[FAMSA_GAP_TERM_EXT];

        curr_row[0].D = -infty;
        curr_row[0].H = -infty;
        matrix.set_dir_all(i, 0, 2);

        if (i < prof1_width) {
            if (i == 1)
                curr_row[0].V = std::max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
            else
                curr_row[0].V = std::max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;
        } else {
            curr_row[0].V = -infty;
        }

        const float* c1_i = profile1.get_counters(i);
        std::array<std::pair<int, float>, 32> col1;
        size_t col1_size = 0;
        float col1_n_non_gaps = 0.0f;
        for (int k = 0; k < FAMSA_NO_SYMBOLS; ++k) {
            float cnt = c1_i[k];
            if (cnt > 0.0f) {
                col1[col1_size++] = {k, cnt};
                if (k < FAMSA_NO_AMINOACIDS)
                    col1_n_non_gaps += cnt;
            }
        }

        float n_gap_prof1_start_open = 0.0f, n_gap_prof1_start_ext = 0.0f;
        float n_gap_prof1_start_term_open = 0.0f, n_gap_prof1_start_term_ext = 0.0f;
        famsa_solve_gaps_start(i, prof1_width, prof1_card, profile1,
            n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext);

        float n_gap_prof1_cont_ext = 0.0f, n_gap_prof1_cont_term_ext = 0.0f;
        famsa_solve_gaps_cont(i, prof1_width, prof1_card, profile1,
            n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext);

        float n_gaps_prof1_to_change      = c1_i[FAMSA_GAP_OPEN];
        float n_gaps_prof1_term_to_change = c1_i[FAMSA_GAP_TERM_OPEN];

        for (size_t j = 1; j <= prof2_width; ++j) {
            const float* scores2_column = profile2.get_scores(j);

            float t = 0.0f;
            for (size_t k = 0; k < col1_size; ++k) {
                t += col1[k].second * scores2_column[col1[k].first];
            }

            // Consistency bonus if available (cVal is per-pair average, so multiply by prof1_card * prof2_card to match FAMSA cardinality)
            if (consistencyTable && i - 1 < consistencyTable->size() && j - 1 < (*consistencyTable)[i - 1].size()) {
                float cVal = (*consistencyTable)[i - 1][j - 1];
                if (cVal > 0.0f) {
                    t += cVal * consistencyWeight * prof1_card * prof2_card;
                }
            }

            // State D
            float t_D = prev_row[j - 1].D + t;

            float t_H = prev_row[j - 1].H;
            if (n_gaps_prof1_to_change > 0.0f || n_gaps_prof1_term_to_change > 0.0f) {
                float delta_t = n_gaps_prof1_to_change * (scores2_column[FAMSA_GAP_EXT] - scores2_column[FAMSA_GAP_OPEN]) +
                                n_gaps_prof1_term_to_change * (scores2_column[FAMSA_GAP_TERM_EXT] - scores2_column[FAMSA_GAP_TERM_OPEN]);
                t_H = t_H + t + delta_t;
            } else {
                t_H = t_H + t;
            }

            float t_V = prev_row[j - 1].V + t + gaps_prof2_change[j] * col1_n_non_gaps;

            if (t_D > t_H && t_D > t_V) {
                curr_row[j].D = t_D;
                matrix.set_dir_D(i, j, 0);
            } else if (t_H > t_V) {
                curr_row[j].D = t_H;
                matrix.set_dir_D(i, j, 1);
            } else {
                curr_row[j].D = t_V;
                matrix.set_dir_D(i, j, 2);
            }

            // State H
            float gap_corr_H = prof2_gaps[j].open * n_gap_prof1_start_open +
                               prof2_gaps[j].ext * n_gap_prof1_start_ext +
                               prof2_gaps[j].term_open * n_gap_prof1_start_term_open +
                               prof2_gaps[j].term_ext * n_gap_prof1_start_term_ext;

            float t_D_H = curr_row[j - 1].D + gap_corr_H;
            float t_H_H = curr_row[j - 1].H + prof2_gaps[j].ext * n_gap_prof1_cont_ext +
                                             prof2_gaps[j].term_ext * n_gap_prof1_cont_term_ext;

            if (i > 1 && j > 1) {
                float t_V_H = curr_row[j - 1].V + gap_corr_H;
                if (t_D_H > t_H_H && t_D_H > t_V_H) {
                    curr_row[j].H = t_D_H;
                    matrix.set_dir_H(i, j, 0);
                } else if (t_V_H > t_H_H) {
                    curr_row[j].H = t_V_H;
                    matrix.set_dir_H(i, j, 2);
                } else {
                    curr_row[j].H = t_H_H;
                    matrix.set_dir_H(i, j, 1);
                }
            } else {
                if (t_D_H > t_H_H) {
                    curr_row[j].H = t_D_H;
                    matrix.set_dir_H(i, j, 0);
                } else {
                    curr_row[j].H = t_H_H;
                    matrix.set_dir_H(i, j, 1);
                }
            }

            // State V
            float gap_corr_V = prof1_gap_open_curr * gap_corrections[j].n_gap_start_open +
                               prof1_gap_ext_curr * gap_corrections[j].n_gap_start_ext +
                               prof1_gap_term_open_curr * gap_corrections[j].n_gap_start_term_open +
                               prof1_gap_term_ext_curr * gap_corrections[j].n_gap_start_term_ext;

            float t_D_V = prev_row[j].D + gap_corr_V;
            float t_V_V = prev_row[j].V + prof1_gap_ext_curr * gap_corrections[j].n_gap_cont_ext +
                                         prof1_gap_term_ext_curr * gap_corrections[j].n_gap_cont_term_ext;

            if (i > 1 && j > 1) {
                float t_H_V = prev_row[j].H + gap_corr_V;
                if (t_D_V > t_H_V && t_D_V > t_V_V) {
                    curr_row[j].V = t_D_V;
                    matrix.set_dir_V(i, j, 0);
                } else if (t_H_V > t_V_V) {
                    curr_row[j].V = t_H_V;
                    matrix.set_dir_V(i, j, 1);
                } else {
                    curr_row[j].V = t_V_V;
                    matrix.set_dir_V(i, j, 2);
                }
            } else {
                if (t_D_V > t_V_V) {
                    curr_row[j].V = t_D_V;
                    matrix.set_dir_V(i, j, 0);
                } else {
                    curr_row[j].V = t_V_V;
                    matrix.set_dir_V(i, j, 2);
                }
            }
        }
        curr_row.swap(prev_row);
    }

    // Traceback
    std::vector<int8_t> path;
    path.reserve(prof1_width + prof2_width);

    int cur_dir = 0;
    if (prev_row[prof2_width].D >= prev_row[prof2_width].H && prev_row[prof2_width].D >= prev_row[prof2_width].V)
        cur_dir = 0;
    else if (prev_row[prof2_width].H > prev_row[prof2_width].V)
        cur_dir = 1;
    else
        cur_dir = 2;

    size_t i = prof1_width, j = prof2_width;
    while (i > 0 || j > 0) {
        if (cur_dir == 0) {
            path.push_back(0); // Match
            cur_dir = matrix.get_dir_D(i, j);
            --i; --j;
        } else if (cur_dir == 1) {
            path.push_back(1); // Gap in Profile 1 (Query residue)
            cur_dir = matrix.get_dir_H(i, j);
            --j;
        } else {
            path.push_back(2); // Gap in Profile 2 (Ref residue)
            cur_dir = matrix.get_dir_V(i, j);
            --i;
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}


