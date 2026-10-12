#include "align-method.hpp"

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

} // namespace accurate
} // namespace msa
