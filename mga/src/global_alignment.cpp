#include "global_alignment.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


// ==========================================
// 1. Global Alignment with Affine Gap
// ==========================================
CigarString runGlobalAlignment(const std::string& ref, const std::string& qry) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;
    
    int n = ref.size();
    int m = qry.size();

    std::vector<std::vector<int>> M(n + 1, std::vector<int>(m + 1, INF));
    std::vector<std::vector<int>> I(n + 1, std::vector<int>(m + 1, INF));
    std::vector<std::vector<int>> D(n + 1, std::vector<int>(m + 1, INF));

    M[0][0] = 0;
    for (int i = 1; i <= n; i++) {
        D[i][0] = GAP_OPEN + i * GAP_EXT;
    }
    for (int j = 1; j <= m; j++) {
        I[0][j] = GAP_OPEN + j * GAP_EXT;
    }

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score = (ref[i - 1] == qry[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
            M[i][j] = std::max({M[i - 1][j - 1], I[i - 1][j - 1], D[i - 1][j - 1]}) + score;
            I[i][j] = std::max(M[i][j - 1] + GAP_OPEN, I[i][j - 1]) + GAP_EXT;
            D[i][j] = std::max(M[i - 1][j] + GAP_OPEN, D[i - 1][j]) + GAP_EXT;
        }
    }

    int finalScore = std::max({M[n][m], I[n][m], D[n][m]});
    int state = 0;
    if (finalScore == M[n][m]) state = 0;
    else if (finalScore == I[n][m]) state = 1;
    else state = 2;
    
    int i = n, j = m;
    int matches = 0;
    CigarString cigar_ops;
    auto add_op = [&](char op) {
        if (!cigar_ops.empty() && cigar_ops.back().second == op) {
            cigar_ops.back().first++;
        } else {
            cigar_ops.push_back({1, op});
        }
    };

    int alignment_length = 0;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            if (state == 0) {
                int score = (ref[i - 1] == qry[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
                if (ref[i - 1] == qry[j - 1]) matches++;
                
                if (M[i][j] == M[i - 1][j - 1] + score) state = 0;
                else if (M[i][j] == I[i - 1][j - 1] + score) state = 1;
                else state = 2;
                
                add_op('M');
                i--; j--;
            } else if (state == 1) {
                if (I[i][j] == M[i][j - 1] + GAP_OPEN + GAP_EXT) state = 0;
                else state = 1;
                add_op('I');
                j--;
            } else if (state == 2) {
                if (D[i][j] == M[i - 1][j] + GAP_OPEN + GAP_EXT) state = 0;
                else state = 2;
                add_op('D');
                i--;
            }
        } else if (i > 0 && j == 0) {
            add_op('D'); i--;
        } else if (j > 0 && i == 0) {
            add_op('I'); j--;
        }
    }

    std::reverse(cigar_ops.begin(), cigar_ops.end());
    std::string cigar_str = "";
    for (auto op : cigar_ops) {
        cigar_str += std::to_string(op.first) + op.second;
        alignment_length += op.first;
    }

    float approxIdentity = static_cast<float>(matches) / static_cast<float>(alignment_length);

    return cigar_ops;
}

// ==========================================
// 2. Semi-Global Alignment with Affine Gap
// 適用於 Query 被完整包含在 Ref 中的某一段 (允許 Ref 首尾有 Free Overhang)
// ==========================================
CigarString runSemiGlobalAlignment(const std::string& ref, const std::string& qry) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;
    
    bool swapped = false;
    std::string ref_seq = ref;
    std::string qry_seq = qry;

    if (ref_seq.size() < qry_seq.size()) {
        std::swap(ref_seq, qry_seq);
        swapped = true;
    }
    
    int n = ref_seq.size(); // Ref 通常較長 (若 ref 較短，已 swap)
    int m = qry_seq.size();

    std::vector<std::vector<int>> M(n + 1, std::vector<int>(m + 1, INF));
    std::vector<std::vector<int>> I(n + 1, std::vector<int>(m + 1, INF));
    std::vector<std::vector<int>> D(n + 1, std::vector<int>(m + 1, INF));

    // 初始化邊界 (Semi-Global 的關鍵差異)
    for (int i = 0; i <= n; i++) {
        M[i][0] = 0; // Query 可以在 Ref 的任何位置無損起步
    }
    for (int j = 1; j <= m; j++) {
        I[0][j] = GAP_OPEN + j * GAP_EXT; // Query 若先被消耗依然要扣分
    }

    int maxScore = INF;
    int endI = 0;
    int endState = 0;

    // 填表
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score = (ref_seq[i - 1] == qry_seq[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;

            M[i][j] = std::max({M[i - 1][j - 1], I[i - 1][j - 1], D[i - 1][j - 1]}) + score;
            I[i][j] = std::max(M[i][j - 1] + GAP_OPEN, I[i][j - 1]) + GAP_EXT;
            D[i][j] = std::max(M[i - 1][j] + GAP_OPEN, D[i - 1][j]) + GAP_EXT;

            // 記錄當 Query 完全對齊完畢時 (j == m)，在 Ref 各個位置的最高分
            if (j == m) {
                // 🌟 修改 1：將 >= 改為 >。如果分數平手，保留最早達到的位置 (強制集中在左側)
                // 🌟 修改 2：如果分數平手，且前面的最高分是處於 Gap 狀態，強制優先選 Match (狀態 0)
                if (M[i][j] > maxScore) { 
                    maxScore = M[i][j]; endI = i; endState = 0; 
                } else if (M[i][j] == maxScore && endState != 0) { 
                    endI = i; endState = 0; 
                }
                
                if (I[i][j] > maxScore) { maxScore = I[i][j]; endI = i; endState = 1; }
                if (D[i][j] > maxScore) { maxScore = D[i][j]; endI = i; endState = 2; }
            }
        }
    }

    std::vector<std::pair<int, char>> cigar_ops;
    auto add_op = [&](char op) {
        if (swapped) {
            if (op == 'I') op = 'D';
            else if (op == 'D') op = 'I';
        }
        if (!cigar_ops.empty() && cigar_ops.back().second == op) {
            cigar_ops.back().first++;
        } else {
            cigar_ops.push_back({1, op});
        }
    };

    // Right overhang of Ref (Query is shorter, Ref has remaining bases)
    for (int i = 0; i < n - endI; i++) add_op('D');

    int i = endI, j = m;
    int state = endState;
    int matches = 0;

    int alignment_length = 0;

    while (i > 0 && j > 0) {
        if (state == 0) {
            int score = (ref_seq[i - 1] == qry_seq[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
            if (ref_seq[i - 1] == qry_seq[j - 1]) matches++;
            
            if (M[i][j] == M[i - 1][j - 1] + score) state = 0;
            else if (M[i][j] == I[i - 1][j - 1] + score) state = 1;
            else state = 2;
            
            add_op('M');
            i--; j--;
            alignment_length++;
        } else if (state == 1) {
            if (I[i][j] == M[i][j - 1] + GAP_OPEN + GAP_EXT) state = 0;
            else state = 1;
            add_op('I'); j--;
            alignment_length++;
        } else if (state == 2) {
            if (D[i][j] == M[i - 1][j] + GAP_OPEN + GAP_EXT) state = 0;
            else state = 2;
            add_op('D'); i--;
            alignment_length++;
        }
    }

    // Left overhang
    while (i > 0) { add_op('D'); i--; }
    while (j > 0) { add_op('I'); j--; }

    std::reverse(cigar_ops.begin(), cigar_ops.end());
    std::string cigar_str = "";
    for (auto op : cigar_ops) cigar_str += std::to_string(op.first) + op.second;

    return cigar_ops;
}

CigarString runGlobalAlignment(const Consensus& refCons, const Consensus& qryCons) {
    return runGlobalAlignment(refCons.getConsensusString(), qryCons.getConsensusString());
}

CigarString runSemiGlobalAlignment(const Consensus& refCons, const Consensus& qryCons) {
    return runSemiGlobalAlignment(refCons.getConsensusString(), qryCons.getConsensusString());
}

// ==========================================
// 3. Tiling Alignment (GACT) CPU Implementation
// ==========================================
CigarString runTilingAlignment(const std::string& ref, const std::string& qry) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;

    // Tile configuration
    const int T = 200;        // Tile size
    const int O = 50;         // Overlap between tiles

    // Scoring scheme
    const int16_t MATCH = 2;
    const int16_t MISMATCH = -1;
    const int16_t GAP = -2;

    // Traceback direction constants
    const uint8_t DIR_DIAG = 1;
    const uint8_t DIR_UP   = 2;
    const uint8_t DIR_LEFT = 3;

    int32_t refTotalLen = static_cast<int32_t>(ref.size());
    int32_t qryTotalLen = static_cast<int32_t>(qry.size());

    bool lastTile = false;
    int16_t maxScore = 0;
    int32_t reference_idx = 0;
    int32_t query_idx = 0;

    std::vector<uint8_t> tbDir(T * T, 0);
    std::vector<int16_t> wf_scores(3 * (T + 1), -9999);
    std::vector<uint8_t> localPath(2 * T, 0);

    CigarString raw_cigar;
    auto add_op = [&](char op) {
        if (!raw_cigar.empty() && raw_cigar.back().second == op) {
            raw_cigar.back().first++;
        } else {
            raw_cigar.push_back({1, op});
        }
    };

    while (!lastTile) {
        int32_t refLen = std::min((int32_t)T, refTotalLen - reference_idx);
        int32_t qryLen = std::min((int32_t)T, qryTotalLen - query_idx);

        if ((reference_idx + refLen == refTotalLen) && (query_idx + qryLen == qryTotalLen)) {
            lastTile = true;
        }

        std::fill(wf_scores.begin(), wf_scores.end(), -9999);

        int32_t best_ti = refLen;
        int32_t best_tj = qryLen;

        // Wavefront Scoring Loop (Diagonal Traversal)
        for (int k = 0; k <= refLen + qryLen; ++k) {
            int curr_k   = (k % 3) * (T + 1);
            int pre_k    = ((k + 2) % 3) * (T + 1);
            int prepre_k = ((k + 1) % 3) * (T + 1);

            int i_start = std::max(0, k - qryLen);
            int i_end   = std::min(refLen, k);

            for (int i = i_start; i <= i_end; ++i) {
                int j = k - i;

                int16_t score = -9999;
                uint8_t direction = DIR_DIAG;

                if (i == 0 && j == 0) {
                    score = maxScore;
                    maxScore = -9999;
                } else if (i == 0) {
                    score = wf_scores[pre_k + i] + GAP;
                    direction = DIR_LEFT;
                } else if (j == 0) {
                    score = wf_scores[pre_k + (i - 1)] + GAP;
                    direction = DIR_UP;
                } else {
                    char r_char = ref[reference_idx + (i - 1)];
                    char q_char = qry[query_idx + (j - 1)];

                    int16_t score_diag = wf_scores[prepre_k + (i - 1)] + (r_char == q_char ? MATCH : MISMATCH);
                    int16_t score_up   = wf_scores[pre_k + (i - 1)] + GAP;
                    int16_t score_left = wf_scores[pre_k + i] + GAP;

                    score = score_diag;
                    direction = DIR_DIAG;

                    if (score_up > score) {
                        score = score_up;
                        direction = DIR_UP;
                    }
                    if (score_left > score) {
                        score = score_left;
                        direction = DIR_LEFT;
                    }
                }

                wf_scores[curr_k + i] = score;

                if (i > 0 && j > 0) {
                    tbDir[(i - 1) * T + (j - 1)] = direction;
                }

                if (!lastTile) {
                    if (i > (refLen - O) && j > (qryLen - O)) {
                        if (score >= maxScore) {
                            maxScore = score;
                            best_ti = i;
                            best_tj = j;
                        }
                    }
                }
            }
        } // End Wavefront Loop

        // Traceback
        int ti = (!lastTile) ? best_ti : refLen;
        int tj = (!lastTile) ? best_tj : qryLen;

        int next_ref_advance = ti;
        int next_qry_advance = tj;

        int localLen = 0;
        while (ti > 0 || tj > 0) {
            uint8_t dir;
            if (ti == 0) {
                dir = DIR_LEFT;
            } else if (tj == 0) {
                dir = DIR_UP;
            } else {
                dir = tbDir[(ti - 1) * T + (tj - 1)];
            }

            localPath[localLen++] = dir;

            if (dir == DIR_DIAG) { ti--; tj--; }
            else if (dir == DIR_UP) { ti--; }
            else { tj--; }
        }

        // Convert reversed local path to CIGAR
        for (int k = localLen - 1; k >= 0; --k) {
            uint8_t dir = localPath[k];
            if (dir == DIR_DIAG) add_op('M');
            else if (dir == DIR_UP) add_op('D');
            else if (dir == DIR_LEFT) add_op('I');
        }

        reference_idx += next_ref_advance;
        query_idx     += next_qry_advance;
    } // End Tile Loop

    return raw_cigar;
}

CigarString runTilingAlignment(const Consensus& refCons, const Consensus& qryCons) {
    return runTilingAlignment(refCons.getConsensusString(), qryCons.getConsensusString());
}

