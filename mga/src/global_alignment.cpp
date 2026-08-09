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
// 3. Tiling Alignment (GACT) CPU Implementation with Affine Gap
// ==========================================
CigarString runTilingAlignment(const std::string& ref, const std::string& qry) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;

    // Tile configuration
    const int T = 200;        // Tile size
    const int O = 50;         // Overlap between tiles

    int32_t refTotalLen = static_cast<int32_t>(ref.size());
    int32_t qryTotalLen = static_cast<int32_t>(qry.size());

    bool lastTile = false;
    int maxScore = 0;
    int32_t reference_idx = 0;
    int32_t query_idx = 0;

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

        std::string ref_sub = ref.substr(reference_idx, refLen);
        std::string qry_sub = qry.substr(query_idx, qryLen);

        int n = refLen;
        int m = qryLen;

        std::vector<std::vector<int>> M(n + 1, std::vector<int>(m + 1, INF));
        std::vector<std::vector<int>> I(n + 1, std::vector<int>(m + 1, INF));
        std::vector<std::vector<int>> D(n + 1, std::vector<int>(m + 1, INF));

        M[0][0] = maxScore;
        for (int i = 1; i <= n; ++i) {
            D[i][0] = maxScore + GAP_OPEN + i * GAP_EXT;
        }
        for (int j = 1; j <= m; ++j) {
            I[0][j] = maxScore + GAP_OPEN + j * GAP_EXT;
        }

        int best_score = INF;
        int best_ti = n;
        int best_tj = m;

        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                int score = (ref_sub[i - 1] == qry_sub[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;

                int prev_max = std::max({M[i - 1][j - 1], I[i - 1][j - 1], D[i - 1][j - 1]});
                if (prev_max > INF / 2) {
                    M[i][j] = prev_max + score;
                }

                int i_from_m = (M[i][j - 1] > INF / 2) ? (M[i][j - 1] + GAP_OPEN + GAP_EXT) : INF;
                int i_from_i = (I[i][j - 1] > INF / 2) ? (I[i][j - 1] + GAP_EXT) : INF;
                I[i][j] = std::max(i_from_m, i_from_i);

                int d_from_m = (M[i - 1][j] > INF / 2) ? (M[i - 1][j] + GAP_OPEN + GAP_EXT) : INF;
                int d_from_d = (D[i - 1][j] > INF / 2) ? (D[i - 1][j] + GAP_EXT) : INF;
                D[i][j] = std::max(d_from_m, d_from_d);

                if (!lastTile) {
                    if (i > (n - O) && j > (m - O)) {
                        int cell_max = std::max({M[i][j], I[i][j], D[i][j]});
                        if (cell_max >= best_score) {
                            best_score = cell_max;
                            best_ti = i;
                            best_tj = j;
                        }
                    }
                }
            }
        }

        int ti = (!lastTile) ? best_ti : n;
        int tj = (!lastTile) ? best_tj : m;

        if (!lastTile) {
            maxScore = best_score;
        }

        int next_ref_advance = ti;
        int next_qry_advance = tj;

        // Traceback inside tile
        int cell_max = std::max({M[ti][tj], I[ti][tj], D[ti][tj]});
        int state = 0;
        if (cell_max == M[ti][tj]) state = 0;
        else if (cell_max == I[ti][tj]) state = 1;
        else state = 2;

        std::vector<char> local_ops;
        while (ti > 0 || tj > 0) {
            if (ti > 0 && tj > 0) {
                if (state == 0) {
                    int score = (ref_sub[ti - 1] == qry_sub[tj - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
                    if (M[ti][tj] == M[ti - 1][tj - 1] + score) state = 0;
                    else if (M[ti][tj] == I[ti - 1][tj - 1] + score) state = 1;
                    else state = 2;
                    local_ops.push_back('M');
                    ti--; tj--;
                } else if (state == 1) {
                    if (I[ti][tj] == M[ti][tj - 1] + GAP_OPEN + GAP_EXT) state = 0;
                    else state = 1;
                    local_ops.push_back('I');
                    tj--;
                } else if (state == 2) {
                    if (D[ti][tj] == M[ti - 1][tj] + GAP_OPEN + GAP_EXT) state = 0;
                    else state = 2;
                    local_ops.push_back('D');
                    ti--;
                }
            } else if (ti > 0 && tj == 0) {
                local_ops.push_back('D');
                ti--;
            } else if (tj > 0 && ti == 0) {
                local_ops.push_back('I');
                tj--;
            }
        }

        for (auto it = local_ops.rbegin(); it != local_ops.rend(); ++it) {
            add_op(*it);
        }

        reference_idx += next_ref_advance;
        query_idx     += next_qry_advance;
    } // End Tile Loop

    return raw_cigar;
}

CigarString runTilingAlignment(const Consensus& refCons, const Consensus& qryCons) {
    return runTilingAlignment(refCons.getConsensusString(), qryCons.getConsensusString());
}

// ==========================================
// 1b. Global Alignment with Linear Gap (No Affine)
// ==========================================
CigarString runGlobalAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;

    int n = ref.size();
    int m = qry.size();

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

    dp[0][0] = 0;
    for (int i = 1; i <= n; ++i) {
        dp[i][0] = i * gapPenalty;
    }
    for (int j = 1; j <= m; ++j) {
        dp[0][j] = j * gapPenalty;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int score = (ref[i - 1] == qry[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
            int match = dp[i - 1][j - 1] + score;
            int del = dp[i - 1][j] + gapPenalty;
            int ins = dp[i][j - 1] + gapPenalty;
            dp[i][j] = std::max({match, del, ins});
        }
    }

    int i = n, j = m;
    CigarString cigar_ops;
    auto add_op = [&](char op) {
        if (!cigar_ops.empty() && cigar_ops.back().second == op) {
            cigar_ops.back().first++;
        } else {
            cigar_ops.push_back({1, op});
        }
    };

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            int score = (ref[i - 1] == qry[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
            if (dp[i][j] == dp[i - 1][j - 1] + score) {
                add_op('M');
                i--; j--;
            } else if (dp[i][j] == dp[i - 1][j] + gapPenalty) {
                add_op('D');
                i--;
            } else {
                add_op('I');
                j--;
            }
        } else if (i > 0) {
            add_op('D');
            i--;
        } else {
            add_op('I');
            j--;
        }
    }

    std::reverse(cigar_ops.begin(), cigar_ops.end());
    return cigar_ops;
}

CigarString runGlobalAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty) {
    return runGlobalAlignmentLinearGap(refCons.getConsensusString(), qryCons.getConsensusString(), gapPenalty);
}

// ==========================================
// 2b. Semi-Global Alignment with Linear Gap (No Affine)
// ==========================================
CigarString runSemiGlobalAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;

    bool swapped = false;
    std::string ref_seq = ref;
    std::string qry_seq = qry;

    if (ref_seq.size() < qry_seq.size()) {
        std::swap(ref_seq, qry_seq);
        swapped = true;
    }

    int n = ref_seq.size();
    int m = qry_seq.size();

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

    for (int i = 0; i <= n; ++i) {
        dp[i][0] = 0;
    }
    for (int j = 1; j <= m; ++j) {
        dp[0][j] = j * gapPenalty;
    }

    int maxScore = INF;
    int endI = 0;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int score = (ref_seq[i - 1] == qry_seq[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
            int match = dp[i - 1][j - 1] + score;
            int del = dp[i - 1][j] + gapPenalty;
            int ins = dp[i][j - 1] + gapPenalty;
            dp[i][j] = std::max({match, del, ins});

            if (j == m) {
                if (dp[i][j] >= maxScore) {
                    maxScore = dp[i][j];
                    endI = i;
                }
            }
        }
    }

    CigarString cigar_ops;
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

    for (int i = 0; i < n - endI; ++i) add_op('D');

    int i = endI, j = m;
    while (i > 0 && j > 0) {
        int score = (ref_seq[i - 1] == qry_seq[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
        if (dp[i][j] == dp[i - 1][j - 1] + score) {
            add_op('M');
            i--; j--;
        } else if (dp[i][j] == dp[i - 1][j] + gapPenalty) {
            add_op('D');
            i--;
        } else {
            add_op('I');
            j--;
        }
    }

    while (i > 0) { add_op('D'); i--; }
    while (j > 0) { add_op('I'); j--; }

    std::reverse(cigar_ops.begin(), cigar_ops.end());
    return cigar_ops;
}

CigarString runSemiGlobalAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty) {
    return runSemiGlobalAlignmentLinearGap(refCons.getConsensusString(), qryCons.getConsensusString(), gapPenalty);
}

// ==========================================
// 3b. Tiling Alignment (GACT) CPU Implementation with Linear Gap
// ==========================================
CigarString runTilingAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty) {
    CigarString result;
    if (ref.empty() || qry.empty()) return result;

    const int T = 200;        // Tile size
    const int O = 50;         // Overlap between tiles

    int32_t refTotalLen = static_cast<int32_t>(ref.size());
    int32_t qryTotalLen = static_cast<int32_t>(qry.size());

    bool lastTile = false;
    int maxScore = 0;
    int32_t reference_idx = 0;
    int32_t query_idx = 0;

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

        std::string ref_sub = ref.substr(reference_idx, refLen);
        std::string qry_sub = qry.substr(query_idx, qryLen);

        int n = refLen;
        int m = qryLen;

        std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, INF));

        dp[0][0] = maxScore;
        for (int i = 1; i <= n; ++i) {
            dp[i][0] = maxScore + i * gapPenalty;
        }
        for (int j = 1; j <= m; ++j) {
            dp[0][j] = maxScore + j * gapPenalty;
        }

        int best_score = INF;
        int best_ti = n;
        int best_tj = m;

        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j) {
                int score = (ref_sub[i - 1] == qry_sub[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
                int match = (dp[i - 1][j - 1] > INF / 2) ? (dp[i - 1][j - 1] + score) : INF;
                int del   = (dp[i - 1][j] > INF / 2) ? (dp[i - 1][j] + gapPenalty) : INF;
                int ins   = (dp[i][j - 1] > INF / 2) ? (dp[i][j - 1] + gapPenalty) : INF;
                dp[i][j]  = std::max({match, del, ins});

                if (!lastTile) {
                    if (i > (n - O) && j > (m - O)) {
                        if (dp[i][j] >= best_score) {
                            best_score = dp[i][j];
                            best_ti = i;
                            best_tj = j;
                        }
                    }
                }
            }
        }

        int ti = (!lastTile) ? best_ti : n;
        int tj = (!lastTile) ? best_tj : m;

        if (!lastTile) {
            maxScore = best_score;
        }

        int next_ref_advance = ti;
        int next_qry_advance = tj;

        // Traceback inside tile
        std::vector<char> local_ops;
        int i = ti, j = tj;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                int score = (ref_sub[i - 1] == qry_sub[j - 1]) ? MATCH_SCORE : MISMATCH_PENALTY;
                if (dp[i][j] == dp[i - 1][j - 1] + score) {
                    local_ops.push_back('M');
                    i--; j--;
                } else if (dp[i][j] == dp[i - 1][j] + gapPenalty) {
                    local_ops.push_back('D');
                    i--;
                } else {
                    local_ops.push_back('I');
                    j--;
                }
            } else if (i > 0) {
                local_ops.push_back('D');
                i--;
            } else {
                local_ops.push_back('I');
                j--;
            }
        }

        for (auto it = local_ops.rbegin(); it != local_ops.rend(); ++it) {
            add_op(*it);
        }

        reference_idx += next_ref_advance;
        query_idx     += next_qry_advance;
    } // End Tile Loop

    return raw_cigar;
}

CigarString runTilingAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty) {
    return runTilingAlignmentLinearGap(refCons.getConsensusString(), qryCons.getConsensusString(), gapPenalty);
}


