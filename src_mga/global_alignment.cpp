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
