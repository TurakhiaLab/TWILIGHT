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

AlignmentResult Aligner::align_affine_local (const std::string& reference, const std::string& query, char type, msa::Params& params)
{
    const int gOpen = static_cast<int>(params.gapOpen);
    const int gExt  = static_cast<int>(params.gapExtend);
    const int MIN_INF = -100000000; 

    const int totalRef = static_cast<int>(reference.size());
    const int totalQry = static_cast<int>(query.size());

    // 建立 Full-Matrix
    std::vector<std::vector<int>> score(totalRef + 1, std::vector<int>(totalQry + 1, 0));
    std::vector<std::vector<uint8_t>> traceback(totalRef + 1, std::vector<uint8_t>(totalQry + 1, 0));
    
    // E 用於記錄 Query 的 gap (對應 Up 方向)
    std::vector<int> E(totalQry + 1, MIN_INF);

    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    // 填表 (DP)
    for (int i = 1; i <= totalRef; ++i) {
        int F_val = MIN_INF; // F 用於記錄 Reference 的 gap (對應 Left 方向)
        for (int j = 1; j <= totalQry; ++j) {
            uint8_t tb = 0;

            // 計算 E (Up - Delete in query)
            int e_open = score[i - 1][j] + gOpen; 
            int e_ext  = E[j] + gExt;
            if (e_open >= e_ext) {
                E[j] = e_open;
            } else {
                E[j] = e_ext;
                tb |= 0x04; // affine gap extension flag for E
            }

            // 計算 F (Left - Insert in query)
            int f_open = score[i][j - 1] + gOpen;
            int f_ext  = F_val + gExt;
            if (f_open >= f_ext) {
                F_val = f_open;
            } else {
                F_val = f_ext;
                tb |= 0x08; // affine gap extension flag for F
            }

            // 計算 Diagonal (Match/Mismatch)
            int diag = score[i - 1][j - 1] + residueScore(reference[i - 1], query[j - 1], type, params);

            // Local alignment: 分數下限為 0
            int best = 0; 
            uint8_t h_src = 0;

            if (diag > best) { best = diag; h_src = 1; }
            if (E[j] > best) { best = E[j]; h_src = 2; }
            if (F_val > best) { best = F_val; h_src = 3; }

            score[i][j] = best;
            tb |= h_src;
            traceback[i][j] = tb;

            // 記錄全域最高分作為 traceback 起點
            if (best > max_score) {
                max_score = best;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback 階段
    int i = max_i;
    int j = max_j;
    int currentState = 0;
    std::vector<std::pair<int, int>> aligned_pairs;

    while (i > 0 && j > 0) {
        // Local alignment 遇到 0 就停止 (只有在一般狀態下才檢查，若正在處理 gap 則需先回到 match 狀態)
        if (score[i][j] <= 0 && currentState == 0) break; 
        
        uint8_t tb = traceback[i][j];
        
        if (currentState == 0) { 
            uint8_t h_src = tb & 0x03;
            if (h_src == 1) { // 來自 Diagonal
                aligned_pairs.push_back({i - 1, j - 1});
                --i; --j;
            } else if (h_src == 2) { // 來自 Up (E)
                currentState = 1; 
            } else if (h_src == 3) { // 來自 Left (F)
                currentState = 2;
            } else { 
                break; // 遇到 0 (沒有來源)
            }
        } 
        else if (currentState == 1) { // 處理 Up (E) 狀態
            bool e_from_e = (tb & 0x04) != 0; 
            --i; 
            if (!e_from_e) currentState = 0; // 回到 Match 狀態
        } 
        else if (currentState == 2) { // 處理 Left (F) 狀態
            bool f_from_f = (tb & 0x08) != 0; 
            --j; 
            if (!f_from_f) currentState = 0; // 回到 Match 狀態
        }
    }

    // 因為是從尾巴往前找，需要反轉
    std::reverse(aligned_pairs.begin(), aligned_pairs.end());

    // 計算結果
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

static constexpr float CONSISTENCY_ALPHA = 100.0f;

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
                consistencyBonus = CONSISTENCY_ALPHA * consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            // M_mat = Dot Product + Consistency Bonus
            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;
            
            
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
    int totalCells = (N + 1) * (M + 1); 
    std::vector<uint8_t> tb_M(totalCells, 0);
    std::vector<uint8_t> tb_I(totalCells, 0);
    std::vector<uint8_t> tb_D(totalCells, 0);

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
            int idx = i * (M + 1) + j;

            // 2. Dynamic Border Penalties for Trailing Gaps
            // If at the last column (j==M), downward movement is a trailing deletion
            // If at the last row (i==N), rightward movement is a trailing insertion
            float pos_gapOpen_ref   = (j == M) ? param.gapTerminal : gapOp[0][i - 1];
            float pos_gapExtend_ref = (j == M) ? param.gapTerminal : gapEx[0][i - 1];

            float pos_gapOpen_qry   = (i == N) ? param.gapTerminal : gapOp[1][j - 1];
            float pos_gapExtend_qry = (i == N) ? param.gapTerminal : gapEx[1][j - 1];
            
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
                consistencyBonus = CONSISTENCY_ALPHA * consistencyWeight * (*consistencyTable)[i - 1][j - 1];
            }

            float match_score = scoreProfileOptimized(reference[i - 1], transQry[j - 1], alphabetSize) + consistencyBonus;
            
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
