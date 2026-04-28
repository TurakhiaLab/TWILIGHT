#ifndef MSA_HPP
#include "../msa.hpp"
#endif

#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <cuda_runtime.h>
#include <tbb/parallel_for.h>
#include <numeric>

#define MIN_INF -30000 // Min value of int16_t

__device__ inline int get_idx(int i, int j, int totalSequence) {
    if (i > j) { int t = i; i = j; j = t; }
    return i * totalSequence - i * (i + 1) / 2 + (j - i - 1);
}

__device__ inline int8_t device_letterIdx(char type, char inChar) {
    char c = inChar;
    if (c >= 'a' && c <= 'z') c -= 32; 
    if (type == 'p') {
        switch(c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'D': return 2;
            case 'E': return 3;
            case 'F': return 4;
            case 'G': return 5;
            case 'H': return 6;
            case 'I': return 7;
            case 'K': return 8;
            case 'L': return 9;
            case 'M': return 10;
            case 'N': return 11;
            case 'P': return 12;
            case 'Q': return 13;
            case 'R': return 14;
            case 'S': return 15;
            case 'T': return 16;
            case 'V': return 17;
            case 'W': return 18;
            case 'Y': return 19;
            default: return 20;
        }
    } else { 
        switch(c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': case 'U': return 3;
            default: return 4;
        }
    }
}

__device__ inline int32_t get_pos_GPU(int u, int v, int posU, int totalSequence, uint64_t* d_pairOffsets, int32_t* d_forward, int32_t* d_backward) {
    if (u == v) return posU;
    
    // get_idx 巨集內建了 if (i > j) 的交換，所以 id 絕對正確
    int id = get_idx(u, v, totalSequence); 
    uint64_t offset = d_pairOffsets[id];
    
    if (u < v) {
        return d_forward[offset + posU];
    } else {
        return d_backward[offset + posU];
    }
}

// 1. 初始化 Kernel：使用 uint64_t 支援超大緊湊陣列
__global__ void initArrays_GPU(
    int32_t totalPairs,
    uint64_t totalElements, 
    int32_t* d_forward,
    int32_t* d_backward,
    float* d_weights
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < totalPairs) {
        d_weights[idx] = -1.0f;
    }
    for (uint64_t i = idx; i < totalElements; i += blockDim.x * gridDim.x) {
        d_forward[i] = -1;
        d_backward[i] = -1;
    }
}

// 3. Build Library Kernel：使用雙重 Offsets
__global__ void buildLibrary_GPU(
    int32_t actualPairs,
    int32_t totalSequenceCount,
    int2* d_pairIndices,
    int* d_seqIDs,
    int* d_offsets,
    char* d_seqs,
    int2* d_alignedPairs,
    int32_t* d_alignedPairsLen,
    uint64_t* d_jobOffsets,
    int32_t* d_forward,
    int32_t* d_backward,
    float* d_weights,
    uint64_t* d_pairOffsets
) {
    for (int pair = blockIdx.x * blockDim.x + threadIdx.x; pair < actualPairs; pair += gridDim.x * blockDim.x) {
        int seq1_idx = d_pairIndices[pair].x;
        int seq2_idx = d_pairIndices[pair].y;
        int seqA = d_seqIDs[seq1_idx];
        int seqB = d_seqIDs[seq2_idx];

        int id = get_idx(seqA, seqB, totalSequenceCount);

        uint64_t offset = d_pairOffsets[id];
        uint64_t jobOffset = d_jobOffsets[pair];

        int len = d_alignedPairsLen[pair];
        int local_match = 0;
        
        for (int k = 0; k < len; ++k) {
            int2 p = d_alignedPairs[jobOffset + k];
            int base = p.x;
            int align = p.y;
            
            char c1 = d_seqs[d_offsets[seq1_idx] + base];
            char c2 = d_seqs[d_offsets[seq2_idx] + align];
            if (c1 >= 'a' && c1 <= 'z') c1 -= 32;
            if (c2 >= 'a' && c2 <= 'z') c2 -= 32;
            if (c1 == c2) local_match++;

            if (seqA < seqB) {
                d_forward[offset + base] = align;
                d_backward[offset + align] = base;
            } else {
                d_forward[offset + align] = base;
                d_backward[offset + base] = align;
            }
        }
        d_weights[id] = (len > 0) ? (float)local_match / len : 0.0f;
    }
}

// 🚀 極速版 Sparse Sub-Matrix：徹底消滅 O(L^2) 的龐大記憶體消耗
__global__ void computePairWeight_GPU_SparseRow(
    int num_jobs,
    int2* d_computeJobs,        
    int* d_target_locals,       
    int num_targets,
    int totalSequence,
    int* d_seqLengths,          
    int32_t* d_forward, 
    int32_t* d_backward, 
    float* d_weights,
    uint64_t* d_pairOffsets,
    int32_t* d_sparse_posU,
    int32_t* d_sparse_posV,
    float* d_sparse_weight,
    int32_t* d_sparse_counts,
    int max_edges_per_pair,
    int max_seq_len
) {

    int bx = blockIdx.x;
    for (int job_idx = bx; job_idx < num_jobs; job_idx += gridDim.x) {

        int u = d_computeJobs[job_idx].x;
        int v = d_computeJobs[job_idx].y;
        int lenU = d_seqLengths[u];
        int lenV = d_seqLengths[v];

        int id_uv = get_idx(u, v, totalSequence);
        float w_uv = d_weights[id_uv];

        // 🔥 神奇魔法：宣告動態 Shared Memory
        extern __shared__ float s_num_V[]; 
        
        __shared__ int edge_count;
        if (threadIdx.x == 0) edge_count = 0;

        // 初始化 s_num_V
        for (int i = threadIdx.x; i < lenV; i += blockDim.x) {
            s_num_V[i] = 0.0f;
        }
        __syncthreads();

        for (int posU = 0; posU < lenU; ++posU) {
            
            // 1. Direct Edge
            int directV = get_pos_GPU(u, v, posU, totalSequence, d_pairOffsets, d_forward, d_backward);
            if (threadIdx.x == 0) {
                if (directV != -1 && directV >= 0 && directV < lenV) {
                    s_num_V[directV] = w_uv; 
                }
            }
            __syncthreads();

            // 2. Highway Extension
            for (int c_idx = threadIdx.x; c_idx < num_targets; c_idx += blockDim.x) {
                int c = d_target_locals[c_idx];
                if (c == u || c == v) continue;
            
                float w_uc = d_weights[get_idx(u, c, totalSequence)];
                if (w_uc > 0.0f) {
                    int posC = get_pos_GPU(u, c, posU, totalSequence, d_pairOffsets, d_forward, d_backward);
                    if (posC != -1) {
                        float w_cv = d_weights[get_idx(c, v, totalSequence)];
                        if (w_cv > 0.0f) {
                            int posV = get_pos_GPU(c, v, posC, totalSequence, d_pairOffsets, d_forward, d_backward);
                            if (posV != -1 && posV >= 0 && posV < lenV) {
                                float min_w = fminf(w_uc, w_cv);
                                // 🔥 這裡的 atomicAdd 是在 Shared Memory 執行，快如閃電！
                                atomicAdd(&s_num_V[posV], min_w);
                            }
                        }
                    }
                }
            }
            __syncthreads();

            // 3. 收集與清理
            for (int i = threadIdx.x; i < lenV; i += blockDim.x) {
                float w = s_num_V[i];
                if (w > 0.0f) {
                    // 這裡的 atomicAdd 是針對 __shared__ 變數，一樣超快
                    int write_idx = atomicAdd(&edge_count, 1);
                    if (write_idx < max_edges_per_pair) {
                        uint64_t base_idx = (uint64_t)job_idx * max_edges_per_pair + write_idx;
                        d_sparse_posU[base_idx] = posU;
                        d_sparse_posV[base_idx] = i;
                        d_sparse_weight[base_idx] = w;
                    }
                    s_num_V[i] = 0.0f; // 重置
                }
            }
            __syncthreads();
        }

        if (threadIdx.x == 0) {
            d_sparse_counts[job_idx] = (edge_count < max_edges_per_pair) ? edge_count : max_edges_per_pair;
        }
    }
}

__global__ void alignmentOnGPU_affine(
    int32_t numPairs,      
    int2* d_pairIndices,
    int* d_offsets,        
    char* d_seqs,          
    int2* d_alignedPairs,  
    int32_t* d_alignedPairsLen, 
    uint64_t* d_jobOffsets, // 【新增】精準的 Traceback 記憶體起點
    float* param,
    char type
) {
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    
    int matrixSize = (type == 'n') ? 5 : 21;
    int paramSize = (type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    int16_t gapOpen = static_cast<int16_t>(param[paramSize - 4]);
    int16_t gapExtend = static_cast<int16_t>(param[paramSize - 3]);

    const int T = 220;
    const int O = 200;

    __shared__ uint8_t tbDir[(T * T) / 2];
    __shared__ int16_t wf_H[3 * (T + 1)]; 
    __shared__ int16_t wf_E[3 * (T + 1)];
    __shared__ int16_t wf_F[3 * (T + 1)];

    __shared__ int16_t scoreMat[21 * 21];

    __shared__ int2 localPairs[2 * T];
    __shared__ int localLen; 

    __shared__ bool lastTile;
    __shared__ int16_t maxScore;
    __shared__ int32_t best_ti;
    __shared__ int32_t best_tj;

    __shared__ int16_t redBestS[512]; 
    __shared__ int32_t redBestI[512];
    __shared__ int32_t redBestJ[512];
    __shared__ uint8_t redBestTB[512]; 

    __shared__ int8_t shared_ref[T];
    __shared__ int8_t shared_qry[T];

    __shared__ int16_t last_max_score;
    __shared__ int16_t last_max_state;
    __shared__ uint8_t best_tb;

    for (int pair = bx; pair < numPairs; pair += gridDim.x) {
        if (tx == 0) { 
            lastTile = false; 
            last_max_score = 0; 
            last_max_state = 0; 
        }
        for (int i = tx; i < matrixSize * matrixSize; i += blockDim.x) {
            scoreMat[i] = static_cast<int16_t>(param[i]);
        }
        __syncthreads();

        int32_t currentPairPathLen = 0; 
        int32_t reference_idx = 0; 
        int32_t query_idx = 0;  

        int seq1_idx = d_pairIndices[pair].x;
        int seq2_idx = d_pairIndices[pair].y;

        int32_t refStart = d_offsets[seq1_idx];
        int32_t qryStart = d_offsets[seq2_idx];
        int32_t refTotalLen = d_offsets[seq1_idx + 1] - refStart;
        int32_t qryTotalLen = d_offsets[seq2_idx + 1] - qryStart;

        uint64_t pairOffset = d_jobOffsets[pair];


        while (true) {
            int32_t refLen = min(T, refTotalLen - reference_idx); 
            int32_t qryLen = min(T, qryTotalLen - query_idx); 
            int16_t localBest  = MIN_INF;
            int32_t localBestI = refLen;
            int32_t localBestJ = qryLen;
            uint8_t localBestTB = 0;

            if(tx == 0) {
                if ((reference_idx + refLen == refTotalLen) && (query_idx + qryLen == qryTotalLen)) {
                    lastTile = true;
                }
                best_ti = refLen; 
                best_tj = qryLen;
            }
            __syncthreads();

            for (int s = tx; s < refLen; s += blockDim.x) {
                shared_ref[s] = device_letterIdx(type, d_seqs[refStart + reference_idx + s]);
            }
            for (int s = tx; s < qryLen; s += blockDim.x) {
                shared_qry[s] = device_letterIdx(type, d_seqs[qryStart + query_idx + s]);
            }
            __syncthreads();

            for (int k = 0; k <= refLen + qryLen; ++k) {
                int curr_k   = (k % 3) * (T + 1);
                int pre_k    = ((k + 2) % 3) * (T + 1); 
                int prepre_k = ((k + 1) % 3) * (T + 1); 

                int i_start = max(0, k - qryLen);
                int i_end   = min(refLen, k);

                for (int i = i_start + tx; i <= i_end; i += blockDim.x) {
                    int j = k - i; 

                    int16_t h_val = MIN_INF;
                    int16_t e_val = MIN_INF;
                    int16_t f_val = MIN_INF;
                    uint8_t tb = 0;

                    if (i == 0 && j == 0) {
                        h_val = 0; tb = 0;
                    } 
                    else if (i == 0) {
                        if (reference_idx == 0 && query_idx == 0) {
                            h_val = 0; tb = 0;
                        } else {
                            int penalty = (last_max_state == 2 && j == 1) ? gapExtend : (j == 1 ? gapOpen : gapExtend);
                            h_val = wf_H[pre_k + i] + penalty;
                            tb = 3; 
                            if (j > 1 || last_max_state == 2) tb |= 0x08; 
                        }
                    } 
                    else if (j == 0) {
                        if (reference_idx == 0 && query_idx == 0) {
                            h_val = 0; tb = 0;
                        } else {
                            int penalty = (last_max_state == 1 && i == 1) ? gapExtend : (i == 1 ? gapOpen : gapExtend);
                            h_val = wf_H[pre_k + (i - 1)] + penalty;
                            tb = 2; 
                            if (i > 1 || last_max_state == 1) tb |= 0x04; 
                        }
                    } 
                    else {
                        char r_char = shared_ref[i - 1];
                        char q_char = shared_qry[j - 1];

                        int16_t H_up = wf_H[pre_k + (i - 1)];
                        int16_t E_up = wf_E[pre_k + (i - 1)];
                        int e_open = H_up + gapOpen;
                        int e_ext  = E_up + gapExtend;
                        if (e_open >= e_ext) {
                            e_val = e_open;
                        } else {
                            e_val = e_ext;
                            tb |= 0x04;
                        }

                        int16_t H_left = wf_H[pre_k + i];
                        int16_t F_left = wf_F[pre_k + i];
                        int f_open = H_left + gapOpen;
                        int f_ext  = F_left + gapExtend;
                        if (f_open >= f_ext) {
                            f_val = f_open;
                        } else {
                            f_val = f_ext;
                            tb |= 0x08;
                        }

                        int16_t H_diag = wf_H[prepre_k + (i - 1)];

                        int ref_idx_char = shared_ref[i - 1]; // 已經是 int 了
                        int qry_idx_char = shared_qry[j - 1]; // 已經是 int 了
                        int16_t sub_score = scoreMat[ref_idx_char * matrixSize + qry_idx_char];
                        // int diag_score = H_diag + sub_score;
                        int diag_score = max((int)MIN_INF, (int)H_diag + sub_score);
                        // Apply the same max() clamping to e_open, e_ext, f_open, and f_ext

                        int best = MIN_INF;
                        uint8_t h_src = 0;
                        if (diag_score > best) { best = diag_score; h_src = 1; }
                        if (e_val > best) { best = e_val; h_src = 2; }
                        if (f_val > best) { best = f_val; h_src = 3; }

                        h_val = best;
                        tb |= h_src;

                        // Use 4 bits to store tb
                        int linear_idx = (i - 1) * T + (j - 1);
                        int byte_idx = linear_idx >> 1;          
                        int shift = (linear_idx & 1) << 2;       

                        uint8_t old_val = tbDir[byte_idx];
                        old_val &= ~(0x0F << shift);             
                        old_val |= (tb << shift);                
                        tbDir[byte_idx] = old_val;         
                    }

                    wf_H[curr_k + i] = h_val;
                    wf_E[curr_k + i] = e_val;
                    wf_F[curr_k + i] = f_val;

                    if (!lastTile) {
                        int start_i = max(1, refLen - O);
                        int start_j = max(1, qryLen - O);
                        bool in_bottom_strip = (i >= start_i && i <= refLen && j >= 1 && j <= qryLen);
                        bool in_right_strip  = (i >= 1 && i < start_i && j >= start_j && j <= qryLen);
                        if (in_bottom_strip || in_right_strip) {
                            if (h_val > localBest) {
                                localBest  = h_val; localBestI = i; localBestJ = j; localBestTB = tb;
                            }
                        }
                    } else {
                        bool is_last_tile_ref = (reference_idx + refLen == refTotalLen);
                        bool is_last_tile_qry = (query_idx + qryLen == qryTotalLen);
                        if ((is_last_tile_ref && i == refLen) || (is_last_tile_qry && j == qryLen)) {
                            if (h_val > localBest) {
                                localBest = h_val; localBestI = i; localBestJ = j; localBestTB = tb;
                            }
                        }
                    }
                }
                __syncthreads();
            } 
            
            redBestS[tx] = localBest;
            redBestI[tx] = localBestI;
            redBestJ[tx] = localBestJ;
            redBestTB[tx] = localBestTB;
            __syncthreads();

            for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
                if (tx < stride) {
                    if (redBestS[tx + stride] > redBestS[tx]) {
                        redBestS[tx] = redBestS[tx + stride];
                        redBestI[tx] = redBestI[tx + stride];
                        redBestJ[tx] = redBestJ[tx + stride];
                        redBestTB[tx] = redBestTB[tx + stride];
                    }
                }
                __syncthreads();
            }

            if (tx == 0) {
                best_ti = redBestI[0];
                best_tj = redBestJ[0];
                best_tb = redBestTB[0];
            }
            __syncthreads();

            int ti = best_ti;
            int tj = best_tj;

            int next_ref_advance = ti;
            int next_qry_advance = tj;

            if(tx == 0) {
                localLen = 0;
                int currentState = 0;

                while (ti > 0 || tj > 0) {
                    if (reference_idx == 0 && query_idx == 0 && ti == 0 && tj == 0) break;
                    if (ti == 0 && tj == 0) break;

                    uint8_t tb;
                    if (ti == 0) { 
                        tb = 3; 
                        if (tj > 1 || last_max_state == 2) tb |= 0x08; 
                    } else if (tj == 0) { 
                        tb = 2; 
                        if (ti > 1 || last_max_state == 1) tb |= 0x04; 
                    } else {
                        int linear_idx = (ti - 1) * T + (tj - 1);
                        int byte_idx = linear_idx >> 1;
                        if (byte_idx < 0 || byte_idx >= (T * T) / 2) {
                            printf("CRASH: byte_idx %d out of bounds! ti=%d, tj=%d\n", byte_idx, ti, tj);
                            break; // Prevent the crash so you can see the log
                        }
                        int shift = (linear_idx & 1) << 2;
                        tb = (tbDir[byte_idx] >> shift) & 0x0F; // 取回 4-bit
                    }

                    if (currentState == 0) { 
                        uint8_t h_src = tb & 0x03;
                        if (h_src == 1) {
                            localPairs[localLen++] = {reference_idx + ti - 1, query_idx + tj - 1};
                            --ti; --tj;
                        } else if (h_src == 2) { 
                            currentState = 1; 
                        } else if (h_src == 3) { 
                            currentState = 2;
                        } else { 
                            if (ti > 0 && tj == 0) currentState = 1;
                            else if (tj > 0 && ti == 0) currentState = 2;
                            else break;
                        }
                    } 
                    else if (currentState == 1) { 
                        // 🔥 防撞牆：只有當 ti > 0 時，才允許繼續往上回溯
                        if (ti > 0) {
                            bool e_from_e = (tb & 0x04) != 0; 
                            --ti; 
                            if (!e_from_e) currentState = 0; 
                        } else {
                            // 撞到頂端邊界了，強制切回狀態 0 讓它轉彎
                            currentState = 0;
                        }
                    } 
                    else if (currentState == 2) { 
                        // 🔥 防撞牆：只有當 tj > 0 時，才允許繼續往左回溯
                        if (tj > 0) {
                            bool f_from_f = (tb & 0x08) != 0; 
                            --tj; 
                            if (!f_from_f) currentState = 0; 
                        } else {
                            // 撞到左側邊界了，強制切回狀態 0 讓它轉彎
                            currentState = 0;
                        }
                    }
                }
            }
            __syncthreads();

            for (int k = tx; k < localLen; k += blockDim.x) {
                int pos = currentPairPathLen + k;
                d_alignedPairs[pairOffset + pos] = localPairs[localLen - 1 - k];
            }
            __syncthreads();


            currentPairPathLen += localLen;
            reference_idx += next_ref_advance;
            query_idx     += next_qry_advance;

            if (tx == 0) {
                last_max_score = redBestS[0];
                uint8_t h_src = best_tb & 0x03;
                if (h_src == 1) last_max_state = 0;
                else if (h_src == 2) last_max_state = 1;
                else if (h_src == 3) last_max_state = 2;
                else last_max_state = 0;
            }
            __syncthreads();

            if (lastTile) break;      
        } 

        if (tx == 0) d_alignedPairsLen[pair] = currentPairPathLen;
        __syncthreads();
    }
}


namespace msa {
namespace accurate {
namespace gpu {

struct GPU_pointers_consistency {
    char* deviceSequences = nullptr;
    int* deviceOffsets = nullptr;
    int2* deviceAlignedPairs = nullptr;
    int32_t* deviceAlignedPairsLen = nullptr;
    float* deviceParam = nullptr;
    int2* devicePairIndices = nullptr;
    int* deviceSeqIDs = nullptr;

    int32_t* deviceForward = nullptr;
    int32_t* deviceBackward = nullptr;
    float* deviceWeights = nullptr;

    uint64_t* deviceJobOffsets = nullptr;
    uint64_t* devicePairOffsets = nullptr;

    char* hostSequences = nullptr;
    int* hostOffsets = nullptr;
    float* hostParam = nullptr;
    int2* hostPairIndices = nullptr;
    int* hostSeqIDs = nullptr;

    int* deviceSeqLengths = nullptr;
    int* hostSeqLengths = nullptr;

    int2* deviceComputeJobs = nullptr;
    int* deviceTargetLocals = nullptr;
    float* deviceChunkMatrix = nullptr;
    uint64_t* deviceChunkMatrixOffsets = nullptr;

    float* hostChunkMatrix = nullptr;
    uint64_t* hostChunkMatrixOffsets = nullptr;
    int2* hostComputeJobs = nullptr;
    
    std::vector<int32_t> hostForward;
    std::vector<int32_t> hostBackward;
    std::vector<float> hostWeights;

    int numSeqs = 0;
    int totalSeqLen = 0;
    int totalSequenceCount = 0; 
    int actualPairs = 0;
    int totalPairs = 0; 
    int paramSize = 0;
    uint64_t totalElements = 0;
    uint64_t totalJobElements = 0;

    void memAllocate(int numSeqs, int totalSeqLen, int sequenceCount, int actualPairs, int totalPairs, Option* option, Params& param, uint64_t elements, uint64_t jobElements) {
        this->numSeqs = numSeqs;
        this->totalSeqLen = totalSeqLen;
        this->totalSequenceCount = sequenceCount;
        this->actualPairs = actualPairs;
        this->totalPairs = totalPairs;
        this->totalElements = elements;
        this->totalJobElements = jobElements;
        this->paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;

        hostSequences = new char[totalSeqLen];
        hostOffsets = new int[numSeqs + 1];
        hostSeqIDs = new int[numSeqs];
        hostParam = new float[paramSize];
        hostPairIndices = new int2[actualPairs];

        for (int i = 0; i < param.matrixSize; ++i) {
            for (int j = 0; j < param.matrixSize; ++j) {
                hostParam[i * param.matrixSize + j] = param.scoringMatrix[i][j];
            }
        }
        hostParam[paramSize - 4] = param.gapOpen;
        hostParam[paramSize - 3] = param.gapExtend;
        hostParam[paramSize - 2] = param.gapBoundary;
        hostParam[paramSize - 1] = param.xdrop;

        cudaMalloc((void**)&deviceSequences, totalSeqLen * sizeof(char));
        cudaMalloc((void**)&deviceOffsets, (numSeqs + 1) * sizeof(int));
        cudaMalloc((void**)&deviceSeqIDs, numSeqs * sizeof(int));
        cudaMalloc((void**)&deviceParam, paramSize * sizeof(float));
        cudaMalloc((void**)&devicePairIndices, actualPairs * sizeof(int2));
        
        cudaMalloc((void**)&deviceJobOffsets, (actualPairs + 1) * sizeof(uint64_t));
        cudaMalloc((void**)&devicePairOffsets, (totalPairs + 1) * sizeof(uint64_t));
        cudaMalloc((void**)&deviceAlignedPairs, totalJobElements * sizeof(int2));
        cudaMalloc((void**)&deviceAlignedPairsLen, actualPairs * sizeof(int32_t));
        
        cudaMalloc((void**)&deviceForward, totalElements * sizeof(int32_t));
        cudaMalloc((void**)&deviceBackward, totalElements * sizeof(int32_t));
        cudaMalloc((void**)&deviceWeights, totalPairs * sizeof(float));

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/memAllocate [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void memHost2Device(const std::vector<std::string>& currentSequences, const std::vector<int>& seqIDs) {
        int currentOffset = 0;
        for (int i = 0; i < numSeqs; ++i) {
            hostOffsets[i] = currentOffset;
            hostSeqIDs[i] = seqIDs[i];
            std::copy(currentSequences[i].begin(), currentSequences[i].end(), hostSequences + currentOffset);
            currentOffset += currentSequences[i].length();
        }
        hostOffsets[numSeqs] = currentOffset;

        cudaMemcpy(deviceSequences, hostSequences, totalSeqLen * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceOffsets, hostOffsets, (numSeqs + 1) * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceSeqIDs, hostSeqIDs, numSeqs * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceParam, hostParam, paramSize * sizeof(float), cudaMemcpyHostToDevice);
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/memHost2Device [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void memAllocateChunk(uint64_t max_floats, int max_pairs, int totalSequenceCount) {
        cudaMalloc((void**)&deviceComputeJobs, max_pairs * sizeof(int2));
        cudaMalloc((void**)&deviceTargetLocals, totalSequenceCount * sizeof(int));
        cudaMalloc((void**)&deviceChunkMatrix, max_floats * sizeof(float));
        cudaMalloc((void**)&deviceChunkMatrixOffsets, (max_pairs + 1) * sizeof(uint64_t));
        cudaMalloc((void**)&deviceSeqLengths, totalSequenceCount * sizeof(int));

        hostChunkMatrix = new float[max_floats];
        hostChunkMatrixOffsets = new uint64_t[max_pairs + 1];
        hostComputeJobs = new int2[max_pairs];
        hostSeqLengths = new int[totalSequenceCount];
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/memAllocateChunk [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void freePhase1Memory() {
        if (deviceSequences) { cudaFree(deviceSequences); deviceSequences = nullptr; }
        if (deviceOffsets) { cudaFree(deviceOffsets); deviceOffsets = nullptr; }
        if (deviceSeqIDs) { cudaFree(deviceSeqIDs); deviceSeqIDs = nullptr; }
        if (deviceParam) { cudaFree(deviceParam); deviceParam = nullptr; }
        if (devicePairIndices) { cudaFree(devicePairIndices); devicePairIndices = nullptr; }
        if (deviceJobOffsets) { cudaFree(deviceJobOffsets); deviceJobOffsets = nullptr; }
        if (deviceAlignedPairs) { cudaFree(deviceAlignedPairs); deviceAlignedPairs = nullptr; }
        if (deviceAlignedPairsLen) { cudaFree(deviceAlignedPairsLen); deviceAlignedPairsLen = nullptr; }
        
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/freePhase1Memory [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (hostSequences) { delete[] hostSequences; hostSequences = nullptr; }
        if (hostOffsets) { delete[] hostOffsets; hostOffsets = nullptr; }
        if (hostSeqIDs) { delete[] hostSeqIDs; hostSeqIDs = nullptr; }
        if (hostParam) { delete[] hostParam; hostParam = nullptr; }
        if (hostPairIndices) { delete[] hostPairIndices; hostPairIndices = nullptr; }
    }

    void freeMemory() {
        delete[] hostSequences; delete[] hostOffsets; delete[] hostSeqIDs; delete[] hostParam; delete[] hostPairIndices;
        cudaFree(deviceSequences); cudaFree(deviceOffsets); cudaFree(deviceSeqIDs);
        cudaFree(deviceParam); cudaFree(devicePairIndices);
        cudaFree(deviceJobOffsets); cudaFree(devicePairOffsets);
        cudaFree(deviceAlignedPairs); cudaFree(deviceAlignedPairsLen);
        cudaFree(deviceForward); cudaFree(deviceBackward); cudaFree(deviceWeights);

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/freeMemory [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (deviceComputeJobs) cudaFree(deviceComputeJobs);
        if (deviceTargetLocals) cudaFree(deviceTargetLocals);
        if (deviceChunkMatrix) cudaFree(deviceChunkMatrix);
        if (deviceChunkMatrixOffsets) cudaFree(deviceChunkMatrixOffsets);
        if (deviceSeqLengths) cudaFree(deviceSeqLengths);

        if (hostChunkMatrix) delete[] hostChunkMatrix;
        if (hostChunkMatrixOffsets) delete[] hostChunkMatrixOffsets;
        if (hostComputeJobs) delete[] hostComputeJobs;
        if (hostSeqLengths) delete[] hostSeqLengths;
    }
};


std::shared_ptr<msa::accurate::SubtreeAccurateState> buildSubtreeAccurateState_GPU(SequenceDB* database, Option* option, Tree* tree, int subtreeIdx, Params& params) {
    auto time0 = std::chrono::high_resolution_clock::now();
    std::cout << "0. Prepare Data and Transfer to Device.\n";

    const auto& sequences = database->sequences;
    std::vector<int> activeSeqIdx;
    activeSeqIdx.reserve(sequences.size());
    for (std::size_t seqIdx = 0; seqIdx < sequences.size(); ++seqIdx) {
        if (!sequences[seqIdx]->lowQuality || option->noFilter) activeSeqIdx.push_back(seqIdx);
    }

    std::vector<std::string> currentSequences(sequences.size());
    std::vector<int> seqIDs(sequences.size());
    int maxSeqID = -1;
    int totalSeqLen = 0;
    std::vector<int> seqLengths(sequences.size(), 0);

    for (std::size_t idx = 0; idx < activeSeqIdx.size(); ++idx) {
        const std::size_t seqIdx = activeSeqIdx[idx];
        auto seq = msa::getCurrentSequence(sequences[seqIdx]);
        currentSequences[seqIdx] = seq;
        seqLengths[seqIdx] = seq.size();
        seqIDs[seqIdx] = static_cast<int>(sequences[seqIdx]->id);
        totalSeqLen += seq.size();
        if (seqIDs[seqIdx] > maxSeqID) maxSeqID = seqIDs[seqIdx];
    }

    int sequenceCount = maxSeqID + 1;
    auto accurateState = std::make_shared<SubtreeAccurateState>(subtreeIdx, sequenceCount, activeSeqIdx, seqLengths);
    auto& ConsistencyLibrary = accurateState->directLib;

    std::vector<std::pair<std::size_t, std::size_t>> pairJobs;
    auto addPair = [&](std::size_t a, std::size_t b) {
        if (ConsistencyLibrary.idx(a, b) == -1) return;
        if (a > b) std::swap(a, b);
        if (a != b) pairJobs.push_back({a, b});
    };

    // --- 【保留你原本完美的 Sparse/Dense 分群與 Reps 產生邏輯】 ---
    if (activeSeqIdx.size() <= ConsistencyLibrary.DENSE_LIMIT) {
        // Dense Mode
        ConsistencyLibrary.cluster_members.push_back({});
        for (std::size_t seqIdx : activeSeqIdx) {
            ConsistencyLibrary.cluster_members[0].push_back(seqIdx);
            int local_id = ConsistencyLibrary.global_to_local[seqIdx];
            ConsistencyLibrary.all_reps.push_back(local_id);
            ConsistencyLibrary.is_rep[local_id] = true;
        }
        ConsistencyLibrary.cluster_reps.push_back(ConsistencyLibrary.cluster_members[0]);
        
        pairJobs.reserve((activeSeqIdx.size() * (activeSeqIdx.size() - 1)) / 2);
        for (std::size_t i = 0; i < activeSeqIdx.size(); ++i) {
            for (std::size_t j = i + 1; j < activeSeqIdx.size(); ++j) {
                addPair(activeSeqIdx[i], activeSeqIdx[j]);
            }
        }
    } else {
        // Sparse Mode
        std::size_t targetClusters = ConsistencyLibrary.TARGET_CLUSTERS;
        std::size_t targetSize = std::max(static_cast<std::size_t>(1), activeSeqIdx.size() / targetClusters);
        
        std::vector<std::vector<int>> clusters;
        auto leftover = gatherClustersFromTree(tree->root, database, targetSize, clusters);
        if (!leftover.empty()) {
            if (clusters.empty()) clusters.push_back(leftover);
            else clusters.back().insert(clusters.back().end(), leftover.begin(), leftover.end());
        }

        int max_reps_per_cluster = ConsistencyLibrary.REPS_PER_CLUSTER;
        std::vector<int> representatives;
        representatives.reserve(clusters.size() * max_reps_per_cluster);
        
        int cluster_id = 0;
        for (const auto& cl : clusters) {
            std::vector<int> current_members;
            std::vector<int> current_reps;
            
            for (std::size_t s : cl) {
                int local_s = ConsistencyLibrary.global_to_local[s];
                ConsistencyLibrary.seq_to_cluster[local_s] = cluster_id;
                current_members.push_back(s);
            }
            ConsistencyLibrary.cluster_members.push_back(current_members);

            if (cl.size() <= static_cast<std::size_t>(max_reps_per_cluster)) {
                for (std::size_t s : cl) current_reps.push_back(s);
            } else {
                std::size_t chunk_size = cl.size() / max_reps_per_cluster;
                for (int i = 0; i < max_reps_per_cluster; ++i) {
                    std::size_t start = i * chunk_size;
                    std::size_t end = (i == max_reps_per_cluster - 1) ? cl.size() : (start + chunk_size);
                    std::size_t best_seq = cl[start];
                    int max_len = -1;
                    for (std::size_t j = start; j < end; ++j) {
                        int l = currentSequences[cl[j]].size();
                        if (l > max_len) { max_len = l; best_seq = cl[j]; }
                    }
                    current_reps.push_back(best_seq);
                }
            }
            
            for (int r : current_reps) {
                representatives.push_back(r);
                int local_r = ConsistencyLibrary.global_to_local[r];
                ConsistencyLibrary.all_reps.push_back(local_r);
                ConsistencyLibrary.is_rep[local_r] = true;
            }
            ConsistencyLibrary.cluster_reps.push_back(current_reps);
            cluster_id++;
            
            // Cluster: All-to-all
            for (std::size_t i = 0; i < cl.size(); ++i) {
                for (std::size_t j = i + 1; j < cl.size(); ++j) {
                    addPair(cl[i], cl[j]);
                }
            }
        }
        
        // Rep: All-to-all
        for (std::size_t i = 0; i < representatives.size(); ++i) {
            for (std::size_t j = i + 1; j < representatives.size(); ++j) {
                addPair(representatives[i], representatives[j]);
            }
        }
        std::sort(pairJobs.begin(), pairJobs.end());
        pairJobs.erase(std::unique(pairJobs.begin(), pairJobs.end()), pairJobs.end());
    }
    // 注意：確保你在這裡有把產生的 pairJobs 推進去！
    // -------------------------------------------------------------

    int actualPairs = pairJobs.size();
    int totalPairs = ConsistencyLibrary.totalPairs;

    // 🔥 1. 計算完美壓縮記憶體 Offset (極限省 VRAM)
    std::vector<uint64_t> hostPairOffsets(totalPairs + 1, 0);
    uint64_t current_element = 0;
    // 只為有實際被選為 Job 的 Pair 分配 forward/backward 空間
    for (int i = 0; i < actualPairs; ++i) {
        int u = pairJobs[i].first;
        int v = pairJobs[i].second;
        int pId = ConsistencyLibrary.idx(u, v);
        hostPairOffsets[pId] = current_element;
        current_element += std::max(seqLengths[u], seqLengths[v]);
    }
    hostPairOffsets[totalPairs] = current_element;

    std::vector<uint64_t> hostJobOffsets(actualPairs + 1, 0);
    uint64_t current_job = 0;
    std::vector<int2> hostPairIndices(actualPairs);
    for (int i = 0; i < actualPairs; ++i) {
        hostJobOffsets[i] = current_job;
        int u = pairJobs[i].first;
        int v = pairJobs[i].second;
        hostPairIndices[i] = {u, v};
        current_job += (seqLengths[u] + seqLengths[v]);
    }
    hostJobOffsets[actualPairs] = current_job;

    // 2. 配置 GPU
    GPU_pointers_consistency gpu_ptrs;
    gpu_ptrs.memAllocate(activeSeqIdx.size(), totalSeqLen, sequenceCount, actualPairs, totalPairs, option, params, current_element, current_job);
    
    // 只傳入 active 的 sequences 給 device
    std::vector<std::string> activeSeqsStrs;
    std::vector<int> activeSeqIDs;
    for(int idx : activeSeqIdx) {
        activeSeqsStrs.push_back(currentSequences[idx]);
        activeSeqIDs.push_back(seqIDs[idx]);
    }
    gpu_ptrs.memHost2Device(activeSeqsStrs, activeSeqIDs);

    cudaMemcpy(gpu_ptrs.deviceJobOffsets, hostJobOffsets.data(), (actualPairs + 1) * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_ptrs.devicePairOffsets, hostPairOffsets.data(), (totalPairs + 1) * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_ptrs.devicePairIndices, hostPairIndices.data(), actualPairs * sizeof(int2), cudaMemcpyHostToDevice);

    const int blockSize = 256;
    const int numBlocks = 4096;
    int grid_dim = (actualPairs < numBlocks) ? actualPairs : numBlocks;
    
    std::cout << "1. All-to-all Pairwise Alignment: [GPU executing]\r" << std::flush;
    cudaError_t err;

    // 3. 執行 Kernels
    initArrays_GPU<<<grid_dim, blockSize>>>(totalPairs, current_element, gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward, gpu_ptrs.deviceWeights);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: consistency.cu/initArrays_GPU [" << cudaGetErrorString(err) << "]" << std::endl;
        exit(EXIT_FAILURE);
    }

    alignmentOnGPU_affine<<<grid_dim, blockSize>>>(
        actualPairs, gpu_ptrs.devicePairIndices, gpu_ptrs.deviceOffsets,
        gpu_ptrs.deviceSequences, gpu_ptrs.deviceAlignedPairs, gpu_ptrs.deviceAlignedPairsLen,
        gpu_ptrs.deviceJobOffsets, gpu_ptrs.deviceParam, option->type
    );
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: consistency.cu/alignmentOnGPU_affine [" << cudaGetErrorString(err) << "]" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "1-2. Build library.\n";

    buildLibrary_GPU<<<grid_dim, blockSize>>>(
        actualPairs, sequenceCount, gpu_ptrs.devicePairIndices, gpu_ptrs.deviceSeqIDs,
        gpu_ptrs.deviceOffsets, gpu_ptrs.deviceSequences, gpu_ptrs.deviceAlignedPairs, gpu_ptrs.deviceAlignedPairsLen,
        gpu_ptrs.deviceJobOffsets, gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward, gpu_ptrs.deviceWeights, gpu_ptrs.devicePairOffsets
    );
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: consistency.cu/buildLibrary_GPU [" << cudaGetErrorString(err) << "]" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto time1 = std::chrono::high_resolution_clock::now();
    std::cerr << "1. All-to-all Pairwise Alignment: [" << actualPairs << "/" << actualPairs << " pairs aligned] (100.0%)\n";

    // 4. 🚀 將精準計算好的資料拉回 CPU
    gpu_ptrs.hostWeights.resize(totalPairs);
    gpu_ptrs.hostForward.resize(current_element);
    gpu_ptrs.hostBackward.resize(current_element);

    cudaMemcpy(gpu_ptrs.hostWeights.data(), gpu_ptrs.deviceWeights, totalPairs * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_ptrs.hostForward.data(), gpu_ptrs.deviceForward, current_element * sizeof(int32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(gpu_ptrs.hostBackward.data(), gpu_ptrs.deviceBackward, current_element * sizeof(int32_t), cudaMemcpyDeviceToHost);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: consistency.cu/Phase1 DtoH [" << cudaGetErrorString(err) << "]" << std::endl;
        exit(EXIT_FAILURE);
    }
    // 5. 填入 CPU 的 Sparse Records 中！
    ConsistencyLibrary.records.resize(actualPairs);
    for (int i = 0; i < actualPairs; ++i) {
        int u = pairJobs[i].first;
        int v = pairJobs[i].second;
        int pId = ConsistencyLibrary.idx(u, v);
        
        ConsistencyLibrary.pair_to_record_idx[pId] = i;
        auto& rec = ConsistencyLibrary.records[i];
        
        rec.weight = gpu_ptrs.hostWeights[pId];
        uint64_t offset = hostPairOffsets[pId];
        
        // 利用 Offset 切割出一維 vector，零記憶體拷貝開銷！
        rec.forward.assign(gpu_ptrs.hostForward.begin() + offset, gpu_ptrs.hostForward.begin() + offset + seqLengths[u]);
        rec.backward.assign(gpu_ptrs.hostBackward.begin() + offset, gpu_ptrs.hostBackward.begin() + offset + seqLengths[v]);
    }

    auto time2 = std::chrono::high_resolution_clock::now();

    // -------------------------------------------------------------
    // [區塊 4]：🚀 GPU 趁熱打鐵！啟動 Sparse Compaction Chunking 運算！
    // -------------------------------------------------------------
    bool is_dense = (activeSeqIdx.size() <= ConsistencyLibrary.DENSE_LIMIT);
    std::vector<int> target_locals;
    if (is_dense) {
        target_locals.resize(sequenceCount);
        std::iota(target_locals.begin(), target_locals.end(), 0);
    } else {
        target_locals = ConsistencyLibrary.all_reps;
        std::sort(target_locals.begin(), target_locals.end());
        target_locals.erase(std::unique(target_locals.begin(), target_locals.end()), target_locals.end());
    }

    std::vector<int2> computeJobs;
    for (size_t i = 0; i < target_locals.size(); ++i) {
        for (size_t j = i + 1; j < target_locals.size(); ++j) {
            computeJobs.push_back({target_locals[i], target_locals[j]});
        }
    }
    int total_compute_jobs = computeJobs.size();

    // 取得最長序列長度，作為 GPU Workspace 的基準
    int max_seq_len = 0;
    for (int i = 0; i < sequenceCount; ++i) {
        if (seqLengths[i] > max_seq_len) max_seq_len = seqLengths[i];
    }
    gpu_ptrs.freePhase1Memory();

    // ====================================================================
    // 🔥 全新 VRAM 配置：再也不用擔心 OOM！
    // ====================================================================
    // 因為改成了 Row-by-row，記憶體負擔極小，我們可以安全地將 Chunk 放大，並調高 Edge 容量上限
    int max_jobs_per_chunk = 10000; 
    int max_edges_per_pair = 20000; // 提升兩倍，避免 Dense 區域被截斷

    float* d_global_num_V;
    cudaMalloc(&d_global_num_V, (size_t)max_jobs_per_chunk * max_seq_len * sizeof(float));

    int32_t *d_sparse_posU, *d_sparse_posV, *d_sparse_counts;
    float *d_sparse_weight;
    cudaMalloc(&d_sparse_posU, (size_t)max_jobs_per_chunk * max_edges_per_pair * sizeof(int32_t));
    cudaMalloc(&d_sparse_posV, (size_t)max_jobs_per_chunk * max_edges_per_pair * sizeof(int32_t));
    cudaMalloc(&d_sparse_weight, (size_t)max_jobs_per_chunk * max_edges_per_pair * sizeof(float));
    cudaMalloc(&d_sparse_counts, max_jobs_per_chunk * sizeof(int32_t));

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: consistency.cu/Phase 4 cudaMalloc [" << cudaGetErrorString(err) << "]" << std::endl;
        exit(EXIT_FAILURE);
    }

    int32_t *h_sparse_posU, *h_sparse_posV, *h_sparse_counts;
    float *h_sparse_weight;
    size_t max_copy_elements = (size_t)max_jobs_per_chunk * max_edges_per_pair;

    cudaMallocHost(&h_sparse_posU, max_copy_elements * sizeof(int32_t));
    cudaMallocHost(&h_sparse_posV, max_copy_elements * sizeof(int32_t));
    cudaMallocHost(&h_sparse_weight, max_copy_elements * sizeof(float));
    cudaMallocHost(&h_sparse_counts, max_jobs_per_chunk * sizeof(int32_t));
    
    // 簡單的 Job 陣列
    int2* d_computeJobs;
    cudaMalloc(&d_computeJobs, max_jobs_per_chunk * sizeof(int2));
    // ====================================================================

    int* d_seqLengths;
    cudaMalloc(&d_seqLengths, sequenceCount * sizeof(int));
    cudaMalloc(&(gpu_ptrs.deviceTargetLocals), target_locals.size() * sizeof(int));
    cudaMemcpy(d_seqLengths, seqLengths.data(), sequenceCount * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_ptrs.deviceTargetLocals, target_locals.data(), target_locals.size() * sizeof(int), cudaMemcpyHostToDevice);

    int current_job_idx = 0;
    std::cerr << "3. Compute Pair Weights (GPU Opt Sparse)...\n";

    while (current_job_idx < total_compute_jobs) {
        auto gpu_st = std::chrono::high_resolution_clock::now();
        int chunk_job_count = std::min(max_jobs_per_chunk, total_compute_jobs - current_job_idx);

        cudaMemcpy(d_computeJobs, &computeJobs[current_job_idx], chunk_job_count * sizeof(int2), cudaMemcpyHostToDevice);

        // 🚀 啟動極速 Row-by-Row Kernel
        // 一個 Block 處理一個 job，所以 Grid Size 必須等於 chunk_job_count
        // 只需要計算 num_V 需要的 bytes 數
        size_t shared_mem_size = max_seq_len * sizeof(float);

        computePairWeight_GPU_SparseRow<<<chunk_job_count, 256, shared_mem_size>>>(
        // computePairWeight_GPU_SparseRow<<<chunk_job_count, 256>>>(
            chunk_job_count, d_computeJobs, gpu_ptrs.deviceTargetLocals, target_locals.size(),
            sequenceCount, d_seqLengths, gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward,
            gpu_ptrs.deviceWeights, gpu_ptrs.devicePairOffsets,
            d_sparse_posU, d_sparse_posV, d_sparse_weight, d_sparse_counts, max_edges_per_pair,
            max_seq_len
        );
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: consistency.cu/computePairWeight_GPU_SparseRow [" << cudaGetErrorString(err) << "]" << std::endl;
            exit(EXIT_FAILURE);
        }

        // 僅將輕量級 Sparse 資料拉回 CPU
        size_t copy_elements = (size_t)chunk_job_count * max_edges_per_pair;
        cudaMemcpy(h_sparse_counts, d_sparse_counts, chunk_job_count * sizeof(int32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_sparse_posU, d_sparse_posU, copy_elements * sizeof(int32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_sparse_posV, d_sparse_posV, copy_elements * sizeof(int32_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_sparse_weight, d_sparse_weight, copy_elements * sizeof(float), cudaMemcpyDeviceToHost);

        auto cpu_st = std::chrono::high_resolution_clock::now();
        // CPU 平行轉譯
        tbb::parallel_for(tbb::blocked_range<int>(0, chunk_job_count), [&](const tbb::blocked_range<int>& range) {
            for (int p = range.begin(); p < range.end(); ++p) {
                int u = computeJobs[current_job_idx + p].x;
                int v = computeJobs[current_job_idx + p].y;
                int lenU = seqLengths[u];
                int lenV = seqLengths[v];
                
                int pId = ConsistencyLibrary.idx(u, v);
                int rId = ConsistencyLibrary.pair_to_record_idx[pId];
                if (rId == -1) continue;
            
                auto& rec = ConsistencyLibrary.records[rId];
                
                int edge_count = h_sparse_counts[p];
                uint64_t base_idx = (uint64_t)p * max_edges_per_pair;
            
                // --- 1. 預先計算容量 (Counting Pass) ---
                std::vector<int> countU(lenU, 0);
                std::vector<int> countV(lenV, 0);
                for (int k = 0; k < edge_count; ++k) {
                    countU[h_sparse_posU[base_idx + k]]++;
                    countV[h_sparse_posV[base_idx + k]]++;
                }
            
                // --- 2. 一次性分配與 Reserve (消滅 Reallocation) ---
                // 使用 assign({}) 速度遠快於 resize() + clear()
                rec.extendedForward.assign(lenU, {});
                rec.extendedBackward.assign(lenV, {});
                
                for(int i = 0; i < lenU; ++i) rec.extendedForward[i].reserve(countU[i]);
                for(int i = 0; i < lenV; ++i) rec.extendedBackward[i].reserve(countV[i]);
            
                // --- 3. 安全且極速的寫入 ---
                for (int k = 0; k < edge_count; ++k) {
                    int posU = h_sparse_posU[base_idx + k];
                    int posV = h_sparse_posV[base_idx + k];
                    float w  = h_sparse_weight[base_idx + k];
                
                    rec.extendedForward[posU].push_back({posV, w});
                    rec.extendedBackward[posV].push_back({static_cast<int32_t>(posU), w});
                }
            }
        });

         auto cpu_en = std::chrono::high_resolution_clock::now();
         std::cout << "GPU time: " << std::chrono::duration_cast<std::chrono::milliseconds>(cpu_st - gpu_st).count() << " ms\n";
         std::cout << "CPU time: " << std::chrono::duration_cast<std::chrono::milliseconds>(cpu_en - cpu_st).count() << " ms\n";


        current_job_idx += chunk_job_count;
        double percent = 100.0 * current_job_idx / total_compute_jobs;
        std::cerr << "  [" << current_job_idx << "/" << total_compute_jobs << " pairs computed] ("
                  << std::fixed << std::setprecision(1) << percent << "%)\r" << std::flush;
    }
    std::cerr << "\n";

    cudaFreeHost(h_sparse_posU); cudaFreeHost(h_sparse_posV);
    cudaFreeHost(h_sparse_weight); cudaFreeHost(h_sparse_counts);

    cudaFree(d_global_num_V);
    cudaFree(d_sparse_posU);
    cudaFree(d_sparse_posV);
    cudaFree(d_sparse_weight);
    cudaFree(d_sparse_counts);
    cudaFree(d_computeJobs);
    
    gpu_ptrs.freeMemory();
    auto time3 = std::chrono::high_resolution_clock::now();
    

    ConsistencyLibrary.computeResidueSupport(); // 呼叫 CPU 原生的
    
    
    auto time4 = std::chrono::high_resolution_clock::now();
    
    if (option->printDetail) {
        auto ms = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1000000; };
        std::cerr << "Built accurate-mode direct library for subtree " << subtreeIdx
                  << " with " << actualPairs << " pairs in " << ms(time0, time4) << " ms.\n"
                  << "Time breakdown: \n"
                  << "  1. Build pairwise library (GPU): " << ms(time0, time1) << " ms\n"
                  << "  2. Transfer to CPU Struct:       " << ms(time1, time2) << " ms\n"
                  << "  3. Compute pair weights (GPU):   " << ms(time2, time3) << " ms\n"
                  << "  4. Compute residue support:      " << ms(time3, time4) << " ms\n";
    }

    return accurateState;
}

} // namespace gpu
} // namespace accurate
} // namespace msa