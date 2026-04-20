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

#define MIN_INF -30000 // Min value of int16_t

__device__ inline int get_idx(int i, int j, int totalSequence) {
    if (i > j) { int t = i; i = j; j = t; }
    return i * totalSequence - i * (i + 1) / 2 + (j - i - 1);
}

__global__ void buildLibrary_GPU(
    int32_t actualPairs,
    int32_t maxSeqLen,
    int32_t totalSequenceCount,
    int2* d_pairIndices,
    int* d_seqIDs,
    int* d_offsets,
    char* d_seqs,
    int2* d_alignedPairs,
    int32_t* d_alignedPairsLen,
    int32_t* d_forward,
    int32_t* d_backward,
    float* d_weights
) {
    for (int pair = blockIdx.x * blockDim.x + threadIdx.x; pair < actualPairs; pair += gridDim.x * blockDim.x) {
        int seq1_idx = d_pairIndices[pair].x;
        int seq2_idx = d_pairIndices[pair].y;
        int seqA = d_seqIDs[seq1_idx];
        int seqB = d_seqIDs[seq2_idx];

        int id = get_idx(seqA, seqB, totalSequenceCount);
        int offset = id * maxSeqLen;

        int len = d_alignedPairsLen[pair];
        
        int local_match = 0;
        for (int k = 0; k < len; ++k) {
            int2 p = d_alignedPairs[pair * 2 * maxSeqLen + k];
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

__global__ void computeResidueSupport_GPU(
    int32_t actualPairs,
    int32_t maxSeqLen,
    int32_t totalSequenceCount,
    int2* d_pairIndices,
    int* d_seqIDs,
    int32_t* d_forward,
    float* d_weights,
    float* d_residueSupport
) {
    for (int pair = blockIdx.x * blockDim.x + threadIdx.x; pair < actualPairs; pair += gridDim.x * blockDim.x) {
        int seq1_idx = d_pairIndices[pair].x;
        int seq2_idx = d_pairIndices[pair].y;
        int seqA = d_seqIDs[seq1_idx];
        int seqB = d_seqIDs[seq2_idx];
        if (seqA > seqB) { int t = seqA; seqA = seqB; seqB = t; }

        int id = get_idx(seqA, seqB, totalSequenceCount);
        float w = d_weights[id];
        if (w <= 0.0f) continue;

        int offset = id * maxSeqLen;

        for (int posI = 0; posI < maxSeqLen; ++posI) {
            int posJ = d_forward[offset + posI];
            if (posJ != -1) {
                atomicAdd(&d_residueSupport[seqA * maxSeqLen + posI], w);
                atomicAdd(&d_residueSupport[seqB * maxSeqLen + posJ], w);
            }
        }
    }
}

__global__ void computePairWeight_GPU(
    int32_t actualPairs,
    int32_t maxSeqLen,
    int32_t totalSequenceCount,
    int2* d_pairIndices,
    int* d_seqIDs,
    int32_t* const __restrict__ d_forward,
    int32_t* const __restrict__ d_backward,
    float* const __restrict__ d_weights,
    float* __restrict__ d_pairWeight
) {
    // 改變策略：讓一個 Block (例如 256 個 Threads) 合作處理一個 Pair
    for (int pair = blockIdx.x; pair < actualPairs; pair += gridDim.x) {
        
        int seq1_idx = d_pairIndices[pair].x;
        int seq2_idx = d_pairIndices[pair].y;
        int seqA = d_seqIDs[seq1_idx];
        int seqB = d_seqIDs[seq2_idx];
        
        if (seqA > seqB) { int t = seqA; seqA = seqB; seqB = t; }

        int id = get_idx(seqA, seqB, totalSequenceCount);
        float w_AB = d_weights[id];
        float base_weight = w_AB > 0.0f ? w_AB : 0.0f;
        int offset = id * maxSeqLen;

        // 1. 利用 Block 內的 Threads 平行填入 Base weight
        for (int posA = threadIdx.x; posA < maxSeqLen; posA += blockDim.x) {
            d_pairWeight[offset + posA] = base_weight;
        }
        
        // __syncthreads() 在這裡其實不需要，因為每個 thread 只會讀寫它專屬的 posA

        // 2. seqC 迴圈依然在外層，Block 內所有 Threads 會同步執行這個迴圈
        for (int seqC = 0; seqC < totalSequenceCount; ++seqC) {
            if (seqC == seqA || seqC == seqB) continue;

            int id_AC = get_idx(seqA, seqC, totalSequenceCount);
            float w_AC = d_weights[id_AC];
            if (w_AC <= 0.0f) continue; 

            int id_BC = get_idx(seqB, seqC, totalSequenceCount);
            float w_BC = d_weights[id_BC];
            if (w_BC <= 0.0f) continue;

            float min_w = (w_AC < w_BC) ? w_AC : w_BC;

            int offset_AC = id_AC * maxSeqLen;
            int offset_BC = id_BC * maxSeqLen;

            int32_t* ptr_C_A = (seqA < seqC) ? (d_forward + offset_AC) : (d_backward + offset_AC);
            int32_t* ptr_C_B = (seqB < seqC) ? (d_forward + offset_BC) : (d_backward + offset_BC);

            // 3. 核心加速點：用 threadIdx.x 來平行化 posA！
            // 此時記憶體存取將會達到「完美連續 (Perfectly Coalesced)」
            for (int posA = threadIdx.x; posA < maxSeqLen; posA += blockDim.x) {
                int posB = d_forward[offset + posA];
                if (posB == -1) continue;

                int posC_A = ptr_C_A[posA];
                if (posC_A == -1) continue;

                int posC_B = ptr_C_B[posB];
                if (posC_B == -1) continue;

                if (posC_A == posC_B) {
                    // 因為不同 Thread 負責不同的 posA，這裡絕對不會有 Race Condition！不需要 Atomic！
                    d_pairWeight[offset + posA] += min_w;
                }
            }
        }
    }
}

__global__ void initArrays_GPU(
    int32_t totalPairs,
    int32_t maxSeqLen,
    int32_t* d_forward,
    int32_t* d_backward,
    float* d_weights
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < totalPairs) {
        d_weights[idx] = -1.0f;
    }
    int totalElements = totalPairs * maxSeqLen;
    for (int i = idx; i < totalElements; i += blockDim.x * gridDim.x) {
        d_forward[i] = -1;
        d_backward[i] = -1;
    }
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

__global__ void alignmentOnGPU_affine(
    int32_t numPairs,      
    int32_t maxSeqLen,     
    int2* d_pairIndices,
    int* d_offsets,        
    char* d_seqs,          
    int2* d_alignedPairs,  
    int32_t* d_alignedPairsLen, 
    float* param,
    char type
) {
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    
    int matrixSize = (type == 'n') ? 5 : 21;
    int paramSize = (type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    int16_t gapOpen = static_cast<int16_t>(param[paramSize - 4]);
    int16_t gapExtend = static_cast<int16_t>(param[paramSize - 3]);

    const int T = 256;
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

        int32_t pairOffset = pair * 2 * maxSeqLen; 


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
                        int diag_score = H_diag + sub_score;

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
                        bool e_from_e = (tb & 0x04) != 0; 
                        --ti; 
                        if (!e_from_e) currentState = 0; 
                    } 
                    else if (currentState == 2) { 
                        bool f_from_f = (tb & 0x08) != 0; 
                        --tj; 
                        if (!f_from_f) currentState = 0; 
                    }
                }
            }
            __syncthreads();

            for (int k = tx; k < localLen; k += blockDim.x) {
                int pos = currentPairPathLen + k;
                if (pos < 2 * maxSeqLen) {
                    d_alignedPairs[pairOffset + pos] = localPairs[localLen - 1 - k];
                }
            }

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

        if (tx == 0) {
            d_alignedPairsLen[pair] = currentPairPathLen;
        }
        __syncthreads();
    } 
    __syncthreads();
}


namespace msa {
namespace accurate {
namespace gpu {

struct GPU_pointers_consistency {
    char* deviceSequences = nullptr;
    int* deviceOffsets = nullptr;
    float* deviceResidueSupport = nullptr;
    float* devicePairWeight = nullptr;
    int2* deviceAlignedPairs = nullptr;
    int32_t* deviceAlignedPairsLen = nullptr;
    float* deviceParam = nullptr;
    int2* devicePairIndices = nullptr;
    int* deviceSeqIDs = nullptr;

    int32_t* deviceForward = nullptr;
    int32_t* deviceBackward = nullptr;
    float* deviceWeights = nullptr;

    char* hostSequences = nullptr;
    int* hostOffsets = nullptr;
    float* hostResidueSupport = nullptr;
    float* hostPairWeight = nullptr;
    float* hostParam = nullptr;
    int2* hostPairIndices = nullptr;
    int* hostSeqIDs = nullptr;

    int numSeqs = 0;
    int totalSeqLen = 0;
    int maxLen = 0;
    int totalSequenceCount = 0; // maxSeqID + 1
    int actualPairs = 0;
    int totalPairs = 0; // based on totalSequenceCount
    int paramSize = 0;

    void memAllocate(int numSeqs, int totalSeqLen, int maxLen, int totalSequenceCount, Option* option, Params& param) {
        this->numSeqs = numSeqs;
        this->totalSeqLen = totalSeqLen;
        this->maxLen = maxLen;
        this->totalSequenceCount = totalSequenceCount;
        this->actualPairs = numSeqs * (numSeqs - 1) / 2;
        this->totalPairs = totalSequenceCount * (totalSequenceCount - 1) / 2;
        this->paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;

        // 1. Allocate host memory
        hostSequences = new char[totalSeqLen];
        hostOffsets = new int[numSeqs + 1];
        hostSeqIDs = new int[numSeqs];
        hostResidueSupport = new float[totalSequenceCount * maxLen];
        hostPairWeight = new float[totalPairs * maxLen];
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

        int p = 0;
        for (int i = 0; i < numSeqs - 1; ++i) {
            for (int j = i + 1; j < numSeqs; ++j) {
                hostPairIndices[p++] = {i, j};
            }
        }

        // 2. Allocate device memory
        cudaMalloc((void**)&deviceSequences, totalSeqLen * sizeof(char));
        cudaMalloc((void**)&deviceOffsets, (numSeqs + 1) * sizeof(int));
        cudaMalloc((void**)&deviceSeqIDs, numSeqs * sizeof(int));
        cudaMalloc((void**)&deviceResidueSupport, totalSequenceCount * maxLen * sizeof(float));
        cudaMalloc((void**)&devicePairWeight, totalPairs * maxLen * sizeof(float));
        cudaMalloc((void**)&deviceAlignedPairs, actualPairs * 2 * maxLen * sizeof(int2));
        cudaMalloc((void**)&deviceAlignedPairsLen, actualPairs * sizeof(int32_t));
        cudaMalloc((void**)&deviceParam, paramSize * sizeof(float));
        cudaMalloc((void**)&devicePairIndices, actualPairs * sizeof(int2));
        
        cudaMalloc((void**)&deviceForward, totalPairs * maxLen * sizeof(int32_t));
        cudaMalloc((void**)&deviceBackward, totalPairs * maxLen * sizeof(int32_t));
        cudaMalloc((void**)&deviceWeights, totalPairs * sizeof(float));

        std::string error = cudaGetErrorString(cudaGetLastError());
        if (error != "no error") {
            fprintf(stderr, "CUDA Error: Failed to allocate consistency memory [%s]\n", error.c_str());
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

        std::fill(hostResidueSupport, hostResidueSupport + (totalSequenceCount * maxLen), 0.0f);
        std::fill(hostPairWeight, hostPairWeight + (totalPairs * maxLen), 0.0f);

        cudaMemcpy(deviceSequences, hostSequences, totalSeqLen * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceOffsets, hostOffsets, (numSeqs + 1) * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceSeqIDs, hostSeqIDs, numSeqs * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceResidueSupport, hostResidueSupport, totalSequenceCount * maxLen * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(devicePairWeight, hostPairWeight, totalPairs * maxLen * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceParam, hostParam, paramSize * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(devicePairIndices, hostPairIndices, actualPairs * sizeof(int2), cudaMemcpyHostToDevice);
        
        std::string error = cudaGetErrorString(cudaGetLastError());
        if (error != "no error") {
            fprintf(stderr, "CUDA Error: Failed to copy consistency memory to device [%s]\n", error.c_str());
        }
    }

    void freeMemory() {
        if (hostSequences) delete[] hostSequences;
        if (hostOffsets) delete[] hostOffsets;
        if (hostSeqIDs) delete[] hostSeqIDs;
        if (hostResidueSupport) delete[] hostResidueSupport;
        if (hostPairWeight) delete[] hostPairWeight;
        if (hostParam) delete[] hostParam;
        if (hostPairIndices) delete[] hostPairIndices;

        if (deviceSequences) cudaFree(deviceSequences);
        if (deviceOffsets) cudaFree(deviceOffsets);
        if (deviceSeqIDs) cudaFree(deviceSeqIDs);
        if (deviceResidueSupport) cudaFree(deviceResidueSupport);
        if (devicePairWeight) cudaFree(devicePairWeight);
        if (deviceAlignedPairs) cudaFree(deviceAlignedPairs);
        if (deviceAlignedPairsLen) cudaFree(deviceAlignedPairsLen);
        if (deviceParam) cudaFree(deviceParam);
        if (devicePairIndices) cudaFree(devicePairIndices);

        if (deviceForward) cudaFree(deviceForward);
        if (deviceBackward) cudaFree(deviceBackward);
        if (deviceWeights) cudaFree(deviceWeights);
    }
};

std::shared_ptr<msa::accurate::SubtreeAccurateState> buildSubtreeAccurateState_GPU(SequenceDB* database, Option* option, int subtreeIdx, Params& params) {
    auto accurateState = std::make_shared<SubtreeAccurateState>();
    accurateState->subtreeIdx = subtreeIdx;

    auto time0 = std::chrono::high_resolution_clock::now();
    std::cout << "0. Prepare Data and Transfer to Device.\n";

    const auto& sequences = database->sequences;
    std::vector<std::string> currentSequences;
    std::vector<int> seqIDs;
    currentSequences.reserve(sequences.size());
    seqIDs.reserve(sequences.size());

    int maxLength = -1;
    int maxSeqID = -1;
    int totalSeqLen = 0;

    for (std::size_t seqIdx = 0; seqIdx < sequences.size(); ++seqIdx) {
        if (!sequences[seqIdx]->lowQuality || option->noFilter) {
            auto seq = std::string(sequences[seqIdx]->alnStorage[sequences[seqIdx]->storage], sequences[seqIdx]->len);
            currentSequences.push_back(seq);
            seqIDs.push_back(static_cast<int>(sequences[seqIdx]->id));

            int seqLen = static_cast<int>(seq.size());
            if (seqLen > maxLength) maxLength = seqLen;
            if (static_cast<int>(sequences[seqIdx]->id) > maxSeqID) maxSeqID = static_cast<int>(sequences[seqIdx]->id);
            totalSeqLen += seqLen;
        }
    }

    int numSeqs = currentSequences.size();
    int sequenceCount = maxSeqID + 1; // Match CPU logic where library uses maxSeqID + 1

    GPU_pointers_consistency gpu_ptrs;
    gpu_ptrs.memAllocate(numSeqs, totalSeqLen, maxLength, sequenceCount, option, params);
    gpu_ptrs.memHost2Device(currentSequences, seqIDs);

    // Call Alignment Kernel
    int blockSize = 256;
    int numBlocks = 4096;
    int actualPairs = gpu_ptrs.totalPairs;
    int grid_dim = (actualPairs < numBlocks) ? actualPairs : numBlocks;
    
    std::cout << "1. All-to-all Pairwise Alignment: [GPU executing]\r" << std::flush;

    alignmentOnGPU_affine<<<grid_dim, blockSize>>>(
        actualPairs,
        maxLength,
        gpu_ptrs.devicePairIndices,
        gpu_ptrs.deviceOffsets,
        gpu_ptrs.deviceSequences,
        gpu_ptrs.deviceAlignedPairs,
        gpu_ptrs.deviceAlignedPairsLen,
        gpu_ptrs.deviceParam,
        option->type
    );
    cudaDeviceSynchronize();

    std::string error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to execute alignmentOnGPU_affine [%s]\n", error.c_str());
    }

    std::cout << "1-1. Initialize library memory.\n";

    int totalPairs = gpu_ptrs.totalPairs;
    
    initArrays_GPU<<<grid_dim, blockSize>>>(
        totalPairs, maxLength,
        gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward, gpu_ptrs.deviceWeights
    );
    cudaDeviceSynchronize();

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to execute initArrays_GPU [%s]\n", error.c_str());
    }

    std::cerr << "1-2. Build library.\n";

    buildLibrary_GPU<<<grid_dim, blockSize>>>(
        actualPairs, maxLength, sequenceCount,
        gpu_ptrs.devicePairIndices, gpu_ptrs.deviceSeqIDs, gpu_ptrs.deviceOffsets,
        gpu_ptrs.deviceSequences, gpu_ptrs.deviceAlignedPairs, gpu_ptrs.deviceAlignedPairsLen,
        gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward, gpu_ptrs.deviceWeights
    );
    cudaDeviceSynchronize();

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to execute buildLibrary_GPU [%s]\n", error.c_str());
    }
    
    std::cerr << "1. All-to-all Pairwise Alignment: [" << actualPairs << "/" << actualPairs << " pairs aligned] (100.0%)\n";

    auto time1 = std::chrono::high_resolution_clock::now();
    
    std::cerr << "2. Compute Residue Support: ";
    computeResidueSupport_GPU<<<grid_dim, blockSize>>>(
        actualPairs, maxLength, sequenceCount,
        gpu_ptrs.devicePairIndices, gpu_ptrs.deviceSeqIDs,
        gpu_ptrs.deviceForward, gpu_ptrs.deviceWeights, gpu_ptrs.deviceResidueSupport
    );
    cudaDeviceSynchronize();

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to execute computeResidueSupport_GPU [%s]\n", error.c_str());
    }
    std::cerr << "Done.\n";

    auto time2 = std::chrono::high_resolution_clock::now();
    
    std::cout << "3. Compute Pair Weights: [" << actualPairs << "/" << actualPairs << "] (100.0%)\n";
    computePairWeight_GPU<<<grid_dim, blockSize>>>(
        actualPairs, maxLength, sequenceCount,
        gpu_ptrs.devicePairIndices, gpu_ptrs.deviceSeqIDs,
        gpu_ptrs.deviceForward, gpu_ptrs.deviceBackward, gpu_ptrs.deviceWeights, gpu_ptrs.devicePairWeight
    );
    cudaDeviceSynchronize();

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to execute computePairWeight_GPU [%s]\n", error.c_str());
    }

    accurateState->directLib = DirectPairLibrary(sequenceCount, maxLength);
    accurateState->directLib.residueSupport.resize(sequenceCount * maxLength, 0.0f);

    cudaMemcpy(accurateState->directLib.weights.data(), gpu_ptrs.deviceWeights, totalPairs * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(accurateState->directLib.forward.data(), gpu_ptrs.deviceForward, totalPairs * maxLength * sizeof(int32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(accurateState->directLib.backward.data(), gpu_ptrs.deviceBackward, totalPairs * maxLength * sizeof(int32_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(accurateState->directLib.pairWeight.data(), gpu_ptrs.devicePairWeight, totalPairs * maxLength * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(accurateState->directLib.residueSupport.data(), gpu_ptrs.deviceResidueSupport, sequenceCount * maxLength * sizeof(float), cudaMemcpyDeviceToHost);

    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to transfer data back to host [%s]\n", error.c_str());
    }
    auto time3 = std::chrono::high_resolution_clock::now();
    
    if (option->printDetail) {
        auto ms = [](auto a, auto b) {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1000000;
        };
        std::cerr << "Built accurate-mode direct library for subtree " << subtreeIdx
                  << " with " << actualPairs 
                  << " pairs in " << ms(time0, time3) << " ms.\n"
                  << "Time breakdown: \n"
                  << "  1. Build pairwise library: " << ms(time0, time1) << " ms\n"
                  << "  2. Compute residue support: " << ms(time1, time2) << " ms\n"
                  << "  3. Compute pair weights: " << ms(time2, time3) << " ms\n";
    }
    
    gpu_ptrs.freeMemory();
    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") {
        fprintf(stderr, "CUDA Error: Failed to free memory on device [%s]\n", error.c_str());
    }

    return accurateState;
}

} // namespace gpu
} // namespace accurate
} // namespace msa