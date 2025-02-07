#ifndef ALIGN_HPP
#include "align.cuh"
#endif

__global__ void alignGrpToGrp_freq(float* freq, int8_t *aln, int32_t* len, int32_t* num, int32_t* alnLen, int32_t *seqInfo,  float* gapOpen, float* gapExtend, float* param)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    const int fLen = FRONT_WAVE_LEN;
    const int p_marker = 128;
    const int maxDiv = FRONT_WAVE_LEN/THREAD_NUM;
    const int tile_aln_size = 2*p_marker;

    int32_t pairNum = seqInfo[0];    

    __syncthreads();
    if (bx < pairNum && alnLen[bx] >= 0) {
        __shared__ int8_t tile_aln[tile_aln_size];
        __shared__ int16_t S  [3*fLen];
        __shared__ int16_t CS [3*fLen];
        __shared__ int16_t I  [2*fLen];
        __shared__ int16_t CI [2*fLen];
        __shared__ int16_t D  [2*fLen];
        __shared__ int16_t CD [2*fLen];
        __shared__ int8_t tb  [33*p_marker];
        __shared__ int32_t idx [2]; 
        // [0] reference_idx
        // [1] query_idx
        __shared__ int16_t max_score_list [maxDiv];
        __shared__ int16_t sub_max_score_list [fLen/maxDiv]; 
        __shared__ int16_t max_score;
        __shared__ int16_t aln_idx;
        __shared__ bool last_tile;
        __shared__ bool converged; 
        __shared__ bool conv_logic;
        __shared__ bool global;
        __shared__ bool includeHead;
        __shared__ bool reachBoundary;
        
        int32_t seqLen = seqInfo[1];
        
        int32_t refLen = len[2*bx];
        int32_t qryLen = len[2*bx+1];
        float refNum =  __int2float_rn(num[2*bx]);
        float qryNum =  __int2float_rn(num[2*bx+1]);

        
        // float p_gapOpen = param[25];
        float p_gapExtend = param[26];
        // float p_gapClose = param[27];
        float p_xdrop = param[28];
        float p_scoreMat [25];
        float denominator = refNum * qryNum;
        int32_t refFreqStart = 6*(2*bx)*seqLen;
        int32_t qryFreqStart = 6*(2*bx+1)*seqLen;
                            

        for (int i = 0; i < 25; ++i) p_scoreMat[i] = param[i];
        
        int32_t tile = 0;
        int32_t last_k = 0;
        int32_t L[3], U[3];

        int32_t ftr_length [p_marker+1];
        int32_t ftr_lower_limit [p_marker+1];
        // initialize global values
        if (tx == 0) {
            last_tile = false;
            reachBoundary = false;
            global = (alnLen[bx] <= 1);
            includeHead = (alnLen[bx] % 2 == 0);
            alnLen[bx] = 0;
        } 
        if (tx < 2) {
            idx[tx] = 0;
        }

        __syncthreads();
        
        while (!last_tile) {
            int16_t inf = 2 * p_xdrop + 1;
            int32_t reference_idx = idx[0];
            int32_t query_idx = idx[1];
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            int16_t score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            int16_t max_score_prime = -(p_xdrop+1);
            
            int32_t tb_idx = 0;
            for (int32_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (int32_t i = 0; i < p_marker+1; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }
            if (tx == 0) {
                max_score = 0;
                converged = false;
                conv_logic = false;
            }
            int32_t conv_value = 0; 
            
            int8_t tb_state = 0;
            int8_t ptr = 0;  // Main pointer
            int8_t state = 0;
            
            int32_t prev_conv_s = -1;

            bool Iptr = false;
            bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
            __syncthreads();
            
            
            // Initialize shared memory
            if (tx < THREAD_NUM) {
                for (int i = tx; i < 3*fLen; i = i+THREAD_NUM) {                     
                    S[i] = -1; CS[i] = -1;
                    if (i < 2*fLen) { I[i] = -1; CI[i] = -1; D[i] = -1; CD[i] = -1;}
                    if (i < tile_aln_size) {tile_aln[i] = -1;}
                }
            }
            aln_idx = 0;
            __syncthreads();

            // Start alignment
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                if (tx < fLen/maxDiv) {
                    sub_max_score_list[tx] = -(p_xdrop+1); 
                }
                __syncthreads();
                if (tx < maxDiv) {
                    max_score_list[tx] = -(p_xdrop+1); 
                }
                __syncthreads();
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d No more cells to compute based on x-drop criteria. Align length = %d, length left (r,q): (%d,%d)\n", bx, alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }
                
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d Anti-diagonal larger than the maximum limit (%d)! align length = %d, length left (r,q): (%d,%d)\n",  bx,fLen, alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }

                
                // __syncthreads();
                ptr = 0;
                int32_t ftrLen = U[k%3] - L[k%3] + 1;
                int32_t tbLen = (((ftrLen-1) >> 1)+1);
                if (k <= p_marker) {
                    ftr_length[ftr_length_idx] = tbLen;
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    ftr_addr += tbLen;
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;      
                }
                __syncthreads();
                int32_t threadRound = ftrLen / THREAD_NUM + 1;
                int32_t lastRound = ftrLen % THREAD_NUM;
                if (lastRound == 0) {
                    threadRound -= 1;
                    lastRound = THREAD_NUM;
                }
                for (int32_t rn = 0; rn < threadRound; rn += 1) {
                    int32_t activeThread = (rn != threadRound-1) ? THREAD_NUM : lastRound;
                    if (tx < activeThread) {
                        int32_t i=L[k%3]+rn*THREAD_NUM+tx;
                        int32_t Lprime = max(0, k-reference_length + 1);
                        int32_t j= min(k, reference_length - 1) - (i-Lprime);
                        if (j < 0) { printf("tx: %d, ERROR: j less than 0.\n", tx);}
                        int16_t match = -inf, insExt = -inf, delExt = -inf;
                        int16_t insOp = -inf, delOp = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            int16_t similarScore = 0;
                            float numerator = 0;
                            int32_t refFreqIdx = reference_idx + j;
                            int32_t qryFreqIdx = query_idx + i;
                            for (int l=0; l<6; l++) {
                                for (int m=0; m<6; m++) {
                                    if (m == 5 && l == 5)      numerator += 0;
                                    else if (m == 5 || l == 5) numerator = __fmaf_rn(__fmul_rn(freq[refFreqStart+6*(refFreqIdx)+l], freq[qryFreqStart+6*(qryFreqIdx)+m]), p_gapExtend, numerator);
                                    else                       numerator = __fmaf_rn(__fmul_rn(freq[refFreqStart+6*(refFreqIdx)+l], freq[qryFreqStart+6*(qryFreqIdx)+m]), p_scoreMat[m*5+l], numerator);
                                }
                            }
                            similarScore = __float2int_rn(__fdiv_rn(numerator, denominator));
                            if (offsetDiag < 0) match = similarScore;
                            else                match = S[(k+1)%3*fLen+offsetDiag] + similarScore;
                        }
                        
                        int16_t pos_gapOpen_ref = __float2int_rn(gapOpen[2*bx*seqLen+reference_idx+j]);
                        int16_t pos_gapOpen_qry = __float2int_rn(gapOpen[(2*bx+1)*seqLen+query_idx+i]);
                        int16_t pos_gapExtend_ref = __float2int_rn(gapExtend[2*bx*seqLen+reference_idx+j]);
                        int16_t pos_gapExtend_qry = __float2int_rn(gapExtend[(2*bx+1)*seqLen+query_idx+i]);
                        
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            delOp =  S[(k+2)%3*fLen+offsetUp] + pos_gapOpen_ref;
                            delExt = D[(k+1)%2*fLen+offsetUp] + pos_gapExtend_ref;
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + pos_gapOpen_qry;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + pos_gapExtend_qry;
                        }
                        int16_t tempI, tempD, tempH;

                        // Use DPX instuctions to calculate Max
                        unsigned Ext = (delExt << 16) | (insExt & 0xFFFF);
                        unsigned Op =  (delOp << 16)  | (insOp & 0xFFFF);
                        unsigned ExtOp = __vibmax_s16x2(Ext, Op, &Dptr, &Iptr);
                        tempD = (ExtOp >> 16) & 0xFFFF;
                        tempI = ExtOp & 0xFFFF;
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;
                        // Use DPX instuctions to calculate Max and argMax
                        int32_t mat32 = (match << 16) | (0 & 0xFFFF);
                        int32_t ins32 = (tempI << 16) | (1 & 0xFFFF);
                        int32_t del32 = (tempD << 16) | (2 & 0xFFFF);
                        int32_t max32 = __vimax3_s32(mat32, ins32, del32);
                        tempH = (max32 >> 16) & 0xFFFF;
                        ptr = (max32) & 0xFF;
                        if (tempH < max_score-p_xdrop) tempH = -inf;
                        S[(k%3)*fLen+offset] = tempH;
                        score = tempH;
                        sub_max_score_list[tx] = score;

                        if (k == p_marker - 1) { // Convergence algorithm
                            CS[(k%3)*fLen+offset] = (3 << 8) | (i & 0xFF); 
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 8) | (i & 0xFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 8) | (i & 0xFF);
                            CD[(k%2)*fLen+offset] = (2 << 8) | (i & 0xFF);
                        } 
                        else if (k >= p_marker + 1){
                            if (Iptr) {
                                CI[(k%2)*fLen+offset] = CI[((k+1)%2)*fLen+offsetLeft]; 
                            } else {
                                CI[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetLeft]; 
                            }
                            if (Dptr) {
                                CD[(k%2)*fLen+offset] = CD[((k+1)%2)*fLen+offsetUp];
                            } else {
                                CD[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetUp];
                            }
                            if (ptr == 0) {
                                CS[(k%3)*fLen+offset] = CS[((k+1)%3)*fLen+offsetDiag];
                            } else if (ptr == 1) {
                                CS[(k%3)*fLen+offset] = CI[(k%2)*fLen+offset];
                            } else {
                                CS[(k%3)*fLen+offset] = CD[(k%2)*fLen+offset];
                            } 
                        }
                        if (Iptr) ptr |= 0x04; 
                        if (Dptr) ptr |= 0x08;
                    }
                    if (k <= p_marker) {
                        if (tx % 2 == 0) ptr = (ptr << 4) & 0xF0;
                        else             ptr = ptr & 0x0F;
                        ptr += __shfl_xor_sync(0xffffffff, ptr, 1);
                        if (tx < activeThread && tx%2 == 0) tb[tb_idx+rn*THREAD_NUM/2+tx/2] = ptr;
                    }
                    __syncthreads();
                    // Calculate Max
                    for (int32_t r = activeThread/2; r > 0; r >>= 1) {
                        if (tx < r) {
                            // sub_max_score_list[tx]  = (sub_max_score_list[tx+r] > sub_max_score_list[tx]) ? sub_max_score_list[tx+r] : sub_max_score_list[tx];
                            sub_max_score_list[tx] = max(sub_max_score_list[tx+r], sub_max_score_list[tx]);
                        }
                        __syncthreads();
                    }
                    if (tx == 0) max_score_list[rn] = sub_max_score_list[0];
                }
                tb_idx += tbLen;
                __syncthreads();
                
                for (int32_t r = threadRound/2; r > 0; r >>= 1) {
                    if (tx < r) {
                        // max_score_list[tx]  = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_list[tx+r] : max_score_list[tx];
                        max_score_list[tx] = max(max_score_list[tx+r], max_score_list[tx]);
                    }
                    __syncthreads();
                }
                if (max_score_prime < max_score_list[0]) {
                    max_score_prime = max_score_list[0];
                }
                if (max_score_prime >= INT16_MAX - __float2int_rn(p_scoreMat[0])) {
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d Max number (%d) larger than the maximum limit (%d)! align length = %d, length left (r,q): (%d,%d)\n",  bx, max_score_prime, (INT16_MAX - __float2int_rn(p_scoreMat[0])), alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }
                __syncthreads();
                
                

                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                while (newL <= U[k%3]) {         
                    int32_t offset = newL - L[k%3];
                    if (S[(k%3)*fLen+offset] <= -inf) newL++;
                    else break;
                }
                while (newU >= L[k%3]) {
                    int32_t offset = newU - L[k%3];
                    if (S[(k%3)*fLen+offset] <= -inf) newU--;
                    else break;
                }
                int32_t v1 = query_length-1;
                int32_t v2 = k+2-reference_length;
                int32_t v3 = newU+1;
                int32_t Lprime = max(0, v2);
                
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 

                if (tx == 0 && !global) {
                    int32_t seqBand = v1-Lprime;
                    if (seqBand <= (U[k%3]-L[k%3]+1)) reachBoundary = true;
                }
                
                
                
                
                
                if (tx == 0 && !reachBoundary) {
                    if ((!converged) && (k < reference_length + query_length - 2)) {
                        int32_t start = newL - L[k%3];
                        int32_t length = newU - newL;
                        int32_t conv_I = CI[(k%2)*fLen+start];
                        int32_t conv_D = CD[(k%2)*fLen+start];
                        int32_t conv_S = CS[(k%3)*fLen+start];
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_I != CI[(k%2)*fLen+i]){
                                conv_I = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_D != CD[(k%2)*fLen+i]){
                                conv_D = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_S != CS[(k%3)*fLen+i]){
                                conv_S = -1;
                                break;
                            } 
                        }
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                            // conv_score = max_score_prime;
                        }
                        prev_conv_s = conv_S;
                    }
                    max_score = max_score_prime;
                    // if ((converged) && (max_score > conv_score)){
                    //     conv_logic = true;
                    // }
                    if ((converged)){
                        conv_logic = true;
                    }
                }
                __syncthreads();
                last_k = k;
                if (conv_logic || reachBoundary) break;
            }
            __syncthreads();
            if (tx == 0) {
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFF;
                    tb_state = (conv_value >> 8) & 0xFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                      tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    // if (conv_logic && bx == 0 && tx == 0) printf("%d, %d, (%d, %d)\n", tile, conv_ref_idx, conv_query_idx);
                
                } 
                else {
                    if (global) {
                        if (last_tile == true) {
                            conv_query_idx = 0;
                            conv_ref_idx = 0;
                            tb_start_addr = 0;
                            tb_start_ftr = -1;
                            tb_state = 0; 
                            alnLen[bx] = -1;
                        }
                        else if (last_k < p_marker) {
                            conv_query_idx = query_length-1;
                            conv_ref_idx = reference_length-1;
                            // tb_start_addr = ftr_addr-1;
                            tb_start_addr = ftr_addr-1;
                            tb_start_ftr = last_k;
                            tb_state = 0;
                            last_tile = true;
                        }
                        else {
                            conv_query_idx = CS[(last_k%3)*fLen] & 0xFF;
                            tb_state = (CS[(last_k%3)*fLen] >> 8) & 0xFF;
                            conv_ref_idx = p_marker - conv_query_idx; 
                            conv_ref_idx -= (tb_state == 3) ? 1: 0;
                            tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                            tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                              tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);

                            tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                        }
                    }
                    else {
                        conv_query_idx = 0;
                        conv_ref_idx = 0;
                        tb_start_addr = 0;
                        tb_start_ftr = -1;
                        tb_state = 0; 
                        last_tile = true;
                        if (tile == 0) {
                            alnLen[bx] = -1;
                            printf("%d: Cannot find any convergence point\n", bx);
                        }
                    }
                }
                int32_t addr = tb_start_addr; 
                int32_t ftr = tb_start_ftr;
                int32_t traceback_idx = conv_query_idx;
                int32_t qry_idx = conv_query_idx;
                int32_t ref_idx = conv_ref_idx;
                int8_t  tb_value = 0;
                state = tb_state%3;
                int8_t  dir = 0;
                // if(xdrop) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, qry_idx , ftr_lower_limit[ftr]);
                
                while (ftr >= 0) {
                    if (addr < 0) {
                        if (bx == 0) printf("ERROR: tb addr < 0 (%d)!\n", addr);
                        break;
                    }
                    
                    int32_t addrOffset = (traceback_idx  - ftr_lower_limit[ftr]);
                    if (addrOffset % 2 == 0)  tb_value = (tb[addr] >> 4) & 0x0F;
                    else                      tb_value = tb[addr] & 0x0F;
                    if (state == 0) { // Current State M
                        state = tb_value & 0x03;
                        if (state == 0) {
                            dir = 0;
                        }
                        else if (state == 1) {
                            dir = 1;
                            if (tb_value & 0x04) {
                                state = 1;
                            } else {
                                state = 0;
                            }   
                        }
                        else {
                            dir = 2;
                            if (tb_value & 0x08) {
                                state = 2;
                            } else {
                                state = 0;
                            }
                        }
                    } 
                    else if (state == 1) { // Current State I
                        dir = 1;
                        if (tb_value & 0x04) {
                            state = 1;
                        } else {
                            state = 0;
                        }
                    } 
                    else { // Current State D
                        dir = 2;
                        if (tb_value & 0x08) {
                            state = 2;
                        } else {
                            state = 0;
                        }
                    }
                    // addr = addr - (traceback_idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);
                    addr = addr - (traceback_idx  - ftr_lower_limit[ftr])/2 - (ftr_length[ftr - 1]);


                    if (dir == 0) {
                        // addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]);
                        addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]-1)/2;
                        ftr -= 2;
                        traceback_idx -= 1;
                        qry_idx--;
                        ref_idx--;
                        
                    }
                    else if (dir == 1) {
                        // addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]);
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]-1)/2;
                        ftr -= 1;
                        traceback_idx -=1;
                        qry_idx--;
                        
                    }
                    else {
                        // addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1);
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1-1)/2;
                        ftr -= 1;
                        ref_idx--;
                    }
                    tile_aln[aln_idx] = dir;
                    ++aln_idx;    
                    // state = next_state;
                    
                }
                state = tb_state % 3;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
            }
            __syncthreads();

            int32_t alnStartIdx = bx * 2 * seqLen + alnLen[bx];
            int txNum = (tile == 0 && includeHead) ? aln_idx : aln_idx - 1;
            // int txNum = (tile > 0) ? aln_idx - 1: aln_idx;
            if (tile == 0 && includeHead) {
                if (tx < THREAD_NUM) {
                    for (int x = tx; x < txNum; x += THREAD_NUM) aln[alnStartIdx+x] = tile_aln[aln_idx-1-x];
                }
            }
            else {
                if (tx < THREAD_NUM) {
                    for (int x = tx; x < txNum; x += THREAD_NUM) aln[alnStartIdx+x] = tile_aln[aln_idx-2-x];
                }
            }
            if ((reference_idx == refLen-1 || query_idx == qryLen-1) && global) {
                alnStartIdx = bx * 2 * seqLen + alnLen[bx];
                if (reference_idx == refLen-1 && query_idx < qryLen-1) { // dir == 1
                    for (int q = 0; q < qryLen-query_idx-1; ++q) {
                        aln[alnStartIdx+txNum+q] = 1;
                    }
                    alnLen[bx] += qryLen-query_idx-1;
                    last_tile = true;
                }
                if (query_idx == qryLen-1 && reference_idx < refLen-1) { // dir == 2
                    for (int r = 0; r < refLen-reference_idx-1; ++r) {
                        aln[alnStartIdx+txNum+r] = 2;
                    }
                    alnLen[bx] += refLen-reference_idx-1;
                    last_tile = true;
                }
            }
            if (tx == 0 && txNum > 0) {
                alnLen[bx] += txNum;
                idx[0] = reference_idx;
                idx[1] = query_idx;
            }
            __syncthreads();
            tile++;
        }
        __syncthreads();
    }
    return;

}

__global__ void alignGrpToGrp_seq(char* seqs, int8_t *aln, int32_t* len, int32_t* num, int32_t* alnLen, int32_t *seqInfo,  float* gapOpen, float* gapExtend, float* weight, float* param)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    const int fLen = FRONT_WAVE_LEN;
    const int p_marker = 128;
    const int maxDiv = FRONT_WAVE_LEN/THREAD_NUM;
    const int tile_aln_size = 2*p_marker;

    int32_t pairNum = seqInfo[0];    

    __syncthreads();
    if (bx < pairNum) {
        // Shared memory
        __shared__ int8_t tile_aln[tile_aln_size];
        __shared__ int16_t S  [3*fLen];
        __shared__ int16_t CS [3*fLen];
        __shared__ int16_t I  [2*fLen];
        __shared__ int16_t CI [2*fLen];
        __shared__ int16_t D  [2*fLen];
        __shared__ int16_t CD [2*fLen];
        __shared__ int8_t tb  [33*p_marker];
        __shared__ int32_t idx [2]; 
        __shared__ int16_t max_score_list [maxDiv];
        __shared__ int16_t sub_max_score_list [fLen/maxDiv]; 
        __shared__ int16_t max_score;
        __shared__ int16_t aln_idx;
        __shared__ bool last_tile;
        __shared__ bool converged; 
        __shared__ bool conv_logic;
        
        int32_t seqLen = seqInfo[1];
        int32_t refLen = len[2*bx];
        int32_t qryLen = len[2*bx+1];
        int32_t refNum =  num[2*bx];
        int32_t qryNum =  num[2*bx+1];

        int32_t refStart = 0;
        int32_t qryStart = 0;
        for (int i = 0; i < 2*bx; ++i) {
            refStart += num[i];
        }
        qryStart = refStart + num[2*bx];

        // float p_gapOpen = param[25];
        float p_gapExtend = param[26];
        // float p_gapClose = param[27];
        float p_xdrop = param[28];
        float p_scoreMat [25];
        float denominator = __int2float_rn(refNum) * __int2float_rn(qryNum);
        float refFreq [6];
        float qryFreq [6]; 
        for (int i = 0; i < 25; ++i) p_scoreMat[i] = param[i];
        
        int32_t tile = 0;
        int32_t last_k = 0;
        int32_t L[3], U[3];

        int32_t ftr_length [p_marker+1];
        int32_t ftr_lower_limit [p_marker+1];
        // initialize global values
        if (tx == 0) {
            last_tile = false;
        } 
        if (tx < 2) {
            idx[tx] = 0;
        }

        __syncthreads();
        
        while (!last_tile) {
            int16_t inf = p_xdrop + 1;
            int32_t reference_idx = idx[0];
            int32_t query_idx = idx[1];
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            int16_t score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            int16_t max_score_prime = -(p_xdrop+1);
            
            int32_t tb_idx = 0;
            for (int32_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (int32_t i = 0; i < p_marker+1; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }
            if (tx == 0) {
                max_score = 0;
                converged = false;
                conv_logic = false;
            }
            int32_t conv_value = 0; 
            int8_t tb_state = 0;
            int8_t ptr = 0;  // Main pointer
            int8_t state = 0;
            int32_t prev_conv_s = -1;

            bool Iptr = false;
            bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
            
            // Initialize shared memory
            if (tx < THREAD_NUM) {
                for (int i = tx; i < 3*fLen; i = i+THREAD_NUM) {                     
                    S[i] = -1; CS[i] = -1;
                    if (i < 2*fLen) { I[i] = -1; CI[i] = -1; D[i] = -1; CD[i] = -1;}
                    if (i < tile_aln_size) {tile_aln[i] = -1;}
                }
            }
            aln_idx = 0;
            __syncthreads();

            // Start alignment
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                if (tx < THREAD_NUM) {
                    sub_max_score_list[tx] = -(p_xdrop+1); 
                }
                if (tx < maxDiv) {
                    max_score_list[tx] = -(p_xdrop+1); 
                }
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d No more cells to compute based on x-drop criteria. Align length = %d, length left (r,q): (%d,%d)\n", bx, alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d Anti-diagonal larger than the maximum limit (%d)! align length = %d, length left (r,q): (%d,%d)\n",  bx,fLen, alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }
                // __syncthreads();
                ptr = 0;
                int32_t ftrLen = U[k%3] - L[k%3] + 1;
                int32_t tbLen = (((ftrLen-1) >> 1)+1);
                if (k <= p_marker) {
                    ftr_length[ftr_length_idx] = tbLen;
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    ftr_addr += tbLen;
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;      
                }
                int32_t threadRound = ftrLen / THREAD_NUM + 1;
                int32_t lastRound = ftrLen % THREAD_NUM;
                if (lastRound == 0) {
                    threadRound -= 1;
                    lastRound = THREAD_NUM;
                }
                for (int32_t rn = 0; rn < threadRound; rn += 1) {
                    int32_t activeThread = (rn != threadRound-1) ? THREAD_NUM : lastRound;
                    __syncthreads();
                    if (tx < activeThread) {
                        int32_t i=L[k%3]+rn*THREAD_NUM+tx;
                        int32_t Lprime = max(0, k-reference_length + 1);
                        int32_t j= min(k, reference_length - 1) - (i-Lprime);
                        if (j < 0) { printf("tx: %d, ERROR: j less than 0.\n", tx);}
                        int16_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            int16_t similarScore = 0;
                            float numerator = 0;
                            int32_t refFreqIdx = reference_idx + j;
                            int32_t qryFreqIdx = query_idx + i;
                            for (int l=0; l<6; l++) {
                                refFreq[l] = 0.0;
                                qryFreq[l] = 0.0;
                            }
                            // Calcualte frequency
                            for (int l = 0; l < refNum; ++l) {
                                if      (seqs[(refStart+l)*seqLen+refFreqIdx] == 'A') refFreq[0] += 1.0 * weight[refStart+l];
                                else if (seqs[(refStart+l)*seqLen+refFreqIdx] == 'C') refFreq[1] += 1.0 * weight[refStart+l];
                                else if (seqs[(refStart+l)*seqLen+refFreqIdx] == 'G') refFreq[2] += 1.0 * weight[refStart+l];
                                else if (seqs[(refStart+l)*seqLen+refFreqIdx] == 'T') refFreq[3] += 1.0 * weight[refStart+l];
                                else if (seqs[(refStart+l)*seqLen+refFreqIdx] == '-') refFreq[5] += 1.0 * weight[refStart+l];
                                else                                                  refFreq[4] += 1.0 * weight[refStart+l];
                            }
                            for (int m = 0; m < qryNum; ++m) {
                                if      (seqs[(qryStart+m)*seqLen+qryFreqIdx] == 'A') qryFreq[0] += 1.0 * weight[qryStart+m];
                                else if (seqs[(qryStart+m)*seqLen+qryFreqIdx] == 'C') qryFreq[1] += 1.0 * weight[qryStart+m];
                                else if (seqs[(qryStart+m)*seqLen+qryFreqIdx] == 'G') qryFreq[2] += 1.0 * weight[qryStart+m];
                                else if (seqs[(qryStart+m)*seqLen+qryFreqIdx] == 'T') qryFreq[3] += 1.0 * weight[qryStart+m];
                                else if (seqs[(qryStart+m)*seqLen+qryFreqIdx] == '-') qryFreq[5] += 1.0 * weight[qryStart+m];
                                else                                                  qryFreq[4] += 1.0 * weight[qryStart+m];
                            }
                            for (int l=0; l<6; l++) {
                                for (int m=0; m<6; m++) {
                                    if (m == 5 && l == 5)      numerator += 0;
                                    else if (m == 5 || l == 5) numerator = __fmaf_rn(__fmul_rn(refFreq[l], qryFreq[m]), p_gapExtend, numerator);
                                    else                       numerator = __fmaf_rn(__fmul_rn(refFreq[l], qryFreq[m]), p_scoreMat[m*5+l], numerator);
                                }
                            }
                            similarScore = __float2int_rn(__fdiv_rn(numerator, denominator));
                            if (offsetDiag < 0) match = similarScore;
                            else                match = S[(k+1)%3*fLen+offsetDiag] + similarScore;
                            
                        }

                        int16_t pos_gapOpen_ref = __float2int_rn(gapOpen[2*bx*seqLen+reference_idx+j]);
                        int16_t pos_gapOpen_qry = __float2int_rn(gapOpen[(2*bx+1)*seqLen+query_idx+i]);
                        int16_t pos_gapExtend_ref = __float2int_rn(gapExtend[2*bx*seqLen+reference_idx+j]);
                        int16_t pos_gapExtend_qry = __float2int_rn(gapExtend[(2*bx+1)*seqLen+query_idx+i]);
                        
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            delOp =  S[(k+2)%3*fLen+offsetUp] + pos_gapOpen_ref;
                            delExt = D[(k+1)%2*fLen+offsetUp] + pos_gapExtend_ref;
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + pos_gapOpen_qry;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + pos_gapExtend_qry;
                        }
                        int16_t tempI, tempD, tempH;

                        unsigned Ext = (delExt << 16) | (insExt & 0xFFFF);
                        unsigned Op =  (delOp << 16)  | (insOp & 0xFFFF);
                        unsigned ExtOp = __vibmax_s16x2(Ext, Op, &Dptr, &Iptr);
                        tempD = (ExtOp >> 16) & 0xFFFF;
                        tempI = ExtOp & 0xFFFF;
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;

                        
                        int32_t mat32 = (match << 16) | (0 & 0xFFFF);
                        int32_t ins32 = (tempI << 16) | (1 & 0xFFFF);
                        int32_t del32 = (tempD << 16) | (2 & 0xFFFF);
                        int32_t max32 = __vimax3_s32(mat32, ins32, del32);
                        tempH = (max32 >> 16) & 0xFFFF;
                        ptr = (max32) & 0xFF;

                        if (tempH < max_score-p_xdrop) {
                            tempH = -inf;
                        }
                        S[(k%3)*fLen+offset] = tempH;
                        score = tempH;
                        sub_max_score_list[tx] = score;

                        // Convergence algorithm
                        if (k == p_marker - 1) { 
                            CS[(k%3)*fLen+offset] = (3 << 8) | (i & 0xFF); 
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 8) | (i & 0xFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 8) | (i & 0xFF);
                            CD[(k%2)*fLen+offset] = (2 << 8) | (i & 0xFF);
                        } 
                        else if (k >= p_marker + 1){
                            if (Iptr) {
                                CI[(k%2)*fLen+offset] = CI[((k+1)%2)*fLen+offsetLeft]; 
                            } else {
                                CI[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetLeft]; 
                            }
                            if (Dptr) {
                                CD[(k%2)*fLen+offset] = CD[((k+1)%2)*fLen+offsetUp];
                            } else {
                                CD[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetUp];
                            }
                            if (ptr == 0) {
                                CS[(k%3)*fLen+offset] = CS[((k+1)%3)*fLen+offsetDiag];
                            } else if (ptr == 1) {
                                CS[(k%3)*fLen+offset] = CI[(k%2)*fLen+offset];
                            } else {
                                CS[(k%3)*fLen+offset] = CD[(k%2)*fLen+offset];
                            } 
                        }
                        if (Iptr) ptr |= 0x04; 
                        if (Dptr) ptr |= 0x08;
                    }
                    // Use 4 bits to store a ptr
                    if (k <= p_marker) {
                        if (tx % 2 == 0) ptr = (ptr << 4) & 0xF0;
                        else             ptr = ptr & 0x0F;
                        ptr += __shfl_xor_sync(0xffffffff, ptr, 1);
                        if (tx < activeThread && tx%2 == 0) tb[tb_idx+rn*THREAD_NUM/2+tx/2] = ptr;
                    }
                    __syncthreads();
                    // Calculate Max
                    for (int32_t r = activeThread/2; r > 0; r >>= 1) {
                        if (tx < r) {
                            sub_max_score_list[tx]  = (sub_max_score_list[tx+r] > sub_max_score_list[tx]) ? sub_max_score_list[tx+r] : sub_max_score_list[tx];
                        }
                        __syncthreads();
                    }
                    if (tx == 0) max_score_list[rn] = sub_max_score_list[0];
                }
                tb_idx += tbLen;
                __syncthreads();
                // Calculate Max
                for (int32_t r = threadRound/2; r > 0; r >>= 1) {
                    if (tx < r) {
                        max_score_list[tx]  = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_list[tx+r] : max_score_list[tx];
                    }
                    __syncthreads();
                }
                if (max_score_prime < max_score_list[0]) {
                    max_score_prime = max_score_list[0];
                }
                if (max_score_prime >= INT16_MAX - __float2int_rn(p_scoreMat[0])) {
                    if (tx == 0) {
                        last_tile = true;
                        // printf("No.%d Max number (%d) larger than the maximum limit (%d)! align length = %d, length left (r,q): (%d,%d)\n",  bx, max_score_prime, (INT16_MAX - __float2int_rn(p_scoreMat[0])), alnLen[bx], reference_length, query_length);
                    }
                    __syncthreads();
                    break;
                }
                
                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                while (newL <= U[k%3]) {         
                    int32_t offset = newL - L[k%3];
                    if (S[(k%3)*fLen+offset] <= -inf) newL++;
                    else break;
                }
                while (newU >= L[k%3]) {
                    int32_t offset = newU - L[k%3];
                    if (S[(k%3)*fLen+offset] <= -inf) newU--;
                    else break;
                }
                int32_t v1 = query_length-1;
                int32_t v2 = k+2-reference_length;
                int32_t v3 = newU+1;
                int32_t Lprime = max(0, v2);
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 
                __syncthreads();
                if (tx == 0) {
                    if ((!converged) && (k < reference_length + query_length - 2)) {
                        int32_t start = newL - L[k%3];
                        int32_t length = newU - newL;
                        int32_t conv_I = CI[(k%2)*fLen+start];
                        int32_t conv_D = CD[(k%2)*fLen+start];
                        int32_t conv_S = CS[(k%3)*fLen+start];
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_I != CI[(k%2)*fLen+i]){
                                conv_I = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_D != CD[(k%2)*fLen+i]){
                                conv_D = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_S != CS[(k%3)*fLen+i]){
                                conv_S = -1;
                                break;
                            } 
                        }
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                        }
                        prev_conv_s = conv_S;
                    }
                    max_score = max_score_prime;
                    if (converged) {
                        conv_logic = true;
                    }
                }
                __syncthreads();
                last_k = k;
                if (conv_logic) break;
            }
            if (tx == 0) {
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFF;
                    tb_state = (conv_value >> 8) & 0xFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                      tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                } 
                else {
                    if (last_tile == true) { // Met the condition of fallback to cpu
                        conv_query_idx = 0;
                        conv_ref_idx = 0;
                        tb_start_addr = 0;
                        tb_start_ftr = -1;
                        tb_state = 0; 
                        alnLen[bx] = -1;
                    }
                    else if (last_k < p_marker) {
                        conv_query_idx = query_length-1;
                        conv_ref_idx = reference_length-1;
                        tb_start_addr = ftr_addr-1;
                        tb_start_ftr = last_k;
                        tb_state = 0;
                        last_tile = true;
                    }
                    else {
                        conv_query_idx = CS[(last_k%3)*fLen] & 0xFF;
                        tb_state = (CS[(last_k%3)*fLen] >> 8) & 0xFF;
                        conv_ref_idx = p_marker - conv_query_idx; 
                        conv_ref_idx -= (tb_state == 3) ? 1: 0;
                        tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                        tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                          tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);

                        tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    }
                }
                int32_t addr = tb_start_addr; 
                int32_t ftr = tb_start_ftr;
                int32_t traceback_idx = conv_query_idx;
                int32_t qry_idx = conv_query_idx;
                int32_t ref_idx = conv_ref_idx;
                int8_t  tb_value = 0;
                state = tb_state%3;
                int8_t  dir = 0;
                while (ftr >= 0) {
                    if (addr < 0) {
                        if (bx == 0) printf("ERROR: tb addr < 0 (%d)!\n", addr);
                        break;
                    }
                    
                    int32_t addrOffset = (traceback_idx  - ftr_lower_limit[ftr]);
                    if (addrOffset % 2 == 0)  tb_value = (tb[addr] >> 4) & 0x0F;
                    else                      tb_value = tb[addr] & 0x0F;
                    if (state == 0) { // Current State M
                        state = tb_value & 0x03;
                        if (state == 0) {
                            dir = 0;
                        }
                        else if (state == 1) {
                            dir = 1;
                            if (tb_value & 0x04) {
                                state = 1;
                            } else {
                                state = 0;
                            }   
                        }
                        else {
                            dir = 2;
                            if (tb_value & 0x08) {
                                state = 2;
                            } else {
                                state = 0;
                            }
                        }
                    } 
                    else if (state == 1) { // Current State I
                        dir = 1;
                        if (tb_value & 0x04) {
                            state = 1;
                        } else {
                            state = 0;
                        }
                    } 
                    else { // Current State D
                        dir = 2;
                        if (tb_value & 0x08) {
                            state = 2;
                        } else {
                            state = 0;
                        }
                    }
                    addr = addr - (traceback_idx  - ftr_lower_limit[ftr])/2 - (ftr_length[ftr - 1]);
                    if (dir == 0) {
                        addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]-1)/2;
                        ftr -= 2;
                        traceback_idx -= 1;
                        qry_idx--;
                        ref_idx--;
                        
                    }
                    else if (dir == 1) {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]-1)/2;
                        ftr -= 1;
                        traceback_idx -=1;
                        qry_idx--;
                        
                    }
                    else {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1-1)/2;
                        ftr -= 1;
                        ref_idx--;
                    }
                    tile_aln[aln_idx] = dir;
                    ++aln_idx;     
                }
                state = tb_state % 3;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
            }
            __syncthreads();

            int32_t alnStartIdx = bx * 2 * seqLen + alnLen[bx];
            int txNum = (tile > 0) ? aln_idx - 1: aln_idx;
            if (tile == 0) {
                if (tx < THREAD_NUM) {
                    for (int x = tx; x < txNum; x += THREAD_NUM) aln[alnStartIdx+x] = tile_aln[aln_idx-1-x];
                }
            }
            else {
                if (tx < THREAD_NUM) {
                    for (int x = tx; x < txNum; x += THREAD_NUM) aln[alnStartIdx+tx] = tile_aln[aln_idx-2-tx];
                }
            }
            if (tx == 0) {
                if ((reference_idx == refLen-1 || query_idx == qryLen-1)) {
                    alnStartIdx = bx * 2 * seqLen + alnLen[bx];
                    if (reference_idx == refLen-1 && query_idx < qryLen-1) { // dir == 1
                        for (int q = 0; q < qryLen-query_idx-1; ++q) {
                            aln[alnStartIdx+txNum+q] = 1;
                        }
                        alnLen[bx] += qryLen-query_idx-1;
                        last_tile = true;
                    }
                    if (query_idx == qryLen-1 && reference_idx < refLen-1) { // dir == 2
                        for (int r = 0; r < refLen-reference_idx-1; ++r) {
                            aln[alnStartIdx+txNum+r] = 2;
                        }
                        alnLen[bx] += refLen-reference_idx-1;
                        last_tile = true;
                    }
                }
                alnLen[bx] += txNum;
                idx[0] = reference_idx;
                idx[1] = query_idx;
            }
            __syncthreads();
            tile++;
        }
        __syncthreads();
    }
    return;

}