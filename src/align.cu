#ifndef ALIGN_HPP
#include "align.cuh"
#define INF (1<<20)
#endif

/*
__global__ void alignGrpToGrp_talco(uint16_t* freq, int8_t *aln, int32_t* len, int32_t* alnLen, int32_t *seqInfo, paramType* param)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    // int bs = blockDim.x;
    // int gs = gridDim.x;
    // int tidx = bx*bs+tx;

    const int fLen = FRONT_WAVE_LEN; // frontier length (assuming anti-diagonal length cannot exceed 1024)
    const int p_marker = 128;
    const int maxDiv = FRONT_WAVE_LEN/THREAD_NUM;
    
    // int32_t threadNum = seqInfo[3];
    int32_t pairNum = seqInfo[2];    

    __syncthreads();
    if (bx < pairNum) {
    // if (bx == 0) {
        __shared__ int8_t tile_aln[13*fLen/10];
        __shared__ int16_t S  [3*fLen];
        __shared__ int32_t CS [3*fLen];
        __shared__ int16_t I  [2*fLen];
        __shared__ int32_t CI [2*fLen];
        __shared__ int16_t D  [2*fLen];
        __shared__ int32_t CD [2*fLen];
        __shared__ int8_t tb  [33*128];
        __shared__ uint32_t idx [2]; 
        // [0] reference_idx
        // [1] query_idx

        // 1024*(2+4)*7=43008 (6144)
        // 33*128 = 4224 (1920)
        // 1.5*1024 = 1536 (384)
        // 4*2 = 8 (376)


        __shared__ int16_t max_score_list [maxDiv];
        __shared__ int16_t sub_max_score_list [fLen/maxDiv]; 
        __shared__ int16_t max_score;
        __shared__ bool last_tile;
        __shared__ bool xdrop;
        __shared__ bool converged; 
        __shared__ bool conv_logic;
        
        int32_t seqLen = seqInfo[0];
        int32_t refLen = len[2*bx];
        int32_t qryLen = len[2*bx+1];

        paramType p_match = param[0];
        paramType p_mismatch = param[1];
        paramType p_gapOpen = param[2];
        paramType p_gapExtend = param[3];
        paramType p_xdrop = param[4];
        // paramType p_marker = param[5];
        paramType p_trans = param[6];
        // if (bx == 0 && tx == 0) printf("%d, %d, %d, %d, %d, %d, %d\n", p_match, p_mismatch, p_trans, p_gapOpen, p_gapExtend, p_xdrop, p_marker);

        int16_t tile = 0;
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
            int16_t reference_idx = idx[0];
            int16_t query_idx = idx[1];
            // if (tx == 0) printf("No. %d, (r, q)=(%d, %d)\n", bx, reference_idx, query_idx);
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            // int32_t score = 0; 
            int16_t score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            // int32_t max_score_prime = -(p_xdrop+1);
            int16_t max_score_prime = -(p_xdrop+1);
            // int32_t max_score_start_addr; 
            // int32_t max_score_start_ftr;
            
            // int32_t max_score_ref_idx;    
            // int32_t max_score_query_idx;
            int16_t tb_idx = 0;
            for (size_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (size_t i = 0; i < 2*p_marker; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }
            if (tx == 0) {
                max_score = 0;
                converged = false;
                conv_logic = false;
                xdrop = false;
            }
            // __syncthreads();
            int32_t conv_score = 0; 
            int32_t conv_value = 0; 
            
            int8_t tb_state = 0;
            int8_t ptr = 0;  // Main pointer
            int8_t state = 0;
            int16_t aln_idx = 0;

            
            int32_t prev_conv_s = -1;

            bool Iptr = false;
            bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
            // if (tx == 0) printf("Tile: %d, refIdx: %d, qryIdx: %d, refLen: %d, qryLen: %d\n", tile, reference_idx, query_idx, reference_length, query_length);
            __syncthreads();
            
            
            // Initialize shared memory
            if (tx < fLen) {
                for (int i = tx; i < 3*fLen; i = i+fLen) {                     
                    S[i] = -1;
                    CS[i] = -1;
                    if (i < 2*fLen) {
                        I[i] = -1;
                        CI[i] = -1;
                        D[i] = -1;
                        CD[i] = -1;
                        tile_aln[i] = -1;
                    }
                }
            }
            __syncthreads();

            
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                if (tx < fLen/maxDiv) {
                    sub_max_score_list[tx] = -(p_xdrop+1); 
                }
                __syncthreads();
                if (tx < maxDiv) {
                    max_score_list[tx] = -(p_xdrop+1); 
                }
                __syncthreads();
                // if (bx == 45 && tx == 0) printf("tile: %d, k: %d, len: %d, L: %d, U: %d, last_k: %d\n", tile, k, alnLen[bx], L[k%3]+1, U[k%3], last_k);
                // if (tidx == 0) printf("ref+qry: %d\n", reference_length + query_length - 1);
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    if (tx == 0) {
                        float threshold = 0.9;
                        int minLen = min(refLen, qryLen);
                        printf("No.%d No more cells to compute based on x-drop critieria, align length = %d, k = %d, threshold: %d\n", bx, alnLen[bx], k, static_cast<int32_t>(threshold*minLen));
                        if (alnLen[bx] < static_cast<int32_t>(threshold*minLen)) {
                            xdrop = true;
                            last_tile = true;
                        }
                    }
                    __syncthreads();
                    break;
                }
                // if (bx == 45 && tx == 0) printf("tile: %d, k: %d, len: %d, L: %d, U: %d, last_k: %d\n", tile, k, alnLen[bx], L[k%3]+1, U[k%3], last_k);
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    if (tx == 0) printf("No.%d ERROR: anti-diagonal larger than the max limit!\n", bx);
                    if (tx == 0) last_tile = true;
                    __syncthreads();
                    break;
                }
                
                // __syncthreads();
                ptr = 0;
                int ftrLen = U[k%3] - L[k%3] + 1;
                int tbLen = ((ftrLen-1)/2+1);
                if (k <= p_marker) {
                    ftr_length[ftr_length_idx] = tbLen;
                    // ftr_length[ftr_length_idx] = (U[k%3] - L[k%3] + 1);
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    // ftr_addr += U[k%3] - L[k%3] + 1;
                    ftr_addr += tbLen;
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;      
                }
                __syncthreads();
                int threadRound = ftrLen / THREAD_NUM + 1;
                int lastRound = ftrLen % THREAD_NUM;
                if (lastRound == 0) {
                    threadRound -= 1;
                    lastRound = THREAD_NUM;
                }
                // if (tx < ftrLen) {
                // for (int activeThread = tx; activeThread < ftrLen; activeThread += threadNum) {
                for (int rn = 0; rn < threadRound; rn += 1) {
                    int activeThread = (rn != threadRound-1) ? THREAD_NUM : lastRound;
                    if (tx < activeThread) {
                        int16_t i=L[k%3]+rn*THREAD_NUM+tx;
                        int16_t Lprime = max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1);
                        int16_t j= min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime);
                        if (j < 0) { printf("tx: %d, ERROR: j less than 0.\n", tx);}
                        int16_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        int16_t score_from_prev_tile = 0;
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            if (k==0 && tile>0) {
                                score_from_prev_tile = tile*10;
                            }
                            int16_t similarScore = 0;
                            float denominator = 0;
                            float numerator = 0;
                            int refFreqStart = 6*(2*bx)*seqLen;
                            int qryFreqStart = 6*(2*bx+1)*seqLen;
                            int refFreqIdx = reference_idx + j;
                            int qryFreqIdx = query_idx + i;
                            for (int l=0; l<6; l++) {
                                for (int m=0; m<6; m++) {
                                    denominator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m];
                                    if ((m == 4 || l == 4) || (m == 5 && l == 5)) numerator += 0;
                                    else if (m == 5 || l == 5)                    numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_gapExtend;
                                    else if (m == l)                              numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_match;
                                    else if (abs(m-l) == 2)                       numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_trans;
                                    else                                          numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_mismatch;
                                }
                            }
                            similarScore = static_cast<int16_t>(roundf(numerator/denominator));
                            if (offsetDiag < 0) match = similarScore + score_from_prev_tile;
                            else                match = S[(k+1)%3*fLen+offsetDiag] + similarScore + score_from_prev_tile;
                        }
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            delOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            delExt = D[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                        }
                        int16_t tempI, tempD, tempH;

                        unsigned Ext = static_cast<uint32_t>((delExt) << 16 | (insExt & 0xFFFF));
                        unsigned Op = static_cast<uint32_t>((delOp) << 16  | (insOp & 0xFFFF));
                        unsigned int ExtOp = __vibmax_s16x2(Ext, Op, &Dptr, &Iptr);
                        tempD = (ExtOp >> 16) & 0xFFFF;
                        tempI = ExtOp & 0xFFFF;

                        // tempI =  insOp;
                        // tempD =  delOp;
                        // Iptr = false;
                        // Dptr = false;
                        // if (insExt >= insOp) {
                        //     tempI = insExt;
                        //     Iptr = true;
                        // }
                        // if (delExt >= delOp) {
                        //     tempD = delExt;
                        //     Dptr = true;
                        // }
                        int32_t mat32 = ((match) << 16 | (0 & 0xFFFF));
                        int32_t ins32 = ((tempI) << 16 | (1 & 0xFFFF));
                        int32_t del32 = ((tempD) << 16 | (2 & 0xFFFF));
                        int32_t max32 = __vimax3_s32(mat32, ins32, del32);
                        tempH = (max32 >> 16) & 0xFFFF;
                        ptr = (max32) & 0xFF;
                        // if (match > tempI) {
                        //     if (match > tempD) {
                        //         tempH = match;
                        //         ptr = 0;
                        //     }
                        //     else {
                        //         tempH = tempD;
                        //         ptr = 2;
                        //     }
                        // }
                        // else if (tempI > tempD) {
                        //     tempH = tempI;
                        //     ptr = 1;
                        // }
                        // else {
                        //     tempH = tempD;
                        //     ptr = 2;
                        // }
                        if (tempH < max_score-p_xdrop) {
                            tempH = -inf;
                        }
                        // if(bx == 0 && tx == (U[k%3]-L[k%3])/2 && tile == 0) printf("k: %d, (%d, %d), idx: %d, H: %d, D: %d, I: %d, ptr: %d\n", k, i, j, tx, tempH, tempD, tempI, ptr);
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;
                        S[(k%3)*fLen+offset] = tempH;
                        score = tempH;


                        sub_max_score_list[tx] = score;

                        if (k == p_marker - 1) { // Convergence algorithm
                            CS[(k%3)*fLen+offset] = (3 << 16) | (i & 0xFFFF); 
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 16) | (i & 0xFFFF);
                            CD[(k%2)*fLen+offset] = (2 << 16) | (i & 0xFFFF);
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
                    for (uint32_t r = activeThread/2; r > 0; r >>= 1) {
                        if (tx < r) {
                            sub_max_score_list[tx]  = (sub_max_score_list[tx+r] > sub_max_score_list[tx]) ? sub_max_score_list[tx+r] : sub_max_score_list[tx];
                        }
                        __syncthreads();
                    }
                    if (tx == 0) max_score_list[rn] = sub_max_score_list[0];
                }
                tb_idx += tbLen;
                __syncthreads();
                
                for (uint32_t r = threadRound/2; r > 0; r >>= 1) {
                    if (tx < r) {
                        max_score_list[tx]  = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_list[tx+r] : max_score_list[tx];
                    }
                    __syncthreads();
                }
                if (max_score_prime < max_score_list[0]) {
                    max_score_prime = max_score_list[0];
                }
                __syncthreads();
                
                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                // if (tx == 0) printf("k: %d,newL: %d, newU: %d\n",k, L[k%3], U[k%3]);
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
                int32_t v1 = static_cast<int32_t>(query_length)-1;
                int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                int32_t v3 = newU+1;
                int32_t Lprime = max(static_cast<int32_t>(0), v2);
                // __syncthreads();
                // if (tx == 0) printf("newU: %d, + %d = %d ", newU, 1, v3);
                // if (tx == 0) printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                // if (tx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
        
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 
                // if (tx == 0 && bx == 42 && tile == 0) printf("DEBIG4\n");
                // if (tidx == 0) printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
                // sL[(k+1)%3] = L[(k+1)%3];
                // sU[(k+1)%3] = U[(k+1)%3];

                if (tx == 0) {
                    if ((!converged) && (k < reference_length + query_length - 2)) {
                        int32_t start = newL - L[k%3];
                        int32_t length = newU - newL;
                        // if(bx == 1 && tx == 0 && tile >= 462 )printf("k: %d, st: %d, l: %d\n",k, start, length);
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
                        // if(bx == 4 && tx == 0 && tile == 100) printf("k: %d, CONV: %d, %d, %d\n", k, conv_I, conv_D, conv_S);
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                            conv_score = max_score_prime;
                            // if (DEBUG)  std::cout << "Converged at: " << conv_value << "\n";
                            // if (bx == 0) printf("Converged at: %d\n", conv_value);
                        }
                        prev_conv_s = conv_S;
                    }
                    // Update max_score
                    max_score = max_score_prime;
                    // printf("Max score: %d\n", max_score);
                    if ((converged) && (max_score > conv_score)){
                        conv_logic = true;
                    }
                }
                __syncthreads();
                // if (tx == 0 && bx == 42 && tile == 0) printf("DEBIG5\n");
                last_k = k;
                if (conv_logic) break;
            }
            __syncthreads();
            if (tx == 0) {
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFFFF;
                    tb_state = (conv_value >> 16) & 0xFFFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    // tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                    //                                   tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                      tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                } 
                else {
                    if (last_tile == true) {
                        conv_query_idx = 0;
                        conv_ref_idx = 0;
                        tb_start_addr = 0;
                        tb_start_ftr = -1;
                        tb_state = 0; 
                        // alnLen[bx] = 0;
                        if (!xdrop) alnLen[bx] = 0;
                        else        alnLen[bx] = -100;
                        // if (bx == 45) printf("op 1\n");
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
                        // if (bx == 45) printf("op 3 last_k %d\n", last_k);
                        conv_query_idx = CS[(last_k%3)*fLen] & 0xFFFF;
                        tb_state = (CS[(last_k%3)*fLen] >> 16) & 0xFFFF;
                        conv_ref_idx = p_marker - conv_query_idx; 
                        conv_ref_idx -= (tb_state == 3) ? 1: 0;
                        tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                        // tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                        //                                   tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                        tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                          tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);

                        tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    }
                }
                int32_t addr = tb_start_addr; 
                int16_t ftr = tb_start_ftr;
                int16_t traceback_idx = conv_query_idx;
                int16_t qry_idx = conv_query_idx;
                int16_t ref_idx = conv_ref_idx;
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
                    // tb_value = tb[addr];
                    // tb_value = (tb[addr] >> 4) & 0x0F;
                    // if (bx == 45) printf("ERROR: After tb addr (%d)!\n", addr);
                    // if (bx == 0 && tx == 0) printf("addr: %d, ftr: %d, dir: %d, (r,q): (%d, %d), tb: %d, off: %d\n", addr, ftr, (dir&0xFFFF), ref_idx, qry_idx, tb[addr] & 0x00FF, addrOffset);
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
                    // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
                }
                
                state = tb_state % 3;

                // Write global memory
                // int32_t refIndex = reference_idx;
                // int32_t qryIndex = query_idx;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
            }
            __syncthreads();
            int32_t numRound = aln_idx / THREAD_NUM + 1;
            int32_t lastRound = aln_idx % THREAD_NUM;
            int32_t alnStartIdx = bx * 2 * seqLen + alnLen[bx];
            for (int round = 0; round < numRound; ++round) {
                int txNum = (round == numRound - 1) ? lastRound : THREAD_NUM;
                if (tile > 0) txNum -= 1;
                // globalAlnIdx += txNum;
                // if (bx == 4 && tx == 0) {
                //     printf("Tile: %d, tbLen: %d, startIdx: %d, %d\n", tile, txNum, alnStartIdx, alnStartIdx-bx*2*seqLen);
                //     for (int x = 0; x < txNum; ++x) {
                //         if (tile == 0) printf("%d,", (tile_aln[aln_idx-1-threadNum*round-x] & 0xFFFF));
                //         else           printf("%d,", (tile_aln[aln_idx-2-threadNum*round-x] & 0xFFFF));
                //     }
                //     printf("\n");
                // }
                if (tx == 0) {
                    if (tile == 0) {
                        for (int x = 0; x < txNum; ++x) {
                            aln[alnStartIdx+THREAD_NUM*round+x] = tile_aln[aln_idx-1-THREAD_NUM*round-x];
                        }
                    }
                    else {
                        for (int x = 0; x < txNum; ++x) {
                            aln[alnStartIdx+THREAD_NUM*round+x] = tile_aln[aln_idx-2-THREAD_NUM*round-x];
                        }
                    }
                }
                // if (tx < txNum) {
                //     if (tile == 0) aln[alnStartIdx+threadNum*round+tx] = tile_aln[aln_idx-1-threadNum*round-tx];
                //     else           aln[alnStartIdx+threadNum*round+tx] = tile_aln[aln_idx-2-threadNum*round-tx];
                // }
                // if (bx == 4 && tx == 0) {
                //     for (int x = 0; x < txNum; ++x) {
                //         printf("%d,", (aln[alnStartIdx+threadNum*round+x] & 0xFFFF));
                //     }
                //     printf("\n");
                // }
                __syncthreads();
                if (tx == 0) {
                    alnLen[bx] += txNum;
                }
            }
            if (tx == 0) {
                idx[0] = reference_idx;
                idx[1] = query_idx;
            }
            __syncthreads();
            tile++;
        }
        __syncthreads();
        // __syncthreads();
        // TODO: Add score to debug
        // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
    }
    return;

}
*/


__global__ void alignGrpToGrp_talco(uint16_t* freq, int8_t *aln, int32_t* len, int32_t* alnLen, int32_t *seqInfo, paramType* param)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;

    const int fLen = FRONT_WAVE_LEN; // frontier length (assuming anti-diagonal length cannot exceed 1024)
    const int p_marker = 128;
    const int maxDiv = FRONT_WAVE_LEN/THREAD_NUM;
    
    // int32_t threadNum = seqInfo[3];
    int32_t pairNum = seqInfo[2];    

    __syncthreads();
    if (bx < pairNum) {
    // if (bx == 0) {
        __shared__ int8_t tile_aln[13*fLen/10];
        __shared__ int16_t S  [3*fLen];
        __shared__ int32_t CS [3*fLen];
        __shared__ int16_t I  [2*fLen];
        __shared__ int32_t CI [2*fLen];
        __shared__ int16_t D  [2*fLen];
        __shared__ int32_t CD [2*fLen];
        __shared__ int8_t tb  [33*128];
        __shared__ uint32_t idx [2]; 
        // [0] reference_idx
        // [1] query_idx
        // 1024*(2+4)*7=43008 (6144)
        // 33*128 = 4224 (1920)
        // 1.5*1024 = 1536 (384)
        // 4*2 = 8 (376)


        __shared__ int16_t max_score_list [maxDiv];
        __shared__ int16_t sub_max_score_list [fLen/maxDiv]; 
        __shared__ int16_t max_score;
        __shared__ bool last_tile;
        __shared__ bool xdrop;
        __shared__ bool converged; 
        __shared__ bool conv_logic;
        
        int32_t seqLen = seqInfo[0];
        // int32_t scoreMode = seqInfo[4];

        int32_t refLen = len[2*bx];
        int32_t qryLen = len[2*bx+1];

        paramType p_gapOpen = param[25];
        paramType p_gapExtend = param[26];
        paramType p_xdrop = param[27];

        paramType p_scoreMat [25];
        for (int i = 0; i < 25; ++i) p_scoreMat[i] = param[i];
        
        int16_t tile = 0;
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
            int16_t reference_idx = idx[0];
            int16_t query_idx = idx[1];
            // if (tx == 0) printf("No. %d, (r, q)=(%d, %d)\n", bx, reference_idx, query_idx);
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            // int32_t score = 0; 
            int16_t score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            // int32_t max_score_prime = -(p_xdrop+1);
            int16_t max_score_prime = -(p_xdrop+1);
            
            int16_t tb_idx = 0;
            for (size_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (size_t i = 0; i < 2*p_marker; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }
            if (tx == 0) {
                max_score = 0;
                converged = false;
                conv_logic = false;
                xdrop = false;
            }
            // __syncthreads();
            int32_t conv_score = 0; 
            int32_t conv_value = 0; 
            
            int8_t tb_state = 0;
            int8_t ptr = 0;  // Main pointer
            int8_t state = 0;
            int16_t aln_idx = 0;

            
            int32_t prev_conv_s = -1;

            bool Iptr = false;
            bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
            // if (tx == 0) printf("Tile: %d, refIdx: %d, qryIdx: %d, refLen: %d, qryLen: %d\n", tile, reference_idx, query_idx, reference_length, query_length);
            __syncthreads();
            
            
            // Initialize shared memory
            if (tx < fLen) {
                for (int i = tx; i < 3*fLen; i = i+fLen) {                     
                    S[i] = -1;
                    CS[i] = -1;
                    if (i < 2*fLen) {
                        I[i] = -1;
                        CI[i] = -1;
                        D[i] = -1;
                        CD[i] = -1;
                        tile_aln[i] = -1;
                    }
                }
            }
            __syncthreads();

            
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
                        float threshold = 0.9;
                        int minLen = min(refLen, qryLen);
                        if (alnLen[bx] < static_cast<int32_t>(threshold*minLen)) {
                            xdrop = true;
                            last_tile = true;
                            printf("No.%d No more cells to compute based on x-drop criteria, align length = %d, k = %d, threshold: %d\n", bx, alnLen[bx], k, static_cast<int32_t>(threshold*minLen));
                        }
                        // last_tile = true;
                        // printf("No.%d No more cells to compute based on x-drop criteria, align length = %d, k = %d, minSeqLen: %d\n", bx, alnLen[bx], k, static_cast<int32_t>(minLen));
                    }
                    __syncthreads();
                    break;
                }
                
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    if (tx == 0) printf("No.%d ERROR: anti-diagonal larger than the max limit! align length = %d, k = %d\n", bx, alnLen[bx], k);
                    if (tx == 0) last_tile = true;
                    __syncthreads();
                    break;
                }
                
                // __syncthreads();
                ptr = 0;
                int ftrLen = U[k%3] - L[k%3] + 1;
                int tbLen = ((ftrLen-1)/2+1);
                if (k <= p_marker) {
                    ftr_length[ftr_length_idx] = tbLen;
                    // ftr_length[ftr_length_idx] = (U[k%3] - L[k%3] + 1);
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    // ftr_addr += U[k%3] - L[k%3] + 1;
                    ftr_addr += tbLen;
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;      
                }
                __syncthreads();
                int threadRound = ftrLen / THREAD_NUM + 1;
                int lastRound = ftrLen % THREAD_NUM;
                if (lastRound == 0) {
                    threadRound -= 1;
                    lastRound = THREAD_NUM;
                }
                // if (tx < ftrLen) {
                // for (int activeThread = tx; activeThread < ftrLen; activeThread += threadNum) {
                for (int rn = 0; rn < threadRound; rn += 1) {
                    int activeThread = (rn != threadRound-1) ? THREAD_NUM : lastRound;
                    if (tx < activeThread) {
                        int16_t i=L[k%3]+rn*THREAD_NUM+tx;
                        int16_t Lprime = max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1);
                        int16_t j= min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime);
                        if (j < 0) { printf("tx: %d, ERROR: j less than 0.\n", tx);}
                        int16_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        int16_t score_from_prev_tile = 0;
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            if (k==0 && tile>0) {
                                score_from_prev_tile = tile*10;
                            }
                            int16_t similarScore = 0;
                            float denominator = 0;
                            float numerator = 0;
                            int refFreqStart = 6*(2*bx)*seqLen;
                            int qryFreqStart = 6*(2*bx+1)*seqLen;
                            int refFreqIdx = reference_idx + j;
                            int qryFreqIdx = query_idx + i;
                            for (int l=0; l<6; l++) {
                                for (int m=0; m<6; m++) {
                                    denominator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m];
                                    if (m == 5 && l == 5)      numerator += 0;
                                    else if (m == 5 || l == 5) numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_gapExtend;
                                    else                       numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_scoreMat[m*5+l];
                                }
                            }
                            similarScore = static_cast<int16_t>(roundf(numerator/denominator));
                            if (offsetDiag < 0) match = similarScore + score_from_prev_tile;
                            else                match = S[(k+1)%3*fLen+offsetDiag] + similarScore + score_from_prev_tile;
                        }
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            delOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            delExt = D[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                        }
                        int16_t tempI, tempD, tempH;

                        unsigned Ext = static_cast<uint32_t>((delExt) << 16 | (insExt & 0xFFFF));
                        unsigned Op = static_cast<uint32_t>((delOp) << 16  | (insOp & 0xFFFF));
                        unsigned int ExtOp = __vibmax_s16x2(Ext, Op, &Dptr, &Iptr);
                        tempD = (ExtOp >> 16) & 0xFFFF;
                        tempI = ExtOp & 0xFFFF;

                        int32_t mat32 = ((match) << 16 | (0 & 0xFFFF));
                        int32_t ins32 = ((tempI) << 16 | (1 & 0xFFFF));
                        int32_t del32 = ((tempD) << 16 | (2 & 0xFFFF));
                        int32_t max32 = __vimax3_s32(mat32, ins32, del32);
                        tempH = (max32 >> 16) & 0xFFFF;
                        ptr = (max32) & 0xFF;
                        if (tempH < max_score-p_xdrop) {
                            tempH = -inf;
                        }
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;
                        S[(k%3)*fLen+offset] = tempH;
                        score = tempH;


                        sub_max_score_list[tx] = score;

                        if (k == p_marker - 1) { // Convergence algorithm
                            CS[(k%3)*fLen+offset] = (3 << 16) | (i & 0xFFFF); 
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 16) | (i & 0xFFFF);
                            CD[(k%2)*fLen+offset] = (2 << 16) | (i & 0xFFFF);
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
                    for (uint32_t r = activeThread/2; r > 0; r >>= 1) {
                        if (tx < r) {
                            sub_max_score_list[tx]  = (sub_max_score_list[tx+r] > sub_max_score_list[tx]) ? sub_max_score_list[tx+r] : sub_max_score_list[tx];
                        }
                        __syncthreads();
                    }
                    if (tx == 0) max_score_list[rn] = sub_max_score_list[0];
                }
                tb_idx += tbLen;
                __syncthreads();
                
                for (uint32_t r = threadRound/2; r > 0; r >>= 1) {
                    if (tx < r) {
                        max_score_list[tx]  = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_list[tx+r] : max_score_list[tx];
                    }
                    __syncthreads();
                }
                if (max_score_prime < max_score_list[0]) {
                    max_score_prime = max_score_list[0];
                }
                __syncthreads();
                
                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                // if (tx == 0) printf("k: %d,newL: %d, newU: %d\n",k, L[k%3], U[k%3]);
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
                int32_t v1 = static_cast<int32_t>(query_length)-1;
                int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                int32_t v3 = newU+1;
                int32_t Lprime = max(static_cast<int32_t>(0), v2);
                // __syncthreads();
                // if (tx == 0) printf("newU: %d, + %d = %d ", newU, 1, v3);
                // if (tx == 0) printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                // if (tx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
        
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 
                // if (tx == 0 && bx == 42 && tile == 0) printf("DEBIG4\n");
                // if (tidx == 0) printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
                // sL[(k+1)%3] = L[(k+1)%3];
                // sU[(k+1)%3] = U[(k+1)%3];

                if (tx == 0) {
                    if ((!converged) && (k < reference_length + query_length - 2)) {
                        int32_t start = newL - L[k%3];
                        int32_t length = newU - newL;
                        // if(bx == 1 && tx == 0 && tile >= 462 )printf("k: %d, st: %d, l: %d\n",k, start, length);
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
                        // if(bx == 4 && tx == 0 && tile == 100) printf("k: %d, CONV: %d, %d, %d\n", k, conv_I, conv_D, conv_S);
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                            conv_score = max_score_prime;
                            // if (DEBUG)  std::cout << "Converged at: " << conv_value << "\n";
                            // if (bx == 0) printf("Converged at: %d\n", conv_value);
                        }
                        prev_conv_s = conv_S;
                    }
                    // Update max_score
                    max_score = max_score_prime;
                    // printf("Max score: %d\n", max_score);
                    if ((converged) && (max_score > conv_score)){
                        conv_logic = true;
                    }
                }
                __syncthreads();
                // if (tx == 0 && bx == 42 && tile == 0) printf("DEBIG5\n");
                last_k = k;
                if (conv_logic) break;
            }
            __syncthreads();
            if (tx == 0) {
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFFFF;
                    tb_state = (conv_value >> 16) & 0xFFFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    // tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                    //                                   tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                      tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                } 
                else {
                    if (last_tile == true) {
                        conv_query_idx = 0;
                        conv_ref_idx = 0;
                        tb_start_addr = 0;
                        tb_start_ftr = -1;
                        tb_state = 0; 
                        // alnLen[bx] = 0;
                        if (!xdrop) alnLen[bx] = 0;
                        else        alnLen[bx] = -100;
                        // if (bx == 45) printf("op 1\n");
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
                        // if (bx == 45) printf("op 3 last_k %d\n", last_k);
                        conv_query_idx = CS[(last_k%3)*fLen] & 0xFFFF;
                        tb_state = (CS[(last_k%3)*fLen] >> 16) & 0xFFFF;
                        conv_ref_idx = p_marker - conv_query_idx; 
                        conv_ref_idx -= (tb_state == 3) ? 1: 0;
                        tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                        // tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                        //                                   tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                        tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2])/2) : 
                                                          tb_start_addr +  ((conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]) / 2);

                        tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    }
                }
                int32_t addr = tb_start_addr; 
                int16_t ftr = tb_start_ftr;
                int16_t traceback_idx = conv_query_idx;
                int16_t qry_idx = conv_query_idx;
                int16_t ref_idx = conv_ref_idx;
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
                    // tb_value = tb[addr];
                    // tb_value = (tb[addr] >> 4) & 0x0F;
                    // if (bx == 45) printf("ERROR: After tb addr (%d)!\n", addr);
                    // if (bx == 0 && tx == 0) printf("addr: %d, ftr: %d, dir: %d, (r,q): (%d, %d), tb: %d, off: %d\n", addr, ftr, (dir&0xFFFF), ref_idx, qry_idx, tb[addr] & 0x00FF, addrOffset);
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
                    // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
                }
                
                state = tb_state % 3;

                // Write global memory
                // int32_t refIndex = reference_idx;
                // int32_t qryIndex = query_idx;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
            }
            __syncthreads();
            int32_t numRound = aln_idx / THREAD_NUM + 1;
            int32_t lastRound = aln_idx % THREAD_NUM;
            int32_t alnStartIdx = bx * 2 * seqLen + alnLen[bx];
            for (int round = 0; round < numRound; ++round) {
                int txNum = (round == numRound - 1) ? lastRound : THREAD_NUM;
                if (tile > 0) txNum -= 1;
                // globalAlnIdx += txNum;
                // if (bx == 4 && tx == 0) {
                //     printf("Tile: %d, tbLen: %d, startIdx: %d, %d\n", tile, txNum, alnStartIdx, alnStartIdx-bx*2*seqLen);
                //     for (int x = 0; x < txNum; ++x) {
                //         if (tile == 0) printf("%d,", (tile_aln[aln_idx-1-threadNum*round-x] & 0xFFFF));
                //         else           printf("%d,", (tile_aln[aln_idx-2-threadNum*round-x] & 0xFFFF));
                //     }
                //     printf("\n");
                // }
                if (tx == 0) {
                    if (tile == 0) {
                        for (int x = 0; x < txNum; ++x) {
                            aln[alnStartIdx+THREAD_NUM*round+x] = tile_aln[aln_idx-1-THREAD_NUM*round-x];
                        }
                    }
                    else {
                        for (int x = 0; x < txNum; ++x) {
                            aln[alnStartIdx+THREAD_NUM*round+x] = tile_aln[aln_idx-2-THREAD_NUM*round-x];
                        }
                    }
                }
                // if (tx < txNum) {
                //     if (tile == 0) aln[alnStartIdx+threadNum*round+tx] = tile_aln[aln_idx-1-threadNum*round-tx];
                //     else           aln[alnStartIdx+threadNum*round+tx] = tile_aln[aln_idx-2-threadNum*round-tx];
                // }
                // if (bx == 4 && tx == 0) {
                //     for (int x = 0; x < txNum; ++x) {
                //         printf("%d,", (aln[alnStartIdx+threadNum*round+x] & 0xFFFF));
                //     }
                //     printf("\n");
                // }
                __syncthreads();
                if (tx == 0) {
                    alnLen[bx] += txNum;
                }
            }
            if (tx == 0) {
                idx[0] = reference_idx;
                idx[1] = query_idx;
            }
            __syncthreads();
            tile++;
        }
        __syncthreads();
        // __syncthreads();
        // TODO: Add score to debug
        // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
    }
    return;

}




int8_t updateState(STATE& currentHState, STATE& currentIState, STATE& currentDState)
{
    int8_t currentState = 0;
    switch (currentHState)
    {
    case STATE::HI:
        currentState = (currentState | 0x01); break;
    case STATE::HD:
        currentState = (currentState | 0x02); break;
    default: 
        currentState = currentState; break;
    }

    switch (currentIState)
    {
    case STATE::II:
        currentState = (currentState | 0x04);  break;
    default:
        currentState = currentState; break;
    }
    
    switch (currentDState)
    {
    case STATE::DD:
        currentState = (currentState | 0x08); break;
    default:
        currentState = currentState; break;
    }
    return currentState;
}

void tracebackGrpToGrp (int8_t state, std::vector<int8_t>& TB, std::vector<int32_t>& wfLL, std::vector<int32_t>& wfLen, std::vector<int8_t>& aln, int32_t refLen, int32_t queryLen) {

    int32_t refIndex = refLen-1;
    int32_t queryIndex = queryLen-1;
    
    int32_t k=wfLL.size()-1;
    int32_t tbIndex = TB.size()-1;
    int8_t currentTB;
    int8_t dir;

    aln.clear();

    while (refIndex>=0 && queryIndex>=0) {
        currentTB = TB[tbIndex];
        if (state == 0) { // Current State M
            state = currentTB & 0x03;
            if (state == 0) dir = 0;
            else if (state == 1)
            {
                dir = 1;
                if (currentTB & 0x04) state = 1;
                else state = 0;   
            }
            else
            {
                dir = 2;
                if (currentTB & 0x08) state = 2;
                else state = 0;
            }
        }
        else if (state == 1) { // Current State I
            dir = 1;
            if (currentTB & 0x04) state = 1;
            else state = 0;
        } 
        else { // Current State D
            dir = 2;
            if (currentTB & 0x08) state = 2;
            else state = 0;
        }

        tbIndex -= (refIndex - wfLL[k] + 1 + wfLen[k-1]);
        if (dir == 0) {
            tbIndex -= wfLen[k-2]; tbIndex += refIndex - wfLL[k-2];
            k--;k--;queryIndex--;refIndex--;
        }
        else if (dir == 1) {
            tbIndex += (refIndex - wfLL[k-1] + 1);
            k--;queryIndex--;
        }
        else {
            tbIndex += (refIndex - wfLL[k-1]);
            k--;refIndex--;
        }
        aln.push_back(dir);   
    }


    while (refIndex>=0)
    {
        aln.push_back(static_cast<int8_t>(2));
        k--;refIndex--;
    }

    while (queryIndex>=0)
    {
        aln.push_back(static_cast<int8_t>(1));
        k--;queryIndex--;
    }
    return;
}

void alignGrpToGrp_traditional (uint16_t* freq, int32_t seqLen, int32_t refLen, int32_t queryLen, Params& param, std::vector<int8_t>& aln) {

    if (refLen<=0 || queryLen<=0) {fprintf(stderr, "Error: Ref/Query length <= 0\n"); exit(1);}
    
    int maxWFLen = 0; //wavefront length
    // if (param.marker != 0) maxWFLen = param.marker+2;
    // else maxWFLen = refLen + queryLen + 2;
    maxWFLen = refLen + queryLen + 2;

    // std::cout << charFreqRef[0].size() << "-" << charFreqQuery[0].size() << "\n";
    // int score = 0;
    int32_t *H[3], *I[2], *D[2];
    int32_t L[3], U[3];
    std::vector<int32_t> wfLL, wfLen;
    paramType scoreMat [25];
    paramType gapOpen = (param.userDefine) ? param.userGapOpen : param.gapOpen;
    paramType gapExtend = (param.userDefine) ? param.userGapExtend : param.gapExtend;
    if (!param.userDefine) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 5 || j == 5)          scoreMat[i*5+j] = 0;
                else if (i == j)               scoreMat[i*5+j] = param.match;
                else if (i-j == 2 || j-i == 2) scoreMat[i*5+j] = param.trans;
                else                           scoreMat[i*5+j] = param.mismatch;
            }
        }
    }
    else if (param.userDefine == 1) {
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) scoreMat[i*5+j] = param.userMatrix[i][j];
    }

    //Output
    int8_t state=0;
    std::vector<int8_t> TB;

    for(size_t i=0; i<3; i++)
    {
        H[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        if (i<2) 
        {
            I[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
            D[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        }
        L[i]=0;
        U[i]=0;
    }

    for (size_t i=0; i<3; i++)
    {
        for (size_t j=0; j<(size_t)maxWFLen; j++)
        {
            H[i][j] = -INF;
            if (i<2) {I[i][j] = 0; D[i][j] = -INF;} 
        }
    }

    for (int32_t k=0; k<refLen+queryLen+1; k++) {
        L[k%3] = (k<=queryLen)?0:k-queryLen;
        U[k%3] = (k<=refLen)?k:refLen;
        wfLL.push_back(L[k%3]);
        wfLen.push_back(U[k%3]-L[k%3]+1);
        for(int32_t i=L[k%3]; i<U[k%3]+1; i++) { // i->Ref Index 
            int32_t j=(k-i); //j->Query Index
            int32_t match = -INF, insOp = -INF, delOp = -INF, insExt = -INF, delExt = -INF;
            int32_t offset = i-L[k%3];
            int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1;
            int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
            int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
            STATE currentHState = STATE::HH;
            STATE currentIState = STATE::IH;
            STATE currentDState = STATE::DH;

            if (k==0) match=0;

            if (offsetDiag>=0 && j>0) {
                int32_t similarScore = 0;
                float denominator = 0;
                float numerator = 0;
                for (int l=0; l<6; l++) {
                    for (int m=0; m<6; m++) {
                        float freqMult = freq[6*(i-1)+l]*freq[6*(seqLen+j-1)+m];
                        denominator += freqMult;
                        if (m == 5 && l == 5)      numerator += 0;
                        else if (m == 5 || l == 5) numerator += freqMult*param.gapExtend;
                        else                       numerator += freqMult*scoreMat[m*5+l];
                    }
                }
                similarScore = static_cast<int32_t>(round(numerator/denominator));
                // if (i == (U[k%3]-L[k%3])/2) printf("k: %d, idx: %d, score: %d\n", k, i, similarScore);
                match = H[(k+1)%3][offsetDiag] + similarScore;
            }
            if (offsetUp >= 0) {
                insOp = H[(k+2)%3][offsetUp] + gapOpen;
                insExt = I[(k+1)%2][offsetUp] + gapExtend;
            }
            if (offsetLeft >=0) {
                delOp = H[(k+2)%3][offsetLeft] + gapOpen;
                delExt = D[(k+1)%2][offsetLeft] + gapExtend;
            }
            I[k%2][offset] =  insOp;
            D[k%2][offset] =  delOp;
            
            if (insExt >= insOp) {
                I[k%2][offset] = insExt;
                currentIState = STATE::II;
            }
            if (delExt >= delOp) {
                D[k%2][offset] = delExt;
                currentDState = STATE::DD;
            }
            if (match > I[k%2][offset]) {
                if (match > D[k%2][offset]) H[k%3][offset] = match;
                else 
                {
                    H[k%3][offset] = D[k%2][offset];
                    currentHState = STATE::HD;
                }
            }
            else if (I[k%2][offset] > D[k%2][offset]) {
                H[k%3][offset] = I[k%2][offset];
                currentHState = STATE::HI;
            }
            else {
                H[k%3][offset] = D[k%2][offset];
                currentHState = STATE::HD;
            }

            TB.push_back(updateState(currentHState, currentIState, currentDState));
            state = currentHState;
        }
    }

    // Deallocate memory for scores
    for (size_t sIndx=0; sIndx<3; sIndx++) {
        std::free(H[sIndx]);
        if (sIndx < 2) {
            std::free(I[sIndx]);
            std::free(D[sIndx]);
        }
    }

    
    tracebackGrpToGrp (state, TB, wfLL, wfLen, aln, refLen, queryLen);
    return;

}