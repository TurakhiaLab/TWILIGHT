#ifndef ALIGN_HPP
#include "align.cuh"
#define INF (1<<20)
#endif

__device__ int8_t updateState_cuda(STATE& currentHState, STATE& currentIState, STATE& currentDState)
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

bool getSimilarityScore(int32_t refIndex, int32_t queryIdx, std::vector<std::vector<float>>& charFreqRef, std::vector<std::vector<float>>& charFreqQuery)
{
    float score=0;
    for (int i=0; i<5; i++) score+= std::sqrt(charFreqRef[refIndex][i]*charFreqQuery[queryIdx][i]);

    return (score>0.95); //ToDo: How to clerverly decide cut-off point
} 

/*
__global__ void alignGrpToGrp_talco(char *ref, char *qry, int16_t* param, char *alignment, int32_t* seqInfo)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    // int gs = gridDim.x;
    int tidx = bx*bs+tx;

    // const size_t threadNum = 128;
    const int fLen = 256; // frontier length (assuming anti-diagonal length cannot exceed 1024)
    const int charFreqLen = 256;

    int32_t threadNum = seqInfo[6];
    int32_t blockNum = seqInfo[5];

    
    
        

    // __syncthreads();
    if (bx == 0) {
        __shared__ int8_t tile_aln[2*fLen];
        // __shared__ int16_t S  [3*fLen];
        __shared__ float S  [3*fLen];
        __shared__ int32_t CS [3*fLen];
        // __shared__ int16_t I  [2*fLen];
        __shared__ float I  [2*fLen];
        __shared__ int32_t CI [2*fLen];
        // __shared__ int16_t D  [2*fLen];
        __shared__ float D  [2*fLen];
        __shared__ int32_t CD [2*fLen];
        __shared__ int8_t tb  [64*fLen]; // May be improved to 4 bit
        __shared__ uint8_t charFreqRef [6*charFreqLen];
        __shared__ uint8_t charFreqQry [6*charFreqLen];
        __shared__ uint32_t idx [3]; 
        // [0] reference_idx
        // [1] query_idx
        // [2] globalAlnIdx
        // __shared__ int32_t max_score_list [fLen]; 
        __shared__ float max_score_list [fLen]; 
        __shared__ uint16_t max_score_ref [fLen]; 
        __shared__ uint16_t max_score_query [fLen];
        __shared__ float max_score;
        __shared__ bool last_tile;
        __shared__ bool converged; 
        __shared__ bool conv_logic;

        int32_t seqLen = seqInfo[0];
        int32_t refLen = seqInfo[1];
        int32_t qryLen = seqInfo[2];
        int32_t refNum = seqInfo[3];
        int32_t qryNum = seqInfo[4];
        
        int16_t p_match = param[0];
        int16_t p_mismatch = param[1];
        int16_t p_gapOpen = param[2];
        int16_t p_gapExtend = param[3];
        int16_t p_xdrop = param[4];
        int16_t p_marker = param[5];

        int16_t tile = 0;
        int32_t L[3], U[3];

        int32_t ftr_length [2*fLen];
        int32_t ftr_lower_limit [2*fLen];
        // initialize global values
        if (tidx == 0) {
            last_tile = false;
        } 
        if (tidx < 3) {
            idx[tidx] = 0;
        }
        
        __syncthreads();
        // if (tidx == 0) printf("refNum: %d, qryNum: %d\n", refNum, qryNum);


        while (!last_tile) {
            int32_t inf = p_xdrop + 1;
            int16_t reference_idx = idx[0];
            int16_t query_idx = idx[1];
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            // int32_t score = 0; 
            float score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            // int32_t max_score_prime = -(p_xdrop+1);
            float max_score_prime = static_cast<float>(-(p_xdrop+1));
            int32_t max_score_start_addr; 
            int32_t max_score_start_ftr;
            int32_t max_score_ref_idx;    
            int32_t max_score_query_idx;
            int16_t tb_idx = 0;
            for (size_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (size_t i = 0; i < 2*fLen; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }

            if (tx == 0) {
                max_score = 0;
                converged = false;
                conv_logic = false;
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
            // if (tidx == 0) printf("Tile: %d, refIdx: %d, qryIdx: %d, refLen: %d, qryLen: %d\n", tile, reference_idx, query_idx, reference_length, query_length);
            __syncthreads();
            
            
            // Initialize shared memory
            if (tidx < fLen) {
                for (int i = tidx; i < 3*fLen; i = i+fLen) {                     
                    S[i] = -1;
                    CS[i] = -1;
                    if (i < 2*fLen) {
                        I[i] = -1;
                        CI[i] = -1;
                        D[i] = -1;
                        CD[i] = -1;
                        tile_aln[i] = -1;
                    }
                    if (i < fLen) {
                        max_score_list [i] = -(p_xdrop+1); 
                        max_score_ref [i] = 0; 
                        max_score_query [i] = 0;
                    }
                }
            }
            // __syncthreads();
            //if (tidx == 0) printf("DEBUG: break 1\n");
            // if (tidx == 0) printf("Tile: %d, refIdx: %d, qryIdx: %d\n", tile, reference_idx, query_idx);
            // if (tidx == 0) {
            //     for (int i = 0; i < 10; ++i) printf("%c", ref[reference_idx+i]);
            //     printf("\n");
            //     for (int i = 0; i < 10; ++i) printf("%c", qry[query_idx+i]);
            //     printf("\n");
            // }
            
            __syncthreads();
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                if (tidx < fLen) {
                    max_score_list [tidx] = -(p_xdrop+1); 
                    max_score_ref [tidx] = 0; 
                    max_score_query [tidx] = 0;
                }
                __syncthreads();
                // if (tidx == 0) printf("k: %d\n", k);
                // if (tidx == 0) printf("ref+qry: %d\n", reference_length + query_length - 1);
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    // std::cout << "No more cells to compute based on x-drop critieria" << std::endl;
                    // if (tidx == 0) printf("No more cells to compute based on x-drop critieria\n");
                    if (tidx == 0) printf("No more cells to compute based on x-drop critieria (L, U) = (%d, %d)\n", L[k%3], U[k%3]+1);
                    __syncthreads();
                    break;
                }
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    // fprintf(stderr, "ERROR: anti-diagonal larger than the max limit!\n");
                    if (tidx == 0) printf("ERROR: anti-diagonal larger than the max limit!\n");
                    __syncthreads();
                    break;
                }
                // __syncthreads();
                // if (tidx == 0) printf("DEBUG: break 1.2\n");
                // __syncthreads();
                // if (tidx == 0) {
                if (k <= p_marker) {
                    // printf("ftr length: %d ftr_addr: %d  ftr lower limit: %d idx0: %d, idx1: %d\n", U[k%3] - L[k%3] + 1, ftr_addr, L[k%3], ftr_length_idx, ftr_lower_limit_idx);
                    // printf("Before idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    ftr_length[ftr_length_idx] = (U[k%3] - L[k%3] + 1);
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    // printf("before ftr_addr: %d\n", ftr_addr);
                    // printf("UL %d\n", sU[k%3] - sL[k%3] + 1);
                    ftr_addr += U[k%3] - L[k%3] + 1;
                    // printf("after ftr_addr: %d\n", ftr_addr);
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;
                    // printf("After idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    // std::cout << "ftr length: " << U[k%3] - L[k%3] + 1 << " ftr_addr: " << ftr_addr << " ftr lower limit: " << L[k%3] << " tb len: " << tb.size() << std::endl;
                    
                }
                // int32_t numRound = (U[k%3]+1-L[k%3]) / threadNum + 1;
                int32_t leftGrid = U[k%3]+1-L[k%3];
                // }   
                //__syncthreads();
                // if (tidx == 0) printf("DEBUG: break 1.4\n");
                // Calculate character frequency
                if (tidx < charFreqLen) {
                    for (int i = tidx; i < charFreqLen; i += charFreqLen) {
                        for (int j = i*6; j < i*6+6; j++) {
                            charFreqRef[j] = 0;
                            charFreqQry[j] = 0;
                        }
                        for (int j = 0; j < refNum; j++) {
                            int k = seqLen*j;
                            if      (ref[k+reference_idx+i]=='A' || ref[k+reference_idx+i]=='a') charFreqRef[i*6]+=1;
                            else if (ref[k+reference_idx+i]=='C' || ref[k+reference_idx+i]=='c') charFreqRef[i*6+1]+=1;
                            else if (ref[k+reference_idx+i]=='G' || ref[k+reference_idx+i]=='g') charFreqRef[i*6+2]+=1;
                            else if (ref[k+reference_idx+i]=='T' || ref[k+reference_idx+i]=='t') charFreqRef[i*6+3]+=1;
                            else if (ref[k+reference_idx+i]=='N' || ref[k+reference_idx+i]=='n') charFreqRef[i*6+4]+=1;
                            else charFreqRef[i*6+5]+=1;
                        }
                        for (int j = 0; j < qryNum; j++) {
                            int k = seqLen*j;
                            if      (qry[k+query_idx+i]=='A' || qry[k+query_idx+i]=='a') charFreqQry[i*6]+=1;
                            else if (qry[k+query_idx+i]=='C' || qry[k+query_idx+i]=='c') charFreqQry[i*6+1]+=1;
                            else if (qry[k+query_idx+i]=='G' || qry[k+query_idx+i]=='g') charFreqQry[i*6+2]+=1;
                            else if (qry[k+query_idx+i]=='T' || qry[k+query_idx+i]=='t') charFreqQry[i*6+3]+=1;
                            else if (qry[k+query_idx+i]=='N' || qry[k+query_idx+i]=='n') charFreqQry[i*6+4]+=1;
                            else charFreqQry[i*6+5]+=1;
                        }
                    }
                }
                // __syncthreads();
                // if (tidx == 0) printf("DEBUG1: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                // if (tidx == 0) if (tile == 462) printf("k: %d, L: %d, U: %d, (%d, %d)\n", k, L[k%3], U[k%3]+1, reference_length, query_length);
                // if (tidx == 0) printf("k: %d, Max Score: %d\n", k, max_score);
                // printf("%d, max: %d\n", tidx, max_score_list[tidx]);
                // for (int round = 0; round < numRound; ++round) {
                    if (tidx < min(leftGrid, (int32_t)threadNum)) {
                        // int16_t i=L[k%3]+round*threadNum+tidx;
                        int16_t i=L[k%3]+tidx;
                        int16_t Lprime = max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1);
                        int16_t j= min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime);
                        if (j < 0) if (tx == 0) { printf("ERROR: j less than 0.\n");}
                        // int32_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        float match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        // int score_from_prev_tile = 0;
                        float score_from_prev_tile = 0;
                        // if (tidx == 0) printf("O: %d, OD: %d OU: %d, OL: %d\n", offset, offsetDiag, offsetUp, offsetLeft);
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            if (k==0 && tile>0) {
                                score_from_prev_tile = tile*10;
                            }
                            float similarScore=0;
                            float denominator = 0;
                            float numerator = 0;
                            for (int l=0; l<6; l++) {
                                for (int m=0; m<6; m++) {
                                    denominator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m];
                                    if (m == 4 || l == 4) numerator += 0;
                                    else if (m == l)      numerator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m]*p_match;
                                    else                  numerator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m]*p_mismatch;
                                }
                                
                                // similarScore += sqrtf(charFreqRef[5*(i)+l]*charFreqQry[5*(j)+l]);
                            }
                            similarScore = numerator/denominator;
                            // if (k == 0 && tidx == 0) printf("Old: %f, New: %f\n", similarScore, numerator/denominator); 
                            // bool groupMatch = (similarScore>0.8);
                            // if (groupMatch) {
                            //     if (offsetDiag < 0) match = p_match + score_from_prev_tile; 
                            //     else match = S[(k+1)%3*fLen+offsetDiag] + p_match + score_from_prev_tile; 
                            // }
                            // else {
                            //     if (offsetDiag < 0) match = p_mismatch + score_from_prev_tile;
                            //     else match = S[(k+1)%3*fLen+offsetDiag] + p_mismatch + score_from_prev_tile;
                            // }
                            if (offsetDiag < 0) match = similarScore + score_from_prev_tile;
                            else                match = S[(k+1)%3*fLen+offsetDiag] + similarScore + score_from_prev_tile;
                            // if (tidx == 0) printf("match: %d\n", match);
                        }
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            // insOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            // insExt = I[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                            delOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            delExt = D[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                            // if (tidx == 0) printf("del: %d, %d\n", delOp, delExt);
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            // delOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            // delExt = D[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                        }
        
                        int16_t tempI, tempD, tempH;
                        tempI =  insOp;
                        tempD =  delOp;
                        Iptr = false;
                        Dptr = false;
                        if (insExt >= insOp) {
                            tempI = insExt;
                            Iptr = true;
                        }
                        if (delExt >= delOp) {
                            tempD = delExt;
                            Dptr = true;
                        }
                        if (match > tempI) {
                            if (match > tempD) {
                                tempH = match;
                                ptr = 0;
                            }
                            else {
                                tempH = tempD;
                                ptr = 2;
                            }
                        }
                        else if (tempI > tempD) {
                            tempH = tempI;
                            ptr = 1;
                        }
                        else {
                            tempH = tempD;
                            ptr = 2;
                        }
                        // printf("%d: m: %d, i: %d, d:%d, h: %d, max: %d\n",tidx, match, tempI, tempD, tempH, max_score);
                        if (tempH < max_score-p_xdrop) {
                            tempH = -inf;
                        }
                        // printf("%d: m: %d, i: %d, d:%d, h: %d\n",tidx, match, tempI, tempD, tempH);
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;
                        S[(k%3)*fLen+offset] = tempH;

                        score = tempH;

                        max_score_query[tidx] = i;
                        max_score_ref[tidx] = j;
                        max_score_list[tidx] = score;

                        // if (tile == 462) printf("No.%d (r,q)=(%d,%d) pre_max:%f, score: %f\n", tidx,  max_score_ref[tidx], max_score_query[tidx], max_score_prime, max_score_list[tidx]);
                        // if (tidx == (U[k%3]-L[k%3])/2) printf("k: %d, idx: %d, H: %d, D: %d, I: %d\n", k, tidx, tempH, tempD, tempI);
                        if (k == p_marker - 1) { // Convergence algorithm
                            CS[(k%3)*fLen+offset] = (3 << 16) | (i & 0xFFFF); 
                            // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  tempCS << "\n";
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 16) | (i & 0xFFFF);
                            CD[(k%2)*fLen+offset] = (2 << 16) | (i & 0xFFFF);
                            // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  tempCS <<  " " << tempCI <<  " " << tempCD << "\n";
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
                        // if (tile == 0 && tidx == 0) printf("k: %d, CS: %d, CI: %d, CD: %d, Iptr: %d, Dptr: %d\n",k, CS[k%3*fLen+offset], CI[k%2*fLen+offset], CD[k%2*fLen+offset], Iptr, Dptr);
                        if (Iptr) {
                            // std::cout << (ptr & 0xFF) << " ";
                            ptr |= 0x04; 
                            // std::cout << (ptr & 0xFF) << "\n";
                        }
                        if (Dptr) {
                            // std::cout << (ptr & 0xFF) << " ";
                            ptr |= 0x08;
                            // std::cout << (ptr & 0xFF) << "\n";
                        }
                        // if (tidx == 0) printf("Before KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);                      
                        if (k <= p_marker){
                            tb[tb_idx+tidx] = ptr;
                        }
                        // if (tidx == 0) printf("After KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    }
                    __syncthreads();
                    // Calculate Max
                    for (uint32_t r = threadNum/2; r > 0; r >>= 1) {
                        if (tidx < r) {
                            // if (tile == 462) printf("Reduction max No.%d & %d, max_score: %f & %f, pre_max: %f\n", tidx, tidx+r, max_score_list[tidx], max_score_list[tidx+r], max_score_prime);
                            max_score_ref[tidx]   = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_ref[tidx+r] : max_score_ref[tidx];
                            max_score_query[tidx] = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_query[tidx+r] : max_score_query[tidx];
                            max_score_list[tidx]  = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_list[tidx+r] : max_score_list[tidx];
                           
                        }
                        __syncthreads();
                    }
                    // if (tidx == 0) if (tile == 462) printf("max (%d, %d) = %f\n",max_score_ref[0], max_score_query[0], max_score_list[0]);
                    // Update Max
                    if (max_score_prime < max_score_list[0]) {
                        max_score_prime = max_score_list[0];
                        // if (tidx == 0) {
                        // max_score_prime = max_score_list[0];
                        if (k <= p_marker) {
                            max_score_ref_idx = max_score_ref[0];
                            max_score_query_idx = max_score_query[0];
                            max_score_start_addr = ftr_addr - (U[k%3] - L[k%3] + 1)  + (max_score_query[0] - L[k%3]);
                            max_score_start_ftr = k;
                        }
                            // printf("max (%d, %d) = %d, %d, %d\n",max_score_ref_idx, max_score_query_idx, max_score_prime, max_score_start_addr, max_score_list[0]);
                        // }
                    }
                    // __syncthreads();
                    
                // }
                // __syncthreads();
                if (k <= p_marker){
                    tb_idx += min((int)leftGrid, (int)threadNum);
                }
                // if (tidx == 0) printf("DEBUG2: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                // if (tidx == 0) printf("k: %d,newL: %d, newU: %d\n",k, L[k%3], U[k%3]);
                while (newL <= U[k%3]) {         
                    int32_t offset = newL - L[k%3];
                    // printf("k: %d,newL: %d, offset: %d\n",S[(k%3)*fLen+offset], L[k%3], offset);
                    if (S[(k%3)*fLen+offset] <= -inf) {
                        newL++;
                    }
                    else {
                        break;
                    }
                }
                while (newU >= L[k%3]) {
                    int32_t offset = newU - L[k%3];
                    // if (tidx == 0) printf("DDDDD: %d, %d, %d, %d\n", newU, L[k%3], offset, S[(k%3)*fLen+offset]);
                    if (S[(k%3)*fLen+offset] <= -inf) {
                        newU--;
                    }
                    else {
                        break;
                    }
                }
                int32_t v1 = static_cast<int32_t>(query_length)-1;
                int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                int32_t v3 = newU+1;
                int32_t Lprime = max(static_cast<int32_t>(0), v2);
                // __syncthreads();
                // if (tidx == 0) printf("newU: %d, + %d = %d ", newU, 1, v3);
                // if (tidx == 0) printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                // if (tidx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
        
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 
        
                // if (tidx == 0) printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
        
                
                // sL[(k+1)%3] = L[(k+1)%3];
                // sU[(k+1)%3] = U[(k+1)%3];

                if (tidx == 0) {
                    // printf("k: %d,newL: %d, newU: %d\n",k, newL, newU);
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
                        // printf("k: %d, CONV: %d, %d, %d\n", k, conv_I, conv_D, conv_S);
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                            conv_score = max_score_prime;
                            // if (DEBUG)  std::cout << "Converged at: " << conv_value << "\n";
                            // printf("Converged at: %d\n", conv_value);
                        }
                        prev_conv_s = conv_S;
                    }
                    // int32_t v1 = static_cast<int32_t>(query_length)-1;
                    // int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                    // int32_t v3 = newU+1;

                    // int32_t Lprime = max(static_cast<int32_t>(0), v2);

                    // printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                    // // if (tidx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
            
                    // L[(k+1)%3] = max(newL, Lprime);
                    // U[(k+1)%3] = min(v1, v3); 
            
                    // printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                    // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
            
                    
                    // sL[(k+1)%3] = L[(k+1)%3];
                    // sU[(k+1)%3] = U[(k+1)%3];
    

                    // printf("DEBUG44: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                    // printf("DEBUG6 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
                    // Update max_score
                    max_score = max_score_prime;
                    // printf("Max score: %d\n", max_score);

                    if ((converged) && (max_score > conv_score)){
                        conv_logic = true;
                        // printf("Tile %d Convergence logic found: %d\n", tile, max_score);
                        // if (DEBUG) std::cout << "Convergence logic found: ";
                        // break;
                    }
                    
                }
                __syncthreads();
                if (conv_logic) break;
            }
            __syncthreads();
            // if (tidx == 0) printf("DEBUG: break 2\n");
            if (tidx == 0) {
                // printf("Frontier addr: %d \ntb_start_ftr: %d\nmarker: %d\n", ftr_addr, ftr_length_idx, p_marker);
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                // printf("ftr_addr: %d\n", ftr_addr);
                // printf("ftr_l_idx: %d\n", ftr_length_idx);
                // printf("ftr[]: %d\n", ftr_length[ftr_length_idx - 1]);
                // printf("limit_idx: %d\n", ftr_lower_limit_idx);
                // printf("limit[]: %d\n", ftr_lower_limit[ftr_lower_limit_idx - 1]);
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFFFF;
                    tb_state = (conv_value >> 16) & 0xFFFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                                                      tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    // if (DEBUG) std::cout <<  " conv query idx: " << conv_query_idx << " " << (tb_state&0xFFFF) << " " << conv_ref_idx << " " << conv_value << std::endl;
                    // printf(" conv query idx: %d, %d, %d, %d\n", conv_query_idx ,(tb_state&0xFFFF) ,conv_ref_idx , conv_value );
                
                } else {
                    conv_query_idx = max_score_query_idx;
                    conv_ref_idx = max_score_ref_idx;
                    tb_start_addr = max_score_start_addr;
                    tb_start_ftr = max_score_start_ftr;
                    tb_state = 0;
                    // printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", max_score_query_idx, max_score_ref_idx, max_score_start_addr, max_score_start_ftr);
                    last_tile = true;
                }
                // printf("Before: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                // printf("After: convref: %d, convqry:%d\n", conv_ref_idx, conv_query_idx);
                // printf("After: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                // Start Traceback
                // printf("DEBUG: break 3\n");
                // Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln);
                int32_t addr = tb_start_addr; 
                int16_t ftr = tb_start_ftr;
                int16_t traceback_idx = conv_query_idx;
                int16_t qry_idx = conv_query_idx;
                int16_t ref_idx = conv_ref_idx;
                int8_t  tb_value = 0;
                state = tb_state%3;

                int8_t  dir = 0;
                // bool checkpoint = false;
                // int count = 0;
                // if(tile == 155) printf("ftr: %d\n", ftr);
                // if(tile == 155) printf("state: %d\n", tb_state);
                // if(tile == 155) printf("addr: %d\n", addr);
                // if(tile == 155) printf("aln_idx: %d\n", aln_idx);
                // if(tile == 155) printf("dir: %d\n", dir);
                // if(tile == 155) printf("traceback_idx: %d\n", traceback_idx);
                while (ftr >= 0) {
                    // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
                    if (addr < 0) {
                        printf("ERROR: tb addr < 0 (%d)!\n", addr);
                        // exit(1);
                        break;
                    }
                    // count++;
                    tb_value = tb[addr];
                    // if (DEBUG)
                    // {
                    // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
                    // printf(" fL[ftr - 1]: %d, ll[ftr-2]: %d\n", ftr_length[ftr - 1], ftr_lower_limit[ftr-2]);
                    // printf(" fL[ftr - 2]: %d", ftr_length[ftr - 2]);
                    // printf(" Tb: %d", ( tb_value&0xFFFF));
                    // }
                    // if (tile < 10) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d, tb_value: %d\n", addr, (state & 0xFFFF), ftr, traceback_idx , ftr_lower_limit[ftr], tb_value);
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
                    // if (tile == 0) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d, tb_value: %d, dir: %d, aln_idx: %d\n", addr, (state & 0xFFFF), ftr, traceback_idx , ftr_lower_limit[ftr], tb_value, dir, aln_idx);
                    addr = addr - (traceback_idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);

                    if (dir == 0) {
                        addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]);
                        ftr -= 2;
                        traceback_idx -= 1;
                        qry_idx--;
                        ref_idx--;
                    }
                    else if (dir == 1) {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]);
                        ftr -= 1;
                        traceback_idx -=1;
                        qry_idx--;
                    }
                    else {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1);
                        ftr -= 1;
                        ref_idx--;
                    }
                    tile_aln[aln_idx] = dir;
                    ++aln_idx;    
                    // state = next_state;
                    // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
                }
                // std::cout << count << std::endl;   

                // if (DEBUG) std::cout << ref_idx << " " << qry_idx << std::endl; 
                
                state = tb_state % 3;
                // if (DEBUG) {
                //     std::cout << "tb_state: " <<  (tb_state & 0xFFFF) << std::endl;
                //     int count = 0;
                //     for (auto &a: aln){
                //         std::cout << count << ": " << (a & 0xFFFF) << "\n";
                //         count += 1;
                //     }
                // }

                // Write global memory
                int32_t refIndex = reference_idx;
                int32_t qryIndex = query_idx;
                // reference_idx += conv_query_idx;
                // query_idx += conv_ref_idx;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
                int16_t globalAlnIdx = idx[2];
                // printf("TB: tile: %d, refidx: %d, qryidx:%d, globalAln: %d, aln_idx: %d\n", tile, reference_idx, query_idx, globalAlnIdx, aln_idx);
                for (int j = aln_idx -1; j >= 0; --j) {
                    // if (tile >= 0 && tile <= 2 ) printf("j: %d, dir:%d, globalAln: %d\n", j, tile_aln[j], globalAlnIdx);
                    // if (tile >= 0 && tile <= 2 ) printf("ridx: %d, qidx:%d, r: %c, q: %c\n", refIndex, qryIndex, ref[refIndex], qry[qryIndex]);
                    if (j == aln_idx-1 && tile>0){
                        // if (tile == 1) {
                            ++refIndex;
                            ++qryIndex;
                        // }
                        continue;
                    }
                    if ( (tile_aln[j] & 0xFFFF) == 0) {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = ref[i*seqLen+refIndex]; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = qry[i*seqLen+qryIndex]; 
                        qryIndex++;refIndex++;
                    }
                    else if ((tile_aln[j] & 0xFFFF) == 1) {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = ref[i*seqLen+refIndex];  
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = '-'; 
                        refIndex++;
                    }
                    else {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = '-'; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = qry[i*seqLen+qryIndex]; 
                        qryIndex++;
                    }
                    globalAlnIdx++;
                }
                // printf("Finish TB of tile %d!!!!\n", tile);
                idx[0] = reference_idx;
                idx[1] = query_idx;
                idx[2] = globalAlnIdx;
                // tile++;
            }
           
            __syncthreads();
            if (tidx < threadNum) tile++;
        }
        __syncthreads();
        
        // __syncthreads();
        // TODO: Add score to debug
        // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
    }


}
*/

__global__ void alignGrpToGrp_talco(/*char *seqs, */ uint16_t* freq, int8_t *aln, /*int32_t* seqIdx,*/ int32_t* len, int32_t* alnLen, int32_t *seqInfo, int16_t* param)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    // int gs = gridDim.x;
    // int tidx = bx*bs+tx;

    // const size_t threadNum = 128;
    const int fLen = 512; // frontier length (assuming anti-diagonal length cannot exceed 1024)
    const int16_t reducedValue = (1 << 14);
    
    int32_t threadNum = seqInfo[3];
    int32_t blockNum = seqInfo[4];
    int32_t pairNum = seqInfo[2];    

    __syncthreads();
    if (bx < pairNum) {
    // if (bx == 0) {
        __shared__ int8_t tile_aln[2*fLen];
        __shared__ int16_t S  [3*fLen];
        __shared__ int32_t CS [3*fLen];
        __shared__ int16_t I  [2*fLen];
        __shared__ int32_t CI [2*fLen];
        __shared__ int16_t D  [2*fLen];
        __shared__ int32_t CD [2*fLen];
        __shared__ int8_t tb  [32*fLen]; // May be improved to 4 bit
        __shared__ uint32_t idx [3]; 
        // [0] reference_idx
        // [1] query_idx
        // [2] globalAlnIdx
        __shared__ int16_t max_score_list [fLen]; 
        __shared__ int16_t max_score;
        __shared__ uint16_t reducedTimes;
        __shared__ bool last_tile;
        __shared__ bool converged; 
        __shared__ bool conv_logic;
        
        int32_t seqLen = seqInfo[0];
        int32_t seqNum = seqInfo[1];
        int32_t refLen = len[2*bx];
        int32_t qryLen = len[2*bx+1];
        // int32_t refStartIdx = seqIdx[2*bx];
        // int32_t qryStartIdx = seqIdx[2*bx+1];
        // int32_t refNum = seqIdx[2*bx+1] - seqIdx[2*bx];
        // int32_t qryNum = 0; 
        // if (bx == pairNum-1) qryNum = seqNum - seqIdx[2*bx+1];
        // else                 qryNum = seqIdx[2*bx+2] - seqIdx[2*bx+1]; 

        // if (tx == 0) printf("refNum: %d, qryNum: %d\n", refNum, qryNum);
        // if (tx == 0 && bx == 4) printf("refLen: %d, qryLen: %d\n", refLen, qryLen);
       
        int16_t p_match = param[0];
        int16_t p_mismatch = param[1];
        int16_t p_gapOpen = param[2];
        int16_t p_gapExtend = param[3];
        int16_t p_xdrop = param[4];
        int16_t p_marker = param[5];

        int16_t tile = 0;
        int32_t last_k = 0;
        int32_t L[3], U[3];

        int32_t ftr_length [2*fLen];
        int32_t ftr_lower_limit [2*fLen];
        // initialize global values
        if (tx == 0) {
            last_tile = false;
            // no_converge = false;
        } 
        if (tx < 3) {
            idx[tx] = 0;
        }

        __syncthreads();
        
        while (!last_tile) {
            int16_t inf = p_xdrop + 1;
            int16_t reference_idx = idx[0];
            int16_t query_idx = idx[1];
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            // int32_t score = 0; 
            int16_t score = 0;
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            // int32_t max_score_prime = -(p_xdrop+1);
            int16_t max_score_prime = -(p_xdrop+1);
            int32_t max_score_start_addr; 
            int32_t max_score_start_ftr;
            // int32_t max_score_ref_idx;    
            // int32_t max_score_query_idx;
            int16_t tb_idx = 0;
            for (size_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (size_t i = 0; i < 2*fLen; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }
            if (tx == 0) {
                // max_score = 0;
                // if (max_score > reducedValue) {
                //     max_score -= reducedValue;
                //     reducedTimes += 1;
                // }
                max_score = 0;
                converged = false;
                conv_logic = false;
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
                    if (i < fLen) {
                        max_score_list [i] = -(p_xdrop+1); 
                        // max_score_ref [i] = 0; 
                        // max_score_query [i] = 0;
                    }
                }
            }
            __syncthreads();

            
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                // Initial shared memory for recording max scores and index
                if (tx < fLen) {
                    max_score_list[tx] = -(p_xdrop+1); 
                    // max_score_ref[tx] = 0; 
                    // max_score_query[tx] = 0;
                }
                __syncthreads();
                // if (tx == 0) printf("k: %d\n", k);
                // if (tidx == 0) printf("ref+qry: %d\n", reference_length + query_length - 1);
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    if (tx == 0) printf("No.%d No more cells to compute based on x-drop critieria (L, U) = (%d, %d)\n", bx, L[k%3], U[k%3]+1);
                    // if (tx == 0) last_tile = true;
                    __syncthreads();
                    break;
                }
                if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                    if (tx == 0) printf("No.%d ERROR: anti-diagonal larger than the max limit!\n", bx);
                    if (tx == 0) last_tile = true;
                    __syncthreads();
                    break;
                }
                
                // __syncthreads();
                if (k <= p_marker) {
                    ftr_length[ftr_length_idx] = (U[k%3] - L[k%3] + 1);
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    ftr_addr += U[k%3] - L[k%3] + 1;
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;      
                }
                __syncthreads();
                // Calculate character frequency (original method)
                // if (tx < fLen) {
                //     for (int j = tx*6; j < tx*6+6; j++) {
                //         charFreqRef[j] = 0;
                //         charFreqQry[j] = 0;
                //     }
                //     for (int j = 0; j < refNum; j++) {
                //         int k = seqLen*(refStartIdx+j);
                //         if      (seqs[k+reference_idx+tx]=='A' || seqs[k+reference_idx+tx]=='a') charFreqRef[tx*6]+=1;
                //         else if (seqs[k+reference_idx+tx]=='C' || seqs[k+reference_idx+tx]=='c') charFreqRef[tx*6+1]+=1;
                //         else if (seqs[k+reference_idx+tx]=='G' || seqs[k+reference_idx+tx]=='g') charFreqRef[tx*6+2]+=1;
                //         else if (seqs[k+reference_idx+tx]=='T' || seqs[k+reference_idx+tx]=='t') charFreqRef[tx*6+3]+=1;
                //         else if (seqs[k+reference_idx+tx]=='N' || seqs[k+reference_idx+tx]=='n') charFreqRef[tx*6+4]+=1;
                //         else charFreqRef[tx*6+5]+=1;
                //     }
                //     for (int j = 0; j < qryNum; j++) {
                //         int k = seqLen*(qryStartIdx+j);
                //         if      (seqs[k+query_idx+tx]=='A' || seqs[k+query_idx+tx]=='a') charFreqQry[tx*6]+=1;
                //         else if (seqs[k+query_idx+tx]=='C' || seqs[k+query_idx+tx]=='c') charFreqQry[tx*6+1]+=1;
                //         else if (seqs[k+query_idx+tx]=='G' || seqs[k+query_idx+tx]=='g') charFreqQry[tx*6+2]+=1;
                //         else if (seqs[k+query_idx+tx]=='T' || seqs[k+query_idx+tx]=='t') charFreqQry[tx*6+3]+=1;
                //         else if (seqs[k+query_idx+tx]=='N' || seqs[k+query_idx+tx]=='n') charFreqQry[tx*6+4]+=1;
                //         else charFreqQry[tx*6+5]+=1;
                //     }
                // }
                __syncthreads();
                
                // if (tx == 0 && bx == 1 && tile >= 462) printf("Tile: %d, k: %d, L: %d, U: %d, (%d, %d)\n",tile, k, L[k%3], U[k%3]+1, reference_length, query_length);
                int32_t leftGrid = U[k%3]+1-L[k%3];
                if (tx < leftGrid) {
                    int16_t i=L[k%3]+tx;
                    int16_t Lprime = max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1);
                    int16_t j= min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime);
                    if (j < 0) if (tx == 0) { printf("ERROR: j less than 0.\n");}
                    int16_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                    int32_t offset = i-L[k%3];
                    int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                    int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                    int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                    int16_t score_from_prev_tile = 0;
                    if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                        if (k==0 && tile>0) {
                            score_from_prev_tile = tile*10;
                            // score_from_prev_tile = max_score;
                        }
                        int16_t similarScore = 0;
                        float denominator = 0;
                        float numerator = 0;
                        // for (int l=0; l<6; l++) {
                        //     for (int m=0; m<6; m++) {
                        //         denominator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m];
                        //         if (m == 4 || l == 4) numerator += 0;
                        //         else if (m == l)      numerator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m]*p_match;
                        //         else                  numerator += charFreqRef[6*(j)+l]*charFreqQry[6*(i)+m]*p_mismatch;
                        //     }
                        // }
                        int refFreqStart = 6*(2*bx)*seqLen;
                        int qryFreqStart = 6*(2*bx+1)*seqLen;
                        int refFreqIdx = reference_idx + j;
                        int qryFreqIdx = query_idx + i;
                        for (int l=0; l<6; l++) {
                            for (int m=0; m<6; m++) {
                                denominator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m];
                                if (m == 4 || l == 4) numerator += 0;
                                else if (m == l)      numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_match;
                                else                  numerator += freq[refFreqStart+6*(refFreqIdx)+l]*freq[qryFreqStart+6*(qryFreqIdx)+m]*p_mismatch;
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
                    // max_score_query[tx] = i;
                    // max_score_ref[tx] = j;
                    max_score_list[tx] = score;
                    // if (bx == 4 && tx == (U[k%3]-L[k%3])/2 && (tile < 10)) printf("k: %d, (%d, %d), idx: %d, H: %d, D: %d, I: %d\n", k, i, j, tx, tempH, tempD, tempI);
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
                    // if(bx == 4 && i == L[k%3] && tile == 100) printf("k: %d, Convergence Unique Id's: %d, %d, %d\n", k, CS[(k%3)*fLen+offset], CI[(k%2)*fLen+offset], CD[(k%2)*fLen+offset]);
                    if (Iptr) {
                        // std::cout << (ptr & 0xFF) << " ";
                        ptr |= 0x04; 
                        // std::cout << (ptr & 0xFF) << "\n";
                    }
                    if (Dptr) {
                        // std::cout << (ptr & 0xFF) << " ";
                        ptr |= 0x08;
                        // std::cout << (ptr & 0xFF) << "\n";
                    }
                    // if (tidx == 0) printf("Before KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);                      
                    if (k <= p_marker){
                        tb[tb_idx+tx] = ptr;
                    }
                    // if (tidx == 0) printf("After KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                }
                __syncthreads();
                // Calculate Max
                for (uint32_t r = fLen/2; r > 0; r >>= 1) {
                    if (tx < r) {
                        // if (tile == 462) printf("Reduction max No.%d & %d, max_score: %f & %f, pre_max: %f\n", tidx, tidx+r, max_score_list[tidx], max_score_list[tidx+r], max_score_prime);
                        // max_score_ref[tx]   = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_ref[tx+r] : max_score_ref[tx];
                        // max_score_query[tx] = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_query[tx+r] : max_score_query[tx];
                        max_score_list[tx]  = (max_score_list[tx+r] > max_score_list[tx]) ? max_score_list[tx+r] : max_score_list[tx];
                    }
                    __syncthreads();
                }
                // if (tidx == 0) if (tile == 462) printf("max (%d, %d) = %f\n",max_score_ref[0], max_score_query[0], max_score_list[0]);
                // Update Max
                if (max_score_prime < max_score_list[0]) {
                    max_score_prime = max_score_list[0];
                    if (k <= p_marker) {
                        // max_score_ref_idx = max_score_ref[0];
                        // max_score_query_idx = max_score_query[0];
                        // max_score_start_addr = ftr_addr - (U[k%3] - L[k%3] + 1)  + (max_score_query[0] - L[k%3]);
                        max_score_start_ftr = k;
                    }
                }
                __syncthreads();


                if (k <= p_marker){
                    tb_idx += (int)leftGrid;
                }
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
                            // if (bx == 4) printf("Converged at: %d\n", conv_value);
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
                last_k = k;
                if (conv_logic) break;
            }
            __syncthreads();
            // if (tx == 0) printf("DEBUG: break 2\n");
            if (tx == 0) {
                // printf("Frontier addr: %d \ntb_start_ftr: %d\nmarker: %d\n", ftr_addr, ftr_length_idx, p_marker);
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                // printf("ftr_addr: %d\n", ftr_addr);
                // printf("ftr_l_idx: %d\n", ftr_length_idx);
                // printf("ftr[]: %d\n", ftr_length[ftr_length_idx - 1]);
                // printf("limit_idx: %d\n", ftr_lower_limit_idx);
                // printf("limit[]: %d\n", ftr_lower_limit[ftr_lower_limit_idx - 1]);
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFFFF;
                    tb_state = (conv_value >> 16) & 0xFFFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                                                      tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    // if (DEBUG) std::cout <<  " conv query idx: " << conv_query_idx << " " << (tb_state&0xFFFF) << " " << conv_ref_idx << " " << conv_value << std::endl;
                    // if (bx == 4 && tile == 100) printf(" conv query idx: %d, %d, %d, %d\n", conv_query_idx ,(tb_state&0xFFFF) ,conv_ref_idx , conv_value );
                
                } 
                else {
                    // int32_t last_k = reference_length + query_length - 2;
                    // printf("k: %d, Convergence Unique Id's: %d, %d, %d\n", last_k, CS[(last_k%3)*fLen], CI[(last_k%3)*fLen], CD[(last_k%3)*fLen]);
                    if (last_tile == true) {
                        conv_query_idx = 0;
                        conv_ref_idx = 0;
                        tb_start_addr = 0;
                        tb_start_ftr = -1;
                        tb_state = 0; 
                        alnLen[bx] = 0;
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
                        conv_query_idx = CS[(last_k%3)*fLen] & 0xFFFF;
                        tb_state = (CS[(last_k%3)*fLen] >> 16) & 0xFFFF;
                        conv_ref_idx = p_marker - conv_query_idx; 
                        conv_ref_idx -= (tb_state == 3) ? 1: 0;
                        tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                        tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                                                          tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                        tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                        // conv_query_idx = max_score_query_idx;
                        // conv_ref_idx = max_score_ref_idx;
                        // tb_start_addr = max_score_start_addr;
                        // tb_start_ftr = max_score_start_ftr;
                    }
                    
                    // if (bx == 1 && tx == 0 && tile >= 462) printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", max_score_query_idx, max_score_ref_idx, max_score_start_addr, max_score_start_ftr);
                    // if (bx == 1 && tx == 0 && tile >= 462) printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", conv_query_idx, conv_ref_idx, tb_start_addr , tb_start_ftr);
                    // if (bx == 1 && tx == 0 && tile >= 462) printf("MAX: qryLen: %d, refLen: %d\n", query_length, reference_length);
                    
                }
                // printf("Before: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                // printf("After: convref: %d, convqry:%d\n", conv_ref_idx, conv_query_idx);
                // printf("After: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                // Start Traceback
                // printf("DEBUG: break 3\n");
                // Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln);
                int32_t addr = tb_start_addr; 
                int16_t ftr = tb_start_ftr;
                int16_t traceback_idx = conv_query_idx;
                int16_t qry_idx = conv_query_idx;
                int16_t ref_idx = conv_ref_idx;
                int8_t  tb_value = 0;
                state = tb_state%3;
                int8_t  dir = 0;
                // bool checkpoint = false;
                // int count = 0;
                // if(tile == 155) printf("ftr: %d\n", ftr);
                // if(tile == 155) printf("state: %d\n", tb_state);
                // if(tile == 155) printf("addr: %d\n", addr);
                // if(tile == 155) printf("aln_idx: %d\n", aln_idx);
                // if(tile == 155) printf("dir: %d\n", dir);
                // if(tile == 155) printf("traceback_idx: %d\n", traceback_idx);
                while (ftr >= 0) {
                    // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
                    if (addr < 0) {
                        printf("ERROR: tb addr < 0 (%d)!\n", addr);
                        // exit(1);
                        break;
                    }
                    // count++;
                    tb_value = tb[addr];
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
                    // if (tile == 0) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d, tb_value: %d, dir: %d, aln_idx: %d\n", addr, (state & 0xFFFF), ftr, traceback_idx , ftr_lower_limit[ftr], tb_value, dir, aln_idx);
                    addr = addr - (traceback_idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);

                    if (dir == 0) {
                        addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]);
                        ftr -= 2;
                        traceback_idx -= 1;
                        qry_idx--;
                        ref_idx--;
                    }
                    else if (dir == 1) {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]);
                        ftr -= 1;
                        traceback_idx -=1;
                        qry_idx--;
                    }
                    else {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1);
                        ftr -= 1;
                        ref_idx--;
                    }
                    tile_aln[aln_idx] = dir;
                    ++aln_idx;    
                    // state = next_state;
                    // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
                }
                // std::cout << count << std::endl;   

                // if (DEBUG) std::cout << ref_idx << " " << qry_idx << std::endl; 
                
                state = tb_state % 3;

                // Write global memory
                int32_t refIndex = reference_idx;
                int32_t qryIndex = query_idx;
                // reference_idx += conv_query_idx;
                // query_idx += conv_ref_idx;
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
                // if (bx == 4 && tx == 0 ) printf("TB: tile: %d, refidx: %d, qryidx:%d, last:%d\n", tile, reference_idx, query_idx, last_tile);
                // printf("TB: tile: %d, refidx: %d, qryidx:%d, globalAln: %d, aln_idx: %d\n", tile, reference_idx, query_idx, globalAlnIdx, aln_idx);
            }
            __syncthreads();
            int32_t numRound = aln_idx / threadNum + 1;
            int32_t lastRound = aln_idx % threadNum;
            // int16_t globalAlnIdx = idx[2];
            int32_t alnStartIdx = bx * 2 * seqLen + alnLen[bx];
            for (int round = 0; round < numRound; ++round) {
                int txNum = (round == numRound - 1) ? lastRound : threadNum;
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
                            aln[alnStartIdx+threadNum*round+x] = tile_aln[aln_idx-1-threadNum*round-x];
                        }
                    }
                    else {
                        for (int x = 0; x < txNum; ++x) {
                            aln[alnStartIdx+threadNum*round+x] = tile_aln[aln_idx-2-threadNum*round-x];
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
                // idx[2] = globalAlnIdx;
            }
            __syncthreads();
            tile++;
        }
        __syncthreads();
        // __syncthreads();
        // TODO: Add score to debug
        // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
    }


}