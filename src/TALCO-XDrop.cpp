#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

void Talco_xdrop::Align_freq (
    Params params,
    const std::vector<std::vector<float>>& freqRef,
    const std::vector<std::vector<float>>& freqQry,
    const std::vector<std::vector<float>>& gapOp,
    const std::vector<std::vector<float>>& gapEx,
    const std::vector<std::vector<float>>& gapCl,
    const std::pair<float, float>& num,
    std::vector<int8_t>& aln,
    int16_t& errorType
) {
        
        int8_t state = 0;
        int32_t reference_idx, query_idx;
        reference_idx = 0; query_idx = 0;
        bool last_tile = false;
        int tile = 0;
        while (!last_tile) {
            std::vector<int8_t> tile_aln;
            Talco_xdrop::Tile(
                freqRef, 
                freqQry, 
                gapOp,
                gapEx,
                gapCl,
                num,
                params, 
                reference_idx, 
                query_idx, 
                tile_aln, 
                state, 
                last_tile, 
                tile, 
                errorType
            );
            if (tile_aln.empty()) {
                aln.clear();
                return;
            }
            for (int i= tile_aln.size()-1; i>=0; i--){
                if (i == tile_aln.size()-1 && tile>0) continue;
                aln.push_back(tile_aln[i]);
            }
            tile_aln.clear();
            tile++;
        }
        return;
}

int32_t Talco_xdrop::Reduction_tree(const int32_t *C, const int32_t start, const int32_t length){
    int32_t conv = C[start];
    for (int32_t i = start + 1; i <= start + length; i++ ){
        if (conv != C[i]){
            conv = -1;
            break;
        } 
    }
    return conv;
}

void Talco_xdrop::Traceback(
    const std::vector<int32_t> &ftr_length, 
    const std::vector<int32_t> &ftr_lower_limit, 
    const int32_t tb_start_addr, 
    const int16_t tb_start_ftr,
    const int8_t tb_state,
    const int16_t tb_start_idx,
    const int16_t ref_start_idx,
    const std::vector<int8_t> &tb,
    std::vector<int8_t> &aln,
    bool firstTile
    ){
    int32_t addr = tb_start_addr; 
    int16_t ftr = tb_start_ftr;
    int16_t idx = tb_start_idx;
    int16_t query_idx = tb_start_idx;
    int16_t ref_idx = ref_start_idx;
    int8_t  state = tb_state;
    int8_t  tb_value = 0;

    int8_t  dir = 0;
    while (ftr >= 0) {
        if (addr < 0) {
            fprintf(stderr, "ERROR: tb addr < 0!\n");
            // exit(1);
        }
        tb_value = tb[addr];
        if (state == 0) { // Current State M
            state = tb_value & 0x03;
            if (state == 0)
            {
                dir = 0;
            }
            else if (state == 1)
            {
                dir = 1;
                if (tb_value & 0x04) {
                    state = 1;
                } else {
                    state = 0;
                }   
            }
            else
            {
                dir = 2;
                if (tb_value & 0x08) {
                    state = 2;
                } else {
                    state = 0;
            }
            }

        } else if (state == 1) { // Current State I
            dir = 1;
            if (tb_value & 0x04) {
                state = 1;
            } else {
                state = 0;
            }
        } else { // Current State D
            dir = 2;
            if (tb_value & 0x08) {
                state = 2;
            } else {
                state = 0;
            }
        }
        addr = addr - (idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);
        if (dir == 0){
            addr = addr - (ftr_length[ftr - 2]) + (idx - ftr_lower_limit[ftr  - 2]);
            ftr -= 2;
            idx -= 1;
            query_idx--;
            ref_idx--;
        }else if (dir == 1){
            addr = addr + (idx - ftr_lower_limit[ftr  - 1]);
            ftr -= 1;
            idx -=1;
            query_idx--;
        }else{
            addr = addr + (idx - ftr_lower_limit[ftr  - 1] + 1);
            ftr -= 1;
            ref_idx--;
        }
        aln.push_back(dir);   
    }
}

void Talco_xdrop::Tile (
    // const std::vector<std::string>& reference,
    // const std::vector<std::string>& query,
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp,
    const std::vector<std::vector<float>>& gapEx,
    const std::vector<std::vector<float>>& gapCl,   
    const std::pair<float, float>& num,
    Params param,
    int32_t &reference_idx,
    int32_t &query_idx,
    std::vector<int8_t> &aln,
    int8_t &state,
    bool &last_tile,
    const int &tile,
    int16_t& errorType // 1: xdrop, 2: exceed anti-doagonal limit
    ) {
        
        // Initialising variables
        int32_t inf = param.xdrop + 1;
        int32_t marker = param.marker;
        int32_t fLen = param.fLen;
        bool converged = false; bool conv_logic = false;
        float refNum = num.first, qryNum = num.second; 
        int32_t reference_length = reference.size() - reference_idx; 
        int32_t query_length = query.size() - query_idx;
        int32_t score = 0; int32_t max_score = 0; int32_t max_score_prime = -inf; // int32_t max_score_ref_idx = 0; int32_t max_score_query_idx = 0;
        int32_t conv_score = 0; int32_t conv_value = 0; int32_t conv_ref_idx = 0; int32_t conv_query_idx = 0; 
        int32_t tb_start_addr = 0; int32_t tb_start_ftr = 0; // int32_t max_score_start_addr = 0; int32_t max_score_start_ftr = 0;
        int8_t tb_state = 0;
        int8_t ptr = 0;  // Main pointer
        bool Iptr = false;
        bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend

        // For xdrop
        // int32_t max_score_marker = -inf; int32_t max_score_marker_ref_idx = 0; int32_t max_score_marker_query_idx = 0;
        // int32_t max_score_marker_start_addr = 0; int32_t max_score_marker_start_ftr = 0;
        // int8_t tb_state_marker = 0;

        float denominator = refNum * qryNum;
                    

        int32_t L[3], U[3];
        int32_t *S[3], *I[2], *D[2];
        int32_t *CS[3], *CI[2], *CD[2];
        std::vector<int8_t> tb;
        std::vector<int32_t> ftr_length;
        std::vector<int32_t> ftr_lower_limit;
        int32_t ftr_addr = 0;
        int32_t last_k = 0;
        // int32_t xdrop = false;
        int32_t prev_conv_s = -1;
        paramType scoreMat [25];
        // paramType gapOpen = param.gapOpen;
        paramType gapExtend = param.gapExtend;
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) scoreMat[i*5+j] = param.scoreMatrix[i][j];
        for (size_t sIndx=0; sIndx<3; sIndx++) { // Allocate memory for S, I, D, and CS array
            S[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
            CS[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
            if (sIndx < 2) {
                I[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
                D[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
                CI[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
                CD[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
            }
            L[sIndx] = sIndx;
            U[sIndx] = -sIndx;
        }
        
        for (int16_t sIndx=0; sIndx<3; sIndx++){ // Initialise memory for S, I, D array
            for (int16_t sLenIndex=0; sLenIndex<fLen; sLenIndex++){
                S[sIndx][sLenIndex] = -1;
                CS[sIndx][sLenIndex] = -1;
                if (sIndx < 2) {  
                    I[sIndx][sLenIndex] = -1;
                    D[sIndx][sLenIndex] = -1;
                    CI[sIndx][sLenIndex] = -1;
                    CD[sIndx][sLenIndex] = -1;
                }
            }
        }

        if ((reference_length < 0) || (query_length < 0)) {
            std::cout << reference_length << " " << query_length << std::endl;
            fprintf(stderr, "ERROR: Reference/Query index exceeded limit!\n");
            errorType = 3;
            aln.clear();
            return;
        }
        for (int32_t k = 0; k < reference_length + query_length - 1; k++){
            // printf("Tile: %d, k: %d, L: %d, U: %d, (%d, %d)\n", tile, k, L[k%3], U[k%3]+1, reference_length, query_length);
            if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                last_tile = true;
                errorType = 1;
                aln.clear();
                return;
            }
            
            if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                // fprintf(stderr, "ERROR: anti-diagonal larger than the max limit!\n");
                last_tile = true;
                errorType = 2;
                aln.clear();
                return;
            }

            if (k <= marker) {
                ftr_length.push_back(U[k%3] - L[k%3] + 1);
                ftr_lower_limit.push_back(L[k%3]);
                ftr_addr += U[k%3] - L[k%3] + 1;
            }

            // if (tile == 0) printf("k:%d, i_st: %d, i_en: %d\n",k, L[k%3], U[k%3]+1);
            for (int32_t i = L[k%3]; i < U[k%3]+1; i++) { // i-> query_idx, j -> reference_idx
                int32_t Lprime = std::max(0, static_cast<int32_t>(k)-static_cast<int32_t>(reference_length) + 1); 
                int32_t j = std::min(static_cast<int32_t>(k), static_cast<int32_t>(reference_length - 1)) - (i-Lprime); 
                
                if (j < 0) {
                    fprintf(stderr, "ERROR: j less than 0.\n");
                    exit(1);
                }

                int32_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                int32_t offset = i-L[k%3];
                int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                

                
                // int score_from_prev_tile = 0;
                if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                    // if (k==0 && tile>0)
                    // {
                    //     score_from_prev_tile = tile*10;
                    // }
                    int32_t similarScore = 0;
                    // float denominator = 0;
                    float numerator = 0;
                    for (int l = 0; l < 6; ++l) {
                        for (int m = 0; m < 6; ++m) {
                            // denominator += reference[reference_idx+j][l]*query[query_idx+i][m];
                            if (m == 5 && l == 5)      numerator += 0;
                            else if (m == 5 || l == 5) numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*gapExtend;
                            else                       numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*scoreMat[m*5+l];
                        }
                    }  

                    similarScore = static_cast<int32_t>(roundeven(numerator/denominator));
                    // if (reference.size() == 1558 && query.size() == 1552 && (refNum > 1 || qryNum > 1)) printf("%d: (%d,%d), %d, %f\n", k, i, j, similarScore, numerator/denominator);
                    if (offsetDiag < 0) match = similarScore;
                    else                match = S[(k+1)%3][offsetDiag] + similarScore;
                    // match = S[(k+1)%3][offsetDiag] + similarScore + score_from_prev_tile;
                }

                int32_t pos_gapOpen_ref =   static_cast<int32_t>(roundeven(gapOp[0][reference_idx+j]));
                int32_t pos_gapOpen_qry =   static_cast<int32_t>(roundeven(gapOp[1][query_idx+i]));
                int32_t pos_gapExtend_ref = static_cast<int32_t>(roundeven(gapEx[0][reference_idx+j]));
                int32_t pos_gapExtend_qry = static_cast<int32_t>(roundeven(gapEx[1][query_idx+i]));
                
                if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                    // delOp = S[(k+2)%3][offsetUp] + gapOpen;
                    delOp = S[(k+2)%3][offsetUp] + pos_gapOpen_ref;
                    // delExt = D[(k+1)%2][offsetUp] + gapExtend;
                    delExt = D[(k+1)%2][offsetUp] + pos_gapExtend_ref;
                }

                if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                    // insOp = S[(k+2)%3][offsetLeft] + gapOpen;
                    insOp = S[(k+2)%3][offsetLeft] + pos_gapOpen_qry;
                    // insExt = I[(k+1)%2][offsetLeft] + gapExtend;
                    insExt = I[(k+1)%2][offsetLeft] + pos_gapExtend_qry;
                }

                
                I[k%2][offset] = insOp;
                D[k%2][offset] = delOp;

                Iptr = false;
                Dptr = false;

                if (insExt >= insOp) {
                    I[k%2][offset] = insExt;
                    Iptr = true;
                }
                if (delExt >= delOp) {
                    D[k%2][offset] = delExt;
                    Dptr = true;
                }

                if (match > I[k%2][offset]) {
                    if (match > D[k%2][offset]) {
                        S[k%3][offset] = match;
                        ptr = 0;
                    }
                    else {
                        S[k%3][offset] = D[k%2][offset];
                        ptr = 2;
                    }
                }
                else if (I[k%2][offset] > D[k%2][offset]) {
                    S[k%3][offset] = I[k%2][offset];
                    ptr = 1;
                }
                else {
                    S[k%3][offset] = D[k%2][offset];
                    ptr = 2;
                }
                
                if (S[k%3][offset] < max_score-param.xdrop) {
                    S[k%3][offset] = -inf;
                }

                score = S[k%3][offset];

                // if (reference.size() == 1558 && query.size() == 1552 && (refNum > 1 || qryNum > 1)) printf("%d: (%d,%d), H: %d, I: %d, D: %d, ptr: %d\n", k, i, j, S[k%3][offset], I[k%2][offset], D[k%2][offset], (ptr&0xFFFF));

                if (max_score_prime < score) {
                    max_score_prime = score;
                }
                

                if (k == marker - 1) { // Convergence algorithm
                    CS[k%3][offset] = (3 << 16) | (i & 0xFFFF); 
                    // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  CS[k%3][offset] << "\n";
                } else if (k == marker) {
                    CS[k%3][offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                    CI[k%2][offset] = (1 << 16) | (i & 0xFFFF);
                    CD[k%2][offset] = (2 << 16) | (i & 0xFFFF);
                    // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  CS[k%3][offset] <<  " " << CI[k%2][offset] <<  " " << CD[k%2][offset] << "\n";
                    // if (score > max_score_marker) {
                    //     max_score_marker = score;
                    //     max_score_marker_ref_idx = j;
                    //     max_score_marker_query_idx = i;
                    //     max_score_marker_start_addr = ftr_addr - (U[k%3] - L[k%3] + 1)  + (i - L[k%3]);
                    //     max_score_marker_start_ftr = k;
                    //     tb_state_marker = ptr;
                    // }
                } 
                else if (k >= marker + 1){
                    if (Iptr) {
                        CI[k%2][offset] = CI[(k+1)%2][offsetLeft]; 
                    } else {
                        CI[k%2][offset] = CS[(k+2)%3][offsetLeft]; 
                    }

                    if (Dptr) {
                        CD[k%2][offset] = CD[(k+1)%2][offsetUp];
                    } else {
                        CD[k%2][offset] = CS[(k+2)%3][offsetUp];
                    }

                    if (ptr == 0) {
                        CS[k%3][offset] = CS[(k+1)%3][offsetDiag];
                    } else if (ptr == 1) {
                        CS[k%3][offset] = CI[k%2][offset];
                    } else {
                        CS[k%3][offset] = CD[k%2][offset];
                    } 
                }
                // if(tile == 0) printf("k: %d, Convergence Unique Id's: %d, %d, %d\n", k, CS[(k%3)*fLen+offset], CI[(k%2)*fLen+offset], CD[(k%2)*fLen+offset]);
                // if (k < 5)
                // {
                //     std::cout << "Index: " << i << ", " << j << " Ins:" << insOp << ", " << insExt << " match:"<< match << " Del:" << delOp << ", " << delExt << " " << S[(k+2)%3][offsetUp] << " ptr:" << (ptr&0xFFFF) << " " << Iptr<< std::endl;
                // }
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
                if (k <= marker){
                    tb.push_back(ptr);
                    // std::cout << (ptr & 0xFFFF) << " ";
                }
            }
            
            // std::cout << "\n";
            
            int32_t newL = L[k%3];
            int32_t newU = U[k%3];


            while (newL <= U[k%3]) {
                int32_t offset = newL - L[k%3];
                if (S[k%3][offset] <= -inf) {
                    newL++;
                }
                else {
                    break;
                }
            }

            while (newU >= L[k%3]) {
                int32_t offset = newU - L[k%3];
                if (S[k%3][offset] <= -inf) {
                    newU--;
                }
                else {
                    break;
                }
            }
            
            if ((!converged) && (k < reference_length + query_length - 2)) {
                int32_t conv_I = Reduction_tree(CI[k%2], newL - L[k%3], newU - newL);
                int32_t conv_D = Reduction_tree(CD[k%2], newL - L[k%3], newU - newL);
                int32_t conv_S = Reduction_tree(CS[k%3], newL - L[k%3], newU - newL);
                // if (tile == 100 ) printf("k: %d, CONV: %d, %d, %d\n", k, conv_I, conv_D, conv_S);
                if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                    converged = true; 
                    conv_value = prev_conv_s;
                    conv_score = max_score_prime;
                    // if (DEBUG)  std::cout << "Converged at: " << conv_value << "\n";
                    // std::cout << "Converged at: " << conv_value << "\n";
                }
                prev_conv_s = conv_S;
            }


            int32_t v1 = static_cast<int32_t>(query_length)-1;
            int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
            int32_t v3 = newU+1;

            int32_t Lprime = std::max(static_cast<int32_t>(0), v2);
            
            L[(k+1)%3] = std::max(newL, Lprime);
            U[(k+1)%3] = std::min(v1, v3); 
            
            // Update max_score
            max_score = max_score_prime;
            // if (tile >= 100 && tile <= 103 ) printf("Max score: %d\n", max_score);
            last_k = k;
            if ((converged) && (max_score > conv_score)){
                conv_logic = true;
                // if (DEBUG) std::cout << "Convergence logic found: ";
                // printf("Tile %d Convergence logic found: %d\n", tile, max_score);
                break;
            }
        }
        
        // if (DEBUG) std::cout <<  "Frontier addr: " << ftr_addr << " \ntb_start_ftr: " << ftr_length.size() << "\nmarker: " << params.marker << std::endl;
        // if (tile == 0) std::cout <<  "Frontier addr: " << ftr_addr << " \ntb_start_ftr: " << ftr_length.size() << "\nmarker: " << params.marker << std::endl;
        if (conv_logic) {
            conv_query_idx = conv_value & 0xFFFF;
            tb_state = (conv_value >> 16) & 0xFFFF;
            conv_ref_idx = marker - conv_query_idx; 
            conv_ref_idx -= (tb_state == 3) ? 1: 0;
            tb_start_addr = ftr_addr - ftr_length[ftr_length.size() - 1];
            tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length.size() - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 2]) : tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 1]);
            tb_start_ftr = (tb_state == 3) ? ftr_length.size() - 2: ftr_length.size() - 1;
            // if (DEBUG) std::cout <<  " conv query idx: " << conv_query_idx << " " << (tb_state&0xFFFF) << " " << conv_ref_idx << " " << conv_value << std::endl;
        } else {
            // Global Alignment
            // int32_t last_k = reference_length + query_length - 2;
            // std::cout << "K:"<< last_k <<" Convergence Unique Id's: " <<  CS[last_k%3][0] <<  " " << CI[last_k%2][0] <<  " " << CD[last_k%2][0] << "\n";
            if (last_k < marker) {
                conv_query_idx = query_length-1;
                conv_ref_idx = reference_length-1;
                tb_start_addr = ftr_addr-1;
                tb_start_ftr = last_k;
                tb_state = 0;
                last_tile = true;
            }
            else {
                conv_query_idx = CS[(last_k%3)][0] & 0xFFFF;
                tb_state = (CS[(last_k%3)][0] >> 16) & 0xFFFF;
                conv_ref_idx = marker - conv_query_idx; 
                conv_ref_idx -= (tb_state == 3) ? 1: 0;
                tb_start_addr = ftr_addr - ftr_length[ftr_length.size() - 1];
                tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length.size() - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 2]) : 
                                                  tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 1]);
                tb_start_ftr = (tb_state == 3) ? ftr_length.size() - 2: ftr_length.size() - 1;
            }
            // if (tile >= 462) printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", max_score_query_idx, max_score_ref_idx, max_score_start_addr, max_score_start_ftr);
            // if (tile >= 462) printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", conv_query_idx, conv_ref_idx, tb_start_addr , tb_start_ftr);
            // if (tile >= 462) printf("MAX: qryLen: %d, refLen: %d\n", query_length, reference_length);
            // conv_query_idx = max_score_query_idx;
            // conv_ref_idx = max_score_ref_idx;
            // tb_start_addr = max_score_start_addr;
            // tb_start_ftr = max_score_start_ftr;
            // tb_state = 0;
            // printf("MAX: qry: %d, ref: %d, addr:%d, ftr:%d\n", max_score_query_idx, max_score_ref_idx, max_score_start_addr, max_score_start_ftr);
            // last_tile = true;
        }

        reference_idx += conv_ref_idx;
        query_idx += conv_query_idx;

        if (reference_idx == reference.size()-1 && query_idx < query.size()-1) { // dir == 1
            for (int q = 0; q < query.size()-query_idx-1; ++q) aln.push_back(static_cast<int16_t>(1));
            last_tile = true;
        }
        if (query_idx == query.size()-1 && reference_idx < reference.size()-1) { // dir == 2
            for (int r = 0; r < reference.size()-reference_idx-1; ++r) aln.push_back(static_cast<int16_t>(2));
            last_tile = true;
        }
                
        // if (conv_query_idx < 10) printf("TB: tile: %d, refidx: %d, qryidx:%d, last:%d\n", tile, reference_idx, query_idx, last_tile);
        // if (conv_query_idx < 10) printf("TB: tile: %d, refLen: %d, qryLen:%d, last:%d\n", tile, reference_length, query_length, last_tile);
        // if (DEBUG) std::cout <<  "tb_start_addr: " << tb_start_addr << " \ntb_start_ftr: " << tb_start_ftr << std::endl;
        bool firstTile = (tile == 0);
        Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln, firstTile);
        state = tb_state%3;
        // printf("conv_idx: %d, %d\n", conv_ref_idx, conv_query_idx);
        
        
        

        // Deallocate memory
        for (size_t sIndx=0; sIndx<3; sIndx++) {
            std::free(S[sIndx]);
            std::free(CS[sIndx]);
            if (sIndx < 2) { 
                std::free(I[sIndx]);
                std::free(D[sIndx]);
                std::free(CI[sIndx]);
                std::free(CD[sIndx]);
            }
        }

    }


/*

int Talco_xdrop::Score (
    Params params,
    const std::vector<int8_t> &aln,
    const std::string &reference,
    const std::string &query,
    const int ref_idx,
    const int query_idx
){  
    int score = 0;
    int ridx = ref_idx, qidx = query_idx;
    int m = 0,mm= 0,go=0,ge=0;
    // printf ("aln: %d %d %d\n",aln.size(), ref_idx, query_idx);
    for (int i = aln.size() - 1; i >= 0; i--){
        int8_t state = aln[i];	
        // printf ("state %d, ri %d, qi %d", state, ridx, qidx);
        if (state == 0) { // Current State M	
            if (reference[ridx] == query[qidx]){	
                score += params.match;
                m++;	
            } else {	
                score += params.mismatch;
                mm++;	
            }	
            ridx--;	
            qidx--;	
        } else if (state == 1) { // Current State I	
            if (i < aln.size() - 1 && aln[i+1] == 1) {	
                score += params.gapExtend;
                ge++;	
            } else {	
                score += params.gapOpen;
                go++;	
            }	
            qidx--;	

        } else { // Current State D	
            if (i < aln.size() - 1 && aln[i+1] == 2) {	
                score += params.gapExtend;	
                ge++;	
            } else {	
                score += params.gapOpen;	
                go++;	
            }	
            ridx--;	
        }	
        // printf(" Score: %d\n", score);
    }
    // printf("count %d %d %d %d, ridx %d, qidx %d\n", m, mm, go ,ge, ridx, qidx);

    return score;
    
}

void Talco_xdrop::Align (
    Params params,
    const std::vector<std::string>& reference,
    const std::vector<std::string>& query,
    std::vector<int8_t>& aln
    // size_t num_alignments
    ) {
        
        clock_t start, end;
        double time = 0;
        // int count = 0;
        // n -> alignment number
        // for (size_t n = 0; n < num_alignments; n++) {
        
            // std::vector<int8_t> aln;
            int8_t state = 0;
            int score;
            // std::cout << reference[n] << " " << query[n] << " " << std::endl;

            // initialising variables
            int32_t reference_idx, query_idx;
            reference_idx = 0; query_idx = 0;
            bool last_tile = false;
            
            std::vector<std::vector<int>> freqRef;
            std::vector<std::vector<int>> freqQry;
            for (int c = 0; c < std::max(reference[0].size(), query[0].size()); ++c) {
                std::vector<int> temp;
                for (int f = 0; f < 6; ++f) temp.push_back(0);
                if (c < reference[0].size()) freqRef.push_back(temp);
                if (c < query[0].size()) freqQry.push_back(temp);
            }
            for (auto r: reference) {
                for (int c = 0; c < r.size(); ++c) {
                    if (r[c] == 'A')      freqRef[c][0] += 1;
                    else if (r[c] == 'C') freqRef[c][1] += 1;
                    else if (r[c] == 'G') freqRef[c][2] += 1;
                    else if (r[c] == 'T') freqRef[c][3] += 1;
                    else if (r[c] == 'N') freqRef[c][4] += 1;
                    else                  freqRef[c][5] += 1;
                }
            }
            for (auto q: query) {
                for (int c = 0; c < q.size(); ++c) {
                    if (q[c] == 'A')      freqQry[c][0] += 1;
                    else if (q[c] == 'C') freqQry[c][1] += 1;
                    else if (q[c] == 'G') freqQry[c][2] += 1;
                    else if (q[c] == 'T') freqQry[c][3] += 1;
                    else if (q[c] == 'N') freqQry[c][4] += 1;
                    else                  freqQry[c][5] += 1;
                }
            } 
            // std::cout << freqRef.size() << ',' << freqQry.size() << '\n';
            int tile = 0;
            while (!last_tile) {
                std::vector<int8_t> tile_aln;
                // Talco_xdrop::Tile(reference, query, params, reference_idx, query_idx, tile_aln, state, last_tile, tile);
                Talco_xdrop::Tile(freqRef, freqQry, params, reference_idx, query_idx, tile_aln, state, last_tile, tile);
                // printf("Tile: %d, (r,q) = (%d,%d), tile_aln_size: %d\n",tile, reference_idx, query_idx, tile_aln.size());
                for (int i= tile_aln.size()-1; i>=0; i--){
                    if (i == tile_aln.size()-1 && tile>0){
                        continue;
                    }
                    // std::cout << i << "\n";
                    aln.push_back(tile_aln[i]);
                }
                
                tile_aln.clear();
                
                tile++;
                // aln.clear();
            }
            // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
            
            // for (int i = 0; i < aln.size(); ++i) {
            //     // std::cout << aln[i];
            //     if (aln[i] & 0xFFFF != 0) 
            // }
        // int count = 0; 
        // for (auto &a: aln){
        //     if ((a & 0xFFFF) != 0) std::cout << "TB: " << (a & 0xFFFF) << " at " << count << '\n';
        //     count += 1;
        // }
        // printf("Alignment (Length, Score): (%d, %d)\n ", aln.size(), score);
        // break;

        // }
        
        // printf("%f\n", time);
        // printf("%d\n", count);

    }
*/