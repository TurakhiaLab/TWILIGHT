#include "TALCO-XDrop.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <ctime>
#include <unordered_map>
#include <cmath>

#define DEBUG false

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
    std::vector<int8_t> &aln
){
    int32_t addr = tb_start_addr; 
    int16_t ftr = tb_start_ftr;
    int16_t idx = tb_start_idx;
    int16_t query_idx = tb_start_idx;
    int16_t ref_idx = ref_start_idx;
    int8_t  state = tb_state;
    int8_t  tb_value = 0;

    int8_t  dir = 0;
    bool checkpoint = false;
    // int count = 0;
    while (ftr >= 0) {
        // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
        if (addr < 0) {
            fprintf(stderr, "ERROR: tb addr < 0!\n");
            // exit(1);
        }
        // count++;
        
        tb_value = tb[addr];
        if (DEBUG)
        {
            std::cout << "Start Addr:" << addr << " state: " << (state & 0xFFFF) << " ,ftr: " << ftr << " ,idx: " << idx << " ,ll[ftr-1]: " <<  ftr_lower_limit[ftr];
            std::cout << " fL[ftr - 1]: " << ftr_length[ftr - 1] << " ,ll[ftr-2]: " <<  ftr_lower_limit[ftr-2];
            std::cout << " fL[ftr - 2]: " << ftr_length[ftr - 2];
            std::cout << " Tb: " << ( tb_value&0xFFFF);
        }

        
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
        // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d, tb_value: %d, dir: %d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr], tb_value, dir);
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

        // if (ref_idx==0 && query_idx==0)
        // {
        //     checkpoint = true;
        // }
        // printf("ftr: %d\n", ftr);
        // printf("addr: %d\n", addr);
        // printf("aln_idx: %d\n", aln.size());
        // printf("dir: %d\n", dir);
        // printf("traceback_idx: %d\n", idx);
        aln.push_back(dir);   
        // std::cout << "aln_idx:" << aln.size() << '\n'; 
        // state = next_state;
        // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
    }
    // std::cout << count << std::endl;   

    if (DEBUG) std::cout << ref_idx << " " << query_idx << std::endl; 
    // std::cout << "aln_idx:" << aln.size() << '\n'; 
    // if (!checkpoint)
    // {

    //     while (ref_idx != -1)
    //     {
    //         state = 2;
    //         aln.push_back(state);
    //         ref_idx--;
    //         std::cout << "Here";
    //     }
    //     while (query_idx != -1)
    //     {
    //         state = 1;
    //         aln.push_back(state);
    //         query_idx--;
    //         std::cout << "Here";
    //     }
    // }

}

void Talco_xdrop::Tile (
    // const std::vector<std::string>& reference,
    // const std::vector<std::string>& query,
    const std::vector<std::vector<int>>& reference,
    const std::vector<std::vector<int>>& query,
    Params params,
    int32_t &reference_idx,
    int32_t &query_idx,
    std::vector<int8_t> &aln,
    int8_t &state,
    bool &last_tile,
    const int &tile
    ) {
        
        // Initialising variables
        int32_t inf = params.xdrop + 1;
        // int32_t fLen = (1 << 10); // frontier length (assuming anti-diagonal length cannot exceed 1024)
        int32_t fLen = (1 << 11); 
        bool converged = false; bool conv_logic = false;
        int32_t reference_length = reference.size() - reference_idx; 
        int32_t query_length = query.size() - query_idx;
        int32_t score = 0; int32_t max_score = 0; int32_t max_score_prime = -inf; int32_t max_score_ref_idx = 0; int32_t max_score_query_idx = 0;
        int32_t conv_score = 0; int32_t conv_value = 0; int32_t conv_ref_idx = 0; int32_t conv_query_idx = 0; 
        int32_t tb_start_addr = 0; int32_t tb_start_ftr = 0; int32_t max_score_start_addr = 0; int32_t max_score_start_ftr = 0;
        int8_t tb_state = 0;
        int8_t ptr = 0;  // Main pointer
        bool Iptr = false;
        bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend

        int32_t L[3], U[3];
        int16_t *S[3], *I[2], *D[2];
        int32_t *CS[3], *CI[2], *CD[2];
        std::vector<int8_t> tb;
        std::vector<int32_t> ftr_length;
        std::vector<int32_t> ftr_lower_limit;
        int32_t ftr_addr = 0;
        int32_t last_k = 0;

        int32_t prev_conv_s = -1;
        
        for (size_t sIndx=0; sIndx<3; sIndx++) { // Allocate memory for S, I, D, and CS array
            S[sIndx] = (int16_t*) std::malloc(fLen*sizeof(int16_t));
            CS[sIndx] = (int32_t*) std::malloc(fLen*sizeof(int32_t));
            if (sIndx < 2) {
                I[sIndx] = (int16_t*) std::malloc(fLen*sizeof(int16_t));
                D[sIndx] = (int16_t*) std::malloc(fLen*sizeof(int16_t));
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
            exit(1); 
        }
        // printf("Tile: %d (%d, %d)\n", tile, reference_length, query_length);
        for (int32_t k = 0; k < reference_length + query_length - 1; k++){
            // printf("Tile: %d, k: %d, L: %d, U: %d, (%d, %d)\n", tile, k, L[k%3], U[k%3]+1, reference_length, query_length);
            if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                std::cout << "No more cells to compute based on x-drop critieria" << std::endl;
                break;
            }
            
            if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                fprintf(stderr, "ERROR: anti-diagonal larger than the max limit!\n");
                last_tile = true;
                break;
            }

            if (k <= params.marker) {
                ftr_length.push_back(U[k%3] - L[k%3] + 1);
                ftr_lower_limit.push_back(L[k%3]);
                ftr_addr += U[k%3] - L[k%3] + 1;
                // std::cout << k << " ftr length: " << U[k%3] - L[k%3] + 1 << " ftr_addr: " << ftr_addr << " ftr lower limit: " << L[k%3] << " tb len: " << tb.size() << std::endl;
            }
            // if (tile == 0) printf("k:%d, i_st: %d, i_en: %d\n",k, L[k%3], U[k%3]+1);
            for (int16_t i = L[k%3]; i < U[k%3]+1; i++) { // i-> query_idx, j -> reference_idx
                int16_t Lprime = std::max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1); 
                int16_t j = std::min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime); 
                
                if (j < 0) {
                    fprintf(stderr, "ERROR: j less than 0.\n");
                    exit(1);
                }

                int32_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                int32_t offset = i-L[k%3];
                int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;

                
                int score_from_prev_tile = 0;
                if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                    if (k==0 && tile>0)
                    {
                        score_from_prev_tile = tile*10;
                    }
                    
                    int32_t similarScore = 0;
                    float denominator = 0;
                    float numerator = 0;
                    for (int l = 0; l < 6; ++l) {
                        for (int m = 0; m < 6; ++m) {
                            denominator += reference[reference_idx+j][l]*query[query_idx+i][m];
                            if (m == 4 || l == 4) numerator += 0;
                            else if (m == l)      numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*params.match;
                            else                  numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*params.mismatch;
                        }
                    }
                    similarScore = static_cast<int32_t>(std::floor(numerator/denominator));
                    // for (int t = 0; t < 6; ++t) std::cout << reference[reference_idx+j][t] << ',';
                    // std::cout << '\n';
                    // for (int t = 0; t < 6; ++t) std::cout << query[query_idx+i][t] << ',';
                    // std::cout << '\n';
                    // printf("(%d, %d) = (%c, %c), score: %d\n",j, i, reference[0][reference_idx+j], query[0][query_idx+i], similarScore);
                    match = S[(k+1)%3][offsetDiag] + similarScore + score_from_prev_tile;
                    // if (reference[0][reference_idx+j] == 'N' || reference[0][reference_idx+j] == 'n' ||
                    //     query[0][query_idx+i] == 'N' || query[0][query_idx+i] == 'n') {
                    //     match = S[(k+1)%3][offsetDiag] + score_from_prev_tile; // match score = 0 
                    // }
                    // else if (reference[0][reference_idx+j] == query[0][query_idx+i]) {
                    //     match = S[(k+1)%3][offsetDiag] + params.match + score_from_prev_tile; 
                    // }
                    // else {
                    //     match = S[(k+1)%3][offsetDiag] + params.mismatch + score_from_prev_tile;
                    // }
                }

                

                if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                    delOp = S[(k+2)%3][offsetUp] + params.gapOpen;
                    // if (k==1 && tile>0 && state==2)
                    // {
                    //     D[(k+1)%2][offsetUp] = val;
                    // }
                    delExt = D[(k+1)%2][offsetUp] + params.gapExtend;
                }

                if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                    insOp = S[(k+2)%3][offsetLeft] + params.gapOpen;
                    // if (k==1 && tile>0 && state==1)
                    // {
                    //     I[(k+1)%2][offsetLeft] = val;
                    // }
                    insExt = I[(k+1)%2][offsetLeft] + params.gapExtend;
                   
                }

                // if (k ==0 && tile > 0)
                // {
                //     if (state == 1)
                //     {
                //         insExt = 10;
                //         val = 10;
                //     }
                //     else if (state == 2)
                //     {
                //         delExt = 10;
                //         val = 10;
                //     }
                // }

                


                I[k%2][offset] =  insOp;
                D[k%2][offset] =  delOp;
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
                
                // if (tile > 0)
                // {
                //     if (i==0)
                // }


                if (S[k%3][offset] < max_score-params.xdrop) {
                    S[k%3][offset] = -inf;
                }

                score = S[k%3][offset];

                // if (tile == 462) printf("(r,q)=(%d,%d) pre_max:%d, score: %d\n", i, j, max_score_prime, score);

                // if (i == 0) printf("k: %d, %d, score: %d\n", k, max_score_prime, score);
                // if (tile == 0 ) if (i == L[k%3]+(U[k%3]-L[k%3])/2) printf("k: %d, (%d, %d), idx: %d, H: %d, D: %d, I: %d, %c, %c\n", k, i, j, i-L[k%3], S[k%3][offset], D[k%2][offset], I[k%2][offset], reference[reference_idx+j], query[query_idx+i]);
                if (max_score_prime < score) {
                    max_score_prime = score;
                    if (k <= params.marker) {
                        max_score_ref_idx = j;
                        max_score_query_idx = i;
                        max_score_start_addr = ftr_addr - (U[k%3] - L[k%3] + 1)  + (i - L[k%3]);
                        max_score_start_ftr = k;
                    }
                    // if (i == U[k%3]) printf("max (%d, %d), %d",max_score_ref_idx, max_score_query_idx, max_score_start_addr);
                }
                

                if (k == params.marker - 1) { // Convergence algorithm
                    CS[k%3][offset] = (3 << 16) | (i & 0xFFFF); 
                    if(DEBUG) std::cout << "Convergence Unique Id's: " <<  CS[k%3][offset] << "\n";
                } else if (k == params.marker) {
                    CS[k%3][offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                    CI[k%2][offset] = (1 << 16) | (i & 0xFFFF);
                    CD[k%2][offset] = (2 << 16) | (i & 0xFFFF);
                    if(DEBUG) std::cout << "Convergence Unique Id's: " <<  CS[k%3][offset] <<  " " << CI[k%2][offset] <<  " " << CD[k%2][offset] << "\n";
                } 
                else if (k >= params.marker + 1){
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
                if (k <= params.marker){
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
            if ((converged) && (max_score > conv_score)){
                conv_logic = true;
                // if (DEBUG) std::cout << "Convergence logic found: ";
                // printf("Tile %d Convergence logic found: %d\n", tile, max_score);
                break;
            }
        }
        // Deallocate memory for scores
        // for (size_t sIndx=0; sIndx<3; sIndx++) {
        //     std::free(S[sIndx]);
        //     std::free(CS[sIndx]);
        //     if (sIndx < 2) {
        //         std::free(I[sIndx]);
        //         std::free(D[sIndx]);
        //         std::free(CI[sIndx]);
        //         std::free(CD[sIndx]);
        //     }
        // }
        if (DEBUG) std::cout <<  "Frontier addr: " << ftr_addr << " \ntb_start_ftr: " << ftr_length.size() << "\nmarker: " << params.marker << std::endl;
        // if (tile == 0) std::cout <<  "Frontier addr: " << ftr_addr << " \ntb_start_ftr: " << ftr_length.size() << "\nmarker: " << params.marker << std::endl;
        if (conv_logic) {
            conv_query_idx = conv_value & 0xFFFF;
            tb_state = (conv_value >> 16) & 0xFFFF;
            conv_ref_idx = params.marker - conv_query_idx; 
            conv_ref_idx -= (tb_state == 3) ? 1: 0;
            tb_start_addr = ftr_addr - ftr_length[ftr_length.size() - 1];
            tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length.size() - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 2]) : tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 1]);
            tb_start_ftr = (tb_state == 3) ? ftr_length.size() - 2: ftr_length.size() - 1;
            if (DEBUG) std::cout <<  " conv query idx: " << conv_query_idx << " " << (tb_state&0xFFFF) << " " << conv_ref_idx << " " << conv_value << std::endl;
        } else {
            // Global Alignment
            // int32_t last_k = reference_length + query_length - 2;
            // std::cout << "K:"<< last_k <<" Convergence Unique Id's: " <<  CS[last_k%3][0] <<  " " << CI[last_k%2][0] <<  " " << CD[last_k%2][0] << "\n";
            if (last_k < params.marker) {
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
                conv_ref_idx = params.marker - conv_query_idx; 
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
        // printf("TB: tile: %d, refidx: %d, qryidx:%d, last:%d\n", tile, reference_idx, query_idx, last_tile);
        if (DEBUG) std::cout <<  "tb_start_addr: " << tb_start_addr << " \ntb_start_ftr: " << tb_start_ftr << std::endl;

        Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln);
        state = tb_state%3;
        
        
        
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