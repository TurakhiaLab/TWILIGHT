/*
    MIT License

    Copyright (c) 2023 Turakhia Lab

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>

#define I_BOUNDARY -2
#define D_BOUNDARY -3

Talco_xdrop::Params::Params(msa::Params& param) {
    this->matrixSize = param.matrixSize;
    this->scoreMatrix = new float * [this->matrixSize];
    for (int i = 0; i < matrixSize; ++i) this->scoreMatrix[i] = new float [matrixSize];
    for (int i = 0; i < matrixSize; ++i) {
        for (int j = 0; j < matrixSize; ++j) {
            this->scoreMatrix[i][j] = param.scoringMatrix[i][j];
        }
    }
    this->gapOpen = param.gapOpen;
    this->gapExtend = param.gapExtend;
    this->gapCharScore = param.gapExtend;
    this->gapBoundary = param.gapBoundary;
    this->xdrop = static_cast<int32_t> (1000 * -1 * this->gapExtend);
    this->fLen = (1 << 12);
    this->marker = (1 << 10); //reduce this value to save memory
    this->alnType = 0;
}

Talco_xdrop::Params::~Params() {
    for (int i = 0; i < this->matrixSize; ++i) delete [] this->scoreMatrix[i];
    delete [] this->scoreMatrix;
}

using scoreType = float;

void Talco_xdrop::Align_freq (
    Params* params,
    const std::vector<std::vector<float>>& freqRef,
    const std::vector<std::vector<float>>& freqQry,
    const std::vector<std::vector<float>>& gapOp,
    const std::vector<std::vector<float>>& gapEx,
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
                // if (freqRef.size() == 1776 && freqQry.size() == 6124) std::cout << (tile_aln[i] & 0xFFFF);
                aln.push_back(tile_aln[i]);
            }
            // if (freqRef.size() == 1776 && freqQry.size() == 6124) std::cout << '\n';
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

inline void freeMem(scoreType** S, scoreType** I, scoreType** D, int32_t** CS, int32_t** CI, int32_t** CD) {
    for (size_t sIndx=0; sIndx<3; sIndx++) {
        delete [] S[sIndx];
        delete [] CS[sIndx];
        if (sIndx < 2) { 
            delete [] I[sIndx];
            delete [] D[sIndx];
            delete [] CI[sIndx];
            delete [] CD[sIndx];
        }
    }
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
        if (ftr > 0) addr = addr - (idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);
        if (dir == 0){
            if (ftr > 1) addr = addr - (ftr_length[ftr - 2]) + (idx - ftr_lower_limit[ftr  - 2]);
            ftr -= 2;
            idx -= 1;
            query_idx--;
            ref_idx--;
        }else if (dir == 1){
            if (ftr > 0) addr = addr + (idx - ftr_lower_limit[ftr  - 1]);
            ftr -= 1;
            idx -=1;
            query_idx--;
        }else{
            if (ftr > 0) addr = addr + (idx - ftr_lower_limit[ftr  - 1] + 1);
            ftr -= 1;
            ref_idx--;
        }
        aln.push_back(dir);   
        if (firstTile &&  (ref_idx < 0 || query_idx < 0)) break;
    }
    if (firstTile) {
        while (ref_idx > -1) {
            aln.push_back(2);
            ref_idx--;
        }
        while (query_idx > -1) {
            aln.push_back(1);
            query_idx--;
        }
    }
}

void Talco_xdrop::Tile (
    // const std::vector<std::string>& reference,
    // const std::vector<std::string>& query,
    const std::vector<std::vector<float>>& reference,
    const std::vector<std::vector<float>>& query,
    const std::vector<std::vector<float>>& gapOp,
    const std::vector<std::vector<float>>& gapEx,
    const std::pair<float, float>& num,
    Params* param,
    int32_t &reference_idx,
    int32_t &query_idx,
    std::vector<int8_t> &aln,
    int8_t &state,
    bool &last_tile,
    const int &tile,
    int16_t& errorType // 1: xdrop, 2: exceed anti-doagonal limit
    ) {
                
        // Initialising variables
        scoreType inf = 2.0*param->xdrop + 1.0;
        int32_t marker = param->marker;
        bool converged = false, conv_logic = false;
        float refNum = num.first, qryNum = num.second; 
        int32_t reference_length = reference.size() - reference_idx; 
        int32_t query_length = query.size() - query_idx;
        int32_t fLen = std::min(param->fLen, std::min(reference_length, query_length));
        scoreType max_score = 0; scoreType max_score_prime = -inf;
        scoreType conv_score = 0; int32_t conv_value = 0; int32_t conv_ref_idx = 0; int32_t conv_query_idx = 0; 
        int32_t tb_start_addr = 0; int32_t tb_start_ftr = 0; // int32_t max_score_start_addr = 0; int32_t max_score_start_ftr = 0;
        int8_t tb_state = 0;
        
        // For xdrop
        // int32_t max_score_marker = -inf; int32_t max_score_marker_ref_idx = 0; int32_t max_score_marker_query_idx = 0;
        // int32_t max_score_marker_start_addr = 0; int32_t max_score_marker_start_ftr = 0;
        // int8_t tb_state_marker = 0;

        float denominator = refNum * qryNum;

        bool type = (reference[0].size() == 6) ? 0 : 1; // 0: dna/rna, 1: protein
        float gapOpen = param->gapOpen;
        float gapExtend = param->gapExtend;
        float gapOpenAtEnds = (param->alnType == 0) ? gapOpen : 0;
        float gapExtendAtEnds = (param->alnType == 0) ? gapExtend : 0;
        
        int32_t L[3], U[3];
        scoreType *S[3], *I[2], *D[2];
        int32_t *CS[3], *CI[2], *CD[2];
        std::vector<int8_t> tb;
        std::vector<int32_t> ftr_length;
        std::vector<int32_t> ftr_lower_limit;
        int32_t ftr_addr = 0;
        int32_t last_k = 0;
        int32_t prev_conv_s = -1;

        for (size_t sIndx=0; sIndx<3; sIndx++) { // Allocate memory for S, I, D, and CS array
            S[sIndx] = new scoreType [fLen];
            CS[sIndx] = new int32_t [fLen];
            if (sIndx < 2) {
                I[sIndx] = new scoreType [fLen];
                D[sIndx] = new scoreType [fLen];
                CI[sIndx] = new int32_t [fLen];
                CD[sIndx] = new int32_t [fLen];
            }
            L[sIndx] = sIndx;
            U[sIndx] = -sIndx;
        }
        
        for (int32_t sIndx=0; sIndx<3; sIndx++){ // Initialise memory for S, I, D array
            for (int32_t sLenIndex=0; sLenIndex<fLen; sLenIndex++){
                S[sIndx][sLenIndex] = -1;
                CS[sIndx][sLenIndex] = -1;
                if (sIndx < 2) {  
                    I[sIndx][sLenIndex] = -1;
                    D[sIndx][sLenIndex] = -1;
                    CI[sIndx][sLenIndex] = I_BOUNDARY;
                    CD[sIndx][sLenIndex] = D_BOUNDARY;
                }
            }
        }

        if ((reference_length < 0) || (query_length < 0)) {
            // std::cout << reference_length << " " << query_length << std::endl;
            fprintf(stderr, "ERROR: Reference/Query index exceeded limit!\n");
            errorType = 3;
            aln.clear();
            freeMem(S, I, D, CS, CI, CD);
            return;
        }
        for (int32_t k = 0; k < reference_length + query_length - 1; k++){
            // if (reference.size() == 1776 && query.size() == 6124) printf("Tile: %d, k: %d, L: %d, U: %d, (%d, %d)\n", tile, k, L[k%3], U[k%3]+1, reference_length, query_length);
            if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                last_tile = true;
                errorType = 1;
                aln.clear();
                freeMem(S, I, D, CS, CI, CD);
                return;
            }
            
            if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                // fprintf(stderr, "ERROR: anti-diagonal larger than the max limit!\n");
                last_tile = true;
                errorType = 2;
                aln.clear();
                freeMem(S, I, D, CS, CI, CD);
                return;
            }

            if (k <= marker) {
                ftr_length.push_back(U[k%3] - L[k%3] + 1);
                ftr_lower_limit.push_back(L[k%3]);
                ftr_addr += U[k%3] - L[k%3] + 1;
            }

            // std::vector<scoreType> max_scores (U[k%3]+1-L[k%3], max_score_prime);
            // int tb_size = tb.size();
            // if (k <= marker) tb.resize(tb_size+U[k%3]+1-L[k%3], -1);

            // tbb::this_task_arena::isolate([&] { 
            // tbb::parallel_for(tbb::blocked_range<int>(L[k%3], U[k%3]+1), [&](tbb::blocked_range<int> range) {
            // for (int32_t i = range.begin(); i < range.end(); ++i) {
            for (int32_t i = L[k%3]; i < U[k%3]+1; i++) { // i-> query_idx, j -> reference_idx
                int8_t ptr = 0;  // Main pointer
                scoreType score = 0;
                bool Iptr = false;
                bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
                int32_t Lprime = std::max(0, static_cast<int32_t>(k)-static_cast<int32_t>(reference_length) + 1); 
                int32_t j = std::min(static_cast<int32_t>(k), static_cast<int32_t>(reference_length - 1)) - (i-Lprime); 
                if (j < 0) {
                    fprintf(stderr, "ERROR: j less than 0.\n");
                    exit(1);
                }
                scoreType match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                int32_t offset = i-L[k%3];
                int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                if ((k==0) || 
                    ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3])) ||
                    (tile == 0 && (i == 0 || j == 0 ))) {
                    scoreType similarScore = 0;
                    float numerator = 0;
                    if (type == 0) {
                        for (int l = 0; l < 6; ++l) {
                            for (int m = 0; m < 6; ++m) {
                                if (m == 5 && l == 5)      numerator += 0;
                                else if (m == 5 || l == 5) numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*param->gapCharScore;
                                else                       numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*param->scoreMatrix[m][l];
                            }
                        }
                    }
                    else {
                        for (int l = 0; l < 22; ++l) {
                            for (int m = 0; m < 22; ++m) {
                                if (m == 21 && l == 21)      numerator += 0;
                                else if (m == 21 || l == 21) numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*gapExtend;
                                else                         numerator += reference[reference_idx+j][l]*query[query_idx+i][m]*param->scoreMatrix[m][l];
                            }
                        }
                    }
                    similarScore = numerator/denominator;
                    if  (tile == 0 && (i == 0 || j == 0 )) {
                        if (i == 0 && j == 0) match = similarScore;
                        else                  match = similarScore + gapOpenAtEnds + gapExtendAtEnds * (std::max(0, std::max(reference_idx + j, query_idx + i) - 1));
                    }
                    else if (offsetDiag < 0)               match = similarScore;
                    else                                   match = S[(k+1)%3][offsetDiag] + similarScore;
                }
                scoreType pos_gapOpen_ref =  gapOp[0][reference_idx+j];
                scoreType pos_gapOpen_qry =  gapOp[1][query_idx+i];
                scoreType pos_gapExtend_ref = gapEx[0][reference_idx+j];
                scoreType pos_gapExtend_qry = gapEx[1][query_idx+i];
                if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                    delOp =  S[(k+2)%3][offsetUp] + pos_gapOpen_ref;
                    delExt = D[(k+1)%2][offsetUp] + pos_gapExtend_ref;
                }
                if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                    insOp =  S[(k+2)%3][offsetLeft] + pos_gapOpen_qry;
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
                // if (tile == 0 && param->alnType == 0 && k < 20) std::cout << '(' << i << ',' << j << ')' << match << ',' << I[k%2][offset] << ',' << D[k%2][offset] << '\n';
                if (match >= I[k%2][offset]) {
                    if (match >= D[k%2][offset]) {
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
                if (S[k%3][offset] < max_score-param->xdrop) {
                    S[k%3][offset] = -inf;
                }
                
                score = S[k%3][offset];
                
                if (max_score_prime < score) {
                    max_score_prime = score;
                }
                // if (max_scores[offset] < score) max_scores[offset] = score;

                // if (i == L[k%3] && reference.size() == 3304 && query.size() == 1948) {
                // if (pos_gapExtend_ref < -100000 || pos_gapExtend_qry < -100000 || pos_gapOpen_ref < -100000 || pos_gapOpen_qry < -100000) {
                //     fprintf(stderr, "ERROR: Extremely large gap penalties!\n");
                //     std::cout << k << ',' << i << ',' << j << ',' << pos_gapOpen_ref << ',' << pos_gapExtend_ref << ',' << pos_gapOpen_qry << ',' << pos_gapExtend_qry << '\n';
                //     // exit(1);
                // }
                // }
                // if (i == L[k%3] && reference.size() == 3304 && query.size() == 1948) {
                //     std::cout << k << ',' << i << ',' << j << ',' << match << ',' << I[k%2][offset] << ',' << D[k%2][offset] << ',' << S[k%3][offset] << ',' << max_score_prime << '\n';
                //     std::cout << pos_gapOpen_ref << ',' << pos_gapExtend_ref << ',' << pos_gapOpen_qry << ',' << pos_gapExtend_qry << '\n';
                //     std::cout << offset << ',' << offsetDiag << ',' << offsetUp << ',' << offsetLeft << '\n';
                // }
                

                if (k == marker - 1) { // Convergence algorithm
                    CS[k%3][offset] = (3 << 16) | (i & 0xFFFF); 
                } else if (k == marker) {
                    CS[k%3][offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                    CI[k%2][offset] = (1 << 16) | (i & 0xFFFF);
                    CD[k%2][offset] = (2 << 16) | (i & 0xFFFF);
                } 
                else if (k >= marker + 1){
                    if (Iptr) {
                        CI[k%2][offset] = (offsetLeft >= 0) ? CI[(k+1)%2][offsetLeft] : I_BOUNDARY; 
                    } else {
                        CI[k%2][offset] = (offsetLeft >= 0 && CS[(k+2)%3][offsetLeft] != -1) ? CS[(k+2)%3][offsetLeft] : I_BOUNDARY; 
                    }

                    if (Dptr) {
                        CD[k%2][offset] = (offsetUp >= 0) ? CD[(k+1)%2][offsetUp] : D_BOUNDARY;
                    } else {
                        CD[k%2][offset] = (offsetUp >= 0 && CS[(k+2)%3][offsetUp] != -1) ? CS[(k+2)%3][offsetUp] : D_BOUNDARY;
                    }

                    if (ptr == 0) {
                        CS[k%3][offset] = CS[(k+1)%3][offsetDiag];
                    } else if (ptr == 1) {
                        CS[k%3][offset] = CI[k%2][offset];
                    } else {
                        CS[k%3][offset] = CD[k%2][offset];
                    } 
                }
                if (Iptr) {
                    ptr |= 0x04; 
                }
                if (Dptr) {
                    ptr |= 0x08;
                }
                if (k <= marker){
                    // tb[tb_size+offset] = ptr;
                    tb.push_back(ptr);
                }
            }
            // });
            // });
            // max_score_prime = *std::max_element(max_scores.begin(), max_scores.end());
            
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
                if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                    converged = true; 
                    conv_value = prev_conv_s;
                    conv_score = max_score_prime;
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
            max_score =(max_score_prime < 0) ? 0 : max_score_prime;
            last_k = k;
            if ((converged) && (max_score > conv_score)){
                conv_logic = true;
                break;
            }
        }
        
        if (conv_logic) {
            conv_query_idx = conv_value & 0xFFFF;
            tb_state = (conv_value >> 16) & 0xFFFF;
            conv_ref_idx = marker - conv_query_idx; 
            conv_ref_idx -= (tb_state == 3) ? 1: 0;
            tb_start_addr = ftr_addr - ftr_length[ftr_length.size() - 1];
            tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length.size() - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 2]) : tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit.size() - 1]);
            tb_start_ftr = (tb_state == 3) ? ftr_length.size() - 2: ftr_length.size() - 1;
        } else {
            // Global Alignment
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
        }

        if (conv_query_idx == (D_BOUNDARY & 0xFFFF)) {
            conv_query_idx = 0;
            conv_ref_idx = param->marker;
        }
        else if (conv_query_idx == (I_BOUNDARY & 0xFFFF)) {
            conv_query_idx = param->marker;
            conv_ref_idx = 0;
        }
            
        reference_idx += conv_ref_idx;
        query_idx += conv_query_idx;

        reference_length = reference.size() - reference_idx; 
        query_length = query.size() - query_idx;
        if ((reference_length < 0) || (query_length < 0)) {
            fprintf(stderr, "ERROR: Reference/Query index exceeded limit!\n");
            std::cerr << reference_length << ',' << reference_idx << ',' << conv_ref_idx << '\n';
            std::cerr << query_length << ',' << query_idx << ',' << conv_query_idx << '\n';
            
            errorType = 3;
            aln.clear();
            freeMem(S, I, D, CS, CI, CD);
            return;
        }
        

        if (reference_idx == reference.size()-1 && query_idx < query.size()-1) { // dir == 1
            for (int q = 0; q < query.size()-query_idx-1; ++q) aln.push_back(static_cast<int16_t>(1));
            last_tile = true;
        }
        if (query_idx == query.size()-1 && reference_idx < reference.size()-1) { // dir == 2
            for (int r = 0; r < reference.size()-reference_idx-1; ++r) aln.push_back(static_cast<int16_t>(2));
            last_tile = true;
        }
        if (reference_idx == reference.size()-1 && query_idx == query.size()-1) last_tile = true;
                
        bool firstTile = (tile == 0);
        Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln, firstTile);
        state = tb_state%3;

        // Deallocate memory
        freeMem(S, I, D, CS, CI, CD);
        return;

    }