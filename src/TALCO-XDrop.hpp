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
#define TALCO_HPP

#include <stdint.h>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cmath>


namespace Talco_xdrop {
    struct Params {
        float** scoreMatrix;
        float gapOpen;
        float gapBoundary;
        float gapExtend;
        int32_t matrixSize;
        int32_t xdrop;
        int32_t fLen;
        int32_t marker;
        int type;

        void updateXDrop(int32_t new_xdrop) { 
            this->xdrop = new_xdrop;
        }
        void updateFLen(int32_t new_flen) { this->fLen = new_flen;}

        Params(float* t_param, int matrixSize) {
            this->matrixSize = matrixSize;
            this->scoreMatrix = new float * [matrixSize];
            for (int i = 0; i < matrixSize; ++i) this->scoreMatrix[i] = new float [matrixSize];
            for (int i = 0; i < matrixSize; ++i) {
                for (int j = 0; j < matrixSize; ++j) {
                    this->scoreMatrix[i][j] = t_param[i*matrixSize+j];
                }
            }
            this->gapOpen = t_param[matrixSize*matrixSize];
            this->gapExtend = t_param[matrixSize*matrixSize+1];
            this->gapBoundary = t_param[matrixSize*matrixSize+2];
            this->xdrop = static_cast<int32_t> (1000 * -1 * this->gapExtend);
            this->fLen = (1 << 12);
            this->marker = (1 << 10); //reduce this value to save memory
        }
        ~Params() {
            for (int i = 0; i < this->matrixSize; ++i) delete [] this->scoreMatrix[i];
            delete [] this->scoreMatrix;
        }
    };
    
    void Align_freq (
        Params* params,
        const std::vector<std::vector<float>>& freqRef,
        const std::vector<std::vector<float>>& freqQry,
        const std::vector<std::vector<float>>& gapOp,
        const std::vector<std::vector<float>>& gapEx,
        const std::pair<float, float>& num,
        std::vector<int8_t>& aln,
        int16_t& errorType
    );

    void Tile (
        const std::vector<std::vector<float>>& reference, 
        const std::vector<std::vector<float>>& query, 
        const std::vector<std::vector<float>>& gapOp,
        const std::vector<std::vector<float>>& gapEx,
        const std::pair<float, float>& num,
        Params* params,
        int32_t &reference_idx,
        int32_t &query_idx,
        std::vector<int8_t> &aln,
        int8_t &state,
        bool &last_tile,
        const int &tile,
        int16_t& errorType
    );

    int32_t Reduction_tree (
        const int32_t *C,
        const int32_t L,
        const int32_t U
    );
    void Traceback(
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
    );
    
}

#endif