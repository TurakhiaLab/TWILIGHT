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

typedef float paramType;

namespace Talco_xdrop {
    struct Params {
        paramType scoreMatrix [5][5];
        paramType gapOpen;
        paramType gapClose;
        paramType gapExtend;
        int32_t xdrop;
        int32_t fLen;
        int32_t marker;

        void updateXDrop(int32_t new_xdrop) { this->xdrop = new_xdrop;}
        void updateFLen(int32_t new_flen) { this->fLen = new_flen;}

        Params(paramType* t_param) {
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    this->scoreMatrix[i][j] = t_param[i*5+j];
                }
            }
            this->gapOpen = t_param[25];
            this->gapExtend = t_param[26];
            this->gapClose = t_param[27];
            this->xdrop = 1000;
            this->fLen = (1 << 12);
            this->marker = (1 << 11); //reduce this value to save memory
        }
    };
    // struct Params 
    // {
    //     int16_t match;
    //     int16_t mismatch;
    //     int16_t gapOpen;
    //     int16_t gapExtend;

    //     int16_t xdrop;
    //     int16_t marker;
    //     int16_t scoreMode;

    //     int16_t hoxd70 [5][5] = { {  91, -114,  -31, -123, -100},
    //                             {-114,  100, -125,  -31, -100},
    //                             { -31, -125,  100, -114, -100},
    //                             {-123,  -31, -114,   91, -100},
    //                             {-100, -100, -100, -100, -100} }; 
    
    //     int16_t hoxd70_gapOpen = -400;
    //     int16_t hoxd70_gapExtend = -30;

    //     Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen, int16_t t_gapExtend, int16_t t_xdrop, int16_t scoreMode) : 
    //         match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop), scoreMode(scoreMode) {}
    // };

    void Align_freq (
        Params params,
        // const std::vector<std::string>& reference,
        // const std::vector<std::string>& query,
        const std::vector<std::vector<float>>& freqRef,
        const std::vector<std::vector<float>>& freqQry,
        const std::vector<std::vector<float>>& gapOp,
        const std::vector<std::vector<float>>& gapEx,
        const std::vector<std::vector<float>>& gapCl,
        std::vector<int8_t>& aln,
        int16_t& errorType
        // size_t num_alignments
    );

    void Tile (
        // const std::vector<std::string>& reference,
        // const std::vector<std::string>& query,
        const std::vector<std::vector<float>>& reference, 
        const std::vector<std::vector<float>>& query, 
        const std::vector<std::vector<float>>& gapOp,
        const std::vector<std::vector<float>>& gapEx,
        const std::vector<std::vector<float>>& gapCl,
        Params params,
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

    int Score (
        Params params,
        const std::vector<int8_t> &aln,
        const std::string &reference,
        const std::string &query,
        const int ref_idx,
        const int query_idx
    );

    
}

#endif