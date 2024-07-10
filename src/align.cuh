#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <string>
#include <vector>
#include <thrust/device_vector.h>
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <tbb/parallel_for.h>

const int FRONT_WAVE_LEN = 1024;
const int THREAD_NUM = 256;

typedef float paramType;

struct Params 
{
    paramType match;
    paramType mismatch;
    paramType gapOpen;
    paramType gapExtend; //for gap-affine
    paramType trans; // transition

    paramType xdrop; //optional for now
    paramType marker; //optional for now
    int scoreMode; //0: match/mismatch, 1: hoxd70

    // paramType hoxd70 [5][5] = { {  91, -114,  -31, -123, -100},
    //                             {-114,  100, -125,  -31, -100},
    //                             { -31, -125,  100, -114, -100},
    //                             {-123,  -31, -114,   91, -100},
    //                             {-100, -100, -100, -100, -100} }; 
    // paramType hoxd70_gapOpen = -400;
    // paramType hoxd70_gapExtend = -30;

    paramType hoxd70 [5][5] = { {  2.22,  -1.86,  -1.46,  -1.39 },  // a
                                { -1.86,   1.16,  -2.48,  -1.05 },  // c
                                { -1.46,  -2.48,   1.03,  -1.74 },  // g
                                { -1.39,  -1.05,  -1.74,   1.65 }}; // t
    paramType hoxd70_gapOpen = -20;
    paramType hoxd70_gapExtend = -2;


    // Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, int mode) : //linear gap
    //     match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), scoreMode(mode) {}

    // Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend, int mode) : //affine gap 
    //     match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend), scoreMode(mode) {}

    Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend, paramType t_xdrop, int mode) : //xdrop with affine gap 
        match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend), scoreMode(mode) {
            // if (mode == 0)      this->xdrop = round((FRONT_WAVE_LEN/3)*(-t_gapExtend));
            // // else if (mode == 1) this->xdrop = round((FRONT_WAVE_LEN/3)*(-hoxd70_gapExtend));
            // else if (mode == 1) this->xdrop = 1000000;
            this->xdrop = round((FRONT_WAVE_LEN/3)*(-t_gapExtend));
        }

    // Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend, paramType t_xdrop, paramType marker, int mode) : //TALCO-Xdrop with affine gap 
    //     match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop), marker(marker), scoreMode(mode) {}
};


// __global__ void alignGrpToGrp_talco
// (
//     char* ref,
//     char* query,
//     int16_t* param,
//     char* alignment,
//     int32_t* seqInfo
    
// );

__global__ void alignGrpToGrp_talco
(
    // char* seqs,
    uint16_t* freq,
    int8_t* aln,
    // int32_t* seqIdx,
    int32_t* len,
    int32_t* alnLen,
    int32_t* seqInfo,
    paramType* param
);

void alignGrpToGrp_traditional
(
    uint16_t* freq,
    int32_t seqLen,
    int32_t refLen,
    int32_t qryLen,
    Params& param,
    std::vector<int8_t>& aln
);


enum STATE 
{
    HH = 0,
    HI = 1,
    HD = 2,
    IH = 3,
    II = 4,
    DH = 5,
    DD = 6
};

void tracebackGrpToGrp
(
    int8_t state,
    std::vector<int8_t>& TB,
    std::vector<int32_t>& wfLL,
    std::vector<int32_t>& wfLen,
    std::vector<int8_t>& aln,
    int32_t ref,
    int32_t query
);

int8_t updateState(STATE& currentHState, STATE& currentIState, STATE& currentDState);


#endif