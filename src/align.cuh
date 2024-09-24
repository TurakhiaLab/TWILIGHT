#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <tbb/parallel_for.h>

const int FRONT_WAVE_LEN = 1024+512;
const int THREAD_NUM = 256;

typedef float paramType;

struct Params 
{
    paramType gapOpen;
    paramType gapExtend; //for gap-affine
    paramType gapClose;
    paramType xdrop; //optional for now
    // paramType marker; //optional for now
    int userDefine; //0: simple match/mismatch, 1: kimura, 2: userdefined
    // hoxd70
    // paramType userMatrix [5][5] = { {  91, -114,  -31, -123, -100},
    //                                 {-114,  100, -125,  -31, -100},
    //                                 { -31, -125,  100, -114, -100},
    //                                 {-123,  -31, -114,   91, -100},
    //                                 {-100, -100, -100, -100, -100} }; 
    // paramType userGapOpen = -400;
    // paramType userGapExtend = -30;

    paramType scoringMatrix [5][5];

    paramType userMatrix [5][5] = { {  2.22,  -1.86,  -1.46,  -1.39 },  // a
                                    { -1.86,   1.16,  -2.48,  -1.05 },  // c
                                    { -1.46,  -2.48,   1.03,  -1.74 },  // g
                                    { -1.39,  -1.05,  -1.74,   1.65 }}; // t
    paramType userGapOpen = -10;
    paramType userGapClose = -10;
    paramType userGapExtend = -2;

    Params() {};
};

__global__ void alignGrpToGrp_talco
(
    // uint16_t* freq,
    // int32_t* freq,
    float* freq,
    int8_t* aln,
    // int32_t* seqIdx,
    int32_t* len,
    int32_t* num,
    int32_t* alnLen,
    int32_t* seqInfo,
    float* gapOpen,
    float* gapCont,
    float* gapEnd,
    paramType* param
);

/*
void alignGrpToGrp_traditional
(
    // uint16_t* freq,
    // int32_t* freq,
    // int32_t seqLen,
    const std::vector<std::vector<int>>& reference,
    const std::vector<std::vector<int>>& query, 
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
*/

#endif