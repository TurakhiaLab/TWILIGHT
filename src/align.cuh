#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <string>
#include <vector>
#include <thrust/device_vector.h>
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <tbb/parallel_for.h>

const int FRONT_WAVE_LEN = 512;

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
    


    Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen) : //linear gap
        match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen) {}

    Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend) : //affine gap 
        match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend) {}

    Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend, paramType t_xdrop) : //xdrop with affine gap 
        match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop){}

    Params(paramType t_match, paramType t_mismatch, paramType t_trans, paramType t_gapOpen, paramType t_gapExtend, paramType t_xdrop, paramType marker) : //TALCO-Xdrop with affine gap 
        match(t_match), mismatch(t_mismatch), trans(t_trans), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop), marker(marker) {}
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