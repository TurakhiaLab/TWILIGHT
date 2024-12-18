#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <bits/stdc++.h>

<<<<<<< HEAD
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
=======
const int FRONT_WAVE_LEN = 1024+512;
const int THREAD_NUM = 256;
const int BLOCKSIZE = 2048;

// General Usage
__global__ void alignGrpToGrp_freq
>>>>>>> 3eea4f208434ec5fef8f9e150c042140efb8ba22
(
    float* freq,
    int8_t* aln,
    int32_t* len,
    int32_t* num,
    int32_t* alnLen,
    int32_t* seqInfo,
<<<<<<< HEAD
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
=======
    float* gapOpen,
    float* gapExtend,
    float* param
>>>>>>> 3eea4f208434ec5fef8f9e150c042140efb8ba22
);

// Specify for small profiles (< 24 sequences)
// The only difference with alignGrpToGrp_freq 
// is that the sequence is copied instead of 
// the frequency and the frequency is calculated on the device
__global__ void alignGrpToGrp_seq
(
    char* seqs,
    int8_t* aln,
    int32_t* len,
    int32_t* num,
    int32_t* alnLen,
    int32_t* seqInfo,
    float* gapOpen,
    float* gapExtend,
    float* weight,
    float* param
);
#endif