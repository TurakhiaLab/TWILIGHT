#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <bits/stdc++.h>

const int _FRONT_WAVE_LEN = 1350;
const int _THREAD_NUM = 256;
const int _BLOCKSIZE = 2048;
const int _MARKER = 200;

// General Usage
__global__ void alignGrpToGrp_freq
(
    float* freq,
    int8_t* aln,
    int32_t* len,
    int32_t* num,
    int32_t* alnLen,
    int32_t* seqInfo,
    float* gapOpen,
    float* gapExtend,
    float* param
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
