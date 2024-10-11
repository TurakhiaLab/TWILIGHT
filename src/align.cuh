#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include <tbb/parallel_for.h>


#ifndef SETTING_HPP
#include "setting.cuh"
#endif

const int FRONT_WAVE_LEN = 1024+512;
const int THREAD_NUM = 256;

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
    float* param
);
#endif