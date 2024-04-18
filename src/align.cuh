#ifndef ALIGN_HPP
#define ALIGN_HPP

#include <string>
#include <vector>
#include <thrust/device_vector.h>
#include <iostream>
#include <cmath>
#include <bits/stdc++.h>

struct Params 
{
    int16_t match;
    int16_t mismatch;
    int16_t gapOpen;
    int16_t gapExtend; //for gap-affine

    int16_t xdrop; //optional for now
    int16_t marker; //optional for now


    Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen) : //linear gap
        match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen) {}

    Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen, int16_t t_gapExtend) : //affine gap 
        match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen), gapExtend(t_gapExtend) {}

    Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen, int16_t t_gapExtend, int16_t t_xdrop) : //xdrop with affine gap 
        match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop){}

    Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen, int16_t t_gapExtend, int16_t t_xdrop, int16_t marker) : //TALCO-Xdrop with affine gap 
        match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop), marker(marker) {}
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
    int16_t* param
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

#endif