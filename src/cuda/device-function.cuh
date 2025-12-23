#ifndef DEVICE_HPP
#define DEVICE_HPP

#include <bits/stdc++.h>



namespace device_function {
    constexpr int _FRONT_WAVE_LEN = 1350;
    constexpr int _THREAD_NUM = 256;
    constexpr int _BLOCKSIZE = 2048;
    constexpr int _MARKER = 200;
    
    __global__ void parallelProfileAlignment (
        float* freqPointers [],
        int8_t* alnPointers [],
        int32_t* len,
        int32_t* num,
        int32_t* alnLen,
        int32_t* seqInfo,
        float* gapOpenPointers [],
        float* gapExtendPointers [],
        float* param
    );

    __global__ void parallelProfileAlignment_Fast(
        float* freq, 
        int8_t *aln, 
        int32_t* len, 
        int32_t* num, 
        int32_t* alnLen, 
        int32_t *seqInfo,  
        float* gapOpen, 
        float* gapExtend, 
        float* param
    );
}


#endif
