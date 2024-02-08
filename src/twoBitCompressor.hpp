#ifndef TWOBITCOMPRESSOR_HPP
#define TWOBITCOMPRESSOR_HPP

#include <string>
#include <iostream>
#include <tbb/parallel_for.h>
#include "kseq.h"

void twoBitCompressor(std::string seq, size_t seqLen, uint32_t* compressedSeq);

#endif