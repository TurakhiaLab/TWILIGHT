#ifndef MSAGPU_HPP
#define MSAGPU_HPP

#ifndef MSA_HPP
#include "msa-cpu.hpp"
#endif

#ifndef ALIGN_HPP
#include "align.cuh"
#endif

void getGpuInfo (po::variables_map& vm, msa::option* option);
void msaGpu(  Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaGpu_s(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaOnSubtreeGpu (Tree* T, msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param);

#endif