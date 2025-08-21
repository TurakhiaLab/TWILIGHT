#ifndef MSAGPU_HPP
#define MSAGPU_HPP

#ifndef MSA_HPP
#include "../msa-cpu.hpp"
#endif

#ifndef ALIGN_HPP
#include "align.hip.hpp"
#endif


const int MAX_PROFILE_LEN = (1 << 16);
const int MAX_PROFILE_LEN_P = (1 << 14);

void getGpuInfo (po::variables_map& vm, msa::option* option);
void msaGpu(  Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaOnSubtreeGpu (Tree* T, msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param);
void placementGpu(Tree *T, msa::utility *util, msa::option *option, partitionInfo_t *partition, Params &param);
bool comp (Node* A, Node* B);

#endif
