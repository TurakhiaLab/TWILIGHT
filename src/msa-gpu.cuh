#ifndef MSAGPU_HPP
#define MSAGPU_HPP

#ifndef MSA_HPP
#include "msa-cpu.hpp"
#endif

#ifndef ALIGN_HPP
#include "align.cuh"
#endif

void createOverlapAlnGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaGpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaOnSubtreeGpu (Tree* T, msa::utility* util, msa::option* option, paritionInfo_t* partition, Params& param);
void alignSubtreesGpu (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param);

#endif