#ifndef MSA_HPP
#define MSA_HPP

#ifndef ALNUTIL_HPP
#include "align-util.hpp"
#endif

void msaOnSubtree  (Tree* T,             msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param);
void placement(     Tree *T,             msa::utility *util, msa::option *option, partitionInfo_t *partition, Params &param);
void alignSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option,                             Params& param);
void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option,                             Params& param);

void createOverlapAlnCpu(Tree* tree,                std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaCpu(             Tree* tree,                std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void transitivityMerge(  Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option);

#endif