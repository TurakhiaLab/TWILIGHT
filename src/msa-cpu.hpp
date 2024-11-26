#ifndef MSA_HPP
#define MSA_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <boost/program_options.hpp> 
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/mutex.h>
#include<tbb/spin_rw_mutex.h>
#include <cstdio>
#include <sys/stat.h>
#include <queue>

#ifndef TREE_HPP
#include "tree.hpp"
#endif

#ifndef SETTING_HPP
#include "setting.hpp"
#endif

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#ifndef PARTITION_HPP
#include "treePartition.hpp"
#endif

void msaOnSubtree  (Tree* T,             msa::utility* util, msa::option* option, partitionInfo_t* partition, Params& param);
void alignSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option,                             Params& param);
void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option,                             Params& param);

void createOverlapAlnCpu(Tree* tree,                std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaCpu(             Tree* tree,                std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void transitivityMerge(  Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option);

void updateNode(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util);
void calculateProfileFreq(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, int32_t seqLen, int32_t offsetf, int32_t offsetg);
void calculatePSGOP(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t seqLen, int32_t offsetf, int32_t offsetg, Params& param);

void removeGappyColumns(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, int32_t& newRef, int32_t& newQry, int32_t seqLen, int32_t offsetf, int32_t offsetg);
void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx, msa::option* option);

void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln);
void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx);

void fallback2cpu(std::vector<int>& fallbackPairs, Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option);
void getPostOrderList(Node* node, std::stack<Node*>& postStack);
void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID);
double calColumnSimilarity(Tree* tree, Node* node, msa::utility* util, Params& param);

#endif