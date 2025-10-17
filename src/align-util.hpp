#ifndef ALNUTIL_HPP
#define ALNUTIL_HPP

#include <fstream>
#include <tbb/parallel_invoke.h>

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


const int32_t PROFILE_LEN_TH = 1000;
const int32_t SEQNUM_TH = 100;
using gappyColumnQueue = std::pair<std::queue<std::tuple<int, int, float*>>, std::queue<std::tuple<int, int, float*>>>;
using gappyColumnList = std::vector<gappyColumnQueue>;
using alnPathList = std::vector<int8_t>;

void updateNode(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util);
void calculateProfileFreq(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t profileLen, int32_t profileSize, std::pair<int, int> startPos);
void calculatePSGOP(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t seqLen, std::pair<int32_t,int32_t> offset, std::pair<int32_t,int32_t>lens, Params& param);

// void removeGappyColumns(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<int,int>& rawLens);
void removeGappyColumns(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, gappyColumnQueue& gappyColumns, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<int,int>& rawLens);

// void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx);
void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, gappyColumnQueue& gappyColumns, std::pair<int, int>& debugIdx, float* hostParam, Params& param);

void mergeInsertions (Tree* tree, Node* nodeRef, std::vector<Node*>& nodes, msa::utility* util, std::vector<std::vector<int8_t>>& alnBad);
void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln);
void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx);

void fallback2cpu(std::vector<int>& fallbackPairs, Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option);
void getPostOrderList(Node* node, std::stack<Node*>& postStack);
void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID);
double calColumnSimilarity(Tree* tree, Node* node, msa::utility* util, Params& param);
void createAlnPairs(Tree* tree, msa::utility* util, msa::option* option, std::vector<std::pair<Node*, Node*>>& alnPairs);
void mergedAlignedSeqs(Tree* tree, msa::utility* util, msa::option* option, const std::vector<std::pair<Node*, Node*>>& alnPairs);
void removeEnds(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<alnPathList,alnPathList>& endsPaths);
void smith_waterman_DNA(const std::string &seq1, const std::string &seq2, int& max_r, int& max_q);

bool checkAlignedDescendant(Node* node, msa::utility* util);
void collectAlignedDescendant(Node* node, msa::utility* util, std::vector<std::string>& descendants);
void createAlnPairs(Tree* tree, msa::utility* util, msa::option* option, std::vector<std::pair<Node*, Node*>>& alnPairs);


void global_alignmentDNA(const std::string &seq1, const std::string &seq2, std::vector<int8_t>& alnPath);
void global_alignmentProtein(const std::string &seq1, const std::string &seq2, std::vector<int8_t>& alnPath);
std::string getConsensusDNA(float* profile, int len);
std::string getConsensusProtein(float* profile, int len);

#endif


