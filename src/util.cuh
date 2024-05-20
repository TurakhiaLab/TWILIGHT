#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>
#include <boost/program_options.hpp> 
#include <tbb/parallel_for.h>


#include "../src/kseq.h"
#include "zlib.h"

#ifndef TREE_HPP
#include "../src/tree.hpp"
#endif

#ifndef MSA_HPP
#include "../src/msa.hpp"
#endif

#ifndef ALIGN_HPP
#include "../src/align.cuh"
#endif

#ifndef ALIGN_TALCO_HPP
#include "../src/TALCO-XDrop.hpp"
#endif

namespace po = boost::program_options;

// KSEQ_INIT2(, gzFile, gzread)

// po::options_description mainDesc("MSA Command Line Arguments");

// void parseArguments(int argc, char** argv);
void printTree(Node* node);
void printLeaves(Node* node);
Tree* readNewick(po::variables_map& vm);
// void readSequences(po::variables_map& vm, msa::utility* util);
void checkAlignment(std::vector<std::string>& ref);
void msaPostOrderTraversal_gpu_org(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void createOverlapMSA(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);

void msaPostOrderTraversal_gpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void msaPostOrderTraversal_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void transitivityMerge_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util);
void transitivityMerge_cpu_mod(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util);
void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& hier, std::stack<Node*> msaStack, int grpID, int mode);
void getPostOrderList(Node* node, std::stack<Node*>& msaStack);
__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result);
double getSPScore_cpu(std::vector<std::string>& alignment, Params& param);
double getSPScore_gpu(std::vector<std::string>& alignment, msa::utility* util, Params& param);


#endif