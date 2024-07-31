#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <boost/program_options.hpp> 
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>


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
// void outputFile(std::string fileName, msa::utility* util, Tree* T, int grpID);
void printTree(Node* node, int grpID);
void printLeaves(Node* node);
Tree* readNewick(po::variables_map& vm);

// void getLongestDescendent(Tree* tree, msa::utility* util);
// void getSubMSANodes(std::vector<Node*>& subMSANodes, Node* startNode, int N);
// void getSubMSAs(std::map<int, std::vector<std::vector<Node*>>>& subMSAs, Tree* T, int N);
// void getAlnPairs(std::vector<std::vector<std::pair<Node*, Node*>>>& alnPairs, std::vector<std::vector<Node*>>& clusters);
// void storeMSA(Tree* T, std::vector<Node*>& nodes, msa::utility* util, int level);
// void resetSeqMem(std::vector<Node*>& nodes, msa::utility* util);
// void merger_ref(Tree* tree, std::map<int, Node*>& refNodes, msa::utility* util, std::string& refString, std::string& qryString, int qLevel);
// void merger_qry(Tree* tree, std::vector<Node*>& qryNodes, msa::utility* util, std::string& refString, std::string& qryString, int qLevel);
// void transitivityMerge_regressive(Tree* tree, std::map<int, std::vector<std::vector<Node*>>>& subMSAs, msa::utility* util);
// void msaPostOrderTraversal_multigpu_regressive(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);

void createOverlapMSA(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param);
void transitivityMerge_cpu_mod(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util);

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& hier, std::stack<Node*> msaStack, int grpID, int mode);
void getPostOrderList(Node* node, std::stack<Node*>& msaStack);

__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result);
double getSPScore_cpu(std::vector<std::string>& alignment, Params& param);
double getSPScore_gpu(std::vector<std::string>& alignment, msa::utility* util, Params& param);


#endif