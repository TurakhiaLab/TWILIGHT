#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <boost/program_options.hpp> 
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <cstdio>
#include <sys/stat.h>



#include "kseq.h"
#include "zlib.h"

#ifndef TREE_HPP
#include "tree.hpp"
#endif

#ifndef MSA_HPP
#include "msa.hpp"
#endif

#ifndef ALIGN_HPP
#include "align.cuh"
#endif

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#ifndef PARTITION_HPP
#include "treePartition.hpp"
#endif


namespace po = boost::program_options;
// po::options_description mainDesc("MSA Command Line Arguments");



// KSEQ_INIT2(, gzFile, gzread)
Params* setParameters(po::variables_map& vm);
void setOptions(po::variables_map& vm, msa::option* option);
// po::options_description mainDesc("MSA Command Line Arguments");

// void parseArguments(po::options_description mainDesc, int argc, char** argv);
void readSequences(po::variables_map& vm, msa::utility* util, Tree* tree);
void readSequences(std::string seqFileName, msa::utility* util, Tree* tree);
void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition);
void outputFinal (std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util, int& totalSeqs);

void printTree(Node* node, int grpID);
void printLeaves(Node* node);
Tree* readNewick(po::variables_map& vm);
Tree* readNewick(std::string treeFileName);
void readFreq(std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util);
void getFreq(Tree* tree, paritionInfo_t* partition, msa::utility* util);

void outputFile(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputSubtree(std::string fileName, Tree* T);
void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string> seqs);
// void readMSANUpdate(std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util, int& totalSeqs);


void createOverlapMSA(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param);
void transitivityMerge(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util);
void msaOnSubtree (Tree* T, msa::utility* util, msa::option* option, paritionInfo_t* partition, Params& param);
void alignSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param);
void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util);

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID);
void getPostOrderList(Node* node, std::stack<Node*>& msaStack);

__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result);
double getSPScore_cpu(std::vector<std::string>& alignment, Params& param);
double getSPScore_gpu(msa::utility* util, Params& param);

// Regressive method
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


#endif