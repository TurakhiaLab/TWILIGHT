#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <boost/program_options.hpp> 
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/mutex.h>
#include <cstdio>
#include <sys/stat.h>
#include <queue>



#include "kseq.h"
#include "zlib.h"

#ifndef TREE_HPP
#include "tree.hpp"
#endif

#ifndef SETTING_HPP
#include "setting.hpp"
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

void readSequences(po::variables_map& vm, msa::utility* util, msa::option* option, Tree* tree);
void readSequences(std::string seqFileName, msa::utility* util, msa::option* option, Tree* tree);
void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option);
Tree* readNewick(po::variables_map& vm);
Tree* readNewick(std::string treeFileName);
void readFreq(std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util);

void printTree(Node* node, int grpID);
void printLeaves(Node* node);

void outputAln(std::string fileName, msa::utility* util, msa::option* option, Tree* T, int grpID);
void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputSubtree(Tree* tree, msa::option* option, int subtreeIdx);
// void outputFinal (std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs);
void outputFinal (po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs);
void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string>& seqs);
void outputSubtreeCIGAR(std::string fileName, std::map<std::string, std::string>& seqs);

// auxiliary
bool cmp(std::string a, std::string b);
void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID);
void getPostOrderList(Node* node, std::stack<Node*>& msaStack);
void getSubtreeNewick(Node* root, std::string& outputString);


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