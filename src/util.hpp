#ifndef UTIL_HPP
#define UTIL_HPP

#include "kseq.h"
#include <zlib.h>

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

// print Tree (Only for debugging)
void printTree(Node* node, int grpID);

// read
void readSequences(msa::utility* util, msa::option* option, Tree* tree);
Tree* readNewick(std::string treeFileName);
void readFrequency(msa::utility* util, msa::option* option);
void readMSA_and_Seqs(msa::utility* util, msa::option* option, Tree* tree);
void readBackboneAln(msa::utility* util, msa::option* option, Tree* tree);
void readNewSequences(msa::utility* util, msa::option* option, Tree* tree);

// output
void outputAln(msa::utility* util, msa::option* option, Tree* T);
void outputSubAln(msa::utility* util, msa::option* option, Tree* T, int subtreeIdx);
void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputSubtreeTrees(Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option);
void outputSubtree(Tree* tree, msa::option* option, int subtreeIdx);
void outputPrunedTree(Tree* T, msa::option* option);
void outputFinal (Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs);
void outputSubtreeSeqs(std::string fileName, std::vector<std::pair<std::string, std::string>>& seqs, bool compressed);

// auxiliary
bool cmp(std::string a, std::string b);
bool cmp2(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b);
void getSubtreeNewick(Node* root, std::string& outputString);
double calSPScore(std::string alnFile, msa::utility* util, Params* param);
void storeFreq(msa::utility* util, msa::option* option, Tree* T, int grpID);
void updateSeqLen(Tree* tree, partitionInfo_t* partition, msa::utility* util);
Tree* pruneTree(Tree* T, std::string& seqFileName);

#endif