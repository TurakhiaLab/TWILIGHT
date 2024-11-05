#ifndef UTIL_HPP
#define UTIL_HPP

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
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#include "kseq.h"
#include "zlib.h"

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

void readSequences(msa::utility* util, msa::option* option, Tree* tree);
void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option);
Tree* readNewick(std::string treeFileName);
void readFrequency(msa::utility* util, msa::option* option, Tree* tree);

void printTree(Node* node, int grpID);
void printLeaves(Node* node);

void outputAln(msa::utility* util, msa::option* option, Tree* T);
void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputSubtree(Tree* tree, msa::option* option, int subtreeIdx);
void outputFinal (Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs);
void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string>& seqs);
void outputSubtreeCIGAR(std::string fileName, std::map<std::string, std::string>& seqs);

void storeFreq(msa::utility* util, Tree* T, int grpID);
void updateSeqLen(Tree* tree, partitionInfo_t* partition, msa::utility* util);
// auxiliary
bool cmp(std::string a, std::string b);
void getSubtreeNewick(Node* root, std::string& outputString);
double calSPScore(std::string alnFile, msa::utility* util, Params* param);

#endif