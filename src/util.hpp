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

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#ifndef PARTITION_HPP
#include "treePartition.hpp"
#endif

void readSequences(po::variables_map& vm, msa::utility* util, msa::option* option, Tree* tree);
void readSequences(std::string seqFileName, msa::utility* util, msa::option* option, Tree* tree);
void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option);
Tree* readNewick(po::variables_map& vm);
Tree* readNewick(std::string treeFileName);
void readFreq(std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util);

void printTree(Node* node, int grpID);
void printLeaves(Node* node);

void outputAln(std::string fileName, msa::utility* util, msa::option* option, Tree* T);
void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID);
void outputSubtree(Tree* tree, msa::option* option, int subtreeIdx);
void outputFinal (po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs);
void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string>& seqs);
void outputSubtreeCIGAR(std::string fileName, std::map<std::string, std::string>& seqs);

// auxiliary
bool cmp(std::string a, std::string b);
void getSubtreeNewick(Node* root, std::string& outputString);
double calSPScore(std::string alnFile, msa::utility* util, Params* param);

#endif