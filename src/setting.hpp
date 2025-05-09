#ifndef SETTING_HPP
#define SETTING_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include <boost/program_options.hpp> 
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

const std::unordered_set<char> protein_only = {'E', 'F', 'I', 'J', 'L', 'P', 'Q', 'Z'};
const std::unordered_set<char> nucleotide_only = {'U'};
const std::map<char, int> PROTEIN = {{'A',0}, {'C',1}, {'D',2}, {'E',3}, {'F',4}, {'G',5}, {'H',6}, {'I',7}, {'K',8}, {'L',9}, {'M',10}, {'N',11}, {'P',12}, {'Q',13}, {'R',14}, {'S',15}, {'T',16}, {'V',17}, {'W',18}, {'Y',19}, {'-',21}}; // 20 for all other characters (ambiguous)
const std::map<char, int> NUCLEOTIDE = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'U', 3}, {'-', 5}}; // 4 for all other characters (ambiguous)

const float BLOSUM62[20][20] = {
//    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    { 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2 },
    { 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2 },
    {-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3 },
    {-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2 },
    {-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3 },
    { 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3 },
    {-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2 },
    {-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1 },
    {-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2 },
    {-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1 },
    {-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1 },
    {-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2 },
    {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3 },
    {-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1 },
    {-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2 },
    { 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2 },
    { 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2 },
    { 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1 },
    {-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2 },
    {-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7 }};

struct Params 
{
    float gapOpen;
    float gapExtend; //for gap-affine
    float gapBoundary; // gap penalty at ends
    float xdrop; //optional for now
    float** scoringMatrix;
    int matrixSize;
    Params(po::variables_map& vm, char type);
    ~Params();
};

namespace msa
{
    struct option
    {
        int gpuNum;
        int cpuNum;
        int maxSubtree;
        int maxSubSubtree;
        int gappyHorizon;
        int alnMode; //0: MSA from raw sequences, 1: merge multiple MSA files
        float gappyVertical;
        float lenDev;
        float maxAmbig;
        std::vector<int> gpuIdx;
        int nowProcess;
        bool debug;
        bool cpuOnly;
        bool noFilter;
        char type; // 'n' for dna/rna, 'p' for protein
        // == Position-specific option control ==
        bool psgop;
        bool psgopAuto;
        bool calSim; // calculate column similarity
        bool redo;
        bool alignGappy;
        float divTH = 0.1;
        // ======================================
        std::string treeFile;
        std::string seqFile;
        std::string outFile;
        std::string msaDir;
        std::string tempDir;
        std::string subtreeDir;
        std::string merger;
        std::string outType;
        bool printDetail;
        bool deleteTemp;
        option(po::variables_map& vm);
    };

    struct utility
    {
        std::unordered_map<std::string, std::string> rawSeqs;
        std::map<int, std::vector<std::string>> badSequences;
        // std::map<int, std::vector<int>> subtrees;
        std::unordered_map<std::string, int>  seqsIdx;
        std::unordered_map<int, std::string> seqsName;
        std::unordered_map<std::string, int>  seqsLen;
        std::unordered_map<int, bool> seqsStorage;
        std::map<int, std::vector<std::vector<float>>> profileFreq;
        std::unordered_map<std::string, std::string> seqsCIGAR;

        std::unordered_map<int, bool> lowQuality;
        
        char** alnStorage[2] = {nullptr, nullptr};
        int nowProcess = 0;
        int subtreeIdx = 0;
        int memLen = 0;
        int* seqMemLen = 0;
        int memNum = 0;
        int seqLen = 0;
        int maxRawSeqLen = 0;
        const float timesBigger = 2.0;
        void changeStorage(int idx); 
        void setSubtreeIdx(int idx); 
        void seqsFree();
        void seqFree(int i);
        void seqMalloc(int seqNum, int seqLen, option* option);
        void seqMalloc(int seqLen, int idx);
        void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap, option* option);
        void storeCIGAR();
        void memCheck(int seqLen, option* option);
        void memCheck(int seqLen, int idx);
        void clearAll();
        void debug(int& debugNum);
    };
} 

#endif