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

#include <boost/program_options.hpp> 
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct Params 
{
    float gapOpen;
    float gapExtend; //for gap-affine
    float gapBoundary; // gap penalty at ends
    float xdrop; //optional for now
    float scoringMatrix [5][5];
    Params(po::variables_map& vm);
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
        std::vector<int> gpuIdx;
        int nowProcess;
        bool debug;
        bool cpuOnly;
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
        void debug();
    };
} 

#endif