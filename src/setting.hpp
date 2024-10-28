#ifndef SETTING_HPP
#define SETTING_HPP

#include <iostream>
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
#include <tbb/mutex.h>


namespace po = boost::program_options;

typedef float paramType;

struct Params 
{
    float gapOpen;
    float gapExtend; //for gap-affine
    float gapClose;
    float xdrop; //optional for now
    int userDefine; //0: simple match/mismatch, 1: kimura, 2: userdefined
    // hoxd70
    // float userMatrix [5][5] = { {  91, -114,  -31, -123, -100},
    //                                 {-114,  100, -125,  -31, -100},
    //                                 { -31, -125,  100, -114, -100},
    //                                 {-123,  -31, -114,   91, -100},
    //                                 {-100, -100, -100, -100, -100} }; 
    // float userGapOpen = -400;
    // float userGapExtend = -30;

    float scoringMatrix [5][5];
    float userMatrix [5][5] = { {  2.22,  -1.86,  -1.46,  -1.39 },  // a
                                { -1.86,   1.16,  -2.48,  -1.05 },  // c
                                { -1.46,  -2.48,   1.03,  -1.74 },  // g
                                { -1.39,  -1.05,  -1.74,   1.65 }}; // t
    float userGapOpen = -10;
    float userGapClose = -10;
    float userGapExtend = -2;
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
        float gappyVertical;
        std::vector<int> gpuIdx;
        int nowProcess;
        bool debug;
        bool cpuOnly;
        bool psgop;
        std::string outType;
        std::string tempDir;
        std::string merger;
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
        
        char** alnStorage[2] = {nullptr, nullptr};
        int nowProcess = 0;
        size_t subtreeIdx = 0;
        size_t memLen;
        size_t* seqMemLen;
        size_t memNum;
        size_t seqLen;
        size_t seqNum;
        int nowStore = 0;
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