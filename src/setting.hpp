#ifndef SETTING_HPP
#define SETTING_HPP

#include <string>
#include <unordered_map>
#include <map>
#include <vector>

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
    };

    struct utility
    {
        std::unordered_map<std::string, std::string> rawSeqs;
        std::unordered_map<std::string, std::string> alnSeqs;
        std::map<int, std::vector<std::string>> badSequences;
        std::map<int, std::vector<int>> subtrees;
        std::unordered_map<std::string, int>  seqsIdx;
        std::unordered_map<int, std::string> seqsName;
        std::unordered_map<std::string, int>  seqsLen;
        std::unordered_map<int, bool> seqsStorage;
        std::map<int, std::vector<std::vector<float>>> profileFreq;
        
        char** alnStorage[2] = {nullptr, nullptr};
        int nowProcess = 0;
        size_t subtreeIdx = 0;
        size_t memLen;
        size_t memNum;
        size_t seqLen;
        size_t seqNum;
        int nowStore = 0;
        const float timesBigger = 2.0;
        void changeStorage(int idx); 
        void setSubtreeIdx(int idx); 
        void seqsFree();
        void seqFree(int i);
        void seqMalloc(int seqNum, int seqLen);
        void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap);
        void memCheck(int seqLen);
        void clearAll();
        void debug();
    };
} 

#endif