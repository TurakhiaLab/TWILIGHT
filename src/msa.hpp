#ifndef MSA_HPP
#define MSA_HPP

#include <string>
#include <unordered_map>

namespace msa
{
    struct utility
    {
        std::unordered_map<std::string, std::string> seqs;
        std::unordered_map<std::string, int>  seqsIdx;
        std::unordered_map<std::string, int>  seqsLen;
        std::unordered_map<int, bool> seqsStorage;
        char* seqBuf[2] = {nullptr, nullptr};
        // char* seqBuf0 = nullptr;
        // char* seqBuf1 = nullptr;
        int memLen;
        int memNum;
        int nowStore = 0;
        void changeStorage(int idx) {
            seqsStorage[idx] = (seqsStorage[idx] == 0) ? 1 : 0;
        }
        void seqFree(){
            free(seqBuf[0]);
            free(seqBuf[1]);
        }
        void seqMalloc(int seqNum, int seqLen){
            if (seqBuf[0] != nullptr && seqBuf[1] != nullptr)
                seqFree();
            const float timesBigger = 1.5;
            int adjustLen = static_cast<int>(seqLen * timesBigger); 
            seqBuf[0] = (char*)malloc(seqNum * adjustLen * sizeof(char));
            seqBuf[1] = (char*)malloc(seqNum * adjustLen * sizeof(char));
            memNum = seqNum;
            memLen = adjustLen;
            printf("Num: %d, Len: %d\n", memNum, memLen);
        }
        void memCheck(int seqLen) {
            if (seqLen > (0.9*memLen)) {
                seqMalloc(memNum, memLen);
            }
        }
    };
} 

#endif