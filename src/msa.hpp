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
            return;
        }
        void seqFree(){
            free(seqBuf[0]);
            free(seqBuf[1]);
            return;
        }
        void seqMalloc(int seqNum, int seqLen){
            if (seqBuf[0] != nullptr && seqBuf[1] != nullptr) {
                char* temp[2] = {nullptr, nullptr};
                const float timesBigger = 1.5;
                int adjustLen = static_cast<int>(seqLen * timesBigger); 
                temp[0] = (char*)malloc(seqNum * adjustLen * sizeof(char));
                temp[1] = (char*)malloc(seqNum * adjustLen * sizeof(char));
                for (int j = 0; j < seqNum; ++j) {
                    for (int i = 0; i < adjustLen; ++i) {
                        if (i < memLen) {
                            temp[0][j*adjustLen+i] = seqBuf[0][j*memLen+i];
                            temp[1][j*adjustLen+i] = seqBuf[1][j*memLen+i];
                        }
                        else {
                            temp[0][j*adjustLen+i] = 0;
                            temp[1][j*adjustLen+i] = 0;
                        }
                    }
                }
                seqFree();
                seqBuf[0] = temp[0];
                seqBuf[1] = temp[1];
                memLen = adjustLen;
            }
            else {
                const float timesBigger = 1.5;
                int adjustLen = static_cast<int>(seqLen * timesBigger); 
                seqBuf[0] = (char*)malloc(seqNum * adjustLen * sizeof(char));
                seqBuf[1] = (char*)malloc(seqNum * adjustLen * sizeof(char));
                memNum = seqNum;
                memLen = adjustLen;
            }    
            printf("Allocate memory... Size = %d x %d\n", memNum, memLen);
            return;
        }
        void memCheck(int seqLen) {
            if (seqLen > memLen) {
                std::cout << "Reallocate Memory...\n";
                seqMalloc(memNum, memLen);
            }
            return;
        }
    };
} 

#endif