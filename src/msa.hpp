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
        char** seqBuf[2] = {nullptr, nullptr};
        int memLen;
        int memNum;
        int nowStore = 0;
        const float timesBigger = 1.3;
        void changeStorage(int idx) {
            seqsStorage[idx] = (seqsStorage[idx] == 0) ? 1 : 0;
            return;
        }
        void seqFree(){
            for (int i = 0; i < this->memNum; ++i) {
                delete [] seqBuf[0][i];
                delete [] seqBuf[1][i];
            }
            delete [] seqBuf[0];
            delete [] seqBuf[1];
            return;
        }
        void seqMalloc(int seqNum, int seqLen){
            if (seqBuf[0] != nullptr && seqBuf[1] != nullptr) {
                char** temp[2] = {nullptr, nullptr};
                int adjustLen = static_cast<int>(seqLen * this->timesBigger);
                temp[0] = new char*[seqNum]; 
                temp[1] = new char*[seqNum]; 
                for (int i = 0; i < seqNum; ++i) {
                    temp[0][i] = new char[adjustLen];
                    temp[1][i] = new char[adjustLen];
                }
                for (int j = 0; j < seqNum; ++j) {
                    for (int i = 0; i < adjustLen; ++i) {
                        if (i < memLen) {
                            temp[0][j][i] = seqBuf[0][j][i];
                            temp[1][j][i] = seqBuf[1][j][i];
                        }
                        else {
                            temp[0][j][i] = 0;
                            temp[1][j][i] = 0;
                        }
                    }
                }
                seqFree();
                seqBuf[0] = temp[0];
                seqBuf[1] = temp[1];
                memLen = adjustLen;
            }
            else {
                int adjustLen = static_cast<int>(seqLen * this->timesBigger); 
                this->seqBuf[0] = new char*[seqNum]; 
                this->seqBuf[1] = new char*[seqNum]; 
                for (int i = 0; i < seqNum; ++i) {
                    this->seqBuf[0][i] = new char[adjustLen];
                    this->seqBuf[1][i] = new char[adjustLen];
                }
                for (int j = 0; j < seqNum; ++j) {
                    for (int i = 0; i < adjustLen; ++i) {
                        this->seqBuf[0][j][i] = 0;
                        this->seqBuf[1][j][i] = 0;
                    }
                }
                this->memNum = seqNum;
                this->memLen = adjustLen;
            }    
            printf("Allocate memory... Size = %d x %d\n", memNum, memLen);
            return;
        }
        
        void memCheck(int seqLen) {
            if (seqLen > this->memLen*0.9) {
                printf("Reallocate Memory. SeqLen (%d) > MemLen (%d)\n", seqLen, this->memLen);
                seqMalloc(this->memNum, seqLen);
            }
            return;
        }
    };
} 

#endif