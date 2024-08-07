#ifndef MSA_HPP
#define MSA_HPP

#include <string>
#include <unordered_map>
#include <map>

namespace msa
{
    struct utility
    {
        // std::unordered_map<std::string, std::string> seqs;
        std::map<int, std::vector<int>> subtrees;
        std::unordered_map<std::string, int>  seqsIdx;
        std::unordered_map<std::string, int>  seqsLen;
        std::unordered_map<int, bool> seqsStorage;
        std::map<int, std::vector<std::vector<uint16_t>>> profileFreq;
        char** alnStorage[2] = {nullptr, nullptr};
        char** seqs = nullptr;
        int nowProcess = 0;
        size_t subtreeIdx = 0;
        size_t gpuNum;
        size_t seqLen;
        size_t seqNum;
        size_t memLen;
        size_t memNum;
        int nowStore = 0;
        const float timesBigger = 1.3;
        void changeStorage(int idx) {
            seqsStorage[idx] = (seqsStorage[idx] == 0) ? 1 : 0;
            return;
        }
        void setSubtreeIdx(int idx) {
            this->subtreeIdx = idx;
            return;
        }
        void seqsFree(){
            for (int i = 0; i < this->seqNum; ++i) {
                delete [] this->seqs[i];
            }
            delete [] this->seqs;
            this->seqs = nullptr;
            return;
        }
        void reset() {
            seqsFree();
            subtrees.clear();
            seqsIdx.clear();
            seqsLen.clear();
            seqsStorage.clear();
            seqLen = 0;
            seqNum = 0;
            memLen = 0;
            memNum = 0;
            
        }

        void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap){
            this->seqNum = seqsMap.size();
            std::cout << "seqsMap size: " << seqsMap.size() << '\n';
            this->seqs = new char*[this->seqNum]; 
            this->seqLen = seqLen;
            
            int s = 0;
            for (auto seq: seqsMap) {
                this->seqs[s] = new char[seqLen];
                std::string seqName = seq.first;
                int subtreeIdx = seq.second.second;
                if (this->subtrees.find(subtreeIdx) == this->subtrees.end()) {
                    std::vector<int> temp;
                    this->subtrees[subtreeIdx] = temp;
                }
                this->subtrees[subtreeIdx].push_back(s);
                for (int j = 0; j < seqLen; ++j) {
                    this->seqs[s][j] =  (j < seq.second.first.size()) ? seq.second.first[j] : 0;
                }
                this->seqsIdx[seqName] = s;
                this->seqsLen[seqName] = seq.second.first.size();
                ++s;
                seqsMap.erase(seqName);
            }
            assert(s == seqNum);
            this->seqNum = seqNum;
            this->seqLen = seqLen;
            printf("Allocate memory for sequence storage, size = %lu x %lu\n", this->seqNum, seqLen);
            return;
        }
        void updateSeqs() {
            size_t len = max(this->memLen, this->seqLen);
            for (int i = 0; i < this->seqNum; ++i) {
                char* temp = new char[len];
                if (std::find(this->subtrees[this->subtreeIdx].begin(), this->subtrees[this->subtreeIdx].end(), i) != this->subtrees[this->subtreeIdx].end()) {
                    int s = this->seqsStorage[i];
                    for (int j = 0; j < len; ++j) {
                        temp[j] = (j < this->memLen) ? this->alnStorage[s][i][j] : 0;
                    }
                }
                else {
                    for (int j = 0; j < len; ++j) {
                        temp[j] = (j < this->seqLen) ? this->seqs[i][j] : 0;
                    }
                    delete [] this->seqs[i];
                    this->seqs[i] = nullptr;
                }
                
                this->seqs[i] = temp;
            }
            this->alnFree();
            this->seqLen = len;
        }
        void alnFree(){
            for (auto i: this->subtrees[this->subtreeIdx]) {
                delete [] this->alnStorage[0][i];
                delete [] this->alnStorage[1][i];
            }
            delete [] this->alnStorage[0];
            delete [] this->alnStorage[1];
            this->alnStorage[0] = nullptr;
            this->alnStorage[1] = nullptr;
            this->memLen = 0;
            this->memNum = 0;
            return;
        }
        void alnMalloc(size_t alnLen){
            int adjustLen = static_cast<int>(alnLen * this->timesBigger);
            if (alnStorage[0] != nullptr && alnStorage[1] != nullptr) {
                char** temp[2] = {nullptr, nullptr};
                temp[0] = new char*[this->seqNum]; 
                temp[1] = new char*[this->seqNum];
                for (int i = 0; i < this->seqNum; ++i) {
                    temp[0][i] = nullptr;
                    temp[1][i] = nullptr;
                }
                for (auto i: this->subtrees[this->subtreeIdx]) {
                    temp[0][i] = new char[adjustLen];
                    temp[1][i] = new char[adjustLen];
                    for (int j = 0; j < adjustLen; ++j) {
                        temp[0][i][j] = (j < memLen) ? this->alnStorage[0][i][j] : 0;
                        temp[1][i][j] = (j < memLen) ? this->alnStorage[1][i][j] : 0;
                    }
                    delete [] this->alnStorage[0][i];
                    delete [] this->alnStorage[1][i];
                    
                }
                delete [] this->alnStorage[0];
                delete [] this->alnStorage[1];
                this->alnStorage[0] = temp[0];
                this->alnStorage[1] = temp[1];
            }
            else {
                this->memNum = this->subtrees[this->subtreeIdx].size();
                this->alnStorage[0] = new char*[this->seqNum]; 
                this->alnStorage[1] = new char*[this->seqNum]; 
                for (int i = 0; i < this->seqNum; ++i) {
                    this->alnStorage[0][i] = nullptr;
                    this->alnStorage[1][i] = nullptr;
                }
                for (auto i: this->subtrees[this->subtreeIdx]) {
                    this->alnStorage[0][i] = new char[adjustLen];
                    this->alnStorage[1][i] = new char[adjustLen];
                    for (int j = 0; j < adjustLen; ++j) {
                        if (j < this->seqLen) {
                            this->alnStorage[0][i][j] = this->seqs[i][j];
                            this->alnStorage[1][i][j] = this->seqs[i][j];
                        }
                        else {
                            this->alnStorage[0][i][j] = 0;
                            this->alnStorage[1][i][j] = 0;
                        }
                        // this->alnStorage[1][i][j] = 0;
                    }
                    delete [] this->seqs[i];
                    this->seqs[i] = nullptr;
                }
            }    
            this->memLen = adjustLen;
            printf("Allocate memory for alignment, size = %lu x %lu\n", memNum, memLen);
            return;
        }
        
        void memCheck(int alnLen) {
            if (alnLen > this->memLen) {
                printf("Reallocate Memory. SeqLen (%d) > MemLen (%lu)\n", alnLen, this->memLen);
                this->alnMalloc(alnLen);
            }
            return;
        }
    };
} 

#endif