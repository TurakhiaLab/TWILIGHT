#ifndef MSA_HPP
#define MSA_HPP

#include <string>
#include <unordered_map>
#include <map>

// namespace msa
// {
//     struct utility
//     {
//         std::map<int, std::vector<Node*>> badSequences;
//         std::map<int, std::vector<int>> subtrees;
//         std::unordered_map<std::string, int>  seqsIdx;
//         std::unordered_map<std::string, int>  seqsLen;
//         std::unordered_map<int, bool> seqsStorage;
//         std::map<int, std::vector<std::vector<int32_t>>> profileFreq;
//         char** alnStorage[2] = {nullptr, nullptr};
//         char** seqs = nullptr;
//         int nowProcess = 0;
//         size_t subtreeIdx = 0;
//         size_t gpuNum;
//         size_t seqLen;
//         size_t seqNum;
//         size_t memLen;
//         size_t memNum;
//         int nowStore = 0;
//         const float timesBigger = 1.3;
//         void changeStorage(int idx) {
//             seqsStorage[idx] = (seqsStorage[idx] == 0) ? 1 : 0;
//             return;
//         }
//         void setSubtreeIdx(int idx) {
//             this->subtreeIdx = idx;
//             return;
//         }
//         void seqsFree(){
//             for (int i = 0; i < this->seqNum; ++i) {
//                 delete [] this->seqs[i];
//             }
//             delete [] this->seqs;
//             this->seqs = nullptr;
//             return;
//         }
//         void reset() {
//             seqsFree();
//             subtrees.clear();
//             seqsIdx.clear();
//             seqsLen.clear();
//             seqsStorage.clear();
//             seqLen = 0;
//             seqNum = 0;
//             memLen = 0;
//             memNum = 0;
            
//         }

//         void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap){
//             this->seqNum = seqsMap.size();
//             this->seqs = new char*[this->seqNum]; 
//             this->seqLen = seqLen;
            
//             int s = 0;
//             for (auto seq: seqsMap) {
//                 this->seqs[s] = new char[seqLen];
//                 std::string seqName = seq.first;
//                 int subtreeIdx = seq.second.second;
//                 if (this->subtrees.find(subtreeIdx) == this->subtrees.end()) {
//                     std::vector<int> temp;
//                     this->subtrees[subtreeIdx] = temp;
//                 }
//                 this->subtrees[subtreeIdx].push_back(s);
//                 for (int j = 0; j < seqLen; ++j) {
//                     this->seqs[s][j] =  (j < seq.second.first.size()) ? seq.second.first[j] : 0;
//                 }
//                 this->seqsIdx[seqName] = s;
//                 this->seqsLen[seqName] = seq.second.first.size();
//                 ++s;
//                 seqsMap.erase(seqName);
//             }
//             assert(s == seqNum);
//             this->seqNum = seqNum;
//             this->seqLen = seqLen;
//             printf("Allocate memory for sequence storage, size = %lu x %lu\n", this->seqNum, seqLen);
//             return;
//         }
//         void updateSeqs() {
//             size_t len = max(this->memLen, this->seqLen);
//             std::cout << "update Len: " << len << '\n';
//             for (int i = 0; i < this->seqNum; ++i) {
//                 char* temp = new char[len];
//                 if (std::find(this->subtrees[this->subtreeIdx].begin(), this->subtrees[this->subtreeIdx].end(), i) != this->subtrees[this->subtreeIdx].end()) {
//                     int s = this->seqsStorage[i];
//                     for (int j = 0; j < len; ++j) {
//                         temp[j] = (j < this->memLen) ? this->alnStorage[s][i][j] : 0;
//                     }
//                 }
//                 else {
//                     for (int j = 0; j < len; ++j) {
//                         temp[j] = (j < this->seqLen) ? this->seqs[i][j] : 0;
//                     }
//                     delete [] this->seqs[i];
//                     this->seqs[i] = nullptr;
//                 }
                
//                 this->seqs[i] = temp;
//             }
//             this->alnFree();
//             this->seqLen = len;
//         }
//         void alnFree(){
//             for (auto i: this->subtrees[this->subtreeIdx]) {
//                 delete [] this->alnStorage[0][i];
//                 delete [] this->alnStorage[1][i];
//             }
//             delete [] this->alnStorage[0];
//             delete [] this->alnStorage[1];
//             this->alnStorage[0] = nullptr;
//             this->alnStorage[1] = nullptr;
//             this->memLen = 0;
//             this->memNum = 0;
//             return;
//         }
//         void alnMalloc(size_t alnLen){
//             int adjustLen = static_cast<int>(alnLen * this->timesBigger);
//             if (alnStorage[0] != nullptr && alnStorage[1] != nullptr) {
//                 char** temp[2] = {nullptr, nullptr};
//                 temp[0] = new char*[this->seqNum]; 
//                 temp[1] = new char*[this->seqNum];
//                 for (int i = 0; i < this->seqNum; ++i) {
//                     temp[0][i] = nullptr;
//                     temp[1][i] = nullptr;
//                 }
//                 for (auto i: this->subtrees[this->subtreeIdx]) {
//                     temp[0][i] = new char[adjustLen];
//                     temp[1][i] = new char[adjustLen];
//                     for (int j = 0; j < adjustLen; ++j) {
//                         temp[0][i][j] = (j < memLen) ? this->alnStorage[0][i][j] : 0;
//                         temp[1][i][j] = (j < memLen) ? this->alnStorage[1][i][j] : 0;
//                     }
//                     delete [] this->alnStorage[0][i];
//                     delete [] this->alnStorage[1][i];
                    
//                 }
//                 delete [] this->alnStorage[0];
//                 delete [] this->alnStorage[1];
//                 this->alnStorage[0] = temp[0];
//                 this->alnStorage[1] = temp[1];
//             }
//             else {
//                 this->memNum = this->subtrees[this->subtreeIdx].size();
//                 this->alnStorage[0] = new char*[this->seqNum]; 
//                 this->alnStorage[1] = new char*[this->seqNum]; 
//                 for (int i = 0; i < this->memNum; ++i) {
//                     this->alnStorage[0][i] = nullptr;
//                     this->alnStorage[1][i] = nullptr;
//                 }
//                 for (auto i: this->subtrees[this->subtreeIdx]) {
//                     this->alnStorage[0][i] = new char[adjustLen];
//                     this->alnStorage[1][i] = new char[adjustLen];
//                     for (int j = 0; j < adjustLen; ++j) {
//                         if (j < this->seqLen) {
//                             this->alnStorage[0][i][j] = this->seqs[i][j];
//                             this->alnStorage[1][i][j] = this->seqs[i][j];
//                         }
//                         else {
//                             this->alnStorage[0][i][j] = 0;
//                             this->alnStorage[1][i][j] = 0;
//                         }
//                         // this->alnStorage[1][i][j] = 0;
//                     }
//                     delete [] this->seqs[i];
//                     this->seqs[i] = nullptr;
//                 }
//             }    
//             this->memLen = adjustLen;
//             printf("Allocate memory for alignment, size = %lu x %lu\n", memNum, memLen);
//             return;
//         }
        
//         void memCheck(int alnLen) {
//             if (alnLen > this->memLen) {
//                 printf("Reallocate Memory. SeqLen (%d) > MemLen (%lu)\n", alnLen, this->memLen);
//                 this->alnMalloc(alnLen);
//             }
//             return;
//         }
//     };
// } 

namespace msa
{
    struct option
    {
        int gpuNum;
        int cpuNum;
        int maxSubtree;
        int maxSubSubtree;
        std::vector<int> gpuIdx;
        int nowProcess;
        bool debug;
        bool readBatches;
    };

    struct utility
    {
        std::map<int, std::vector<std::string>> badSequences;
        std::map<int, std::vector<int>> subtrees;
        std::unordered_map<std::string, int>  seqsIdx;
        std::unordered_map<std::string, int>  seqsLen;
        std::unordered_map<int, bool> seqsStorage;
        std::map<int, std::vector<std::vector<int32_t>>> profileFreq;
        
        char** alnStorage[2] = {nullptr, nullptr};
        int nowProcess = 0;
        size_t subtreeIdx = 0;
        size_t gpuNum;
        size_t memLen;
        size_t memNum;
        size_t seqLen;
        size_t seqNum;
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
            for (int i = 0; i < this->memNum; ++i) {
                delete [] alnStorage[0][i];
                delete [] alnStorage[1][i];
            }
            delete [] alnStorage[0];
            delete [] alnStorage[1];
            alnStorage[0] = nullptr;
            alnStorage[1] = nullptr;
            return;
        }
        
        void seqMalloc(int seqNum, int seqLen){
            if (alnStorage[0] != nullptr && alnStorage[1] != nullptr) {
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
                            temp[0][j][i] = alnStorage[0][j][i];
                            temp[1][j][i] = alnStorage[1][j][i];
                        }
                        else {
                            temp[0][j][i] = 0;
                            temp[1][j][i] = 0;
                        }
                    }
                }
                seqsFree();
                alnStorage[0] = temp[0];
                alnStorage[1] = temp[1];
                memLen = adjustLen;
            }
            else {
                int adjustLen = static_cast<int>(seqLen * this->timesBigger); 
                this->alnStorage[0] = new char*[seqNum]; 
                this->alnStorage[1] = new char*[seqNum]; 
                for (int i = 0; i < seqNum; ++i) {
                    this->alnStorage[0][i] = new char[adjustLen];
                    this->alnStorage[1][i] = new char[adjustLen];
                }
                for (int j = 0; j < seqNum; ++j) {
                    for (int i = 0; i < adjustLen; ++i) {
                        this->alnStorage[0][j][i] = 0;
                        this->alnStorage[1][j][i] = 0;
                    }
                }
                this->memNum = seqNum;
                this->memLen = adjustLen;
            }    
            printf("Allocate memory... Size = %d x %d\n", memNum, memLen);
            return;
        }
        void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap){
            this->seqMalloc(seqsMap.size(), seqLen);
            int s = 0;
            for (auto seq: seqsMap) {
                for (int i = 0; i < this->memLen; ++i) {
                    if (i < seq.second.first.size()) this->alnStorage[0][s][i] = seq.second.first[i];
                    else                             this->alnStorage[0][s][i] = 0; 
                    this->alnStorage[1][s][i] = 0;  
                }
                // std::cout << seq.first << ':' << s << ':' << seq.second.size() << '\n';
                this->seqsIdx[seq.first] = s;
                this->seqsLen[seq.first] = seq.second.first.size();
                this->seqsStorage[s] = 0;
                ++s;
            }
            printf("Allocate memory for sequence storage, size = %lu x %lu\n", this->memNum, this->memLen);
            return;
        }
        
        void memCheck(int seqLen) {
            if (seqLen > this->memLen) {
                printf("Reallocate Memory. SeqLen (%d) > MemLen (%d)\n", seqLen, this->memLen);
                seqMalloc(this->memNum, seqLen);
            }
            return;
        }
    };
} 

#endif