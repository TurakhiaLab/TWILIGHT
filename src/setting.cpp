#ifndef MSA_HPP
#include "setting.hpp"
#endif


void msa::utility::changeStorage(int idx) {
    this->seqsStorage[idx] = (seqsStorage[idx] == 0) ? 1 : 0;
    return;
}

void msa::utility::setSubtreeIdx(int idx) {
    this->subtreeIdx = idx;
    return;
}

void msa::utility::seqsFree(){
    for (size_t i = 0; i < this->memNum; ++i) {
        if (this->alnStorage[0][i] != nullptr) delete [] this->alnStorage[0][i];
        if (this->alnStorage[1][i] != nullptr) delete [] this->alnStorage[1][i];
    }
    delete [] this->alnStorage[0];
    delete [] this->alnStorage[1];
    this->alnStorage[0] = nullptr;
    this->alnStorage[1] = nullptr;
    return;
}

void msa::utility::seqFree(int i) {
    delete [] this->alnStorage[0][i];
    delete [] this->alnStorage[1][i];
    this->alnStorage[0][i] = nullptr;
    this->alnStorage[1][i] = nullptr;
    return;
}

void msa::utility::seqMalloc(int seqNum, int seqLen){
    if (this->alnStorage[0] != nullptr && this->alnStorage[1] != nullptr) {
        char** temp[2] = {nullptr, nullptr};
        int adjustLen = static_cast<int>(seqLen * this->timesBigger);
        temp[0] = new char*[seqNum]; 
        temp[1] = new char*[seqNum]; 
        // for (int i = 0; i < seqNum; ++i) {
        //     temp[0][i] = new char[adjustLen];
        //     temp[1][i] = new char[adjustLen];
        // }
        for (int i = 0; i < seqNum; ++i) {
            temp[0][i] = new char[adjustLen];
            temp[1][i] = new char[adjustLen];
            for (int j = 0; j < adjustLen; ++j) {
                temp[0][i][j] = (j < memLen) ? alnStorage[0][i][j] : 0;
                temp[1][i][j] = (j < memLen) ? alnStorage[1][i][j] : 0;
            }
            delete [] this->alnStorage[0][i];
            delete [] this->alnStorage[1][i];
        }
        delete [] this->alnStorage[0];
        delete [] this->alnStorage[1];
        this->alnStorage[0] = temp[0];
        this->alnStorage[1] = temp[1];
        this->memLen = adjustLen;
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
    printf("Allocate memory... Size = %lu x %lu\n", memNum, memLen);
    return;
}

void msa::utility::seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap){
    this->seqMalloc(seqsMap.size(), seqLen);
    int s = 0;
    for (auto seq: seqsMap) {
        for (int i = 0; i < this->memLen; ++i) {
            if (i < seq.second.first.size()) this->alnStorage[0][s][i] = seq.second.first[i];
            else                             this->alnStorage[0][s][i] = 0; 
            this->alnStorage[1][s][i] = 0;  
        }
        // if (s == 0) std::cout << seq.first << '\n' << s << '\n' << seq.second.first << '\n';
        this->seqsIdx[seq.first] = s;
        this->seqsName[s] = seq.first;
        this->seqsLen[seq.first] = seq.second.first.size();
        this->seqsStorage[s] = 0;
        ++s;
    }
    printf("Allocate memory for sequence storage, size = %lu x %lu\n", this->memNum, this->memLen);
    return;
}

void msa::utility::memCheck(int seqLen) {
    if (seqLen > this->memLen) {
        printf("Reallocate Memory. SeqLen (%d) > MemLen (%ld)\n", seqLen, this->memLen);
        seqMalloc(this->memNum, seqLen);
    }
    return;
}

void msa::utility::clearAll() {
    this->rawSeqs.clear();
    this->badSequences.clear();
    this->subtrees.clear();
    this->seqsIdx.clear();
    this->seqsName.clear();
    this->seqsLen.clear();
    this->seqsStorage.clear();
}

void msa::utility::debug() {
    int alnLen = 0, offset = 0;
    bool theFirst = true;
    for (auto s: this->rawSeqs) {
        std::string r = "";
        int sIdx = this->seqsIdx[s.first];
        int storage = this->seqsStorage[sIdx];
        offset = 0;
        while (this->alnStorage[storage][sIdx][offset] != 0) {
            if (this->alnStorage[storage][sIdx][offset] != '-') {
                r += this->alnStorage[storage][sIdx][offset];
            }
            ++offset;
        }
        if (theFirst) {alnLen = offset; theFirst = false;}
        else {
            if (alnLen != offset) printf("seq: %s, the sequence length (%d) did not match (%d)\n", s.first.c_str(), offset, alnLen);
        }
        if (r != s.second) {
            printf("seq: %s, the sequence did not match\n", s.first.c_str());
            if (r.size() != s.second.size()) printf("Wrong length. %lu/%lu.\n",r.size(), s.second.size());
            for (int i = 0; i < s.second.size(); ++i) {
                if (r[i] != s.second[i]) {
                    printf("Mismatch at position %d.\n", i);
                    break;
                }
            }                
        }
        
    }
    this->rawSeqs.clear();
        
}