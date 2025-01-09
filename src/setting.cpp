#ifndef SETTING_HPP
#include "setting.hpp"
#endif

Params::Params(po::variables_map& vm) {
    bool userDefine = vm.count("matrix");
    float gapOp = vm["gap-open"].as<float>();
    float gapEx = vm["gap-extend"].as<float>();
    float mat = vm["match"].as<float>();
    float mis = vm["mismatch"].as<float>();
    float trans = vm["transition"].as<float>();
    float xdrop = round(vm["xdrop"].as<float>());
    float gapCl = gapOp;
    
    if (gapOp > 0 || gapEx > 0) {
        std::cerr << "ERROR: Gap penalties should be less or equal than 0.\n";
        exit(1);
    }
    
    if (!userDefine) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4)        this->scoringMatrix[i][j] = vm.count("wildcard") ? mat : 0.0;
                else if (i == j)             this->scoringMatrix[i][j] = mat;
                else if (std::abs(i-j) == 2) this->scoringMatrix[i][j] = trans;
                else                         this->scoringMatrix[i][j] = mis;
            }
        }
        this->gapOpen = gapOp; this->gapExtend = gapEx; this->gapClose = gapCl;
        this->xdrop = (this->gapExtend == 0) ? xdrop : -1*xdrop*this->gapExtend;
    }
    else {
        std::string matrixFileName = vm["matrix"].as<std::string>();
        std::ifstream matrixFile(matrixFileName);
        if (!matrixFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", matrixFileName.c_str());
            exit(1);
        }
        std::string dna;
        std::vector<int> dnaVec;
        int readCount = 0;
        while (matrixFile >> dna) {
            if (readCount < 4) {
                if      (dna == "A" || dna == "a") dnaVec.push_back(0);
                else if (dna == "C" || dna == "c") dnaVec.push_back(1);
                else if (dna == "G" || dna == "g") dnaVec.push_back(2);
                else if (dna == "T" || dna == "t" || dna == "U" || dna == "u") dnaVec.push_back(3);
                else { std::cerr << "ERROR: Illegal nucleotide type \"" << dna << "\"\n"; exit(1);}
                readCount++;
            }
            else {
                int x = (readCount-4) / 4;
                int y = (readCount-4) % 4;
                int i = dnaVec[x];
                int j = dnaVec[y];
                this->scoringMatrix[i][j] = std::stof(dna.c_str());
                readCount++;
            }
        }
        matrixFile.close();
        float Nscore = vm.count("wildcard") ? (this->scoringMatrix[0][0] + this->scoringMatrix[1][1] + this->scoringMatrix[2][2] + this->scoringMatrix[3][3]) / 4 : 0.0;
        for (int i = 0; i < 5; ++i) {
            this->scoringMatrix[i][4] = Nscore;
            this->scoringMatrix[4][i] = Nscore;
        }
        this->gapOpen = gapOp; this->gapExtend = gapEx; this->gapClose = gapCl;
        this->xdrop =  (this->gapExtend == 0) ? xdrop : -1*xdrop*this->gapExtend;
    }
    if (vm.count("verbose")) {
        std::cout << "======== Parameters ========\n";
        std::cout << "    A    C    G    T    N   \n";
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                std::cout << std::setw(5) << this->scoringMatrix[i][j];
            }
            std::cout << '\n';
        }
        std::cout << "GapOpen: " << this->gapOpen << " / GapExtend: " << this->gapExtend << " / Xdrop: " << this->xdrop << '\n';
    }
    std::cout << "============================\n";
}

msa::option::option(po::variables_map& vm) {
    if (!vm.count("tree") && !vm.count("sequences") && !vm.count("files")) {
        std::cerr << "ERROR: No input file.\n";
        exit(1);
    }
    if (!vm.count("files") && (!vm.count("sequences") || !vm.count("tree"))) {
        std::cerr << "ERROR: A tree file and a raw sequeunce file are required for building MSA from raw seqeunces.\n";
        exit(1);
    }
    if (vm.count("sequences") && vm.count("files")) {
        std::cerr << "ERROR: Both modes cannot run at the same time.\n";
        exit(1);
    }
    if (vm.count("tree")) this->treeFile = vm["tree"].as<std::string>();
    if (vm.count("sequences")) {this->seqFile = vm["sequences"].as<std::string>(); this->alnMode = 0;}
    if (vm.count("files")) {this->msaDir = vm["files"].as<std::string>(); this->alnMode = 1;}
    if (!vm.count("output")) {
        std::cerr << "ERROR: An output file name is required.\n";
        exit(1);
    }
    
    this->outFile = vm["output"].as<std::string>();
    if (fs::exists(this->outFile)) {
        std::cerr << "ERROR: " << this->outFile << " already exists. Please use another file name.\n";
        exit(1);
    }
    std::ofstream outFile(this->outFile);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", this->outFile.c_str());
        exit(1);
    }
    outFile.close();

    int maxSubtreeSize = (vm.count("max-subtree")) ? vm["max-subtree"].as<int>() : INT32_MAX;
    int maxSubSubtreeSize = (vm.count("max-group")) ? vm["max-group"].as<int>() : INT32_MAX;
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    int cpuNum = (vm.count("cpu")) ? vm["cpu"].as<int>() : maxCpuThreads;
    if (cpuNum <= 0) {
        std::cerr << "ERROR: Requested cpu cores <= 0.\n";
        exit(1);
    }
    if (cpuNum > maxCpuThreads) {
        std::cerr << "ERROR: Requested cpu cores more available threads.\n";
        exit(1);
    }
    float gappyVertical = vm["rgc"].as<float>();
    float gappyHorizon;
    if (gappyVertical > 1 || gappyVertical <= 0) {
        std::cerr << "ERROR: The value of rgc should be in (0,1]\n";
        exit(1);
    }
    if (vm.count("gappy-horizon")) {
        gappyHorizon = vm["gappy-horizon"].as<float>();
        if (gappyHorizon - round(gappyHorizon) != 0 || gappyHorizon < 1) {
            std::cerr << "ERROR: The value of gappy-horizon should be an positive integer\n";
            exit(1);
        }
    }
    else {
        gappyHorizon = 1;
    }
    
    if (maxSubtreeSize < INT32_MAX || vm.count("files")) {
        std::string tempDir;
        if (!vm.count("temp-dir")) {
            int idx = 1;
            std::string tempDir_org = "./twilight_temp";
            tempDir = tempDir_org;
            while (true) {
                if (mkdir(tempDir.c_str(), 0777) == -1) {
                    if( errno == EEXIST ) {
                        tempDir = tempDir_org + '_' + std::to_string(idx);
                        ++idx;
                    }
                    else { fprintf(stderr, "ERROR: Can't create directory: %s\n", tempDir.c_str()); exit(1); }
                }
                else break;
            }
            
        }
        else {
            tempDir = vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            if (mkdir(tempDir.c_str(), 0777) == -1) {
                if( errno == EEXIST ) {
                    std::cerr << "ERROR: " << tempDir << " already exists. In order to prevent your file from being overwritten, please delete this folder or use another folder name.\n";
                    exit(1);
                }
                else { fprintf(stderr, "ERROR: Can't create directory: %s\n", tempDir.c_str()); exit(1); }
            }
        }
        std::cout << tempDir << " created for storing temporary alignments\n";
        this->tempDir = tempDir;
    }
    std::string outType = vm["output-type"].as<std::string>();
    if (outType != "FASTA" && outType != "CIGAR") {
        std::cerr << "ERROR: Unrecognized output type \"" << outType <<"\"\n";
        exit(1);
    }
    std::string merger = vm["merge"].as<std::string>();
    if (merger == "T" || merger == "t") merger = "transitivity";
    else if (merger == "P" || merger == "p") merger = "profile";
    else {
        std::cerr << "ERROR: Unrecognized method to merge subtrees \"" << merger <<"\"\n";
        exit(1);
    }
    
    if (vm.count("psgop")) {
        std::string psgop = vm["psgop"].as<std::string>();
        if      (psgop == "y" || psgop == "Y" || psgop == "yes" || psgop == "YES" || psgop == "Yes") {this->psgopAuto = false; this->psgop = true;}
        else if (psgop == "n" || psgop == "N" || psgop == "no" || psgop == "NO" || psgop == "No")    {this->psgopAuto = false; this->psgop = false;}
        else if (psgop == "auto" || psgop == "Auto" || psgop == "AUTO")                              {this->psgopAuto = true;  this->psgop = false;}
        else {
            std::cerr << "ERROR: Unrecognized option \"" << psgop <<"\" for position-specific gap open penalty.\n";
            exit(1);
        } 
    }
    else {
        this->psgop = false;
        this->psgopAuto = true;
    }
    this->redo = false;
    this->calSim = false;
    this->cpuNum = cpuNum; 
    this->gpuNum = 0;
    this->gpuIdx = std::vector<int> (0);
    this->maxSubtree = maxSubtreeSize;
    this->maxSubSubtree = maxSubSubtreeSize;
    this->debug = vm.count("check");
    this->gappyHorizon = gappyHorizon;
    this->gappyVertical = gappyVertical;
    this->cpuOnly = vm.count("cpu-only");
    this->outType = outType;
    this->merger = merger;
    this->printDetail = vm.count("verbose");
    this->deleteTemp = !vm.count("keep-temp");

    
    std::cout << "====== Configuration =======\n";
    if (this->maxSubtree != INT32_MAX) 
    std::cout << "Max-subtree: " << this->maxSubtree << '\n';
    if (this->maxSubSubtree != INT32_MAX) 
    std::cout << "Max-group: " << this->maxSubSubtree << '\n';
    if (this->gappyVertical == 1) 
    std::cout << "Disable removing gappy columns.\n";
    else
    std::cout << "Removing gappy columns: " << this->gappyVertical << '\n';
    printf("Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);
    
}

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
    delete [] this->seqMemLen;
    this->alnStorage[0] = nullptr;
    this->alnStorage[1] = nullptr;
    this->memNum = 0;
    this->memLen = 0;
    return;
}

void msa::utility::seqFree(int i) {
    delete [] this->alnStorage[0][i];
    delete [] this->alnStorage[1][i];
    this->alnStorage[0][i] = nullptr;
    this->alnStorage[1][i] = nullptr;
    return;
}

void msa::utility::seqMalloc(int seqNum, int seqLen, option* option){
    if (this->alnStorage[0] != nullptr && this->alnStorage[1] != nullptr) {
        char** temp[2] = {nullptr, nullptr};
        int adjustLen = static_cast<int>(seqLen * this->timesBigger);
        temp[0] = new char*[seqNum]; 
        temp[1] = new char*[seqNum]; 
        for (int i = 0; i < seqNum; ++i) {
            temp[0][i] = new char[adjustLen];
            temp[1][i] = new char[adjustLen];
            for (int j = 0; j < adjustLen; ++j) {
                temp[0][i][j] = (j < this->seqMemLen[i]) ? alnStorage[0][i][j] : 0;
                temp[1][i][j] = (j < this->seqMemLen[i]) ? alnStorage[1][i][j] : 0;
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
        this->seqMemLen = new int[seqNum];
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
    if (option->printDetail) printf("Allocate memory... Size = %d x %d\n", memNum, memLen);
    return;
}

void msa::utility::seqMalloc(int seqLen, int idx){
    int adjustLen = static_cast<int>(seqLen * this->timesBigger);
    char* temp_0 = new char[adjustLen];
    char* temp_1 = new char[adjustLen];
    for (int j = 0; j < adjustLen; ++j) {
        temp_0[j] = (j < this->seqMemLen[idx]) ? alnStorage[0][idx][j] : 0;
        temp_1[j] = (j < this->seqMemLen[idx]) ? alnStorage[1][idx][j] : 0;
    }
    delete [] this->alnStorage[0][idx];
    delete [] this->alnStorage[1][idx];
    this->alnStorage[0][idx] = temp_0;
    this->alnStorage[1][idx] = temp_1;
    this->seqMemLen[idx] = adjustLen;
    // if (adjustLen > this->memLen) this->memLen = adjustLen;
    return;
}

void msa::utility::seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap, option* option){
    this->seqMalloc(seqsMap.size(), seqLen, option);
    int s = 0;
    for (auto seq: seqsMap) {
        for (int i = 0; i < this->memLen; ++i) {
            if (i < seq.second.first.size()) this->alnStorage[0][s][i] = seq.second.first[i];
            else                             this->alnStorage[0][s][i] = 0; 
            this->alnStorage[1][s][i] = 0;  
        }
        this->seqsIdx[seq.first] = s;
        this->seqsName[s] = seq.first;
        this->seqsLen[seq.first] = seq.second.first.size();
        this->seqsStorage[s] = 0;
        this->seqMemLen[s] = this->memLen;
        
        ++s;
    }
    // printf("Allocate memory for sequence storage, size = %lu x %lu\n", this->memNum, this->memLen);
    return;
}

void msa::utility::memCheck(int seqLen, option* option) {
    if (seqLen > this->memLen) {
        if (option->printDetail) printf("Reallocate Memory. SeqLen (%d) > MemLen (%d)\n", seqLen, this->memLen);
        seqMalloc(this->memNum, seqLen, option);
    }
    return;
}

void msa::utility::memCheck(int seqLen, int idx) {
    if (seqLen > this->seqMemLen[idx]) {
        // printf("No: %d Reallocate Memory. SeqLen (%d) > MemLen (%ld)\n",idx, seqLen, this->seqMemLen[idx]);
        seqMalloc(seqLen, idx);
    }
    return;
}

void msa::utility::storeCIGAR() {
    tbb::spin_rw_mutex writeMutex;
    tbb::parallel_for(tbb::blocked_range<int>(0, this->seqsIdx.size()), [&](tbb::blocked_range<int> range){ 
    for (int n = range.begin(); n < range.end(); ++n) {
        auto seq = std::make_pair(this->seqsName[n], n);
    // for (auto seq: this->seqsIdx) {
        std::string seqName = seq.first;
        int sIdx = seq.second;
        int storage = this->seqsStorage[sIdx];
        std::string cigar = "";
        char type = (this->alnStorage[storage][sIdx][0] == '-') ? 'D' : 'M';
        int num = 0;
        int i = 0;
        while (this->alnStorage[storage][sIdx][i] != 0) {
            if (this->alnStorage[storage][sIdx][i] == '-') {
                if (type == 'M') {
                    cigar += (std::to_string(num) + type);
                    type = 'D'; num = 1;
                }
                else if (type == 'D') ++num;
            }
            else {
                if (type == 'M') ++num;
                else {
                    cigar += (std::to_string(num) + type);
                    type = 'M'; num = 1;
                }
            }
            ++i;
        }
        cigar += (std::to_string(num) + type);
        {
            tbb::spin_rw_mutex::scoped_lock lock(writeMutex);
            this->seqsCIGAR[seqName] = cigar;
        }
        // this->seqsCIGAR[seqName] = cigar;
    }
    });
    return;
}

void msa::utility::clearAll() {
    this->rawSeqs.clear();
    this->badSequences.clear();
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
            if (alnLen != offset) printf("%s: the sequence length (%d) did not match the MSA length(%d)\n", s.first.c_str(), offset, alnLen);
        }
        if (r != s.second) {
            printf("%s: after removing the gaps, the alignment did not match the original sequence.\n", s.first.c_str());
            // if (r.size() != s.second.size()) printf("Wrong length. %lu/%lu.\n",r.size(), s.second.size());
            // for (int i = 0; i < s.second.size(); ++i) {
            //     if (r[i] != s.second[i]) {
            //         printf("Mismatch at position %d.\n", i);
            //         break;
            //     }
            // }                
        }
        
    }
    this->rawSeqs.clear();
        
}