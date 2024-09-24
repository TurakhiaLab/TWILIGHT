#ifndef UTIL_HPP
#include "util.cuh"
#endif

KSEQ_INIT2(, gzFile, gzread);
// po::options_description mainDesc("MSA Command Line Arguments");

Params* setParameters(po::variables_map& vm) {
    
    int userDefine = vm["scoring-matrix"].as<int>();
    Params* param = new Params();
    param->userDefine = userDefine;
    // Kimura
    paramType gapOp = vm["gap-open"].as<paramType>();
    paramType gapCl = vm["gap-close"].as<paramType>();
    paramType gapEx = vm["gap-extend"].as<paramType>();
    paramType mat = vm["match"].as<paramType>();
    paramType mis = vm["mismatch"].as<paramType>();
    paramType trans = vm["trans"].as<paramType>();
    
    paramType xdrop = round(vm["xdrop"].as<paramType>());
    
    if (userDefine == 0) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4)        param->scoringMatrix[i][j] = 0;
                else if (i == j)             param->scoringMatrix[i][j] = mat;
                else if (std::abs(i-j) == 2) param->scoringMatrix[i][j] = mis+(trans-1)*(mat-mis)/(trans);
                else                         param->scoringMatrix[i][j] = mis;
                // std::cout << std::setw(10) << param->scoringMatrix[i][j] << ' ';
            }
            // std::cout << '\n';
        }
        param->gapOpen = gapOp; param->gapExtend = gapEx; param->gapClose = gapCl;
        param->xdrop = (param->gapExtend == 0) ? xdrop : -1*xdrop*param->gapExtend;
    }
    else if (userDefine == 1) {
        int i, j;
        // generatenuc1pam from MAFFT
        float R[4][4], mut[4], total = 0.0, temp;
        float pam1[4][4];
        float freq [4] = {0.25, 0.25, 0.25, 0.25};
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j) {
                if (i == j)                  R[i][j] = 0.0;
                else if (std::abs(i-j) == 2) R[i][j] = trans;
                else                         R[i][j] = 1.0;
            }
        }
	    for(i = 0; i < 4; ++i ) {
            temp = 0.0;
            for(j = 0; j < 4; ++j ) temp += R[i][j] * freq[j];
            mut[i] = temp;
	    	total += temp * freq[i];
        }
	    for (i = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j) {
                if (i != j) pam1[i][j] = 0.01 / total * R[i][j] * freq[j];
                else        pam1[i][j] = 1.0 - 0.01 / total * mut[i];
            }
        }
        // MtxmltDouble from MAFFT
        float pamx [4][4];
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] = (i == j) ? 1.0 : 0.0;
        int PAM_N = vm["pam-n"].as<int>();
        for (int t = 0; t < PAM_N; ++t) {
            int k;
            float temp [4], s;
            for (i = 0; i < 4; ++i) {
                for (k = 0; k < 4; ++k) temp[k] = pamx[i][k];
                for (j = 0; j < 4; ++j) {
                    s = 0;
                    for (k = 0; k < 4; ++k) s += temp[k] * pam1[k][j];
                    pamx[i][j] = s;
                }
            }
        }
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] /= freq[j];
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] = (pamx[i][j] == 0.0) ? std::log10(0.00001) * 1000.0 : std::log10( pamx[i][j] ) * 1000.0;
        double avg1 = 0, avg2 = 0;
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) avg2 += freq[i]*freq[j]*pamx[i][j];
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] -= avg2;
        for (i = 0; i < 4; ++i) avg1 += 0.25*pamx[i][i];
        int offset = 59;
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] *= 600 / avg1;
        for (i = 0; i < 4; ++i) for (j = 0; j < 4; ++j) pamx[i][j] -= offset;
        
        for (i = 0; i < 5; ++i) for (j = 0; j < 5; ++j) 
            param->scoringMatrix[i][j] = (i == 4 || j == 4) ? 0.0 : std::round(pamx[i][j] / 10.0);
        
        
        float gop = 1.53; float ge = 0.123;
        param->gapOpen = std::round(-gop*100); param->gapExtend = std::round(-ge*100);
        param->gapClose = param->gapOpen;
        param->xdrop = (param->gapExtend == 0) ? xdrop : -1*xdrop*param->gapExtend;
    }
    else if (userDefine == 2) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                param->scoringMatrix[i][j] = param->userMatrix[i][j];
            }
        }
        param->gapOpen = param->userGapOpen; param->gapExtend = param->userGapExtend;
        param->gapOpen = param->userGapClose;
        param->xdrop = (param->gapExtend == 0) ? xdrop : -1*xdrop*param->gapExtend;
    }
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            std::cout << std::setw(5) << param->scoringMatrix[i][j];
        }
        std::cout << '\n';
    }
    
    std::cout << "GapOpen: " << param->gapOpen << " / GapExtend: " << param->gapExtend << " / GapClose: " << param->gapClose << " / Xdrop: " << param->xdrop << '\n';
    // exit(1);
    return param;
}

void setOptions(po::variables_map& vm, msa::option* option) {
    int maxSubtreeSize = vm["max-subtree-size"].as<int>();
    int maxSubSubtreeSize = (vm.count("max-leaves")) ? vm["max-leaves"].as<int>() : INT_MAX;
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    int cpuNum = (vm.count("cpu-num")) ? vm["cpu-num"].as<int>() : maxCpuThreads;
    if (cpuNum <= 0) {
        std::cerr << "ERROR: requested cpu threads <= 0.\n";
        exit(1);
    }
    if (cpuNum > maxCpuThreads) {
        std::cerr << "ERROR: requested cpu threads more available threads.\n";
        exit(1);
    }
    printf("Maximum available CPU threads: %d. Using %d CPU threads.\n", maxCpuThreads, cpuNum);
    tbb::task_scheduler_init init(cpuNum);
    int maxGpuNum;
    cudaGetDeviceCount(&maxGpuNum);
    int gpuNum = (vm.count("gpu-num")) ? vm["gpu-num"].as<int>() : maxGpuNum;
    if (gpuNum <= 0) {
        std::cerr << "ERROR: requested number of GPU <= 0.\n";
        exit(1);
    }
    if (gpuNum > maxGpuNum) {
        std::cerr << "ERROR: requested number of GPU more than available GPUs.\n";
        exit(1);
    }
    std::vector<int> gpuIdx;
    if (vm.count("gpu-index")) {
        std::string gpuIdxString = vm["gpu-index"].as<std::string>();
        std::string id = "";
        for (int i = 0; i < gpuIdxString.size(); ++i) {
            if (gpuIdxString[i] == ',') {
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else if (i == gpuIdxString.size()-1) {
                id += gpuIdxString[i];
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else id += gpuIdxString[i];
        }
    }
    else {
        for (int i = 0; i < gpuNum; ++i) gpuIdx.push_back(i);
    }
    if (gpuIdx.size() != gpuNum) {
        std::cerr << "ERROR: the number of requested GPUs does not match the number of specified gpu indexes.\n";
        exit(1);
    }
    for (auto id: gpuIdx) {
        if (id >= maxGpuNum) {
            std::cerr << "ERROR: specified gpu index >= the number of GPUs\n";
            exit(1);
        }
    }
    printf("Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);
    float gappyVertical = vm["gappy-column"].as<float>();
    float gappyHorizon;
    if (gappyVertical > 1 || gappyVertical <= 0) {
        std::cerr << "ERROR: the value of gappy-column should be (0,1]\n";
        exit(1);
    }
    if (vm.count("gappy-length")) {
        gappyHorizon = vm["gappy-length"].as<float>();
        if (gappyHorizon - round(gappyHorizon) != 0 || gappyHorizon < 1) {
            std::cerr << "ERROR: the value of gappy-length should be an positive integer\n";
            exit(1);
        }
    }
    else {
        gappyHorizon = 0;
    }
    
    option->cpuNum = cpuNum; 
    option->gpuNum = gpuNum;
    option->gpuIdx = gpuIdx;
    option->maxSubtree = maxSubtreeSize;
    option->maxSubSubtree = maxSubSubtreeSize;
    option->debug = vm.count("debug");
    option->readBatches = vm.count("read-batches");
    option->gappyHorizon = gappyHorizon;
    option->gappyVertical = gappyVertical;
}

void readSequences(po::variables_map& vm, msa::utility* util, Tree* tree)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();

    std::string seqFileName = vm["sequences"].as<std::string>();
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    int seqNum = 0, maxLen = 0, minLen = INT_MAX;
    uint64_t totalLen = 0;
    
    std::map<std::string, std::pair<std::string, int>> seqs;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        std::string seqName = kseq_rd->name.s;
        int subtreeIdx = tree->allNodes[seqName]->grpID;
        seqs[seqName] = std::make_pair(std::string(kseq_rd->seq.s, seqLen), subtreeIdx);
        if (seqLen > maxLen) maxLen = seqLen;
        if (seqLen < minLen) minLen = seqLen;
        totalLen += seqLen;
    }
    
    seqNum = seqs.size();
    uint32_t avgLen = totalLen/seqNum;
    std::cout << "=== Sequence information ===\n";
    std::cout << "Number : " << seqNum << '\n';
    std::cout << "Max. Length: " << maxLen << '\n';
    std::cout << "Min. Length: " << minLen << '\n';
    std::cout << "Avg. Length: " << avgLen << '\n';
    std::cout << "============================\n";

    util->seqsMallocNStore(maxLen, seqs);
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

void readSequences(std::string seqFileName, msa::utility* util, Tree* tree)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    int seqNum = 0, maxLen = 0;
    uint64_t totalLen = 0;
    
    std::map<std::string, std::pair<std::string, int>> seqs;
    
    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        std::string seqName = kseq_rd->name.s;
        int subtreeIdx = tree->allNodes[seqName]->grpID;
        seqs[seqName] = std::make_pair(std::string(kseq_rd->seq.s, seqLen), subtreeIdx);
        if (seqLen > maxLen) maxLen = seqLen;
        totalLen += seqLen;
    }
    
    seqNum = seqs.size();
    uint32_t avgLen = totalLen/seqNum;
    std::cout << "(Num, MaxLen, AvgLen) = (" << seqNum << ", " << maxLen << ", " << avgLen << ")\n";

    util->seqsMallocNStore(maxLen, seqs);
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string tempDir;
    if (!vm.count("temp-dir")) {
        tempDir =  "./temp";
        if (mkdir(tempDir.c_str(), 0777) == -1) {
            if( errno == EEXIST ) {
                std::cout << tempDir << " already exists.\n";
            }
            else {
                fprintf(stderr, "ERROR: cant create directory: %s\n", tempDir.c_str());
                exit(1);
            }
            
        }
        else {
            std::cout << tempDir << " created\n";
        }
    }
    else tempDir = vm["temp-dir"].as<std::string>();
            
    if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);

    std::string seqFileName = vm["sequences"].as<std::string>();
    size_t maxLen = 0, totalLen = 0, seqNum = 0;
    std::cout << "Total " <<  partition->partitionsRoot.size() << " subtrees.\n";
    for (auto subroot: partition->partitionsRoot) {
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        Tree* subT = new Tree(subroot.second.first);
        gzFile f_rd = gzopen(seqFileName.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
            exit(1);
        }
        kseq_t* kseq_rd = kseq_init(f_rd);
        std::map<std::string, std::string> seqs;
        while (kseq_read(kseq_rd) >= 0) {
            size_t seqLen = kseq_rd->seq.l;
            std::string seqName = kseq_rd->name.s;
            if (subT->allNodes.find(seqName) != subT->allNodes.end()) {
                seqs[seqName] = std::string(kseq_rd->seq.s, seqLen);
                if (seqLen > maxLen) maxLen = seqLen;
                totalLen += seqLen;
                ++seqNum;
            }
        }
        std::string subtreeFileName = "subtree-" + std::to_string(subtreeIdx);
        std::string subtreeTreeFile = tempDir + '/' + subtreeFileName + ".nwk";
        outputSubtree(subtreeTreeFile, subT);
        std::string subtreeSeqFile = tempDir + '/' + subtreeFileName + ".raw.fa";
        outputSubtreeSeqs(subtreeSeqFile, seqs);
        seqs.clear();
        // subtreeIdx += 1;
        delete subT;
        
        kseq_destroy(kseq_rd);
        gzclose(f_rd);

    }

    uint32_t avgLen = totalLen/seqNum;
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Seqences: (Num, MaxLen, AvgLen) = (" << seqNum << ", " << maxLen << ", " << avgLen << ")\n";
    std::cout << "Created " << partition->partitionsRoot.size() << " subtree files in " << tempDir << " in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

void outputFinal (std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util, int& totalSeqs) {
    for (auto id: tree->root->msa) {
        // std::cout << "ID: " << id << '\n';
        int subtree = tree->allNodes[id]->grpID;
        std::string tempAlnFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".temp.aln";
        gzFile f_rd = gzopen(tempAlnFile.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", tempAlnFile.c_str());
            exit(1);
        }

        kseq_t* kseq_rd = kseq_init(f_rd);
        std::map<std::string, std::string> seqs;
        std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
        size_t seqLen = tree->allNodes[id]->msaAln.size();
        while (kseq_read(kseq_rd) >= 0) {
            size_t seqLen = kseq_rd->seq.l;
            std::string seqName = kseq_rd->name.s;
            std::string finalAln = "";
            std::string seq = std::string(kseq_rd->seq.s, seqLen);
            int orgIdx = 0;
            for (int j = 0; j < aln.size(); ++j) {
                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                    finalAln += seq[orgIdx];
                    ++orgIdx;
                }
                else {
                    finalAln += '-';
                }
            }
            seqs[seqName] = finalAln;
        }
        std::string subtreeSeqFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".final.aln";
        outputSubtreeSeqs(subtreeSeqFile, seqs);    
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        // std::remove(tempAlnFile.c_str());
        totalSeqs += seqs.size();
        uint16_t** freq = new uint16_t* [6];
        for (int i = 0; i < 6; ++i) {
            freq[i] = new uint16_t [seqLen];
            for (int j = 0; j < seqLen; ++j) freq[i][j] = 0;
        }
        for (auto sIdx = seqs.begin(); sIdx != seqs.end(); ++sIdx) {
            for (int j = 0; j < seqLen; ++j) {
                if      (sIdx->second[j] == 'A' || sIdx->second[j] == 'a') freq[0][j]+=1;
                else if (sIdx->second[j] == 'C' || sIdx->second[j] == 'c') freq[1][j]+=1;
                else if (sIdx->second[j] == 'G' || sIdx->second[j] == 'g') freq[2][j]+=1;
                else if (sIdx->second[j] == 'T' || sIdx->second[j] == 't' ||
                         sIdx->second[j] == 'U' || sIdx->second[j] == 'u') freq[3][j]+=1;
                else if (sIdx->second[j] == 'N' || sIdx->second[j] == 'n') freq[4][j]+=1;
                else                                                       freq[5][j]+=1;
            }
        }

        std::string subtreeFreqFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".final.freq.txt";
        std::ofstream outFile(subtreeFreqFile);
        if (!outFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", subtreeFreqFile.c_str());
            exit(1);
        }
        // Info subtreeIdx, seqNum, seqLen
        outFile << subtree << ',' << seqs.size() << ',' << seqLen << '\n';
        seqs.clear();
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < seqLen-1; ++j) {
                outFile << freq[i][j] << ',';
            }
            outFile << freq[i][seqLen-1] << '\n';
        }

        for (int i = 0; i < 6; ++i) delete [] freq[i];
        delete [] freq;
        outFile.close();
    }
    
    tree->root->msa.clear();
    return;
}

void printTree(Node* node, int grpID)
{
    if (grpID == -1) {
        if (node->parent == nullptr)
            std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << "ROOT\t" << node->numLeaves << '\t' << node->weight <<std::endl;
        else
            std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << "\t" << node->numLeaves << '\t' << node->weight << std::endl;  

        if (node->children.size() == 0) return;
        // std::cout << "Print children\n";
        for (auto &c: node->children) printTree(c, -1);
    }
    else {
        if (node->grpID == grpID) {
            if (node->parent == nullptr)
                std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << "ROOT\t" << node->grpID << std::endl;
            else
                std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << "\t" << node->grpID << std::endl;  
            if (node->children.size() == 0) return;
            // std::cout << "Print children\n";
        }
        for (auto &c: node->children) printTree(c, grpID);
    }
    
}

void printLeaves(Node* node)
{
    if (node->children.size() == 0) {
        std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << '\t' << node->grpID << std::endl;
        return;
    }
    for (auto &c: node->children) printLeaves(c);
}

Tree* readNewick(po::variables_map& vm)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::string treeFileName = vm["tree"].as<std::string>();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; inputStream >> newick;
    Tree *T = new Tree(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Newick string read in: " <<  treeBuiltTime.count() << " ns\n";

    // printTree(T->root);

    return T;
}

Tree* readNewick(std::string treeFileName)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; inputStream >> newick;
    Tree *T = new Tree(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Newick string read in: " <<  treeBuiltTime.count() << " ns\n";

    return T;
}

void readFreq(std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util) {
    for (auto subroot:  partition->partitionsRoot) {
        int subtree = tree->allNodes[subroot.first]->grpID;
        std::string freqFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".freq.txt";
        std::ifstream inputStream(freqFile);
        if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", freqFile.c_str()); exit(1); }
        std::vector<std::vector<int32_t>> freq;
        util->profileFreq[subtree] = freq;
        std::string rawInput;
        int idx; size_t seqNum, seqLen; 
        getline(inputStream, rawInput);
        std::string num = "";
        std::vector<int32_t> numbers;
        for (int i = 0; i < rawInput.size(); ++i) {
            if (rawInput[i] == ',') {
                numbers.push_back(std::atoi(num.c_str()));
                num = "";
            }
            else if (i == rawInput.size()-1) {
                num += rawInput[i];
                numbers.push_back(std::atoi(num.c_str()));
                num = "";
            }
            else num += rawInput[i];
        }
        assert(numbers.size() == 3);
        idx = numbers[0];
        seqNum = numbers[1];
        seqLen = numbers[2]; 
        assert(idx == subtree);
        numbers.clear();
        util->seqsLen[subroot.first] = seqLen;
        for (int i = 0; i < 6; ++i) {
            std::vector<int32_t> charFreq;
            util->profileFreq[subtree].push_back(charFreq);
            getline(inputStream, rawInput);
            for (int j = 0; j < rawInput.size(); ++j) {
                if (rawInput[j] == ',') {
                    util->profileFreq[subtree][i].push_back(std::atoi(num.c_str()));
                    num = "";
                }
                else if (j == rawInput.size()-1) {
                    num += rawInput[j];
                    util->profileFreq[subtree][i].push_back(std::atoi(num.c_str()));
                    num = "";
                }
                else num += rawInput[j];
            }
            if (util->profileFreq[subtree][i].size() > util->seqLen) util->seqLen = util->profileFreq[subtree][i].size();
        }
        inputStream.close();
        // std::remove(freqFile.c_str());
    }
    return;
}


void getFreq(Tree* tree, paritionInfo_t* partition, msa::utility* util) {
    for (auto subroot:  partition->partitionsRoot) {
        int subtree = tree->allNodes[subroot.first]->grpID;
        size_t seqLen = util->seqsLen[subroot.first]; 
        std::vector<std::vector<int32_t>> freq;
        util->profileFreq[subtree] = freq;
        for (int i = 0; i < 6; ++i) {
            std::vector<int32_t> temp;
            util->profileFreq[subtree].push_back(temp);
            for (int j = 0; j < seqLen; ++j) {
                util->profileFreq[subtree][i].push_back(0);
            }
        }
        
        for (auto sIdx: tree->allNodes[subroot.first]->msaIdx) { 
            tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
            int storage = util->seqsStorage[sIdx];
            for (int s = r.begin(); s < r.end(); ++s) {
            // for (int s = 0; s < refLen; ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') util->profileFreq[subtree][0][s]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') util->profileFreq[subtree][1][s]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') util->profileFreq[subtree][2][s]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') util->profileFreq[subtree][3][s]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') util->profileFreq[subtree][4][s]+=1;
                else                                                                                             util->profileFreq[subtree][5][s]+=1;
            }
            });
        }
    }
    // std::cout << "Finish getting frequency,\n";
    return;
}

void getPostOrderList(Node* node, std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(node); 
    Node* current; 
  
    while (!s1.empty()) {
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (auto ch: current->children) {
            if (ch->grpID == current->grpID) {
                s1.push(ch);
            }      
        }
    } 
    return;
}

bool cmp(std::string a, std::string b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
}

bool cmp_msaIdx(Node* a, Node* b) {
    return a->msaIdx.size() > b->msaIdx.size();
}

void outputFile(std::string fileName, msa::utility* util, Tree* T, int grpID) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    std::vector<std::string> seqs;
    for (auto seq: util->seqsIdx) {
        seqs.push_back(seq.first);
    }
    size_t seqLen = util->seqsLen[T->root->identifier];
    std::sort(seqs.begin(), seqs.end(), cmp);
    for (int s = 0; s < seqs.size(); ++s) {
        int sIdx = util->seqsIdx[seqs[s]];
        int storage = util->seqsStorage[sIdx];
        if (std::find(T->root->msaIdx.begin(), T->root->msaIdx.end(), sIdx) != T->root->msaIdx.end()) {
            outFile << '>' << seqs[s] << "\n";
            outFile.write(&util->alnStorage[storage][sIdx][0], seqLen);
            // int i = 0;
            // while (util->alnStorage[storage][sIdx][i] != 0) {
            //     outFile << util->alnStorage[storage][sIdx][i];
            //     ++i;
            // }
            outFile << '\n';
        }
        util->seqFree(sIdx);
    }
    outFile.close();
    util->seqsFree();
    return;
}

void outputFreq(std::string fileName, msa::utility* util, Tree* T, int grpID) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    // Info subtreeIdx, seqNum, seqLen
    outFile << grpID << ',' << T->root->msaIdx.size() << ',' << util->seqsLen[T->root->identifier] << '\n';

    size_t seqLen = util->seqsLen[T->root->identifier];
    std::cout << "seqLen: " << seqLen << '\n';
    uint32_t** freq = new uint32_t* [6];
    for (int i = 0; i < 6; ++i) {
        freq[i] = new uint32_t [seqLen];
        for (int j = 0; j <  seqLen; ++j) freq[i][j] = 0;
    }
    
    for (int sIdx: T->root->msaIdx) {
        int storage = util->seqsStorage[sIdx];
        for (int j = 0; j <  seqLen; ++j) {
            if      (util->alnStorage[storage][sIdx][j] == 'A' || util->alnStorage[storage][sIdx][j] == 'a') freq[0][j]+=1;
            else if (util->alnStorage[storage][sIdx][j] == 'C' || util->alnStorage[storage][sIdx][j] == 'c') freq[1][j]+=1;
            else if (util->alnStorage[storage][sIdx][j] == 'G' || util->alnStorage[storage][sIdx][j] == 'g') freq[2][j]+=1;
            else if (util->alnStorage[storage][sIdx][j] == 'T' || util->alnStorage[storage][sIdx][j] == 't' ||
                     util->alnStorage[storage][sIdx][j] == 'U' || util->alnStorage[storage][sIdx][j] == 'u') freq[3][j]+=1;
            else if (util->alnStorage[storage][sIdx][j] == 'N' || util->alnStorage[storage][sIdx][j] == 'n') freq[4][j]+=1;
            else                                                                                             freq[5][j]+=1;
        }
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < seqLen-1; ++j) {
            outFile << freq[i][j] << ',';
        }
        outFile << freq[i][seqLen-1] << '\n';
    }
    outFile.close();
    for (int i = 0; i < 6; ++i) delete [] freq[i];
    delete [] freq;
}

void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string> seqs) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        outFile << '>' << it->first << '\n';
        outFile << it->second << '\n';
    }

    outFile.close();
}

void getSubtreeNewick(Node* root, std::string& outputString) {
	if(root->children.size() != 0) {
		outputString += "(";
		for(int n = 0; n < root->children.size(); ++n) {
			if(n != 0) outputString += ",";
			getSubtreeNewick(root->children[n], outputString);
		}
		outputString += ")";
	}
	else {
		outputString += (root->identifier + ':' + std::to_string(root->branchLength));
    }
}

void outputSubtree(std::string fileName, Tree* T) {
	std::string out_str = "";
	getSubtreeNewick(T->root, out_str);
	out_str += ";\n";
	std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
	outFile << out_str;
	outFile.close();
}

void msaOnSubtree (Tree* T, msa::utility* util, msa::option* option, paritionInfo_t* partition, Params& param) {
    
    auto progressiveStart = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    for (auto &p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
        int grpID = p.second.first->grpID;
        getMsaHierachy(subhier, msaStack, grpID);
        
        for (auto h: subhier) {
            while (hier.size() < h.second+1) {
                std::vector<std::pair<Node*, Node*>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
    }


    // std::vector<std::pair<Node*, Node*>> debug;
    // for (int m = 0; m < hier.size(); ++m) {
    // // for (int m = 0; m < 2; ++m) {
    //     std::cout << m << '\n';
    //     for (int n = 0; n < hier[m].size(); ++n) {
    //         if (hier[m][n].first->identifier == "node_1943" || 
    //             hier[m][n].second->identifier == "node_1943" ) {
    //         std::cout << m << ',' << n << ':' << hier[m][n].first->identifier <<'(' <<hier[m][n].first->grpID<< ')' << ',' << hier[m][n].second->identifier<<'(' <<hier[m][n].second->grpID<< ')' << '\n';
    //         debug.push_back(std::make_pair(hier[m][n].first, hier[m][n].second));
    //         }
    //     }
    //     if (!debug.empty()) break;
    // }
    // Node* ch = T->allNodes[debug[0].first->identifier]->children[0];
    // std::cout << debug.size() << '\n';
    // std::string before = "", after = "";
    // std::string r = "", q = "";
    // int sIdx = util->seqsIdx[ch->identifier];
    // int storage = util->seqsStorage[sIdx];
    // // std::string r = "";
    // int j = 0;
    // while (util->alnStorage[storage][sIdx][j] != 0) {
    //     if (util->alnStorage[storage][sIdx][j] != '-') {
    //         before += util->alnStorage[storage][sIdx][j];
    //     }
    //     ++j;
    // }
    // std::cout << before.size() << '/' << util->seqsLen[ch->identifier] << '\n';
    // // std::cout << before << '\n';
    // // exit(1);
    // msaPostOrderTraversal_multigpu(T, debug, util, param);
    // j = 0;
    // storage = util->seqsStorage[sIdx];
    // while (util->alnStorage[storage][sIdx][j] != 0) {
    //     if (util->alnStorage[storage][sIdx][j] != '-') {
    //         after += util->alnStorage[storage][sIdx][j];
    //     }
    //     r += util->alnStorage[storage][sIdx][j];
    //     ++j;
    // }
    // sIdx = util->seqsIdx[debug[0].second->identifier];
    // storage = util->seqsStorage[sIdx];
    // j = 0;
    // while (util->alnStorage[storage][sIdx][j] != 0) {
    //     // if (util->alnStorage[storage][sIdx][j] != '-') {
    //         q += util->alnStorage[storage][sIdx][j];
    //     // }
    //     ++j;
    // }
    // std::cout << before.size() << '/' << after.size() << '\n';
    // for (int i = 0; i < before.size(); ++i) {
    //     if (before[i] != after[i]) std::cout << i << '\t' << before[i] << '\t' << after[i] << '\n';
    // }
    // std::cout << r.size() << ':' << q.size() << '\n';
    // std::cout << r << '\n' << q << '\n';
    // std::cout << util->seqsLen[debug[0].first->identifier] << ':' << util->seqsLen[debug[0].second->identifier] << '\n';
    // exit(1); 

    std::unordered_map<std::string, std::string> beforeAln;
    
    // for (auto s: T->allNodes) {
    //     if (s.second->is_leaf()) {
    //         std::string seqName = s.first;
    //         int sIdx = util->seqsIdx[s.first];
    //         int storage = util->seqsStorage[sIdx];
    //         std::string r = "";
    //         int j = 0;
    //         while (util->alnStorage[storage][sIdx][j] != 0) {
    //             if (util->alnStorage[storage][sIdx][j] != '-') {
    //                 r += util->alnStorage[storage][sIdx][j];
    //             }
    //             ++j;
    //         }
    //         beforeAln[seqName] = r;
    //     }
    // }
    // std::cout << "beforeAln Size: " << beforeAln.size() << '\n';
    int level = 0;
    for (auto m: hier) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        msaPostOrderTraversal_multigpu(T, m, util, option, param);
        // int err = 0;
        // for (auto s: beforeAln) {
        //     int sIdx = util->seqsIdx[s.first];
        //     int storage = util->seqsStorage[sIdx];
        //     std::string r = "", q = "";
        //     int offset = 0;
        //     while (util->alnStorage[storage][sIdx][offset] != 0) {
        //         if (util->alnStorage[storage][sIdx][offset] != '-') {
        //             r += util->alnStorage[storage][sIdx][offset];
        //         }
        //         q += util->alnStorage[storage][sIdx][offset];
        //         ++offset;
        //     }
        //     if (r != s.second) {
        //         err += 1;
        //         printf("seq: %s, the sequence did not match\n", s.first.c_str());
        //         if (r.size() != s.second.size()) {
        //             std::cout << "Wrong length. " << r.size() << '/' << s.second.size() << ".\n";
        //         }
        //         for (int i = 0; i < s.second.size(); ++i) {
        //             if (r[i] != s.second[i]) {
        //                 std::cout << "Mismatch at position " << i << '\n';
        //                 break;
        //             }
        //         }                
        //     }
        //     // if (s.first == "ON811175.1") std::cout << "ON811175.1: " << q.size() << '\n';
        //     if (err>0) exit(1);
        // }



        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
        else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
        ++level;
    }
    // Push msa results to roots of sub-subtrees
    for (auto p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<Node*> msaArray;
        while (!msaStack.empty()) {
            msaArray.push_back(msaStack.top());
            msaStack.pop();
        }
        if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
            if (msaArray.size() == 2) {
                T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
                util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
                break;
            }
            for (int m = msaArray.size()-2; m >=0; --m) {
                if (msaArray[m]->msaIdx.size()>0) {
                    T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
                    util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
                    break;
                }
            }
        }
    }
    auto progressiveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds progessiveTime = progressiveEnd - progressiveStart;
    std::cout << "Progressive alignment in " <<  progessiveTime.count() / 1000000000 << " s\n";
    
    if (util->badSequences.empty()) return;
    auto badStart = std::chrono::high_resolution_clock::now();
    std::cout << "Adding bad profiles back.\n";
    
    int levelThreshold = 0;
    int maxIteration = 2;
    int iteration = 0;
    bool lastIter = false;
    std::map<std::string, int> NodeAlnOrder;
    std::vector<std::pair<std::pair<Node*, Node*>, int>> alnOrder;
        
    
    while (!util->badSequences.empty() && iteration < maxIteration) {
        // if (iteration == maxIteration-1) lastIter = true;
        if (lastIter || iteration == maxIteration-1) util->nowProcess = 1;
        ++iteration;
        hier.clear();
        NodeAlnOrder.clear();
        alnOrder.clear();
        // for (auto bbbb: util->badSequences) {
        //     std::cout << "GrpID: " << bbbb.first  << " Size: " << bbbb.second.size() << '\n';
        //     for (int n = 0; n < bbbb.second.size(); ++n) std::cout << n << ':' << bbbb.second[n] << '\n';
        // }

        int badSeqBefore = 0, badSeqAfter = 0;
        int badProfileBefore = 0, badProfileAfter = 0;
        for (auto p: partition->partitionsRoot) {
            if (util->badSequences.find(p.second.first->grpID) != util->badSequences.end()) {
                auto badSeqName = util->badSequences[p.second.first->grpID];
                std::vector<Node*> badSeqNode;
                for (auto name: badSeqName) badSeqNode.push_back(T->allNodes[name]);
                // Sort bad profiles based on number of sequences
                std::sort(badSeqNode.begin(), badSeqNode.end(), cmp_msaIdx);
                
                if (badSeqNode.size() > levelThreshold && levelThreshold != 0) {
                    int tempSize = badSeqNode.size();
                    for (int i = 0; i < tempSize-levelThreshold; ++i) {
                        badSeqNode.pop_back();
                    }
                }
                
                badSeqName.clear();
                for (auto node: badSeqNode) badSeqName.push_back(node->identifier);
                
                badProfileBefore += badSeqNode.size();
                for (auto n: badSeqName) badSeqBefore += T->allNodes[n]->msaIdx.size();
                
                std::vector<std::string> nodeLeft;
                while (badSeqName.size() > 1) {
                    nodeLeft.clear();
                    for (int i = 0; i < badSeqName.size()-1; i+=2) {
                        int firstIdx  = (NodeAlnOrder.find(badSeqName[i]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i]]+1 : 0;
                        int secondIdx = (NodeAlnOrder.find(badSeqName[i+1]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[i+1]]+1 : 0;
                        int maxIdx = max(firstIdx, secondIdx);
                        NodeAlnOrder[badSeqName[i]] = maxIdx;
                        NodeAlnOrder[badSeqName[i+1]] = maxIdx;
                        alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[badSeqName[i]], T->allNodes[badSeqName[i+1]]), maxIdx));
                        nodeLeft.push_back(badSeqName[i]);
                    }
                    if (badSeqName.size()%2 == 1) nodeLeft.push_back(badSeqName.back());
                    badSeqName = nodeLeft;
                }
                assert(badSeqName.size() == 1);
                int idx  = (NodeAlnOrder.find(badSeqName[0]) != NodeAlnOrder.end()) ? NodeAlnOrder[badSeqName[0]]+1 : 0;
                alnOrder.push_back(std::make_pair(std::make_pair(T->allNodes[p.second.first->identifier], T->allNodes[badSeqName[0]]), idx));
            }
        }
        for (auto h: alnOrder) {
            while (hier.size() < h.second+1) {
                std::vector<std::pair<Node*, Node*>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
        util->badSequences.clear();
        level = 0;
        if (!hier.empty()) {
            std::cout << "Iteraton " << iteration-1 << ". Total " << hier.size() << " levels.\n";
            for (auto m: hier) {
                auto alnStart = std::chrono::high_resolution_clock::now();
                msaPostOrderTraversal_multigpu(T, m, util, option, param);
                auto alnEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                if (m.size() > 1) std::cout << "Level "<< level << ", aligned " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
                else              std::cout << "Level "<< level << ", aligned " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
                ++level;
            }
        }

        for (auto bad: util->badSequences) badProfileAfter += bad.second.size();
        for (auto bad: util->badSequences) for (auto name: bad.second) badSeqAfter += T->allNodes[name]->msaIdx.size();
        if (badSeqBefore == badSeqAfter) lastIter = true;
        std::cout << "== The number of Bad profiles/sequences ==\n";
        std::cout << "Before: " << std::setw(5) << badProfileBefore << " / " << badSeqBefore << '\n';
        std::cout << "After:  " << std::setw(5) << badProfileAfter << " / " << badSeqAfter << '\n';
        std::cout << "==========================================\n";
    }

    auto badEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds badTime = badEnd - badStart;
    std::cout << "Added bad profiles in " <<  badTime.count() / 1000000000 << " s\n";
    return;
}

void alignSubtrees (Tree* T, Tree* newT, msa::utility* util, msa::option* option, Params& param) {
    std::vector<std::pair<Node*, Node*>> type1Aln;
    for (auto n: newT->allNodes) {
        for (auto m: n.second->children) {
            if (newT->allNodes[m->identifier]->grpID == newT->allNodes[n.second->identifier]->grpID) {
                type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
            }
        }
    }
    // for (auto n: type1Aln) std::cout << n.first->identifier << ':' << util->seqsLen[n.first->identifier] << 
    //                              ',' << n.second->identifier << ':' << util->seqsLen[n.second->identifier] <<'\n';
    std::cout << "Align sub-subtrees, total " << type1Aln.size() << " pairs.\n"; 
    createOverlapMSA(T, type1Aln, util, option, param);
    return;
}

void mergeSubtrees (Tree* T, Tree* newT, msa::utility* util) {
    for (auto n: newT->allNodes) {
        T->allNodes[n.first]->msa.clear();
        T->allNodes[n.first]->msa.push_back(n.first);
    }
    
    std::vector<std::pair<Node*, Node*>> mergePairs;
    for (auto n: newT->allNodes) {
        if (n.second->children.size() > 1) {
            mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
            for (int i = 1; i < n.second->children.size(); ++i) {
                mergePairs.push_back(std::make_pair(n.second->children[0], n.second->children[i]));
            }
        }
        else if (n.second->children.size() == 1) {
            mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
        }
    }
    
    std::vector<std::pair<Node*, Node*>> singleLevel;
    std::map<std::string, char> addedNodes;
    
    int totalEdges = 0;
    int totalLevels = 0;
    while (true) {
        auto roundStart = std::chrono::high_resolution_clock::now();
        ++totalLevels;
        addedNodes.clear();
        singleLevel.clear();
        for (auto it = mergePairs.begin(); it != mergePairs.end();) {
            Node* a = it->first;
            Node* b = it->second;
            if ((a->parent != nullptr && b->parent != nullptr) &&  
                (addedNodes.find(a->identifier) == addedNodes.end() && addedNodes.find(b->identifier) == addedNodes.end())) {
                singleLevel.push_back(std::make_pair(a, b));
                for (auto id: T->allNodes[a->identifier]->msa) addedNodes[id] = 0;
                for (auto id: T->allNodes[b->identifier]->msa) addedNodes[id] = 0;
                mergePairs.erase(it);
            }
            else {
                ++it;
            }
        }
        bool breakLoop = false;
        if (singleLevel.empty()) {
            for (auto mp: mergePairs) {
                singleLevel.push_back(mp);
            }
            breakLoop = true;
        }
        transitivityMerge(T, newT, singleLevel, util);
        for (auto n: singleLevel) {
            std::vector<std::string> temp = T->allNodes[n.first->identifier]->msa;

            for (int r = 1; r < T->allNodes[n.first->identifier]->msa.size(); ++r) {
                std::string grpNode = T->allNodes[n.first->identifier]->msa[r];
                for (auto id: T->allNodes[n.second->identifier]->msa) T->allNodes[grpNode]->msa.push_back(id);
            }
            for (auto id: T->allNodes[n.second->identifier]->msa) T->allNodes[n.first->identifier]->msa.push_back(id);
            for (int r = 1; r < T->allNodes[n.second->identifier]->msa.size(); ++r) {
                std::string grpNode = T->allNodes[n.second->identifier]->msa[r];
                for (auto id: temp) T->allNodes[grpNode]->msa.push_back(id);
            }
            for (auto id: temp) T->allNodes[n.second->identifier]->msa.push_back(id);
        }
        auto roundEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds roundTime = roundEnd - roundStart;
        if (singleLevel.size() > 1) {
            std::cout << "Merged "<< singleLevel.size() << " edges in " << roundTime.count() / 1000000 << " ms\n";
        }
        else {
            std::cout << "Merged "<< singleLevel.size() << " edge in " << roundTime.count() / 1000000 << " ms\n";
        }
        totalEdges += singleLevel.size();
        if (breakLoop) break;
        if (totalLevels % 100 == 0) std::cout << "=============" << totalLevels << "=================\n";
    }
    std::cout << "Total Edges/Levels: " << totalEdges << '/' << totalLevels << '\n';
    return;
}


void createOverlapMSA(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;

    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 

    
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;

    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    float** hostGapOp = new float* [gpuNum];
    float** hostGapEx = new float* [gpuNum];
    float** hostGapCl = new float* [gpuNum];

    

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceGapCl = new float* [gpuNum];
  
    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] = (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2 * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapCl[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            
            
            cudaSetDevice(gn);
            // cudaSetDevice(1);
            // cudaError_t error;
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            cudaMalloc((void**)&deviceGapOp[gn], 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn], 2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapCl[gn], 2 * seqLen * numBlocks * sizeof(float));


            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            // error = cudaGetLastError(); printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            std::vector<std::pair<int, int>> seqIdx;
            
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEx[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCl[gn][n] = 0;
                // Calculate Frequency
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    
                    if (util->nowProcess < 2) {
                        float refWeight = 0.0; float qryWeight = 0.0;
                        // get sum of weights
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                            int storage = util->seqsStorage[sIdx];
                            std::string name = util->seqsName[sIdx];
                            float w = tree->allNodes[name]->weight / refWeight * refNum;
                            // float w = 1.0;
                            tbb::this_task_arena::isolate( [&]{
                            tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                            for (int s = r.begin(); s < r.end(); ++s) {
                            // for (int s = 0; s < refLen; ++s) {
                                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                                if (s > 0) {
                                    if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+s] += 1.0;
                                    if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+s] += 1.0;
                                    // if (util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+s] += 1.0;
                                    if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+s] += 1.0;
                                }

                            }
                            });
                            });
                        }
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                            int storage = util->seqsStorage[sIdx];
                            std::string name = util->seqsName[sIdx];
                            // float w = 1.0;
                            float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                            tbb::this_task_arena::isolate( [&]{
                            tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                            for (int s = r.begin(); s < r.end(); ++s) {
                            // for (int s = 0; s < qryLen; ++s) {
                                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                                else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                                if (s > 0) {
                                    if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+seqLen+s] += 1.0;
                                    if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+seqLen+s] += 1.0;
                                    // if (util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+seqLen+s] += 1.0;
                                    if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+seqLen+s] += 1.0;
                                }
                            }
                            });
                            });

                        }
                        for (int s = 0; s < refLen; ++s) { hostGapOp[gn][2*seqLen*n+s] /= refNum; hostGapEx[gn][2*seqLen*n+s] /= refNum; hostGapCl[gn][2*seqLen*n+s] /= refNum;} 
                        for (int s = 0; s < qryLen; ++s) { hostGapOp[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapEx[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapCl[gn][2*seqLen*n+seqLen+s] /= qryNum;}
                        // Clustalw's method
                        for (int s = 0; s < refLen; ++s) {
                            if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                                hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5])*1.0 / refNum);
                            }
                            else {
                                int backSites = 8;
                                int distance_from_gap = 0;
                                int increPenalty = false;
                                for (int d = s-1; d > max(s-backSites, 0); --d) {
                                    ++distance_from_gap;
                                    if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                        increPenalty = true;
                                        break;
                                    }
                                }
                                hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                            }
                        }
                        for (int s = 0; s < qryLen; ++s) {
                            if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                                hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                            }
                            else {
                                int backSites = 8;
                                int distance_from_gap = 0;
                                int increPenalty = false;
                                for (int d = s-1; d > max(s-backSites, 0); --d) {
                                    ++distance_from_gap;
                                    if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                        increPenalty = true;
                                        break;
                                    }
                                }
                                hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                            }
                        }

                    }
                    else {
                        int subtreeRef = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
                        for (int i = 0; i < 6; ++i) {
                            for (int s = 0; s < refLen; ++s) {
                                hostFreq[gn][12*seqLen*n+6*s+i] = util->profileFreq[subtreeRef][i][s]; 
                            }
                        }
                        int subtreeQry = tree->allNodes[nodes[nIdx].second->identifier]->grpID;
                        for (int i = 0; i < 6; ++i) {
                            for (int s = 0; s < qryLen; ++s) {
                                hostFreq[gn][12*seqLen*n+6*(seqLen+s)+i] = util->profileFreq[subtreeQry][i][s]; 
                            }
                        }
                    }
                    hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;    
                }


                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);


                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceGapCl[gn],
                    deviceParam[gn]
                );
                // cudaDeviceSynchronize();
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
                
                tbb::this_task_arena::isolate([&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    // if (nIdx % 100 == 0) {
                    if (hostAlnLen[gn][n] <= 0) {
                        std::vector<int8_t> aln;
                        std::vector<std::vector<float>> freqRef;
                        std::vector<std::vector<float>> freqQry;
                        std::vector<std::vector<float>> gapOp, gapEx, gapCl;
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            
                        for (int r = 0; r < refLen; r++) {
                            std::vector<float> temp;
                            for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*r+f]);
                            freqRef.push_back(temp);
                        }
                        for (int q = 0; q < qryLen; q++) {
                            std::vector<float> temp;
                            for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*(seqLen+q)+f]);
                            freqQry.push_back(temp);
                        }
                        for (int i = 0; i < 2; ++i) {
                            std::vector<float> temp;
                            gapOp.push_back(temp); gapEx.push_back(temp); gapCl.push_back(temp);
                        }
                        for (int r = 0; r < hostLen[gn][2*n]; ++r) {
                            gapOp[0].push_back(hostGapOp[gn][2*n*seqLen+r]);
                            gapEx[0].push_back(hostGapEx[gn][2*n*seqLen+r]);
                            gapCl[0].push_back(hostGapCl[gn][2*n*seqLen+r]);
                        }
                        for (int q = 0; q < hostLen[gn][2*n+1]; ++q) {
                            gapOp[1].push_back(hostGapOp[gn][2*n*seqLen+seqLen+q]);
                            gapEx[1].push_back(hostGapEx[gn][2*n*seqLen+seqLen+q]);
                            gapCl[1].push_back(hostGapCl[gn][2*n*seqLen+seqLen+q]);
                        }
                        // for (int s = 0; s < refLen; ++s) for (int j = 0; j < 6; ++j) freqRef[s][j] = hostFreq[gn][12*seqLen*n+6*s+j];
                        // for (int s = 0; s < qryLen; ++s) for (int j = 0; j < 6; ++j) freqQry[s][j] = hostFreq[gn][12*seqLen*n+6*(seqLen+s)+j];
                        Talco_xdrop::Params talco_params(hostParam);
                        while (aln.empty()) {
                            int16_t errorType = 0;
                            Talco_xdrop::Align_freq (
                                talco_params,
                                freqRef,
                                freqQry,
                                gapOp,
                                gapEx,
                                gapCl,
                                aln,
                                errorType
                            );
                            // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                            if (errorType == 1) {
                                std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                talco_params.updateXDrop(2*talco_params.xdrop);
                            }
                            if (errorType == 2) {
                                std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                talco_params.updateFLen(talco_params.fLen << 1);
                            }
                        }
        
                        // Talco_xdrop::Params talco_params(hostParam);
                        // int16_t errorType = 0;
                        // Talco_xdrop::Align_freq (
                        //     talco_params,
                        //     freqRef,
                        //     freqQry,
                        //     aln,
                        //     errorType
                        // );
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                        // std::cout << "CPU fallback on No. " << n << " (" << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ")\n";
                        printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
                    }
                    else {
                        std::vector<int8_t> aln;
                        for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            aln.push_back(hostAln[gn][n*2*seqLen+j]);
                        }
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                    }
                } 
                });   
                });
            }  
        }
    });
    
    // for (auto n: nodes) {
    //     tree->allNodes[n.first->identifier]->msa.clear();
    //     tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
    //     tree->allNodes[n.second->identifier]->msa.clear();
    //     tree->allNodes[n.second->identifier]->msa.push_back(n.second->identifier);    
    // }
    
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        // cudaSetDevice(1);
        cudaSetDevice(gn);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEx[gn]);
        cudaFree(deviceGapCl[gn]);
        cudaDeviceSynchronize();  
        free(hostFreq[gn]);
        free(hostAlnLen[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAln[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEx[gn]);
        free(hostGapCl[gn]);
    }
    
    free(hostParam);

    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceGapOp;
    delete [] deviceGapEx;
    delete [] deviceGapCl;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostGapCl;
    return;
}

void transitivityMerge(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util) {
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range) {
    for (int s = range.begin(); s < range.end(); ++s) {
        auto n = nodes[s];
    // for (auto n: nodes) {
        if (newtree->allNodes[n.first->identifier]->parent == nullptr) {
            
            for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
            std::vector<int8_t> rootAln;
            for (int i = 0; i < tree->allNodes[n.second->identifier]->msaAln.size(); ++i) {
                if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 2) rootAln.push_back(1);
                else if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 1) rootAln.push_back(2);
                else rootAln.push_back(tree->allNodes[n.second->identifier]->msaAln[i]);
            }
            tree->allNodes[n.first->identifier]->msaAln = rootAln;

            int seqLen = tree->allNodes[n.first->identifier]->msaAln.size();
            util->seqsLen[n.first->identifier] = seqLen;
            if (util->nowProcess < 2) {
                util->memCheck(seqLen);
                for (auto id: tree->allNodes[n.first->identifier]->msa) {
                    // auto id = tree->root->msa[k];
                    std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
                    for (auto sIdx: tree->allNodes[id]->msaIdx) {
                        int orgIdx = 0;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        for (int j = 0; j < aln.size(); ++j) {
                            if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                                util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->alnStorage[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[id] = aln.size();
                        util->changeStorage(sIdx);
                    }
                }
                // });
                for (auto id: tree->allNodes[n.first->identifier]->msa) {
                    if (id != tree->allNodes[n.first->identifier]->identifier) {
                        for (auto Idx: tree->allNodes[id]->msaIdx) {
                            tree->allNodes[n.first->identifier]->msaIdx.push_back(Idx);
                        }
                    }
                }
                tree->allNodes[n.first->identifier]->msa.clear();
                continue;
            }
            // else {
            //     util->seqLen = seqLen;
            //     if (util->seqNum != 0) {
            //         for (auto id: tree->allNodes[n.first->identifier]->msa) {
            //             std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
            //             for (auto sIdx: tree->allNodes[id]->msaIdx) {
            //                 char* seq = new char[seqLen+1];
            //                 int storage = util->seqsStorage[sIdx];
            //                 int orgIdx = 0;
            //                 for (int j = 0; j < seqLen+1; ++j) {
            //                     if (j < aln.size()) {
            //                         if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
            //                             seq[j] = util->seqs[sIdx][orgIdx];
            //                             orgIdx++;
            //                         }
            //                         else {
            //                             seq[j] = '-';
            //                         }
            //                     }
            //                     else {
            //                         seq[j] = 0;
            //                     }

            //                 }
            //                 util->seqsLen[id] = aln.size();
            //                 // util->changeStorage(sIdx);
            //                 delete [] util->seqs[sIdx];
            //                 util->seqs[sIdx] = seq;
            //             }
            //         }
            //         // });
            //         for (auto id: tree->allNodes[n.first->identifier]->msa) {
            //             if (id != tree->allNodes[n.first->identifier]->identifier) {
            //                 for (auto Idx: tree->allNodes[id]->msaIdx) {
            //                     tree->allNodes[n.first->identifier]->msaIdx.push_back(Idx);
            //                 }
            //             }
            //         }
            //         tree->allNodes[n.first->identifier]->msa.clear();
                
            //     }
            //     continue;
            // }
            
        }
        int8_t refGap = (newtree->allNodes[n.first->identifier]->level == newtree->allNodes[n.second->identifier]->level) ? 2 : 1; 
        std::vector<std::vector<int8_t>> refAln, qryAln;
        std::vector<std::vector<int8_t>> refNewAln, qryNewAln;
        // refAln.push_back(tree->allNodes[n.first->identifier]->msaAln);
        // qryAln.push_back(tree->allNodes[n.second->identifier]->msaAln);
        for (auto id: tree->allNodes[n.first->identifier]->msa)  refAln.push_back(tree->allNodes[id]->msaAln);
        for (auto id: tree->allNodes[n.second->identifier]->msa) qryAln.push_back(tree->allNodes[id]->msaAln);
        
        std::vector<int8_t> refOpAln = tree->allNodes[n.first->identifier]->msaAln;
        std::vector<int8_t> qryOpAln = tree->allNodes[n.second->identifier]->msaAln;
        int32_t refLen = refOpAln.size();
        int32_t qryLen = qryOpAln.size();
        int32_t refNum = refAln.size();
        int32_t qryNum = qryAln.size();
        int32_t seqLen = max(refLen, qryLen);
        for (int i = 0; i < refNum; ++i) {
            std::vector<int8_t> temp;
            refNewAln.push_back(temp);
        }
        for (int i = 0; i < qryNum; ++i) {
            std::vector<int8_t> temp;
            qryNewAln.push_back(temp);
        }
        // std::cout << refLen << ':';
        // for (int i = 0; i < refLen; ++i) {
        //     if ((refOpAln[i] & 0XFFFF) == refGap || (refOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // std::cout << qryLen << ':';
        // for (int i = 0; i < qryLen; ++i) {
        //     if ((qryOpAln[i] & 0XFFFF) == 2 || (qryOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // aln = 0: match
        // aln = 1: ref gap
        // aln = 2: qry gap
        // aln = 3: permenant gap
        int32_t rIdx = 0, qIdx = 0;
        while (rIdx < refLen && qIdx < qryLen) {
            if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                qIdx++;rIdx++;
            }
            else if (((refOpAln[rIdx] & 0xFFFF) == refGap || (refOpAln[rIdx] & 0xFFFF) == 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                int consecGap = 0;
                int k = rIdx;
                while (((refOpAln[k] & 0xFFFF) == refGap || (refOpAln[k] & 0xFFFF) == 3) && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
            }
            else if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) == 2 || (qryOpAln[qIdx] & 0xFFFF) == 3)) {
                int consecGap = 0;
                int k = qIdx;
                while (((qryOpAln[k] & 0xFFFF) == 2 || (qryOpAln[k] & 0xFFFF) == 3) && k < qryLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (((refOpAln[kr] & 0xFFFF) == refGap || (refOpAln[kr] & 0xFFFF) == 3) && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
                consecGap = 0;
                while (((qryOpAln[kq] & 0xFFFF) == 2 || (qryOpAln[kq] & 0xFFFF) == 3) && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
        }
        // printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][g]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);    
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][g]);
            }
        }
        assert (refNewAln[0].size() == qryNewAln[0].size());
        for (int i = 0; i < tree->allNodes[n.first->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.first->identifier]->msa[i];
            tree->allNodes[id]->msaAln = refNewAln[i];
        } 
        for (int i = 0; i < tree->allNodes[n.second->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.second->identifier]->msa[i];
            tree->allNodes[id]->msaAln = qryNewAln[i];
        } 
        
        // std::vector<std::string> temp = tree->allNodes[n.first->identifier]->msa;

        // for (int r = 1; r < tree->allNodes[n.first->identifier]->msa.size(); ++r) {
        //     std::string grpNode = tree->allNodes[n.first->identifier]->msa[r];
        //     for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[grpNode]->msa.push_back(id);
        // }
        // for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
        // for (int r = 1; r < tree->allNodes[n.second->identifier]->msa.size(); ++r) {
        //     std::string grpNode = tree->allNodes[n.second->identifier]->msa[r];
        //     for (auto id: temp) tree->allNodes[grpNode]->msa.push_back(id);
        // }
        // for (auto id: temp) tree->allNodes[n.second->identifier]->msa.push_back(id);
    }
    });
    
    return;
}

/*
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, Params& param)
{

    // auto freqStart = std::chrono::high_resolution_clock::now();
    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msaIdx.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                        util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                        break;
                    }
                }
            }
        }
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msaIdx.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                        util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                        break;
                    }
                }
            }
        }
    }

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = util->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
    
    paramType* hostParam = (paramType*)malloc(28 * sizeof(paramType)); 

    if (!param.userDefine) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 5 || j == 5)          hostParam[i*5+j] = 0;
                else if (i == j)               hostParam[i*5+j] = param.match;
                else if (i-j == 2 || j-i == 2) hostParam[i*5+j] = param.trans;
                else                           hostParam[i*5+j] = param.mismatch;
            }
        }
        hostParam[25] = param.gapOpen;
        hostParam[26] = param.gapExtend;
        hostParam[27] = param.xdrop;
    }
    else {
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.userMatrix[i][j];
        hostParam[25] = param.userGapOpen;
        hostParam[26] = param.userGapExtend;
        hostParam[27] = param.xdrop;
    }
    

    std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    // allocate memory on host and device
    uint16_t** hostFreq = new uint16_t* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];

    uint16_t** deviceFreq = new uint16_t* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    paramType**  deviceParam = new paramType* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
    // nowRound.store(roundGPU-1);


    int maxThreads = tbb::this_task_arena::max_concurrency();
    
    // int ThreadsPerGPU = maxThreads / gpuNum;
    bool* cpuFallback = new bool[nodes.size()];
    for (int i = 0; i < nodes.size(); ++i) cpuFallback[i] = false;

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] = (uint16_t*)malloc(12 * seqLen * numBlocks * sizeof(uint16_t));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(5                     * sizeof(int32_t));
            
            cudaSetDevice(gn);
            // cudaSetDevice(1);
            // cudaError_t error;
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(uint16_t));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 5 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn],  28 * sizeof(paramType));

            cudaMemcpy(deviceParam[gn], hostParam, 28 * sizeof(paramType), cudaMemcpyHostToDevice);
            // error = cudaGetLastError(); printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            std::vector<std::pair<int, int>> seqIdx;
            
            
            while (nowRound < roundGPU) {
            // while (nowRound >= 0) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                int seqNum = 0;
                // std::cout << "GPU: " << gn << " Rn: " << rn << " Pairs: " << alnPairs << '\n';

                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                seqIdx.clear();

                // Calculate Frequency
                // tbb::task_arena ta {maxThreads - gpuNum};
                // tbb::task_group tg;
                // tg.run ([&] () {
                // tbb::this_task_arena::isolate ([&] () {
                // std::cout << "SSS\n";
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> aln_range) { 
                // for (int n = aln_range.begin(); n < aln_range.end(); ++n) {
                
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t qryIdx = 0;
                    int32_t refIdx = 0;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    refIdx = seqNum;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        int maxLen = max(refLen, qryLen);
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < refLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1;
                            else                                                                                             hostFreq[gn][12*seqLen*n+6*s+5]+=1;
                        }
                        });
                        });
                        seqNum += 1;
                    }
                    qryIdx = seqNum;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                    int storage = util->seqsStorage[sIdx];
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                    for (int s = r.begin(); s < r.end(); ++s) {
                    // for (int s = 0; s < qryLen; ++s) {
                        if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1;
                        else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1;
                        else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1;
                        else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                 util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1;
                        else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1;
                        else                                                                                             hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1;
                    }
                    });
                    });
                    seqNum += 1;
                }
                    hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    seqIdx.push_back(std::make_pair(refIdx, qryIdx));
                }
                
                hostSeqInfo[gn][0] = seqLen;
                hostSeqInfo[gn][1] = seqNum;
                hostSeqInfo[gn][2] = alnPairs;
                hostSeqInfo[gn][3] = numBlocks;
                hostSeqInfo[gn][4] = param.userDefine;
        
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(uint16_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 5 * sizeof(int32_t), cudaMemcpyHostToDevice);
                
                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceParam[gn]
                );
                cudaDeviceSynchronize();
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
                int maxAlnLen = 0;
                for (int n = 0; n <  alnPairs; ++n) {
                    if (hostAlnLen[gn][n] > maxAlnLen) maxAlnLen = hostAlnLen[gn][n];
                }
                util->memCheck(maxAlnLen);
                if (rn % 10 == 0 && rn > 0) std::cout << rn*numBlocks << " pairs have been processed.\n";
                
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    
                    int32_t refNum = seqIdx[n].second - seqIdx[n].first;
                    int32_t qryNum = (n !=  alnPairs-1) ? seqIdx[n+1].first - seqIdx[n].second : seqNum - seqIdx[n].second;
                    int32_t nIdx = n + rn*numBlocks;

                    // if (nIdx % 400 == 399) {
                    if (hostAlnLen[gn][n] <= 0) {
                        cpuFallback[nIdx] = true;
                        // int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        // int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        // uint16_t *freq = new uint16_t[12*seqLen]; 
                        // for (int i = 0; i < 12*seqLen; ++i) freq[i] = hostFreq[gn][12*seqLen*n+i];
                        // std::vector<int8_t> aln;
                        // alignGrpToGrp_traditional (
                        //     freq,
                        //     seqLen,
                        //     refLen,
                        //     qryLen,
                        //     param,
                        //     aln
                        // );
                        // delete [] freq;
                        // int32_t alnLen = aln.size();
                        // util->memCheck(alnLen);
                        // std::reverse(aln.begin(), aln.end());
                        // for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                        //     int storeFrom = util->seqsStorage[sIdx];
                        //     int storeTo = 1 - util->seqsStorage[sIdx];
                        //     int orgIdx = 0;
                        //     for (int j = 0; j < aln.size(); ++j) {
                        //         if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                        //             util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        //             orgIdx++;
                        //         }
                        //         else {
                        //             util->alnStorage[storeTo][sIdx][j] = '-';
                        //         }
                        //     }
                        //     util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                        //     util->changeStorage(sIdx);
                        // }
                        // for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                        //     int storeFrom = util->seqsStorage[sIdx];
                        //     int storeTo = 1 - util->seqsStorage[sIdx];
                        //     int orgIdx = 0;
                        //     for (int j = 0; j < aln.size(); ++j) {
                        //         if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                        //             util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        //             orgIdx++;
                        //         }
                        //         else {
                        //             util->alnStorage[storeTo][sIdx][j] = '-';
                        //         }
                        //     }
                        //     util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                        //     util->changeStorage(sIdx);
                        // }
                        // std::cout << "CPU fallback (traditional global alignment) on No. " << n << " (" << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ")\n";
                    }
                    else {
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                            int orgIdx = 0;
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 2) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                            }
                            util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][n];
                            util->changeStorage(sIdx);
                        }
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            int orgIdx = 0;
                            for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                            }
                            util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
                            util->changeStorage(sIdx);
                        }
                        for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                        }
                        tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                    }
                    
                }
                });
                });
                
            }  

            
        }
    });
    
    // free memory  
    free(hostParam);
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(gn);
        // cudaSetDevice(1);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaDeviceSynchronize();  
        free(hostFreq[gn]);
        free(hostAln[gn]);
        free(hostLen[gn]);
        free(hostAlnLen[gn]);
        free(hostSeqInfo[gn]);
    }
    
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceLen;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostAln;
    delete [] hostSeqInfo;
    
    
    // CPU Fallback
    std::vector<int> fallbackPairs;
    for (int i = 0; i < nodes.size(); ++i) if (cpuFallback[i]) fallbackPairs.push_back(i);
    delete [] cpuFallback;
    if (fallbackPairs.size() > 0) std::cout << "CPU Fallback. Num of pairs: " << fallbackPairs.size() << '\n';
    else return;
    tbb::parallel_for(tbb::blocked_range<int>(0, fallbackPairs.size()), [&](tbb::blocked_range<int> range) {
    for (int n = range.begin(); n < range.end(); ++n) {
                
    // for (auto nIdx: fallbackPairs) {
        int nIdx = fallbackPairs[n];
        uint16_t *freq = new uint16_t[12*seqLen]; 
        for (int i = 0; i < 12*seqLen; ++i) freq[i] = 0;
        // int32_t qryIdx = 0;
        // int32_t refIdx = 0;
        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            int maxLen = max(refLen, qryLen);
            // tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
            // for (int s = r.begin(); s < r.end(); ++s) {
            for (int s = 0; s < refLen; ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freq[6*s+0]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freq[6*s+1]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freq[6*s+2]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freq[6*s+3]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') freq[6*s+4]+=1;
                else                                                                                             freq[6*s+5]+=1;
            }
            // });
        }
        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            // tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
            // for (int s = r.begin(); s < r.end(); ++s) {
            for (int s = 0; s < qryLen; ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freq[6*(seqLen+s)+0]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freq[6*(seqLen+s)+1]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freq[6*(seqLen+s)+2]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freq[6*(seqLen+s)+3]+=1;
                else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') freq[6*(seqLen+s)+4]+=1;
                else                                                                                             freq[6*(seqLen+s)+5]+=1;
            }
            // });
        }
        std::vector<int8_t> aln;
        alignGrpToGrp_traditional (
            freq,
            seqLen,
            refLen,
            qryLen,
            param,
            aln
        );
        delete [] freq;
        int32_t alnLen = aln.size();
        util->memCheck(alnLen);
        std::reverse(aln.begin(), aln.end());
        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            int orgIdx = 0;
            for (int j = 0; j < aln.size(); ++j) {
                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                    orgIdx++;
                }
                else {
                    util->alnStorage[storeTo][sIdx][j] = '-';
                }
            }
            util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
            util->changeStorage(sIdx);
        }
        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            int orgIdx = 0;
            for (int j = 0; j < aln.size(); ++j) {
                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                    orgIdx++;
                }
                else {
                    util->alnStorage[storeTo][sIdx][j] = '-';
                }
            }
            util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
            util->changeStorage(sIdx);
        }
        for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
        }
        tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
        printf("CPU fallback (traditional global alignment) on No. %d (%s), Alignment Length: %d\n", nIdx, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
    }
    });
    
    
    
    return;
}
*/

/*
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{

    
    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        // is leaf
        if (n.first->is_leaf() && n.first->msaIdx.empty()) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.first->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->msaIdx.empty()) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.second->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        // if (nodes.size() < 5) {
        //     std::cout << n.first->identifier << '(' << util->seqsLen[n.first->identifier] << ')' <<  tree->allNodes[n.first->identifier]->msaIdx.size() << '\t';
        //     std::cout << n.second->identifier << '(' << util->seqsLen[n.second->identifier] << ')'<< tree->allNodes[n.second->identifier]->msaIdx.size() << '\n';
        // }
    }

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
            
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
            
    
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
    
            

    // std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    // allocate memory on host and device
    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    
    paramType** hostGapOp = new paramType* [gpuNum]; // gap open
    paramType** hostGapCt = new paramType* [gpuNum]; // gap continous
    paramType** hostGapEn = new paramType* [gpuNum]; // gap end

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    paramType** deviceGapOp = new paramType* [gpuNum];
    paramType** deviceGapCt = new paramType* [gpuNum];
    paramType** deviceGapEn = new paramType* [gpuNum];
    
    paramType**  deviceParam = new paramType* [gpuNum];

    std::atomic<int> nowRound;
    nowRound.store(0);
            
    // nowRound.store(roundGPU-1);

    // int ThreadsPerGPU = maxThreads / gpuNum;
    bool* cpuFallback = new bool[nodes.size()];
    std::vector<std::vector<int8_t>> alnCpu;
    for (int i = 0; i < nodes.size(); ++i) {
        cpuFallback[i] = false;
        if (util->nowProcess == 1) {
            std::vector<int8_t> aln;
            alnCpu.push_back(aln);
        }
        
    }
    
            

    // tbb::task_arena limited_arena(2);
    // limited_arena.execute([&]{

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] =  (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2                     * sizeof(int32_t));
            hostGapOp[gn] = (paramType*)malloc(2 * seqLen * numBlocks * sizeof(paramType));
            hostGapCt[gn] = (paramType*)malloc(2 * seqLen * numBlocks * sizeof(paramType));
            hostGapEn[gn] = (paramType*)malloc(2 * seqLen * numBlocks * sizeof(paramType));
            cudaSetDevice(option->gpuIdx[gn]);
            // cudaSetDevice(1);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],  2 * seqLen * numBlocks * sizeof(paramType));
            cudaMalloc((void**)&deviceGapCt[gn],  2 * seqLen * numBlocks * sizeof(paramType));
            cudaMalloc((void**)&deviceGapEn[gn],  2 * seqLen * numBlocks * sizeof(paramType));

            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            
            while (nowRound < roundGPU) {
            // while (nowRound >= 0) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                int seqNum = 0;
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCt[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEn[gn][n] = 0;

                // Calculate Frequency
                // tbb::task_arena ta {maxThreads - gpuNum};
                // tbb::task_group tg;
                // tg.run ([&] () {
                // tbb::this_task_arena::isolate ([&] () {
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> aln_range) { 
                // for (int n = aln_range.begin(); n < aln_range.end(); ++n) {
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    float refWeight = 0.0; float qryWeight = 0.0;
                    // get sum of weights
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        float w = tree->allNodes[name]->weight / refWeight * refNum;
                        // float w = 1.0;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < refLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                            if (s > 0) {
                                if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+s] += 1.0;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapEn[gn][2*seqLen*n+s] += 1.0;
                            }
                            
                        }
                        });
                        });
                        seqNum += 1;
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        // float w = 1.0;
                        float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < qryLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                            if (s > 0) {
                                if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+seqLen+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+seqLen+s] += 1.0;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+seqLen+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapEn[gn][2*seqLen*n+seqLen+s] += 1.0;
                            }
                        }
                        });
                        });
                        seqNum += 1;
                    }
                    hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                    for (int s = 0; s < refLen; ++s) { hostGapOp[gn][2*seqLen*n+s] /= refNum; hostGapCt[gn][2*seqLen*n+s] /= refNum; hostGapEn[gn][2*seqLen*n+s] /= refNum;} 
                    for (int s = 0; s < qryLen; ++s) { hostGapOp[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapCt[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapEn[gn][2*seqLen*n+seqLen+s] /= qryNum;}
                    
                    // ClustalW's method
                    
                    for (int s = 0; s < refLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5]) * 1.0 / refNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    for (int s = 0; s < qryLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    
                }
                for (int i = 0; i < 50; ++i) std::cout << hostGapOp[gn][i] << ',';
                std::cout << '\n';
                for (int i = 0; i < 50; ++i) std::cout << hostGapCt[gn][i] << ',';
                std::cout << '\n';
                for (int i = 0; i < 50; ++i) std::cout << hostGapEn[gn][i] << ',';
                std::cout << '\n';
                // for (int i = 0; i < 100; ++i) {
                //     float sumFreq = 0;
                //     for (int j = 0; j < 6; ++j) sumFreq += hostFreq[gn][6*i+j]; 
                //     std::cout << sumFreq << ',';
                // }
                // std::cout << '\n';
                
                
                // if (nodes[0].first->identifier == "node_1") {
                // std::ofstream outFile("ref.freq");
                // outFile << tree->allNodes[nodes[0].first->identifier]->grpID << ',' << tree->allNodes[nodes[0].first->identifier]->msaIdx.size() << ',' << util->seqsLen[nodes[0].first->identifier] << '\n';
                // uint32_t** freq = new uint32_t* [6];
                // int seq_len = util->seqsLen[nodes[0].first->identifier];
                // for (int i = 0; i < 6; ++i) {
                //     freq[i] = new uint32_t [seq_len];
                //     for (int j = 0; j <  seq_len; ++j) freq[i][j] = 0;
                // }
                // for (auto sIdx: tree->allNodes[nodes[0].first->identifier]->msaIdx) {
                //     int storage = util->seqsStorage[sIdx];
                //     for (int j = 0; j <  seq_len; ++j) {
                //         if      (util->alnStorage[storage][sIdx][j] == 'A' || util->alnStorage[storage][sIdx][j] == 'a') freq[0][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'C' || util->alnStorage[storage][sIdx][j] == 'c') freq[1][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'G' || util->alnStorage[storage][sIdx][j] == 'g') freq[2][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'T' || util->alnStorage[storage][sIdx][j] == 't' ||
                //                  util->alnStorage[storage][sIdx][j] == 'U' || util->alnStorage[storage][sIdx][j] == 'u') freq[3][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'N' || util->alnStorage[storage][sIdx][j] == 'n') freq[4][j]+=1;
                //         else                                                               freq[5][j]+=1;
                //     }
                // }
                // for (int i = 0; i < 6; ++i) {
                //     for (int j = 0; j < seq_len-1; ++j) {
                //         outFile << freq[i][j] << ',';
                //     }
                //     outFile << freq[i][seq_len-1] << '\n';
                // }
                // outFile.close();
                // for (int i = 0; i < 6; ++i) delete [] freq[i];
                // delete [] freq;
                // std::ofstream outFile2("qry.freq");
                // outFile2 << tree->allNodes[nodes[0].second->identifier]->grpID << ',' << tree->allNodes[nodes[0].second->identifier]->msaIdx.size() << ',' << util->seqsLen[nodes[0].second->identifier] << '\n';
                // freq = new uint32_t* [6];
                // seq_len = util->seqsLen[nodes[0].second->identifier];
                // for (int i = 0; i < 6; ++i) {
                //     freq[i] = new uint32_t [seq_len];
                //     for (int j = 0; j <  seq_len; ++j) freq[i][j] = 0;
                // }
                // for (auto sIdx: tree->allNodes[nodes[0].second->identifier]->msaIdx) {
                //     int storage = util->seqsStorage[sIdx];
                //     for (int j = 0; j <  seq_len; ++j) {
                //         if      (util->alnStorage[storage][sIdx][j] == 'A' || util->alnStorage[storage][sIdx][j] == 'a') freq[0][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'C' || util->alnStorage[storage][sIdx][j] == 'c') freq[1][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'G' || util->alnStorage[storage][sIdx][j] == 'g') freq[2][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'T' || util->alnStorage[storage][sIdx][j] == 't' ||
                //                  util->alnStorage[storage][sIdx][j] == 'U' || util->alnStorage[storage][sIdx][j] == 'u') freq[3][j]+=1;
                //         else if (util->alnStorage[storage][sIdx][j] == 'N' || util->alnStorage[storage][sIdx][j] == 'n') freq[4][j]+=1;
                //         else                                                               freq[5][j]+=1;
                //     }
                // }
                // for (int i = 0; i < 6; ++i) {
                //     for (int j = 0; j < seq_len-1; ++j) {
                //         outFile2 << freq[i][j] << ',';
                //     }
                //     outFile2 << freq[i][seq_len-1] << '\n';
                // }
                // outFile2.close();
                // for (int i = 0; i < 6; ++i) delete [] freq[i];
                // delete [] freq;
                // }
                
                
                // if (nodes[0].first->identifier == "node_1") {
                // std::string freqFile = "ref.freq";
                // std::ifstream inputStream(freqFile);
                // if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", freqFile.c_str()); exit(1); }
                // std::vector<std::vector<uint16_t>> freqRef;
                // std::string rawInput;
                // getline(inputStream, rawInput);
                // std::string num = "";
                // std::vector<uint16_t> numbers;
                // for (int i = 0; i < rawInput.size(); ++i) {
                //     if (rawInput[i] == ',') {
                //         numbers.push_back(std::atoi(num.c_str()));
                //         num = "";
                //     }
                //     else if (i == rawInput.size()-1) {
                //         num += rawInput[i];
                //         numbers.push_back(std::atoi(num.c_str()));
                //         num = "";
                //     }
                //     else num += rawInput[i];
                // }
                // assert(numbers.size() == 3);
                // int refLen = numbers.back();
                // hostLen[gn][2*0] = refLen;
                // numbers.clear();
                // for (int i = 0; i < 6; ++i) {
                //     std::vector<uint16_t> charFreq;
                //     freqRef.push_back(charFreq);
                //     getline(inputStream, rawInput);
                //     for (int j = 0; j < rawInput.size(); ++j) {
                //         if (rawInput[j] == ',') {
                //             freqRef[i].push_back(std::atoi(num.c_str()));
                //             num = "";
                //         }
                //         else if (j == rawInput.size()-1) {
                //             num += rawInput[j];
                //             freqRef[i].push_back(std::atoi(num.c_str()));
                //             num = "";
                //         }
                //         else num += rawInput[j];
                //     }
                // }
                // inputStream.close();
                // // std::remove(freqFile.c_str());
                // for (int i = 0; i < 6; ++i) {
                //     for (int s = 0; s < refLen; ++s) {
                //         hostFreq[gn][12*seqLen*0+6*s+i] = freqRef[i][s]; 
                //     }
                // }
                // freqFile = "qry.freq";
                // std::ifstream inputStream2(freqFile);
                // if (!inputStream2) { fprintf(stderr, "Error: Can't open file: %s\n", freqFile.c_str()); exit(1); }
                // std::vector<std::vector<uint16_t>> freqQry;
                // getline(inputStream2, rawInput);
                // num = "";
                // numbers.clear();
                // for (int i = 0; i < rawInput.size(); ++i) {
                //     if (rawInput[i] == ',') {
                //         numbers.push_back(std::atoi(num.c_str()));
                //         num = "";
                //     }
                //     else if (i == rawInput.size()-1) {
                //         num += rawInput[i];
                //         numbers.push_back(std::atoi(num.c_str()));
                //         num = "";
                //     }
                //     else num += rawInput[i];
                // }
                // assert(numbers.size() == 3);
                // int qryLen = numbers.back();
                // hostLen[gn][2*0+1] = qryLen;
                // numbers.clear();
                // for (int i = 0; i < 6; ++i) {
                //     std::vector<uint16_t> charFreq;
                //     freqQry.push_back(charFreq);
                //     getline(inputStream2, rawInput);
                //     for (int j = 0; j < rawInput.size(); ++j) {
                //         if (rawInput[j] == ',') {
                //             freqQry[i].push_back(std::atoi(num.c_str()));
                //             num = "";
                //         }
                //         else if (j == rawInput.size()-1) {
                //             num += rawInput[j];
                //             freqQry[i].push_back(std::atoi(num.c_str()));
                //             num = "";
                //         }
                //         else num += rawInput[j];
                //     }
                // }
                // inputStream2.close();
                // // std::remove(freqFile.c_str());
                // for (int i = 0; i < 6; ++i) {
                //     for (int s = 0; s < qryLen; ++s) {
                //         hostFreq[gn][12*seqLen*0+6*(seqLen+s)+i] = freqQry[i][s]; 
                //     }
                // }

                // }
                

                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
        
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(paramType), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEn[gn], hostGapEn[gn], 2 * seqLen * numBlocks * sizeof(paramType), cudaMemcpyHostToDevice);

                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapCt[gn],
                    deviceGapEn[gn],
                    deviceParam[gn]
                );
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();

                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                int maxAlnLen = 0;
                for (int n = 0; n <  alnPairs; ++n) {
                    if (hostAlnLen[gn][n] > maxAlnLen) maxAlnLen = hostAlnLen[gn][n];
                    // std::cout << n << ':' <<  hostAlnLen[gn][n] << '\n'; 
                }
                util->memCheck(maxAlnLen);
                if (rn % 10 == 0 && rn > 0) std::cout << rn*numBlocks << " pairs have been processed.\n";
                
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                for (int n = range.begin(); n < range.end(); ++n) {
                // for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    
                    
                    // if (nIdx % 1000 == 0 || hostAlnLen[gn][n] <= 0) {
                    if (hostAlnLen[gn][n] <= 0) {
                        cpuFallback[nIdx] = true;
                        if (util->nowProcess == 1) {
                            std::vector<int8_t> aln;
                            std::vector<std::vector<int32_t>> freqRef;
                            std::vector<std::vector<int32_t>> freqQry;
                            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                            for (int r = 0; r < refLen; r++) {
                                std::vector<int32_t> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*r+f]);
                                freqRef.push_back(temp);
                            }
                            for (int q = 0; q < qryLen; q++) {
                                std::vector<int32_t> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*(seqLen+q)+f]);
                                freqQry.push_back(temp);
                            }

                            Talco_xdrop::Params talco_params(hostParam);
                            while (aln.empty()) {
                                int16_t errorType = 0;
                                Talco_xdrop::Align_freq (
                                    talco_params,
                                    freqRef,
                                    freqQry,
                                    aln,
                                    errorType
                                );
                                // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                                if (errorType == 1) {
                                    std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                    talco_params.updateXDrop(2*talco_params.xdrop);
                                }
                                if (errorType == 2) {
                                    std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                    talco_params.updateFLen(talco_params.fLen << 1);
                                }
                            }
                            alnCpu[nIdx] = aln;
                        }
                        
                    }
                    else {
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                            int orgIdx = 0;
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 2) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                            }
                            // std::cout << "orgIdx: " << orgIdx << '\n';
                            util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][n];
                            util->changeStorage(sIdx);
                        }
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            int orgIdx = 0;
                            for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                    
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                                
                            }
                            util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
                            util->changeStorage(sIdx);
                        }
                        for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                        }
                        tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                    }
                    
                }
                });
                });
                // cudaFree(deviceFreq[gn]);
                // cudaFree(deviceAln[gn]);
                // cudaFree(deviceLen[gn]);
                // cudaFree(deviceAlnLen[gn]);
                // cudaFree(deviceSeqInfo[gn]);
                // cudaFree(deviceParam[gn]);
                // cudaDeviceSynchronize();  
                // free(hostFreq[gn]);
                // free(hostAln[gn]);
                // free(hostLen[gn]);
                // free(hostAlnLen[gn]);
                // free(hostSeqInfo[gn]);
            }  

            
        }
    });
    // });
    
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(option->gpuIdx[gn]);
        // cudaSetDevice(gn);
        // cudaSetDevice(1);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEn[gn]);
        cudaDeviceSynchronize();  
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
        if (freeErr != "no error") {
            printf("ERROR: Free memory Last%s!\n", freeErr.c_str());
            exit(1);
        }
        free(hostFreq[gn]);
        free(hostAln[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAlnLen[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEn[gn]);
    }
    
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceGapOp;
    delete [] deviceGapEn;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEn;
    
    
    // CPU Fallback
    std::vector<int> fallbackPairs;
    for (int i = 0; i < nodes.size(); ++i) if (cpuFallback[i]) fallbackPairs.push_back(i);
    delete [] cpuFallback;
    
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    if (util->nowProcess == 0) {
        std::cout << "Bad alignments. Num of pairs: " << fallbackPairs.size() << '\n';
        for (int i = 0; i < fallbackPairs.size(); ++i) {
            int nIdx = fallbackPairs[i];
            int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
            if (util->badSequences.find(grpID) == util->badSequences.end()) {
                std::vector<std::string> temp;
                util->badSequences[grpID] = temp;
            }
            if (refNum < qryNum) {
                // proceed with qry and store ref as bad sequences
                int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                util->seqsLen[nodes[nIdx].second->identifier] = refLen;
                util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
                auto temp = nodes[nIdx].second->msaIdx;
                tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << refNum <<" seqeuences).\n";
            }
            else {
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << qryNum <<" seqeuences).\n";
            }
        }
    }
    else {
        std::cout << "CPU Fallback. Num of pairs: " << fallbackPairs.size() << '\n';
        // std::vector<std::vector<int8_t>> alns;
        // for (int i = 0; i < fallbackPairs.size(); ++i) {
        //     std::vector<int8_t> aln;
        //     alns.push_back(aln);
        // }
        // tbb::this_task_arena::isolate( [&]{
        // tbb::parallel_for(tbb::blocked_range<int>(0, fallbackPairs.size()), [&](tbb::blocked_range<int> range) {
        // for (int n = range.begin(); n < range.end(); ++n) {
        //     int nIdx = fallbackPairs[n];
        //     std::vector<int8_t> aln;
        //     std::vector<std::vector<int32_t>> freqRef;
        //     std::vector<std::vector<int32_t>> freqQry;
        //     int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
        //     int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
        //     for (int r = 0; r < refLen; r++) {
        //         std::vector<int32_t> temp;
        //         for (int f = 0; f < 6; ++f) temp.push_back(0);
        //         freqRef.push_back(temp);
        //     }
        //     for (int q = 0; q < qryLen; q++) {
        //         std::vector<int32_t> temp;
        //         for (int f = 0; f < 6; ++f) temp.push_back(0);
        //         freqQry.push_back(temp);
        //     }           
        //     for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
        //         int storage = util->seqsStorage[sIdx];
        //         // int maxLen = max(refLen, qryLen);
        //         tbb::this_task_arena::isolate( [&]{
        //         tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
        //         for (int s = r.begin(); s < r.end(); ++s) {
        //         // for (int s = 0; s < refLen; ++s) {
        //             if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freqRef[s][0]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freqRef[s][1]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freqRef[s][2]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
        //                      util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freqRef[s][3]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') freqRef[s][4]+=1;
        //             else                                                                                             freqRef[s][5]+=1;
        //         }
        //         });
        //         });
        //     }
        //     for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
        //         int storage = util->seqsStorage[sIdx];
        //         tbb::this_task_arena::isolate( [&]{
        //         tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
        //         for (int s = r.begin(); s < r.end(); ++s) {
        //         // for (int s = 0; s < qryLen; ++s) {
        //             if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freqQry[s][0]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freqQry[s][1]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freqQry[s][2]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
        //                      util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freqQry[s][3]+=1;
        //             else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') freqQry[s][4]+=1;
        //             else                                                                                             freqQry[s][5]+=1;
        //         }
        //         });
        //         });
        //     }
        //     Talco_xdrop::Params talco_params(hostParam);
        //     while (aln.empty()) {
        //         int16_t errorType = 0;
        //         Talco_xdrop::Align_freq (
        //             talco_params,
        //             freqRef,
        //             freqQry,
        //             aln,
        //             errorType
        //         );
        //         // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
        //         if (errorType == 1) {
        //             std::cout << "Updated x-drop value on No. " << nIdx << '\n';
        //             talco_params.updateXDrop(2*talco_params.xdrop);
        //         }
        //         if (errorType == 2) {
        //             std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
        //             talco_params.updateFLen(talco_params.fLen << 1);
        //         }
        //     }
        //     alns[n] = aln;
        //     // if (aln.empty()) {
        //     //     alignGrpToGrp_traditional (
        //     //         freqRef,
        //     //         freqQry,
        //     //         refLen,
        //     //         qryLen,
        //     //         param,
        //     //         aln
        //     //     );
        //     // }
        // }
        // });
        // });


        std::vector<std::vector<int8_t>> alns;
        for (int i = 0; i < fallbackPairs.size(); ++i) alns.push_back(alnCpu[fallbackPairs[i]]);
        alnCpu.clear();
        int maxAlnLen = 0;
        for (auto aln: alns) maxAlnLen = (aln.size() > maxAlnLen) ? aln.size() : maxAlnLen;
        util->memCheck(maxAlnLen);

        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, fallbackPairs.size()), [&](tbb::blocked_range<int> range) {
        for (int n = range.begin(); n < range.end(); ++n) {
            int nIdx = fallbackPairs[n];
            auto aln = alns[n];
            // if (nodes[nIdx].first->identifier == "node_2") for (auto a: aln) std::cout << (a&0xFFFF);
            // if (nodes[nIdx].first->identifier == "node_2") std::cout << '\n';
            
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
            }
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
            // if (nodes[nIdx].first->identifier == "node_2") printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
            printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
        }
        });
        });
    }
    free(hostParam);
            

    
    
    
    return;
}
*/


void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option, Params& param)
{

    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        // is leaf
        if (n.first->is_leaf() && n.first->msaIdx.empty()) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.first->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->msaIdx.empty()) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.second->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        // if (nodes.size() < 5) {
        //     std::cout << n.first->identifier << '(' << util->seqsLen[n.first->identifier] << ')' <<  tree->allNodes[n.first->identifier]->msaIdx.size() << '\t';
        //     std::cout << n.second->identifier << '(' << util->seqsLen[n.second->identifier] << ')'<< tree->allNodes[n.second->identifier]->msaIdx.size() << '\n';
        // }
    }

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = option->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
            
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
            
    
    paramType* hostParam = (paramType*)malloc(29 * sizeof(paramType)); 
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.scoringMatrix[i][j];
    hostParam[25] = param.gapOpen;
    hostParam[26] = param.gapExtend;
    hostParam[27] = param.gapClose;
    hostParam[28] = param.xdrop;
    
    

    // std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    // allocate memory on host and device
    float** hostFreq = new float* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostNum = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];
    
    float** hostGapOp = new float* [gpuNum]; // gap open
    float** hostGapEx = new float* [gpuNum]; // gap extend
    float** hostGapCl = new float* [gpuNum]; // gap close

    float** deviceFreq = new float* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceNum = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    
    float** deviceGapOp = new float* [gpuNum];
    float** deviceGapEx = new float* [gpuNum];
    float** deviceGapCl = new float* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];

    // float** rawFreq =  new float* [gpuNum];
    // float** rawGapOp = new float* [gpuNum];
    // float** rawGapCt = new float* [gpuNum];
    // float** rawGapEn = new float* [gpuNum];
    // for (int gn = 0; gn < gpuNum; ++gn) {
    //     rawFreq[gn] = (float*)malloc(12 * seqLen * sizeof(float));    
    //     rawGapOp[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    //     rawGapCt[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    //     rawGapEn[gn] = (float*)malloc(2 * seqLen * sizeof(float));
    // }

    std::atomic<int> nowRound;
    nowRound.store(0);
    // nowRound.store(roundGPU-1);

    // int ThreadsPerGPU = maxThreads / gpuNum;
    bool* cpuFallback = new bool[nodes.size()];
    std::vector<std::vector<int8_t>> alnCpu;
    for (int i = 0; i < nodes.size(); ++i) {
        cpuFallback[i] = false;
        if (util->nowProcess == 1) {
            std::vector<int8_t> aln;
            alnCpu.push_back(aln);
        }
        
    }
    
            

    // tbb::task_arena limited_arena(2);
    // limited_arena.execute([&]{

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] =  (float*)malloc(12 * seqLen * numBlocks * sizeof(float));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostNum[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(2                     * sizeof(int32_t));
            hostGapOp[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapEx[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            hostGapCl[gn] = (float*)malloc(2 * seqLen * numBlocks * sizeof(float));
            
            cudaSetDevice(option->gpuIdx[gn]);
            // cudaSetDevice(1);
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceNum[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 2 * sizeof(int32_t));
            cudaMalloc((void**)&deviceGapOp[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapEx[gn],  2 * seqLen * numBlocks * sizeof(float));
            cudaMalloc((void**)&deviceGapCl[gn],  2 * seqLen * numBlocks * sizeof(float));

            cudaMalloc((void**)&deviceParam[gn],  29 * sizeof(paramType));
            
            std::string error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda malloc %s!\n", error.c_str());
            cudaMemcpy(deviceParam[gn], hostParam, 29 * sizeof(paramType), cudaMemcpyHostToDevice);
            error = cudaGetErrorString(cudaGetLastError()); // printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            if (error != "no error") printf("ERROR: Cuda copy param %s!\n", error.c_str());
            std::vector<std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>> gappyColumns;
            
            while (nowRound < roundGPU) {
            // while (nowRound >= 0) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                gappyColumns.clear();
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostNum[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapOp[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapEx[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostGapCl[gn][n] = 0;

                // Calculate Frequency
                // tbb::task_arena ta {maxThreads - gpuNum};
                // tbb::task_group tg;
                // tg.run ([&] () {
                // tbb::this_task_arena::isolate ([&] () {
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> aln_range) { 
                // for (int n = aln_range.begin(); n < aln_range.end(); ++n) {
                
                
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    float refWeight = 0.0; float qryWeight = 0.0;

                    // for (int s = 0; s < 12*seqLen; ++s) rawFreq[gn][s] = 0.0;
                    // for (int s = 0; s < 2*seqLen ; ++s) {rawGapOp[gn][s] = 0.0; rawGapCt[gn][s] = 0.0; rawGapEn[gn][s] = 0.0;}
                    // get sum of weights
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        float w = tree->allNodes[name]->weight / refWeight * refNum;
                        // float w = 1.0;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < refLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*s+5]+=1.0*w;
                            if (s > 0) {
                                if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+s] += 1.0;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+s] += 1.0;
                            }
                            
                        }
                        });
                        });
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        std::string name = util->seqsName[sIdx];
                        // float w = 1.0;
                        float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < qryLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1.0*w;
                            else if (util->alnStorage[storage][sIdx][s] == '-')                                              hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1.0*w;
                            if (s > 0) {
                                if (util->alnStorage[storage][sIdx][s-1] != '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapOp[gn][2*seqLen*n+seqLen+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] == '-') hostGapEx[gn][2*seqLen*n+seqLen+s] += 1.0;
                                // if (                                               util->alnStorage[storage][sIdx][s] == '-') hostGapCt[gn][2*seqLen*n+seqLen+s] += 1.0;
                                if (util->alnStorage[storage][sIdx][s-1] == '-' && util->alnStorage[storage][sIdx][s] != '-') hostGapCl[gn][2*seqLen*n+seqLen+s] += 1.0;
                            }
                        }
                        });
                        });
                    }
                    for (int s = 0; s < refLen; ++s) { hostGapOp[gn][2*seqLen*n+s] /= refNum; hostGapEx[gn][2*seqLen*n+s] /= refNum; hostGapCl[gn][2*seqLen*n+s] /= refNum;} 
                    for (int s = 0; s < qryLen; ++s) { hostGapOp[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapEx[gn][2*seqLen*n+seqLen+s] /= qryNum; hostGapCl[gn][2*seqLen*n+seqLen+s] /= qryNum;}
                    
                    // ClustalW's method
                    
                    for (int s = 0; s < refLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*s+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+s] = param.gapOpen * 0.3 * ((refNum-hostFreq[gn][12*seqLen*n+6*s+5]) * 1.0 / refNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*d+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    for (int s = 0; s < qryLen; ++s) {
                        if (hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5] > 0) {
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = param.gapOpen * 0.3 * ((qryNum-hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]) * 1.0 / qryNum);
                        }
                        else {
                            int backSites = 8;
                            int distance_from_gap = 0;
                            int increPenalty = false;
                            for (int d = s-1; d > max(s-backSites, 0); --d) {
                                ++distance_from_gap;
                                if (hostFreq[gn][12*seqLen*n+6*(seqLen+d)+5] > 0) {
                                    increPenalty = true;
                                    break;
                                }
                            }
                            hostGapOp[gn][2*seqLen*n+seqLen+s] = (increPenalty) ? param.gapOpen * static_cast<paramType>((2 + ((backSites - distance_from_gap)*2.0)/backSites)) : param.gapOpen;
                        }
                    }
                    
                    std::queue<std::pair<int,int>> gappyRef, gappyQry;
                    if (option->gappyHorizon > 0) {
                        float gappyVertical = option->gappyVertical;
                        int gappyHorizon = option->gappyHorizon, gappyLength;
                        int rawIdx = 0, newIdx = 0;
                        int gapRef = 0, gapQry = 0;
                        while (true) {
                            if (rawIdx >= refLen) {
                                for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*i+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*gappyLength+5]/refNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyRef.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHR: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapRef += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, refLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*newIdx+t] = hostFreq[gn][12*seqLen*n+6*rawIdx+t];
                                hostGapOp[gn][2*seqLen*n+newIdx] = hostGapOp[gn][2*seqLen*n+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+newIdx] = hostGapEx[gn][2*seqLen*n+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+newIdx] = hostGapCl[gn][2*seqLen*n+rawIdx];
                                ++newIdx;
                            } 
                        }
                        hostLen[gn][2*n] = newIdx; 
                        // std::cout << n << '#' << rawIdx << '/' << newIdx << '\t';
                        rawIdx = 0; newIdx = 0;
                        while (rawIdx < seqLen) {
                            if (rawIdx == qryLen) {
                                for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < 6; ++j) hostFreq[gn][12*seqLen*n+6*(seqLen+i)+j] = 0;
                                break;
                            }
                            for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) if (hostFreq[gn][12*seqLen*n+6*(seqLen+gappyLength)+5]/qryNum <= gappyVertical) break;
                            if (gappyLength - rawIdx >= gappyHorizon) {
                                gappyQry.push(std::make_pair(rawIdx, gappyLength-rawIdx)); 
                                // if (alnPairs < 2) std::cout << "PUSHQ: " << rawIdx << '-' << gappyLength-rawIdx << '\n';
                                gapQry += gappyLength-rawIdx;
                                rawIdx += gappyLength-rawIdx;
                            }
                            // if (n == 0 && alnPairs == 3) std::cout << newIdx << ':' << newIdx + gapQry << ':' << rawIdx << ':' << qryLen << ':' << seqLen << '\n';
                            for (rawIdx = rawIdx; rawIdx < min(gappyLength+1, qryLen); ++rawIdx) {
                                for (int t = 0; t < 6; ++t) hostFreq[gn][12*seqLen*n+6*(seqLen+newIdx)+t] = hostFreq[gn][12*seqLen*n+6*(seqLen+rawIdx)+t];
                                hostGapOp[gn][2*seqLen*n+seqLen+newIdx] = hostGapOp[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapEx[gn][2*seqLen*n+seqLen+newIdx] = hostGapEx[gn][2*seqLen*n+seqLen+rawIdx]; 
                                hostGapCl[gn][2*seqLen*n+seqLen+newIdx] = hostGapCl[gn][2*seqLen*n+seqLen+rawIdx];
                                ++newIdx;
                            } 
                        }
                        if (!gappyRef.empty()) std::cout << "No. " << n << ": Removing " << gapRef << " columns from refernce.\n";
                        if (!gappyQry.empty()) std::cout << "No. " << n << ": Removing " << gapQry << " columns from query.\n";
                    
                        hostLen[gn][2*n+1] = newIdx;
                        assert(hostLen[gn][2*n] + gapRef == refLen);
                        assert(hostLen[gn][2*n+1] + gapQry == qryLen);
                    }
                    else {
                        hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    }
                    // std::cout << rawIdx << '/' << newIdx << '\n';
                    hostNum[gn][2*n] = refNum; hostNum[gn][2*n+1] = qryNum;
                    gappyColumns.push_back(std::make_pair(gappyRef, gappyQry));
                }
                    
                // for (int i = 0; i < 50; ++i) std::cout << hostGapOp[gn][i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapEx[gn][i] << ',';
                // std::cout << '\n';
                // for (int i = 0; i < 50; ++i) std::cout << hostGapCl[gn][i] << ',';
                // std::cout << '\n';
                
                hostSeqInfo[gn][0] = alnPairs;
                hostSeqInfo[gn][1] = seqLen;
                
        
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceNum[gn], hostNum[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 2 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapOp[gn], hostGapOp[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapEx[gn], hostGapEx[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceGapCl[gn], hostGapCl[gn], 2 * seqLen * numBlocks * sizeof(float), cudaMemcpyHostToDevice);

                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceNum[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceGapOp[gn],
                    deviceGapEx[gn],
                    deviceGapCl[gn],
                    deviceParam[gn]
                );
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();

                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") {
                    printf("ERROR: After kernel %s!\n", aerr.c_str());
                    exit(1);
                }

                int maxAlnLen = 0;
                for (int n = 0; n <  alnPairs; ++n) {
                    if (hostAlnLen[gn][n] > maxAlnLen) maxAlnLen = hostAlnLen[gn][n];
                    // std::cout << n << ':' <<  hostAlnLen[gn][n] << '\n'; 
                }
                util->memCheck(maxAlnLen);
                if (rn % 10 == 0 && rn > 0) std::cout << rn*numBlocks << " pairs have been processed.\n";
                
                // tbb::this_task_arena::isolate( [&]{
                // tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
                // for (int n = range.begin(); n < range.end(); ++n) {
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
                    int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
                    if (hostAlnLen[gn][n] <= 0) {
                        cpuFallback[nIdx] = true;
                        if (util->nowProcess == 1) {
                            std::vector<int8_t> aln_old, aln;
                            std::vector<std::vector<float>> freqRef;
                            std::vector<std::vector<float>> freqQry;
                            std::vector<std::vector<float>> gapOp, gapEx, gapCl;
                            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                            for (int r = 0; r < hostLen[gn][2*n]; r++) {
                                std::vector<float> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*r+f]);
                                freqRef.push_back(temp);
                            }
                            for (int q = 0; q < hostLen[gn][2*n+1]; q++) {
                                std::vector<float> temp;
                                for (int f = 0; f < 6; ++f) temp.push_back(hostFreq[gn][12*seqLen*n+6*(seqLen+q)+f]);
                                freqQry.push_back(temp);
                            }
                            for (int i = 0; i < 2; ++i) {
                                std::vector<float> temp;
                                gapOp.push_back(temp); gapEx.push_back(temp); gapCl.push_back(temp);
                            }
                            for (int r = 0; r < hostLen[gn][2*n]; ++r) {
                                gapOp[0].push_back(hostGapOp[gn][2*n*seqLen+r]);
                                gapEx[0].push_back(hostGapEx[gn][2*n*seqLen+r]);
                                gapCl[0].push_back(hostGapCl[gn][2*n*seqLen+r]);
                            }
                            for (int q = 0; q < hostLen[gn][2*n+1]; ++q) {
                                gapOp[1].push_back(hostGapOp[gn][2*n*seqLen+seqLen+q]);
                                gapEx[1].push_back(hostGapEx[gn][2*n*seqLen+seqLen+q]);
                                gapCl[1].push_back(hostGapCl[gn][2*n*seqLen+seqLen+q]);
                            }
                            
                            Talco_xdrop::Params talco_params(hostParam);
                            while (aln_old.empty()) {
                                int16_t errorType = 0;
                                Talco_xdrop::Align_freq (
                                    talco_params,
                                    freqRef,
                                    freqQry,
                                    gapOp,
                                    gapEx,
                                    gapCl,
                                    aln_old,
                                    errorType
                                );
                                // assert((aln.empty() && errorType != 0) || (!aln.empty() && errorType == 0));
                                if (errorType == 1) {
                                    std::cout << "Updated x-drop value on No. " << nIdx << '\n';
                                    talco_params.updateXDrop(talco_params.xdrop << 1);
                                }
                                if (errorType == 2) {
                                    std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                                    talco_params.updateFLen(talco_params.fLen << 1);
                                }
                            }
                            int rIdx = 0, qIdx = 0, j = 0; 
                            // for (int j = 0; j < aln_old.size(); ++j) {
                            while (j < aln_old.size() || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty())) {
                            // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                                bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                                bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                                int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                                int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                                if (gapR || gapQ) {
                                    // std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                    if (gapRLen >= gapQLen) {
                                        for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                        for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                    }
                                    else {
                                        for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                        for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                    }
                                    if (gapR) gappyColumns[n].first.pop();
                                    if (gapQ) gappyColumns[n].second.pop();
                                }
                                else {
                                    switch ((aln_old[j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                    ++j;
                                }
                                if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                    for (j = j; j < aln_old.size(); ++j) {
                                        switch ((aln_old[j] & 0xFFFF)) {
                                            case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                            case 1: ++qIdx;         aln.push_back(1); break;
                                            case 2: ++rIdx;         aln.push_back(2); break;
                                        }
                                    }
                                }
                            }
                            // std::cout << n << '#' << rIdx << '-' << refLen << '^' << qIdx << '-' << qryLen << '\n';
                            assert(rIdx == refLen); assert(qIdx == qryLen);
                            alnCpu[nIdx] = aln;
                        }
                    }
                    else {
                        std::vector<int8_t> aln;
                        int rIdx = 0, qIdx = 0; 
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        int j = 0;
                        int oneInAln = 0;
                        while (j < hostAlnLen[gn][n] || (!gappyColumns[n].first.empty() || !gappyColumns[n].second.empty()) ) {
                        // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            bool gapR = (gappyColumns[n].first.empty())  ? false : (rIdx == gappyColumns[n].first.front().first);
                            bool gapQ = (gappyColumns[n].second.empty()) ? false : (qIdx == gappyColumns[n].second.front().first);
                            int gapRLen = (!gapR) ? 0 : gappyColumns[n].first.front().second;
                            int gapQLen = (!gapQ) ? 0 : gappyColumns[n].second.front().second;
                            if (gapR || gapQ) {
                                // if (n == 0 ) std::cout << "No." << n << '\t' << rIdx << ':' << gapR << '/' << gapRLen << '\t' << qIdx << ':' << gapQ << '/' << gapQLen << '$' << (gappyColumns[n].second.empty() ? 0: gappyColumns[n].second.front().first) <<'\n';
                                if (gapRLen >= gapQLen) {
                                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                                }
                                else {
                                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                                }
                                if (gapR) gappyColumns[n].first.pop();
                                if (gapQ) gappyColumns[n].second.pop();
                            }
                            else {
                                switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                    case 0: ++rIdx; ++qIdx; aln.push_back(0); ++oneInAln; break;
                                    case 1: ++qIdx;         aln.push_back(1); ++oneInAln; break;
                                    case 2: ++rIdx;         aln.push_back(2); break;
                                }
                                ++j;
                            }
                            if (gappyColumns[n].first.empty() && gappyColumns[n].second.empty()) {
                                for (j = j; j < hostAlnLen[gn][n]; ++j) {
                                    switch ((hostAln[gn][n*2*seqLen+j] & 0xFFFF)) {
                                        case 0: ++rIdx; ++qIdx; aln.push_back(0); break;
                                        case 1: ++qIdx;         aln.push_back(1); break;
                                        case 2: ++rIdx;         aln.push_back(2); break;
                                    }
                                }
                            }
                        }
                        util->memCheck(aln.size());
                        if(rIdx != refLen || qIdx != qryLen) std::cout << n << '#' << rIdx << '-' << refLen << '^' << qIdx <<':' << oneInAln << '-' << qryLen << '\n';
                        // if (n == 0) std::cout << (aln[0] & 0xFFFF) << '\n';
                        // if (n == 0 && alnPairs == 3) {
                        //     for (int i = 0; i < aln.size(); ++i) {
                        //         if (i % 10 == 0) std::cout << '\n' << i/10 << '\t';
                        //         std::cout << (aln[i] & 0xFFFF);
                        //         
                        //     }
                        // }
                        assert(rIdx == refLen); assert(qIdx == qryLen);
                        for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                            int orgIdx = 0;
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {   
                            for (int j = 0; j < aln.size(); ++j) {
                                // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 2) {
                                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                            }
                            // util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][n];
                            util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                            util->changeStorage(sIdx);
                        }
                        for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            int storeFrom = util->seqsStorage[sIdx];
                            int storeTo = 1 - util->seqsStorage[sIdx];
                            int orgIdx = 0;
                            // for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            for (int j = 0; j < aln.size(); ++j) {
                                // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
                                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                                    util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                                    orgIdx++;
                                    
                                }
                                else {
                                    util->alnStorage[storeTo][sIdx][j] = '-';
                                }
                                
                            }
                            // util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
                            util->seqsLen[nodes[nIdx].second->identifier] = aln.size();  
                            util->changeStorage(sIdx);
                        }
                        for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                        }
                        tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                    }
                    
                }
                // });
                // });
            }  

            
        }
    });
    // });
    
    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(option->gpuIdx[gn]);
        // cudaSetDevice(gn);
        // cudaSetDevice(1);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceNum[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceGapOp[gn]);
        cudaFree(deviceGapEx[gn]);
        cudaFree(deviceGapCl[gn]);
        cudaDeviceSynchronize();  
        std::string freeErr = cudaGetErrorString(cudaGetLastError());
        if (freeErr != "no error") {
            printf("ERROR: Free memory Last%s!\n", freeErr.c_str());
            exit(1);
        }
        free(hostFreq[gn]);
        free(hostAln[gn]);
        free(hostLen[gn]);
        free(hostNum[gn]);
        free(hostAlnLen[gn]);
        free(hostSeqInfo[gn]);
        free(hostGapOp[gn]);
        free(hostGapEx[gn]);
        free(hostGapCl[gn]);
        // free(rawGapOp[gn]);
        // free(rawGapCt[gn]);
        // free(rawGapEn[gn]); 
        // free(rawFreq[gn]);    
    }
    
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] deviceLen;
    delete [] deviceNum;
    delete [] deviceGapOp;
    delete [] deviceGapEx;
    delete [] deviceGapCl;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostNum;
    delete [] hostAln;
    delete [] hostSeqInfo;
    delete [] hostGapOp;
    delete [] hostGapEx;
    delete [] hostGapCl;
    // delete [] rawFreq;
    // delete [] rawGapOp;  
    // delete [] rawGapCt;
    // delete [] rawGapEn;    
    // CPU Fallback
    std::vector<int> fallbackPairs;
    for (int i = 0; i < nodes.size(); ++i) if (cpuFallback[i]) fallbackPairs.push_back(i);
    delete [] cpuFallback;
    
    if (fallbackPairs.empty()) {
        free(hostParam);
        return;
    }
    if (util->nowProcess == 0) {
        std::cout << "Bad alignments. Num of pairs: " << fallbackPairs.size() << '\n';
        for (int i = 0; i < fallbackPairs.size(); ++i) {
            int nIdx = fallbackPairs[i];
            int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
            int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
            int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
            if (util->badSequences.find(grpID) == util->badSequences.end()) {
                std::vector<std::string> temp;
                util->badSequences[grpID] = temp;
            }
            if (refNum < qryNum) {
                // proceed with qry and store ref as bad sequences
                int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];

                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                util->seqsLen[nodes[nIdx].second->identifier] = refLen;
                util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
                auto temp = nodes[nIdx].second->msaIdx;
                tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << refNum <<" seqeuences).\n";
            }
            else {
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                std::cout << "Deferring the profile on " << nodes[nIdx].second->identifier << " (" << qryNum <<" seqeuences).\n";
            }
        }
    }
    else {
        std::cout << "CPU Fallback. Num of pairs: " << fallbackPairs.size() << '\n';
        std::vector<std::vector<int8_t>> alns;
        for (int i = 0; i < fallbackPairs.size(); ++i) alns.push_back(alnCpu[fallbackPairs[i]]);
        alnCpu.clear();
        int maxAlnLen = 0;
        for (auto aln: alns) maxAlnLen = (aln.size() > maxAlnLen) ? aln.size() : maxAlnLen;
        util->memCheck(maxAlnLen);

        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, fallbackPairs.size()), [&](tbb::blocked_range<int> range) {
        for (int n = range.begin(); n < range.end(); ++n) {
            int nIdx = fallbackPairs[n];
            auto aln = alns[n];
            // if (nodes[nIdx].first->identifier == "node_2") for (auto a: aln) std::cout << (a&0xFFFF);
            // if (nodes[nIdx].first->identifier == "node_2") std::cout << '\n';
            
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int j = 0; j < aln.size(); ++j) {
                    if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                        util->alnStorage[storeTo][sIdx][j] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][j] = '-';
                    }
                }
                util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                util->changeStorage(sIdx);
            }
            for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
            }
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
            // if (nodes[nIdx].first->identifier == "node_2") printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
            printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", nIdx, nodes[nIdx].first->identifier.c_str(), aln.size());
        }
        });
        });
    }
    free(hostParam);
            
    return;
}


void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID) {
    
    std::map<std::string, int> NodeAlnOrder;

    while(!postOrder.empty()) {
        Node* node = postOrder.top();
        // Leaf or not belongs to the subtree
        if (!(node->grpID==-1 || node->grpID==grpID) || node->is_leaf()) {
            postOrder.pop();
            continue;
        }
        // Collect children that belong to the subtree
        std::vector<Node*> children;
        for (int i = 0; i < node->children.size(); ++i) {
            if (node->children[i]->grpID == grpID) children.push_back(node->children[i]);
        }
        // Useless node, remove from the subtree
        if (children.empty()) {
            node->grpID = -2;
            postOrder.pop();
            continue;
        }
        // Only one child, merge the child and the parent 
        if (children.size() == 1 && node->parent != nullptr) {
            if (node->parent->grpID == grpID) {
                for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                    if (node->parent->children[chIdx]->identifier == node->identifier) {
                        node->parent->children[chIdx] = children[0];
                    }
                }
                postOrder.pop();
                continue;
            }
            
        }
        // Replace the first child with node
        NodeAlnOrder[node->identifier] = (NodeAlnOrder.find(children[0]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[0]->identifier] : -1;
        children[0] = node;
        while (children.size() > 1) {
            std::vector<Node*> nodeLeft;
            for (int i = 0; i < children.size()-1; i+=2) {
                int firstIdx  = (NodeAlnOrder.find(children[i]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i]->identifier]+1 : 0;
                int secondIdx = (NodeAlnOrder.find(children[i+1]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i+1]->identifier]+1 : 0;
                int maxIdx = max(firstIdx, secondIdx);
                NodeAlnOrder[children[i]->identifier] = maxIdx;
                NodeAlnOrder[children[i+1]->identifier] = maxIdx;
                alnOrder.push_back(std::make_pair(std::make_pair(children[i], children[i+1]), maxIdx));
                nodeLeft.push_back(children[i]);
            }
            if (children.size()%2 == 1) nodeLeft.push_back(children.back());
            children = nodeLeft;
        }
        postOrder.pop();
    }
}



__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result) {
    int32_t seqNum     = seqInfo[0]; 
    int32_t seqLen     = seqInfo[1]; 
    int32_t numBlocks  = seqInfo[2]; 
    int32_t numThreads = seqInfo[3]; 
    int32_t match      = seqInfo[4]; 
    int32_t mismatch   = seqInfo[5]; 
    int32_t gap        = seqInfo[6]; 

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;
    const int threadNum = 512;
    // int tidx = bx*bs+tx;
    if (bx < numBlocks) {
        __shared__ int64_t columnScore [threadNum];
        if (tx < numThreads) {
            columnScore[tx] = 0;
        }
        __syncthreads();
        for (int l = bx; l < seqLen; l=l+gs) {
            for (int i = 0; i < seqNum-1; ++i) {
                if (tx < numThreads) {
                    for (int j = i + 1 + tx; j < seqNum; j=j+bs) {
                        if (seqs[i*seqLen+l] == 'N' || seqs[j*seqLen+l] == 'N' || 
                            seqs[i*seqLen+l] == 'n' || seqs[j*seqLen+l] == 'n' ||
                           (seqs[i*seqLen+l] == '-' && seqs[j*seqLen+l] == '-' )) {
                            columnScore[tx] += 0;
                        }
                        else if (seqs[i*seqLen+l] == '-' || seqs[j*seqLen+l] == '-') {
                            columnScore[tx] += gap;
                        }
                        else if (seqs[i*seqLen+l] == seqs[j*seqLen+l]) {
                            columnScore[tx] += match;
                        }
                        else {
                            columnScore[tx] += mismatch;
                        }
                    }
                }
                __syncthreads();
            }
        }
        for (uint32_t r = threadNum/2; r > 0; r >>= 1) {
            if (tx < r) {
                columnScore[tx]   = columnScore[tx] + columnScore[tx+r];
            }
            __syncthreads();
        }
        if (tx == 0) result[bx] = columnScore[0];
    }

}

double getSPScore_gpu(msa::utility* util, Params& param) {
    // auto scoreStart = std::chrono::high_resolution_clock::now();
    int numBlocks = 1024;
    int blockSize = 512;
    // size_t seqNum = util->memNum;
    // size_t seqLen = util->seqsLen["node_1"];
    size_t seqNum = util->memNum;
    size_t seqLen = 0;

    int storage = util->seqsStorage[0];
    while (util->alnStorage[storage][0][seqLen] != 0) {
        ++seqLen;
    }
    printf("(Num, Len) - (%lu, %lu)\n", seqNum, seqLen);
    
    char*    hostSeqs = (char*)malloc(seqLen * seqNum * sizeof(char));
    int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
    int64_t* hostResult = (int64_t*)malloc(numBlocks * sizeof(int64_t));
    // int seqCount = 0;
    for (int i = 0; i < seqNum; ++i) {
        int storage = util->seqsStorage[i];
        for (int j = 0; j < seqLen; ++j) {
            hostSeqs[i*seqLen+j] = util->alnStorage[storage][i][j];
        }
    }
    // for (int j = 0; j < seqLen*seqNum; ++j) { 
    //     if (j%seqLen < alignment[seqCount].size()) {
    //         hostSeqs[j] = util->seqs[seqCount][j%seqLen];
    //     }
    //     else hostSeqs[j] = 0;
    //     if (j%seqLen == seqLen-1) ++seqCount;
    // }
    

    hostSeqInfo[0] = seqNum;
    hostSeqInfo[1] = seqLen;
    hostSeqInfo[2] = numBlocks;
    hostSeqInfo[3] = blockSize;
    hostSeqInfo[4] = param.scoringMatrix[0][0];
    hostSeqInfo[5] = param.scoringMatrix[0][1];
    hostSeqInfo[6] = param.gapExtend;
    for (int i = 0; i < numBlocks; ++i) hostResult[i] = 0;

    char*    deviceSeqs;
    int32_t* deviceSeqInfo;
    int64_t* deviceResult;

    cudaMalloc((void**)&deviceSeqs, seqLen * seqNum * sizeof(char));
    cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
    cudaMalloc((void**)&deviceResult, numBlocks * sizeof(int64_t));

    cudaMemcpy(deviceSeqs, hostSeqs, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceResult, hostResult, numBlocks * sizeof(int64_t), cudaMemcpyHostToDevice);

    calSPScore<<<numBlocks, blockSize>>>(
        deviceSeqs, 
        deviceSeqInfo,
        deviceResult
    );
    cudaDeviceSynchronize();
    cudaMemcpy(hostResult, deviceResult, numBlocks * sizeof(int64_t), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    int64_t gpuScore = 0;
    for (int i = 0; i < numBlocks; ++i) {
        // std::cout << "Block: "<<i<<", score: " << hostResult[i] << '\n';
        gpuScore += hostResult[i];
    }
    double score = static_cast<double>(gpuScore)/static_cast<double>(seqNum*(seqNum-1)/2);
    // std::cout << "GPU score: " << score << '\n';
    // auto scoreEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds scoreTime = scoreEnd - scoreStart;
    
    // printf("GPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
    
    cudaFree(deviceSeqs);
    cudaFree(deviceSeqInfo);
    cudaFree(deviceResult);
    free(hostSeqs);
    free(hostSeqInfo);
    free(hostResult);
    return score;
} 



// Regressive Method (not maintained)

/*
void getLongestDescendent(Tree* tree, msa::utility* util) {
    std::stack<Node*> postOrder;
    getPostOrderList(tree->root, postOrder);
    while (!postOrder.empty()) {
        Node* n = postOrder.top();
        postOrder.pop();
        if (n->children.empty()) n->setLongestDescendant(n);
        else if (n->children.size() == 1) n->setLongestDescendant(n->children[0]->longestDescendant);
        else {
            Node* longest = n->children[0]->longestDescendant;
            for (int i = 1; i < n->children.size(); ++i) {
                if (util->seqs[n->children[i]->longestDescendant->identifier].size() > util->seqs[longest->identifier].size()) {
                    longest = n->children[i]->longestDescendant;
                }
            }
            n->setLongestDescendant(longest);
        }
    }
    return;

}

void getSubMSANodes(std::vector<Node*>& subMSANodes, Node* startNode, int N) {
    if (startNode->is_leaf()) return;
    std::queue<Node*> visited;
    std::queue<Node*> added;
    visited.push(startNode);
    int minNumNodes = N;
    while (visited.size() < minNumNodes && !visited.empty()) {
        
        Node* currentNode = visited.front();
        visited.pop();
        // std::cout << "start: " << currentNode->identifier << '\n';
        if (!currentNode->children.empty()) {
            for (auto ch: currentNode->children) {
                // std::cout << "ch: " << ch->identifier << '\n';
                visited.push(ch);
            }
        }
        else {
            added.push(currentNode);
            minNumNodes -= 1;
        }
    }
    while (!visited.empty()) {
        added.push(visited.front());
        visited.pop();
    }
    while (!added.empty()) {
        subMSANodes.push_back(added.front());
        added.pop();
    }
    // for (auto a: subMSANodes) std::cout << a->longestDescendant->identifier << ',';
    // std::cout << '\n';
    return;
    
}

void getSubMSAs(std::map<int, std::vector<std::vector<Node*>>>& subMSAs, Tree* T, int N) {
    // subMSAs (Level, clusters)
    subMSAs.clear();
    int level = 0;
    Node* subRoot = T->root;
    std::vector<Node*> parentMSA;
    std::vector<std::pair<int, std::vector<Node*>>> subMSAList;
    getSubMSANodes(parentMSA, subRoot, N);
    subMSAList.push_back(std::make_pair(level, parentMSA));
    int start = 0, end = subMSAList.size();
    level++;
    while (true) {
        for (int idx = start; idx < end; ++idx) {
            for (auto node: subMSAList[idx].second) {
                // std::cout << node->identifier << '\n';
                std::vector<Node*> childMSA;
                getSubMSANodes(childMSA, node, N);
                if (!childMSA.empty()) {
                    subMSAList.push_back(std::make_pair(level, childMSA));
                }
            }
        }
        start = end;
        end = subMSAList.size();
        if (start == end) break;
        ++level;
    }
    for (int i = 0; i < level; ++i) {
        std::vector<std::vector<Node*>> temp;
        subMSAs[i] = temp;
    }
    for (auto msa: subMSAList) {
        subMSAs[msa.first].push_back(msa.second);
    }
    // for (auto it = subMSAs.begin(); it != subMSAs.end(); ++it) {
    //     std::cout << "Level: " << it->first << '\n';
    //     for (auto s: it->second) {
    //         for (auto ss: s) std::cout << ss->longestDescendant->identifier << ',';
    //         std::cout << '\n';
    //     }
    // }
    return;
}

void getAlnPairs(std::vector<std::vector<std::pair<Node*, Node*>>>& alnPairs, std::vector<std::vector<Node*>>& clusters) {
    std::vector<std::vector<Node*>> currentLevel, nextLevel;
    currentLevel = clusters;
    while (!currentLevel.empty()) {
        std::vector<std::pair<Node*, Node*>> temp;
        int cidx = 0;
        for (auto cluster: currentLevel) {
            // std::cout << cidx << ':' << cluster.size() <<'\n';
            if (cluster.size() == 1) continue;
            std::vector<Node*> nextLevelTemp;
            for (int i = 0; i < cluster.size()-1; i+=2) {
                nextLevelTemp.push_back(cluster[i]);
                temp.push_back(std::make_pair(cluster[i], cluster[i+1]));
            }
            if (cluster.size()%2 == 1) nextLevelTemp.push_back(cluster.back());
            nextLevel.push_back(nextLevelTemp);
        }
        alnPairs.push_back(temp);
        // std::cout << nextLevel.size() << '\n';
        currentLevel = nextLevel;
        nextLevel.clear();
    }
    return;
}

void storeMSA(Tree* T, std::vector<Node*>& nodes, msa::utility* util, int level) {
    for (auto n: nodes) {
        assert(n->msa.size() == n->msaIdx.size());
        for (int i = 0; i < n->msa.size(); ++i) {
            int sIdx = n->msaIdx[i];
            std::string sIdentifier = n->msa[i];
            int storage = util->seqsStorage[sIdx];
            int j = 0;
            std::string msa = "";
            while (util->seqBuf[storage][sIdx][j] != 0) {
                msa += util->seqBuf[storage][sIdx][j];
                ++j;
            }
            T->allNodes[sIdentifier]->msaSeq[level] = msa;
        }
    }
    return;
}

void resetSeqMem(std::vector<Node*>& nodes, msa::utility* util) {
    for (auto n: nodes) {
        assert(n->msa.size() == n->msaIdx.size());
        for (int i = 0; i < n->msa.size(); ++i) {
            int sIdx = n->msaIdx[i];
            std::string sIdentifier = n->msa[i];
            int storage = util->seqsStorage[sIdx];
            std::string rawSeq = util->seqs[sIdentifier];
            for (int j = 0; j < util->memLen; ++j) {
                if (j < rawSeq.size()) util->seqBuf[storage][sIdx][j] = rawSeq[j];
                else util->seqBuf[storage][sIdx][j] = 0;
            }
            storage = 1 - storage;
            for (int j = 0; j < util->memLen; ++j) {
                util->seqBuf[storage][sIdx][j] = 0;
            }
        }
    }
    return;
}

std::string removeGap(std::string s) {
    std::string s_noGap = "";
    for (int i = 0; i < s.size(); ++i) if (s[i] != '-') s_noGap += s[i];
    return s_noGap;
}

void merger_ref(Tree* tree, std::map<int, Node*>& refNodes, msa::utility* util, std::string& refString, std::string& qryString, int qLevel) {
    int rIdx = 0, qIdx = 0, alnIdx = 0;
    int refLen = refString.size(), qryLen = qryString.size();
    if (removeGap(refString) != removeGap(qryString)) {
        std::cout << "Unmatched Seq.\n"; exit(1);
    }
    
    while (rIdx < refLen && qIdx < qryLen) {
        if (alnIdx > refLen && alnIdx > qryLen) util->memCheck(alnIdx);
        if (refString[rIdx] == qryString[qIdx] && refString[rIdx] != '-') {
            for (auto n: refNodes) {
                int sIdx = n.first;
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - storeFrom;
                util->seqBuf[storeTo][sIdx][alnIdx] = util->seqBuf[storeFrom][sIdx][rIdx];
            }
            ++alnIdx;++qIdx;++rIdx;
        }
        else if (refString[rIdx] == '-' && qryString[qIdx] != '-') {
            int consecGap = 0;
            int k = rIdx;
            while (refString[k] == '-' && k < refLen) {
                ++consecGap;
                ++k;
            }
            for (size_t g = 0; g < consecGap; ++g) {
                for (auto n: refNodes) {
                    int sIdx = n.first;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - storeFrom;
                    util->seqBuf[storeTo][sIdx][alnIdx] = util->seqBuf[storeFrom][sIdx][rIdx];
                }
                ++alnIdx;++rIdx;
            }
        }
        else if (refString[rIdx] != '-' && qryString[qIdx] == '-') {
            int consecGap = 0;
            int k = qIdx;
            while (qryString[k] == '-' && k < qryLen) {
                ++consecGap;
                ++k;
            }
            for (size_t g = 0; g < consecGap; ++g) {
                for (auto n: refNodes) {
                    int sIdx = n.first;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - storeFrom;
                    util->seqBuf[storeTo][sIdx][alnIdx] = '-';
                }
                ++alnIdx;++qIdx;
            }
        }
        else {
            int consecGap = 0;
            int kr = rIdx, kq = qIdx;
            while (refString[kr] == '-' && kr < refLen) {
                ++consecGap;
                ++kr;
            }
            for (size_t g = 0; g < consecGap; ++g) {
                for (auto n: refNodes) {
                    int sIdx = n.first;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - storeFrom;
                    util->seqBuf[storeTo][sIdx][alnIdx] = util->seqBuf[storeFrom][sIdx][rIdx];
                }
                ++alnIdx;++rIdx;
            }
            consecGap = 0;
            while (qryString[kq] == '-' && kq < qryLen) {
                ++consecGap;
                ++kq;
            }
            for (size_t g = 0; g < consecGap; ++g) {
                for (auto n: refNodes) {
                    int sIdx = n.first;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - storeFrom;
                    util->seqBuf[storeTo][sIdx][alnIdx] = '-';
                }
                ++alnIdx;++qIdx;
            }
        }
    }
    while (rIdx < refLen) {
        for (auto n: refNodes) {
            int sIdx = n.first;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - storeFrom;
            util->seqBuf[storeTo][sIdx][alnIdx] = util->seqBuf[storeFrom][sIdx][rIdx];
        }
        ++alnIdx;++rIdx;
    }
    while (qIdx < qryLen) {
        for (auto n: refNodes) {
            int sIdx = n.first;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - storeFrom;
            util->seqBuf[storeTo][sIdx][alnIdx] = '-';
        }
        ++alnIdx;++qIdx;
    }
    for (auto n: refNodes) {
        int sIdx = n.first;
        util->changeStorage(sIdx);
    }
    return;
}

void merger_qry(Tree* tree, std::vector<Node*>& qryNodes, msa::utility* util, std::string& refString, std::string& qryString, int qLevel) {
    int rIdx = 0, qIdx = 0, alnIdx = 0;
    int refLen = refString.size(), qryLen = qryString.size();
    // std::cout << refString << '\n' << qryString << '\n';
    if (removeGap(refString) != removeGap(qryString)) {
        std::cout << "XXXXXX\n";
    }
    std::map<int, int> qrySeqIdx;
    for (auto n: qryNodes) qrySeqIdx[util->seqsIdx[n->longestDescendant->identifier]] = util->seqsStorage[util->seqsIdx[n->longestDescendant->identifier]];
    // for (auto n: qrySeqIdx) std::cout << n.first << ',';
    // std::cout << '\n';
    assert (refLen >= qryLen);
    for (int i = 0; i < refLen; ++i) {
        if (refString[i] == qryString[qIdx]) {
            for (auto n: qryNodes) {
                int sIdx = util->seqsIdx[n->longestDescendant->identifier];
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - storeFrom;
                // util->seqBuf[storeTo][sIdx][i] = util->seqBuf[storeFrom][sIdx][qIdx];
                util->seqBuf[storeTo][sIdx][i] = n->longestDescendant->msaSeq[qLevel][qIdx];
            }
            ++qIdx;
        }
        else {
            for (auto n: qryNodes) {
                int sIdx = util->seqsIdx[n->longestDescendant->identifier];
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - storeFrom;
                util->seqBuf[storeTo][sIdx][i] = '-';
            }
        }
    }   
    for (auto n: qryNodes) {
        int sIdx = util->seqsIdx[n->longestDescendant->identifier];
        util->changeStorage(sIdx);
    }
    return;
}

void transitivityMerge_regressive(Tree* tree, std::map<int, std::vector<std::vector<Node*>>>& subMSAs, msa::utility* util) {

    auto mergeStart = std::chrono::high_resolution_clock::now();
    // clear seq storage
    for (int sIdx = 0; sIdx < util->memNum; ++sIdx) {
        util->seqsStorage[sIdx] = 0;
        for (int i = 0; i < util->memLen; ++i) {
            util->seqBuf[0][sIdx][i] = 0;
            util->seqBuf[1][sIdx][i] = 0;
        }
    }
    // store MSA to storage
    assert(subMSAs[0].size() == 1);
    for (auto n: tree->allNodes) {
        if (n.second->is_leaf()) {
            int sIdx = util->seqsIdx[n.second->longestDescendant->identifier];
            std::string seq = tree->allNodes[n.second->longestDescendant->identifier]->msaSeq.begin()->second;
            for (int i = 0; i < seq.size(); ++i) {
                util->seqBuf[0][sIdx][i] = seq[i];
            }
        }
    }
    
    // for (int i = 0; i < util->memNum; ++i) {
    //     int storage = util->seqsStorage[i];
    //     int s = 0;
    //     while (true) {
    //         if (util->seqBuf[storage][i][s] == 0) break;
    //         ++s;
    //     }
    //     std::cout << i << ':' << s << '\n';
    // }
    std::map<int, Node*> refNodes;
    
    for (auto n: subMSAs[0][0]) {
        int longestIdx = util->seqsIdx[n->longestDescendant->identifier];
        if (refNodes.find(longestIdx) == refNodes.end()) {
            refNodes[longestIdx] = n->longestDescendant;
        }
    }
    for (int i = 1; i < subMSAs.size(); ++i) {
        auto parentMSAs = subMSAs[i-1];
        int childIdx = 0;
        int parentIdx = 0;
        
        for (auto pMSA: parentMSAs) {
            for (auto node: pMSA) {
                auto longestD = node->longestDescendant;
                while (childIdx < subMSAs[i].size()) {
                    auto cMSA = subMSAs[i][childIdx];
                    bool is_child = false;
                    for (auto cNode: cMSA) {
                        if (cNode->longestDescendant->identifier == longestD->identifier) {
                            is_child = true;
                            break;
                        }
                    }
                    if (is_child) {
                        printf("Merge Ref: %d-%d (qry) to %d-%d (ref).\n", i, childIdx, i-1, parentIdx); 
                        std::string refString = "";
                        int sIdx = util->seqsIdx[longestD->identifier];
                        int storage = util->seqsStorage[sIdx];
                        int s = 0;
                        while (true) {
                            if (util->seqBuf[storage][sIdx][s] == 0) break;
                            refString += util->seqBuf[storage][sIdx][s];
                            ++s;
                        }
                        std::string qryString = tree->allNodes[longestD->identifier]->msaSeq[i];
                        merger_ref(tree, refNodes, util, refString, qryString, i);
                        ++childIdx;
                        
                        // debug
                        // std::string refString_post = "";
                        // sIdx = util->seqsIdx[longestD->identifier];
                        // storage = util->seqsStorage[sIdx];
                        // s = 0;
                        // while (true) {
                        //     if (util->seqBuf[storage][sIdx][s] == 0) break;
                        //     refString_post += util->seqBuf[storage][sIdx][s];
                        //     ++s;
                        // }
                        // std::string rawr = "", rawq = "";
                        // for (int k = 0; k < refString.size(); ++k) if (refString[k] != '-') rawr += refString[k];
                        // for (int k = 0; k < refString_post.size(); ++k) if (refString_post[k] != '-') rawq += refString_post[k];
                        // if (rawr != rawq) {
                            
                        //     std::cout << "Post: Unmatched Seq.\n";
                        //     std::cout << refString << '\n' << qryString << '\n';
                        //     std::cout << refString_post << '\n';
                        //     std::cout << rawr << '\n' << rawq << '\n';
                        //     exit(1);
                        // }
                    }
                    else break;
                }
            }
            ++parentIdx;
        }

        childIdx = 0; parentIdx = 0;
        
        for (auto pMSA: parentMSAs) {
            for (auto node: pMSA) {
                auto longestD = node->longestDescendant;
                while (childIdx < subMSAs[i].size()) {
                    auto cMSA = subMSAs[i][childIdx];
                    bool is_child = false;
                    for (auto cNode: cMSA) {
                        if (cNode->longestDescendant->identifier == longestD->identifier) {
                            is_child = true;
                            break;
                        }
                    }
                    if (is_child) {
                        printf("Merge Qry: %d-%d (qry) to %d-%d (ref).\n", i, childIdx, i-1, parentIdx); 
                        std::string refString = "";
                        int sIdx = util->seqsIdx[longestD->identifier];
                        int storage = util->seqsStorage[sIdx];
                        int s = 0;
                        while (true) {
                            if (util->seqBuf[storage][sIdx][s] == 0) break;
                            refString += util->seqBuf[storage][sIdx][s];
                            ++s;
                        }
                        std::string qryString = tree->allNodes[longestD->identifier]->msaSeq[i];
                        merger_qry(tree, cMSA, util, refString, qryString, i);
                        ++childIdx;
                        // debug
                        // std::string refString_post = "";
                        // sIdx = util->seqsIdx[longestD->identifier];
                        // storage = util->seqsStorage[sIdx];
                        // s = 0;
                        // while (true) {
                        //     if (util->seqBuf[storage][sIdx][s] == 0) break;
                        //     refString_post += util->seqBuf[storage][sIdx][s];
                        //     ++s;
                        // }
                        // std::string rawr = "", rawq = "";
                        // for (int k = 0; k < refString.size(); ++k) if (refString[k] != '-') rawr += refString[k];
                        // for (int k = 0; k < refString_post.size(); ++k) if (refString_post[k] != '-') rawq += refString_post[k];
                        // if (rawr != rawq) {
                            
                        //     std::cout << "PostMergeQ: Unmatched Seq.\n";
                        //     std::cout << refString << '\n' << qryString << '\n';
                        //     std::cout << refString_post << '\n';
                        //     exit(1);
                        // } 
                        for (auto n: cMSA) {
                            int longestIdx = util->seqsIdx[n->longestDescendant->identifier];
                            if (refNodes.find(longestIdx) == refNodes.end()) {
                                refNodes[longestIdx] = n->longestDescendant;
                            }
                        }
                    }
                    else break;
                }
            }
            ++parentIdx;
        }
    }
    auto mergeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
    printf("Merge time : %d us\n",  (mergeTime.count() / 1000));      
    return;
}

void msaPostOrderTraversal_multigpu_regressive(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{
    auto freqStart = std::chrono::high_resolution_clock::now();
    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        if (n.first->msaIdx.size() == 0) {
            n.first->msaIdx.push_back(util->seqsIdx[n.first->longestDescendant->identifier]);
            n.first->msa.push_back(n.first->longestDescendant->identifier);
            util->seqsLen[n.first->identifier] = util->seqsLen[n.first->longestDescendant->identifier];
        }
        if (n.second->msaIdx.size() == 0) {
            n.second->msaIdx.push_back(util->seqsIdx[n.second->longestDescendant->identifier]);
            n.second->msa.push_back(n.second->longestDescendant->identifier);
            util->seqsLen[n.second->identifier] = util->seqsLen[n.second->longestDescendant->identifier];
        }
        // std::cout << n.first->identifier << ':' << n.second->identifier << '\n';
        // for (auto id: n.first->msaIdx) std::cout << id << ',';
        // std::cout << '\n';
        // for (auto id: n.second->msaIdx) std::cout << id << ',';
        // std::cout << '\n';
    }

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum;
    cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;

    int* alignSize = new int[roundGPU];
    int32_t* seqNum = new int32_t[roundGPU];
    uint16_t** hostFreq = new uint16_t* [roundGPU];
    int8_t**   hostAln = new int8_t* [roundGPU];
    int32_t**  hostLen = new int32_t* [roundGPU];
    int32_t**  hostAlnLen = new int32_t* [roundGPU];
    int32_t**  hostSeqInfo = new int32_t* [roundGPU];
    
    paramType* hostParam = (paramType*)malloc(28 * sizeof(paramType)); 

    if (param.scoreMode == 0) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 5 || j == 5)          hostParam[i*5+j] = 0;
                else if (i == j)               hostParam[i*5+j] = param.match;
                else if (i-j == 2 || j-i == 2) hostParam[i*5+j] = param.trans;
                else                           hostParam[i*5+j] = param.mismatch;
            }
        }
        hostParam[25] = param.gapOpen;
        hostParam[26] = param.gapExtend;
        hostParam[27] = param.xdrop;
    }
    else if (param.scoreMode == 1) {
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.hoxd70[i][j];
        hostParam[25] = param.hoxd70_gapOpen;
        hostParam[26] = param.hoxd70_gapExtend;
        hostParam[27] = param.xdrop;
    }
    
    
    std::vector<std::vector<uint16_t*>> freq;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> len;
    for (int rn = 0; rn < roundGPU; ++rn) {
        int pairsLeft = nodes.size() - rn*numBlocks;
        if (pairsLeft < numBlocks) alignSize[rn] = pairsLeft;
        else alignSize[rn] = numBlocks;
        seqNum[rn] = 0;
        hostFreq[rn] = (uint16_t*)malloc(12*seqLen * alignSize[rn] * sizeof(uint16_t));
        hostAln[rn] = (int8_t*)malloc(2*seqLen * alignSize[rn] * sizeof(int8_t));
        hostLen[rn] = (int32_t*)malloc(2*alignSize[rn] * sizeof(int32_t));
        hostAlnLen[rn] = (int32_t*)malloc(alignSize[rn] * sizeof(int32_t));
        hostSeqInfo[rn] = (int32_t*)malloc(5 * sizeof(int32_t));
        // store all sequences to array
        std::vector<uint16_t*> freqTemp;
        std::vector<std::pair<int32_t, int32_t>> seqIdxTemp;
        std::vector<std::pair<int32_t, int32_t>> lenTemp;
        for (int n = 0; n < alignSize[rn]; ++n) {
            int32_t nIdx = n + rn*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            refIdx = seqNum[rn];
            uint16_t *temp = new uint16_t[12*seqLen]; 
            for (int i = 0; i < 12*seqLen; ++i) temp[i]=0;
            // assert(temp.size() == 12*seqLen);
            // tbb::blocked_range<int> rangeRef(0, refLen);
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][sIdx][s] == 'A' || util->seqBuf[storage][sIdx][s] == 'a') temp[6*s+0]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'C' || util->seqBuf[storage][sIdx][s] == 'c') temp[6*s+1]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'G' || util->seqBuf[storage][sIdx][s] == 'g') temp[6*s+2]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'T' || util->seqBuf[storage][sIdx][s] == 't' ||
                             util->seqBuf[storage][sIdx][s] == 'U' || util->seqBuf[storage][sIdx][s] == 'u') temp[6*s+3]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'N' || util->seqBuf[storage][sIdx][s] == 'n') temp[6*s+4]+=1;
                    else                                                                                     temp[6*s+5]+=1;
                }
                });
                seqNum[rn] += 1;
            }
            qryIdx = seqNum[rn];
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][sIdx][s] == 'A' || util->seqBuf[storage][sIdx][s] == 'a') temp[6*(seqLen+s)+0]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'C' || util->seqBuf[storage][sIdx][s] == 'c') temp[6*(seqLen+s)+1]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'G' || util->seqBuf[storage][sIdx][s] == 'g') temp[6*(seqLen+s)+2]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'T' || util->seqBuf[storage][sIdx][s] == 't' ||
                             util->seqBuf[storage][sIdx][s] == 'U' || util->seqBuf[storage][sIdx][s] == 'u') temp[6*(seqLen+s)+3]+=1;
                    else if (util->seqBuf[storage][sIdx][s] == 'N' || util->seqBuf[storage][sIdx][s] == 'n') temp[6*(seqLen+s)+4]+=1;
                    else                                                                     temp[6*(seqLen+s)+5]+=1;
                }
                });
                seqNum[rn] += 1;
            }
            seqIdxTemp.push_back(std::make_pair(refIdx, qryIdx));
            lenTemp.push_back(std::make_pair(refLen, qryLen));
            freqTemp.push_back(temp);
        }
        freq.push_back(freqTemp);
        len.push_back(lenTemp);
        seqIdx.push_back(seqIdxTemp);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, roundGPU), [&](tbb::blocked_range<int> range){
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            for (int j = 0; j < 2*alignSize[gn]; ++j) { 
                if (j%2 == 0) hostLen[gn][j] = len[gn][j/2].first;
                else          hostLen[gn][j] = len[gn][j/2].second;
            }
            for (int j = 0; j < alignSize[gn]; ++j) {
                for (int l = 0; l < 12*seqLen; ++l) {
                    hostFreq[gn][12*seqLen*j+l] = freq[gn][j][l];
                }
            }
            for (int j = 0; j < 2*seqLen*alignSize[gn]; ++j) { 
                hostAln[gn][j] = 0;
            }
            for (int j = 0; j < alignSize[gn]; ++j) { 
                hostAlnLen[gn][j] = 0;
            }
            hostSeqInfo[gn][0] = seqLen;
            hostSeqInfo[gn][1] = seqNum[gn];
            hostSeqInfo[gn][2] = alignSize[gn];
            hostSeqInfo[gn][3] = numBlocks;
            hostSeqInfo[gn][4] = param.scoreMode;
        }
    });

    auto freqEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds freqTime = freqEnd -freqStart;
    printf("Preprocessing time : %d ms\n",  (freqTime.count() / 1000000));        

    auto kernelStart = std::chrono::high_resolution_clock::now();
    uint16_t** deviceFreq = new uint16_t* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];
    // int nowRound = 0;
    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            cudaSetDevice(gn);
            int nowMemSize = alignSize[gn];
            // cudaError_t error;
            cudaMalloc((void**)&deviceFreq[gn], 12*seqLen * alignSize[gn] * sizeof(uint16_t));
            // error = cudaGetLastError(); printf("CUDA error Freq1: %s, %d\n, E",cudaGetErrorString(error), gn); 
            cudaMalloc((void**)&deviceAln[gn], 2*seqLen * alignSize[gn] * sizeof(int8_t));
            // error = cudaGetLastError(); printf("CUDA error Freq2: %s\n",cudaGetErrorString(error)); 
            cudaMalloc((void**)&deviceLen[gn], 2*alignSize[gn] * sizeof(int32_t));
            // error = cudaGetLastError(); printf("CUDA error Freq3: %s\n",cudaGetErrorString(error)); 
            cudaMalloc((void**)&deviceAlnLen[gn], alignSize[gn] * sizeof(int32_t));
            // error = cudaGetLastError(); printf("CUDA error Freq4: %s\n",cudaGetErrorString(error)); 
            cudaMalloc((void**)&deviceSeqInfo[gn], 5 * sizeof(int32_t));
            // error = cudaGetLastError(); printf("CUDA error Freq5: %s\n",cudaGetErrorString(error)); 
            cudaMalloc((void**)&deviceParam[gn], 28 * sizeof(paramType));
            cudaMemcpy(deviceParam[gn], hostParam, 28 * sizeof(paramType), cudaMemcpyHostToDevice);
            // error = cudaGetLastError(); printf("CUDA error Freq6: %s\n",cudaGetErrorString(error)); 
                
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                if (alignSize[rn] != nowMemSize) {
                    // cudaSetDevice(gn);
                    cudaFree(deviceFreq[gn]);
                    cudaFree(deviceAln[gn]);
                    cudaFree(deviceLen[gn]);
                    cudaFree(deviceAlnLen[gn]);
                    // cudaFree(deviceParam[gn]);
                    // cudaFree(deviceSeqInfo[gn]);
                    // error = cudaGetLastError(); printf("CUDA error Free: %s\n",cudaGetErrorString(error)); 
                    cudaDeviceSynchronize();
                    cudaMalloc((void**)&deviceFreq[gn], 12*seqLen*alignSize[rn] * sizeof(uint16_t));
                    cudaMalloc((void**)&deviceAln[gn], 2*seqLen*alignSize[rn] * sizeof(int8_t));
                    cudaMalloc((void**)&deviceLen[gn], 2*alignSize[rn] * sizeof(int32_t));
                    cudaMalloc((void**)&deviceAlnLen[gn], alignSize[rn] * sizeof(int32_t));
                    // cudaMalloc((void**)&deviceSeqInfo[gn], 5 * sizeof(int32_t));
                    // error = cudaGetLastError(); printf("CUDA error Alloc: %s\n",cudaGetErrorString(error)); 
                }
                
                cudaMemcpy(deviceFreq[gn], hostFreq[rn], 12*seqLen * alignSize[rn] * sizeof(uint16_t), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error Freq: %s\n",cudaGetErrorString(error)); 
                cudaMemcpy(deviceAln[gn], hostAln[rn], 2*seqLen * alignSize[rn] * sizeof(int8_t), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error Aln: %s\n",cudaGetErrorString(error)); 
                cudaMemcpy(deviceLen[gn], hostLen[rn], 2*alignSize[rn] * sizeof(int32_t), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error Len: %s\n",cudaGetErrorString(error)); 
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[rn], alignSize[rn] * sizeof(int32_t), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error AlnLen: %s\n",cudaGetErrorString(error)); 
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[rn], 5 * sizeof(int32_t), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error SeqInfo: %s\n",cudaGetErrorString(error)); 
                // cudaMemcpy(deviceParam[gn], hostParam[rn], 7 * sizeof(paramType), cudaMemcpyHostToDevice);
                // error = cudaGetLastError(); printf("CUDA error Param: %s\n",cudaGetErrorString(error)); 
                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceParam[gn]
                );
                cudaDeviceSynchronize();
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
                cudaMemcpy(hostAln[rn], deviceAln[gn], 2*seqLen * alignSize[rn] * sizeof(int8_t), cudaMemcpyDeviceToHost);
                // error = cudaGetLastError(); printf("CUDA error rAln: %s\n",cudaGetErrorString(error)); 
                cudaMemcpy(hostAlnLen[rn], deviceAlnLen[gn], alignSize[rn] * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();  
            }
        }
    });

    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(gn);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaDeviceSynchronize();  
    }

    auto kernelEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    int totalPairs = 0;
    for (int gn = 0; gn < roundGPU; ++gn) totalPairs += alignSize[gn];
    std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";

    auto reAlnStart = std::chrono::high_resolution_clock::now();
    
    int maxAlnLen = 0;
    for (int gn = 0; gn < roundGPU; ++gn) {
       for (int k = 0; k <  alignSize[gn]; ++k) {
            if (hostAlnLen[gn][k] > maxAlnLen) maxAlnLen = hostAlnLen[gn][k];
        }
    }
    util->memCheck(maxAlnLen);
        
    for (int gn = 0; gn < roundGPU; ++gn) {
        if (alignSize[gn] == 0) break;
        tbb::parallel_for(tbb::blocked_range<int>(0, alignSize[gn]), [&](tbb::blocked_range<int> range) {
            // for (int k = 0; k < alignSize[gn]; ++k) {
            for (int k = range.begin(); k < range.end(); ++k) {
                // std::vector<std::string> alignment;
                int32_t refNum = seqIdx[gn][k].second - seqIdx[gn][k].first;
                int32_t qryNum = (k !=  alignSize[gn]-1) ? seqIdx[gn][k+1].first - seqIdx[gn][k].second : seqNum[gn] - seqIdx[gn][k].second;
                int32_t refStart = seqIdx[gn][k].first;
                int32_t qryStart = seqIdx[gn][k].second;
                int32_t nIdx = k + gn*numBlocks;
                if (hostAlnLen[gn][k] <= 0) {
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    std::vector<int8_t> aln;
                    alignGrpToGrp_traditional
                    (
                        freq[gn][k],
                        seqLen,
                        refLen,
                        qryLen,
                        param,
                        aln
                    );
                    int32_t alnLen = aln.size();
                    util->memCheck(alnLen);
                    std::reverse(aln.begin(), aln.end());
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        int orgIdx = 0;
                        for (int j = 0; j < aln.size(); ++j) {
                            if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                                util->seqBuf[storeTo][sIdx][j] = util->seqBuf[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                        util->changeStorage(sIdx);
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        int orgIdx = 0;
                        for (int j = 0; j < aln.size(); ++j) {
                            if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                                util->seqBuf[storeTo][sIdx][j] = util->seqBuf[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                        util->changeStorage(sIdx);
                    }
                    printf("CPU fallback (traditional global alignment) on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                    // printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                }
                // else if (hostAlnLen[gn][k] <= 0) {
                //     std::vector<int8_t> aln;
                //     std::vector<std::vector<int>> freqRef;
                //     std::vector<std::vector<int>> freqQry;
                //     int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                //     int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                //     for (int r = 0; r < refLen; r++) {
                //         std::vector<int> temp;
                //         for (int f = 0; f < 6; ++f) temp.push_back(freq[gn][k][6*r+f]);
                //         freqRef.push_back(temp);
                //     }
                //     for (int q = 0; q < qryLen; q++) {
                //         std::vector<int> temp;
                //         for (int f = 0; f < 6; ++f) temp.push_back(freq[gn][k][6*(seqLen+q)+f]);
                //         freqQry.push_back(temp);
                //     }
                //     Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, 1000, param.marker);
                //     Talco_xdrop::Align_freq (
                //         talco_params,
                //         freqRef,
                //         freqQry,
                //         aln
                //     );
                //     util->memCheck(aln.size());
                //     for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                //         int64_t start = sIdx*util->memLen;
                //         int storeFrom = util->seqsStorage[sIdx];
                //         int storeTo = 1 - util->seqsStorage[sIdx];
                //         int orgIdx = 0;
                //         for (int j = 0; j < aln.size(); ++j) {
                //             if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                //                 util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                //                 orgIdx++;
                //             }
                //             else {
                //                 util->seqBuf[storeTo][start+j] = '-';
                //             }
                //         }
                //         util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                //         util->changeStorage(sIdx);
                //     }
                //     for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                //         int64_t start = sIdx*util->memLen;
                //         int storeFrom = util->seqsStorage[sIdx];
                //         int storeTo = 1 - util->seqsStorage[sIdx];
                //         int orgIdx = 0;
                //         for (int j = 0; j < aln.size(); ++j) {
                //             if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                //                 util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                //                 orgIdx++;
                //             }
                //             else {
                //                 util->seqBuf[storeTo][start+j] = '-';
                //             }
                //         }
                //         util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                //         util->changeStorage(sIdx);
                //     }
                //     printf("CPU fallback (TALCO-Xdrop) on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                // }
                else {
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                        int orgIdx = 0;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 2) {
                                util->seqBuf[storeTo][sIdx][j] = util->seqBuf[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        int orgIdx = 0;
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 1) {
                                util->seqBuf[storeTo][sIdx][j] = util->seqBuf[storeFrom][sIdx][orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][sIdx][j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                }
                // std::cout << "LenB : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
                //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
                for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) 
                    tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msa) 
                    tree->allNodes[nodes[nIdx].first->identifier]->msa.push_back(q);
            }  
        });
        for (int i = 0; i < alignSize[gn]; ++i) delete [] freq[gn][i];
    } 
    auto reAlnEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds reAlnTime = reAlnEnd - kernelEnd;
    printf("Alignment Time: %d us\n", reAlnTime.count() / 1000);

    

    for (int rn = 0; rn < roundGPU; ++rn) {
        free(hostFreq[rn]);
        free(hostAlnLen[rn]);
        free(hostLen[rn]);
        free(hostAln[rn]);
        free(hostSeqInfo[rn]);
    }  
    free(hostParam);

    delete [] alignSize;
    delete [] seqNum;
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostAln;
    // delete [] hostParam;
    delete [] hostSeqInfo;
    return;
}

void createOverlapMSA_subtree(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{

    int numBlocks = 1024; 
    int blockSize = THREAD_NUM;
    int gpuNum = util->gpuNum;
    // cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    // get maximum sequence/profile length 
    int32_t seqLen = util->memLen;
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;

    paramType* hostParam = (paramType*)malloc(28 * sizeof(paramType)); 

    if (!param.userDefine) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 5 || j == 5)          hostParam[i*5+j] = 0;
                else if (i == j)               hostParam[i*5+j] = param.match;
                else if (i-j == 2 || j-i == 2) hostParam[i*5+j] = param.trans;
                else                           hostParam[i*5+j] = param.mismatch;
            }
        }
        hostParam[25] = param.gapOpen;
        hostParam[26] = param.gapExtend;
        hostParam[27] = param.xdrop;
    }
    else {
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) hostParam[i*5+j] = param.userMatrix[i][j];
        hostParam[25] = param.userGapOpen;
        hostParam[26] = param.userGapExtend;
        hostParam[27] = param.xdrop;
    }

    std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    
    uint16_t** hostFreq = new uint16_t* [gpuNum];
    int8_t**   hostAln = new int8_t* [gpuNum];
    int32_t**  hostLen = new int32_t* [gpuNum];
    int32_t**  hostAlnLen = new int32_t* [gpuNum];
    int32_t**  hostSeqInfo = new int32_t* [gpuNum];

    uint16_t** deviceFreq = new uint16_t* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    paramType**  deviceParam = new paramType* [gpuNum];
  
    std::atomic<int> nowRound;
    nowRound.store(0);

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            hostFreq[gn] = (uint16_t*)malloc(12 * seqLen * numBlocks * sizeof(uint16_t));
            hostAln[gn] = (int8_t*)malloc(    2 * seqLen * numBlocks * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(   2 *          numBlocks * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(             numBlocks * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(5 * sizeof(int32_t));
            
            cudaSetDevice(gn);
            // cudaError_t error;
            cudaMalloc((void**)&deviceFreq[gn],  12 * seqLen * numBlocks * sizeof(uint16_t));
            cudaMalloc((void**)&deviceAln[gn],    2 * seqLen * numBlocks * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn],    2 *          numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn],              numBlocks * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 5 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn],  28 * sizeof(paramType));

            cudaMemcpy(deviceParam[gn], hostParam, 28 * sizeof(paramType), cudaMemcpyHostToDevice);
            // error = cudaGetLastError(); printf("CUDA error Malloc: %s\n",cudaGetErrorString(error)); 
            std::vector<std::pair<int, int>> seqIdx;
            
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                int alnPairs = (nodes.size() - rn*numBlocks > numBlocks) ? numBlocks : nodes.size() - rn*numBlocks;
                int seqNum = 0;
                // std::cout << "GPU: " << gn << " Rn: " << rn << " Pairs: " << alnPairs << '\n';
                
                // Initailize 
                for (int n = 0; n < 12*seqLen * numBlocks; ++n) hostFreq[gn][n] = 0;
                for (int n = 0; n <  2*seqLen * numBlocks; ++n) hostAln[gn][n] = 0;
                for (int n = 0; n <  2*         numBlocks; ++n) hostLen[gn][n] = 0;
                for (int n = 0; n <             numBlocks; ++n) hostAlnLen[gn][n] = 0;
                seqIdx.clear();

                // Calculate Frequency
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t nIdx = n + rn*numBlocks;
                    int32_t qryIdx = 0;
                    int32_t refIdx = 0;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    // std::cout << n << "Len: " << refLen << ',' << qryLen << '\n';
                    refIdx = seqNum;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < refLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*s+0]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*s+1]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*s+2]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*s+3]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*s+4]+=1;
                            else                                                                                             hostFreq[gn][12*seqLen*n+6*s+5]+=1;
                        }
                        });
                        seqNum += 1;
                    }
                    qryIdx = seqNum;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                        int storage = util->seqsStorage[sIdx];
                        tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                        for (int s = r.begin(); s < r.end(); ++s) {
                        // for (int s = 0; s < qryLen; ++s) {
                            if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+0]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+1]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+2]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                                     util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+3]+=1;
                            else if (util->alnStorage[storage][sIdx][s] == 'N' || util->alnStorage[storage][sIdx][s] == 'n') hostFreq[gn][12*seqLen*n+6*(seqLen+s)+4]+=1;
                            else                                                                                             hostFreq[gn][12*seqLen*n+6*(seqLen+s)+5]+=1;
                        }
                        });
                        seqNum += 1;
                    }
                    hostLen[gn][2*n] = refLen; hostLen[gn][2*n+1] = qryLen;
                    seqIdx.push_back(std::make_pair(refIdx, qryIdx));
                }

                hostSeqInfo[gn][0] = seqLen;
                hostSeqInfo[gn][1] = seqNum;
                hostSeqInfo[gn][2] = alnPairs;
                hostSeqInfo[gn][3] = numBlocks;
                hostSeqInfo[gn][4] = param.userDefine;
        
                cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * numBlocks * sizeof(uint16_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[gn], 2*numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 5 * sizeof(int32_t), cudaMemcpyHostToDevice);
                
                std::string berr = cudaGetErrorString(cudaGetLastError());
                if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceParam[gn]
                );
                cudaDeviceSynchronize();
                std::string aerr = cudaGetErrorString(cudaGetLastError());
                if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
                
                cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * numBlocks * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], numBlocks * sizeof(int32_t), cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
                // int maxAlnLen = 0;
                // for (int n = 0; n <  alnPairs; ++n) {
                //     if (hostAlnLen[gn][n] > maxAlnLen) maxAlnLen = hostAlnLen[gn][n];
                // }
                // util->memCheck(maxAlnLen);
                
                
                // tbb::parallel_for(tbb::blocked_range<int>(0, alignSize[gn]), [&](tbb::blocked_range<int> range) {
                // for (int k = range.begin(); k < range.end(); ++k) {
                for (int n = 0; n < alnPairs; ++n) {
                    int32_t refNum = seqIdx[n].second - seqIdx[n].first;
                    int32_t qryNum = (n !=  alnPairs-1) ? seqIdx[n+1].first - seqIdx[n].second : seqNum - seqIdx[n].second;
                    int32_t nIdx = n + rn*numBlocks;

                    if (hostAlnLen[gn][n] <= 0) {
                        int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                        int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                        uint16_t *freq = new uint16_t[12*seqLen]; 
                        for (int i = 0; i < 12*seqLen; ++i) freq[i] = hostFreq[gn][12*seqLen*n+i];
                        std::vector<int8_t> aln;
                        alignGrpToGrp_traditional (
                            freq,
                            seqLen,
                            refLen,
                            qryLen,
                            param,
                            aln
                        );
                        delete [] freq;
                        int32_t alnLen = aln.size();
                        util->memCheck(alnLen);
                        std::reverse(aln.begin(), aln.end());
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                        std::cout << "CPU fallback (traditional global alignment) on No. " << n << " (" << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ")\n";
                    }
                    else {
                        std::vector<int8_t> aln;
                        for (int j = 0; j < hostAlnLen[gn][n]; ++j) {
                            aln.push_back(hostAln[gn][n*2*seqLen+j]);
                        }
                        tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                    }
                }    
            }  
        }
    });
    
    for (auto n: nodes) {
        tree->allNodes[n.first->identifier]->msa.clear();
        tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
        tree->allNodes[n.second->identifier]->msa.clear();
        tree->allNodes[n.second->identifier]->msa.push_back(n.second->identifier);    
    }

    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(gn);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceLen[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaFree(deviceParam[gn]);
        cudaDeviceSynchronize();  
        free(hostFreq[gn]);
        free(hostAlnLen[gn]);
        free(hostLen[gn]);
        free(hostAln[gn]);
        free(hostSeqInfo[gn]);
    }
    
    free(hostParam);

    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostLen;
    delete [] hostAln;
    delete [] hostSeqInfo;
    return;
}

*/
