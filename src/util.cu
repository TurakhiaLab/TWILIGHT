#ifndef UTIL_HPP
#include "util.cuh"
#endif

KSEQ_INIT2(, gzFile, gzread);
// po::options_description mainDesc("MSA Command Line Arguments");


// Set variables
Params* setParameters(po::variables_map& vm) {
    
    int userDefine = vm["scoring-matrix"].as<int>();
    Params* param = new Params();
    param->userDefine = userDefine;
    // Kimura
    paramType gapOp = vm["gap-open"].as<paramType>();
    paramType gapCl = gapOp;
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
    int maxSubtreeSize = (vm.count("max-subtree-size")) ? vm["max-subtree-size"].as<int>() : INT32_MAX;
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
    float gappyVertical = vm["gappy-vertical"].as<float>();
    float gappyHorizon;
    if (gappyVertical > 1 || gappyVertical <= 0) {
        std::cerr << "ERROR: the value of gappy-vertical should be (0,1]\n";
        exit(1);
    }
    if (vm.count("gappy-horizon")) {
        gappyHorizon = vm["gappy-horizon"].as<float>();
        if (gappyHorizon - round(gappyHorizon) != 0 || gappyHorizon < 1) {
            std::cerr << "ERROR: the value of gappy-horizon should be an positive integer\n";
            exit(1);
        }
    }
    else {
        gappyHorizon = 0;
    }
    std::string tempDir;
    if (!vm.count("temp-dir")) tempDir =  "./temp";
    else tempDir = vm["temp-dir"].as<std::string>();
    if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);

    std::string outType = vm["output-type"].as<std::string>();
    if (outType != "FASTA" && outType != "CIGAR") {
        std::cerr << "ERROR: Unrecognized output type \"" << outType <<"\"\n";
        exit(1);
    }
    
    option->cpuNum = cpuNum; 
    option->gpuNum = gpuNum;
    option->gpuIdx = gpuIdx;
    option->maxSubtree = maxSubtreeSize;
    option->maxSubSubtree = maxSubSubtreeSize;
    option->debug = vm.count("debug");
    option->gappyHorizon = gappyHorizon;
    option->gappyVertical = gappyVertical;
    option->tempDir = tempDir;
    option->cpuOnly = vm.count("cpu-only");
    option->outType = outType;
}

// print
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


// read
void readSequences(po::variables_map& vm, msa::utility* util, msa::option* option, Tree* tree)
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
        if (tree->allNodes.find(seqName) == tree->allNodes.end()) {
            fprintf(stderr, "Sequence %s is not presented in the tree file.\n", seqName.c_str());
            exit(1);
        }
        if (seqs.find(seqName) != seqs.end()) {
            fprintf(stderr, "Sequence %s already exists.\n", seqName.c_str());
            exit(1);
        }
        int subtreeIdx = tree->allNodes[seqName]->grpID;
        std::string seq = std::string(kseq_rd->seq.s, seqLen);
        seqs[seqName] = std::make_pair(seq, subtreeIdx);
        if (option->debug) util->rawSeqs[seqName] = seq;
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

void readSequences(std::string seqFileName, msa::utility* util, msa::option* option, Tree* tree)
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
        std::string seq = std::string(kseq_rd->seq.s, seqLen);
        seqs[seqName] = std::make_pair(seq, subtreeIdx);
        if (option->debug) util->rawSeqs[seqName] = seq;
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

void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string tempDir;
    if (!vm.count("temp-dir")) tempDir =  "./temp";
    else tempDir = vm["temp-dir"].as<std::string>();
    if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
    if (mkdir(tempDir.c_str(), 0777) == -1) {
        if( errno == EEXIST ) std::cout << tempDir << " already exists.\n";
        else { fprintf(stderr, "ERROR: cant create directory: %s\n", tempDir.c_str()); exit(1); }
    }
    else std::cout << tempDir << " created\n";
    // std::cout << "Total " <<  partition->partitionsRoot.size() << " subtrees.\n";
    std::string seqFileName = vm["sequences"].as<std::string>();
    int maxLen = 0, totalLen = 0, seqNum = 0, minLen = INT_MAX;
    bool storeRaw = false;
    for (auto subroot: partition->partitionsRoot) {
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        Tree* subT = new Tree(subroot.second.first);
        gzFile f_rd = gzopen(seqFileName.c_str(), "r");
        if (!f_rd) { fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str()); exit(1);}
        kseq_t* kseq_rd = kseq_init(f_rd);
        std::map<std::string, std::string> seqs;
        while (kseq_read(kseq_rd) >= 0) {
            size_t seqLen = kseq_rd->seq.l;
            std::string seqName = kseq_rd->name.s;
            if (!storeRaw && option->debug) util->rawSeqs[seqName] = std::string(kseq_rd->seq.s, seqLen);
            if (subT->allNodes.find(seqName) != subT->allNodes.end()) {
                seqs[seqName] = std::string(kseq_rd->seq.s, seqLen);
                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
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
        delete subT;
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        if (!util->rawSeqs.empty()) storeRaw = true;
    }

    uint32_t avgLen = totalLen/seqNum;
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    // std::cout << "Seqences: (Num, MaxLen, AvgLen) = (" << seqNum << ", " << maxLen << ", " << avgLen << ")\n";
    std::cout << "=== Sequence information ===\n";
    std::cout << "Number : " << seqNum << '\n';
    std::cout << "Max. Length: " << maxLen << '\n';
    std::cout << "Min. Length: " << minLen << '\n';
    std::cout << "Avg. Length: " << avgLen << '\n';
    std::cout << "============================\n";
    std::cout << "Created " << partition->partitionsRoot.size() << " subtree files in " << tempDir << " in " <<  seqReadTime.count() / 1000000 << " ms\n";
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
        std::string rawInput;
        int idx, seqNum, seqLen; 
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
        idx = numbers[0]; seqNum = numbers[1]; seqLen = numbers[2]; 
        assert(idx == subtree);
        numbers.clear();
        util->seqsLen[subroot.first] = seqLen;
        if (seqLen > util->seqLen) util->seqLen = seqLen;
        util->profileFreq[subtree] = std::vector<std::vector<float>> (seqLen, std::vector<float> (6, 0.0));
        for (int t = 0; t < 6; ++t) {
            int s = 0;
            getline(inputStream, rawInput);
            for (int j = 0; j < rawInput.size(); ++j) {
                if (rawInput[j] == ',') {
                    util->profileFreq[subtree][s][t] = std::atof(num.c_str());
                    num = "";
                    ++s;
                }
                else if (j == rawInput.size()-1) {
                    num += rawInput[j];
                    util->profileFreq[subtree][s][t] = std::atof(num.c_str());
                    num = "";
                    ++s;
                }
                else num += rawInput[j];
            }
            assert(s == seqLen);
        }
        inputStream.close();
    }
    return;
}

// write (output)

void outputFinal (std::string tempDir, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs) {
    util->seqsIdx.clear();
    util->alnSeqs.clear();
    for (auto id: tree->root->msa) {
        int subtree = tree->allNodes[id]->grpID;
        std::string tempAlnFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".temp.aln";
        gzFile f_rd = gzopen(tempAlnFile.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", tempAlnFile.c_str());
            exit(1);
        }
        // std::cout << tempAlnFile << '\n';
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
            // std::cout << util->alnSeqs.size() << ':' << seqName << '\n' << finalAln << '\n';
            util->alnSeqs[seqName] = finalAln;
        }
        std::string subtreeSeqFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".final.aln";
        outputSubtreeSeqs(subtreeSeqFile, seqs);    
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        // std::remove(tempAlnFile.c_str());
        totalSeqs += seqs.size();
        // uint16_t** freq = new uint16_t* [6];
        // for (int i = 0; i < 6; ++i) {
        //     freq[i] = new uint16_t [seqLen];
        //     for (int j = 0; j < seqLen; ++j) freq[i][j] = 0;
        // }
        // for (auto sIdx = seqs.begin(); sIdx != seqs.end(); ++sIdx) {
        //     for (int j = 0; j < seqLen; ++j) {
        //         if      (sIdx->second[j] == 'A' || sIdx->second[j] == 'a') freq[0][j]+=1;
        //         else if (sIdx->second[j] == 'C' || sIdx->second[j] == 'c') freq[1][j]+=1;
        //         else if (sIdx->second[j] == 'G' || sIdx->second[j] == 'g') freq[2][j]+=1;
        //         else if (sIdx->second[j] == 'T' || sIdx->second[j] == 't' ||
        //                  sIdx->second[j] == 'U' || sIdx->second[j] == 'u') freq[3][j]+=1;
        //         else if (sIdx->second[j] == 'N' || sIdx->second[j] == 'n') freq[4][j]+=1;
        //         else                                                       freq[5][j]+=1;
        //     }
        // }
        // std::string subtreeFreqFile = tempDir + '/' + "subtree-" + std::to_string(subtree) + ".final.freq.txt";
        // std::ofstream outFile(subtreeFreqFile);
        // if (!outFile) {
        //     fprintf(stderr, "ERROR: cant open file: %s\n", subtreeFreqFile.c_str());
        //     exit(1);
        // }
        // // Info subtreeIdx, seqNum, seqLen
        // outFile << subtree << ',' << seqs.size() << ',' << seqLen << '\n';
        // seqs.clear();
        // for (int i = 0; i < 6; ++i) {
        //     for (int j = 0; j < seqLen-1; ++j) {
        //         outFile << freq[i][j] << ',';
        //     }
        //     outFile << freq[i][seqLen-1] << '\n';
        // }
        // for (int i = 0; i < 6; ++i) delete [] freq[i];
        // delete [] freq;
        // outFile.close();
    }
    
    tree->root->msa.clear();
    return;
}

void outputAln(std::string fileName, msa::utility* util, msa::option* option, Tree* T, int grpID) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    std::vector<std::string> seqs;
    if (util->alnSeqs.empty())
        for (auto seq: util->seqsIdx)
            seqs.push_back(seq.first);
    else
        for (auto seq: util->alnSeqs)
            seqs.push_back(seq.first);
    size_t seqLen = util->seqsLen[T->root->identifier];
    std::cout << "seqLen: " << seqLen << '\n';
    std::sort(seqs.begin(), seqs.end(), cmp);
    if (util->alnSeqs.empty()) {
        for (int s = 0; s < seqs.size(); ++s) {
            int sIdx = util->seqsIdx[seqs[s]];
            int storage = util->seqsStorage[sIdx];
            if (std::find(T->root->msaIdx.begin(), T->root->msaIdx.end(), sIdx) != T->root->msaIdx.end()) {
                outFile << '>' << seqs[s] << "\n";
                if (option->outType == "FASTA") {
                    outFile.write(&util->alnStorage[storage][sIdx][0], seqLen);
                    outFile << '\n';
                }
                else if (option->outType == "CIGAR") {
                    std::string cigar = "";
                    char type = 0;
                    int num = 0;
                    int i = 0;
                    while (util->alnStorage[storage][sIdx][i] != 0) {
                        outFile << util->alnStorage[storage][sIdx][i];
                        ++i;
                    }
                }
                outFile << '\n';
            }
            util->seqFree(sIdx);
        }
        util->seqsFree();
    }
    else {
        for (int s = 0; s < seqs.size(); ++s) {
            outFile << ('>' + seqs[s] + '\n');
            outFile << (util->alnSeqs[seqs[s]] + '\n');
        }
    }
    outFile.close();
    
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
    // std::cout << "seqLen: " << seqLen << '\n';
    float** freq = new float* [6];
    for (int i = 0; i < 6; ++i) {
        freq[i] = new float [seqLen];
        for (int j = 0; j <  seqLen; ++j) freq[i][j] = 0;
    }
    float totalWeight = 0;
    for (auto sIdx: T->root->msaIdx) totalWeight += T->allNodes[util->seqsName[sIdx]]->weight;
    for (int sIdx: T->root->msaIdx) {
        int storage = util->seqsStorage[sIdx];
        std::string name = util->seqsName[sIdx];
        float w = T->allNodes[name]->weight / totalWeight * T->root->numLeaves;
        for (int j = 0; j <  seqLen; ++j) {
            if      (util->alnStorage[storage][sIdx][j] == 'A' || util->alnStorage[storage][sIdx][j] == 'a') freq[0][j]+=1*w;
            else if (util->alnStorage[storage][sIdx][j] == 'C' || util->alnStorage[storage][sIdx][j] == 'c') freq[1][j]+=1*w;
            else if (util->alnStorage[storage][sIdx][j] == 'G' || util->alnStorage[storage][sIdx][j] == 'g') freq[2][j]+=1*w;
            else if (util->alnStorage[storage][sIdx][j] == 'T' || util->alnStorage[storage][sIdx][j] == 't' ||
                     util->alnStorage[storage][sIdx][j] == 'U' || util->alnStorage[storage][sIdx][j] == 'u') freq[3][j]+=1*w;
            else if (util->alnStorage[storage][sIdx][j] == 'N' || util->alnStorage[storage][sIdx][j] == 'n') freq[4][j]+=1*w;
            else                                                                                             freq[5][j]+=1*w;
        }
    }
    for (int i = 0; i < 6; ++i) {
        std::string outString = "";
        for (int j = 0; j < seqLen-1; ++j) {
            outString += (std::to_string(freq[i][j]) + ',');
            // outFile << freq[i][j] << ',';
        }
        outString += (std::to_string(freq[i][seqLen-1]) + '\n');
        outFile << outString;
        // outFile << freq[i][seqLen-1] << '\n';
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

// auxiliary
bool cmp(std::string a, std::string b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
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

void getSubtreeNewick(Node* root, std::string& outputString) {
	if(root->children.size() != 0) {
		outputString += "(";
		for(int n = 0; n < root->children.size(); ++n) {
			if(n != 0) outputString += ",";
			getSubtreeNewick(root->children[n], outputString);
		}
		if (root->parent != nullptr) outputString += ("):" + std::to_string(root->branchLength));
        else outputString += ")";
    }
	else {
		outputString += (root->identifier + ':' + std::to_string(root->branchLength));
    }
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
