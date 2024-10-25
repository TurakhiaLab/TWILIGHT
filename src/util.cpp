#ifndef UTIL_HPP
#include "util.hpp"
#endif

KSEQ_INIT2(, gzFile, gzread);

// print
void printTree(Node* node, int grpID)
{
    if (grpID == -1) {
        if (node->parent == nullptr)
            std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << "ROOT\t" << node->grpID << '\t' << node->weight <<std::endl;
        else
            std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << "\t" << node->grpID << '\t' << node->weight << std::endl;  

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
        int seqLen = kseq_rd->seq.l;
        std::string seqName = kseq_rd->name.s;
        if (tree->allNodes.find(seqName) != tree->allNodes.end()) {
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
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);
        
    
    seqNum = seqs.size();
    uint32_t avgLen = totalLen/seqNum;
    if (option->printDetail) {
        std::cout << "=== Sequence information ===\n";
        std::cout << "Number : " << seqNum << '\n';
        std::cout << "Max. Length: " << maxLen << '\n';
        std::cout << "Min. Length: " << minLen << '\n';
        std::cout << "Avg. Length: " << avgLen << '\n';
        std::cout << "============================\n";
    }

    if (tree->m_numLeaves != seqNum) {
        fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %d)\n", tree->m_numLeaves, seqNum); 
        exit(1);
    }
    util->seqsMallocNStore(maxLen, seqs, option);
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
        int seqLen = kseq_rd->seq.l;
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

    util->seqsMallocNStore(maxLen, seqs, option);
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

/*
void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, paritionInfo_t* partition, msa::utility* util, msa::option* option)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string tempDir = option->tempDir;
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
            int seqLen = kseq_rd->seq.l;
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
        outputSubtree(subT, option, subtreeIdx);
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
*/

void readSequencesNoutputTemp(po::variables_map& vm, Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string tempDir = option->tempDir;
    for (auto subroot: partition->partitionsRoot) {
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        Tree* subT = new Tree(subroot.second.first);
        outputSubtree(subT, option, subtreeIdx);
    }
    return;
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

void readFreq(std::string tempDir, Tree* tree, partitionInfo_t* partition, msa::utility* util) {
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
        for (size_t i = 0; i < rawInput.size(); ++i) {
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
        util->seqsIdx[subroot.first] = subtree;
        util->seqsName[subtree] = subroot.first;
        tree->allNodes[subroot.first]->msaAln = std::vector<int8_t> (seqLen, 0);
        if (seqLen > util->seqLen) util->seqLen = seqLen;
        util->profileFreq[subtree] = std::vector<std::vector<float>> (seqLen, std::vector<float> (6, 0.0));
        for (int t = 0; t < 6; ++t) {
            int s = 0;
            getline(inputStream, rawInput);
            for (size_t j = 0; j < rawInput.size(); ++j) {
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

void updateSeqLen(Tree* tree, partitionInfo_t* partition, msa::utility* util) {
    for (auto subroot:  partition->partitionsRoot) {
        int subtree = tree->allNodes[subroot.first]->grpID;
        int seqLen = util->profileFreq[subtree].size();
        util->seqsLen[subroot.first] = seqLen;
        util->seqsIdx[subroot.first] = subtree;
        util->seqsName[subtree] = subroot.first;
        tree->allNodes[subroot.first]->msaAln = std::vector<int8_t> (seqLen, 0);
        if (seqLen > util->seqLen) util->seqLen = seqLen;
    }
    return;
}

void outputFinal (po::variables_map& vm, Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs) {
    util->seqsIdx.clear();
    std::string seqFileName = vm["sequences"].as<std::string>();
    std::cout << "Final Alignment Length: " << tree->allNodes[tree->root->identifier]->msaAln.size() << '\n';
    std::map<std::string, std::string> seqs;
    std::map<int, std::pair<std::string, std::string>> rawSeqs;
    tbb::spin_rw_mutex  writeMutex;
                            
    int proceeded = 0;    
    for (auto subroot: partition->partitionsRoot) {
        gzFile f_rd = gzopen(seqFileName.c_str(), "r");
        if (!f_rd) { fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str()); exit(1);}
        kseq_t* kseq_rd = kseq_init(f_rd);
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        ++proceeded;
        std::cout << "Start writing alignment of subtree No. " << subtreeIdx << ". (" << proceeded << '/' << partition->partitionsRoot.size() << ")\n";
        Tree* subT = new Tree(subroot.second.first);
        int s = 0;
        while (kseq_read(kseq_rd) >= 0) {
            std::string seqName = kseq_rd->name.s;
            if (subT->allNodes.find(seqName) != subT->allNodes.end()) {
                size_t seqLen = kseq_rd->seq.l;
                rawSeqs[s] = std::make_pair(seqName, std::string(kseq_rd->seq.s, seqLen));
                ++s;
            }
        }
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        assert(subT->m_numLeaves == rawSeqs.size());
        tbb::parallel_for(tbb::blocked_range<int>(0, rawSeqs.size()), [&](tbb::blocked_range<int> range){ 
        for (int n = range.begin(); n < range.end(); ++n) {
            std::string seqName = rawSeqs[n].first;
            std::string rawSeq = rawSeqs[n].second;
            std::string updatedSeq = "", alnSeq = "";
            size_t r = 0;
            std::string num_s = "";
            for (size_t i = 0; i < util->seqsCIGAR[seqName].size(); ++i) {
                if (util->seqsCIGAR[seqName][i] == 'M') {
                    int num = std::stoi(num_s.c_str());
                    num_s = "";
                    for (int k = 0; k < num; ++k) {updatedSeq += rawSeq[r]; ++r;}
                }
                else if (util->seqsCIGAR[seqName][i] == 'D') {
                    int num = std::stoi(num_s.c_str());
                    num_s = "";
                    for (int k = 0; k < num; ++k) updatedSeq += '-';
                }
                else {
                    num_s += util->seqsCIGAR[seqName][i];
                }
            }
            assert(r == rawSeq.size());            
            r = 0;
            for (size_t j = 0; j < tree->allNodes[subroot.first]->msaAln.size(); ++j) {
                if ((tree->allNodes[subroot.first]->msaAln[j] & 0xFFFF) == 0 || (tree->allNodes[subroot.first]->msaAln[j] & 0xFFFF) == 2) {
                    alnSeq += updatedSeq[r];
                    ++r;
                }
                else {
                    alnSeq += '-';
                }
            }
            {
                tbb::spin_rw_mutex::scoped_lock lock(writeMutex);
                if (alnSeq.size() != tree->allNodes[subroot.first]->msaAln.size()) {
                    std::cout << "ERROR: length not match. alnSeq (" << alnSeq.size() << ") != msaAln (" << tree->allNodes[subroot.first]->msaAln.size() << '\n'; 
                }
                assert(alnSeq.size() == tree->allNodes[subroot.first]->msaAln.size());
                seqs[seqName] = alnSeq;
            }
        }
        });
        tree->allNodes[subroot.first]->msaAln.clear();
        totalSeqs += seqs.size();
        std::string subtreeSeqFile = option->tempDir + '/' + "subtree-" + std::to_string(subtreeIdx) + ".final.aln";
        if (option->outType == "FASTA") outputSubtreeSeqs(subtreeSeqFile, seqs);   
        else if (option->outType == "CIGAR") outputSubtreeCIGAR(subtreeSeqFile, seqs);   
        delete subT;
        seqs.clear();
        rawSeqs.clear();
    }
    tree->root->msa.clear();
    return;
}

void outputAln(std::string fileName, msa::utility* util, msa::option* option, Tree* T) {
    if (util->nowProcess == 2) {
        std::string command = "cat " + option->tempDir + "/*.final.aln > " + fileName;
        system(command.c_str());
        if (option->deleteTemp) {
            command = "rm -rf " + option->tempDir;
            system(command.c_str());
        }
        return;
    }
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    std::vector<std::string> seqs;
    for (auto seq: util->seqsIdx)
        seqs.push_back(seq.first);
    size_t seqLen = util->seqsLen[T->root->identifier];
    std::cout << "Final Alignment Length: " << seqLen << '\n';
    std::sort(seqs.begin(), seqs.end(), cmp);
    for (int s = 0; s < seqs.size(); ++s) {
        int sIdx = util->seqsIdx[seqs[s]];
        int storage = util->seqsStorage[sIdx];
        if (std::find(T->root->msaIdx.begin(), T->root->msaIdx.end(), sIdx) != T->root->msaIdx.end()) {
            outFile << '>' << seqs[s] << "\n";
            if (option->outType == "FASTA") {
                outFile.write(&util->alnStorage[storage][sIdx][0], seqLen);
            }
            else if (option->outType == "CIGAR") {
                std::string cigar = "";
                char type = (util->alnStorage[storage][sIdx][0] == '-') ? 'D' : 'M';
                int num = 0;
                int i = 0;
                while (util->alnStorage[storage][sIdx][i] != 0) {
                    if (util->alnStorage[storage][sIdx][i] == '-') {
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
                outFile << cigar;
            }
            outFile << '\n';
        }
        util->seqFree(sIdx);
    }
    util->seqsFree();
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
    int seqLen = util->seqsLen[T->root->identifier];
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

void outputSubtree(Tree* tree, msa::option* option, int subtreeIdx)
{
    std::string subtreeFileName = "subtree-" + std::to_string(subtreeIdx);
    std::string subtreeTreeFile = option->tempDir + '/' + subtreeFileName + ".nwk";
    std::string out_str = "";
	getSubtreeNewick(tree->root, out_str);
	out_str += ";\n";
	std::ofstream outFile(subtreeTreeFile);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", subtreeTreeFile.c_str());
        exit(1);
    }
	outFile << out_str;
	outFile.close();
    return;
}

void outputSubtreeSeqs(std::string fileName, std::map<std::string, std::string>& seqs) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    int alnLen = seqs.begin()->second.size();
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        // assert(it->second.size() == alnLen);
        outFile << ('>' + it->first + '\n');
        outFile << (it->second + '\n');
    }
    outFile.close();
}

void outputSubtreeCIGAR(std::string fileName, std::map<std::string, std::string>& seqs) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        outFile << ('>' + it->first + '\n');
        std::string cigar = "";
        char type = (it->second[0] == '-') ? 'D' : 'M';
        int num = 0;
        int i = 0;
        while (i < it->second.size()) {
            if (it->second[i] == '-') {
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
        cigar += (std::to_string(num) + type + '\n');
        outFile << cigar;
    }

    outFile.close();
}

void storeFreq(msa::utility* util, Tree* T, int grpID) {
    int seqLen = util->seqsLen[T->root->identifier];
    std::cout << "Alignment Length: " << seqLen << '\n';
    util->profileFreq[grpID] = std::vector<std::vector<float>> (seqLen, std::vector<float> (6, 0.0));
    float totalWeight = 0;
    for (auto sIdx: T->root->msaIdx) totalWeight += T->allNodes[util->seqsName[sIdx]]->weight;
    
    for (int sIdx: T->root->msaIdx) {
        int storage = util->seqsStorage[sIdx];
        std::string name = util->seqsName[sIdx];
        float w = T->allNodes[name]->weight / totalWeight * T->root->numLeaves;
        tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
        for (int j = r.begin(); j < r.end(); ++j) {
        // for (int j = 0; j <  seqLen; ++j) {
            if      (util->alnStorage[storage][sIdx][j] == 'A' || util->alnStorage[storage][sIdx][j] == 'a') util->profileFreq[grpID][j][0]+=1.0*w;
            else if (util->alnStorage[storage][sIdx][j] == 'C' || util->alnStorage[storage][sIdx][j] == 'c') util->profileFreq[grpID][j][1]+=1.0*w;
            else if (util->alnStorage[storage][sIdx][j] == 'G' || util->alnStorage[storage][sIdx][j] == 'g') util->profileFreq[grpID][j][2]+=1.0*w;
            else if (util->alnStorage[storage][sIdx][j] == 'T' || util->alnStorage[storage][sIdx][j] == 't' ||
                     util->alnStorage[storage][sIdx][j] == 'U' || util->alnStorage[storage][sIdx][j] == 'u') util->profileFreq[grpID][j][3]+=1.0*w;
            else if (util->alnStorage[storage][sIdx][j] == 'N' || util->alnStorage[storage][sIdx][j] == 'n') util->profileFreq[grpID][j][4]+=1.0*w;
            else                                                                                             util->profileFreq[grpID][j][5]+=1.0*w;
        }
        });
    }
    return;
}

// auxiliary
bool cmp(std::string a, std::string b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
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

double calSPScore(std::string alnFile, msa::utility* util, Params* param) {
    
    gzFile f_rd = gzopen(alnFile.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", alnFile.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    std::vector<std::vector<int>> freq;
    int seqNum = 0, alnLen;
    uint64_t totalLen = 0;
    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        if (seqNum == 0) {
            alnLen = seqLen;
            freq = std::vector<std::vector<int>> (alnLen, std::vector<int>(6,0));
        }
        else {
            if (seqLen != alnLen) {
                std::cerr << "ERROR: The alignment lengths are not equal.\n";
                exit(1);
            }
        }
        std::string seqName = kseq_rd->name.s;
        std::string seq = std::string(kseq_rd->seq.s, seqLen);
        tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            if      (seq[s] == 'A' || seq[s] == 'a') freq[s][0]+=1;
            else if (seq[s] == 'C' || seq[s] == 'c') freq[s][1]+=1;
            else if (seq[s] == 'G' || seq[s] == 'g') freq[s][2]+=1;
            else if (seq[s] == 'T' || seq[s] == 't' ||
                     seq[s] == 'U' || seq[s] == 'u') freq[s][3]+=1;
            else if (seq[s] == 'N' || seq[s] == 'n') freq[s][4]+=1;
            else if (seq[s] == '-')                  freq[s][5]+=1;
        }
        });
        ++seqNum;
    }
    
    std::vector<double> spscore (alnLen, 0);
    tbb::parallel_for(tbb::blocked_range<int>(0, alnLen), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        for (int l = 0; l < 5; l++) spscore[s] += param->scoringMatrix[l][l] * (freq[s][l] * (freq[s][l] - 1) / 2);
        for (int l = 0; l < 4; l++) spscore[s] += param->gapExtend * (freq[s][l] * (freq[s][5]));
        for (int l = 0; l < 4; l++) {
            for (int m = l + 1; m < 4; m++) {
                spscore[s] += param->scoringMatrix[l][m] * (freq[s][m] * freq[s][l]);
            }
        }  
    }
    });
    double normScore = 0;
    for (int i = 0; i < alnLen; i++) normScore += spscore[i];
    normScore /= ((seqNum * (seqNum-1)) / 2);
    return normScore;
}   
