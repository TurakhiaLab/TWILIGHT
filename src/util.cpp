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

// read
void readSequences(msa::utility* util, msa::option* option, Tree* tree)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = option->seqFile;
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
        int kk = 0;
        for (auto node: tree->allNodes) {
            if (node.second->is_leaf()) {
                if (seqs.find(node.second->identifier) == seqs.end()) {
                    std::cout << "Missing " << node.second->identifier << '\n';
                    ++kk;
                }
            }
        }
        exit(1);
    }
    util->seqsMallocNStore(maxLen, seqs, option);
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

void readFrequency(msa::utility* util, msa::option* option) {
    auto freqReadStart = std::chrono::high_resolution_clock::now();
    std::string path = option->msaDir;
    std::vector<std::string> seqs (0);
    std::vector<std::string> subroots;
    std::vector<std::string> files;
    std::vector<int> msaLens;
    int subtreeIdx = 0;
    int totalFile = 0;
    for (const auto & msaFile : fs::directory_iterator(path)) {
        totalFile += 1;
        files.push_back(msaFile.path());
    }
    // for (const auto & msaFile : fs::directory_iterator(path)) {
    for (auto msaFileName: files) {
        // std::string msaFileName = msaFile.path();
        gzFile f_rd = gzopen(msaFileName.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
            exit(1);
        }
        std::cout << "Start reading " << msaFileName << " (" << subtreeIdx+1 << '/' << totalFile << ")\n";
        kseq_t* kseq_rd = kseq_init(f_rd);
        int seqNum = 0, msaLen;
        std::string seqName;
        while (kseq_read(kseq_rd) >= 0) {
            int seqLen = kseq_rd->seq.l;
            if (seqNum == 0) {
                msaLen = seqLen;
                seqName = kseq_rd->name.s;
                subroots.push_back(seqName);
                msaLens.push_back(msaLen);
            }
            else {
                if (seqLen != msaLen) {
                    fprintf(stderr, "ERROR: seqeunce length does not match in %s\n", msaFileName.c_str());
                    exit(1);
                }
            }
            std::string seq = std::string(kseq_rd->seq.s, seqLen);
            seqs.push_back(seq);
            ++seqNum;
        }
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        util->profileFreq[subtreeIdx] = std::vector<std::vector<float>> (msaLen, std::vector<float> (6, 0.0));
        for (auto seq: seqs) {
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, msaLen), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
                if      (seq[j] == 'A' || seq[j] == 'a') util->profileFreq[subtreeIdx][j][0]+=1.0;
                else if (seq[j] == 'C' || seq[j] == 'c') util->profileFreq[subtreeIdx][j][1]+=1.0;
                else if (seq[j] == 'G' || seq[j] == 'g') util->profileFreq[subtreeIdx][j][2]+=1.0;
                else if (seq[j] == 'T' || seq[j] == 't' ||
                         seq[j] == 'U' || seq[j] == 'u') util->profileFreq[subtreeIdx][j][3]+=1.0;
                else if (seq[j] == '-')                  util->profileFreq[subtreeIdx][j][5]+=1.0;
                else                                     util->profileFreq[subtreeIdx][j][4]+=1.0;
            }
            });
            });
        }
        if (option->printDetail) {
            std::cout << "File " << subtreeIdx << '(' << msaFileName << ") Num: " << seqNum << ", Length: " << msaLen << '\n';
        }
        util->seqsLen[seqName] = msaLen;
        util->seqsIdx[seqName] = subtreeIdx;
        util->seqsName[subtreeIdx] = seqName;
        if (msaLen > util->seqLen) util->seqLen = msaLen;
        seqs.clear();
        ++subtreeIdx;
    }
    // for (int i = 0; i < subroots.size(); ++i) {
    //     tree->allNodes[subroots[i]]->msaAln = std::vector<int8_t> (msaLens[i], 0);
    // }
    auto freqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = freqReadEnd - freqReadStart;
    std::cout << "MSA files read in " <<  seqReadTime.count() / 1000000 << " ms\n";
    return;
}

Tree* readNewick(std::string treeFileName)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; 
    // inputStream >> std::noskipws >> newick;
    std::getline(inputStream, newick);
    Tree *T = new Tree(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Newick string read in: " <<  treeBuiltTime.count() / 1000000 << " ms\n";

    return T;
}

// Output 
void outputSubtreeTrees(Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option)
{
    std::string tempDir = option->tempDir;
    for (auto subroot: partition->partitionsRoot) {
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        Tree* subT = new Tree(subroot.second.first);
        std::string subtreeName = "subtree-" + std::to_string(subtreeIdx);
        std::string subtreeTreeFile = option->tempDir + '/' + subtreeName + ".nwk";
        std::string out_str = "";
	    getSubtreeNewick(subT->root, out_str);
	    out_str += ";\n";
        std::ofstream outFile(subtreeTreeFile);
        if (!outFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", subtreeTreeFile.c_str());
            exit(1);
        }
	    outFile << out_str;
	    outFile.close();
    }
    return;
}

void outputFinal (Tree* tree, partitionInfo_t* partition, msa::utility* util, msa::option* option, int& totalSeqs) {
    std::vector<std::pair<std::string, std::string>> seqs;    

    if (option->alnMode == 0) {
        util->seqsIdx.clear();
        std::string seqFileName = option->seqFile;
        std::map<int, std::pair<std::string, std::string>> rawSeqs;
        tbb::spin_rw_mutex  writeMutex;
        std::cout << "Final Alignment Length: " << tree->allNodes[tree->root->identifier]->msaAln.size() << '\n';
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
                    seqs.push_back(std::make_pair(seqName, alnSeq));
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
    }
    else if (option->alnMode == 1) {
        std::string path = option->msaDir;
        std::cout << "Final Alignment Length: " << tree->allNodes[tree->root->children[0]->identifier]->msaAln.size() << '\n';
        int subtreeIdx = 0;
        int totalFile = 0;
        for (const auto & msaFile : fs::directory_iterator(path)) totalFile += 1;
        for (const auto & msaFile : fs::directory_iterator(path)) {
            std::string msaFileName = msaFile.path();
            gzFile f_rd = gzopen(msaFileName.c_str(), "r");
            if (!f_rd) {
                fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
                exit(1);
            }
            msaFileName = msaFileName.substr(option->msaDir.size()+1, msaFileName.size());
            std::string subtreeSeqFile = option->tempDir + '/' + msaFileName + ".final.aln";
            std::cout << "Start writing " << subtreeSeqFile << " (" << subtreeIdx+1 << '/' << totalFile << ")\n";
            kseq_t* kseq_rd = kseq_init(f_rd);
            int seqNum = 0, msaLen;
            std::string seqName;
            while (kseq_read(kseq_rd) >= 0) {
                int seqLen = kseq_rd->seq.l;
                if (seqNum == 0) {
                    msaLen = seqLen;
                    seqName = kseq_rd->name.s;
                }
                else {
                    if (seqLen != msaLen) {
                        fprintf(stderr, "ERROR: seqeunce length does not match in %s\n", msaFileName.c_str());
                        exit(1);
                    }
                }
                std::string seq = std::string(kseq_rd->seq.s, seqLen);
                seqs.push_back(std::make_pair(kseq_rd->name.s, seq));
                ++seqNum;
            }
            kseq_destroy(kseq_rd);
            gzclose(f_rd);
            tbb::parallel_for(tbb::blocked_range<int>(0, seqs.size()), [&](tbb::blocked_range<int> range){ 
            for (int n = range.begin(); n < range.end(); ++n) {
                std::string seq = seqs[n].second, alnSeq = "";
                int r = 0;
                for (size_t j = 0; j < tree->allNodes[seqName]->msaAln.size(); ++j) {
                    if ((tree->allNodes[seqName]->msaAln[j] & 0xFFFF) == 0 || (tree->allNodes[seqName]->msaAln[j] & 0xFFFF) == 2) {
                        alnSeq += seq[r];
                        ++r;
                    }
                    else {
                        alnSeq += '-';
                    }
                }
                if (alnSeq.size() != tree->allNodes[seqName]->msaAln.size()) {
                    std::cout << "ERROR: length not match. alnSeq (" << alnSeq.size() << ") != msaAln (" << tree->allNodes[seqName]->msaAln.size() << '\n'; 
                }
                assert(alnSeq.size() == tree->allNodes[seqName]->msaAln.size());
                seqs[n].second = alnSeq;
            }
            });
            
            if (option->outType == "FASTA") outputSubtreeSeqs(subtreeSeqFile, seqs);   
            else if (option->outType == "CIGAR") outputSubtreeCIGAR(subtreeSeqFile, seqs);  
            totalSeqs += seqs.size(); 
            seqs.clear();
            ++subtreeIdx;
        }
    }
    return;
}

void outputAln(msa::utility* util, msa::option* option, Tree* T) {
    std::string fileName = option->outFile;
    if (util->nowProcess == 2) {
        std::string command = "cat " + option->tempDir + "/*.final.aln > " + fileName;
        int catResult = system(command.c_str());
        if (catResult != 0) {
            if (option->deleteTemp) std::cerr << "ERROR: Unable to concatenate subtree alignments. Please check the available space on your device and merge files with \"cat <temp-dir>/*.final.aln <output-file>\". Temporary files will not be deleted.\n";
            else                    std::cerr << "ERROR: Unable to concatenate subtree alignments. Please check the available space on your device and merge files with \"cat <temp-dir>/*.final.aln <output-file>\".\n";
        }
        if (catResult == 0 && option->deleteTemp) {
            command = "rm -rf " + option->tempDir;
            int delResult = system(command.c_str());
            if (delResult != 0) std::cerr << "ERROR: Unable to delete temporary files.\n";
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

void outputSubtreeSeqs(std::string fileName, std::vector<std::pair<std::string, std::string>>& seqs) {
    std::ofstream outFile(fileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    // int alnLen = seqs.begin()->second.size();
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        // assert(it->second.size() == alnLen);
        outFile << ('>' + it->first + '\n');
        outFile << (it->second + '\n');
    }
    outFile.close();
}

void outputSubtreeCIGAR(std::string fileName, std::vector<std::pair<std::string, std::string>>& seqs) {
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
        auto name = root->identifier;
        bool specialChar = (name.find(',') != std::string::npos)  || \
                           (name.find(':') != std::string::npos)  || \
                           (name.find('(') != std::string::npos)  || \
                           (name.find(')') != std::string::npos);
        if (specialChar) name = '\'' + name + '\'';
        outputString += (name + ':' + std::to_string(root->branchLength));
    }
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
            else if (util->alnStorage[storage][sIdx][j] == '-')                                              util->profileFreq[grpID][j][5]+=1.0*w;
            else                                                                                             util->profileFreq[grpID][j][4]+=1.0*w;
        }
        });
    }
    return;
}

double calSPScore(std::string alnFile, msa::utility* util, Params* param) {
    
    gzFile f_rd = gzopen(alnFile.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", alnFile.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    std::vector<std::vector<int>> freq;
    int seqNum = 0, alnLen = 0;
    // uint64_t totalLen = 0;
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
        for (int l = 0; l < 5; l++) spscore[s] += param->gapExtend * (freq[s][l] * (freq[s][5]));
        for (int l = 0; l < 5; l++) {
            for (int m = 0; m < 5; m++) {
                if (l != m) spscore[s] += param->scoringMatrix[l][m] * (freq[s][m] * freq[s][l]);
            }
        }  
    }
    });
    double normScore = 0;
    for (int i = 0; i < alnLen; i++) normScore += spscore[i];
    float matchAvg = 0;
    for (int i = 0; i < 4; ++i) matchAvg += param->scoringMatrix[i][i];
    matchAvg /= 4.0;
    normScore /= ((seqNum * (seqNum-1)) / 2);
    // normScore /= alnLen;
    normScore /= matchAvg;
    return normScore;
}   

