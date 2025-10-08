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
void readSequences(msa::utility* util, msa::option* option, Tree*& tree)
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

    std::vector<std::pair<std::string, std::string>> seqsOut;
    
    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        std::string seqName_full = kseq_rd->name.s;
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";
        if (tree->allNodes.find(seqName_full) != tree->allNodes.end()) seqName = seqName_full;
        else if (tree->allNodes.find(seqName_noblank) != tree->allNodes.end()) seqName = seqName_noblank;
        if (seqName != "") {
            if (seqs.find(seqName) != seqs.end()) {
                printf("WARNING: duplicate leaf names found in the sequence file! Leaf name: %s. Only the first occurrence will be kept.\n", seqName.c_str());
                // exit(1);
            }
            else {
                int subtreeIdx = tree->allNodes[seqName]->grpID;
                std::string seq = std::string(kseq_rd->seq.s, seqLen);
                seqs[seqName] = std::make_pair(seq, subtreeIdx);
                if (option->debug) util->rawSeqs[seqName] = seq;
                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
                if (seqLen == 0) std::cout << "Null sequences, " << seqName << '\n';
                totalLen += seqLen;
            }
        }
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);


    seqNum = seqs.size();
    uint32_t avgLen = totalLen/seqNum;
    util->maxRawSeqLen = maxLen;
    if (tree->m_numLeaves != seqNum && option->alnMode != 2) {
        printf("Warning: Mismatch between the number of leaves and the number of sequences, (%lu != %d)\n", tree->m_numLeaves, seqNum); 
        int kk = 0;
        for (auto node: tree->allNodes) {
            if (node.second->is_leaf()) {
                if (seqs.find(node.second->identifier) == seqs.end()) {
                    std::cout << "Missing " << node.second->identifier << '\n';
                    ++kk;
                }
            }
        }
        std::cout << "Prune the tree according to the existing sequences.\n";
        tree = pruneTree(tree, seqs);
    }
    util->seqsMallocNStore(maxLen, seqs, option);
    for (int i = 0; i < seqNum; ++i) {
        util->lowQuality[i] = false;
    }
    float minLenTh = avgLen * (1-option->lenDev), maxLenTh = avgLen * (1+option->lenDev);
    std::atomic<int> numLowQ;
    numLowQ.store(0);
    int ambig = (option->type == 'n') ? 4 : 20;
    tbb::parallel_for(tbb::blocked_range<int>(0, seqNum), [&](tbb::blocked_range<int> range){ 
    for (int i = range.begin(); i < range.end(); ++i) {
    // for (int i = 0; i < seqNum; ++i) {
        std::string name = util->seqsName[i];
        int len = util->seqsLen[name];
        if (option->lenDev > 0) util->lowQuality[i] = (len > maxLenTh || len < minLenTh);
        if (!util->lowQuality[i]) {
            int countN = 0;
            for (int j = 0; j < len; ++j) {
                if (letterIdx(option->type, toupper(util->alnStorage[0][i][j])) == ambig) ++countN;
            }
            util->lowQuality[i] = (countN > (len * option->maxAmbig));
        }
        if (util->lowQuality[i]) numLowQ.fetch_add(1);
    }
    });

    bool outputLowQ = false;
    if (outputLowQ){
        std::string lowQfile;
        if (option->tempDir != "") lowQfile = option->tempDir+"/lowQ-"+std::to_string(numLowQ)+"-"+std::to_string(seqNum)+".fasta.gz";
        else lowQfile = option->treeFile+".lowQ-"+std::to_string(numLowQ)+"-"+std::to_string(seqNum)+".fasta.gz";
        gzFile outLow = gzopen(lowQfile.c_str(), "wb"); // "wb" = write binary
        if (!outLow) {
            std::cerr << "Failed to open " << outLow << std::endl;
            exit(1);
        }
        for (int sIdx = 0; sIdx < seqNum; ++sIdx) {
            if (util->lowQuality[sIdx]) {
                int storage = util->seqsStorage[sIdx];
                std::string name = util->seqsName[sIdx];
                int sLen = util->seqsLen[name];
                gzprintf(outLow, ">%s\n", name.c_str());
                gzwrite(outLow, &util->alnStorage[storage][sIdx][0], sLen);
                gzprintf(outLow, "\n");
            }
        }
        gzclose(outLow);
    }
    


    if (option->printDetail) {
        std::cout << "===== Sequence Summary =====\n";
        std::cout << "Number : " << seqNum << '\n';
        std::cout << "Max. Length: " << maxLen << '\n';
        std::cout << "Min. Length: " << minLen << '\n';
        std::cout << "Avg. Length: " << avgLen << '\n';
        if (option->noFilter) 
        std::cout << "Deferred sequences: " << numLowQ << '\n';
        else 
        std::cout << "Excluded sequences: " << numLowQ << '\n';
        std::cout << "=============================\n";
    }
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
}

void readFrequency(msa::utility* util, msa::option* option) {
    auto freqReadStart = std::chrono::high_resolution_clock::now();
    std::string path = option->msaDir;
    std::vector<std::string> files;
    int subtreeIdx = 0;
    int totalFile = 0;
    boost::system::error_code ec;
    fs::path msaFolder(path);
    for (fs::recursive_directory_iterator it(msaFolder, ec), eit; it != eit; it.increment(ec)) {
        if (ec) {
            it.pop();
            continue;
        }
        if (!fs::is_directory(it->path())) {
            files.push_back(it->path().string());
            totalFile += 1;
        }
    }
    // for (const auto & msaFile : fs::directory_iterator(path)) {
    std::sort(files.begin(), files.end());
    int profileSize = option->type == 'n' ? 6 : 22;
    for (auto msaFileName: files) {
        // std::string msaFileName = msaFile.path();
        gzFile f_rd = gzopen(msaFileName.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
            exit(1);
        }
        std::cout << "Start reading " << msaFileName << " (" << subtreeIdx+1 << '/' << totalFile << ")\n";
        kseq_t* kseq_rd = kseq_init(f_rd);
        int seqNum = 0, msaLen = 0;
        std::string seqName;
        while (kseq_read(kseq_rd) >= 0) {
            int seqLen = kseq_rd->seq.l;
            if (seqNum == 0) {
                msaLen = seqLen;
                seqName = kseq_rd->name.s;
                util->profileFreq[subtreeIdx] = std::vector<std::vector<float>> (msaLen, std::vector<float> (profileSize, 0.0));
        
            }
            else {
                if (seqLen != msaLen) {
                    fprintf(stderr, "ERROR: seqeunce length does not match in %s\n", msaFileName.c_str());
                    exit(1);
                }
            }
            std::string seq = std::string(kseq_rd->seq.s, seqLen);
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
                int letterIndex = letterIdx(option->type, toupper(seq[j]));
                util->profileFreq[subtreeIdx][j][letterIndex] += 1.0;
            }
            });
            });
            // seqs.push_back(seq);
            ++seqNum;
        }
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
        if (option->printDetail) {
            std::cout << "File " << subtreeIdx << '(' << msaFileName << ") Num: " << seqNum << ", Length: " << msaLen << '\n';
        }
        util->seqsLen[seqName] = msaLen;
        util->seqsIdx[seqName] = subtreeIdx;
        util->seqsName[subtreeIdx] = seqName;
        if (msaLen > util->seqLen) util->seqLen = msaLen;
        ++subtreeIdx;
    }
    auto freqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = freqReadEnd - freqReadStart;
    std::cout << "MSA files read in " <<  seqReadTime.count() / 1000000 << " ms\n";
    return;
}

void readMSA_and_Seqs(msa::utility* util, msa::option* option, Tree* tree) {
    auto freqReadStart = std::chrono::high_resolution_clock::now();
    int profileSize = (option->type == 'n') ? 6 : 22;
    // Construct Tree
    {
        Node* treeRoot = new Node("node_1", 0.0);
        tree->root = treeRoot;
        tree->root->grpID = -1;
        tree->allNodes["node_1"] = treeRoot;
    }
    // Read sequence file
    {
        std::string seqFileName = option->seqFile;
        gzFile seq_rd = gzopen(seqFileName.c_str(), "r");
        if (!seq_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
            exit(1);
        }
        kseq_t* kseq_rd_2 = kseq_init(seq_rd);
        while (kseq_read(kseq_rd_2) >= 0) {
            int seqLen = kseq_rd_2->seq.l;
            util->placedSeqs[kseq_rd_2->name.s] = std::string(kseq_rd_2->seq.s, seqLen);
        }
        kseq_destroy(kseq_rd_2);
        gzclose(seq_rd);
        for (auto seq: util->placedSeqs) {
            Node* newNode = new Node (seq.first, tree->root, 1.0);
            newNode->grpID = tree->root->grpID;
            tree->allNodes[seq.first] = newNode;
        }
    }
    // Read Backbone MSA
    {
        std::string msaFileName = option->backboneAlnFile;
        gzFile f_rd = gzopen(msaFileName.c_str(), "r");
        if (!f_rd) {
            fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
            exit(1);
        }
        kseq_t* kseq_rd = kseq_init(f_rd);
        int seqNum = 0, msaLen = 0;
        while (kseq_read(kseq_rd) >= 0) {
            std::string seqName = kseq_rd->name.s;
            int seqLen = kseq_rd->seq.l;
            if (seqNum == 0) {
                msaLen = seqLen;
                tree->root->msaFreq = std::vector<std::vector<float>> (seqLen, std::vector<float> (profileSize, 0.0));
            }
            else {
                if (seqLen != msaLen) {
                    fprintf(stderr, "ERROR: sequence lengths do not match in %s\n", msaFileName.c_str());
                    exit(1);
                }
            }
            std::string seq = std::string(kseq_rd->seq.s, seqLen);
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, seq.size()), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
                int letterIndex = letterIdx(option->type, toupper(seq[j]));
                tree->root->msaFreq[j][letterIndex] += 1.0;
            }
            });
            });
            tree->root->msa.push_back(seqName);
            ++seqNum;
        }
        kseq_destroy(kseq_rd);
        gzclose(f_rd);
    }

    auto freqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = freqReadEnd - freqReadStart;
    std::cout << "Backbone MSA and unaligned sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
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

void readBackboneAln(msa::utility* util, msa::option* option, Tree* tree)
{
    auto alnReadStart = std::chrono::high_resolution_clock::now();
    std::string fileName = option->backboneAlnFile;
    gzFile f_rd = gzopen(fileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);

    int seqNum = 0, alnLen = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        std::string seqName_full = kseq_rd->name.s;
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";
        if (tree->allNodes.find(seqName_full) != tree->allNodes.end()) seqName = seqName_full;
        else if (tree->allNodes.find(seqName_noblank) != tree->allNodes.end()) seqName = seqName_noblank;
        if (seqName != "") {
            if (alnLen == 0) alnLen = kseq_rd->seq.l;
            else if (alnLen != kseq_rd->seq.l) {
                std::cerr << "ERROR: Backbone MSA contains sequences of unequal lengths." << std::endl;
                exit(1);
            }
            if (util->backboneAln.find(seqName) != util->backboneAln.end()) {
                std::cerr << "ERROR: Sequence " << seqName << " already exists.\n";
                exit(1);
            }
            util->backboneAln[seqName] = std::string(kseq_rd->seq.s, kseq_rd->seq.l);
        }
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    auto alnReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnReadTime = alnReadEnd - alnReadStart;
    std::cout << "Backbone MSA read in " <<  alnReadTime.count() / 1000000 << " ms\n";
    return;
}

void readNewSequences(msa::utility* util, msa::option* option, Tree* tree)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string fileName = option->seqFile;
    gzFile f_rd = gzopen(fileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);

    int seqNum = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        std::string seqName_full = kseq_rd->name.s;
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";
        if (tree->allNodes.find(seqName_full) != tree->allNodes.end()) seqName = seqName_full;
        else if (tree->allNodes.find(seqName_noblank) != tree->allNodes.end()) seqName = seqName_noblank;
        if (seqName != "") {
            if (util->placedSeqs.find(seqName) != util->placedSeqs.end()) {
                std::cerr << "ERROR: Sequence " << seqName << " already exists.\n";
                exit(1);
            }
            util->placedSeqs[seqName] = std::string(kseq_rd->seq.s, kseq_rd->seq.l);
        }
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "New sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";
    return;
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
        // std::string seqFileName = option->seqFile;
        std::map<int, std::pair<std::string, std::string>> rawSeqs;
        tbb::spin_rw_mutex  writeMutex;
        std::cout << "Final alignment length: " << tree->allNodes[tree->root->identifier]->msaAln.size() << '\n';
        int proceeded = 0;    
        for (auto subroot: partition->partitionsRoot) {
            int subtreeIdx = tree->allNodes[subroot.first]->grpID;
            std::string subalnFileName = option->tempDir + "/subtree-" + std::to_string(subtreeIdx) + ".subalignment.aln";
            if (option->compressed) subalnFileName += ".gz";
            gzFile f_rd = gzopen(subalnFileName.c_str(), "r");
            if (!f_rd) { fprintf(stderr, "ERROR: fail to open file: %s\n", subalnFileName.c_str()); exit(1);}
            kseq_t* kseq_rd = kseq_init(f_rd);
            ++proceeded;
            std::cout << "Start writing alignment of subtree No. " << subtreeIdx << ". (" << proceeded << '/' << partition->partitionsRoot.size() << ")\n";
            Tree* subT = new Tree(subroot.second.first);
            int s = 0;
            while (kseq_read(kseq_rd) >= 0) {
                std::string seqName = kseq_rd->name.s;
                size_t seqLen = kseq_rd->seq.l;
                rawSeqs[s] = std::make_pair(seqName, std::string(kseq_rd->seq.s, seqLen));
                ++s;
            }
            kseq_destroy(kseq_rd);
            gzclose(f_rd);
            // assert(subT->m_numLeaves == rawSeqs.size());
            if (subT->m_numLeaves != rawSeqs.size()) std::cout << "Mismatch between tree and seqs: " << subT->m_numLeaves << '/' << rawSeqs.size() << '\n';
            seqs = std::vector<std::pair<std::string, std::string>>(rawSeqs.size(), std::make_pair("",""));
            tbb::parallel_for(tbb::blocked_range<int>(0, rawSeqs.size()), [&](tbb::blocked_range<int> range){ 
            for (int n = range.begin(); n < range.end(); ++n) {
                std::string seqName = rawSeqs[n].first;
                std::string updatedSeq = rawSeqs[n].second, alnSeq = "";
                size_t r = 0;         
                for (size_t j = 0; j < tree->allNodes[subroot.first]->msaAln.size(); ++j) {
                    if ((tree->allNodes[subroot.first]->msaAln[j] & 0xFFFF) == 0 || (tree->allNodes[subroot.first]->msaAln[j] & 0xFFFF) == 2) {
                        alnSeq += updatedSeq[r];
                        ++r;
                    }
                    else {
                        alnSeq += '-';
                    }
                }
                if (alnSeq.size() != tree->allNodes[subroot.first]->msaAln.size()) {
                    std::cout << "ERROR: " << seqName << " length not match. alnSeq (" << alnSeq.size() << ") != msaAln (" << tree->allNodes[subroot.first]->msaAln.size() << '\n'; 
                }
                // assert(alnSeq.size() == tree->allNodes[subroot.first]->msaAln.size());
                seqs[n].first = seqName;
                seqs[n].second = alnSeq;
            }
            });
            tree->allNodes[subroot.first]->msaAln.clear();
            totalSeqs += seqs.size();
            std::string subtreeSeqFile = option->tempDir + '/' + "subtree-" + std::to_string(subtreeIdx) + ".final.aln";
            outputAlignment(subtreeSeqFile, seqs, option->compressed);   
            delete subT;
            seqs.clear();
            rawSeqs.clear();
            if (option->deleteTemp) {
                std::string command = "rm " + subalnFileName;
                int delResult = system(command.c_str());
                if (delResult != 0) std::cerr << "ERROR: Unable to delete subalignment file.\n";
            }
        }
        tree->root->msa.clear();
    }
    else if (option->alnMode == 1) {
        std::string path = option->msaDir;
        std::cout << "Final Alignment Length: " << tree->allNodes[tree->root->children[0]->identifier]->msaAln.size() << '\n';
        int subtreeIdx = 0;
        int totalFile = 0;
        std::vector<std::string> files;
        boost::system::error_code ec;
        fs::path msaFolder(path);
        for (fs::recursive_directory_iterator it(msaFolder, ec), eit; it != eit; it.increment(ec)) {
            if (ec) {
                it.pop();
                continue;
            }
            if (!fs::is_directory(it->path())) {
                totalFile += 1;
                files.push_back(it->path().string());
            }
        }
        // for (const auto & msaFile : fs::directory_iterator(path)) totalFile += 1;
        std::sort(files.begin(), files.end());
        for (auto msaFile : files) {
            std::string msaFileName = msaFile;
            gzFile f_rd = gzopen(msaFileName.c_str(), "r");
            if (!f_rd) {
                fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
                exit(1);
            }
            msaFileName = msaFileName.substr(option->msaDir.size(), msaFileName.size());
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
            
            outputAlignment(subtreeSeqFile, seqs, option->compressed);   
            totalSeqs += seqs.size(); 
            seqs.clear();
            ++subtreeIdx;
        }
    }
    else if (option->alnMode == 2) {
        std::cout << "Final Alignment Length: " << tree->allNodes[tree->root->children[0]->identifier]->msaAln.size() << '\n';
        std::string refName = "";
        {
            // Output updated backbone alignment
            gzFile f_rd = gzopen(option->backboneAlnFile.c_str(), "r");
            if (!f_rd) {
                fprintf(stderr, "ERROR: cant open file: %s\n", option->backboneAlnFile.c_str());
                exit(1);
            }
            fs::path p_backbone(option->backboneAlnFile);
            std::string subtreeSeqFile = option->tempDir + '/' + p_backbone.filename().string().c_str() + ".final.aln";
            if (option->printDetail) std::cout << "Start writing updated backbone alignment\n";
            kseq_t* kseq_rd = kseq_init(f_rd);
            int seqNum = 0;
            while (kseq_read(kseq_rd) >= 0) {
                int seqLen = kseq_rd->seq.l;
                if (seqNum == 0) refName = kseq_rd->name.s;
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
                for (size_t j = 0; j < tree->allNodes[refName]->msaAln.size(); ++j) {
                    if ((tree->allNodes[refName]->msaAln[j] & 0xFFFF) == 0 || (tree->allNodes[refName]->msaAln[j] & 0xFFFF) == 2) {
                        alnSeq += seq[r];
                        ++r;
                    }
                    else {
                        alnSeq += '-';
                    }
                }
                if (alnSeq.size() != tree->allNodes[refName]->msaAln.size()) {
                    std::cout << "ERROR: length not match. alnSeq (" << alnSeq.size() << ") != msaAln (" << tree->allNodes[refName]->msaAln.size() << '\n'; 
                }
                assert(alnSeq.size() == tree->allNodes[refName]->msaAln.size());
                seqs[n].second = alnSeq;
            }
            });
            outputAlignment(subtreeSeqFile, seqs, option->compressed);   
            seqs.clear();
        }
        // Output alignment of newly added sequences
        {
            // Output updated backbone alignment
            gzFile f_rd = gzopen(option->seqFile.c_str(), "r");
            if (!f_rd) {
                fprintf(stderr, "ERROR: cant open file: %s\n", option->backboneAlnFile.c_str());
                exit(1);
            }
            fs::path p_new(option->seqFile);
            std::string subtreeSeqFile = option->tempDir + '/' + p_new.filename().string() + ".final.aln";
            if (option->printDetail) std::cout << "Start writing alignment of newly added sequences\n";
            kseq_t* kseq_rd = kseq_init(f_rd);
            while (kseq_read(kseq_rd) >= 0) {
                int seqLen = kseq_rd->seq.l;
                std::string seq = std::string(kseq_rd->seq.s, seqLen);
                seqs.push_back(std::make_pair(kseq_rd->name.s, seq));
            }
            kseq_destroy(kseq_rd);
            gzclose(f_rd);
            tbb::parallel_for(tbb::blocked_range<int>(0, seqs.size()), [&](tbb::blocked_range<int> range){ 
            for (int n = range.begin(); n < range.end(); ++n) {
                std::string seq = seqs[n].second, alnSeq = "";
                int r = 0;
                for (size_t j = 0; j < tree->allNodes[seqs[n].first]->msaAln.size(); ++j) {
                    if ((tree->allNodes[seqs[n].first]->msaAln[j] & 0xFFFF) == 0 || (tree->allNodes[seqs[n].first]->msaAln[j] & 0xFFFF) == 2) {
                        alnSeq += seq[r];
                        ++r;
                    }
                    else {
                        alnSeq += '-';
                    }
                }
                if (alnSeq.size() != tree->allNodes[refName]->msaAln.size()) {
                    std::cout << "ERROR: length not match. alnSeq (" << alnSeq.size() << ") != msaAln (" << tree->allNodes[refName]->msaAln.size() << '\n'; 
                }
                assert(alnSeq.size() == tree->allNodes[refName]->msaAln.size());
                seqs[n].second = alnSeq;
            }
            });
            outputAlignment(subtreeSeqFile, seqs, option->compressed);   
            seqs.clear();
        }
    }
    return;
}

void outputSubAln(msa::utility* util, msa::option* option, Tree* T, int subtreeIdx) {
    std::string fileName = option->tempDir + "/subtree-" + std::to_string(subtreeIdx) + ".subalignment.aln";
    if (option->alnMode != 2) {
        std::vector<std::pair<std::string, std::string>> seqs;
        size_t seqLen = util->seqsLen[T->root->identifier];
        std::cout << "Subalignment Length: " << seqLen << '\n';
        outputAlignment(fileName, util, T, option->compressed);
    }
    else {
        size_t seqLen = util->backboneAln.begin()->second.size();
        std::cout << "Subalignment Length: " << seqLen << '\n';
        std::vector<std::pair<std::string, std::string>> seqs;
        for (auto seq: util->backboneAln) {
            seqs.push_back(std::make_pair(seq.first, seq.second));
        }
        util->backboneAln.clear();
        for (auto seq: util->placedSeqs) {
            seqs.push_back(std::make_pair(seq.first, seq.second));
        }
        util->placedSeqs.clear();
        outputAlignment(fileName, seqs, option->compressed);
    }
    return;
}

void outputWholeAln(msa::utility* util, msa::option* option, Tree* T) {
    std::string fileName = option->outFile;
    if (util->nowProcess == 2) {
        if (option->compressed) fileName += ".gz";
        std::string command = "cat " + option->tempDir + "/*.final.aln* > " + fileName;
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
    if (option->alnMode != 2) {
        std::vector<std::pair<std::string,std::string>> seqs;
        size_t seqLen = util->seqsLen[T->root->identifier];
        std::cout << "Final Alignment Length: " << seqLen << '\n';
        outputAlignment(fileName, util, T, option->compressed);
    }
    else {
        size_t seqLen = (option->treeFile != "") ? util->backboneAln.begin()->second.size() : T->root->msaAln.size();
        std::cout << "Final Alignment Length: " << seqLen << '\n';
        int catResult = 0;
        // Write to backbone file
        fs::path p_backbone(option->backboneAlnFile);
        std::string backboneFile = option->tempDir + '/' + p_backbone.filename().string() + ".final.aln";
        {
            if (option->treeFile != "") {
                std::vector<std::pair<std::string, std::string>> seqs;
                for (auto seq: util->backboneAln) {
                    if (seq.second.size() != seqLen) std::cout << "ERROR: " << seq.first << " sequence length not match\n";
                    seqs.push_back(std::make_pair(seq.first, seq.second));
                }
                outputAlignment(backboneFile, seqs, option->compressed);
                seqs.clear();
            }
            else {
                std::ofstream refFile;
                if (option->compressed) {
                    refFile.open((backboneFile+".gz"), std::ios::binary);
                }
                else {
                    refFile.open(backboneFile);
                }
                if (!refFile) {
                    fprintf(stderr, "ERROR: fail to open file: %s\n", backboneFile.c_str());
                    exit(1);
                }
                std::vector<std::pair<std::string, std::string>> refSeq;
                std::string msaFileName = option->backboneAlnFile;
                int seqLen = T->root->msaAln.size();
                gzFile f_rd = gzopen(msaFileName.c_str(), "r");
                if (!f_rd) {
                    fprintf(stderr, "ERROR: cant open file: %s\n", msaFileName.c_str());
                    exit(1);
                }
                kseq_t* kseq_rd = kseq_init(f_rd);
                int batchSize = 8000000000 / seqLen; // 8GB memory limit
                while (kseq_read(kseq_rd) >= 0) {
                    std::string seqName = kseq_rd->name.s;
                    int seqLen = kseq_rd->seq.l;
                    refSeq.push_back({seqName, std::string(kseq_rd->seq.s, seqLen)});
                    if (refSeq.size() % batchSize == 0) {
                        tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<int>(0, refSeq.size()), [&](tbb::blocked_range<int> r) {
                        for (int j = r.begin(); j < r.end(); ++j) {
                            std::string updatedSeq = "";
                            int rIdx = 0;
                            for (auto a: T->root->msaAln) {
                                if (a == 0) {
                                    updatedSeq.push_back(refSeq[j].second[rIdx]);
                                    rIdx++;
                                }
                                else updatedSeq.push_back('-');
                            }
                            refSeq[j].second = updatedSeq;
                        }
                        });
                        });

                        std::vector<std::string> chunks(refSeq.size());
                        if (option->compressed) {
                            tbb::parallel_for(size_t(0), refSeq.size(), [&](size_t i) {
                                std::string fasta_chunk = ">" + refSeq[i].first + "\n" + refSeq[i].second + "\n";
                                chunks[i] = gzip_compress(fasta_chunk);
                            });
                        }
                        else {
                            tbb::parallel_for(size_t(0), refSeq.size(), [&](size_t i) {
                                std::string fasta_chunk = ">" + refSeq[i].first + "\n" + refSeq[i].second + "\n";
                                chunks[i] = fasta_chunk;
                            });
                        }
                        refSeq.clear();
                        for (auto &c : chunks) {
                            refFile.write(c.data(), c.size());
                        }

                    }
                }
                kseq_destroy(kseq_rd);
                gzclose(f_rd);
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(0, refSeq.size()), [&](tbb::blocked_range<int> r) {
                for (int j = r.begin(); j < r.end(); ++j) {
                    std::string updatedSeq = "";
                    int rIdx = 0;
                    for (auto a: T->root->msaAln) {
                        if (a == 0) {
                            updatedSeq.push_back(refSeq[j].second[rIdx]);
                            rIdx++;
                        }
                        else updatedSeq.push_back('-');
                    }
                    refSeq[j].second = updatedSeq;
                }
                });
                });

                std::vector<std::string> chunks(refSeq.size());
                if (option->compressed) {
                    tbb::parallel_for(size_t(0), refSeq.size(), [&](size_t i) {
                        std::string fasta_chunk = ">" + refSeq[i].first + "\n" + refSeq[i].second + "\n";
                        chunks[i] = gzip_compress(fasta_chunk);
                    });
                }
                else {
                    tbb::parallel_for(size_t(0), refSeq.size(), [&](size_t i) {
                        std::string fasta_chunk = ">" + refSeq[i].first + "\n" + refSeq[i].second + "\n";
                        chunks[i] = fasta_chunk;
                    });
                }
                refSeq.clear();
                for (auto &c : chunks) {
                    refFile.write(c.data(), c.size());
                }
                refFile.close();
            }
        }
        // Write to placed seq file
        fs::path p_new(option->seqFile);
        std::string newFile = option->tempDir + '/' + p_new.filename().string() + ".final.aln";
        {
            std::vector<std::pair<std::string, std::string>> seqs;
            for (auto seq: util->placedSeqs) {
                if (seq.second.size() != seqLen) std::cout << "ERROR: " << seq.first << " sequence length not match\n";
                seqs.push_back(std::make_pair(seq.first, seq.second));
            }
            outputAlignment(newFile, seqs, option->compressed);
            seqs.clear();
        }

            
        if (option->compressed) {
            fileName += ".gz";
            backboneFile += ".gz";
            newFile += ".gz";
        }
        std::string command = "cat " + backboneFile + " " + newFile + " > " + fileName;
        catResult = system(command.c_str());
        if (catResult != 0) {
            std::cerr << "ERROR: Unable to concatenate alignments.\n";
        }
        if (catResult == 0 && option->deleteTemp) {
            std::string command = "rm -rf " + option->tempDir;
            int delResult = system(command.c_str());
            if (delResult != 0) std::cerr << "ERROR: Unable to delete temporary files.\n";
        }
    }
    
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

void outputPrunedTree(Tree* T, msa::option* option)
{
    fs::path o(option->outFile);
    fs::path t(option->treeFile);
    std::string prunedTreeFileName = (o.parent_path().string() == "") ? t.filename().string() + ".pruned.nwk" : o.parent_path().string() + "/" + t.filename().string() + ".pruned.nwk";
    std::string out_str = "";
	getSubtreeNewick(T->root, out_str);
	out_str += ";\n";
	std::ofstream outFile(prunedTreeFileName);
    if (!outFile) {
        fprintf(stderr, "ERROR: cant open file: %s\n", prunedTreeFileName.c_str());
        exit(1);
    }
	outFile << out_str;
	outFile.close();
    return;
}

void outputSubtreeSeqs(std::string fileName, std::vector<std::pair<std::string, std::string>>& seqs, bool compressed) {
    std::sort(seqs.begin(), seqs.end(), cmp2);
    if (compressed) {
        fileName += ".gz";
        gzFile outFile = gzopen(fileName.c_str(), "wb"); // "wb" = write binary
        if (!outFile) {
            std::cerr << "Failed to open " << fileName << std::endl;
            exit(1);
        }
        for (auto it = seqs.begin(); it != seqs.end(); ++it) {
            std::string header = ">" + it->first + "\n";
            std::string body = it->second + "\n";
            gzwrite(outFile, header.c_str(), header.size());
            gzwrite(outFile, body.c_str(), body.size());
        }
        gzclose(outFile);
    }
    else {
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
}

void outputAlignment(std::string fileName, std::vector<std::pair<std::string, std::string>>& seqs, bool compressed) {
    std::sort(seqs.begin(), seqs.end(), cmp2);
    if (compressed) {
        fileName += ".gz";
        std::vector<std::string> compressed_chunks(seqs.size());
        tbb::parallel_for(size_t(0), seqs.size(), [&](size_t i) {
            std::string fasta_chunk = ">" + seqs[i].first + "\n" + seqs[i].second + "\n";
            compressed_chunks[i] = gzip_compress(fasta_chunk);
        });
        std::ofstream outFile(fileName, std::ios::binary);
        if (!outFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
            exit(1);
        }
        for (auto &c : compressed_chunks) {
            outFile.write(c.data(), c.size());
        }
        outFile.close();
    }
    else {
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
}

void outputAlignment(std::string fileName, msa::utility* util, Tree* T, bool compressed) {
    std::vector<std::string> seqs;
    for (auto seq: util->seqsIdx)
        seqs.push_back(seq.first);
    std::sort(seqs.begin(), seqs.end(), cmp);
    size_t seqLen = util->seqsLen[T->root->identifier];
    if (compressed) {
        fileName += ".gz";
        std::vector<std::string> compressed_chunks(seqs.size());
        tbb::parallel_for(size_t(0), seqs.size(), [&](size_t s) {
            int sIdx = util->seqsIdx[seqs[s]];
            int storage = util->seqsStorage[sIdx];
            if (std::find(T->root->msaIdx.begin(), T->root->msaIdx.end(), sIdx) != T->root->msaIdx.end()) {
                std::string fasta_chunk = ">" + seqs[s] + "\n" + std::string(&util->alnStorage[storage][sIdx][0], seqLen) + "\n";
                compressed_chunks[s] = gzip_compress(fasta_chunk);
            }
            util->seqFree(sIdx);
        });
        std::ofstream outFile(fileName, std::ios::binary);
        if (!outFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
            exit(1);
        }
        for (auto &c : compressed_chunks) {
            outFile.write(c.data(), c.size());
        }
        outFile.close();
    }
    else {
        std::ofstream outFile(fileName);
        if (!outFile) {
            fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
            exit(1);
        }
        for (int s = 0; s < seqs.size(); ++s) {
            int sIdx = util->seqsIdx[seqs[s]];
            int storage = util->seqsStorage[sIdx];
            if (std::find(T->root->msaIdx.begin(), T->root->msaIdx.end(), sIdx) != T->root->msaIdx.end()) {
                outFile << '>' << seqs[s] << "\n";
                outFile.write(&util->alnStorage[storage][sIdx][0], seqLen);
                outFile << '\n';
            }
            util->seqFree(sIdx);
        }
        outFile.close();
    }
    util->seqsFree();
}

// auxiliary
bool cmp(std::string a, std::string b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
}

bool cmp2(std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) {
    if (a.first.size() != b.first.size()) return a.first.size() < b.first.size();
    return a.first < b.first;
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

void storeFreq(msa::utility* util,  msa::option* option, Tree* T, int grpID) {
    int profileSize = option->type == 'n' ? 6 : 22; 
    if (option->alnMode != 2) {
        int seqLen = util->seqsLen[T->root->identifier];
        // std::cout << "Alignment Length: " << seqLen << '\n';
        util->profileFreq[grpID] = std::vector<std::vector<float>> (seqLen, std::vector<float> (profileSize, 0.0));
        float totalWeight = 0;
        for (auto sIdx: T->root->msaIdx) totalWeight += T->allNodes[util->seqsName[sIdx]]->weight;
        for (int sIdx: T->root->msaIdx) {
            int storage = util->seqsStorage[sIdx];
            std::string name = util->seqsName[sIdx];
            float w = T->allNodes[name]->weight / totalWeight * T->root->numLeaves;
            tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
            // for (int j = 0; j <  seqLen; ++j) {
                int letterIndex = letterIdx(option->type, toupper(util->alnStorage[storage][sIdx][j]));
                util->profileFreq[grpID][j][letterIndex] += 1.0 * w;                                                                                       util->profileFreq[grpID][j][4]+=1.0*w;
            }
            });
        }
    }
    else {
        int seqLen = util->backboneAln.begin()->second.size();
        util->profileFreq[grpID] = std::vector<std::vector<float>> (seqLen, std::vector<float> (profileSize, 0.0));
        float totalWeight = 0;
        for (auto nodePair: T->allNodes) {
            if (nodePair.second->is_leaf()) totalWeight += nodePair.second->weight;
        }
        for (auto seq: util->placedSeqs) {
            std::string name = seq.first;
            float w = T->allNodes[name]->weight / totalWeight * T->m_numLeaves;
            tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
            // for (int j = 0; j <  seqLen; ++j) {
                int letterIndex = letterIdx(option->type, toupper(seq.second[j]));
                util->profileFreq[grpID][j][letterIndex] += 1.0 * w; 
            }
            });
        }
        for (auto seq: util->backboneAln) {
            std::string name = seq.first;
            float w = T->allNodes[name]->weight / totalWeight * T->m_numLeaves;
            tbb::parallel_for(tbb::blocked_range<int>(0, seqLen), [&](tbb::blocked_range<int> r) {
            for (int j = r.begin(); j < r.end(); ++j) {
            // for (int j = 0; j <  seqLen; ++j) {
                int letterIndex = letterIdx(option->type, toupper(seq.second[j]));
                util->profileFreq[grpID][j][letterIndex] += 1.0 * w; 
            }
            });
        }
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
    // normScore /= matchAvg;
    return normScore;
}   

Tree* pruneTree(Tree* T, std::string& seqFileName) {
    
    auto treePruneStart = std::chrono::high_resolution_clock::now();
    
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    std::unordered_set<std::string> keepSequences;

    size_t seqNum = 0;
    
    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        std::string seqName_full = kseq_rd->name.s;
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        keepSequences.insert(seqName_full);
        if (seqName_noblank != seqName_full) keepSequences.insert(seqName_noblank);
        seqNum++;
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    Tree* pT = new Tree();
    pT->root = new Node(T->root->identifier, T->root->branchLength);
    pT->root->grpID = -1;
    pT->allNodes[pT->root->identifier] = pT->root;
    
    std::unordered_map<std::string, bool> keepNode;
    
    for (const auto& n : T->allNodes) {
        if (n.second->is_leaf()) {
            keepNode[n.second->identifier] = (keepSequences.find(n.second->identifier) != keepSequences.end());
        }
    }

    std::function<bool(Node*)> hasKeepDescendant = [&](Node* node) -> bool {
        if (node->is_leaf()) return keepNode[node->identifier];
        bool keep = false;
        for (Node* child : node->children) {
            if (hasKeepDescendant(child)) {
                keep = true;
            }
        }
        keepNode[node->identifier] = keep;
        return keep;
    };
    hasKeepDescendant(T->root);
    
    std::function<void(Node*, Node*)> buildPrunedTree = [&](Node* origNode, Node* newParent) {
        if (!keepNode[origNode->identifier]) return;
        
        if (origNode->identifier == T->root->identifier) {
            for (Node* child : T->root->children) {
                buildPrunedTree(child, T->root);
            }
            return;
        } 
        else {
            // Process children
            std::vector<Node*> keepChildren;
            for (Node* child : T->allNodes[origNode->identifier]->children) {
                if (keepNode[child->identifier]) {
                    keepChildren.push_back(child);
                }
            }

            // Handle nodes with only one child (merge them)
            if (keepChildren.size() == 0) {
                if (origNode->is_leaf()) {
                    Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                    newNode->grpID = -1;
                    pT->allNodes[newNode->identifier] = newNode;
                }
                else return;
            }
            else if (keepChildren.size() == 1) {
                Node* onlyChild = keepChildren[0];
                float combinedLength = origNode->branchLength;
                while (true) {
                    std::vector<Node*> temp;
                    combinedLength += onlyChild->branchLength;
                    for (Node* child : onlyChild->children) {
                        if (keepNode[child->identifier]) {
                            temp.push_back(child);
                        }
                    }
                    if (temp.size() > 1) {
                        Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                        newChild->grpID = -1;
                        pT->allNodes[newChild->identifier] = newChild;
                        break;
                    }
                    else if (temp.empty()) {
                        if (onlyChild->is_leaf()) {
                            Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                            newChild->grpID = -1;
                            pT->allNodes[newChild->identifier] = newChild;
                            break;
                        }
                        else return;
                    }
                    onlyChild = temp[0];
                }
                // Process the child's children
                for (Node* grandchild : T->allNodes[onlyChild->identifier]->children) {
                    buildPrunedTree(grandchild, T->allNodes[onlyChild->identifier]);
                }
            } 
            else {
                Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                newNode->grpID = -1;
                pT->allNodes[newNode->identifier] = newNode;
                for (Node* child : T->allNodes[origNode->identifier]->children) {
                    buildPrunedTree(child, T->allNodes[origNode->identifier]);
                }
            }
        }
    };
    
    buildPrunedTree(pT->root, nullptr);
    
    // Update tree properties
    pT->m_numLeaves = 0;
    for (const auto& nodePair : pT->allNodes) {
        if (nodePair.second->is_leaf()) {
            pT->m_numLeaves++;
        }
    }
    
    pT->calLeafNum();
    pT->calSeqWeight();

    std::cout << "Number of Leaves: " << T->m_numLeaves << " (before pruning) -> " << pT->m_numLeaves << " (after pruning)\n";
    if (pT->m_numLeaves == 0) {
        std::cerr << "ERROR: No sequences from the input sequence file are found in the tree file.\n";
        exit(1);
    }
    if (pT->m_numLeaves != seqNum) std::cerr << "WARNING: " << (seqNum - pT->m_numLeaves) << " sequences are missing from the tree and will be ignored.\n";
    delete T;
    auto treePruneEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treePruneTime = treePruneEnd - treePruneStart;
    std::cout << "Tree pruned in: " <<  treePruneTime.count() / 1000000 << " ms\n";
    return pT;
}

Tree* pruneTree(Tree* T, std::map<std::string, std::pair<std::string, int>>& seqs) {
    
    auto treePruneStart = std::chrono::high_resolution_clock::now();

    std::unordered_set<std::string> keepSequences;

    size_t seqNum = 0;

    for (auto seq: seqs) {
        keepSequences.insert(seq.first);
        seqNum++;
    }

    Tree* pT = new Tree();
    pT->root = new Node(T->root->identifier, T->root->branchLength);
    pT->root->grpID = T->root->grpID;
    pT->allNodes[pT->root->identifier] = pT->root;
    
    std::unordered_map<std::string, bool> keepNode;
    
    for (const auto& n : T->allNodes) {
        if (n.second->is_leaf()) {
            keepNode[n.second->identifier] = (keepSequences.find(n.second->identifier) != keepSequences.end());
        }
    }

    std::function<bool(Node*)> hasKeepDescendant = [&](Node* node) -> bool {
        if (node->is_leaf()) return keepNode[node->identifier];
        bool keep = false;
        for (Node* child : node->children) {
            if (hasKeepDescendant(child)) {
                keep = true;
            }
        }
        keepNode[node->identifier] = keep;
        return keep;
    };
    hasKeepDescendant(T->root);
    
    std::function<void(Node*, Node*)> buildPrunedTree = [&](Node* origNode, Node* newParent) {
        if (!keepNode[origNode->identifier]) return;
        
        if (origNode->identifier == T->root->identifier) {
            for (Node* child : T->root->children) {
                buildPrunedTree(child, T->root);
            }
            return;
        } 
        else {
            // Process children
            std::vector<Node*> keepChildren;
            for (Node* child : T->allNodes[origNode->identifier]->children) {
                if (keepNode[child->identifier]) {
                    keepChildren.push_back(child);
                }
            }

            // Handle nodes with only one child (merge them)
            if (keepChildren.size() == 0) {
                if (origNode->is_leaf()) {
                    Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                    newNode->grpID = origNode->grpID;
                    pT->allNodes[newNode->identifier] = newNode;
                }
                else return;
            }
            else if (keepChildren.size() == 1) {
                Node* onlyChild = keepChildren[0];
                float combinedLength = origNode->branchLength;
                while (true) {
                    std::vector<Node*> temp;
                    combinedLength += onlyChild->branchLength;
                    for (Node* child : onlyChild->children) {
                        if (keepNode[child->identifier]) {
                            temp.push_back(child);
                        }
                    }
                    if (temp.size() > 1) {
                        Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                        newChild->grpID = onlyChild->grpID;
                        pT->allNodes[newChild->identifier] = newChild;
                        break;
                    }
                    else if (temp.empty()) {
                        if (onlyChild->is_leaf()) {
                            Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                            newChild->grpID = onlyChild->grpID;
                            pT->allNodes[newChild->identifier] = newChild;
                            break;
                        }
                        else return;
                    }
                    onlyChild = temp[0];
                }
                // Process the child's children
                for (Node* grandchild : T->allNodes[onlyChild->identifier]->children) {
                    buildPrunedTree(grandchild, T->allNodes[onlyChild->identifier]);
                }
            } 
            else {
                Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                newNode->grpID = origNode->grpID;
                pT->allNodes[newNode->identifier] = newNode;
                for (Node* child : T->allNodes[origNode->identifier]->children) {
                    buildPrunedTree(child, T->allNodes[origNode->identifier]);
                }
            }
        }
    };
    
    buildPrunedTree(pT->root, nullptr);
    
    // Update tree properties
    pT->m_numLeaves = 0;
    for (const auto& nodePair : pT->allNodes) {
        if (nodePair.second->is_leaf()) {
            pT->m_numLeaves++;
        }
    }
    
    pT->calLeafNum();
    pT->calSeqWeight();

    std::cout << "Number of Leaves: " << T->m_numLeaves << " (before pruning) -> " << pT->m_numLeaves << " (after pruning)\n";
    if (pT->m_numLeaves == 0) {
        std::cerr << "ERROR: No sequences from the input sequence file are found in the tree file.\n";
        exit(1);
    }
    if (pT->m_numLeaves != seqNum) std::cerr << "WARNING: " << (seqNum - pT->m_numLeaves) << " sequences are missing from the tree and will be ignored.\n";
    delete T;
    auto treePruneEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treePruneTime = treePruneEnd - treePruneStart;
    std::cout << "Tree pruned in: " <<  treePruneTime.count() / 1000000 << " ms\n";
    return pT;
}

// Gzip compression function
std::string gzip_compress(const std::string &input) {
    z_stream zs{};
    deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY); // 15+16 = gzip header

    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    zs.avail_in = input.size();

    std::string output;
    output.resize(compressBound(input.size())); 

    zs.next_out = reinterpret_cast<Bytef*>(&output[0]);
    zs.avail_out = output.size();

    int ret = deflate(&zs, Z_FINISH);
    if (ret != Z_STREAM_END) {
        std::cerr << "Compression error: " << ret << "\n";
        deflateEnd(&zs);
        return {};
    }

    output.resize(zs.total_out);
    deflateEnd(&zs);
    return output;
}