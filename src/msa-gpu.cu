#ifndef UTIL_HPP
#include "../src/util.cuh"
#endif

#include <sys/stat.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>

KSEQ_INIT2(, gzFile, gzread)

namespace po = boost::program_options;

po::options_description mainDesc("MSA Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,i", po::value<std::string>()->required(), "Input tip sequences - Fasta format (required)")
        ("cpu-num,c",  po::value<int>(), "Number of CPU threads")
        ("gpu-num,g",  po::value<int>(), "Number of GPUs")
        ("max-leaves,l",  po::value<int>()->default_value(0), "Maximum number of leaves per sub-subtree")
        ("max-subtree-size,m", po::value<int>()->default_value(1000000), "Maximum number of leaves per subtree")
        
        ("output,o", po::value<std::string>()->default_value(""), "Output file name")
        ("match",      po::value<paramType>()->default_value(2), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(0), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-3), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-1), "Gap extend penalty")
        ("trans",      po::value<paramType>()->default_value(0), "Transition score")
        ("xdrop",      po::value<paramType>()->default_value(0), "X-drop value")
        ("temp-dir", po::value<std::string>(), "Directory for storing temporary files")
        ("read-batches", "Read sequences in batches and create temporary files")
        ("user-define-parameters", "Using user defined parameters. Please modify align.cuh")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment")
        ("debug", "Enable debug and print debug messages")
        ("help,h", "Print help messages");

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
        std::remove(tempAlnFile.c_str());
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


int main(int argc, char** argv) {

    auto mainStart = std::chrono::high_resolution_clock::now();
    
    parseArguments(argc, argv);
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << mainDesc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }



    // Partition tree into subtrees
    Tree* T = readNewick(vm);
    int maxSubtreeSize = vm["max-subtree-size"].as<int>();
    paritionInfo_t* P = new paritionInfo_t(maxSubtreeSize, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);
    // Define MSA utility
    msa::utility* util = new msa::utility;

    // Print hardware information
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    int cpuNum = (vm.count("cpu-num")) ? vm["cpu-num"].as<int>() : maxCpuThreads;
    printf("Maximum available CPU threads: %d. Using %d CPU threads.\n", maxCpuThreads, cpuNum);
    tbb::task_scheduler_init init(cpuNum);
    int maxGpuNum;
    cudaGetDeviceCount(&maxGpuNum);
    int gpuNum = (vm.count("gpu-num")) ? vm["gpu-num"].as<int>() : maxGpuNum;
    gpuNum = (gpuNum == 0) ? maxGpuNum : gpuNum;
    util->gpuNum = gpuNum;
    printf("Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);


    // Set parameters
    paramType mat = vm["match"].as<paramType>();
    paramType mis = vm["mismatch"].as<paramType>();
    paramType gapOp = vm["gap-open"].as<paramType>();
    paramType gapEx = vm["gap-extend"].as<paramType>();
    paramType trans = vm["trans"].as<paramType>();
    paramType xdrop = vm["xdrop"].as<paramType>();
    float gapOp_int =  round(gapOp);
    float gapEx_int =  round(gapEx);
    float xdrop_int =  round(xdrop);
    float gapOp_diff = gapOp-gapOp_int;
    float gapEx_diff = gapEx-gapEx_int;
    float xdrop_diff = xdrop-xdrop_int;
    if (gapOp_diff != 0) printf("WARNING: Floating point gap open penalty is not allowed, %.0f is used intead of %f.\n", gapOp_int, gapOp);
    if (gapEx_diff != 0) printf("WARNING: Floating point gap extend penalty is not allowed, %.0f is used intead of %f.\n", gapEx_int, gapEx);
    if (xdrop_diff != 0) printf("WARNING: Floating point xdrop is not allowed, %.0f is used intead of %f.\n", xdrop_int, xdrop);
    if (xdrop_int == 0) xdrop_int = round((FRONT_WAVE_LEN/3)*(-gapEx_int));
    if (trans == 0) trans = mis + (mat-mis)/2;
    bool userDefine = vm.count("user-define-parameters");
    Params param(mat,mis,trans,gapOp_int,gapEx_int,xdrop_int,userDefine);


    int maxSubSubtreeSize = vm["max-leaves"].as<int>();
    if (maxSubSubtreeSize == 0) maxSubSubtreeSize = INT_MAX;
    bool debug = vm.count("debug");
    bool batches = vm.count("read-batches");
    std::vector<std::string> beforeAln, afterAln;
        
    // read sequences
    if (!batches) {
        readSequences(vm, util, T);
        if (T->m_numLeaves != util->seqsIdx.size()) {
            fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %lu)\n", T->m_numLeaves, util->seqsIdx.size()); 
            exit(1);
        }
        for (auto it = T->allNodes.begin(); it != T->allNodes.end(); ++it) {
            if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                printf("Missing Sequence %s.\n", it->first.c_str());
                exit(1);
            }
        }
        if (debug) {
            for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
                std::string r = "";
                int j = 0;
                while (util->seqs[sIdx][j] != 0) {
                    if (util->seqs[sIdx][j] != '-') {
                        r += util->seqs[sIdx][j];
                    }
                    ++j;
                }
                beforeAln.push_back(r);
            }   
        }
    }
    else {
        readSequencesNoutputTemp(vm, T, P);
    }

    for (auto subRoot: P->partitionsRoot) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        int subtree = T->allNodes[subRoot.first]->grpID;
        std::cout << "Start processing subtree No. " << subtree << '\n';
        
        Tree* subT;
        if (!batches) {
            subT = new Tree(subRoot.second.first);
            util->setSubtreeIdx(subtree);
        }
        else {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeTreeFile = tempDir + '/' + subtreeFileName + ".nwk";
            std::string subtreeSeqFile = tempDir + '/' + subtreeFileName + ".raw.fa";
            subT = readNewick(subtreeTreeFile);
            paritionInfo_t* tempP = new paritionInfo_t(INT_MAX, 0, 0, "centroid"); 
            partitionTree(subT->root, tempP);
            delete tempP;
            readSequences(subtreeSeqFile, util, subT);
            util->setSubtreeIdx(0);
            if (subT->m_numLeaves != util->seqsIdx.size()) {
                fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %lu)\n", subT->m_numLeaves, util->seqsIdx.size()); 
                exit(1);
            }
            for (auto it = subT->allNodes.begin(); it != subT->allNodes.end(); ++it) {
                if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                    printf("Missing Sequence %s.\n", it->first.c_str());
                    exit(1);
                }
            }
        }
        // printTree(subT->root, -1);
        std::cout << "Subtree No." << subtree << " contains "<< subT->m_numLeaves << " sequences.\n";
        // Partition subtree into sub-subtrees
        
        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        
        util->alnMalloc(util->seqLen);
        paritionInfo_t * subP = new paritionInfo_t(maxSubSubtreeSize, 0, 0, "centroid");
        partitionTree(subT->root, subP);
       
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Partition the subtree in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        auto msaStart = std::chrono::high_resolution_clock::now();
        Tree* newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        // std::cout << newSubT->root->getNumLeaves() << '\n';
        msaOnSubtree(subT, util, subP, param);
        auto msaEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        std::cout << "Progressive alignment in " <<  msaTime.count() / 1000000000 << " s\n";
        
        if (subP->partitionsRoot.size() > 1) {
            // Align adjacent sub-subtrees to create overlap alignment
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees(subT, newSubT, util, param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            std::cout << "Aligned adjacent sub-subtrees in " <<  alnTime.count() / 1000000 << " ms\n";

            auto mergeStart = std::chrono::high_resolution_clock::now();
            mergeSubtrees (subT, newSubT, util);
            auto mergeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
            int totalSeqs = subT->root->msaIdx.size();
            std::cout << "Merge " << newSubT->allNodes.size() << " sub-subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            std::cout << "Finished the alignment on subtree No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        util->updateSeqs();
        for (auto sIdx: subT->root->msaIdx) T->allNodes[subT->root->identifier]->msaIdx.push_back(sIdx);
        // for (auto node: subT->allNodes) delete node.second;
        if (batches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeAlnFile = tempDir + '/' + subtreeFileName + ".temp.aln";
            std::string subtreeFreqFile = tempDir + '/' + subtreeFileName + ".freq.txt";
            outputFile(subtreeAlnFile, util, subT, 0);
            outputFreq(subtreeFreqFile, util, subT, subtree);
            util->reset();
        }
        delete subT;
        delete subP;
        delete newSubT;
    }

    if (P->partitionsRoot.size() > 1) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        util->nowProcess = 1; // merge subtrees
        if (batches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            readFreq(tempDir, T, P, util);
        }
        else {
            getFreq(T, P, util);
        }
        alignSubtrees(T, newT, util, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Aligned adjacent subtrees in " <<  alnTime.count() / 1000000 << " ms\n";
        auto mergeStart = std::chrono::high_resolution_clock::now();
        mergeSubtrees (T, newT, util);
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        int totalSeqs = 0;
        if (batches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            outputFinal (tempDir, T, P, util, totalSeqs);
        }
        else totalSeqs = T->root->msaIdx.size();
        std::cout << "Merge " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    }


    // post-alignment debugging
    if (debug) {
        auto dbgStart = std::chrono::high_resolution_clock::now();
        int alnLen = 0;
        for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
            std::string r = "";
            int offset = 0;
            while (util->seqs[sIdx][offset] != 0) {
                if (util->seqs[sIdx][offset] != '-') {
                    r += util->seqs[sIdx][offset];
                }
                ++offset;
            }
            if (sIdx == 0) alnLen = offset;
            else {
                if (alnLen != offset) printf("No: %d, the sequence length (%d) did not match (%d)\n", sIdx, offset, alnLen);
            }
            afterAln.push_back(r);
        }
        for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
            if (beforeAln[sIdx] != afterAln[sIdx]) {
                printf("No: %d, the sequence did not match\n", sIdx);
            }
        }
        auto dbgEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
        std::cout << "Completed checking " << T->m_numLeaves << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
    }
    
    // Calculate sum-of-pairs score
    if (vm.count("sum-of-pairs-score")) {
        auto spStart = std::chrono::high_resolution_clock::now();
        double score = getSPScore_gpu(util, param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
    }
    
    // output MSA
    if (vm["output"].as<std::string>() != "") {
        std::string outFile = vm["output"].as<std::string>();
        auto outStart = std::chrono::high_resolution_clock::now();
        std::string subtreeFreqFile = outFile + ".freq.txt";
        outputFile(outFile, util, T, -1);
        outputFreq(subtreeFreqFile, util, T, -1);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    }
    // release memory
    util->seqsFree();
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return;
}

    // paritionInfo_t * partition = new paritionInfo_t(maxSubtreeSize, 0, 0, "centroid");
    
    // partitionTree(T->root, partition);
    
   
    
    // // int totalLeaves = 0;
    // // for (auto p: partition->partitionsRoot) {
    // //     printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     totalLeaves += p.second.second;
    // // }
    // // std::cout << "Total Leaves: " << totalLeaves << '\n';

    // Tree* newT = reconsturctTree(T->root, partition->partitionsRoot);    
    
    // // Start MSA on subtrees
    // auto msaStart = std::chrono::high_resolution_clock::now();
    // std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    // for (auto &p: partition->partitionsRoot) {
    //     std::stack<Node*> msaStack;
    //     getPostOrderList(p.second.first, msaStack);
        
    //     std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
    //     int grpID = p.second.first->grpID;
    //     getMsaHierachy(subhier, msaStack, grpID, 0);
    //     for (auto h: subhier) {
    //         while (hier.size() < h.second+1) {
    //             std::vector<std::pair<Node*, Node*>> temp;
    //             hier.push_back(temp);
    //         }
    //         hier[h.second].push_back(h.first);
    //     }
    // }
    // int level = 0;
    
    // for (auto m: hier) {
    //     // std::cout << "Aln level: " << level << '\n';
    //     auto alnStart = std::chrono::high_resolution_clock::now();
    //     msaPostOrderTraversal_multigpu(T, m, util, param);
    //     auto alnEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    //     if (m.size() > 1) std::cout << "Level "<< level << ", " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
    //     else              std::cout << "Level "<< level << ", " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
    //     ++level;
    // }
    // auto msaEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds msaTime = msaEnd - msaStart;
    // std::cout << "MSA in " <<  msaTime.count() / 1000000000 << " s\n";

    // // Push MSA to partition roots
    // for (auto p: partition->partitionsRoot) {
    //     std::stack<Node*> msaStack;
    //     getPostOrderList(p.second.first, msaStack);
    //     std::vector<Node*> msaArray;
    //     while (!msaStack.empty()) {
    //         msaArray.push_back(msaStack.top());
    //         msaStack.pop();
    //     }
    //     if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
    //         if (msaArray.size() == 2) {
    //             T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
    //             util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
    //             break;
    //         }
    //         for (int m = msaArray.size()-2; m >=0; --m) {
    //             if (msaArray[m]->msaIdx.size()>0) {
    //                 T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
    //                 util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
    //                 break;
    //             }
    //         }
    //     }
    // }
    

    // // type 1 alignment
    // // int totalLeaves = 0;
    // // for (auto n: newT->allNodes) {
    // //     std::cout << n.first << ':' << T->allNodes[n.first]->msaIdx.size() << '\n';
    // //     totalLeaves += T->allNodes[n.first]->msaIdx.size();
    // //     // std::cout << '\n';
    // // }
    // // std::cout << "totalLeaves: " << totalLeaves << '\n';
    
    // // int totalLeaves = 0;
    // // for (auto p: partition->partitionsRoot) {
    // //     // std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\n';
    // //     // printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     if (T->allNodes[p.first]->msaIdx.size() != p.second.second) {
    // //         std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\t';
    // //         printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     }
    // //     totalLeaves += T->allNodes[p.first]->msaIdx.size();
    // // }
    // // std::cout << "Total Leaves: " << totalLeaves << '\n';
    

    // // Create overlapped alignment (MSA between parent and children)
    // if (partition->partitionsRoot.size() > 1) {
    //     auto alnStart = std::chrono::high_resolution_clock::now();
    //     std::vector<std::pair<Node*, Node*>> type1Aln;
    //     // Copy from util->seqBuf to Node->msa
    //     for (auto n: newT->allNodes) {
    //         for (auto m: n.second->children) {
    //             if (newT->allNodes[m->identifier]->grpID == newT->allNodes[n.second->identifier]->grpID) {
    //                 type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
    //             }
    //         }
    //         T->allNodes[n.second->identifier]->msa.clear();
    //         for (int sIdx: T->allNodes[n.second->identifier]->msaIdx) {
    //             std::string r = "";
    //             int storage = util->seqsStorage[sIdx];
    //             int offset = 0;
    //             while (util->seqBuf[storage][sIdx][offset] != 0) {
    //                 r += util->seqBuf[storage][sIdx][offset];
    //                 ++offset;
    //             }
    //             T->allNodes[n.second->identifier]->msa.push_back(r);
    //         }
            
    //     }
    //     // for (auto n: type1Aln) {
    //     //     std::cout << n.first->identifier << '(' << T->allNodes[n.first->identifier]->msa.size() << ')'
    //     //               << n.second->identifier << '(' << T->allNodes[n.second->identifier]->msa.size() << ")\n";
    //     // }
    //     createOverlapMSA(T, type1Aln, util, param);
    //     auto alnEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    //     if (type1Aln.size() > 1) std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pairs in " <<  alnTime.count() / 1000000 << " ms\n";
    //     else                     std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pair in "  <<  alnTime.count() / 1000000 << " ms\n";
    //     hier.clear();

    //     // Create a new tree with all roots of subtrees
    //     auto mergeStart = std::chrono::high_resolution_clock::now();
    //     paritionInfo_t * newPartition = new paritionInfo_t(std::numeric_limits<size_t>::max(), 0, 0, "centroid");
    //     partitionTree(newT->root, newPartition); 
    //     // printTree(newT->root);
    //     // newPartition->partitionsRoot[midNode->identifier] = std::make_pair(newT->allNodes[midNode->identifier], newT->allNodes[midNode->identifier]->getNumNodes());
    //     // newPartition->partitionsRoot[newT->root->identifier] = std::make_pair(newT->root, newT->root->getNumNodes());
    //     std::vector<std::pair<Node*, Node*>> mergePairs;
    //     std::vector<std::pair<std::pair<Node*, Node*>, float>> sortedMergePairs;
    //     for (auto n: newT->allNodes) {
    //         if (n.second->children.size() > 1) {
    //             mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
    //             for (int i = 1; i < n.second->children.size(); ++i) {
    //                 mergePairs.push_back(std::make_pair(n.second->children[0], n.second->children[i]));
    //             }
    //         }
    //         else if (n.second->children.size() == 1) {
    //             mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
    //         }
    //     }
        
    //     std::vector<std::pair<Node*, Node*>> singleLevel;
    //     std::map<std::string, char> addedNodes;
    //     // std::cout << "\nDEBUG:MERGE PAIRS\n";
    //     // for (auto n: mergePairs) {
    //     //     std::cout << n.first->identifier << ',' << n.second->identifier << ':' << n.second->branchLength << '\n';
    //     // }
    //     while (true) {
    //         auto roundStart = std::chrono::high_resolution_clock::now();
    //         addedNodes.clear();
    //         singleLevel.clear();
    //         for (auto it = mergePairs.begin(); it != mergePairs.end();) {
    //             Node* a = it->first;
    //             Node* b = it->second;
    //             if ((a->parent != nullptr && b->parent != nullptr) &&  
    //                 (addedNodes.find(a->identifier) == addedNodes.end() && addedNodes.find(b->identifier) == addedNodes.end())) {
    //                 singleLevel.push_back(std::make_pair(a, b));
    //                 for (auto id: T->allNodes[a->identifier]->msa) addedNodes[id] = 0;
    //                 for (auto id: T->allNodes[b->identifier]->msa) addedNodes[id] = 0;
    //                 mergePairs.erase(it);
    //             }
    //             else {
    //                 ++it;
    //             }
    //         }
    //         bool breakLoop = false;
    //         if (singleLevel.empty()) {
    //             for (auto mp: mergePairs) {
    //                 singleLevel.push_back(mp);
    //                 // std::cout << mp.first->identifier << ':' << mp.second->identifier << '\n';
    //             }
    //             breakLoop = true;
    //         }
    //         transitivityMerge_cpu_mod(T, newT, singleLevel, util);
    //         // mergeLevels.push_back(singleLevel);
    //         auto roundEnd = std::chrono::high_resolution_clock::now();
    //         std::chrono::nanoseconds roundTime = roundEnd - roundStart;
    //         std::cout << "Merge "<< singleLevel.size() << " edges in " << roundTime.count() / 1000000 << " ms\n";
    //         if (breakLoop) break;
    //     }
    //     // transitivityMerge_cpu(T, mergePairs, util);
    //     // transitivityMerge_cpu_mod(T, newT, mergePairs, util);
    //     auto mergeEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
    //     int totalSeqs = 0;
    //     for (auto &p: newPartition->partitionsRoot) totalSeqs += T->allNodes[p.first]->msaIdx.size();
    //     // std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << T->root->msaIdx.size() << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    //     std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";

    //     T->root->msa.clear();
    //     for (int sIdx: T->root->msaIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }
    // else {
    //     for (int sIdx: T->root->msaIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }   

    // }
    
    // bool regressive = false;
    // if (regressive) {
    //     int clusterSize = vm["max-leaves"].as<int>();
    //     if (clusterSize == 0) clusterSize = 1000;
    //     std::cout << "Cluster Size: " << clusterSize << '\n';
    //     getLongestDescendent(T, util);
    //     std::cout << "Finish longest descendent\n";
    //     std::map<int, std::vector<std::vector<Node*>>> subMSAs;
    //     getSubMSAs(subMSAs, T, clusterSize);
    //     std::cout << "Finish subMSAs, subMSA size: " << subMSAs.size() << "\n";
    //     for (auto level = subMSAs.begin(); level != subMSAs.end(); ++level) {
    //         std::cout << "Level: " << level->first << " Size: " << level->second.size() << '\n';
    //         std::vector<std::vector<std::pair<Node*, Node*>>> alnPairs;
    //         getAlnPairs(alnPairs, level->second);
    //         for (auto pairs: alnPairs) {
    //             msaPostOrderTraversal_multigpu_regressive(T, pairs, util, param);
    //         }
    //         std::vector<Node*> nodes;
    //         for (auto msa: level->second) {
    //             nodes.push_back(msa[0]);
    //         }
    //         storeMSA(T, nodes, util, level->first);
    //         resetSeqMem(nodes, util);
    //         for (auto n: T->allNodes) {
    //             n.second->msaIdx.clear();
    //             n.second->msa.clear();
    //         }
    //     }
    //     std::cout << "Finish MSA\n";
    //     transitivityMerge_regressive(T, subMSAs, util);
    //     std::cout << "Finish Merger\n";
    //     T->root->msa.clear();
    //     for (int sIdx = 0; sIdx < util->memNum; ++sIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }

    


    // Calculate Sum-of-pairs score
    
    // auto spStart = std::chrono::high_resolution_clock::now();
    // double score;
    // if (calSP == "True" || calSP == "TRUE" || calSP == "true" || calSP == "T" || calSP == "t") {
    //     score = getSPScore_gpu(T->root->msa, util, param);
    //     // getSPScore_cpu(T->root->msa, param);
    //     auto spEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds spTime = spEnd - spStart;
    //     std::cout << "SP-score: " << score << ". Runtime " << spTime.count() / 1000000 << " ms\n";
    // }


    // // output file
    // std::string outFile = vm["output"].as<std::string>();
    // // if (outFile == "") outFile = "output.aln"; // default output ifle name
    // if (outFile != "") {
    //     auto outStart = std::chrono::high_resolution_clock::now();
    //     outputFile(outFile, util, T, -1);
    //     auto outEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds outTime = outEnd - outStart;
    //     std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    // }
    
    
    // util->seqFree();
    // auto mainEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    // std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    // return 0;
// }