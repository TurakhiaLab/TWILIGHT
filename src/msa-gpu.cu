#ifndef UTIL_HPP
#include "../src/util.cuh"
#endif

po::options_description mainDesc("MSA Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,i", po::value<std::string>()->required(), "Input tip sequences - Fasta format (required)")
        ("cpu-num,c",  po::value<int>(), "Number of CPU threads")
        ("gpu-num,g",  po::value<int>(), "Number of GPUs")
        ("max-leaves,l",  po::value<int>(), "Maximum number of leaves per sub-subtree")
        ("max-subtree-size,m", po::value<int>()->default_value(1000000), "Maximum number of leaves per subtree")
        
        ("output,o", po::value<std::string>(), "Output file name")
        ("match",      po::value<paramType>()->default_value(18), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(-8), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-100), "Gap open penalty")
        ("gap-close",  po::value<paramType>()->default_value(-100), "Gap close penalty")
        ("gap-extend", po::value<paramType>()->default_value(-5), "Gap extend penalty")
        ("trans",      po::value<paramType>()->default_value(2), "The ratio of transitions to transversions")
        ("pam-n",      po::value<int>()->default_value(200), "'n' in the PAM-n matrix")
        ("xdrop",      po::value<paramType>()->default_value(600), "X-drop value")
        ("scoring-matrix", po::value<int>()->default_value(0), "0: simple, 1: kiruma, 2: user defined")
        ("temp-dir", po::value<std::string>(), "Directory for storing temporary files")
        ("gpu-index", po::value<std::string>(), "Specify the GPU index, separated by commas. Ex. 0,2,3")
        ("gappy-column", po::value<float>()->default_value(1), "If the proportion of gaps in a column exceeds this value, the column will be defined as a gappy column.")
        ("gappy-length", po::value<float>(), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("read-batches", "Read sequences in batches and create temporary files")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment")
        ("debug", "Enable debug on the final alignment")
        ("help,h", "Print help messages");

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

    msa::option* option = new msa::option;
    setOptions(vm, option);

    Params* param = setParameters(vm);

    // Partition tree into subtrees
    Tree* T = readNewick(vm);
    // printTree(T->root, -1);
    paritionInfo_t* P = new paritionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);
    // exit(1);
    // int subtreeCount = 0;
    // for (auto p: P->partitionsRoot) {
    //     ++subtreeCount;
    //     std::cout << subtreeCount << '\t' << T->allNodes[p.first]->getNumLeaves() << '\n';
    // }
    // return;
    
    // Define MSA utility
    msa::utility* util = new msa::utility;

    
    std::unordered_map<std::string, std::string> beforeAln;
        
    // read sequences
    if (!option->readBatches) {
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
        if (option->debug) {
            for (auto s: T->allNodes) {
                if (s.second->is_leaf()) {
                    std::string seqName = s.first;
                    int sIdx = util->seqsIdx[s.first];
                    int storage = util->seqsStorage[sIdx];
                    std::string r = "";
                    int j = 0;
                    while (util->alnStorage[storage][sIdx][j] != 0) {
                        if (util->alnStorage[storage][sIdx][j] != '-') {
                            r += util->alnStorage[storage][sIdx][j];
                        }
                        ++j;
                    }
                    beforeAln[seqName] = r;
                }
            }   
        }
    }
    else {
        readSequencesNoutputTemp(vm, T, P);
    }
    // return;
    for (auto subRoot: P->partitionsRoot) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        int subtree = T->allNodes[subRoot.first]->grpID;
        std::cout << "Start processing subtree No. " << subtree << '\n';
        
        Tree* subT;
        if (!option->readBatches) {
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
        
        paritionInfo_t * subP = new paritionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(subT->root, subP);
       
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " trees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        
        auto msaStart = std::chrono::high_resolution_clock::now();
        Tree* newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        // std::cout << newSubT->root->getNumLeaves() << '\n';
        msaOnSubtree(subT, util, option, subP, *param);
        auto msaEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        std::cout << "MSA on sub-subtree in " <<  msaTime.count() / 1000000000 << " s\n";
        
        if (subP->partitionsRoot.size() > 1) {
            // Align adjacent sub-subtrees to create overlap alignment
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees(subT, newSubT, util, option, *param);
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
        
        for (auto sIdx: subT->root->msaIdx) T->allNodes[subT->root->identifier]->msaIdx.push_back(sIdx);
        // for (auto node: subT->allNodes) delete node.second;
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeAlnFile = tempDir + '/' + subtreeFileName + ".temp.aln";
            std::string subtreeFreqFile = tempDir + '/' + subtreeFileName + ".freq.txt";
            outputFile(subtreeAlnFile, util, subT, 0);
            outputFreq(subtreeFreqFile, util, subT, subtree);
        }
        delete subT;
        delete subP;
        delete newSubT;
    }

    if (P->partitionsRoot.size() > 1) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        util->nowProcess = 2; // merge subtrees
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            readFreq(tempDir, T, P, util);
        }
        else {
            getFreq(T, P, util);
        }
        alignSubtrees(T, newT, util, option, *param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Aligned adjacent subtrees in " <<  alnTime.count() / 1000000 << " ms\n";
        auto mergeStart = std::chrono::high_resolution_clock::now();
        mergeSubtrees (T, newT, util);
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        int totalSeqs = 0;
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            outputFinal (tempDir, T, P, util, totalSeqs);
        }
        else totalSeqs = T->root->msaIdx.size();
        std::cout << "Merge " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    }


    // post-alignment debugging
    if (option->debug) {
        auto dbgStart = std::chrono::high_resolution_clock::now();
        int alnLen = 0;
        // for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
        bool theFirst = true;
        for (auto s: beforeAln) {
            int sIdx = util->seqsIdx[s.first];
            int storage = util->seqsStorage[sIdx];
            std::string r = "";
            int offset = 0;
            while (util->alnStorage[storage][sIdx][offset] != 0) {
                if (util->alnStorage[storage][sIdx][offset] != '-') {
                    r += util->alnStorage[storage][sIdx][offset];
                }
                ++offset;
            }
            if (theFirst) {alnLen = offset; theFirst = false;}
            else {
                if (alnLen != offset) printf("seq: %s, the sequence length (%d) did not match (%d)\n", s.first.c_str(), offset, alnLen);
            }
            if (r != s.second) {
                printf("seq: %s, the sequence did not match\n", s.first.c_str());
                if (r.size() != s.second.size()) {
                    std::cout << "Wrong length. " << r.size() << '/' << s.second.size() << ".\n";
                }
                for (int i = 0; i < s.second.size(); ++i) {
                    if (r[i] != s.second[i]) {
                        std::cout << "Mismatch at position " << i << '\n';
                        break;
                    }
                }                
            }
            
        }
        beforeAln.clear();
        auto dbgEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
        std::cout << "Completed checking " << T->m_numLeaves << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
    }
    
    // Calculate sum-of-pairs score
    if (vm.count("sum-of-pairs-score")) {
        auto spStart = std::chrono::high_resolution_clock::now();
        double score = getSPScore_gpu(util, *param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
    }
    
    // output MSA
    if (vm.count("output")) {
        std::string outFile = vm["output"].as<std::string>();
        if (outFile == "") outFile = "output.aln";
        auto outStart = std::chrono::high_resolution_clock::now();
        std::string subtreeFreqFile = outFile + ".freq.txt";
        outputFreq(subtreeFreqFile, util, T, -1);
        outputFile(outFile, util, T, -1);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    }
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return;
}
