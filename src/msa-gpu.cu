#ifndef UTIL_HPP
#include "util.cuh"
#endif

#ifndef PROCESS_HPP
#include "process.cuh"
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
        ("max-leaves,l",  po::value<int>(), "Maximum number of leaves per sub-subtree, used for transitivity merger")
        ("max-subtree-size,m", po::value<int>(), "Maximum number of leaves per subtree")
        ("output,o", po::value<std::string>(), "Output file name")
        ("match",      po::value<paramType>()->default_value(18), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(-8), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-50), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-5), "Gap extend penalty")
        ("trans",      po::value<paramType>()->default_value(2), "The ratio of transitions to transversions")
        ("pam-n",      po::value<int>()->default_value(200), "'n' in the PAM-n matrix")
        ("xdrop",      po::value<paramType>()->default_value(600), "X-drop value")
        ("scoring-matrix", po::value<int>()->default_value(0), "0: simple, 1: kiruma, 2: user defined")
        ("output-type", po::value<std::string>()->default_value("FASTA"), "FASTA or CIGAR")
        ("temp-dir", po::value<std::string>(), "Directory for storing temporary files")
        ("gpu-index", po::value<std::string>(), "Specify the GPU index, separated by commas. Ex. 0,2,3")
        ("gappy-vertical", po::value<float>()->default_value(1), "If the proportion of gaps in a column exceeds this value, the column will be defined as a gappy column.")
        ("gappy-horizon", po::value<float>(), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment")
        ("debug", "Enable debug on the final alignment")
        ("cpu-only", "Only using CPU to run the program")
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
    msa::utility* util = new msa::utility;

    Params* param = setParameters(vm);

    // Partition tree into subtrees
    Tree* T = readNewick(vm);
    // printTree(T->root, -1);
    paritionInfo_t* P = new paritionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);
        
    // read sequences
    if (P->partitionsRoot.size() == 1) readSequences(vm, util, option, T);
    else readSequencesNoutputTemp(vm, T, P, util, option);
    // return;
    for (auto subRoot: P->partitionsRoot) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        int subtree = T->allNodes[subRoot.first]->grpID;
        std::cout << "Start processing subtree No. " << subtree << '\n';
        Tree* subT;
        if (P->partitionsRoot.size() == 1) {
            subT = new Tree(subRoot.second.first);
            util->setSubtreeIdx(subtree);
        }
        else {
            std::string tempDir = option->tempDir;
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeTreeFile = tempDir + '/' + subtreeFileName + ".nwk";
            std::string subtreeSeqFile = tempDir + '/' + subtreeFileName + ".raw.fa";
            subT = readNewick(subtreeTreeFile);
            paritionInfo_t* tempP = new paritionInfo_t(INT_MAX, 0, 0, "centroid"); 
            partitionTree(subT->root, tempP);
            // printTree(subT->root, -1);
            delete tempP;
            util->clearAll();
            readSequences(subtreeSeqFile, util, option, subT);
            util->setSubtreeIdx(subtree);
            if (subT->root->numLeaves != util->seqsIdx.size()) {
                fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %lu)\n", subT->m_numLeaves, util->seqsIdx.size()); 
                exit(1);
            }
            // int s = 0;
            for (auto it = subT->allNodes.begin(); it != subT->allNodes.end(); ++it) {
                if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                    printf("Missing Sequence %s.\n", it->first.c_str());
                    exit(1);
                }
            }
        }
        std::cout << "Subtree No." << subtree << " contains "<< subT->m_numLeaves << " sequences.\n";
        
        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        
        paritionInfo_t * subP = new paritionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(subT->root, subP);
       
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " sub-subtrees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        
        auto msaStart = std::chrono::high_resolution_clock::now();
        // getSeqsFreq(subT, util);
        // auto freqEnd = std::chrono::high_resolution_clock::now();
        // std::chrono::nanoseconds freqTime = freqEnd - msaStart;
        // std::cout << "Finish getting sequence frequency in " << freqTime.count() / 1000000 << " ms.\n";
        Tree* newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        msaOnSubtree(subT, util, option, subP, *param);
        auto msaEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        std::cout << "MSA on sub-subtree " << subtree << " in " <<  msaTime.count() / 1000000000 << " s\n";
        
        if (subP->partitionsRoot.size() > 1) {
            // Align adjacent sub-subtrees to create overlap alignment
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees(subT, newSubT, util, option, *param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            std::cout << "Aligned adjacent sub-subtree in " <<  alnTime.count() / 1000000 << " ms\n";

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
        // post-alignment debugging
        if (option->debug) {
            auto dbgStart = std::chrono::high_resolution_clock::now();
            util->debug();
            auto dbgEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
            std::cout << "Completed checking " << subT->m_numLeaves << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
        }
        // Calculate sum-of-pairs score
        if (vm.count("sum-of-pairs-score")) {
            auto spStart = std::chrono::high_resolution_clock::now();
            double score = getSPScore_gpu(util, *param);
            auto spEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds spTime = spEnd - spStart;
            std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
        }
        // for (auto node: subT->allNodes) delete node.second;
        if (P->partitionsRoot.size() > 1) {
            std::string tempDir = option->tempDir; 
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeAlnFile = tempDir + '/' + subtreeFileName + ".temp.aln";
            std::string subtreeFreqFile = tempDir + '/' + subtreeFileName + ".freq.txt";
            outputFreq(subtreeFreqFile, util, subT, subtree);
            outputAln(subtreeAlnFile, util, option, subT, subtree);
        }
        delete subT;
        delete subP;
        delete newSubT;
    }

    if (P->partitionsRoot.size() > 1) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        util->nowProcess = 2; // merge subtrees
        readFreq(option->tempDir, T, P, util);
        alignSubtrees(T, newT, util, option, *param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Aligned adjacent subtrees in " <<  alnTime.count() / 1000000 << " ms\n";
        auto mergeStart = std::chrono::high_resolution_clock::now();
        mergeSubtrees (T, newT, util);
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        int totalSeqs = 0;
        outputFinal (option->tempDir, T, P, util, option, totalSeqs);
        std::cout << "Merge " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    }
    
    
    // output MSA
    if (vm.count("output")) {
        if (vm.count("max-subtree-size")) std::cout << "Output files are already in " << option->tempDir << ".\n";
        else {
            std::string outFile = vm["output"].as<std::string>();
            if (outFile == "") outFile = "output.aln";
            auto outStart = std::chrono::high_resolution_clock::now();
            outputAln(outFile, util, option, T, -1);
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
        }
    }
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return;
}
