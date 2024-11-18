#ifndef UTIL_HPP
#include "util.hpp"
#endif

#ifndef MSA_HPP
#include "msa-cpu.hpp"
#endif

po::options_description mainDesc("TWILIGHT Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>(), "Initial Tree - Newick format (required if building MSA from raw sequences)")
        ("sequences,i", po::value<std::string>(), "Input tip sequences - Fasta format")
        ("files,f", po::value<std::string>(), "Path to the directory containing all MSA files. MSA files should be in Fasta format")
        ("cpu-num,c",  po::value<int>(), "Number of CPU cores")
        ("gpu-num,g",  po::value<int>(), "Number of GPUs")
        ("max-leaves,l",  po::value<int>(), "Maximum number of leaves per sub-subtree, used for transitivity merger")
        ("max-subtree,m", po::value<int>(), "Maximum number of leaves per subtree")
        ("temp-dir,d", po::value<std::string>(), "Directory for storing temporary files")
        ("output,o", po::value<std::string>(), "Output file name")
        ("gappy-vertical,v", po::value<float>()->default_value(0.95), "If the proportion of gaps in a column exceeds this value, the column will be defined as a gappy column. Set to 1 to disable this feature.")
        ("gappy-horizon,z", po::value<float>()->default_value(1), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment, have to be used with -o option")
        ("wildcard,w", "Treat unknown 'N' as a wildcard and aligned to usual letters")
        ("print-detail", "Print out every detail")
        ("match",      po::value<paramType>()->default_value(18), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(-8), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-50), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-5), "Gap extend penalty")
        ("xdrop",      po::value<paramType>()->default_value(600), "X-drop value")
        ("scoring-matrix", po::value<int>()->default_value(0), "0: default, 1: kiruma, 2: user defined")
        ("pam-n",      po::value<int>()->default_value(200), "'n' in the PAM-n matrix. Should be used with kimura scoring matrix")
        ("trans",      po::value<paramType>()->default_value(2), "The ratio of transitions to transversions")
        ("output-type", po::value<std::string>()->default_value("FASTA"), "FASTA or CIGAR, CIGAR stands for CIGAR-like compressed format")
        ("merge-subtrees", po::value<std::string>()->default_value("p"), "t: transitivity merger, p: profile alignment")
        ("gpu-index", po::value<std::string>(), "Specify the GPU index, separated by commas. Ex. 0,2,3")
        ("delete-temp", "delete temporary subtree folder and files.")
        ("psgop", po::value<std::string>(), "'y' for enabling and 'n' for disabling position-specific gap open penalty. If not specified, it will be detected automatically.")
        ("output-subtrees", "writing the tree file of each subtree.")
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
    if(vm.count("help")) {
        std::cerr << mainDesc << std::endl;
        return 0;
    }

    msa::option* option = new msa::option(vm);
    tbb::task_scheduler_init init(option->cpuNum);
    msa::utility* util = new msa::utility;
    Params* param = new Params(vm);

    Tree* T, * newT;
    partitionInfo_t* P;


    if (option->alnMode == 0) { // Twilight
        // Partition tree into subtrees
        T = readNewick(option->treeFile);
        P = new partitionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
        partitionTree(T->root, P);
        newT = reconsturctTree(T->root, P->partitionsRoot);
        if (P->partitionsRoot.size() > 1) {
            std::cout << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
            if (vm.count("output-subtrees")) outputSubtreeTrees(T, P, util, option);
        }
        // Start alignmnet on subtrees
        int proceeded = 0;
        auto alnSubtreeStart = std::chrono::high_resolution_clock::now();
        for (auto subRoot: P->partitionsRoot) {
            bool redo = true;
            auto subtreeStart = std::chrono::high_resolution_clock::now();
            ++proceeded;
            int subtree = T->allNodes[subRoot.first]->grpID;
            if (P->partitionsRoot.size() > 1) std::cout << "Start processing subtree No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
            while (redo) {
                Tree* subT = new Tree(subRoot.second.first);
                readSequences(util, option, subT);
                // // Generate MAFFT's files
                // std::vector<std::pair<std::string, std::string>> seqs;
                // for (auto s: util->rawSeqs) {
                //     int sIdx = util->seqsIdx[s.first];
                //     seqs.push_back(std::make_pair(std::to_string(sIdx), s.second));
                //     T->allNodes[s.first]->identifier = std::to_string(sIdx);
                // }
                // outputSubtreeTrees(T, P, util, option);
                // outputSubtreeSeqs("./subtrees/rnasim.raw.mafft.fa", seqs);
                // // printTree(subT->root, -1);
                // exit(1);

                auto treeBuiltStart = std::chrono::high_resolution_clock::now();
                partitionInfo_t * subP = new partitionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
                partitionTree(subT->root, subP);
                auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
                if (subP->partitionsRoot.size() > 1) std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " sub-subtrees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
                // Progressive alignment on each sub-subtree
                Tree* newSubT;
                if (subP->partitionsRoot.size() > 1) newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
                msaOnSubtree(subT, util, option, subP, *param);
                redo = option->redo;
                if (redo) {
                    util->seqsFree();
                    util->clearAll();
                    delete subP;
                    delete subT;
                    if (subP->partitionsRoot.size() > 1) delete newSubT;
                    std::cout << "Switch to position-specific gap penalty and realign the subtree.\n";
                    option->redo = false;
                    continue;
                }
                // If -l option, merge sub-subtree with transitivity merger
                if (subP->partitionsRoot.size() > 1) {
                    // Align adjacent sub-subtrees to create overlap alignment
                    auto alnStart = std::chrono::high_resolution_clock::now();
                    alignSubtrees(subT, newSubT, util, option, *param);
                    auto alnEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                    std::cout << "Aligned adjacent sub-subtree in " <<  alnTime.count() / 1000000 << " ms\n";
                    auto mergeStart = std::chrono::high_resolution_clock::now();
                    mergeSubtrees (subT, newSubT, util, option, *param);
                    auto mergeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
                    int totalSeqs = subT->root->msaIdx.size();
                    std::cout << "Merge " << newSubT->allNodes.size() << " sub-subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
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
                if (P->partitionsRoot.size() > 1) {
                    auto storeStart = std::chrono::high_resolution_clock::now();
                    storeFreq(util, subT, subtree);
                    util->storeCIGAR();
                    util->seqsFree();
                    util->clearAll();
                    auto storeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                    std::cout << "Stored the subtree alignment in " << storeTime.count() / 1000000 << " ms.\n";
                }
                delete subP;
                delete subT;
                if (subP->partitionsRoot.size() > 1) delete newSubT;
            }
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            if (P->partitionsRoot.size() > 1) std::cout << "Finished the alignment on subtree No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
            else                              std::cout << "Finished the alignment in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        if (P->partitionsRoot.size() > 1) {
            auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
            std::cout << "Finsihed alignment on all subtrees in " << alnSubtreeTime.count() / 1000000000 << " s.\n";
        }
        if (P->partitionsRoot.size() > 1) {
            util->nowProcess = 2; // merge subtrees
            updateSeqLen(T, P, util);
            if (option->merger == "transitivity") {
                auto alnStart = std::chrono::high_resolution_clock::now();
                alignSubtrees (T, newT, util, option, *param);
                mergeSubtrees (T, newT, util, option, *param);
                auto alnEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                std::cout << "Merge " << newT->allNodes.size() << " subtrees in " << alnTime.count() / 1000000 << " ms\n";
            }
            else if (option->merger == "profile") {
                auto alnStart = std::chrono::high_resolution_clock::now();
                mergeSubtrees (T, newT, util, option, *param);
                auto alnEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                std::cout << "Profile alignment on " << newT->allNodes.size() << " subtrees in " << alnTime.count() / 1000000 << " ms\n";
            }
            int totalSeqs = 0;
            auto outStart = std::chrono::high_resolution_clock::now();
            outputFinal (T, P, util, option, totalSeqs);
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::cout << "Output " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
        }
    }
    else { // Twilight-Mer
        readFrequency(util, option);
        T = new Tree(util->seqsLen, util->seqsIdx);
        util->nowProcess = 2; // merge subtrees
        auto alnStart = std::chrono::high_resolution_clock::now();
        mergeSubtrees (T, T, util, option, *param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Profile alignment on " << T->allNodes.size() << " subtrees in " << alnTime.count() / 1000000 << " ms\n";
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        outputFinal (T, P, util, option, totalSeqs);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output " << T->allNodes.size() - 1 << " subtrees (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
    }
    
    
    
    // output MSA
    if (vm.count("output")) {
        auto outStart = std::chrono::high_resolution_clock::now();
        outputAln(util, option, T);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
        // Calculate sum-of-pairs score
        if (vm.count("sum-of-pairs-score")) {
            auto spStart = std::chrono::high_resolution_clock::now();
            double score = calSPScore(option->outFile, util, param);
            auto spEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds spTime = spEnd - spStart;
            std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
        }
    }
    delete T;
    
    if (option->alnMode == 0) {
        delete newT;
        delete P;
    }
    
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
