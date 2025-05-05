#ifndef UTIL_HPP
#include "util.hpp"
#endif

#ifndef MSAGPU_HPP
#include "msa-gpu.cuh"
#endif

#include "version.hpp"

po::options_description mainDesc("TWILIGHT Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>(), "Guide Tree - Newick format (required if building MSA from raw sequences).")
        ("sequences,i", po::value<std::string>(), "Input tip sequences - Fasta format (required if building MSA from raw sequences).")
        ("files,f", po::value<std::string>(), "Path to the directory containing all MSA files. MSA files - Fasta format.")
        ("output,o", po::value<std::string>(), "Output file name (required).")
        ("cpu,C",  po::value<int>(), "Number of CPU cores. Default: all available cores.")
        ("gpu,G",  po::value<int>(), "Number of GPUs. Default: all available GPUs.")
        ("gpu-index", po::value<std::string>(), "Specify the GPU to be used, separated by commas. Ex. 0,2,3")
        ("max-subtree,m", po::value<int>(), "Maximum number of leaves in a subtree.")
        ("max-group,g",  po::value<int>(), "Maximum number of leaves in a group within a subtree, used for performing transitivity merger in a subtree.")
        ("temp-dir,d", po::value<std::string>(), "Directory for storing temporary files.")
        ("remove-gappy,r", po::value<float>()->default_value(0.95), "Threshold for removing gappy columns. Set to 1 to disable this feature.")
        // ("gappy-horizon,z", po::value<float>()->default_value(1), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("wildcard,w", "Treat unknown or ambiguous bases as wildcards and align them to usual letters.")
        ("verbose,v", "Print out every detail process.")
        ("psgop,p", po::value<std::string>()->default_value("y"), "y: Enable, n: Disable position-specific gap open penalty.")
        ("match",      po::value<float>()->default_value(18), "Match score.")
        ("mismatch",   po::value<float>()->default_value(-8), "Mismatch penalty for transversions.")
        ("transition", po::value<float>()->default_value(5), "Score for transitions")
        ("gap-open",   po::value<float>()->default_value(-50), "Gap-Open penalty")
        ("gap-extend", po::value<float>()->default_value(-5), "Gap-Extend penalty")
        ("gap-ends",   po::value<float>(), "Gap penalty at ends, default set to the same as the gap extension penalty.")
        ("xdrop",      po::value<float>()->default_value(600), "X-drop value (scale). The actual X-drop will be multiplied by the gap-extend penalty.")
        ("matrix,x", po::value<std::string>(), "Use a user-defined substitution matrix.")
        ("output-type", po::value<std::string>()->default_value("FASTA"), "FASTA or CIGAR, CIGAR stands for CIGAR-like compressed format.")
        ("merge", po::value<std::string>()->default_value("p"), "t: Transitivity merger, p: Progressive alignment.")
        ("keep-temp,k", "Keep the temporary directory.")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment.")
        ("no-align-gappy", "Do not align gappy columns. This will create a longer MSA (larger file).")
        ("length-deviation", po::value<float>(), "Sequences whose lengths deviate from the average by more than the specified fraction will be deferred or excluded.")
        ("max-ambig", po::value<float>()->default_value(0.1), "Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded.")
        ("filter", "Exclude sequences with high ambiguity or length deviation.")
        ("check", "Check the final alignment. Sequences with no legal alignment will be displayed.")
        ("cpu-only", "Run the program only on CPU.")
        ("help,h", "Print help messages")
        ("version", "Show program version");
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
    if(vm.count("help") || argc == 1) {
        std::cerr << mainDesc << std::endl;
        return 0;
    }
    if(vm.count("version")) {
        std::cerr << "TWILIGHT Version " << PROJECT_VERSION << std::endl;
        return 0;
    }

    msa::option* option = new msa::option(vm);
    tbb::task_scheduler_init init(option->cpuNum);
    getGpuInfo(vm, option);
    msa::utility* util = new msa::utility;
    Params* param = new Params(vm);

    Tree* T = nullptr, * newT = nullptr;
    partitionInfo_t* P = nullptr;


    if (option->alnMode == 0) { // Twilight
        // Partition tree into subtrees
        T = readNewick(option->treeFile);
        P = new partitionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
        partitionTree(T->root, P);
        newT = reconsturctTree(T->root, P->partitionsRoot);
        if (P->partitionsRoot.size() > 1) {
            std::cout << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
            outputSubtreeTrees(T, P, util, option);
        }
        // Start alignmnet on subtrees
        int proceeded = 0;
        auto alnSubtreeStart = std::chrono::high_resolution_clock::now();
        for (auto subRoot: P->partitionsRoot) {
            bool redo = true;
            auto subtreeStart = std::chrono::high_resolution_clock::now();
            ++proceeded;
            int subtree = T->allNodes[subRoot.first]->grpID;
            if (P->partitionsRoot.size() > 1) std::cout << "Start processing subalignment No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
            while (redo) {
                Tree* subT = new Tree(subRoot.second.first);
                readSequences(util, option, subT);
                auto treeBuiltStart = std::chrono::high_resolution_clock::now();
                partitionInfo_t * subP = new partitionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
                partitionTree(subT->root, subP);
                auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
                if (subP->partitionsRoot.size() > 1) std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " groups in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
                // Progressive alignment on each sub-subtree
                Tree* newSubT = nullptr;
                if (subP->partitionsRoot.size() > 1) newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
                msaOnSubtreeGpu(subT, util, option, subP, *param);
                redo = option->redo;
                if (redo) {
                    util->seqsFree();
                    util->clearAll();
                    delete subP;
                    delete subT;
                    if (subP->partitionsRoot.size() > 1) delete newSubT;
                    std::cout << "Switch to position-specific gap penalty and realign.\n";
                    option->redo = false;
                    continue;
                }
                // If -g option, merge groups with transitivity merger
                if (subP->partitionsRoot.size() > 1) {
                    auto alnStart = std::chrono::high_resolution_clock::now();
                    alignSubtrees(subT, newSubT, util, option, *param);
                    auto alnEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                    std::cout << "Aligned adjacent groups in " <<  alnTime.count() / 1000000 << " ms\n";
                    auto mergeStart = std::chrono::high_resolution_clock::now();
                    mergeSubtrees (subT, newSubT, util, option, *param);
                    auto mergeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
                    int totalSeqs = subT->root->msaIdx.size();
                    std::cout << "Merge " << newSubT->allNodes.size() << " groups (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
                }
                for (auto sIdx: subT->root->msaIdx) T->allNodes[subT->root->identifier]->msaIdx.push_back(sIdx);
                // post-alignment debugging
                if (option->debug) {
                    auto dbgStart = std::chrono::high_resolution_clock::now();
                    int debugNum = 0;
                    util->debug(debugNum);
                    auto dbgEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
                    std::cout << "Completed checking " << debugNum << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
                }
                if (P->partitionsRoot.size() > 1) {
                    auto storeStart = std::chrono::high_resolution_clock::now();
                    storeFreq(util, subT, subtree);
                    outputSubAln(util, option, subT, subtree);
                    // util->storeCIGAR();
                    util->seqsFree();
                    util->clearAll();
                    auto storeEnd = std::chrono::high_resolution_clock::now();
                    std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                    std::cout << "Stored the subalignments in " << storeTime.count() / 1000000 << " ms.\n";
                }
                delete subP;
                delete subT;
                if (subP->partitionsRoot.size() > 1) delete newSubT;
            }
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            if (P->partitionsRoot.size() > 1) std::cout << "Finished subalignment No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
            else                              std::cout << "Finished the alignment in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        if (P->partitionsRoot.size() > 1) {
            auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
            std::cout << "Finsihed all subalignments in " << alnSubtreeTime.count() / 1000000000 << " s.\n";
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
                std::cout << "Merge " << newT->allNodes.size() << " subalignments in " << alnTime.count() / 1000000 << " ms\n";
            }
            else if (option->merger == "profile") {
                auto alnStart = std::chrono::high_resolution_clock::now();
                mergeSubtrees (T, newT, util, option, *param);
                auto alnEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds alnTime = alnEnd - alnStart;
                std::cout << "Progressive alignment on " << newT->allNodes.size() << " subalignments in " << alnTime.count() / 1000000 << " ms\n";
            }
            int totalSeqs = 0;
            auto outStart = std::chrono::high_resolution_clock::now();
            outputFinal (T, P, util, option, totalSeqs);
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::cout << "Output " << newT->allNodes.size() << " alignments (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
        }
    }
    else { // Twilight-Mer
        readFrequency(util, option);
        T = new Tree(util->seqsLen, util->seqsIdx);
        util->nowProcess = 2; // merge subtrees
        auto alnStart = std::chrono::high_resolution_clock::now();
        P = new partitionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(T->root, P);
        msaOnSubtree(T, util, option, P, *param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Progressive alignment on " << T->allNodes.size() << " MSAs in " << alnTime.count() / 1000000 << " ms\n";
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        outputFinal (T, P, util, option, totalSeqs);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output " << T->allNodes.size() - 1 << " MSAs (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
    }
    
    
    
    // output MSA
    auto outStart = std::chrono::high_resolution_clock::now();
    outputAln(util, option, T);
    auto outEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds outTime = outEnd - outStart;
    std::cout << "Output file to " << option->outFile << " in " <<  outTime.count() / 1000000 << " ms\n";
    // Calculate sum-of-pairs score
    if (vm.count("sum-of-pairs-score")) {
        auto spStart = std::chrono::high_resolution_clock::now();
        double score = calSPScore(option->outFile, util, param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
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
