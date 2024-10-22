#ifndef UTIL_HPP
#include "util.hpp"
#endif

#ifndef MSAGPU_HPP
#include "msa-gpu.cuh"
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
        ("merge-subtrees", po::value<std::string>()->default_value("p"), "t: transitivity merger, p: progressive alignment")
        ("gpu-index", po::value<std::string>(), "Specify the GPU index, separated by commas. Ex. 0,2,3")
        ("gappy-vertical", po::value<float>()->default_value(0.9), "If the proportion of gaps in a column exceeds this value, the column will be defined as a gappy column. Set to 1 to disable this feature.")
        ("gappy-horizon", po::value<float>()->default_value(1), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment, have to be used with -o option")
        ("print-detail", "Print out every detail")
        ("keep-temp", "Keep temporary subtree files. If not specified, temporary files will be deleted to save storage space.")
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

    msa::option* option = new msa::option(vm);
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
    option->gpuNum = gpuNum;
    option->gpuIdx = gpuIdx;
    // setOptions(vm, option);
    msa::utility* util = new msa::utility;

    Params* param = new Params(vm);

    // Partition tree into subtrees
    Tree* T = readNewick(vm);
    std::cout << "Total leaves: " << T->m_numLeaves << '\n';
    // printTree(T->root, -1);
    partitionInfo_t* P = new partitionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);
    if (P->partitionsRoot.size() > 1) std::cout << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
    // if (P->partitionsRoot.size() > 1) readSequencesNoutputTemp(vm, T, P, util, option);
    // exit(1);
    // int totalL = 0;
    // for (auto root: P->partitionsRoot) {
    //     Tree* subT = new Tree(root.second.first);
    //     totalL += subT->m_numLeaves;
    //     std::cout << subT->m_numLeaves << '\t' << totalL << '\n';
    // }

    int proceeded = 0;
    auto alnSubtreeStart = std::chrono::high_resolution_clock::now();
    for (auto subRoot: P->partitionsRoot) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        ++proceeded;
        int subtree = T->allNodes[subRoot.first]->grpID;
        if (P->partitionsRoot.size() > 1) std::cout << "Start processing subtree No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
        Tree* subT = new Tree(subRoot.second.first);
        if (P->partitionsRoot.size() > 1 && !option->deleteTemp) outputSubtree(subT, option, subtree);
        readSequences(vm, util, option, subT);
        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        partitionInfo_t * subP = new partitionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(subT->root, subP);
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        if (subP->partitionsRoot.size() > 1) std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " sub-subtrees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        // auto msaStart = std::chrono::high_resolution_clock::now();
        Tree* newSubT;
        if (subP->partitionsRoot.size() > 1) newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        msaOnSubtreeGpu(subT, util, option, subP, *param);
        // auto msaEnd = std::chrono::high_resolution_clock::now();
        // std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        // std::cout << "MSA on sub-subtree " << subtree << " in " <<  msaTime.count() / 1000000000 << " s\n";
        
        if (subP->partitionsRoot.size() > 1) {
            // Align adjacent sub-subtrees to create overlap alignment
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtreesGpu(subT, newSubT, util, option, *param);
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
        // for (auto node: subT->allNodes) delete node.second;
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
        // readFreq(option->tempDir, T, P, util);
        updateSeqLen(T, P, util);
        if (option->merger == "transitivity") {
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees (T, newT, util, option, *param);
            mergeSubtrees (T, newT, util, option, *param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            std::cout << "Merge " << newT->allNodes.size() << " subtrees in " << alnTime.count() / 1000000 << " ms\n";
            // auto outStart = std::chrono::high_resolution_clock::now();
            // outputFinal (vm, T, P, util, option, totalSeqs);
            // auto outEnd = std::chrono::high_resolution_clock::now();
            // std::chrono::nanoseconds outTime = outEnd - outStart;
            // std::cout << "Output " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
        }
        else if (option->merger == "progressive") {
            auto alnStart = std::chrono::high_resolution_clock::now();
            mergeSubtrees (T, newT, util, option, *param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            std::cout << "Progressive alignment on " << newT->allNodes.size() << " subtrees in " << alnTime.count() / 1000000 << " ms\n";
        }
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        outputFinal (vm, T, P, util, option, totalSeqs);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
    }
    
    
    // output MSA
    if (vm.count("output")) {
        std::string outFile = vm["output"].as<std::string>();
        if (outFile == "") outFile = "output.aln";
        auto outStart = std::chrono::high_resolution_clock::now();
        outputAln(outFile, util, option, T);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
        // Calculate sum-of-pairs score
        if (vm.count("sum-of-pairs-score")) {
            auto spStart = std::chrono::high_resolution_clock::now();
            double score = calSPScore(outFile, util, param);
            auto spEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds spTime = spEnd - spStart;
            std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
        }
    }
    delete T;
    delete P;
    delete newT;
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
