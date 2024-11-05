#ifndef UTIL_HPP
#include "util.hpp"
#endif

#ifndef PROCESS_HPP
#include "msa-cpu.hpp"
#endif

po::options_description mainDesc("MSA Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>(), "Initial Tree - Newick format (required if building MSA from raw sequences)")
        ("sequences,i", po::value<std::string>(), "Input tip sequences - Fasta format")
        ("files, f", po::value<std::string>(), "Path to the directory containing all MSA files. MSA files should be in Fasta format")
        ("cpu-num,c",  po::value<int>(), "Number of CPU cores")
        ("max-leaves,l",  po::value<int>(), "Maximum number of leaves per sub-subtree, used for transitivity merger")
        ("max-subtree,m", po::value<int>(), "Maximum number of leaves per subtree")
        ("temp-dir,d", po::value<std::string>(), "Directory for storing temporary files")
        ("output,o", po::value<std::string>(), "Output file name")
        ("gappy-vertical,v", po::value<float>()->default_value(0.9), "If the proportion of gaps in a column exceeds this value, the column will be defined as a gappy column. Set to 1 to disable this feature.")
        ("gappy-horizon,z", po::value<float>()->default_value(1), "Minimum number of consecutive gappy columns, which will be removed during alignment.")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment, have to be used with -o option")
        ("wildcard,w", "Treat unknown 'N' as a wildcard and aligned to usual letters")
        ("print-detail", "Print out every detail")
        ("match",      po::value<paramType>()->default_value(18), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(-8), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-50), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-5), "Gap extend penalty")
        ("trans",      po::value<paramType>()->default_value(2), "The ratio of transitions to transversions")
        ("xdrop",      po::value<paramType>()->default_value(600), "X-drop value")
        ("scoring-matrix", po::value<int>()->default_value(0), "0: default, 1: kiruma, 2: user defined")
        ("pam-n",      po::value<int>()->default_value(200), "'n' in the PAM-n matrix. Should be used with kimura scoring matrix")
        ("output-type", po::value<std::string>()->default_value("FASTA"), "FASTA or CIGAR, CIGAR stands for CIGAR-like compressed format")
        ("merge-subtrees", po::value<std::string>()->default_value("p"), "t: transitivity merger, p: progressive alignment")
        ("delete-temp", "delete temporary subtree folder and files.")
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

    msa::option* option = new msa::option(vm);
    tbb::task_scheduler_init init(option->cpuNum);
    // setOptions(vm, option);
    msa::utility* util = new msa::utility;

    Params* param = new Params(vm);

    // Partition tree into subtrees
    Tree* T = readNewick(option->treeFile);
    std::cout << "Total leaves: " << T->m_numLeaves << '\n';
    // printTree(T->root, -1);
    partitionInfo_t* P = new partitionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);
    // readSequences(vm, util, option, T);
    if (P->partitionsRoot.size() > 1) std::cout << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
    
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
        if (P->partitionsRoot.size() > 1) outputSubtree(subT, option, subtree);
        util->setSubtreeIdx(subtree);
        readSequences(util, option, subT);
        // std::cout << "Subtree No." << subtree << " contains "<< subT->m_numLeaves << " sequences.\n";

        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        partitionInfo_t * subP = new partitionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(subT->root, subP);
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        if (subP->partitionsRoot.size() > 1) std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " sub-subtrees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        auto msaStart = std::chrono::high_resolution_clock::now();
        Tree* newSubT;
        if (subP->partitionsRoot.size() > 1) newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        msaOnSubtree(subT, util, option, subP, *param);
        auto msaEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        // std::cout << "MSA on sub-subtree " << subtree << " in " <<  msaTime.count() / 1000000000 << " s\n";
        
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
        // for (auto node: subT->allNodes) delete node.second;
        if (P->partitionsRoot.size() > 1) {
            // std::string tempDir = option->tempDir; 
            // std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            // std::string subtreeAlnFile = tempDir + '/' + subtreeFileName + ".temp.aln";
            // std::string subtreeFreqFile = tempDir + '/' + subtreeFileName + ".freq.txt";
            // outputFreq(subtreeFreqFile, util, subT, subtree);
            storeFreq(util, subT, subtree);
            util->storeCIGAR();
            util->seqsFree();
            util->clearAll();
            // outputAln(subtreeAlnFile, util, option, subT, subtree);
            
        }
        delete subT;
        delete subP;
        if (subP->partitionsRoot.size() > 1) delete newSubT;
        auto subtreeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
        std::cout << "Finished the alignment on subtree No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
    }
    if (P->partitionsRoot.size() > 1) {
        auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
        std::cout << "Finsihed alignment on all subtrees in " << alnSubtreeTime.count() / 1000000 << " ms.\n";
    }
    
    if (P->partitionsRoot.size() > 1) {
        util->nowProcess = 2; // merge subtrees
        // readFreq(option->tempDir, T, P, util);
        updateSeqLen(T, P, util);
        if (option->merger == "transitivity") {
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees(T, newT, util, option, *param);
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
        outputFinal (T, P, util, option, totalSeqs);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << outTime.count() / 1000000 << " ms\n";
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
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
