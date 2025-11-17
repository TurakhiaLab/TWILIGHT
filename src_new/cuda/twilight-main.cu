#ifndef MSA_HPP
#include "../msa.hpp"
#endif

#ifndef PHYLO_HPP
#include "../phylogeny.hpp"
#endif

#include "../version.hpp"

#include <tbb/global_control.h>
#include <boost/filesystem.hpp>
#include <chrono>

#define DEFAULT_ALN 0
#define MERGE_MSA 1
#define PLACE_WO_TREE 2
#define PLACE_W_TREE 3

po::options_description mainDesc("TWILIGHT Command Line Arguments", 120);

void parseArguments(int argc, char** argv)
{
    
    // Section: I/O
    po::options_description inputDesc("Inputs\n [1] Build MSA From Unaligned Sequences\n [2] Merge Multiple MSAs\n [3] Add New Sequences to Existing MSA");
    inputDesc.add_options()
        ("tree,t", po::value<std::string>(), "Guide tree (Newick format): required for [1]; optional for [3].")
        ("sequences,i", po::value<std::string>(), "Unaligned sequences file (FASTA format): required for [1] and [3].")
        ("alignment,a", po::value<std::string>(), "Backbone alignments (FASTA format): required for [3].")
        ("files,f", po::value<std::string>(), "Directory containing all MSA files. MSA files (FASTA format): required for [2].");
   
    po::options_description outputDesc("Output/File Options");
    outputDesc.add_options()
        ("output,o", po::value<std::string>(), "Output file name (required).")
        ("temp-dir,d", po::value<std::string>(), "Directory for storing temporary files.")
        ("keep-temp,k", "Keep the temporary directory.")
        ("compress,c", "Write output files in compressed (.gz) format")
        ("overwrite", "Force overwriting the output file.")
        ("write-prune", "Write the pruned tree to the output directory.");

    // Section: Hardware
    po::options_description hardwareDesc("Hardware Usage");
    hardwareDesc.add_options()
        ("cpu,C", po::value<int>(), "Number of CPU cores. Default: all available cores.")
        ("gpu,G", po::value<int>(), "Number of GPUs. Default: all available GPUs.")
        ("gpu-index", po::value<std::string>(), "Specify the GPU to be used, separated by commas. Ex. 0,2,3")
        ("cpu-only", "Run the program only on CPU.");

    // Section: Alignment Parameters
    po::options_description alignDesc("Alignment Options and Parameters");
    alignDesc.add_options()
        ("type", po::value<std::string>(), "Data type. n: nucleotide, p: protein. Will be automatically inferred if not provided.")
        ("max-subtree,m", po::value<int>(), "Maximum number of leaves in a subtree.")
        ("remove-gappy,r", po::value<float>()->default_value(0.95), "Threshold for removing gappy columns. Set to 1 to disable this feature.")
        ("prune", "Prune the input guide tree based on the presence of unaligned sequences.")
        ("wildcard,w", "Treat unknown or ambiguous bases as wildcards and align them to usual letters.")
        ("no-align-gappy", "Do not align gappy columns. This will create a longer MSA (larger file).")
        // ("psgop", po::value<std::string>()->default_value("y"), "y: Enable, n: Disable position-specific gap open penalty.")
        ("length-deviation", po::value<float>(), "Sequences whose lengths deviate from the average by more than the specified fraction will be deferred or excluded.")
        ("max-ambig", po::value<float>()->default_value(0.1), "Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded.")
        ("filter", "Exclude sequences with high ambiguity or length deviation.")
        ("rooted", "Keep the original tree root (disable automatic re-rooting for parallelism)");
        
    // Section: Scoring Parameters
    po::options_description scoringDesc("Scoring Parameters");
    scoringDesc.add_options()
        ("match", po::value<float>()->default_value(18), "Match score.")
        ("mismatch", po::value<float>()->default_value(-8), "Mismatch penalty for transversions.")
        ("transition", po::value<float>()->default_value(-4), "Score for transitions.")
        ("gap-open", po::value<float>()->default_value(-50), "Gap-Open penalty.")
        ("gap-extend", po::value<float>()->default_value(-5), "Gap-Extend penalty.")
        ("gap-ends", po::value<float>(), "Gap penalty at ends, default set to the same as the gap extension penalty.")
        ("xdrop", po::value<float>()->default_value(600), "X-drop value (scale). The actual X-drop will be multiplied by the gap-extend penalty.")
        ("matrix,x", po::value<std::string>(), "Use a user-defined substitution matrix.");        

    po::options_description generalDesc("General");
    generalDesc.add_options()
        ("check", "Check the final alignment. Sequences with no legal alignment will be displayed.")
        ("verbose,v", "Print out every detail process.")
        ("help,h", "Print help messages.")
        ("version", "Show program version.");

    // Setup boost::program_options
    mainDesc.add(inputDesc).add(outputDesc).add(hardwareDesc).add(alignDesc).add(scoringDesc).add(generalDesc);
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
        std::cerr << "Error: " << e.what() << std::endl;
        // std::cerr << mainDesc << std::endl;
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

    msa::Option* option = new msa::Option(vm);
    msa::SequenceDB* database = new msa::SequenceDB();
    tbb::global_control init(tbb::global_control::max_allowed_parallelism, option->cpuNum);
    option->getGpuInfo(vm);
    msa::Params* param = new msa::Params(vm, option->type);

    if (option->alnMode == DEFAULT_ALN) { // Twilight
        phylogeny::Tree* T = new phylogeny::Tree(option->treeFile, option->reroot);
        if (vm.count("prune")) {
            std::unordered_set<std::string> seqNames;
            msa::io::readSequenceNames(option->seqFile, seqNames);
            pruneTree(T, seqNames);
            if (vm.count("write-prune")) msa::io::writePrunedTree(T, option);
        }
        phylogeny::PartitionInfo * P = new phylogeny::PartitionInfo(option->maxSubtree, 0, 0); 
        P->partitionTree(T->root);
        phylogeny::Tree* subRoot_T = phylogeny::constructTreeFromPartitions(T->root, P);
        if (P->partitionsRoot.size() > 1) {
            std::cerr << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
            msa::io::writeSubtrees(T, P, option);
        }
        // Start alignment on subtrees
        int proceeded = 0;
        auto alnSubtreeStart = std::chrono::high_resolution_clock::now();
        for (auto subRoot: P->partitionsRoot) {
            auto subtreeStart = std::chrono::high_resolution_clock::now();
            ++proceeded;
            int subtree = T->allNodes[subRoot.first]->grpID;
            if (P->partitionsRoot.size() > 1) std::cerr << "Start processing subalignment No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
            // Read sequences for each subtree
            phylogeny::Tree* subT = new phylogeny::Tree(subRoot.second.first);
            msa::io::readSequences(option->seqFile, database, option, subT);
            // Progressive alignment on each subtree
            msa::progressive::msaOnSubtree(subT, database, option, *param, msa::progressive::gpu::alignmentKernel_GPU);
            // post-alignment debugging
            if (option->debug) database->debug();
            if (P->partitionsRoot.size() > 1) {
                // Store subtree profile
                auto storeStart = std::chrono::high_resolution_clock::now();
                database->storeSubtreeProfile(subT, option->type, subtree);
                msa::io::writeSubAlignments(database, option, subtree, subT->root->getAlnLen(database->currentTask));
                phylogeny::updateSubrootInfo(subRoot_T->allNodes[subT->root->identifier], subT, subtree);
                database->cleanSubtreeDB();
                auto storeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                std::cerr << "Stored the subalignments in " << storeTime.count() / 1000000 << " ms.\n";
            }
            else {
                // Output final alignment
                auto outStart = std::chrono::high_resolution_clock::now();
                msa::io::writeWholeAlignment(database, option, subT->root->getAlnLen(database->currentTask));
                auto outEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds outTime = outEnd - outStart;
                std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
                std::cerr << "Wrote alignment to " << outFileName << " in " <<  outTime.count() / 1000000 << " ms\n";
            }
            delete subT;
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            if (P->partitionsRoot.size() > 1) std::cerr << "Finished subalignment No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
            else                              std::cerr << "Finished the alignment in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        if (P->partitionsRoot.size() > 1) {
            auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
            std::cerr << "Finished all subalignments in " << alnSubtreeTime.count() / 1000000000 << " s.\n";
            // Merge subalignment
            database->currentTask = 2;
            subRoot_T->showTree();
            msa::progressive::msaOnSubtree(subRoot_T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
            int totalSeqs = 0;
            auto outStart = std::chrono::high_resolution_clock::now();
            msa::io::writeFinalAlignments(database, option, totalSeqs);
            msa::io::writeWholeAlignment(database, option, subRoot_T->root->getAlnLen(database->currentTask));
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
            std::cerr << "Wrote " << subRoot_T->allNodes.size() << " subalignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        }
        delete T;
        delete subRoot_T;
        delete P;
    }
    else if (option->alnMode == MERGE_MSA) { // Twilight Merging alignments
        phylogeny::Tree* T = msa::io::readAlignments_and_buildTree(database, option);
        database->currentTask = 2;
        T->showTree();
        msa::progressive::msaOnSubtree(T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        msa::io::writeFinalAlignments(database, option, totalSeqs);
        msa::io::writeWholeAlignment(database, option, T->root->getAlnLen(database->currentTask));
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
        std::cerr << "Wrote " << T->allNodes.size() << " Alignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        delete T;
    }
    else if (option->alnMode == PLACE_WO_TREE) {
        database->currentTask = 2;
        std::unordered_set<std::string> seqNames;
        msa::io::readSequenceNames(option->seqFile, seqNames);
        phylogeny::Tree* T = new phylogeny::Tree(seqNames);
        msa::io::readSequences(option->seqFile, database, option, T);
        msa::io::readBackboneAlignment(T, database, option);
        // T->showTree();
        msa::progressive::msaOnSubtree(T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
        if (option->debug) database->debug();
        msa::io::update_and_writeAlignment(database, option, option->backboneAlnFile, -1);
        boost::filesystem::path seqPath(option->seqFile);
        std::string placedSeqFile = option->tempDir + "/" + seqPath.stem().string() + ".final.aln";
        msa::io::writeAlignment(placedSeqFile, database, T->root->getAlnLen(database->currentTask), option->compressed);
        msa::io::writeWholeAlignment(database, option, T->root->getAlnLen(database->currentTask));
        delete T;
    }
    else if (option->alnMode == PLACE_W_TREE) {
        phylogeny::Tree* T = new phylogeny::Tree(option->treeFile, option->reroot);
        phylogeny::PartitionInfo * P = new phylogeny::PartitionInfo(option->maxSubtree, 0, 0); 
        P->partitionTree(T->root);
        phylogeny::Tree* subRoot_T = phylogeny::constructTreeFromPartitions(T->root, P);
        if (P->partitionsRoot.size() > 1) {
            std::cerr << "Decomposed the tree into " << P->partitionsRoot.size() << " subtrees.\n";
            msa::io::writeSubtrees(T, P, option);
        }
        // Start alignment on subtrees
        int proceeded = 0;
        auto alnSubtreeStart = std::chrono::high_resolution_clock::now();
        for (auto subRoot: P->partitionsRoot) {
            auto subtreeStart = std::chrono::high_resolution_clock::now();
            ++proceeded;
            int subtree = T->allNodes[subRoot.first]->grpID;
            if (P->partitionsRoot.size() > 1) std::cerr << "Start processing subalignment No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
            // Read backbone alignment and sequences for each subtree
            phylogeny::Tree* subT = new phylogeny::Tree(subRoot.second.first);
            msa::io::readSequences(option->backboneAlnFile, database, option, subT);
            msa::io::readSequences(option->seqFile, database, option, subT);
            phylogeny::Tree* placementT = database->getPlacementTree(subT);
            // Progressive alignment on each subtree
            msa::progressive::msaOnSubtree(placementT, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
            subT->extractResult(placementT);
            delete placementT;
            // post-alignment debugging
            if (option->debug) database->debug();
            if (P->partitionsRoot.size() > 1) {
                // Store subtree profile
                auto storeStart = std::chrono::high_resolution_clock::now();
                database->storeSubtreeProfile(subT, option->type, subtree);
                msa::io::writeSubAlignments(database, option, subtree, subT->root->getAlnLen(database->currentTask));
                phylogeny::updateSubrootInfo(subRoot_T->allNodes[subT->root->identifier], subT, subtree);
                database->cleanSubtreeDB();
                auto storeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                std::cerr << "Stored the subalignments in " << storeTime.count() / 1000000 << " ms.\n";
            }
            else {
                // Output final alignment
                auto outStart = std::chrono::high_resolution_clock::now();
                msa::io::writeWholeAlignment(database, option, subT->root->getAlnLen(database->currentTask));
                auto outEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds outTime = outEnd - outStart;
                std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
                std::cerr << "Wrote alignment to " << outFileName << " in " <<  outTime.count() / 1000000 << " ms\n";
            }
            delete subT;
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            if (P->partitionsRoot.size() > 1) std::cerr << "Finished subalignment No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
            else                              std::cerr << "Finished the alignment in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        if (P->partitionsRoot.size() > 1) {
            auto alnSubtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnSubtreeTime = alnSubtreeEnd - alnSubtreeStart;
            std::cerr << "Finished all subalignments in " << alnSubtreeTime.count() / 1000000000 << " s.\n";
            // Merge subalignment
            database->currentTask = 2;
            subRoot_T->showTree();
            msa::progressive::msaOnSubtree(subRoot_T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
            int totalSeqs = 0;
            auto outStart = std::chrono::high_resolution_clock::now();
            msa::io::writeFinalAlignments(database, option, totalSeqs);
            msa::io::writeWholeAlignment(database, option, subRoot_T->root->getAlnLen(database->currentTask));
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
            std::cerr << "Wrote " << subRoot_T->allNodes.size() << " subalignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        }
        delete T;
        delete subRoot_T;
        delete P;
    }
    
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cerr << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
