#ifndef MSA_HPP
#include "../msa.hpp"
#endif

#include "../version.hpp"

// CLI11 command line parser (https://github.com/CLIUtils/CLI11)
// Licensed under the 3-Clause BSD License.
// Copyright (c) 2017-2024 University of Cincinnati, developed by Henry Schreiner.
#include "../CLI11.hpp"

#include <tbb/global_control.h>
#include <filesystem>
#include <chrono>
#include <hip/hip_runtime.h>

int main(int argc, char** argv) {

    auto mainStart = std::chrono::high_resolution_clock::now();
    
    CLI::App app{"TWILIGHT: Tree-guided Multiple Sequence Alignment", "twilight"};
    app.set_help_all_flag("--help-all", "Expand all help");
    app.set_version_flag("-V,--version", std::string("TWILIGHT Version ") + PROJECT_VERSION, "Show program version.");

    msa::CliOptions opt;

    // Section: I/O
    auto* inputDesc = app.add_option_group("Inputs",
        " [1] Build MSA From Unaligned Sequences\n [2] Merge Multiple MSAs\n [3] Add New Sequences to Existing MSA");
    inputDesc->add_option("-t,--tree", opt.treeFile, "Guide tree (Newick format): required for [1]; optional for [3].");
    inputDesc->add_option("-i,--sequences", opt.seqFile, "Unaligned sequences file (FASTA format): required for [1] and [3].");
    inputDesc->add_option("-a,--alignment", opt.backboneAlnFile, "Backbone alignments (FASTA format): required for [3].");
    inputDesc->add_option("-f,--files", opt.msaDir, "Directory containing all MSA files. MSA files (FASTA format): required for [2].");

    auto* outputDesc = app.add_option_group("Output/File Options");
    outputDesc->add_option("-o,--output", opt.outFile, "Output file name (required).");
    outputDesc->add_option("-d,--temp-dir", opt.tempDir, "Directory for storing temporary files.");
    outputDesc->add_flag("-k,--keep-temp", opt.keepTemp, "Keep the temporary directory.");
    outputDesc->add_flag("-c,--compress", opt.compress, "Write output files in compressed (.gz) format");
    outputDesc->add_flag("--overwrite", opt.overwrite, "Force overwriting the output file.");
    auto* opt_score = outputDesc->add_option("--score", opt.score, "Compute alignment score: 1: Consistency objective score (default), 2: WSP (Weighted Profile Sum-of-Pairs).")->expected(0, 1)->default_val(1);
    outputDesc->add_option("--score-file", opt.scoreFile, "Output file path for score. Default: <output-file>.score.");

    // Section: Hardware
    auto* hardwareDesc = app.add_option_group("Hardware Usage");
    auto* opt_cpu = hardwareDesc->add_option("-C,--cpu", opt.cpuNum, "Number of CPU cores. Default: all available cores.");
    auto* opt_gpu = hardwareDesc->add_option("-G,--gpu", opt.gpuNum, "Number of GPUs. Default: all available GPUs.");
    hardwareDesc->add_option("--gpu-index", opt.gpuIndex, "Specify the GPU to be used, separated by commas. Ex. 0,2,3");
    hardwareDesc->add_flag("--cpu-only", opt.cpuOnly, "Run the program only on CPU.");

    // Section: Alignment Parameters
    auto* alignDesc = app.add_option_group("Alignment Options and Parameters");
    alignDesc->add_option("--type", opt.typeStr, "Data type. n: nucleotide, p: protein. Will be automatically inferred if not provided.");
    auto* opt_max_subtree = alignDesc->add_option("-m,--max-subtree", opt.maxSubtree, "Maximum number of leaves in a subtree.");
    alignDesc->add_option("-r,--remove-gappy", opt.removeGappy, "Threshold for removing gappy columns. Set to 1 to disable this feature.")->default_val(0.95f);
    alignDesc->add_flag("-w,--wildcard", opt.wildcard, "Treat unknown or ambiguous bases as wildcards and align them to usual letters.");
    alignDesc->add_flag("--rooted", opt.rooted, "Keep the original tree root (disable automatic re-rooting for parallelism)");
    alignDesc->add_flag("--prune", opt.prune, "Prune the input guide tree based on the presence of unaligned sequences.");
    alignDesc->add_flag("--write-prune", opt.writePrune, "Write the pruned tree to the output directory.");

    auto* seqFilterDesc = app.add_option_group("Sequence Filtering Options");
    auto* opt_len_dev = seqFilterDesc->add_option("--length-deviation", opt.lenDev, "Sequences whose lengths deviate from the median by more than the specified fraction will be deferred or excluded.");
    seqFilterDesc->add_option("--max-ambig", opt.maxAmbig, "Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded.")->default_val(0.1f);
    auto* opt_max_len = seqFilterDesc->add_option("--max-len", opt.maxLen, "Sequences longer than max-len will be deferred or excluded.");
    auto* opt_min_len = seqFilterDesc->add_option("--min-len", opt.minLen, "Sequences shorter than min-len will be deferred or excluded.");
    seqFilterDesc->add_flag("--filter", opt.filter, "Exclude sequences with high ambiguity or length deviation.");
    seqFilterDesc->add_flag("--write-filtered", opt.writeFiltered, "Write the filtered sequences in FASTA format to the output directory.");

    // Section: Scoring Parameters
    auto* scoringDesc = app.add_option_group("Scoring Parameters");
    scoringDesc->add_option("--match", opt.match, "Match score.")->default_val(18.0f);
    scoringDesc->add_option("--mismatch", opt.mismatch, "Mismatch penalty for transversions.")->default_val(-8.0f);
    scoringDesc->add_option("--transition", opt.transition, "Score for transitions.")->default_val(-4.0f);
    auto* opt_gap_open = scoringDesc->add_option("--gap-open", opt.gapOpen, "Gap-Open penalty in normalized units (relative to average match = 1.0).");
    auto* opt_gap_extend = scoringDesc->add_option("--gap-extend", opt.gapExtend, "Gap-Extend penalty in normalized units (relative to average match = 1.0).");
    auto* opt_gap_ends = scoringDesc->add_option("--gap-ends", opt.gapEnds, "Gap penalty at ends, default set to the same as the gap extension penalty.");
    auto* opt_scale = scoringDesc->add_option("--scale", opt.scale, "Matrix and gap normalization scale constant (default: 600.0 for protein, 20.0 for nucleotide).");
    scoringDesc->add_option("--xdrop", opt.xdrop, "X-drop value (scale). The actual X-drop will be multiplied by the gap-extend penalty.")->default_val(600.0f);
    scoringDesc->add_option("-x,--matrix", opt.matrixFile, "Use a user-defined substitution matrix (only for nucleotide).");
    scoringDesc->add_option("-b,--blosum", opt.blosum, "Substitution matrix to use for protein sequences: 43 (PFASUM43), 45, 62, 80, or 99 (MIQS).")->default_val(62);

    auto* accurateDesc = app.add_option_group("Accuracy & Sampling Options");
    auto* opt_acc_level = accurateDesc->add_option("-A,--accuracy-level", opt.accuracyLevelStr, "Accuracy level: 1: traditional sum-of-pairs scoring (default), 2: + pairwise consistency, 3: + transitive consistency.");
    auto* opt_sample_mode = accurateDesc->add_option("-S,--sample-mode", opt.sampleMode, "Sampling preset mode for consistency: 0: auto (default), 1: all-to-all (slow), 2: high sample rate (med speed), 3: low sample rate (fast).");
    auto* opt_sample_rate = accurateDesc->add_option("--sample-rate", opt.sampleRate, "Sample rate for hierarchical consistency representative selection (must be >0, <1).");
    auto* opt_sample_top = accurateDesc->add_option("--sample-top-centrality", opt.sampleTopCentrality, "Fraction of representative quota to select by highest centrality before column coverage (must be >0, <=1).");
    auto* opt_sample_max_group = accurateDesc->add_option("--sample-max-group", opt.sampleMaxGroup, "Maximum number of sequences in a consistency group. Controls all-to-all pairwise complexity.");
    accurateDesc->add_option("--consistency-weight", opt.consistencyWeight, "Weight multiplier for consistency bonus in accurate mode (default: 5.0).")->default_val(5.0f);
    accurateDesc->add_flag("--no-refine-subtree", opt.noRefineSubtree, "Disable Layer 1 guide subtree refinement in accurate mode (enabled by default).");
    accurateDesc->add_option("--refine-method", opt.refineMethod, "Subtree refinement algorithm: \"upgma\" (default) or \"nj\".")->default_val("upgma");

    // Hidden legacy option
    app.add_flag("--refine-subtree", "Legacy option to enable Layer 1 guide subtree refinement");

    auto* generalDesc = app.add_option_group("General");
    generalDesc->add_flag("--check", opt.check, "Check the final alignment. Sequences with no legal alignment will be displayed.");
    generalDesc->add_flag("-v,--verbose", opt.verbose, "Print out every detail process.");

    if (argc == 1) {
        std::cerr << app.help() << std::endl;
        return 0;
    }

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    opt.hasScore = (opt_score->count() > 0);
    opt.hasCpu = (opt_cpu->count() > 0);
    opt.hasGpu = (opt_gpu->count() > 0);
    opt.hasMaxSubtree = (opt_max_subtree->count() > 0);
    opt.hasLenDev = (opt_len_dev->count() > 0);
    opt.hasMaxLen = (opt_max_len->count() > 0);
    opt.hasMinLen = (opt_min_len->count() > 0);
    opt.hasGapOpen = (opt_gap_open->count() > 0);
    opt.hasGapExtend = (opt_gap_extend->count() > 0);
    opt.hasGapEnds = (opt_gap_ends->count() > 0);
    opt.hasScale = (opt_scale->count() > 0);
    opt.hasAccuracyLevel = (opt_acc_level->count() > 0);
    opt.hasSampleMode = (opt_sample_mode->count() > 0);
    opt.hasSampleRate = (opt_sample_rate->count() > 0);
    opt.hasSampleTopCentrality = (opt_sample_top->count() > 0);
    opt.hasSampleMaxGroup = (opt_sample_max_group->count() > 0);

    msa::Option* option = new msa::Option(opt);
    msa::SequenceDB* database = new msa::SequenceDB();
    tbb::global_control init(tbb::global_control::max_allowed_parallelism, option->cpuNum);
    option->getGpuInfo(opt);
    msa::Params* param = new msa::Params(opt, option->type);

    if (option->alnMode == msa::DEFAULT_ALN || option->alnMode == msa::ACCURATE_ALN) { // Twilight Progressive Alignment
        phylogeny::Tree* T = new phylogeny::Tree(option->treeFile);
        if (opt.prune) {
            std::unordered_set<std::string> seqNames;
            msa::io::readSequenceNames(option->seqFile, seqNames);
            pruneTree(T, seqNames);
            if (opt.writePrune) msa::io::writePrunedTree(T, option);
        }
            msa::io::readSequenceNames(option->seqFile, seqNames);
            pruneTree(T, seqNames);
            if (vm.count("write-prune")) msa::io::writePrunedTree(T, option);
        }
        T->convert2binaryTree();
        T->calLeafNum();
        phylogeny::PartitionInfo * P = new phylogeny::PartitionInfo(option->maxSubtree, 0, 0); 
        P->partitionTree(T->root);
        phylogeny::Tree* subRoot_T = phylogeny::constructTreeFromPartitions_new(T->root, P);
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
            int subtree = (P->partitionsRoot.size() > 1) ? T->allNodes[subRoot.first]->grpID : -1;
            if (P->partitionsRoot.size() > 1) std::cerr << "Start processing subalignment No. " << subtree << ". (" << proceeded << '/' << P->partitionsRoot.size() << ")\n";
            // Read sequences for each subtree
            phylogeny::Tree* subT = new phylogeny::Tree(subRoot.second.first, option->reroot);
            msa::io::readSequences(option->seqFile, database, option, subT, subtree);
            // Progressive alignment on each subtree
            if (!option->accurate) {
                msa::progressive::msaOnSubtree(subT, database, option, *param, nullptr, subtree);
            } else {
                msa::accurate::msaOnSubtree_accurate(subT, database, option, *param, nullptr, subtree);
            }
            // post-alignment debugging
            if (option->debug) database->debug();
            if (P->partitionsRoot.size() > 1) {
                // Store subtree profile
                auto storeStart = std::chrono::high_resolution_clock::now();
                database->storeSubtreeProfile(subT, option->type, subtree);
                msa::io::writeSubAlignments(database, option, subtree, subT->root->getAlnLen(database->currentTask));
                std::string subNodeID = (subRoot_T && subRoot_T->allNodes.find(subRoot.first) != subRoot_T->allNodes.end()) ?
                                        subRoot.first : subT->root->identifier;
                phylogeny::updateSubrootInfo(subRoot_T->allNodes[subNodeID], subT, subtree);
                if (database->subtreeRepresentatives.count(subT->root->identifier)) {
                    database->subtreeRepresentatives[subNodeID] = database->subtreeRepresentatives[subT->root->identifier];
                }
                database->cleanSubtreeDB();
                auto storeEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds storeTime = storeEnd - storeStart;
                std::cerr << "Stored the subalignments in " << storeTime.count() / 1000000 << " ms.\n";
            }
            else {
                // Output final alignment
                auto outStart = std::chrono::high_resolution_clock::now();
                msa::io::writeFinalMSA(database, option, subT->root->getAlnLen(database->currentTask), param);
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
            if (!option->accurate) {
                msa::progressive::msaOnSubtree(subRoot_T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
            } else {
                msa::accurate::mergeSubtreesAccurate(subRoot_T, database, option, *param);
            }
            int totalSeqs = 0;
            auto outStart = std::chrono::high_resolution_clock::now();
            msa::io::update_and_writeAlignments(database, option, totalSeqs);
            msa::io::writeFinalMSA(database, option, subRoot_T->root->getAlnLen(database->currentTask), param);
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
            std::cerr << "Wrote " << subRoot_T->allNodes.size() << " subalignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        }
        delete T;
        delete subRoot_T;
        delete P;
    }
    else if (option->alnMode == msa::MERGE_MSA) { // Twilight Merging alignments
        phylogeny::Tree* T = msa::io::readAlignments_and_buildTree(database, option);
        database->currentTask = 2;
        // T->showTree();
        msa::progressive::msaOnSubtree(T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
        int totalSeqs = 0;
        auto outStart = std::chrono::high_resolution_clock::now();
        msa::io::update_and_writeAlignments(database, option, totalSeqs);
        msa::io::writeFinalMSA(database, option, T->root->getAlnLen(database->currentTask), param);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
        std::cerr << "Wrote " << T->allNodes.size() << " Alignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        delete T;
    }
    else if (option->alnMode == msa::PLACE_WO_TREE) {
        database->currentTask = 2;
        std::unordered_set<std::string> seqNames;
        msa::io::readSequenceNames(option->seqFile, seqNames);
        phylogeny::Tree* T = new phylogeny::Tree(seqNames);
        msa::io::readSequences(option->seqFile, database, option, T);
        msa::io::readBackboneAlignment(T, database, option);
        // T->showTree();
        if (!option->accurate) {
            msa::progressive::msaOnSubtree(T, database, option, *param, msa::progressive::cpu::alignmentKernel_CPU);
        } else {
            msa::accurate::placeWithoutTree_Accurate(T, database, option, *param);
        }
        if (option->debug) database->debug();
        auto outStart = std::chrono::high_resolution_clock::now();
        msa::io::update_and_writeAlignment(database, option, option->backboneAlnFile, -1);
        auto backboneEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds backboneTime = backboneEnd - outStart;
        std::cerr << "Updated and wrote backbone alignment in " << backboneTime.count() / 1000000 << " ms\n";
        std::filesystem::path seqPath(option->seqFile);
        std::string placedSeqFile = option->tempDir + "/" + seqPath.stem().string() + ".final.aln";
        msa::io::writeAlignment(placedSeqFile, database, T->root->getAlnLen(database->currentTask), option->compressed);
        msa::io::writeFinalMSA(database, option, T->root->getAlnLen(database->currentTask), param);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
        std::cerr << "Wrote placed sequences and backbone alignment (total " << T->root->alnNum << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        delete T;
    }
    else if (option->alnMode == msa::PLACE_W_TREE) {
        phylogeny::Tree* T = new phylogeny::Tree(option->treeFile);
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
            phylogeny::Tree* subT = new phylogeny::Tree(subRoot.second.first, false);
            msa::io::readSequences(option->backboneAlnFile, database, option, subT);
            msa::io::readSequences(option->seqFile, database, option, subT);
            if (option->reroot) subT->reroot(true);
            phylogeny::Tree* placementT = database->getPlacementTree(subT);
            // Progressive alignment on each subtree
            msa::progressive::msaOnSubtree(placementT, database, option, *param, msa::progressive::gpu::alignmentKernel_GPU);
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
                msa::io::writeFinalMSA(database, option, subT->root->getAlnLen(database->currentTask), param);
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
            msa::io::update_and_writeAlignments(database, option, totalSeqs);
            msa::io::writeFinalMSA(database, option, subRoot_T->root->getAlnLen(database->currentTask), param);
            auto outEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds outTime = outEnd - outStart;
            std::string outFileName = (option->compressed) ? (option->outFile + ".gz") : option->outFile;
            std::cerr << "Wrote " << subRoot_T->allNodes.size() << " subalignments (total " << totalSeqs << " sequences) to " << outFileName << " in " << outTime.count() / 1000000 << " ms\n";
        }
        delete T;
        delete subRoot_T;
        delete P;
    }

    delete option;
    delete database;
    delete param;
    
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cerr << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
