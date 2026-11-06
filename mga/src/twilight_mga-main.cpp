
#include "mga.hpp"
#include "option.hpp"
#include "phylogeny.hpp"
#include "block_manager.hpp"
#include "panmanWriter.hpp"
#include "message_manager.hpp"
#include "timer.hpp"

#include <tbb/global_control.h>
#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>
#include <chrono>

po::options_description mainDesc("TWILIGHT Command Line Arguments", 120);

void parseArguments(int argc, char** argv)
{
    
    // Section: I/O
    po::options_description inputDesc("Inputs\n [1]  Build MGA From Scatch\n [2]  Add New Genomes to Existing MGA");
    inputDesc.add_options()
        ("tree,t", po::value<std::string>(), "Guide tree (Newick format): required for [1]; optional for [3].")
        ("sequences,i", po::value<std::string>(), "Unaligned sequences file (FASTA format): required for [1] and [3].");
        
    po::options_description outputDesc("Output/File Options");
    outputDesc.add_options()
        ("output,o", po::value<std::string>(), "Output file format: maf or panman. Default: output both formats.")
        ("msa-threshold", po::value<double>()->default_value(0.95), "Minimum sequence occupancy threshold (0.0 - 1.0) for blocks in concatenated MSA. Default: 0.95.")
        ("msa-min-len", po::value<int>()->default_value(500), "Minimum block consensus length (bp) for blocks in concatenated MSA. Default: 500.")
        ("dir,d", po::value<std::string>(), "Directory for storing intermediate and final output files.")
        ("prefix,p", po::value<std::string>(), "Prefix for output files. Default: derived from input sequence file name.")
        ("log", po::value<std::string>(), "Log file path for output messages.")
        ("overwrite", "Force overwriting. Remove the current dir.");

    // Section: Hardware
    po::options_description hardwareDesc("Hardware Usage");
    hardwareDesc.add_options()
        ("cpu,C", po::value<int>(), "Number of CPU cores. Default: all available cores.");

    // Section: Alignment Parameters
    po::options_description alignDesc("Alignment Options and Parameters");
    alignDesc.add_options()
        ("rooted", "Keep the original tree root (disable automatic re-rooting for parallelism)")
        ("circular", "Treat input genomes as circular and orient them relative to the reference genome before alignment.")
        ("ref,r", po::value<std::string>(), "Reference sequence name for circular orientation (default: longest sequence).")
        ("reference", po::value<std::string>(), "Alias for --ref.")
        ("output-oriented", "Write the oriented genome sequences in FASTA format.")
        ("mm-asm", po::value<int>()->default_value(5), "Minimap2 asm preset divergence level: 5, 10, or 20 (corresponding to asm5, asm10, asm20). Default: 5.")
        ("use-medoid", "Use medoid sequences for intermediate alignment to reduce repetition and chimeric noise, then realign to consensus.");

    po::options_description seqFilterDesc("Sequence Filtering Options");
    seqFilterDesc.add_options()
        ("length-deviation", po::value<float>(), "Sequences whose lengths deviate from the median by more than the specified fraction will be deferred or excluded.")
        ("max-ambig", po::value<float>()->default_value(0.1), "Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded.")
        ("max-len", po::value<int>(), "Sequences longer than max-len will be deferred or excluded.")
        ("min-len", po::value<int>(), "Sequences shorter than min-len will be deferred or excluded.")
        ("filter", "Exclude sequences with high ambiguity or length deviation.")
        ("write-filtered", "Write the filtered sequences in FASTA format to the output directory.");
        
    po::options_description generalDesc("General");
    generalDesc.add_options()
        ("check", "Check the final alignment. Sequences with no legal alignment will be displayed.")
        ("debug", "Run internal debugging checks and validation routines.")
        ("verbose,v", "Print out every detail process.")
        ("quiet,q", "Quiet mode")
        ("help,h", "Print help messages.")
        ("version,V", "Show program version.");

    // Setup boost::program_options
    mainDesc.add(inputDesc).add(outputDesc).add(hardwareDesc).add(alignDesc).add(seqFilterDesc).add(generalDesc);
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
        if(vm.count("help")) {
            return 0;
        }
        else {
            mga::msg::error("%s\nFor more details, please use the --help option.\n", e.what());
            return 1;
        }
    }
    if(vm.count("help") || argc == 1) {
        std::ostringstream ss;
        ss << mainDesc << "\n";
        mga::msg::print(ss.str());
        return 0;
    }

    Option option(vm);

    tbb::global_control init(tbb::global_control::max_allowed_parallelism, option.cpuNum);


    // Read Tree
    mga::msg::debug("Reading phylogenetic tree from %s...\n", option.treeFile.c_str());
    auto treeStart = std::chrono::high_resolution_clock::now();
    phylogeny::Tree T(option.treeFile);
    phylogeny::Tree subT(T.root.get(), true);
    auto treeEnd = std::chrono::high_resolution_clock::now();
    auto treeDuration = std::chrono::duration_cast<std::chrono::microseconds>(treeEnd - treeStart);
    mga::msg::info("Tree (Depth - Original: %d, Binary: %d, Reroot: %d) read in %ld us\n",
                   subT.depthOriginal, subT.depthBinary, subT.depthReroot, treeDuration.count());

    // Read Sequences
    auto manager = mga::io::readSequences(option.seqFile, option, subT);
    if (option.circular) manager->orientCircularGenomes(option, option.refSeq);

    mga::progressive::buildMGA(subT, option, manager.get(), 0);

    auto final_set = manager.get()->getBlockSet(manager.get()->getBlockSets().begin()->first);
    if (final_set) {
        // Rename final blockSet to "Root" instead of internal intermediate node ID (e.g. node_3)
        manager->changeBlockSetId(final_set->getId(), "Root");
    }

    if (option.verbose && final_set) {
        final_set->debugValidateQuality(false);
    }

    // Output MAF and Concatenated MSA
    auto outputStart = std::chrono::high_resolution_clock::now();
    auto seqLengths = manager->getSequenceLengths();
    if (option.outputMAF) {
        mga::io::writeMAF(final_set, option.tempDir + "/" + option.prefix + ".maf", seqLengths);
        mga::io::writeLinearMAF(final_set, option.tempDir + "/" + option.prefix + ".linear.maf", seqLengths);
    }
    mga::io::writeConcatMSA(final_set, option.tempDir + "/" + option.prefix + ".concat_msa.fa", option.msaThreshold, option.msaMinLen);

    std::chrono::nanoseconds refineTime(0);
    
    if (final_set && option.outputPanman) {
        std::string panmanOutputPath = option.tempDir + "/" + option.prefix + ".panman";
        mga::panman::PanmanWriter::writePanMAN(final_set, &subT, "", panmanOutputPath);
    }

    auto outputEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds outputTime = outputEnd - outputStart;
    
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    mga::msg::info("Total execution completed in %.6f s\n", static_cast<double>(mainTime.count()) / 1e9);
    return 0;
}
