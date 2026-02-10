#ifndef MGA_HPP
#include "mga.hpp"
#endif

#ifndef OPTION_HPP
#include "option.hpp"
#endif

#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

#include <tbb/global_control.h>
#include <boost/filesystem.hpp>
#include <chrono>
namespace fs = boost::filesystem;

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
        ("overwrite", "Force overwriting the output file.");

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
        ("wildcard,w", "Treat unknown or ambiguous bases as wildcards and align them to usual letters.")
        ("rooted", "Keep the original tree root (disable automatic re-rooting for parallelism)")
        ("prune", "Prune the input guide tree based on the presence of unaligned sequences.")
        ("write-prune", "Write the pruned tree to the output directory.");

    po::options_description seqFilterDesc("Sequence Filtering Options");
    seqFilterDesc.add_options()
        // ("no-align-gappy", "Do not align gappy columns. This will create a longer MSA (larger file).")
        // ("psgop", po::value<std::string>()->default_value("y"), "y: Enable, n: Disable position-specific gap open penalty.")
        ("length-deviation", po::value<float>(), "Sequences whose lengths deviate from the median by more than the specified fraction will be deferred or excluded.")
        ("max-ambig", po::value<float>()->default_value(0.1), "Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded.")
        ("max-len", po::value<int>(), "Sequences longer than max-len will be deferred or excluded.")
        ("min-len", po::value<int>(), "Sequences shorter than min-len will be deferred or excluded.")
        ("filter", "Exclude sequences with high ambiguity or length deviation.")
        ("write-filtered", "Write the filtered sequences in FASTA format to the output directory.");
        
    po::options_description generalDesc("General");
    generalDesc.add_options()
        ("check", "Check the final alignment. Sequences with no legal alignment will be displayed.")
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
            std::cerr << "Error: " << e.what() << std::endl;
            std::cerr << "For more details, please use the --help option.\n";
            return 1;
        }
    }
    if(vm.count("help") || argc == 1) {
        std::cerr << mainDesc << std::endl;
        return 0;
    }

    Option option(vm);

    // Read Tree
    phylogeny::Tree T(option.treeFile);
    phylogeny::Tree subT(T.root.get(), true);

    // Read Sequences
    BlockManager& manager = mga::io::readSequences(option.seqFile, option, subT);
    
    manager.printStatus();

    mga::progressive::msaOnSubtree(subT, option, manager, 0);

    
    

    // mga::identifyPrimaryAlignments(alignments, chains);
    // mga::detectDuplications(alignments);
    // mga::parser::parsePAF(vm["alignment"].as<std::string>());
    
    
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cerr << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
