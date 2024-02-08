#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"

#ifndef TREE_HPP
#include "../src/tree.hpp"
#endif

#ifndef MSA_HPP
#include "../src/msa.hpp"
#endif

#ifndef ALIGN_HPP
#include "../src/align.hpp"
#endif

namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("MSA Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,s", po::value<std::string>()->required(), "Tip sequences - Fasta format (required)")
        ("help,h", "Print help messages");

}

Tree* readNewick(po::variables_map& vm)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::string treeFileName = vm["tree"].as<std::string>();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; inputStream >> newick;
    Tree *T = new Tree(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Newick string read in: " <<  treeBuiltTime.count() << " ns\n";

    return T;
    
}

void readSequences(po::variables_map& vm, msa::utility* util)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["sequences"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        util->seqs[kseq_rd->name.s] = std::string(kseq_rd->seq.s, seqLen);
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}


int main(int argc, char** argv) {

    std::string seqFileName;
    

    parseArguments(argc, argv);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << mainDesc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    // Define MSA utility
    msa::utility* util = new msa::utility;

    // Read Tree (Newick String)
    Tree* T = readNewick(vm);

    // Read Input Sequences (Fasta format)
    readSequences(vm, util);

    // Test alignSeqToSeq
    std::string ref = util->seqs["Vic"];
    std::string query = util->seqs["Lox"];
    Params param(2,1,-2,-1,0,0);
    // std::pair<std::string, std::string> alignment;
    // alignSeqToSeq(ref, query, param, alignment);
    // std::cout<< alignment.first << "\n" << alignment.second << "\n";

    std::vector<std::string> refs, querys;
    refs.push_back(ref);
    querys.push_back(query);
    std::pair<std::vector<std::string>, std::vector<std::string>> alignments;
    alignGrpToGrp(refs, querys, param, alignments);
    for (auto s: alignments.first) std::cout << s << "\n";
    for (auto s: alignments.second) std::cout << s << "\n";

    return 0;
}

