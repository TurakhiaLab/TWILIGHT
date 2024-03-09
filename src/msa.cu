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

#include "../src/treePartition.cpp"

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

void printTree(Node* node)
{
    std::cout << node->identifier << ": " << node->branchLength << "\t" << node->grpID << std::endl; 

    if (node->children.size() == 0) return;
    for (auto &c: node->children) printTree(c);
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

    // printTree(T->root);

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

void checkAlignment(std::vector<std::string>& ref)
{
    size_t len = 0;
    bool set = false;
    // std::cout << "Checking alignment....\n";
    for (auto &r: ref)
    {
        if (!set) len = r.size();
        if (r.size() != len)
        {
            fprintf(stderr, "Error: Alignment Size do not match\n");
        }
    }
}

void msaPostOrderTraversal(Node* node, msa::utility* util, Params& param, int grpID)
{
    // if (!(node->grpID==-1 || node->grpID==grpID)) return;

    if (node->children.size()==0) 
    {
        node->msa.push_back(util->seqs[node->identifier]);
        return;
    }

    for (auto& child: node->children) msaPostOrderTraversal(child, util, param, grpID);

    std::pair<std::vector<std::string>, std::vector<std::string>> alignments;
    std::vector<std::string> ref;
    
    size_t childIndex = 0;
    
    for (childIndex=0; childIndex<node->children.size(); childIndex++)
    {
        // if (node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID)
            ref = node->children[childIndex]->msa;
            break;
    }
    
    
    for (size_t i=childIndex+1; i<node->children.size(); i++)
    {
        /*
        if (!(node->children[i]->grpID == -1 || node->children[i]->grpID == grpID))
        {
            continue;
        }
        */
        std::vector<std::string> query = node->children[i]->msa;
        alignGrpToGrp(ref, query, param, alignments);
        ref.clear();
        for (auto &s: alignments.first) ref.push_back(s);
        for (auto &s: alignments.second) ref.push_back(s);
        alignments.first.clear(); alignments.second.clear();
    }
    node->msa = ref;

    return;
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
    Node* root = T->root;

    // Read Input Sequences (Fasta format)
    readSequences(vm, util);
    
    // Test alignSeqToSeq
    Params param(2,1,-2,-1,0,0);


    auto msaStart = std::chrono::high_resolution_clock::now();
    msaPostOrderTraversal(root, util, param, 0);
    auto msaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds msaTime = msaEnd - msaStart;
    std::cout << "MSA in: " <<  msaTime.count()/1000000 << " ms\n";

    return 0;
}

