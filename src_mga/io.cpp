#ifndef MGA_HPP
#include "mga.hpp"
#endif


#include "kseq.h"

#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <climits>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

KSEQ_INIT2(, gzFile, gzread);


BlockManager& mga::io::readSequences(std::string& fileName, Option& option, Tree& tree) {
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    
    std::string seqFileName = fileName;
    
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);

    // Initiate a Block
    BlockManager& blockManager = BlockManager::instance();

    int seqNum = 0;
    uint64_t totalLen = 0;
    int maxLen = 0;
    int minLen = INT_MAX;

    std::unordered_set<std::string> loadedNames;

    // =========================================================
    // Read sequences and turn into Blocks
    // =========================================================
    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        std::string seqName_full = kseq_rd->name.s;
        
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";

        // Check whether the sequence presented in the guide tree
        if (tree.allNodes.find(seqName_full) != tree.allNodes.end()) {
            seqName = seqName_full;
        } else if (tree.allNodes.find(seqName_noblank) != tree.allNodes.end()) {
            seqName = seqName_noblank;
        }

        if (seqName != "") {
            // Check sequence duplication
            if (loadedNames.find(seqName) != loadedNames.end()) {
                printf("WARNING: duplicate leaf names found in the sequence file! Leaf name: %s. Only the first occurrence will be kept.\n", seqName.c_str());
            } else {
                std::string seqContent = std::string(kseq_rd->seq.s, seqLen);

                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
                if (seqLen == 0) std::cerr << "Null sequences, " << seqName << '\n';
                totalLen += seqLen;
                seqNum++;

                // Placed or not
                tree.allNodes[seqName]->placed = false;
                
                // -----------------------------------------------------------
                // Initialize a block. 1 Block = 1 Sequence
                // -----------------------------------------------------------
                
                std::shared_ptr<Block> newBlock = blockManager.createBlock(seqContent);

                SequenceInfo info;
                info.sequence_id = seqName;
                info.start_coordinate = 0;
                info.end_coordinate = seqLen;

                newBlock->addSequence(info);

                loadedNames.insert(seqName);
            }
        }
    }
    
    kseq_destroy(kseq_rd);
    gzclose(f_rd);


    if (seqNum == 0) {
        std::cerr << "Error: no sequences were read from the input.\n";
        exit(1);
    }

    // =========================================================
    // Print Summary
    // =========================================================
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;

    uint32_t avgLen = (seqNum > 0) ? totalLen / seqNum : 0;

    std::cerr << "===== Sequence Summary =====\n";
    std::cerr << "Number : " << seqNum << '\n';
    std::cerr << "Max. Length: " << maxLen << '\n';
    std::cerr << "Min. Length: " << minLen << '\n';
    std::cerr << "Avg. Length: " << avgLen << '\n';
    std::cerr << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";

    return blockManager;
}