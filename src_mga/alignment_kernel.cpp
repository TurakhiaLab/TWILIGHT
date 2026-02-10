#ifndef MGA_HPP
#include "mga.hpp"
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>


void mga::progressive::alignmentKernel(NodePairVec& alnPairs, BlockManager& blockManager, Option& option) {
    for (auto& pair : alnPairs) {
        Node* node1 = pair.first;
        Node* node2 = pair.second;

        if (node1->blockId == 0 || node2->blockId == 0) {
            std::cerr << "Warning: Skipping alignment pair due to missing blockId." << std::endl;
            continue;
        }

        Block::ID block1_id = node1->blockId;
        Block::ID block2_id = node2->blockId;

        auto block1 = blockManager.getBlock(block1_id);
        auto block2 = blockManager.getBlock(block2_id);

        if (!block1 || !block2) {
            std::cerr << "Error: Could not find blocks for alignment." << std::endl;
            continue;
        }

        // 1. Prepare temporary files for minimap2
        std::string temp_dir = option.tempDir;
        std::string ref_fasta_path = temp_dir + "/ref_" + std::to_string(block1_id) + ".fa";
        std::string qry_fasta_path = temp_dir + "/qry_" + std::to_string(block2_id) + ".fa";
        std::string sam_output_path = temp_dir + "/output_" + std::to_string(block1_id) + "_" + std::to_string(block2_id) + ".sam";

        std::ofstream ref_fasta_file(ref_fasta_path);
        ref_fasta_file << block1->getConsensusAsFasta("ref");
        ref_fasta_file.close();

        std::ofstream qry_fasta_file(qry_fasta_path);
        qry_fasta_file << block2->getConsensusAsFasta("qry");
        qry_fasta_file.close();

        // 2. Run minimap2
        const char* home_dir = getenv("HOME");
        if (home_dir == nullptr) {
            std::cerr << "Error: Could not get HOME directory." << std::endl;
            continue;
        }
        std::string minimap2_path = std::string(home_dir) + "/bin/minimap2";
        std::string command = minimap2_path + " -a " + ref_fasta_path + " " + qry_fasta_path + " > " + sam_output_path;
        
        int system_ret = system(command.c_str());
        if (system_ret != 0) {
            std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
            // Cleanup and continue
            remove(ref_fasta_path.c_str());
            remove(qry_fasta_path.c_str());
            remove(sam_output_path.c_str());
            continue;
        }

        // 3. Parse SAM output
        alnVec alignments = parser::parseMinimap2(sam_output_path, true);
        if (alignments.empty()) {
            std::cerr << "Warning: No alignments produced by minimap2." << std::endl;
        } else {
            // 4. Process alignments
            chainVec chains = getAlignmentChains(alignments);
            identifyPrimaryAlignments(alignments, chains);
            detectDuplications(alignments);
        }

        // 5. Merge blocks
        bool merged = blockManager.mergeBlocks(block1_id, block2_id);
        if (!merged) {
            std::cerr << "Error: Failed to merge blocks " << block1_id << " and " << block2_id << std::endl;
        } else {
            // 6. Update parent node (node1)
            node1->blockId = block1_id; // It now represents the merged block
            
            // Update other properties of node1 if necessary
            auto mergedBlock = blockManager.getBlock(block1_id);
            if(mergedBlock) {
                 node1->alnLen = mergedBlock->getConsensus().length();
                 node1->alnNum += node2->alnNum;
            }
        }

        // 7. Cleanup
        if (option.deleteTemp) {
            remove(ref_fasta_path.c_str());
            remove(qry_fasta_path.c_str());
            remove(sam_output_path.c_str());
        }
    }
}
