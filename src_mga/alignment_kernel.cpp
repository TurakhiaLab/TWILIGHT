#ifndef MGA_HPP
#include "mga.hpp"
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>




void mga::progressive::alignmentKernel(NodePairVec& alnPairs, BlockManager* blockManager, Option& option) {
    std::cout << "Total Pairs: " << alnPairs.size() << '\n'; 
    for (auto& pair : alnPairs) {
        
        Node* node1 = pair.first;
        Node* node2 = pair.second;

        auto BlockSet1 = blockManager->getBlockSet(node1->identifier);
        auto BlockSet2 = blockManager->getBlockSet(node2->identifier);

        // // Test Progressive
        // std::cout << "BlockSet 1 [" << BlockSet1->getId() << "] Longest: " << BlockSet1->getLongestSeq() << "\n";
        // for (auto& seq: BlockSet1->getSequences()) std::cout << seq << ':' << blockManager->getSequenceLength(seq) << '\n';
        // std::cout << "BlockSet 2 [" << BlockSet2->getId() << "] Longest: " << BlockSet2->getLongestSeq() << "\n";
        // for (auto& seq: BlockSet2->getSequences()) std::cout << seq << ':' << blockManager->getSequenceLength(seq) << '\n';
        // BlockSet* resultSet = blockManager->createBlockSet("temp"); 
        // for (auto& seq: BlockSet1->getSequences()) resultSet->addSequence(seq);
        // for (auto& seq: BlockSet2->getSequences()) resultSet->addSequence(seq);
        // blockManager->removeBlockSet(node1->identifier);
        // blockManager->removeBlockSet(node2->identifier);
        // blockManager->changeBlockSetId(resultSet->getId(), node1->identifier); 
        // continue;

        if (!BlockSet1 || !BlockSet2 ) {
            std::cerr << "Warning: Skipping alignment pair due to missing blockId." << std::endl;
            continue;
        }

        stringPairVec consensus1, consensus2;
        stringPairVec remaining1, remaining2;

        if (BlockSet1->getSequenceCount() == 1) BlockSet1->selfMapping(option);
        BlockSet1->getRepresentativeAndRemaining(consensus1, remaining1);
        if (BlockSet2->getSequenceCount() == 1) BlockSet2->selfMapping(option);
        BlockSet2->getRepresentativeAndRemaining(consensus2, remaining2);

        // 1. Prepare temporary files for minimap2
        std::string temp_dir = option.tempDir;
        std::string consensus1_path = temp_dir + "/" + BlockSet1->getId() + ".fa";
        std::string consensus2_path = temp_dir + "/" + BlockSet2->getId() + ".fa";
        

        io::writeAlignment(consensus1_path, consensus1, false, false);
        // io::writeAlignment(consensus1_path, remaining1, false, true);
        io::writeAlignment(consensus2_path, consensus2, false, false);
        // io::writeAlignment(consensus2_path, remaining2, false, true);
        
        // whole-genome alignment between 1 and 2
        std::string c1c2_path = temp_dir + "/output_" + BlockSet1->getId() + "_" + BlockSet2->getId() + ".paf";

        // 2. Run minimap2
        const char* home_dir = getenv("HOME");
        if (home_dir == nullptr) {
            std::cerr << "Error: Could not get HOME directory." << std::endl;
            continue;
        }
        
        // std::string minimap2_path = std::string(home_dir) + "/bin/minimap2";
        std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
        std::string command;
        int system_ret; 

        // minimap2 -cx asm5 -g 500 -r 500 -n 5 -m 50 -N 20 -p 0.8 asm1.fa asm2.fa > aln.paf
        // whole-genome alignment between 1 and 2
        command = minimap2_path + " -cx asm10 -g 500 -r 500,500 -N 20 " + consensus1_path + " " + consensus2_path + " > " + c1c2_path;
        system_ret = system(command.c_str());
        if (system_ret != 0) {
            std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
        }
        
        // Clear sequence file
        // std::remove(consensus1_path.c_str());
        // std::remove(consensus2_path.c_str());
        // std::remove(remaining1_path.c_str());
        // std::remove(remaining2_path.c_str());
        // 3. Parse PAF output
        alnVec mainAlignments = parser::parseMinimap2PAF(c1c2_path);
        if (mainAlignments.empty()) {
            std::cerr << "Warning: No alignments produced by minimap2." << std::endl;
        } 
        else {
            // 4. Process alignments
            chainVec chains = getAlignmentChains(mainAlignments);
            identifyPrimaryAlignments(mainAlignments, chains);
            
            // detectDuplications(mainAlignments);
            // fillUnalignedRegions(mainAlignments, consensus1.begin()->second.size(), consensus2.begin()->second.size());
            // validateCoverage(mainAlignments, consensus1.begin()->second.size(), consensus2.begin()->second.size());
            
            std::cout << "Number of Alignments: " << mainAlignments.size() << '\n';
            
            BlockSet* mergeBlockSet = blockManager->merge(BlockSet1, BlockSet2, mainAlignments);

            mergeBlockSet->debugValidateSegments(false);
            mergeBlockSet->debugValidateLinkages(false);
            mergeBlockSet->refine();
            mergeBlockSet->debugValidateSegments(false);
            mergeBlockSet->debugValidateLinkages(false);
            mergeBlockSet->debugValidateQuality(false);

            blockManager->removeBlockSet(node1->identifier);
            blockManager->removeBlockSet(node2->identifier);
            blockManager->changeBlockSetId(mergeBlockSet->getId(), node1->identifier); 

            // // Align remaining blocks to the main blocksets

            // stringPairVec remaining = std::move(remaining1);
            // remaining.insert(remaining.end(), remaining2.begin(), remaining2.end());
            // stringPairVec mergedBlocks;
            // std::cout << remaining.size() << std::endl;
            // std::string main_path = temp_dir + "/" + mergeBlockSet->getId() + ".fa";
            // std::string remaining_path = temp_dir + "/" + mergeBlockSet->getId() + "_remaining.fa";
            // std::string addremaining_path = temp_dir + "/output_" + mergeBlockSet->getId() + "_addremaining.paf";

            // for (auto& block: mergeBlockSet->getAllBlocks()) {
            //     mergedBlocks.push_back({mergeBlockSet->getId()+"_"+std::to_string(block->getId()), block->getConsensus()});
            // }
            // io::writeAlignment(main_path, mergedBlocks, false, false);
            // io::writeAlignment(remaining_path, remaining, false, false);
        
            // if (!remaining.empty()) {
            //     command = minimap2_path + " -cx asm5 " + main_path + " " + remaining_path + " > " + addremaining_path;
            //     system_ret = system(command.c_str());
            //     if (system_ret != 0) {
            //         std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
            //     }
            //     alnVec addAlignments = parser::parseMinimap2PAF(addremaining_path);
            // }

            // std::cout << "Modify ID of " << mergeBlockSet->getId() << " to " << pair.first->identifier << std::endl;
            // blockManager->changeBlockSetId(mergeBlockSet->getId(), pair.first->identifier);
        }
    }
    blockManager->updateLongestSequences();
    
}
