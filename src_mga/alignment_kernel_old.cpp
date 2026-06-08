#ifndef MGA_HPP
#include "mga.hpp"
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <chrono>




void mga::progressive::alignmentKernel(NodePairVec& alnPairs, BlockManager* blockManager, Option& option) {
    std::cout << "Total Pairs: " << alnPairs.size() << '\n'; 
    
    for (auto& pair : alnPairs) {    

        Node* node1 = pair.first;
        Node* node2 = pair.second;

        auto BlockSet1 = blockManager->getBlockSet(node1->identifier);
        auto BlockSet2 = blockManager->getBlockSet(node2->identifier);

        if (!BlockSet1 || !BlockSet2 ) {
            std::cerr << "Warning: Skipping alignment pair due to missing blockId." << std::endl;
            continue;
        }

        stringPairVec consensus1, consensus2;
        stringPairVec remaining1, remaining2;

        /*
        if (BlockSet1->getSequenceCount() == 1) {
            std::string seqName = BlockSet1->getAllBlocks().front()->getSequences().begin()->first;
            std::string rawSeq = BlockSet1->getAllBlocks().front()->getConsensus();
            
            BlockSet1->selfMapping(option);
            
            auto debug_start = std::chrono::high_resolution_clock::now();
            std::string reconstructedSeq = BlockSet1->reconstructSequence(seqName);
            auto debug_end = std::chrono::high_resolution_clock::now();
            option.debug_time += std::chrono::duration_cast<std::chrono::milliseconds>(debug_end - debug_start).count();

            if (rawSeq != reconstructedSeq) {
                std::cerr << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch. [" << seqName << "]\n";
                exit(1);  
            }
            // if (rawSeq == reconstructedSeq) std::cout << "  ✅ [PERFECT] Validation Passed! Reconstructed sequence perfectly matches the raw sequence. (Len: " << rawSeq.length() << " bp)\n";
            // else                            std::cerr << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch.\n";
        }
        
        
        if (BlockSet2->getSequenceCount() == 1) {
            std::string seqName = BlockSet2->getAllBlocks().front()->getSequences().begin()->first;
            std::string rawSeq = BlockSet2->getAllBlocks().front()->getConsensus();
            
            BlockSet2->selfMapping(option);
            
            auto debug_start = std::chrono::high_resolution_clock::now();
            std::string reconstructedSeq = BlockSet2->reconstructSequence(seqName);
            auto debug_end = std::chrono::high_resolution_clock::now();
            option.debug_time += std::chrono::duration_cast<std::chrono::milliseconds>(debug_end - debug_start).count();

            if (rawSeq != reconstructedSeq) {
                std::cerr << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch. [" << seqName << "]\n";
                exit(1);  
            }

            // if (rawSeq == reconstructedSeq) std::cout << "  ✅ [PERFECT] Validation Passed! Reconstructed sequence perfectly matches the raw sequence. (Len: " << rawSeq.length() << " bp)\n";
            // else                            std::cerr << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch.\n";
            
        }
        */

        BlockSet1->getRepresentativeAndRemaining(consensus1, remaining1);
        BlockSet2->getRepresentativeAndRemaining(consensus2, remaining2);

        // 1. Prepare temporary files for minimap2
        std::string temp_dir = option.tempDir;
        std::string consensus1_path = temp_dir + "/" + BlockSet1->getId() + ".fa";
        std::string consensus2_path = temp_dir + "/" + BlockSet2->getId() + ".fa";
        
        // if (!remaining1.empty()) consensus1.insert(consensus1.end(), remaining1.begin(), remaining1.end());
        // if (!remaining2.empty()) consensus2.insert(consensus2.end(), remaining2.begin(), remaining2.end());

        io::writeAlignment(consensus1_path, consensus1, false, false);
        io::writeAlignment(consensus2_path, consensus2, false, false);
        
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
        command = minimap2_path + " -cx asm5 " + consensus1_path + " " + consensus2_path + " > " + c1c2_path + " 2> /dev/null";
        
        auto minimap2_start = std::chrono::high_resolution_clock::now();
        system_ret = system(command.c_str());
        auto minimap2_end = std::chrono::high_resolution_clock::now();
        option.minimap2_time += std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count();
        
        if (system_ret != 0) {
            std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
        }
        
        // Clear sequence file
        std::remove(consensus1_path.c_str());
        std::remove(consensus2_path.c_str());
        // std::remove(remaining1_path.c_str());
        // std::remove(remaining2_path.c_str());
        // 3. Parse PAF output
        alnVec mainAlignments = parser::parseMinimap2PAF(c1c2_path);
        std::cout << "Number of Minimap2 Alignments: " << mainAlignments.size() << '\n';
        if (mainAlignments.empty()) {
            std::cerr << "Warning: No alignments produced by minimap2." << std::endl;
        } 
        else {
            // 4. Process alignments
            stringMap ref_seqs, qry_seqs;
            for (auto& seqs : consensus1) ref_seqs[seqs.first] = seqs.second;
            for (auto& seqs : consensus2) qry_seqs[seqs.first] = seqs.second;

            std::cout << node1->identifier << '\n';
            identifyPrimaryAlignments(mainAlignments, BlockSet1, BlockSet2, ref_seqs, qry_seqs);
            
            // detectDuplications(mainAlignments);
            // fillUnalignedRegions(mainAlignments, consensus1.begin()->second.size(), consensus2.begin()->second.size());
            // validateCoverage(mainAlignments, consensus1.begin()->second.size(), consensus2.begin()->second.size());
            
            // ==========================================
            // 【新增】：呼叫邊界吸附邏輯 (Boundary Snapping)
            // ==========================================
            std::string refSeq = consensus1.begin()->second;
            std::string qrySeq = consensus2.begin()->second;
            // int snap_threshold = 100; // 你設計的 100bp 容許範圍
            // snapAlignmentsToBlockBoundaries(mainAlignments, BlockSet1, BlockSet2, refSeq, qrySeq, snap_threshold);
            // ==========================================


            std::cout << "Number of Alignments: " << mainAlignments.size() << '\n';
            
            
            auto merge_start = std::chrono::high_resolution_clock::now();
            BlockSet* mergeBlockSet = blockManager->merge(BlockSet1, BlockSet2, mainAlignments);
            auto merge_end = std::chrono::high_resolution_clock::now();
            option.merge_time += std::chrono::duration_cast<std::chrono::milliseconds>(merge_end - merge_start).count();
            
            // mergeBlockSet->debugValidateSegments(false);
            // mergeBlockSet->debugValidateLinkages(false);

            std::cout << BlockSet1->getId() << '\t' << BlockSet2->getId() << '\n';
            std::cout << "Merged Block Count: " << mergeBlockSet->getAllBlocks().size() << '\n';
            std::cout << mergeBlockSet->getId() << '\n';
            
            auto refine_1 = std::chrono::high_resolution_clock::now();
            mergeBlockSet->refineGraph();
            // if (mergeBlockSet->getSequenceCount() >= 10) mergeBlockSet->realignAllToAll(option.tempDir);
            // mergeBlockSet->writeMAF(option.tempDir+"/"+mergeBlockSet->getId()+".maf");
            auto refine_2 = std::chrono::high_resolution_clock::now();
            // option.refineBlockSet_time += std::chrono::duration_cast<std::chrono::milliseconds>(refine_2 - refine_1).count();
            option.refineBlock_time += std::chrono::duration_cast<std::chrono::milliseconds>(refine_2 - refine_1).count();
            
            auto debug_start2 = std::chrono::high_resolution_clock::now();
            mergeBlockSet->debugValidateSegments(false);
            mergeBlockSet->debugValidateLinkages(false);
            mergeBlockSet->debugValidateQualityNew(false);
            // mergeBlockSet->debugValidateBubble(false);
            auto debug_end2 = std::chrono::high_resolution_clock::now();
            option.debug_time += std::chrono::duration_cast<std::chrono::milliseconds>(debug_end2 - debug_start2).count();

            
            blockManager->removeBlockSet(node1->identifier);
            blockManager->removeBlockSet(node2->identifier);
            blockManager->changeBlockSetId(mergeBlockSet->getId(), node1->identifier); 

        }
        // std::remove(c1c2_path.c_str());
    }
    blockManager->updateLongestSequences();
    
    std::cout << "--- Profiling Results (ms) ---\n";
    std::cout << "Minimap2 Time:         " << option.minimap2_time << " ms\n";
    std::cout << "Merge Time:            " << option.merge_time << " ms\n";
    std::cout << "Refine(BlockSet) Time: " << option.refineBlockSet_time << " ms\n";
    std::cout << "Refine(Block) Time:    " << option.refineBlock_time << " ms\n";
    std::cout << "Debug Time:            " << option.debug_time << " ms\n";
    std::cout << "------------------------------\n";
}