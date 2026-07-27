#include "alignment.hpp"
#include "mga.hpp"
#include "option.hpp"
#include "timer.hpp"
#include "type.hpp"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

void mga::progressive::alignmentKernel(NodePairs &alnPairs,
                                       BlockManager *blockManager,
                                       Option &option, Tree &tree) {
  std::cout << "Total Pairs: " << alnPairs.size() << '\n';

  for (auto &pair : alnPairs) {

    Node *node1 = pair.first;
    Node *node2 = pair.second;

    auto BlockSet1 = blockManager->getBlockSet(node1->identifier);
    auto BlockSet2 = blockManager->getBlockSet(node2->identifier);

    if (!BlockSet1 || !BlockSet2) {
      std::cerr << "Warning: Skipping alignment pair due to missing blockId."
                << std::endl;
      continue;
    }

    BlockSet1->setTree(&tree);
    BlockSet2->setTree(&tree);

    const auto& consensus1 = BlockSet1->getAncestralSequence();
    const auto& consensus2 = BlockSet2->getAncestralSequence();

    auto minimap2_start = std::chrono::high_resolution_clock::now();
    auto mainAlignments = runMinimap2(consensus1, consensus2, BlockSet1->getId(), BlockSet2->getId(), option, true, true);
    auto minimap2_end = std::chrono::high_resolution_clock::now();
    std::cout << "Minimap2 [Len: (" << consensus1.size() << "," << consensus2.size() << "), Count: " << mainAlignments.size() << ", Runtime: " << (std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count()) << " us]\n";
    
    if (mainAlignments.empty()) {
      std::cerr << "Warning: No alignments produced by minimap2." << std::endl;
    } else {

      AlignmentCollection alnCollection(mainAlignments, BlockSet1, BlockSet2);

      std::cout << node1->identifier << '\n';
      // std::cout << "Number of Alignments: " << mainAlignments.size() << '\n';
      // std::cout << BlockSet1->getAllBlocks().size() << '\n';
      // std::cout << BlockSet2->getAllBlocks().size() << '\n';

      BlockSetID parentID = (!tree.allNodes[BlockSet1->getId()]->parent) ? "Root" : tree.allNodes[BlockSet1->getId()]->parent->identifier;

      auto merge_start = std::chrono::high_resolution_clock::now();
      BlockSet *mergeBlockSet =
          blockManager->merge(BlockSet1, BlockSet2, alnCollection, parentID,
                              &tree);
      auto merge_end = std::chrono::high_resolution_clock::now();
      option.merge_time += std::chrono::duration_cast<std::chrono::milliseconds>(merge_end - merge_start).count();

      // mergeBlockSet->debugValidateSegments(false);
      // mergeBlockSet->debugValidateLinkages(false);

      std::cout << BlockSet1->getId() << '\t' << BlockSet2->getId() << " -> " << mergeBlockSet->getId() << '\n';
      std::cout << "Merged Block Count: "
                << mergeBlockSet->getAllBlocks().size() << '\n';
      std::cout << mergeBlockSet->getId() << '\n';

      auto refine_1 = std::chrono::high_resolution_clock::now();
      // mergeBlockSet->reconnectBlocks();
      // mergeBlockSet->refineGraph();
      // if (mergeBlockSet->getSequenceCount() >= 10)
      // mergeBlockSet->realignAllToAll(option.tempDir);
      // mergeBlockSet->writeMAF(option.tempDir+"/"+mergeBlockSet->getId()+".maf");
      auto refine_2 = std::chrono::high_resolution_clock::now();
      // option.refineBlockSet_time +=
      // std::chrono::duration_cast<std::chrono::milliseconds>(refine_2 -
      // refine_1).count();
      option.refineBlock_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(refine_2 -
                                                                refine_1)
              .count();

      mga::io::writeMAF(mergeBlockSet,
                        option.tempDir + "/" + mergeBlockSet->getId() + ".maf",
                        true);

      // if (mergeBlockSet->getSequenceCount() > 2)
      // mergeBlockSet->print(std::cout);

      auto debug_start2 = std::chrono::high_resolution_clock::now();
      mergeBlockSet->debugValidateSegments(false);
      // mergeBlockSet->debugValidateLinkages(false);
      mergeBlockSet->debugValidateQuality(false);
      // mergeBlockSet->debugValidateLinearizedBlocks(true);
      mergeBlockSet->debugValidateBlocks(true);
      mergeBlockSet->debugValidateSequences(blockManager);
      // mergeBlockSet->debugValidateBubble(false);
      auto debug_end2 = std::chrono::high_resolution_clock::now();
      option.debug_time +=
          std::chrono::duration_cast<std::chrono::milliseconds>(debug_end2 -
                                                                debug_start2)
              .count();

      blockManager->removeBlockSet(node1->identifier);
      blockManager->removeBlockSet(node2->identifier);
      // blockManager->changeBlockSetId(mergeBlockSet->getId(), node1->identifier); 
      // mergeBlockSet->setDistantBlocks(tree);
      // mergeBlockSet->print(std::cout);
      global_timer.print();
    }

  

    // std::remove(c1c2_path.c_str());
  }

  // blockManager->updateLongestSequences();

  std::cout << "--- Profiling Results (ms) ---\n";
  std::cout << "Minimap2 Time:         " << option.minimap2_time << " ms\n";
  std::cout << "Merge Time:            " << option.merge_time << " ms\n";
  std::cout << "Refine(BlockSet) Time: " << option.refineBlockSet_time
            << " ms\n";
  std::cout << "Refine(Block) Time:    " << option.refineBlock_time << " ms\n";
  std::cout << "Debug Time:            " << option.debug_time << " ms\n";
  std::cout << "------------------------------\n";

  // if (alnPairs.size() < 3) exit(1);
}