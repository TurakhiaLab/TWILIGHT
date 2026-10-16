#include "alignment.hpp"
#include "mga.hpp"
#include "option.hpp"
#include "timer.hpp"
#include "type.hpp"
#include "message_manager.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <time.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Helper function to get monotonic time in milliseconds using POSIX time.h
static inline double getTimeMs() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return static_cast<double>(ts.tv_sec) * 1000.0 + static_cast<double>(ts.tv_nsec) / 1000000.0;
}

void mga::progressive::alignmentKernel(NodePairs &alnPairs, BlockManager *blockManager, Option &option, Tree &tree) {
  
  mga::msg::debug("Total Pairs: %zu\n", alnPairs.size());

  // Execute alignment pairs in parallel across independent tree branches
  tbb::parallel_for(tbb::blocked_range<size_t>(0, alnPairs.size()),
                    [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      auto &pair = alnPairs[i];

      Node *node1 = pair.first;
      Node *node2 = pair.second;

      auto BlockSet1 = blockManager->getBlockSet(node1->identifier);
      auto BlockSet2 = blockManager->getBlockSet(node2->identifier);

      if (!BlockSet1 || !BlockSet2) {
        mga::msg::warn("Skipping alignment pair due to missing blockId.\n");
        continue;
      }

      BlockSet1->setTree(&tree);
      BlockSet2->setTree(&tree);

      const auto& consensus1 = BlockSet1->getAncestralSequence();
      const auto& consensus2 = BlockSet2->getAncestralSequence();

      double minimap2_start = getTimeMs();
      auto mainAlignments = runMinimap2(consensus1, consensus2, BlockSet1->getId(), BlockSet2->getId(), option, mm_config_whole_genome(option), true);
      double minimap2_end = getTimeMs();
      double minimap2_dur = minimap2_end - minimap2_start;
      global_timer.add("Minimap2", minimap2_dur);

      mga::msg::info("Minimap2 [Len: (%zu, %zu), Count: %zu, Runtime: %.2f ms]\n",
                     consensus1.size(), consensus2.size(), mainAlignments.size(), minimap2_dur);
      
      if (mainAlignments.empty()) {
        mga::msg::warn("No alignments produced by minimap2.\n");
      } else {

        AlignmentCollection alnCollection(mainAlignments, BlockSet1, BlockSet2);

        mga::msg::info("%s\n", node1->identifier.c_str());

        // Thread-safe lookup of parent node ID from tree
        auto itNode = tree.allNodes.find(BlockSet1->getId());
        BlockSetID parentID = (itNode == tree.allNodes.end() || !itNode->second->parent)
                                  ? "Root"
                                  : itNode->second->parent->identifier;

        double merge_start = getTimeMs();
        BlockSet *mergeBlockSet =
            blockManager->merge(BlockSet1, BlockSet2, alnCollection, parentID,
                                &tree, &option);
        double merge_end = getTimeMs();
        double merge_dur = merge_end - merge_start;
        global_timer.add("Merge", merge_dur);

        mga::msg::info("%s\t%s -> %s\n", BlockSet1->getId().c_str(), BlockSet2->getId().c_str(), mergeBlockSet->getId().c_str());
        mga::msg::info("Merged Block Count: %zu\n", mergeBlockSet->getAllBlocks().size());
        mga::msg::info("%s\n", mergeBlockSet->getId().c_str());

        double dbg_start = getTimeMs();
        mergeBlockSet->debugValidateQuality(false);
        double dbg_end = getTimeMs();
        double dbg_dur = dbg_end - dbg_start;
        global_timer.add("Debug", dbg_dur);

        blockManager->removeBlockSet(node1->identifier);
        blockManager->removeBlockSet(node2->identifier);
      }
    }
  });

  global_timer.print();

  mga::msg::debug("--- Profiling Results (ms) ---\n");
  mga::msg::debug("Minimap2 Time:         %.2f ms\n", global_timer.getTime("Minimap2"));
  mga::msg::debug("Merge Time:            %.2f ms\n", global_timer.getTime("Merge"));
  mga::msg::debug("Refine(BlockSet) Time: %.2f ms\n", global_timer.getTime("Refine(BlockSet)"));
  mga::msg::debug("Refine(Block) Time:    %.2f ms\n", global_timer.getTime("Refine(Block)"));
  mga::msg::debug("Debug Time:            %.2f ms\n", global_timer.getTime("Debug"));
  mga::msg::debug("------------------------------\n");
}