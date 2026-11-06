#include "alignment.hpp"
#include "mga.hpp"
#include "option.hpp"
#include "timer.hpp"
#include "type.hpp"
#include "message_manager.hpp"
#include "medoid_realigner.hpp"

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
  auto processPair = [&](size_t i) {
    auto &pair = alnPairs[i];

    Node *node1 = pair.first;
    Node *node2 = pair.second;

    auto BlockSet1 = blockManager->getBlockSet(node1->identifier);
    auto BlockSet2 = blockManager->getBlockSet(node2->identifier);

    if (!BlockSet1 || !BlockSet2) {
      mga::msg::warn("Skipping alignment pair due to missing blockId.\n");
      return;
    }

    BlockSet1->setTree(&tree);
    BlockSet2->setTree(&tree);

    // Adaptive Minimap2 preset fallback: start at configured mmAsm, fallback to asm10 and asm20 if 0 alignments
    std::vector<int> asmLevels;
    asmLevels.push_back(option.mmAsm);
    if (option.mmAsm < 10) asmLevels.push_back(10);
    if (option.mmAsm < 20) asmLevels.push_back(20);

    Alignments mainAlignments;

    for (size_t asmIdx = 0; asmIdx < asmLevels.size(); ++asmIdx) {
      int currentAsm = asmLevels[asmIdx];
      double minimap2_start = getTimeMs();

      if (option.useMedoid) {
        const auto& medoid1 = BlockSet1->getAncestralMedoidSequence();
        const auto& medoid2 = BlockSet2->getAncestralMedoidSequence();

        auto medAlignments = runMinimap2(medoid1, medoid2, BlockSet1->getId(), BlockSet2->getId(), option, mm_config_whole_genome(currentAsm));
        double minimap2_end = getTimeMs();
        double minimap2_dur = minimap2_end - minimap2_start;
        global_timer.add("Minimap2", minimap2_dur);

        double realign_dur = 0.0;
        if (!medAlignments.empty()) {
          double realign_start = getTimeMs();
          mainAlignments = realignMedoidAlignmentsToConsensus(medAlignments, BlockSet1, BlockSet2, option);
          double realign_end = getTimeMs();
          realign_dur = realign_end - realign_start;
          global_timer.add("Realignment", realign_dur);
        } else {
          mainAlignments.clear();
        }

        if (currentAsm != asmLevels.front()) {
          mga::msg::debug("Minimap2(Medoid, Fallback asm%d) [Len: (%zu, %zu), Count: %zu, MM2: %.2f ms, Realign: %.2f ms]\n",
                         currentAsm, medoid1.size(), medoid2.size(), mainAlignments.size(), minimap2_dur, realign_dur);
        } else {
          mga::msg::debug("Minimap2(Medoid) [Len: (%zu, %zu), Count: %zu, MM2: %.2f ms, Realign: %.2f ms]\n",
                         medoid1.size(), medoid2.size(), mainAlignments.size(), minimap2_dur, realign_dur);
        }
      } else {
        const auto& consensus1 = BlockSet1->getAncestralSequence();
        const auto& consensus2 = BlockSet2->getAncestralSequence();

        mainAlignments = runMinimap2(consensus1, consensus2, BlockSet1->getId(), BlockSet2->getId(), option, mm_config_whole_genome(currentAsm));
        double minimap2_end = getTimeMs();
        double minimap2_dur = minimap2_end - minimap2_start;
        global_timer.add("Minimap2", minimap2_dur);
      }

      // If alignments were found, exit fallback search
      if (!mainAlignments.empty()) {
        break;
      }

      if (asmIdx + 1 < asmLevels.size()) {
        mga::msg::debug("Minimap2 produced 0 alignments with asm%d, attempting fallback to asm%d...\n",
                       currentAsm, asmLevels[asmIdx + 1]);
      }
    }

    if (mainAlignments.empty()) {
      mga::msg::warn("No alignments produced by minimap2 across all presets. Proceeding to merge as distant blocks.\n");
    }

    AlignmentCollection alnCollection(mainAlignments, BlockSet1, BlockSet2);

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

    double dbg_start = getTimeMs();
    if (option.debug) mergeBlockSet->debugValidateQuality(false);
    double dbg_end = getTimeMs();
    double dbg_dur = dbg_end - dbg_start;
    global_timer.add("Debug", dbg_dur);

    blockManager->removeBlockSet(node1->identifier);
    blockManager->removeBlockSet(node2->identifier);
  };

  // In debug mode, run pairs sequentially to avoid interleaved debug prints.
  // In normal mode, run pairs in parallel across independent tree branches.
  if (option.debug) {
    for (size_t i = 0; i < alnPairs.size(); ++i) {
      processPair(i);
    }
  } else {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, alnPairs.size()),
                      [&](const tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i < range.end(); ++i) {
        processPair(i);
      }
    });
  }

  // global_timer.print();
}