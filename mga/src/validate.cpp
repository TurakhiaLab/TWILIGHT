
#include "block_manager.hpp"
#include "consensus.hpp"
#include "timer.hpp"

#include <cmath>
#include <iomanip>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

void BlockSet::debugValidateSegments(bool verbose) {
  global_timer.start("debug_validate_segments");

  if (verbose) {
    std::cout
        << "\n============================================================\n"
        << "=== BlockSet Debug Validation: " << ID << " ===\n"
        << "============================================================\n";
  }

  int totalBlocks = 0;
  int totalSegments = 0;
  int errorCount = 0;

  for (const auto &blockPair : blocks) {
    std::shared_ptr<Block> blk = blockPair.second;
    int consLen = blk->getConsensus().length();
    totalBlocks++;

    // How many segments are in the block
    int segmentsInBlock = 0;
    for (auto &seqPair : blk->getSequences()) {
      segmentsInBlock += seqPair.second.getSegments().size();
    }

    if (verbose)
      std::cout << "[Block ID: " << blk->getId() << "] "
                << "Consensus Len: " << consLen
                << " | Sequences: " << blk->getSequences().size()
                << " | Segments: " << segmentsInBlock << "\n";

    for (auto &seqPair : blk->getSequences()) {
      std::string seqName = seqPair.first;
      auto &segments = seqPair.second.getSegments();

      for (auto &segPairInner : segments) {
        Segment &seg = segPairInner.second;
        totalSegments++;

        int start = seg.getStart();
        int end = seg.getEnd();
        int coordDiff = std::abs(end - start);

        int gapLen = 0;
        for (auto &var : seg.getVariants()) {
          if (var.getType() == VariantType::GAP) {
            gapLen += (var.getEnd() - var.getStart());
          }
        }

        int calculatedLen = coordDiff + gapLen;

        if (calculatedLen != consLen) {
          std::cerr << "  ❌ [ERROR] At Block " << blk->getId()
                    << " | Seq: " << seqName << " | Seg: [" << start << ", "
                    << end << "] " << (seg.isReverse() ? "(-)" : "(+)")
                    << "\n      => CoordDiff (" << coordDiff << ") + Gaps ("
                    << gapLen << ") = " << calculatedLen << " != Consensus ("
                    << consLen << ")\n";
          errorCount++;
        }
      }
    }
  }

  if (verbose || errorCount > 0) {
    std::cout
        << "------------------------------------------------------------\n";
    std::cout << "Validation Complete. Checked " << totalBlocks
              << " blocks and " << totalSegments << " segments.\n";
  }
  if (errorCount == 0) {
    if (verbose)
      std::cout << "🎉 PERFECT! All segment lengths match their block "
                   "consensus length perfectly.\n";
  } else {
    std::cout << "🚨 CRITICAL: FOUND " << errorCount
              << " LENGTH MISMATCH ERROR(S)!\n";
    global_timer.stop("debug_validate_segments");
    exit(1);
  }
  if (verbose || errorCount > 0) {
    std::cout
        << "============================================================\n\n";
  }
  global_timer.stop("debug_validate_segments");
}

void BlockSet::debugValidateLinkages(bool verbose) {
  global_timer.start("debug_validate_linkages");
  if (verbose) {
    std::cout
        << "\n============================================================\n"
        << "=== Topology & Linkage Debug Validation: " << ID << " ===\n"
        << "============================================================\n";
  }

  // 🌟 1. 加入 copy 欄位，完美契合 VBlock 架構
  struct SegRef {
    Segment *seg;
    std::shared_ptr<Block> blk;
    int copy;
  };
  std::map<std::string, std::vector<SegRef>> seqTracks;

  // 1. Collect all segment
  for (const auto &blockPair : blocks) {
    for (auto &seqPair : blockPair.second->getSequences()) {
      for (auto &segPairInner : seqPair.second.getSegments()) {
        seqTracks[seqPair.first].push_back({
            &segPairInner.second, blockPair.second,
            segPairInner.second.getCopyCount() // 🌟 抓出 Copy ID
        });
      }
    }
  }

  int pointerErrorCount = 0;
  int coordinateGapCount = 0;

  // 🌟 輔助函數：安全地印出 VBlock 格式 (BlockID, cp: X)，防止 nullptr 崩潰
  auto getVBlockStr = [](std::shared_ptr<Block> b, int cp = -1) -> std::string {
    if (!b)
      return "NULL";
    if (cp != -1)
      return std::to_string(b->getId()) + "(cp:" + std::to_string(cp) + ")";
    return std::to_string(b->getId());
  };

  // 2. Check all sequences
  for (auto &trackPair : seqTracks) {
    auto &track = trackPair.second;

    // 🌟 2. 終極穩定排序 (Stable Sort)，防止重疊 Segment 順序錯亂
    std::sort(track.begin(), track.end(), [](const SegRef &a, const SegRef &b) {
      int min_a = std::min(a.seg->getStart(), a.seg->getEnd());
      int min_b = std::min(b.seg->getStart(), b.seg->getEnd());
      if (min_a != min_b)
        return min_a < min_b;

      int max_a = std::max(a.seg->getStart(), a.seg->getEnd());
      int max_b = std::max(b.seg->getStart(), b.seg->getEnd());
      if (max_a != max_b)
        return max_a < max_b;

      if (a.blk->getId() != b.blk->getId())
        return a.blk->getId() < b.blk->getId();
      return a.copy < b.copy;
    });

    if (verbose)
      std::cout << "Checking Sequence Track: " << trackPair.first << " ("
                << track.size() << " segments)\n";

    for (size_t i = 0; i < track.size(); ++i) {
      auto &currRef = track[i];

      if (i > 0) {
        auto &prevRef = track[i - 1];

        int prevEnd = std::max(prevRef.seg->getStart(), prevRef.seg->getEnd());
        int currStart =
            std::min(currRef.seg->getStart(), currRef.seg->getEnd());

        if (prevEnd != currStart) {
          std::cerr << "  ⚠️ [WARNING] Coordinate Gap: VBlock "
                    << getVBlockStr(prevRef.blk, prevRef.copy) << " end("
                    << prevEnd << ") != VBlock "
                    << getVBlockStr(currRef.blk, currRef.copy) << " start("
                    << currStart << ")\n";
          coordinateGapCount++;
        }

        std::shared_ptr<Block> expectedNextForPrev = currRef.blk;
        std::shared_ptr<Block> expectedPrevForCurr = prevRef.blk;
        std::shared_ptr<Block> actualNextForPrev;
        std::shared_ptr<Block> actualPrevForCurr;

        if (this->sequence_names.size() > 1) { // Pangenome graph
          actualNextForPrev = (!prevRef.seg->isReverse())
                                  ? prevRef.seg->getNextBlock().lock()
                                  : prevRef.seg->getPrevBlock().lock();
          actualPrevForCurr = (!currRef.seg->isReverse())
                                  ? currRef.seg->getPrevBlock().lock()
                                  : currRef.seg->getNextBlock().lock();
        } else { // Self-mapping
          actualNextForPrev = prevRef.seg->getNextBlock().lock();
          actualPrevForCurr = currRef.seg->getPrevBlock().lock();
        }

        // 🌟 3. 使用安全印出，避免 actualNextForPrev 為 NULL 時引發 Segfault
        if (actualNextForPrev != expectedNextForPrev) {
          std::cerr << "  ❌ [ERROR] Broken Pointer (Forward): VBlock "
                    << getVBlockStr(prevRef.blk, prevRef.copy)
                    << " does NOT point to VBlock "
                    << getVBlockStr(currRef.blk, currRef.copy)
                    << " (instead: " << getVBlockStr(actualNextForPrev)
                    << ")\n";
          pointerErrorCount++;
        }

        if (actualPrevForCurr != expectedPrevForCurr) {
          std::cerr << "  ❌ [ERROR] Broken Pointer (Backward): VBlock "
                    << getVBlockStr(currRef.blk, currRef.copy)
                    << " does NOT point back to VBlock "
                    << getVBlockStr(prevRef.blk, prevRef.copy)
                    << " (instead: " << getVBlockStr(actualPrevForCurr)
                    << ")\n";
          pointerErrorCount++;
        }
      }
    }
  }

  if (pointerErrorCount == 0 && coordinateGapCount == 0) {
    if (verbose)
      std::cout << "  🎉 PERFECT! All Graph linkages and coordinates are "
                   "contiguous and sound.\n";
  } else {
    std::cout
        << "------------------------------------------------------------\n";
    std::cout << "🚨 SUMMARY: Found " << pointerErrorCount
              << " Pointer Error(s) and " << coordinateGapCount
              << " Coordinate Gap(s)!\n";
    std::cout
        << "============================================================\n\n";
    global_timer.stop("debug_validate_linkages");
    exit(1);
  }
  global_timer.stop("debug_validate_linkages");
}

void BlockSet::debugValidateQuality(bool verbose) {
  global_timer.start("debug_validate_quality");
  std::cout
      << "\n============================================================\n"
      << "=== Pangenome Graph Quality Report (NEW): " << ID << " ===\n"
      << "============================================================\n";

  if (blocks.empty()) {
    std::cout << "  [Warning] Graph is empty. No metrics to calculate.\n";
    global_timer.stop("debug_validate_quality");
    return;
  }

  std::map<std::string, int> seqRealLengths;
  int totalSequenceCount = 0;

  for (auto &blkPair : blocks) {
    for (auto &seqPair : blkPair.second->getSequences()) {
      for (auto &segPair : seqPair.second.getSegments()) {
        seqRealLengths[seqPair.first] +=
            std::abs(segPair.second.getEnd() - segPair.second.getStart());
      }
    }
  }
  totalSequenceCount = seqRealLengths.size();

  int maxSeqLen = 0;
  std::string longestSeqName = "";
  for (const auto &kv : seqRealLengths) {
    if (kv.second > maxSeqLen) {
      maxSeqLen = kv.second;
      longestSeqName = kv.first;
    }
  }

  int threshold95 = std::ceil(totalSequenceCount * 0.95);
  int threshold90 = std::ceil(totalSequenceCount * 0.90);

  uint64_t totalConsensusLen = 0;
  uint64_t singletonLenSum = 0;
  int singletonCount = 0;

  uint64_t globalVarLen = 0;
  uint64_t globalDenominator = 0;

  uint64_t softCore90VarLen = 0;
  uint64_t softCore90Denominator = 0;
  uint64_t narrowCoreVarLen = 0;
  uint64_t narrowCoreDenominator = 0;

  std::vector<int> allBlockLengths;

  int coreBlocksCount = 0;
  int softCore95BlocksCount = 0;
  int softCore90BlocksCount = 0;
  int accessoryBlocksCount = 0;
  int narrowCoreBlocksCount = 0;
  int distantBlocksCount = 0;
  int totalVBlocksCount = 0;

  uint64_t coreLenSum = 0;
  uint64_t softCore95LenSum = 0;
  uint64_t softCore90LenSum = 0;
  uint64_t accessoryLenSum = 0;
  uint64_t narrowCoreLenSum = 0;
  uint64_t distantLenSum = 0;
  int distantMaxLen = 0;

  if (verbose)
    std::cout << "[Block-level Identity Info]\n";

  for (const auto &blkPair : blocks) {
    std::shared_ptr<Block> blk = blkPair.second;
    int consLen = blk->getConsensus().length();
    totalConsensusLen += consLen;
    allBlockLengths.push_back(consLen);

    int blockSegCount = 0;
    uint64_t blockVarLen = 0;
    std::set<std::string> uniqueSeqsInBlock;
    std::set<int> uniqueCopiesInBlock;

    for (auto &seqPair : blk->getSequences()) {
      uniqueSeqsInBlock.insert(seqPair.first);
      for (auto &segPair : seqPair.second.getSegments()) {
        blockSegCount++;
        uniqueCopiesInBlock.insert(segPair.second.getCopyCount());
        for (auto &var : segPair.second.getVariants()) {
          if (var.getType() == VariantType::SNV) {
            blockVarLen += 1;
          } else if (var.getType() == VariantType::GAP) {
            blockVarLen += (var.getEnd() - var.getStart());
          }
        }
      }
    }
    totalVBlocksCount += (uniqueCopiesInBlock.empty() ? 1 : uniqueCopiesInBlock.size());


    // 統計 Distant Block (若該 block 所有 copy 都是 distant)
    if (blk->isAllDistant()) {
      distantBlocksCount++;
      distantLenSum += consLen;
      if (consLen > distantMaxLen) {
        distantMaxLen = consLen;
      }
    }

    // 統計 Core, Soft-core (95%, 90%), Accessory
    int uniqueSeqCount = uniqueSeqsInBlock.size();

    if (uniqueSeqCount == totalSequenceCount &&
        blockSegCount == totalSequenceCount && totalSequenceCount > 1) {
      narrowCoreBlocksCount++;
      narrowCoreLenSum += consLen;
    }

    if (uniqueSeqCount == totalSequenceCount && totalSequenceCount > 1) {
      coreBlocksCount++;
      coreLenSum += consLen;
    }

    if (uniqueSeqCount >= threshold95 && totalSequenceCount > 1) {
      softCore95BlocksCount++;
      softCore95LenSum += consLen;
    }

    if (uniqueSeqCount >= threshold90 && totalSequenceCount > 1) {
      softCore90BlocksCount++;
      softCore90LenSum += consLen;
    }

    if (uniqueSeqCount < totalSequenceCount || totalSequenceCount <= 1) {
      accessoryBlocksCount++;
      accessoryLenSum += consLen;
    }

    // 判斷 Singleton
    if (blockSegCount == 1) {
      singletonLenSum += consLen;
      singletonCount++;
      // Singleton 不印出 Identity，也不列入 Global 計算
    } else if (blockSegCount > 1 && consLen > 0) {
      // 計算並印出非 Singleton 的 Identity
      uint64_t blockDenominator = (uint64_t)blockSegCount * consLen;
      double blockIdentity = 1.0 - ((double)blockVarLen / blockDenominator);

      if (verbose) {
        std::cout << "  ├─ Block ID: " << std::setw(6) << std::left
                  << blk->getId() << " | Len: " << std::setw(7) << consLen
                  << " | Segs: " << std::setw(3) << blockSegCount
                  << " | Identity: " << std::fixed << std::setprecision(4)
                  << blockIdentity << "\n";
      }
      // 累積至全局計算
      globalVarLen += blockVarLen;
      globalDenominator += blockDenominator;

      if (uniqueSeqCount >= threshold90 && totalSequenceCount > 1) {
        softCore90VarLen += blockVarLen;
        softCore90Denominator += blockDenominator;
      }
      if (uniqueSeqCount == totalSequenceCount &&
          blockSegCount == totalSequenceCount && totalSequenceCount > 1) {
        narrowCoreVarLen += blockVarLen;
        narrowCoreDenominator += blockDenominator;
      }
    }
  }

  // ---------------------------------------------------------
  // 3. 計算 N50
  // ---------------------------------------------------------
  std::sort(allBlockLengths.rbegin(), allBlockLengths.rend()); // 降序排列
  uint64_t runningSum = 0;
  int n50 = 0;
  for (int len : allBlockLengths) {
    runningSum += len;
    if (runningSum >= totalConsensusLen / 2) {
      n50 = len;
      break;
    }
  }

  // ---------------------------------------------------------
  // 4. 計算衍生指標與輸出
  // ---------------------------------------------------------
  double lenIncreaseRatio =
      (maxSeqLen > 0) ? ((double)totalConsensusLen / maxSeqLen - 1.0) * 100.0
                      : 0.0;
  double singletonRatio =
      (totalConsensusLen > 0)
          ? ((double)singletonLenSum / totalConsensusLen) * 100.0
          : 0.0;
  double globalIdentity =
      (globalDenominator > 0)
          ? (1.0 - ((double)globalVarLen / globalDenominator))
          : 0.0;
  double softCore90Identity =
      (softCore90Denominator > 0)
          ? (1.0 - ((double)softCore90VarLen / softCore90Denominator))
          : 0.0;
  double narrowCoreIdentity =
      (narrowCoreDenominator > 0)
          ? (1.0 - ((double)narrowCoreVarLen / narrowCoreDenominator))
          : 0.0;
  double distantAvgLen = (distantBlocksCount > 0)
                             ? (double)distantLenSum / distantBlocksCount
                             : 0.0;

  std::cout
      << "\n------------------------------------------------------------\n";
  std::cout << ">>> GRAPH METRICS SUMMARY <<<\n\n";

  std::cout << "[1. Sequence & Graph Size]\n";
  std::cout << "  - Total Sequences        : " << totalSequenceCount << "\n";
  std::cout << "  - Total Blocks           : " << blocks.size() << "\n";
  std::cout << "  - Total VBlocks          : " << totalVBlocksCount << "\n";
  std::cout << "  - Longest Input Sequence : " << maxSeqLen << " bp ("
            << longestSeqName << ")\n";
  std::cout << "  - Total Graph Length     : " << totalConsensusLen << " bp\n";
  std::cout << "  - Graph Size Inflation   : +" << std::fixed
            << std::setprecision(2) << lenIncreaseRatio << " %\n";

  std::cout << "\n[2. Fragmentation & Contiguity]\n";
  std::cout << "  - Block N50              : " << n50 << " bp\n";
  std::cout << "  - Singleton Blocks       : " << singletonCount << " blocks\n";
  std::cout << "  - Singleton Length       : " << singletonLenSum << " bp\n";
  std::cout << "  - Singleton Length Ratio : " << std::fixed
            << std::setprecision(2) << singletonRatio << " %\n";

  std::cout << "\n[3. Evolution & Conservation]\n";
  std::cout << "  - Narrow Core (1-to-1)   : " << narrowCoreBlocksCount
            << " blocks (" << narrowCoreLenSum << " bp)\n";
  std::cout << "  - Strict Core (100%)     : " << coreBlocksCount << " blocks ("
            << coreLenSum << " bp)\n";
  std::cout << "  - Soft Core (>= 95%)     : " << softCore95BlocksCount
            << " blocks (" << softCore95LenSum << " bp)\n";
  std::cout << "  - Soft Core (>= 90%)     : " << softCore90BlocksCount
            << " blocks (" << softCore90LenSum << " bp)\n";
  std::cout << "  - Distant Blocks         : " << distantBlocksCount
            << " blocks (" << distantLenSum << " bp, Max: " << distantMaxLen
            << " bp, Avg: " << std::fixed << std::setprecision(2)
            << distantAvgLen << " bp)\n";
  // std::cout << "  - Accessory (< 100%)     : " << accessoryBlocksCount << "
  // blocks (" << accessoryLenSum << " bp)\n";

  std::cout << "\n[4. Alignment Quality]\n";
  if (globalDenominator > 0) {
    std::cout << "  - Global Average Identity: " << std::fixed
              << std::setprecision(4) << globalIdentity
              << " (Excluded singletons)\n";
  } else {
    std::cout << "  - Global Average Identity: N/A (No valid multi-segment "
                 "blocks found)\n";
  }
  if (softCore90Denominator > 0) {
    std::cout << "  - Broad Core (>=90%) Identity: " << std::fixed
              << std::setprecision(4) << softCore90Identity << "\n";
  } else {
    std::cout << "  - Broad Core (>=90%) Identity: N/A\n";
  }
  if (narrowCoreDenominator > 0) {
    std::cout << "  - Narrow Core (1-to-1) Identity: " << std::fixed
              << std::setprecision(4) << narrowCoreIdentity << "\n";
  } else {
    std::cout << "  - Narrow Core (1-to-1) Identity: N/A\n";
  }

  std::cout
      << "============================================================\n\n";
  global_timer.stop("debug_validate_quality");
}

void BlockSet::debugValidateSequences(BlockManager *manager, bool verbose) {
  global_timer.start("debug_validate_sequences");
  if (verbose) {
    std::cout
        << "\n============================================================\n"
        << "=== Sequence Debug Validation: " << ID << " ===\n"
        << "============================================================\n";
  }

  int err = 0;

  auto seqSet = this->getSequences();
  std::vector<std::string> seqNames(seqSet.begin(), seqSet.end());
  std::sort(seqNames.begin(), seqNames.end());
  size_t numSeqs = seqNames.size();

  struct MismatchLog {
    size_t pos;
    char expected;
    char actual;
    std::string context_orig;
    std::string context_recon;
  };

  struct SequenceResult {
    std::string seqName;
    bool pass = true;
    size_t len_before = 0;
    size_t len_after = 0;
    int mismatch_count = 0;
    std::vector<MismatchLog> error_logs;
  };

  std::vector<SequenceResult> results(numSeqs);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, numSeqs),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (size_t idx = r.begin(); idx < r.end(); ++idx) {
                        const auto &seqName = seqNames[idx];
                        auto &res = results[idx];
                        res.seqName = seqName;

                        auto seq_after = this->reconstructSequence(seqName);
                        auto seq_before = manager->getSequence(seqName);

                        res.len_before = seq_before.length();
                        res.len_after = seq_after.length();
                        res.pass = true;

                        if (res.len_before != res.len_after) {
                          res.pass = false;
                        }

                        size_t min_len = std::min(res.len_before, res.len_after);
                        res.mismatch_count = 0;
                        const int MAX_MISMATCH_PRINT = 1;

                        for (size_t i = 0; i < min_len; ++i) {
                          if (seq_before[i] != seq_after[i]) {
                            res.pass = false;
                            if (res.mismatch_count < MAX_MISMATCH_PRINT) {
                              int ctx_start = std::max(0, (int)i - 5);
                              int ctx_end = std::min((int)min_len, (int)i + 6);
                              int ctx_len = ctx_end - ctx_start;

                              res.error_logs.push_back(
                                  {i, seq_before[i], seq_after[i],
                                   seq_before.substr(ctx_start, ctx_len),
                                   seq_after.substr(ctx_start, ctx_len)});
                            }
                            res.mismatch_count++;
                          }
                        }
                      }
                    });

  for (size_t idx = 0; idx < numSeqs; ++idx) {
    const auto &res = results[idx];
    std::cout << "Validate Sequence " << res.seqName << '\n';
    if (res.pass) {
      if (verbose) {
        std::cout << "  ✅ [PERFECT] Validation Passed! Reconstructed sequence "
                     "perfectly matches the raw sequence. (Len: "
                  << res.len_after << " bp)\n";
      }
    } else {
      err++;
      std::cerr
          << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch.\n";

      // 報告長度差異
      if (res.len_before != res.len_after) {
        std::cerr << "     ├─ Length mismatch: Original = " << res.len_before
                  << " bp, Reconstructed = " << res.len_after << " bp (Diff: "
                  << (long long)res.len_after - (long long)res.len_before << " bp)\n";
      } else {
        std::cerr << "     ├─ Length matches perfectly: " << res.len_before
                  << " bp\n";
      }

      // 報告內容差異
      if (res.mismatch_count > 0) {
        std::cerr
            << "     ├─ Total base mismatches found in overlapping region: "
            << res.mismatch_count << "\n";
        std::cerr << "     └─ Showing first " << res.error_logs.size()
                  << " mismatch(es):\n";

        for (const auto &log : res.error_logs) {
          std::cerr << "          ▶ Pos " << log.pos << ": Expected '"
                    << log.expected << "', but got '" << log.actual << "'\n";
          std::cerr << "            - Orig context : " << log.context_orig
                    << "\n";
          std::cerr << "            - Recon context: " << log.context_recon
                    << "\n";

          // 標示出錯誤位置的對齊標記 (箭頭 ^)
          int relative_pos = log.pos - std::max(0, (int)log.pos - 5);
          std::cerr << "                           "
                    << std::string(relative_pos, ' ') << "^\n";
        }
      } else if (res.len_before != res.len_after) {
        std::cerr << "     └─ Sequences match perfectly up to "
                  << std::min(res.len_before, res.len_after)
                  << " bp, but then one is abruptly truncated.\n";
      }

      // 🔍 詳細診斷：找出該序列在全圖所有 Block 內的 Segments，檢查是否有重疊或斷裂
      struct AuditSeg {
        int start, end;
        bool isRev;
        int copy;
        BlockID blkId;
        int consLen;
      };
      std::vector<AuditSeg> segList;
      for (const auto &blkPair : blocks) {
        auto blk = blkPair.second;
        auto it = blk->getSequences().find(res.seqName);
        if (it != blk->getSequences().end()) {
          for (const auto &sp : it->second.getSegments()) {
            segList.push_back({sp.second.getStart(), sp.second.getEnd(), sp.second.isReverse(),
                               sp.second.getCopyCount(), blk->getId(), (int)blk->getConsensus().length()});
          }
        }
      }
      std::sort(segList.begin(), segList.end(), [](const AuditSeg &a, const AuditSeg &b) {
        return a.start < b.start;
      });

      std::cerr << "     ├─ [SEGMENT AUDIT] Total segments in graph: " << segList.size() << "\n";
      int mismatch_pos = (res.error_logs.empty()) ? 0 : (int)res.error_logs[0].pos;
      for (size_t si = 0; si < segList.size(); ++si) {
        bool is_overlap = (si > 0 && segList[si].start < segList[si - 1].end);
        bool is_gap = (si > 0 && segList[si].start > segList[si - 1].end);
        bool near_mismatch = (std::abs(segList[si].start - mismatch_pos) < 10000 || std::abs(segList[si].end - mismatch_pos) < 10000);

        if (is_overlap || is_gap || near_mismatch) {
          std::cerr << "          " << (is_overlap ? "🚨 [OVERLAP] " : (is_gap ? "⚠️ [GAP] " : "   [SEG] "))
                    << "Block " << segList[si].blkId << " (ConsLen " << segList[si].consLen << ", Copy " << segList[si].copy << ", Strand " << (segList[si].isRev ? "-" : "+") << "): "
                    << "[" << segList[si].start << ", " << segList[si].end << ") len=" << (segList[si].end - segList[si].start) << " bp\n";
          if (is_overlap) {
            std::cerr << "               ↳ Overlaps with previous Block " << segList[si - 1].blkId << " by "
                      << (segList[si - 1].end - segList[si].start) << " bp!\n";
          }
        }
      }
    }
  }

  if (err > 0) {
    std::cerr << "\n❌ [FATAL] Sequence validation failed for " << err
              << " sequences. Exiting.\n";
    global_timer.stop("debug_validate_sequences");
    exit(1);
  }
  std::cout << "  ✅ [PERFECT] All Sequences Validation Passed!\n";
  global_timer.stop("debug_validate_sequences");
}

void BlockSet::debugValidateCopies(bool verbose) {
  global_timer.start("debug_validate_copies");
  bool has_error = false;
  int total_error_blocks = 0;
  int total_error_vblocks = 0;

  for (auto& weak_blk : this->getAllBlocks()) {
    auto blk = weak_blk.lock();
    if (!blk) continue;

    // 依據 copyCount 整理出每個 Copy 底下包含的 Sequence 與其 Segments
    std::map<int, std::map<std::string, std::vector<const Segment*>>> copy_seq_segs;

    for (const auto& seqPair : blk->getSequences()) {
      const std::string& seqName = seqPair.first;
      for (const auto& segPair : seqPair.second.getSegments()) {
        const Segment& seg = segPair.second;
        copy_seq_segs[seg.getCopyCount()][seqName].push_back(&seg);
      }
    }

    bool block_has_error = false;
    for (const auto& copyPair : copy_seq_segs) {
      int copy = copyPair.first;
      for (const auto& seqSegsPair : copyPair.second) {
        const std::string& seqName = seqSegsPair.first;
        const auto& segs = seqSegsPair.second;

        // 🚨 檢查：同一個 (Block, Copy) 底下，同一條 Sequence 絕不能包含 >= 2 個 Segments
        if (segs.size() > 1) {
          has_error = true;
          if (!block_has_error) {
            total_error_blocks++;
            block_has_error = true;
            std::cerr << "\n❌ [COPY VALIDATION ERROR] Block " << blk->getId() 
                      << " (Consensus Len: " << blk->getConsensus().length() << " bp)\n";
          }
          total_error_vblocks++;
          std::cerr << "   ├─ VBlock (Block_" << blk->getId() << ", Copy " << copy << ")\n"
                    << "   │  └─ Sequence '" << seqName << "' has " << segs.size() 
                    << " duplicate segments sharing copyCount=" << copy << ":\n";
          for (size_t sIdx = 0; sIdx < segs.size(); ++sIdx) {
            const Segment* s = segs[sIdx];
            std::cerr << "   │     [" << sIdx + 1 << "] start=" << s->getStart() 
                      << ", end=" << s->getEnd()
                      << ", len=" << std::abs(s->getEnd() - s->getStart()) << " bp"
                      << ", strand=" << (s->isReverse() ? "-" : "+") << "\n";
          }
        }
      }
    }
  }

  if (has_error) {
    std::cerr << "\n💥 [FATAL] debugValidateCopies failed! Found " << total_error_vblocks 
              << " invalid VBlocks with duplicate segments across " << total_error_blocks 
              << " blocks. Terminating execution.\n";
    global_timer.stop("debug_validate_copies");
    exit(1);
  } else {
    if (verbose) {
      std::cout << "  ✅ [PERFECT] debugValidateCopies passed! All VBlocks contain at most one segment per sequence.\n";
    }
  }
  global_timer.stop("debug_validate_copies");
}

void BlockSet::debugValidateLinearizedBlocks(bool verbose) {
  global_timer.start("debug_validate_linearized_blocks");
  // 🌟 1. 改讀取最新的虛擬積木 1D 骨幹
  auto vblocks = this->linear_block_cache;
  if (vblocks.empty()) {
    global_timer.stop("debug_validate_linearized_blocks");
    return;
  }

  // 建立 VBlockID 到 Index 的查表 (這裡用 std::map 比較方便)
  std::map<VBlockID, int> vid_to_index;
  int N = (int)vblocks.size();
  for (int i = 0; i < N; ++i) {
    vid_to_index[vblocks[i]] = i;
  }

  bool all_valid = true;

  // 🌟 2. 收集每一條 Sequence "真實的物理路徑"
  struct SegInfo {
    int start_coord; // 物理座標
    VBlockID vid;    // 所屬的 VBlock
    bool is_reverse; // 是否為反鏈 (inversion)
  };
  std::unordered_map<std::string, std::vector<SegInfo>> true_seq_paths;

  auto allBlocks = this->getAllBlocks();
  for (auto &weak_blk : allBlocks) {
    auto blk = weak_blk.lock();
    if (!blk)
      continue;
    BlockID b_id = blk->getId();

    for (auto &[seqName, seqData] : blk->getSequences()) {
      for (auto &segPair : seqData.getSegments()) {
        int copy = segPair.second.getCopyCount(); // 取得 Copy
        VBlockID vid = {b_id, copy};

        // 確保這個 VBlock 真的有被收錄到骨幹裡
        if (vid_to_index.find(vid) != vid_to_index.end()) {
          true_seq_paths[seqName].push_back(
              {std::min(segPair.second.getStart(), segPair.second.getEnd()),
               vid, segPair.second.isReverse()});
        }
      }
    }
  }

  // 🌟 3. 排序並嚴格檢查 VBlock 的拓撲順序
  for (auto &[seqName, path] : true_seq_paths) {
    if (path.size() <= 1)
      continue;

    // 依照物理座標由小到大排序，還原生物學上真實的走訪順序
    std::sort(path.begin(), path.end(), [](const SegInfo &a, const SegInfo &b) {
      return a.start_coord < b.start_coord;
    });

    for (size_t k = 0; k < path.size() - 1; ++k) {
      int curr_idx = vid_to_index[path[k].vid];
      int next_idx = vid_to_index[path[k + 1].vid];
      bool is_inversion_transition =
          path[k].is_reverse && path[k + 1].is_reverse;

      // ✅ 正常情況
      if (!is_inversion_transition) {
        // 非 inversion 轉換：在 1D 骨幹中嚴格遞增
        if (next_idx > curr_idx) {
          continue;
        }
        // 🔄 環狀情況：從最後面跳回最前面
        else if (curr_idx == N - 1 && next_idx == 0) {
          if (verbose)
            std::cout << "  [Info] Seq " << seqName
                      << " circular wrap-around detected.\n";
          continue;
        }
      } else {
        // inversion 內部轉換：在 1D 骨幹中嚴格遞減 (因為物理順序與轉錄順序相反)
        if (next_idx < curr_idx) {
          continue;
        }
        // 🔄 環狀情況 (反向)：從最前面跳回最後面
        else if (curr_idx == 0 && next_idx == N - 1) {
          if (verbose)
            std::cout << "  [Info] Seq " << seqName
                      << " circular wrap-around (reverse) detected.\n";
          continue;
        }
      }

      // ❌ 異常情況：Self-loop 或是 拓撲倒退！
      {
        std::cerr << "  🚨 [ERROR] Topological violation in Seq " << seqName
                  << "\n"
                  << "      -> VBlock {B:" << path[k].vid.first
                  << ", C:" << path[k].vid.second << "} (index " << curr_idx
                  << ")\n"
                  << "      -> VBlock {B:" << path[k + 1].vid.first
                  << ", C:" << path[k + 1].vid.second << "} (index " << next_idx
                  << ")\n";

        if (curr_idx == next_idx) {
          std::cerr << "      🔥 Reason: Self-loop! Two segments of the same "
                       "sequence share the exact same VBlock. (Copy "
                       "distribution failed!)\n";
        } else {
          std::cerr
              << "      🔥 Reason: Backbone order reversed! DFS topological "
                 "sort failed to maintain physical direction (is_inversion="
              << is_inversion_transition << ").\n";
        }
        all_valid = false;
      }
    }
  }

  if (all_valid) {
    if (verbose)
      std::cout << "  ✅ [Validation] All sequences perfectly align with the "
                   "1D VBlock backbone.\n";
  } else {
    std::cerr
        << "  ❌ [Validation] Found topological violations in VBlock graph!\n";
  }
  global_timer.stop("debug_validate_linearized_blocks");
}

void BlockSet::debugValidateBlocks(bool verbose) {
  global_timer.start("debug_validate_blocks");
  auto linearBlocks = this->getLinearizeBlocks();
  if (linearBlocks.empty()) {
    std::cout << "\n  [Warning] Graph is empty. No blocks to validate.\n";
    global_timer.stop("debug_validate_blocks");
    return;
  }

  // 稍微加寬分隔線以容納新增的欄位
  std::cout << "\n============================================================="
               "=======================================================\n"
            << "=== Linearized Block Distribution & Quality Report ===\n"
            << "==============================================================="
               "=====================================================\n";

  // 定義統計級距 (Bins)
  struct SizeBin {
    std::string label;
    int min_val;
    int max_val;
    int block_count = 0;
    double sum_identity = 0.0;
    int sum_sequences = 0;
    int sum_segments = 0; // 🌟 新增：紀錄該級距的總 Segments
    uint64_t sum_snvs = 0;
    uint64_t sum_gaps = 0;
  };

  std::vector<SizeBin> bins = {{"1-10", 1, 10},
                               {"11-50", 11, 50},
                               {"51-100", 51, 100},
                               {"101-200", 101, 200},
                               {"201-500", 201, 500},
                               {"501-1000", 501, 1000},
                               {"1001-10000", 1001, 10000},
                               {"10000+", 10001, INT32_MAX}};

  if (verbose) {
    std::cout << "[Block-level Details (Linear Order)]\n";
  }

  std::unordered_set<int>
      processed_blocks; // 🌟 核心修改 2：用來記錄已經計算過的 BlockID

  for (VBlockID vbid : linearBlocks) {
    int blkId = vbid.first;

    // 🌟 核心修改 2：如果這個 BlockID 已經被算過，就跳過 (只計算實體積木一次)
    if (processed_blocks.count(blkId))
      continue;
    processed_blocks.insert(blkId);

    auto blk = this->getBlock(blkId);
    if (!blk)
      continue;

    int consLen = blk->getConsensus().length();
    int seqCount =
        blk->getSequences().size(); // 🌟 核心修改 1：移除了 FamID 邏輯

    uint64_t snv_count = 0;
    uint64_t total_gap_len = 0;
    int segCount = 0;

    for (auto &[seqName, seqData] : blk->getSequences()) {
      for (auto &[segId, segNode] : seqData.getSegments()) {
        segCount++; // 🌟 核心修改 3：計算實際的 Segments 數量
        for (auto &var : segNode.getVariants()) {
          if (var.getType() == VariantType::SNV) {
            snv_count += 1;
          } else if (var.getType() == VariantType::GAP) {
            total_gap_len += (var.getEnd() - var.getStart());
          }
        }
      }
    }

    // 計算總變異長度供 Identity 使用
    uint64_t varLen = snv_count + total_gap_len;

    double identity = 1.0; // 預設 1.0 (包含 Singleton)
    if (segCount > 1 && consLen > 0) {
      uint64_t denom = (uint64_t)segCount * consLen;
      identity = 1.0 - ((double)varLen / denom);
    }

    // 單行詳細輸出
    bool low_identity = (identity < 0.995);
    if (verbose && low_identity) {
      std::cout << "  ├─ Block ID: " << std::setw(6) << std::left << blkId
                << " | Len: " << std::setw(7) << consLen
                << " | Seqs: " << std::setw(3) << seqCount
                << " | Segs: " << std::setw(4)
                << segCount // 🌟 新增：印出 Segments 數量
                << " | Identity: " << std::fixed << std::setprecision(4)
                << identity << " | SNVs: " << std::setw(5) << snv_count
                << " | TotalGapLen: " << std::setw(6) << total_gap_len;

      if (segCount == 1)
        std::cout << " (Singleton)";
      std::cout << "\n";
    }

    // 分類到對應的級距中
    for (auto &bin : bins) {
      if (consLen >= bin.min_val && consLen <= bin.max_val) {
        bin.block_count++;
        bin.sum_sequences += seqCount;
        bin.sum_segments += segCount; // 🌟 累積 Segments
        bin.sum_identity += identity;
        bin.sum_snvs += snv_count;
        bin.sum_gaps += total_gap_len;
        break;
      }
    }
  }

  // ==========================================
  // 印出統計結果表格
  // ==========================================
  std::cout << "\n[Block Size Distribution & Statistics]\n";
  std::cout << std::string(116, '-') << "\n"; // 配合新增欄位拉寬
  std::cout << std::setw(15) << std::left << "Size Range"
            << "| " << std::setw(10) << "Blocks"
            << "| " << std::setw(10) << "Avg Seqs"
            << "| " << std::setw(10) << "Avg Segs" // 🌟 表格新增欄位
            << "| " << std::setw(15) << "Avg Identity"
            << "| " << std::setw(15) << "Total SNVs"
            << "| " << std::setw(15) << "Total GapLen"
            << "\n";
  std::cout << std::string(116, '-') << "\n";

  int total_blocks = 0;
  uint64_t grand_total_snvs = 0;
  uint64_t grand_total_gaps = 0;

  for (const auto &bin : bins) {
    if (bin.block_count == 0)
      continue;

    double avg_seqs = (double)bin.sum_sequences / bin.block_count;
    double avg_segs =
        (double)bin.sum_segments / bin.block_count; // 🌟 計算平均 Segments
    double avg_identity = bin.sum_identity / bin.block_count;

    std::cout << std::setw(15) << std::left << bin.label << "| "
              << std::setw(10) << bin.block_count << "| " << std::setw(10)
              << std::fixed << std::setprecision(2) << avg_seqs << "| "
              << std::setw(10) << std::fixed << std::setprecision(2)
              << avg_segs // 🌟 印出
              << "| " << std::setw(15) << std::fixed << std::setprecision(4)
              << avg_identity << "| " << std::setw(15) << bin.sum_snvs << "| "
              << std::setw(15) << bin.sum_gaps << "\n";

    total_blocks += bin.block_count;
    grand_total_snvs += bin.sum_snvs;
    grand_total_gaps += bin.sum_gaps;
  }

  // 全域最終數據 Summary
  std::cout << std::string(116, '-') << "\n";
  std::cout << "Total Unique Blocks Evaluated : " << total_blocks
            << "\n"; // 🌟 改名提示這是 Unique
  std::cout << "Total SNVs in Graph           : " << grand_total_snvs << "\n";
  std::cout << "Total Gap Length in Graph     : " << grand_total_gaps << "\n";
  std::cout << "==============================================================="
               "=====================================================\n\n";
  global_timer.stop("debug_validate_blocks");
}

/*
void BlockSet::debugValidateBubble(bool verbose) {
    std::cout <<
"\n============================================================\n"; std::cout <<
"=== BlockSet: Bubble Topology Validator (Raw Consensus) ===\n"; std::cout <<
"============================================================\n";

    // 尋找氣泡：使用 Genomic 左右兩端的 Hub Block ID 作為 Key
    std::map<std::pair<Block::ID, Block::ID>, std::vector<Block::ID>>
bubble_groups;

    for (auto blk : this->getAllBlocks()) {
        if (!blk) continue;

        Block::ID left_id = 0, right_id = 0;
        bool is_consistent = true;
        bool first = true;

        for (auto& seqPair : blk->getSequences()) {
            for (auto& segPairInner : seqPair.second.getSegments()) {
                Segment& seg = segPairInner.second;

                // 根據方向性，精準取得 Genomic 上的左邊與右邊積木
                auto l_ptr = !seg.isReverse() ? seg.getPrevBlock().lock() :
seg.getNextBlock().lock(); auto r_ptr = !seg.isReverse() ?
seg.getNextBlock().lock() : seg.getPrevBlock().lock();

                Block::ID cur_l = l_ptr ? l_ptr->getId() : 0;
                Block::ID cur_r = r_ptr ? r_ptr->getId() : 0;

                if (first) {
                    left_id = cur_l; right_id = cur_r; first = false;
                } else {
                    if (cur_l != left_id || cur_r != right_id) {
                        is_consistent = false; break;
                    }
                }
            }
            if (!is_consistent) break;
        }

        // 如果這個積木完美夾在兩個 Hub 之間，就把它加入該氣泡群組
        if (is_consistent && left_id != 0 && right_id != 0 && left_id !=
right_id) { Block::ID min_id = std::min(left_id, right_id); Block::ID max_id =
std::max(left_id, right_id); bubble_groups[{min_id,
max_id}].push_back(blk->getId());
        }
    }

    int singletonBubbleLen = 0;

    int bubble_count = 0;
    for (auto& kv : bubble_groups) {
        if (kv.second.size() < 2) continue; // 只有一條路徑，不是氣泡

        bubble_count++;
        Block::ID left_id = kv.first.first;
        Block::ID right_id = kv.first.second;
        auto leftBlk = this->getBlock(left_id);
        auto rightBlk = this->getBlock(right_id);

        int left_seg = 0, right_seg = 0;
        for (auto& seq: leftBlk->getSequences())
            for (auto& seg: seq.second.getSegments()) left_seg++;

        for (auto& seq: rightBlk->getSequences())
            for (auto& seg: seq.second.getSegments()) right_seg++;

        std::cout << "\n[Bubble #" << bubble_count << "]\n";
        std::cout << "  ├─ Flanking Hub A (ID: " << left_id << ") | Len: "
                  << (leftBlk ? leftBlk->getConsensus().length() : 0) << " bp |
Seqs: "
                  << (leftBlk ? left_seg : 0) << "\n";
        std::cout << "  ├─ Flanking Hub B (ID: " << right_id << ") | Len: "
                  << (rightBlk ? rightBlk->getConsensus().length() : 0) << " bp
| Seqs: "
                  << (rightBlk ? right_seg : 0) << "\n";
        std::cout << "  └─ Branches (" << kv.second.size() << " parallel
paths):\n\n";

        for (Block::ID bId : kv.second) {
            auto bBlk = this->getBlock(bId);
            if (!bBlk) continue;

            int bSeg = 0;
            for (auto& seq: bBlk->getSequences())
                for (auto& seg: seq.second.getSegments()) bSeg++;

            // 判斷方向：如果這條分支的序列在 Genomic
上是反向，我們印出來前先把它轉正 bool is_rev = false; if
(!bBlk->getSequences().empty() &&
!bBlk->getSequences().begin()->second.getSegments().empty()) { is_rev =
bBlk->getSequences().begin()->second.getSegments().begin()->second.isReverse();
            }

            std::string seqStr = bBlk->getConsensus();
            if (is_rev) {
                std::reverse(seqStr.begin(), seqStr.end());
                for (char& c : seqStr) {
                    switch (c) {
                        case 'A': c='T'; break; case 'T': c='A'; break;
                        case 'C': c='G'; break; case 'G': c='C'; break;
                        case 'a': c='t'; break; case 't': c='a'; break;
                        case 'c': c='g'; break; case 'g': c='c'; break;
                    }
                }
            }

            // 直接印出這塊分支的 ID、長度、包含序列數，以及最純粹的 Consensus
            std::cout << "       [Block " << std::left << std::setw(6) << bId <<
"] "
                      << "(Len: " << std::setw(4) << seqStr.length() << " bp,
Seqs: " << std::setw(2) << bSeg << ")\n"
                      << "       Seq: " << seqStr << "\n\n";

            if (bSeg == 1) singletonBubbleLen+= seqStr.length();
        }
        std::cout <<
"------------------------------------------------------------\n";
    }

    if (bubble_count == 0) {
        std::cout << "  -> No simple bubbles found!\n";
    } else {
        std::cout << "  -> Total " << bubble_count << " bubbles identified and
validated.\n"; std::cout << "  -> Total singleton bubble length: " <<
singletonBubbleLen << " bp\n";
    }
    std::cout <<
"============================================================\n\n";
}
*/