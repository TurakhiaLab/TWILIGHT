
#include "block_set.hpp"
#include "block_manager.hpp"
#include "phylogeny.hpp"
#include "timer.hpp"
#include "type.hpp"
#include "cigar_util.hpp"
#include "coordinate_manager.hpp"

#include <boost/filesystem.hpp>
#include <string>
#include <queue>
#include <limits>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>


void Alignment::updateEnergy(BlockSet *refBlockSet, BlockSet *qryBlockSet, double beta) {
  double aln_len = static_cast<double>(alnLength);
  double aln_score = alnScore;
  energy = -aln_score;
  // energy = -aln_len + beta * total_variation_count;
}

void BlockSet::buildCaches() {
  // Ancestral
  ancestral_block_cache.clear();
  for (auto &vblkID : linear_block_cache) {
    auto blk = this->getBlock(vblkID.first);
    if (blk->isDistant(vblkID.second))
      continue;
    ancestral_block_cache.push_back(vblkID);
  }

  // Representative
  for (auto vblkID : ancestral_block_cache) {
    auto blk = this->getBlock(vblkID.first);
    if (blk)
      ancestral_seq_cache += blk->getConsensus().getConsensusString();
  }

  // Offset
  ancestral_offset_cache.clear();
  int current_offset = 0;
  for (const auto &vblkID : ancestral_block_cache) {
    auto blk = this->getBlock(vblkID.first);
    if (blk) {
      ancestral_offset_cache[current_offset] = vblkID;
      current_offset += blk->getConsensus().length();
    }
  }

  is_cached = true;
}

void BlockSet::rebuildAllPointers() {
  bool debug_mode = false; // 🌟 開關：控制 Debug 訊息的輸出

  if (debug_mode)
    std::cout << "\n[RebuildPointers] === Start Rebuilding All Pointers for Node '"
              << this->getId() << "' (tree_ptr=" << (tree_ptr ? "SET" : "NULL") << ") ===\n";

  auto allBlocks = this->getAllBlocks();
  if (allBlocks.empty()) {
    if (debug_mode)
      std::cout << "[RebuildPointers] No blocks found. Exiting.\n";
    return;
  }

  // Reset distant status for all blocks and copies before rebuilding the DAG
  for (auto &blk : allBlocks) {
    auto sharedBlk = blk.lock();
    if (!sharedBlk)
      continue;
    int maxCopy = sharedBlk->getMaxCopy();
    for (int c = 0; c <= maxCopy; ++c) {
      sharedBlk->setDistant(false, c);
    }
  }

  // ==========================================
  // 1. 清理舊指標並收集所有 Segments (加入 Copy 資訊)
  // ==========================================
  struct SegRef {
    Segment *seg;
    std::shared_ptr<Block> blk;
    int copy;
  };
  std::map<std::string, std::vector<SegRef>> seqTracks;
  int total_segments = 0;

  for (auto &blk : allBlocks) {
    auto sharedBlk = blk.lock();
    if (!sharedBlk)
      continue;

    sharedBlk->clearPrevBlocksMap();
    sharedBlk->clearNextBlocksMap();

    for (auto &seqPair : sharedBlk->getSequences()) {
      for (auto &segPairInner : seqPair.second.getSegments()) {
        Segment *currentSeg = &segPairInner.second;
        int currentCopy = currentSeg->getCopyCount();

        currentSeg->setPrevBlock(std::shared_ptr<Block>(nullptr));
        currentSeg->setNextBlock(std::shared_ptr<Block>(nullptr));

        seqTracks[seqPair.first].push_back(
            {currentSeg, sharedBlk, currentCopy});
        total_segments++;
      }
    }
  }

  if (debug_mode) {
    std::cout << "[RebuildPointers] Step 1: Extracted " << total_segments
              << " segments across " << seqTracks.size() << " sequences.\n";
  }

  // ==========================================
  // 2. 重新接線 Segment 層級的指標 (TBB 平行加速)
  // ==========================================
  std::vector<std::vector<SegRef> *> trackPtrs;
  trackPtrs.reserve(seqTracks.size());
  for (auto &trackPair : seqTracks) {
    trackPtrs.push_back(&trackPair.second);
  }

  std::mutex debug_mutex; // 用於 TBB 中的 Thread-safe print

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, trackPtrs.size()),
      [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          auto &track = *(trackPtrs[i]);

          // 🌟 依照絕對物理座標排序 (建立由左至右的骨幹，完全不考慮 strand)
          std::sort(track.begin(), track.end(),
                    [](const SegRef &a, const SegRef &b) {
                      return std::min(a.seg->getStart(), a.seg->getEnd()) <
                             std::min(b.seg->getStart(), b.seg->getEnd());
                    });

          // 進行雙向綁定
          for (size_t j = 1; j < track.size(); ++j) {
            auto &prevRef = track[j - 1];
            auto &currRef = track[j];

            if (!prevRef.seg->isReverse())
              prevRef.seg->setNextBlock(currRef.blk);
            else
              prevRef.seg->setPrevBlock(currRef.blk);

            if (!currRef.seg->isReverse())
              currRef.seg->setPrevBlock(prevRef.blk);
            else
              currRef.seg->setNextBlock(prevRef.blk);
          }

          // (可選) 打印每條 sequence 的軌跡狀態
          if (debug_mode && track.size() > 0) {
            std::lock_guard<std::mutex> lock(debug_mutex);
            // std::cout << "  -> Seq Track Sorted: " << track.size() << "
            // segments.\n";
          }
        }
      });

  if (debug_mode)
    std::cout << "[RebuildPointers] Step 2: Segment pointers re-linked "
                 "successfully.\n";

  // ==========================================
  // 3. 建立虛擬積木 (VBlock) 層級的「嚴格線性」指標
  // ==========================================

  std::vector<VBlockID> uniqueNodes;
  for (auto &trackPair : seqTracks) {
    for (auto &ref : trackPair.second) {
      uniqueNodes.push_back({ref.blk->getId(), ref.copy});
    }
  }
  std::sort(uniqueNodes.begin(), uniqueNodes.end());
  uniqueNodes.erase(std::unique(uniqueNodes.begin(), uniqueNodes.end()),
                    uniqueNodes.end());

  size_t N = uniqueNodes.size();
  if (debug_mode)
    std::cout << "[RebuildPointers] Step 3: Generating DAG for " << N
              << " unique VBlocks.\n";

  auto getIdx = [&](const VBlockID &vid) -> size_t {
    auto it = std::lower_bound(uniqueNodes.begin(), uniqueNodes.end(), vid);
    return std::distance(uniqueNodes.begin(), it);
  };

  std::vector<std::vector<size_t>> adj(N);
  std::vector<int> inDegree(N, 0);
  std::vector<int> minCoord(N, std::numeric_limits<int>::max());
  std::vector<std::shared_ptr<Block>> idToBlock(N);
  std::vector<int> state(N, 0); // 0=未走訪, 1=走訪中, 2=已走完

  for (auto &trackPair : seqTracks) {
    for (auto &ref : trackPair.second) {
      size_t u = getIdx({ref.blk->getId(), ref.copy});
      idToBlock[u] = ref.blk;
      int currentMin = std::min(ref.seg->getStart(), ref.seg->getEnd());
      minCoord[u] = std::min(minCoord[u], currentMin);
    }
  }

  // 建立有向圖 (1D 物理相鄰關係)
  for (auto &trackPair : seqTracks) {
    auto &track = trackPair.second;
    for (size_t j = 1; j < track.size(); ++j) {
      size_t u = getIdx({track[j - 1].blk->getId(), track[j - 1].copy});
      size_t v = getIdx({track[j].blk->getId(), track[j].copy});

      if (u != v) {
        if (std::find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
          adj[u].push_back(v);
          inDegree[v]++;
        }
      }
    }
  }

  std::vector<size_t> nodes(N);
  for (size_t i = 0; i < N; ++i)
    nodes[i] = i;

  // 起點排序：入度小優先 -> 座標小優先
  std::sort(nodes.begin(), nodes.end(), [&](size_t a, size_t b) {
    if (inDegree[a] != inDegree[b])
      return inDegree[a] < inDegree[b];
    if (minCoord[a] != minCoord[b])
      return minCoord[a] < minCoord[b];
    return a < b;
  });

  // 分支排序：物理座標小優先
  for (size_t i = 0; i < N; ++i) {
    std::sort(adj[i].begin(), adj[i].end(), [&](size_t a, size_t b) {
      if (minCoord[a] != minCoord[b])
        return minCoord[a] < minCoord[b];
      return a < b;
    });
  }

  std::vector<size_t> postOrder;
  postOrder.reserve(N);

  // 🌟 極速 DFS 遞迴引擎 + Cycle Detection (圖環偵測)
  bool cycle_detected = false;
  std::function<void(size_t)> dfs = [&](size_t u) {
    state[u] = 1; // 標記為走訪中
    for (size_t v : adj[u]) {
      if (state[v] == 0) {
        dfs(v);
      } else if (state[v] == 1) {
        // 🚨 偵測到 Cycle (Back-edge)
        cycle_detected = true;
        if (debug_mode) {
          std::cerr << "[Warning] Cycle detected in Block Graph: VBlock("
                    << uniqueNodes[u].first << ", cp:" << uniqueNodes[u].second
                    << ") -> VBlock(" << uniqueNodes[v].first
                    << ", cp:" << uniqueNodes[v].second << ")\n";
        }
      }
    }
    state[u] = 2; // 標記為走訪完成
    postOrder.push_back(u);
  };

  for (size_t u : nodes) {
    if (state[u] == 0)
      dfs(u);
  }

  std::reverse(postOrder.begin(), postOrder.end());

  if (debug_mode && cycle_detected) {
    std::cout << "[RebuildPointers] Notice: Graph cyclic edges were ignored to "
                 "enforce linearity.\n";
  }

  // ==========================================
  // 4. 將結果寫回實體 Block 並快取 Linearized Backbone
  // ==========================================

  this->linear_block_cache.clear();
  this->linear_block_cache.reserve(postOrder.size());

  if (debug_mode)
    std::cout
        << "[RebuildPointers] Step 4: Constructing Linear Block Order...\n";

  for (size_t i = 0; i < postOrder.size(); ++i) {
    size_t currIdx = postOrder[i];
    VBlockID currVId = uniqueNodes[currIdx];
    auto currBlk = idToBlock[currIdx];
    int currCopy = currVId.second;

    this->linear_block_cache.push_back(currVId);

    if (i > 0) {
      size_t prevIdx = postOrder[i - 1];
      currBlk->setPrevBlock(currCopy, idToBlock[prevIdx]);
    } else {
      currBlk->setPrevBlock(currCopy, std::shared_ptr<Block>(nullptr));
    }

    if (i < postOrder.size() - 1) {
      size_t nextIdx = postOrder[i + 1];
      currBlk->setNextBlock(currCopy, idToBlock[nextIdx]);
    } else {
      currBlk->setNextBlock(currCopy, std::shared_ptr<Block>(nullptr));
    }

    // (可選) 打印最後排出來的線性骨幹順序
    if (debug_mode) {
      std::cout << "  - VBlock(" << currVId.first << ", cp:" << currCopy << ")";
      if (i < postOrder.size() - 1)
        std::cout << " -> ";
      if ((i + 1) % 5 == 0)
        std::cout << "\n";
    }
  }

  // if (tree_ptr) {
  //   this->setDistantBlocks(*tree_ptr);
  // }
  this->buildCaches();
}

bool BlockSet::detectVBlockCycle(const CoordinateManager &coordMgr, bool verbose) const {
  const auto &refIntervals = coordMgr.getRefIntervals();
  const auto &qryIntervals = coordMgr.getQryIntervals();

  std::vector<VBlockID> uniqueVBlocks;
  uniqueVBlocks.reserve(refIntervals.size() + qryIntervals.size());

  for (const auto &kv : refIntervals) {
    uniqueVBlocks.push_back({kv.second.blkId, coordMgr.getCopyId(kv.first, true)});
  }
  for (const auto &kv : qryIntervals) {
    uniqueVBlocks.push_back({kv.second.blkId, coordMgr.getCopyId(kv.first, false)});
  }

  std::sort(uniqueVBlocks.begin(), uniqueVBlocks.end());
  uniqueVBlocks.erase(std::unique(uniqueVBlocks.begin(), uniqueVBlocks.end()), uniqueVBlocks.end());

  size_t N = uniqueVBlocks.size();
  if (N == 0) return false;

  auto getIdx = [&](const VBlockID &vid) -> size_t {
    auto it = std::lower_bound(uniqueVBlocks.begin(), uniqueVBlocks.end(), vid);
    return std::distance(uniqueVBlocks.begin(), it);
  };

  std::vector<std::vector<size_t>> adj(N);

  if (!refIntervals.empty()) {
    auto prevIt = refIntervals.begin();
    auto currIt = std::next(prevIt);
    while (currIt != refIntervals.end()) {
      VBlockID uV = {prevIt->second.blkId, coordMgr.getCopyId(prevIt->first, true)};
      VBlockID vV = {currIt->second.blkId, coordMgr.getCopyId(currIt->first, true)};
      size_t u = getIdx(uV);
      size_t v = getIdx(vV);
      if (u != v && std::find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
        adj[u].push_back(v);
      }
      prevIt = currIt;
      ++currIt;
    }
  }

  if (!qryIntervals.empty()) {
    auto prevIt = qryIntervals.begin();
    auto currIt = std::next(prevIt);
    while (currIt != qryIntervals.end()) {
      VBlockID uV = {prevIt->second.blkId, coordMgr.getCopyId(prevIt->first, false)};
      VBlockID vV = {currIt->second.blkId, coordMgr.getCopyId(currIt->first, false)};
      size_t u = getIdx(uV);
      size_t v = getIdx(vV);
      if (u != v && std::find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
        adj[u].push_back(v);
      }
      prevIt = currIt;
      ++currIt;
    }
  }

  std::vector<int> state(N, 0); // 0=unvisited, 1=visiting, 2=visited
  bool cycle_detected = false;

  std::function<void(size_t)> dfs = [&](size_t u) {
    state[u] = 1;
    for (size_t v : adj[u]) {
      if (state[v] == 0) {
        dfs(v);
      } else if (state[v] == 1) {
        cycle_detected = true;
        if (verbose) {
          std::cerr << "🚨 [VBlock Cycle Detected] Cycle edge: VBlock("
                    << uniqueVBlocks[u].first << ", cp:" << uniqueVBlocks[u].second
                    << ") -> VBlock(" << uniqueVBlocks[v].first
                    << ", cp:" << uniqueVBlocks[v].second << ")\n";
        }
      }
    }
    state[u] = 2;
  };

  for (size_t i = 0; i < N; ++i) {
    if (state[i] == 0) dfs(i);
  }

  return cycle_detected;
}

void BlockSet::rebuildLinearGraph(const CoordinateManager &coordMgr) {
  const auto &refIntervals = coordMgr.getRefIntervals();
  const auto &qryIntervals = coordMgr.getQryIntervals();

  std::vector<VBlockID> refList;
  refList.reserve(refIntervals.size());
  for (const auto &kv : refIntervals) {
    refList.push_back({kv.second.blkId, coordMgr.getCopyId(kv.first, true)});
  }

  std::vector<VBlockID> qryList;
  qryList.reserve(qryIntervals.size());
  for (const auto &kv : qryIntervals) {
    qryList.push_back({kv.second.blkId, coordMgr.getCopyId(kv.first, false)});
  }

  if (refList.empty() && qryList.empty()) {
    this->linear_block_cache.clear();
    this->buildCaches();
    return;
  }

  // 1. 識別同時間存在於 Ref 與 Qry 的 Mode 1 共享 VBlocks (inBothSet)
  std::set<VBlockID> inBothSet;
  for (const auto &vid : refList) {
    if (coordMgr.isVBlockInBoth(vid)) {
      inBothSet.insert(vid);
    }
  }
  // 萬一 coordMgr 標籤無資訊，做集合交集保險補強
  std::set<VBlockID> qryVSet(qryList.begin(), qryList.end());
  for (const auto &vid : refList) {
    if (qryVSet.count(vid)) {
      inBothSet.insert(vid);
    }
  }

  // 2. 雙指標 (pRef, pQry) 依序動態交錯合併 Ref/Qry 1D 順序
  std::vector<VBlockID> finalLinearVBlocks;
  finalLinearVBlocks.reserve(refList.size() + qryList.size());

  size_t pRef = 0;
  size_t pQry = 0;
  std::set<VBlockID> processed;

  while (pRef < refList.size() || pQry < qryList.size()) {
    if (pRef >= refList.size()) {
      if (processed.insert(qryList[pQry]).second) {
        finalLinearVBlocks.push_back(qryList[pQry]);
      }
      pQry++;
      continue;
    }

    if (pQry >= qryList.size()) {
      if (processed.insert(refList[pRef]).second) {
        finalLinearVBlocks.push_back(refList[pRef]);
      }
      pRef++;
      continue;
    }

    const auto &rVid = refList[pRef];
    const auto &qVid = qryList[pQry];

    if (processed.count(rVid)) {
      pRef++;
      continue;
    }

    if (processed.count(qVid)) {
      pQry++;
      continue;
    }

    // 兩邊剛好走到同一個 is_both / 共享 VBlock
    if (rVid == qVid) {
      finalLinearVBlocks.push_back(rVid);
      processed.insert(rVid);
      pRef++;
      pQry++;
      continue;
    }

    bool rIsShared = (inBothSet.count(rVid) > 0);
    bool qIsShared = (inBothSet.count(qVid) > 0);

    if (!rIsShared) {
      // Ref 側目前是 unique VBlock，直接 append 並推進 Ref
      finalLinearVBlocks.push_back(rVid);
      processed.insert(rVid);
      pRef++;
    } else if (!qIsShared) {
      // Qry 側目前是 unique VBlock，直接 append 並推進 Qry
      finalLinearVBlocks.push_back(qVid);
      processed.insert(qVid);
      pQry++;
    } else {
      // 兩邊都是 shared 錨點但點不一樣，尋找哪一邊的指標離對方錨點較近
      size_t rFindQ = refList.size();
      for (size_t k = pRef; k < refList.size(); ++k) {
        if (refList[k] == qVid) { rFindQ = k; break; }
      }

      size_t qFindR = qryList.size();
      for (size_t k = pQry; k < qryList.size(); ++k) {
        if (qryList[k] == rVid) { qFindR = k; break; }
      }

      if (qFindR < qryList.size() && (rFindQ == refList.size() || qFindR - pQry <= rFindQ - pRef)) {
        finalLinearVBlocks.push_back(qVid);
        processed.insert(qVid);
        pQry++;
      } else if (rFindQ < refList.size()) {
        finalLinearVBlocks.push_back(rVid);
        processed.insert(rVid);
        pRef++;
      } else {
        finalLinearVBlocks.push_back(rVid);
        processed.insert(rVid);
        pRef++;
      }
    }
  }

  this->linear_block_cache.clear();

  // 輔助 Lambda：取得特定 VBlock 在給定 Sequence 上的起點座標
  auto getVBlockSegStart = [](const std::shared_ptr<Block>& b, int copy, const std::string& targetSeqName) -> int {
    if (!b) return -1;
    auto& seqs = b->getSequences();
    auto it = seqs.find(targetSeqName);
    if (it != seqs.end()) {
      for (auto& segPair : it->second.getSegments()) {
        if (segPair.second.getCopyCount() == copy) {
          return std::min(segPair.second.getStart(), segPair.second.getEnd());
        }
      }
      if (!it->second.getSegments().empty()) {
        auto& seg = it->second.getSegments().begin()->second;
        return std::min(seg.getStart(), seg.getEnd());
      }
    }
    return -1;
  };

  // 4. 🌟 將所有 Distant VBlocks 根據實體 Sequence 上最近的錨點 Block 插進線性順序中
  for (auto &blkWeak : this->getAllBlocks()) {
    auto blk = blkWeak.lock();
    if (!blk) continue;

    int maxCopy = blk->getMaxCopy();
    for (int c = 0; c <= maxCopy; ++c) {
      if (blk->isDistant(c)) {
        VBlockID distVId = {blk->getId(), c};

        // 尋找此 Distant Block 出現的實體序列名稱與實體座標
        std::string sampleSeq = "";
        int distStart = -1;

        for (auto &seqPair : blk->getSequences()) {
          for (auto &segPair : seqPair.second.getSegments()) {
            if (segPair.second.getCopyCount() == c) {
              sampleSeq = seqPair.first;
              distStart = std::min(segPair.second.getStart(), segPair.second.getEnd());
              break;
            }
          }
          if (!sampleSeq.empty()) break;
        }

        if (sampleSeq.empty() && !blk->getSequences().empty()) {
          sampleSeq = blk->getSequences().begin()->first;
          auto &seg = blk->getSequences().begin()->second.getSegments().begin()->second;
          distStart = std::min(seg.getStart(), seg.getEnd());
        }

        // 在已排好的 finalLinearVBlocks 中尋找同序列上最接近的非-Distant 錨點
        int bestInsertIdx = -1;
        int maxPrevStart = -1;
        int minNextStart = std::numeric_limits<int>::max();
        int minNextIdx = -1;

        for (int i = 0; i < (int)finalLinearVBlocks.size(); ++i) {
          auto anchorBlk = this->getBlock(finalLinearVBlocks[i].first);
          int anchorCopy = finalLinearVBlocks[i].second;
          int anchorStart = getVBlockSegStart(anchorBlk, anchorCopy, sampleSeq);

          if (anchorStart != -1) {
            if (anchorStart <= distStart) {
              if (anchorStart > maxPrevStart) {
                maxPrevStart = anchorStart;
                bestInsertIdx = i;
              }
            } else {
              if (anchorStart < minNextStart) {
                minNextStart = anchorStart;
                minNextIdx = i;
              }
            }
          }
        }

        // 根據找到的鄰居錨點插入
        if (bestInsertIdx != -1) {
          finalLinearVBlocks.insert(finalLinearVBlocks.begin() + bestInsertIdx + 1, distVId);
        } else if (minNextIdx != -1) {
          finalLinearVBlocks.insert(finalLinearVBlocks.begin() + minNextIdx, distVId);
        } else {
          finalLinearVBlocks.push_back(distVId);
        }
      }
    }
  }

  // 5. 寫回實體 Block 與 Segment 的 prev/next 指標以及 linear_block_cache
  this->linear_block_cache.reserve(finalLinearVBlocks.size());

  for (size_t i = 0; i < finalLinearVBlocks.size(); ++i) {
    VBlockID currVId = finalLinearVBlocks[i];
    auto currBlk = this->getBlock(currVId.first);
    int currCopy = currVId.second;

    this->linear_block_cache.push_back(currVId);

    if (currBlk) {
      auto prevBlkPtr = (i > 0) ? this->getBlock(finalLinearVBlocks[i - 1].first) : nullptr;
      auto nextBlkPtr = (i + 1 < finalLinearVBlocks.size()) ? this->getBlock(finalLinearVBlocks[i + 1].first) : nullptr;

      // 更新 Block 層級指標
      currBlk->setPrevBlock(currCopy, prevBlkPtr);
      currBlk->setNextBlock(currCopy, nextBlkPtr);

      // 🌟 同步更新 Segment 層級指標（避免二層以上 Progressive Alignment 時 Segment 連結失效）
      for (auto &seqPair : currBlk->getSequences()) {
        for (auto &segPairInner : seqPair.second.getSegments()) {
          Segment &seg = segPairInner.second;
          if (seg.getCopyCount() == currCopy) {
            if (!seg.isReverse()) {
              seg.setPrevBlock(prevBlkPtr);
              seg.setNextBlock(nextBlkPtr);
            } else {
              seg.setPrevBlock(nextBlkPtr);
              seg.setNextBlock(prevBlkPtr);
            }
          }
        }
      }
    }
  }

  // if (tree_ptr) {
  //   this->setDistantBlocks(*tree_ptr);
  // }
  this->buildCaches();

  std::cout << "  📊 [LinearGraph] Node '" << this->getId()
            << "' -> Linear Graph VBlocks: " << linear_block_cache.size()
            << " | Ancestral VBlocks: " << ancestral_block_cache.size()
            << " | Ancestral Seq Length: " << ancestral_seq_cache.length() << " bp\n";
}

const VBlockIDs &BlockSet::getLinearizeBlocks() {
  if (is_cached)
    return linear_block_cache;
  if (this->getSequenceCount() > 1) {
    std::cerr << "[ERROR] getLinearizeBlocks cache miss for non-leaf node (sequence count = "
              << this->getSequenceCount() << " > 1)! Non-leaf nodes must be rebuilt via rebuildLinearGraph.\n";
    exit(1);
  }
  this->rebuildAllPointers();
  return linear_block_cache;
}

const VBlockIDs &BlockSet::getAncestralBlocks() {
  if (is_cached)
    return ancestral_block_cache;
  if (this->getSequenceCount() > 1) {
    std::cerr << "[ERROR] getAncestralBlocks cache miss for non-leaf node (sequence count = "
              << this->getSequenceCount() << " > 1)! Non-leaf nodes must be rebuilt via rebuildLinearGraph.\n";
    exit(1);
  }
  this->rebuildAllPointers();
  return ancestral_block_cache;
}

const std::string &BlockSet::getAncestralSequence() {
  if (is_cached)
    return ancestral_seq_cache;
  if (this->getSequenceCount() > 1) {
    std::cerr << "[ERROR] getAncestralSequence cache miss for non-leaf node (sequence count = "
              << this->getSequenceCount() << " > 1)! Non-leaf nodes must be rebuilt via rebuildLinearGraph.\n";
    exit(1);
  }
  this->rebuildAllPointers();
  return ancestral_seq_cache;
}

const std::map<int, VBlockID> &BlockSet::getAncestralBlocksOffsets() {
  if (is_cached)
    return ancestral_offset_cache;
  if (this->getSequenceCount() > 1) {
    std::cerr << "[ERROR] getAncestralBlocksOffsets cache miss for non-leaf node (sequence count = "
              << this->getSequenceCount() << " > 1)! Non-leaf nodes must be rebuilt via rebuildLinearGraph.\n";
    exit(1);
  }
  this->rebuildAllPointers();
  return ancestral_offset_cache;
}

std::string BlockSet::reconstructSequence(const std::string &seqName) {
  struct SegNode {
    int start;
    int end;
    bool isRev;
    Variants vars;
    BlockPtr blk;
  };
  std::vector<SegNode> ordered_segments;

  // 1. 收集該序列散落在全圖的所有 Segments
  for (auto &blkPair : blocks) {
    std::shared_ptr<Block> blk = blkPair.second;
    auto &seqs = blk->getSequences();

    auto it = seqs.find(seqName);
    if (it != seqs.end()) {
      for (auto &segPair : it->second.getSegments()) {
        ordered_segments.push_back(
            {segPair.second.getStart(), segPair.second.getEnd(),
             segPair.second.isReverse(), segPair.second.getVariants(), blk});
      }
    }
  }

  if (ordered_segments.empty()) {
    std::cerr << "[Warning] Sequence '" << seqName << "' not found in BlockSet "
              << ID << ".\n";
    return "";
  }

  // 2. 依照基因體真實座標排序 (保證從 5' 端一路拼到 3' 端)
  std::sort(ordered_segments.begin(), ordered_segments.end(),
            [](const SegNode &a, const SegNode &b) {
              return std::min(a.start, a.end) < std::min(b.start, b.end);
            });

  // 3. 依序重建序列
  std::string reconstructedSeq = "";
  for (auto &node : ordered_segments) {
    std::string block_seq = node.blk->getConsensus().getConsensusString();

    // 步驟 A: 應用 SNV
    for (auto &v : node.vars) {
      if (v.getType() == VariantType::SNV) {
        // 安全檢查，避免越界
        if (v.getStart() < block_seq.length()) {
          block_seq[v.getStart()] = v.getAlt();
        }
      }
    }

    // 步驟 B: 應用 GAP
    std::string seg_seq = "";
    int cur = 0;

    Variants sorted_vars = node.vars;
    std::sort(
        sorted_vars.begin(), sorted_vars.end(),
        [](Variant &a, Variant &b) { return a.getStart() < b.getStart(); });

    for (auto &v : sorted_vars) {
      if (v.getType() == VariantType::GAP) {
        if (v.getStart() > cur) {
          seg_seq += block_seq.substr(cur, v.getStart() - cur);
        }
        cur = v.getEnd(); // 直接跳過 GAP
      }
    }
    if (cur < block_seq.length()) {
      seg_seq += block_seq.substr(cur);
    }

    // 步驟 C: 如果是反股 (-)，進行 Reverse Complement
    if (node.isRev) {
      seg_seq = getReverseComplement(seg_seq);
    }

    reconstructedSeq += seg_seq;
  }

  return reconstructedSeq;
}

BlockPtr BlockSet::concatBlocks(const std::vector<BlockID> &target_blocks, bool dryRun) {
  if (target_blocks.empty())
    return nullptr;

  // ==========================================
  // 步驟 1: 計算每個 Block 的長度與偏移量 (Offset)
  // ==========================================
  std::string super_consensus = "";
  std::unordered_map<BlockID, int> block_super_offsets;
  std::unordered_map<BlockID, int> block_lengths;

  int current_offset = 0;
  for (BlockID blkId : target_blocks) {
    auto blk = this->getBlock(blkId);
    if (!blk)
      continue;

    block_super_offsets[blkId] = current_offset;
    int blk_len = blk->getConsensus().length();
    block_lengths[blkId] = blk_len;

    super_consensus += blk->getConsensus().getConsensusString();
    current_offset += blk_len;
  }
  int total_super_len = super_consensus.length();

  // ==========================================
  // 步驟 2: 提取所有 Segments 並加上 Offset
  // ==========================================
  struct SegEntry {
    int offset;
    int len;
    Segment seg;
  };
  std::unordered_map<std::string, std::vector<SegEntry>> seq_to_segments;

  for (BlockID blkId : target_blocks) {
    auto blk = this->getBlock(blkId);
    if (!blk)
      continue;

    int block_offset = block_super_offsets[blkId];
    int block_len = block_lengths[blkId];

    // 不再區分 Copy，直接把所有 Sequences 與 Segments 掃出來
    for (auto &seqPair : blk->getSequences()) {
      const std::string &seqID = seqPair.first;
      for (auto &segPairInner : seqPair.second.getSegments()) {
        seq_to_segments[seqID].push_back(
            {block_offset, block_len, segPairInner.second});
      }
    }
  }

  std::vector<std::string> unique_seqs;
  std::vector<std::vector<SegEntry> *> seq_entries_ptrs;
  unique_seqs.reserve(seq_to_segments.size());
  seq_entries_ptrs.reserve(seq_to_segments.size());

  for (auto &kv : seq_to_segments) {
    unique_seqs.push_back(kv.first);
    seq_entries_ptrs.push_back(&kv.second);
  }

  struct Track {
    std::string seqID;
    Segment seg;
  };
  std::vector<std::vector<Track>> final_tracks_results(unique_seqs.size());

  // ==========================================
  // 步驟 3: TBB 平行處理 (Strand 感知與合併)
  // ==========================================
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, unique_seqs.size()),
      [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          const std::string &seqID = unique_seqs[i];
          std::vector<SegEntry> &entries = *(seq_entries_ptrs[i]);
          std::vector<Track> &final_tracks = final_tracks_results[i];

          if (entries.empty())
            continue;

          std::sort(entries.begin(), entries.end(),
                    [](const SegEntry &a, const SegEntry &b) {
                      if (a.offset != b.offset) return a.offset < b.offset;
                      return a.seg.getStart() < b.seg.getStart();
                    });

          Track current_track;
          current_track.seqID = seqID;
          bool has_active_track = false;

          int current_track_super_end = 0;
          int current_track_super_start = 0;

          // --- 關閉 Track ---
          auto close_active_track = [&]() {
            auto &vars = current_track.seg.getVariants();
            if (!current_track.seg.isReverse()) {
              if (current_track_super_end < total_super_len) {
                if (!vars.empty() &&
                    vars.back().getType() == VariantType::GAP &&
                    vars.back().getEnd() == current_track_super_end) {
                  int old_start = vars.back().getStart();
                  vars.pop_back();
                  vars.push_back(
                      Variant::createGap(old_start, total_super_len));
                } else {
                  vars.push_back(Variant::createGap(current_track_super_end,
                                                    total_super_len));
                }
              }
            } else {
              if (current_track_super_start > 0) {
                if (!vars.empty() &&
                    vars.front().getType() == VariantType::GAP &&
                    vars.front().getStart() == current_track_super_start) {
                  int old_end = vars.front().getEnd();
                  vars.erase(vars.begin());
                  Variants new_vars;
                  new_vars.push_back(Variant::createGap(0, old_end));
                  new_vars.insert(new_vars.end(),
                                  std::make_move_iterator(vars.begin()),
                                  std::make_move_iterator(vars.end()));
                  vars = std::move(new_vars);
                } else {
                  Variants new_vars;
                  new_vars.push_back(
                      Variant::createGap(0, current_track_super_start));
                  new_vars.insert(new_vars.end(),
                                  std::make_move_iterator(vars.begin()),
                                  std::make_move_iterator(vars.end()));
                  vars = std::move(new_vars);
                }
              }
            }
            final_tracks.push_back(std::move(current_track));
            has_active_track = false;
          };

          // --- 開啟新 Track ---
          auto open_new_track = [&](Segment seg, int offset, int len) {
            current_track.seg = std::move(seg);
            if (!current_track.seg.isReverse()) {
              if (offset > 0) {
                auto &vars = current_track.seg.getVariants();
                Variants new_vars;
                new_vars.reserve(vars.size() + 1);
                new_vars.push_back(Variant::createGap(0, offset));
                new_vars.insert(new_vars.end(),
                                std::make_move_iterator(vars.begin()),
                                std::make_move_iterator(vars.end()));
                vars = std::move(new_vars);
              }
              current_track_super_end = offset + len;
            } else {
              if (offset + len < total_super_len) {
                current_track.seg.getVariants().push_back(
                    Variant::createGap(offset + len, total_super_len));
              }
              current_track_super_start = offset;
            }
            has_active_track = true;
          };

          // --- 掃描與合併 ---
          for (size_t j = 0; j < entries.size(); ++j) {
            SegEntry &entry = entries[j];
            int block_offset = entry.offset;
            int block_len = entry.len;

            Segment processed_seg = entry.seg;
            for (auto &var : processed_seg.getVariants())
              var.shift(block_offset);

            if (!has_active_track) {
              open_new_track(std::move(processed_seg), block_offset, block_len);
              continue;
            }

            bool is_same_strand =
                (current_track.seg.isReverse() == processed_seg.isReverse());
            bool can_merge = (is_same_strand && current_track.seg.getEnd() ==
                                                    processed_seg.getStart());

            if (can_merge) {
              if (!current_track.seg.isReverse()) {
                if (block_offset >= current_track_super_end) {
                  current_track.seg.setEnd(processed_seg.getEnd());
                  auto &vars = current_track.seg.getVariants();
                  if (block_offset > current_track_super_end) {
                    if (!vars.empty() &&
                        vars.back().getType() == VariantType::GAP &&
                        vars.back().getEnd() == current_track_super_end) {
                      int old_start = vars.back().getStart();
                      vars.pop_back();
                      vars.push_back(
                          Variant::createGap(old_start, block_offset));
                    } else {
                      vars.push_back(Variant::createGap(current_track_super_end,
                                                        block_offset));
                    }
                  }
                  vars.insert(vars.end(),
                              std::make_move_iterator(
                                  processed_seg.getVariants().begin()),
                              std::make_move_iterator(
                                  processed_seg.getVariants().end()));
                  current_track_super_end = block_offset + block_len;
                } else
                  can_merge = false;
              } else {
                if (block_offset + block_len <= current_track_super_start) {
                  current_track.seg.setEnd(processed_seg.getEnd());
                  auto &vars = current_track.seg.getVariants();
                  Variants new_vars;
                  new_vars.reserve(processed_seg.getVariants().size() + 1 +
                                   vars.size());
                  new_vars.insert(new_vars.end(),
                                  std::make_move_iterator(
                                      processed_seg.getVariants().begin()),
                                  std::make_move_iterator(
                                      processed_seg.getVariants().end()));

                  if (block_offset + block_len < current_track_super_start) {
                    if (!vars.empty() &&
                        vars.front().getType() == VariantType::GAP &&
                        vars.front().getStart() == current_track_super_start) {
                      int old_end = vars.front().getEnd();
                      vars.erase(vars.begin());
                      new_vars.push_back(Variant::createGap(
                          block_offset + block_len, old_end));
                    } else {
                      new_vars.push_back(Variant::createGap(
                          block_offset + block_len, current_track_super_start));
                    }
                  }
                  new_vars.insert(new_vars.end(),
                                  std::make_move_iterator(vars.begin()),
                                  std::make_move_iterator(vars.end()));
                  vars = std::move(new_vars);
                  current_track_super_start = block_offset;
                } else
                  can_merge = false;
              }
            }

            if (!can_merge) {
              close_active_track();
              open_new_track(std::move(processed_seg), block_offset, block_len);
            }
          }
          if (has_active_track)
            close_active_track();
        }
      });

  // ==========================================
  // 步驟 4: 組裝 Final Super Block
  // ==========================================
  BlockPtr super_block;

  if (dryRun) {
    super_block = std::make_shared<Block>(999999, "");
  } else {
    super_block = this->createBlock(super_consensus);
  }
  for (size_t i = 0; i < unique_seqs.size(); ++i) {
    if (final_tracks_results[i].empty())
      continue;
    Sequence seq_info(unique_seqs[i]);
    for (auto &track : final_tracks_results[i])
      seq_info.addSegment(track.seg);
    super_block->addSequence(seq_info);
  }

  return super_block;
}

// ==========================================
// 將外部的 Block 加入此 BlockSet 並賦予新 ID
// ==========================================
BlockPtr BlockSet::addBlock(BlockPtr oldBlock) {
  if (!oldBlock)
    return nullptr;

  BlockID newId = next_block_id_++;

  // 3. 利用舊 Block 的 Consensus 建立全新的 Block
  auto newBlock = std::make_shared<Block>(newId, oldBlock->getConsensus());

  // 4. 深拷貝：將所有的 SequenceInfo 複製過去
  for (const auto &seqPair : oldBlock->getSequences()) {
    newBlock->addSequence(seqPair.second);
  }

  // 5. 註冊進這個 BlockSet 的 Dictionary 中
  blocks[newId] = newBlock;

  invalidateRepCache();

  return newBlock; // 回傳新建立的 Block 智慧指標
}

BlockPtr BlockSet::addBlockReferencing(BlockPtr oldBlock, BlockSet* sourceSet) {
  if (!oldBlock)
    return nullptr;

  BlockID newId = next_block_id_++;

  // 1. 尋找 oldBlock 在 sourceSet 的 ancestral block offsets 裡的 start offset
  int start_offset = -1;
  int consensus_len = oldBlock->getConsensus().length();

  if (sourceSet) {
    const auto& offsets = sourceSet->getAncestralBlocksOffsets();
    for (const auto& pair : offsets) {
      if (pair.second.first == oldBlock->getId()) {
        start_offset = pair.first;
        break;
      }
    }
  }

  // 2. 建立新 Consensus 物件，不拷貝 string，而是用 index 參照
  Consensus newConsensus;
  if (start_offset != -1 && sourceSet) {
    newConsensus = Consensus(sourceSet, start_offset, start_offset + consensus_len, {});
  } else {
    // 防呆：如果找不到 offset，則 fallback 回原本的 consensus
    newConsensus = oldBlock->getConsensus();
  }

  // 3. 利用新建立的 Consensus 建立全新的 Block
  auto newBlock = std::make_shared<Block>(newId, std::move(newConsensus));

  // 4. 深拷貝：將所有的 SequenceInfo 複製過去
  for (const auto &seqPair : oldBlock->getSequences()) {
    newBlock->addSequence(seqPair.second);
  }

  // 5. 註冊進這個 BlockSet 的 Dictionary 中
  blocks[newId] = newBlock;

  invalidateRepCache();

  return newBlock; // 回傳新建立的 Block 智慧指標
}

std::pair<BlockID, BlockID> BlockSet::splitSingleBlock(int parentID, int localCut) {

  auto parent = this->getBlock(parentID);

  if (!parent)
    return {(uint64_t)-1, (uint64_t)-1};

  bool debug = false;

  // 1. 建立新 Block (左右半部)
  auto left = this->createBlock(parent->getConsensus().substr(0, localCut));
  auto right = this->createBlock(parent->getConsensus().substr(localCut));

  // ==========================================
  // 2. TBB 平行化：切割 Sequence 與 Segment
  // ==========================================
  std::vector<std::string> seqIDs;
  seqIDs.reserve(parent->getSequences().size());
  for (auto &kv : parent->getSequences()) {
    seqIDs.push_back(kv.first);
  }

  // 用來儲存平行切割結果的暫存結構，避免 Thread Contention
  struct SplitResult {
    Sequence leftSeq;
    Sequence rightSeq;
    bool hasLeft = false;
    bool hasRight = false;
  };
  std::vector<SplitResult> splitResults(seqIDs.size());

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, seqIDs.size()),
      [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          const std::string &seqID = seqIDs[i];
          auto &parentSeqInfo = parent->getSequences().at(seqID);

          Sequence leftSeqInfo(seqID);
          Sequence rightSeqInfo(seqID);

          for (auto &segPair : parentSeqInfo.getSegments()) {
            Segment oldSeg = segPair.second; // 拷貝出來處理
            auto splitSegs = oldSeg.split(localCut);
            Segment &leftSeg = splitSegs.first;
            Segment &rightSeg = splitSegs.second;

            // 1. 先判斷這個 Segment 切出來後，是否真實擁有物理序列 (非純 Gap)
            bool validLeft = (leftSeg.getStart() != leftSeg.getEnd());
            bool validRight = (rightSeg.getStart() != rightSeg.getEnd());

            // 2. 內部接線：嚴格遵守正反股走向，且「只對真正存在的 Segment
            // 接線」
            if (validLeft && validRight) {
              // 兩邊都有肉：互相連接，並對外連接
              if (!oldSeg.isReverse()) {
                leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                leftSeg.setNextBlock(right);
                rightSeg.setPrevBlock(left);
                rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
              } else {
                rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                rightSeg.setNextBlock(left);
                leftSeg.setPrevBlock(right);
                leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
              }
            } else if (validLeft && !validRight) {
              // 只有左邊有肉：左邊直接繼承原 Segment 的所有對外連接
              leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
              leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
            } else if (!validLeft && validRight) {
              // 只有右邊有肉：右邊直接繼承原 Segment 的所有對外連接
              rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
              rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
            }
            // 如果兩邊都沒肉 (!validLeft &&
            // !validRight)，那就什麼都不用接，直接丟棄

            // 3. 將真實存在的 Segment 放入 Map 裡
            if (validLeft) {
              leftSeqInfo.getSegments()[leftSeg.getStart()] = leftSeg;
            }
            if (validRight) {
              rightSeqInfo.getSegments()[rightSeg.getStart()] = rightSeg;
            }
          }

          // 將結果存入專屬的 index，確保 Thread Safe
          splitResults[i].leftSeq = std::move(leftSeqInfo);
          splitResults[i].rightSeq = std::move(rightSeqInfo);
          splitResults[i].hasLeft =
              !splitResults[i].leftSeq.getSegments().empty();
          splitResults[i].hasRight =
              !splitResults[i].rightSeq.getSegments().empty();
        }
      });

  // 主執行緒快速合併結果 (將 Map 搬進 left/right block)
  for (size_t i = 0; i < seqIDs.size(); ++i) {
    if (splitResults[i].hasLeft)
      left->addSequence(std::move(splitResults[i].leftSeq));
    if (splitResults[i].hasRight)
      right->addSequence(std::move(splitResults[i].rightSeq));
  }

  // ==========================================
  // 3. 【修復核心】：消滅全圖掃描，改為「鄰居局部掃描」+ TBB 平行接線
  // ==========================================
  std::unordered_set<std::shared_ptr<Block>> neighbors;

  // 必須加入 left 和 right 來解開 Self-loop
  neighbors.insert(left);
  neighbors.insert(right);

  // 收集真正有牽連的鄰居 (只看 parent 原本的連線)
  for (auto &seqPair : parent->getSequences()) {
    for (auto &segPair : seqPair.second.getSegments()) {
      if (auto p = segPair.second.getPrevBlock().lock())
        neighbors.insert(p);
      if (auto n = segPair.second.getNextBlock().lock())
        neighbors.insert(n);
    }
  }
  neighbors.erase(parent); // parent 即將被刪除，不用幫它接線

  // 將鄰居轉為 Vector 以供 TBB 平行處理
  std::vector<std::shared_ptr<Block>> neighbor_vec(neighbors.begin(),
                                                   neighbors.end());
  std::atomic<int> rewiredCount{0};

  // TBB 平行接線：因為每條 Thread 處理不同的 Neighbor Block，
  // 其內部的 Segment 也是獨立的，因此絕對 Thread Safe！
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, neighbor_vec.size()),
      [&](const tbb::blocked_range<size_t> &r) {
        int local_rewired = 0;
        for (size_t i = r.begin(); i != r.end(); ++i) {
          auto currentBlock = neighbor_vec[i];

          for (auto &seqPair : currentBlock->getSequences()) {
            std::string seqID = seqPair.first;
            for (auto &segPair : seqPair.second.getSegments()) {
              Segment &seg = segPair.second;

              // 檢查 Prev：如果這段序列是從 parent 來的
              if (seg.getPrevBlock().lock() == parent) {
                bool found = false;
                if (left->getSequences().count(seqID)) {
                  for (auto &lSeg :
                       left->getSequences().at(seqID).getSegments()) {
                    if (lSeg.second.getEnd() == seg.getStart()) {
                      seg.setPrevBlock(left);
                      found = true;
                      local_rewired++;
                      break;
                    }
                  }
                }
                if (!found && right->getSequences().count(seqID)) {
                  for (auto &rSeg :
                       right->getSequences().at(seqID).getSegments()) {
                    if (rSeg.second.getEnd() == seg.getStart()) {
                      seg.setPrevBlock(right);
                      local_rewired++;
                      break;
                    }
                  }
                }
              }

              // 檢查 Next：如果這段序列下一步要走到 parent
              if (seg.getNextBlock().lock() == parent) {
                bool found = false;
                if (left->getSequences().count(seqID)) {
                  for (auto &lSeg :
                       left->getSequences().at(seqID).getSegments()) {
                    if (lSeg.second.getStart() == seg.getEnd()) {
                      seg.setNextBlock(left);
                      found = true;
                      local_rewired++;
                      break;
                    }
                  }
                }
                if (!found && right->getSequences().count(seqID)) {
                  for (auto &rSeg :
                       right->getSequences().at(seqID).getSegments()) {
                    if (rSeg.second.getStart() == seg.getEnd()) {
                      seg.setNextBlock(right);
                      local_rewired++;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
        rewiredCount += local_rewired; // Atomic 累加
      });

  if (debug)
    std::cout << "[DEBUG-SPLIT] Block " << parentID
              << " (Len: " << parent->getConsensus().size() << ") cut at "
              << localCut << " -> L: " << left->getId()
              << " (Len: " << left->getConsensus().size()
              << "), R: " << right->getId()
              << " (Len: " << right->getConsensus().size()
              << " | Rewired pointers: " << rewiredCount.load() << "\n";

  // 4. 安全刪除舊 Block
  this->deleteBlock(parent->getId());

  return {left->getId(), right->getId()};
}

void BlockSet::print(std::ostream &os) const {
  os << "\n===================================================================="
        "======\n";
  os << " 🌐 BLOCK SET ID: " << ID
     << " | Total Blocks In Map: " << blocks.size() << "\n";
  os << "======================================================================"
        "====\n";

  // 因為 getLinearizeBlocks() 在宣告中是非 const，這裡我們透過 const_cast
  // 來安全調用
  auto &mutableSet = const_cast<BlockSet &>(*this);
  VBlockIDs linearBlocks = mutableSet.getLinearizeBlocks();

  if (linearBlocks.empty()) {
    os << "  ⚠️  [Warning] Graph is empty or contains no linearized backbone "
          "blocks.\n";
    os << "===================================================================="
          "======\n\n";
    return;
  }

  // 依照 Linearized 順序逐一印出 Block
  for (size_t i = 0; i < linearBlocks.size(); ++i) {
    VBlockID vbid = linearBlocks[i];
    auto blk = mutableSet.getBlock(vbid.first);

    if (blk) {
      blk->print(os, vbid.second);
      // 如果後面還有 Block，印出一個漂亮的拓撲流向箭頭
      if (i + 1 < linearBlocks.size()) {
        os << "                                   │\n";
        os << "                                   ▼\n";
      }
    } else {
      os << "  ❌ [ERROR] Block ID " << vbid.first
         << " listed in linear backbone but missing from map!\n";
    }
  }
  os << "======================================================================"
        "====\n\n";
}

void BlockSet::setDistantBlocks(Tree &tree, int lookdownDepth) {
  bool debug = true;
  BlockSetID targetNodeId = this->getId();

  auto it = tree.allNodes.find(targetNodeId);
  if (it == tree.allNodes.end())
    return;
  Node *currentNode = it->second;

  if (currentNode->is_leaf() || currentNode->children.size() < 2)
    return;

  std::vector<std::unordered_set<std::string>> sequenceSets;
  tree.getSubLineages(currentNode, lookdownDepth, 0, sequenceSets);

  uint64_t minSequenceSets = 1ULL << lookdownDepth;
  if (sequenceSets.size() < minSequenceSets) {
    if (debug) {
      std::cout << "[DISTANT-SKIP] Node " << targetNodeId
                << ": sequenceSets.size()=" << sequenceSets.size()
                << " < minSequenceSets=" << minSequenceSets
                << ", skipping distant evaluation.\n";
    }
    return;
  }

  // 🔍 Debug: print all sequence sets
  if (debug) {
    std::cout << "[DISTANT-DEBUG] Node " << targetNodeId
              << ": " << sequenceSets.size() << " sequence sets:\n";
    for (size_t i = 0; i < sequenceSets.size(); ++i) {
      std::cout << "  Set " << i << ": {";
      bool first = true;
      for (const auto &s : sequenceSets[i]) {
        if (!first) std::cout << ", ";
        std::cout << s;
        first = false;
      }
      std::cout << "}\n";
    }
  }

  auto all_blocks = this->getAllBlocks();
  std::atomic<int> distantCount{0};
  std::atomic<int> coreCount{0};

  // 🔍 Use sequential loop for debug to avoid interleaved output
  for (const auto &weak_blk : all_blocks) {
        auto blk = weak_blk.lock();
        if (!blk)
          continue;

        // 🌟 步驟 1：將 Block 內的 Sequence 依照 Copy 進行分組
        // 資料結構：Map<Copy號碼, Set<Sequence名稱>>
        std::unordered_map<int, std::unordered_set<std::string>> copyToSeqNames;

        for (auto &seqPair : blk->getSequences()) {
          const std::string &seqName = seqPair.first;
          // 遍歷該 Sequence 底下的 Segments 來取得它們的 Copy
          for (auto &segPairInner : seqPair.second.getSegments()) {
            int copyNum = segPairInner.second
                              .getCopyCount();
            copyToSeqNames[copyNum].insert(seqName);
          }
        }

        // 🌟 步驟 2：針對每一個 Copy 獨立進行投票與判定
        for (const auto &copyPair : copyToSeqNames) {
          int currentCopy = copyPair.first;
          const auto &copySeqNames = copyPair.second;

          int supportedSets = 0;
          for (const auto &leafSet : sequenceSets) {
            for (const auto &seqName : copySeqNames) {
              if (leafSet.find(seqName) != leafSet.end()) {
                supportedSets++;
                break;
              }
            }
          }

          bool isDistant = (supportedSets < 2);

          // 🔍 Debug: print per-block per-copy classification
          if (debug && isDistant) {
            std::cout << "  [DISTANT-DETAIL] Block " << blk->getId()
                      << " Copy " << currentCopy
                      << " -> supportedSets=" << supportedSets
                      << " -> DISTANT | seqs={";
            bool first = true;
            for (const auto &s : copySeqNames) {
              if (!first) std::cout << ", ";
              std::cout << s;
              first = false;
            }
            std::cout << "}\n";
          }

          // 🌟 步驟 3：使用包含 Copy 參數的新 Setter
          if (!isDistant) {
            blk->setDistant(false, currentCopy);
            coreCount++;
          } else {
            blk->setDistant(true, currentCopy);
            distantCount++;
          }
        }
  }

  if (debug) {
    std::cout << "[INFO] Node " << targetNodeId << " (Lookdown "
              << lookdownDepth << " levels -> " << sequenceSets.size()
              << " sets)"
              << " | Core VBlocks=" << coreCount
              << ", Distant VBlocks=" << distantCount << "\n";
  }
}


std::shared_ptr<Block> BlockSet::extractBlock(int extract_start, int extract_end) {
  if (extract_start >= extract_end)
    return nullptr;

  global_timer.start("0 Loop");

  const auto &offset_map = this->getAncestralBlocksOffsets();
  if (offset_map.empty())
    return nullptr;

  auto it_start = offset_map.upper_bound(extract_start);
  if (it_start != offset_map.begin()) {
    --it_start;
  }

  auto it_end = offset_map.upper_bound(extract_end);

  std::vector<VBlockID> target_vblocks_to_concat;
  std::string combined_consensus = "";
  int current_concat_offset = 0;

  struct SegEntry {
    int offset;
    int len;
    Segment seg;
  };
  std::unordered_map<std::string, std::vector<SegEntry>> seq_to_segments;
  global_timer.stop("0 Loop");
  // ==========================================
  // 🌟 只針對「命中區間」的這幾個積木進行走訪
  // ==========================================
  global_timer.start("1st Loop");
  for (auto it = it_start; it != it_end; ++it) {
    int vblk_global_start = it->first;
    VBlockID current_vid = it->second;
    int vblk_len = this->getBlock(current_vid.first)->getConsensus().length();
    int vblk_global_end = vblk_global_start + vblk_len;

    auto blk = this->getBlock(current_vid.first);
    if (!blk)
      continue;

    int target_copy = current_vid.second;
    std::shared_ptr<Block> block_to_process = blk;

    bool need_split_left =
        (extract_start > vblk_global_start && extract_start < vblk_global_end);
    bool need_split_right =
        (extract_end > vblk_global_start && extract_end < vblk_global_end);

    global_timer.start("1st Loop - Split");
    if (need_split_left && need_split_right) {
      int cut1 = extract_start - vblk_global_start;
      auto temp_right = blk->split(cut1, target_copy).second;
      int cut2 = extract_end - extract_start;
      block_to_process = temp_right->split(cut2, target_copy).first;
    } else if (need_split_left) {
      int cut = extract_start - vblk_global_start;
      block_to_process = blk->split(cut, target_copy).second;
    } else if (need_split_right) {
      int cut = extract_end - vblk_global_start;
      block_to_process = blk->split(cut, target_copy).first;
    }
    global_timer.stop("1st Loop - Split");

    if (!block_to_process)
      continue;

    int processed_len = block_to_process->getConsensus().length();
    for (auto &[seqID, seqInfo] : block_to_process->getSequences()) {
      for (auto &[segStart, seg] : seqInfo.getSegments()) {
        if (seg.getCopyCount() == target_copy) {
          seq_to_segments[seqID].push_back(
              {current_concat_offset, processed_len, seg});
        }
      }
    }

    combined_consensus += block_to_process->getConsensus().getConsensusString();
    current_concat_offset += processed_len;
  }
  global_timer.stop("1st Loop");

  if (current_concat_offset == 0)
    return nullptr;

  // 建立全新的獨立巨型 Block (也可以考慮傳入自訂的 sandbox ID)
  int total_len = current_concat_offset;
  BlockID dummy_id = 0;
  auto newBlock = std::make_shared<Block>(dummy_id, combined_consensus);

  // ==========================================
  // 3. 針對每一條 Sequence 進行物理片段縫合與補 Gap
  // ==========================================
  global_timer.start("2nd Loop");
  for (auto &[seqID, entries] : seq_to_segments) {
    if (entries.empty())
      continue;

    // 嚴格按照生物真實座標 (Start) 排序
    std::sort(entries.begin(), entries.end(), [](SegEntry &a, SegEntry &b) {
      return a.seg.getStart() < b.seg.getStart();
    });

    struct Track {
      Segment seg;
    };
    std::vector<Track> final_tracks;
    Track current_track;
    bool has_active_track = false;

    int current_track_super_end = 0;
    int current_track_super_start = 0;

    auto close_active_track = [&]() {
      auto &vars = current_track.seg.getVariants();
      if (!current_track.seg.isReverse()) {
        if (current_track_super_end < total_len) {
          if (!vars.empty() && vars.back().getType() == VariantType::GAP &&
              vars.back().getEnd() == current_track_super_end) {
            int old_start = vars.back().getStart();
            vars.pop_back();
            vars.push_back(Variant::createGap(old_start, total_len));
          } else {
            vars.push_back(
                Variant::createGap(current_track_super_end, total_len));
          }
        }
      } else {
        if (current_track_super_start > 0) {
          if (!vars.empty() && vars.front().getType() == VariantType::GAP &&
              vars.front().getStart() == current_track_super_start) {
            int old_end = vars.front().getEnd();
            vars.erase(vars.begin());
            std::vector<Variant> new_vars;
            new_vars.push_back(Variant::createGap(0, old_end));
            new_vars.insert(new_vars.end(),
                            std::make_move_iterator(vars.begin()),
                            std::make_move_iterator(vars.end()));
            vars = std::move(new_vars);
          } else {
            std::vector<Variant> new_vars;
            new_vars.push_back(
                Variant::createGap(0, current_track_super_start));
            new_vars.insert(new_vars.end(),
                            std::make_move_iterator(vars.begin()),
                            std::make_move_iterator(vars.end()));
            vars = std::move(new_vars);
          }
        }
      }
      final_tracks.push_back(std::move(current_track));
      has_active_track = false;
    };

    auto open_new_track = [&](Segment seg, int offset, int len) {
      current_track.seg = std::move(seg);
      if (!current_track.seg.isReverse()) {
        if (offset > 0) {
          auto &vars = current_track.seg.getVariants();
          std::vector<Variant> new_vars;
          new_vars.reserve(vars.size() + 1);
          new_vars.push_back(Variant::createGap(0, offset));
          new_vars.insert(new_vars.end(), std::make_move_iterator(vars.begin()),
                          std::make_move_iterator(vars.end()));
          vars = std::move(new_vars);
        }
        current_track_super_end = offset + len;
      } else {
        if (offset + len < total_len) {
          current_track.seg.getVariants().push_back(
              Variant::createGap(offset + len, total_len));
        }
        current_track_super_start = offset;
      }
      has_active_track = true;
    };

    for (auto &entry : entries) {
      int block_offset = entry.offset;
      int block_len = entry.len;

      Segment processed_seg = entry.seg;
      for (auto &var : processed_seg.getVariants()) {
        var.shift(block_offset);
      }

      if (!has_active_track) {
        open_new_track(std::move(processed_seg), block_offset, block_len);
        continue;
      }

      bool is_same_strand =
          (current_track.seg.isReverse() == processed_seg.isReverse());
      bool can_merge = is_same_strand &&
                       (current_track.seg.getEnd() == processed_seg.getStart());

      if (can_merge) {
        if (!current_track.seg.isReverse()) {
          if (block_offset >= current_track_super_end) {
            current_track.seg.setEnd(processed_seg.getEnd());
            auto &vars = current_track.seg.getVariants();

            if (block_offset > current_track_super_end) {
              if (!vars.empty() && vars.back().getType() == VariantType::GAP &&
                  vars.back().getEnd() == current_track_super_end) {
                int old_start = vars.back().getStart();
                vars.pop_back();
                vars.push_back(Variant::createGap(old_start, block_offset));
              } else {
                vars.push_back(
                    Variant::createGap(current_track_super_end, block_offset));
              }
            }
            vars.insert(
                vars.end(),
                std::make_move_iterator(processed_seg.getVariants().begin()),
                std::make_move_iterator(processed_seg.getVariants().end()));
            current_track_super_end = block_offset + block_len;
          } else
            can_merge = false;
        } else {
          if (block_offset + block_len <= current_track_super_start) {
            current_track.seg.setEnd(processed_seg.getEnd());
            auto &vars = current_track.seg.getVariants();
            std::vector<Variant> new_vars;
            new_vars.reserve(processed_seg.getVariants().size() + 1 +
                             vars.size());

            new_vars.insert(
                new_vars.end(),
                std::make_move_iterator(processed_seg.getVariants().begin()),
                std::make_move_iterator(processed_seg.getVariants().end()));

            if (block_offset + block_len < current_track_super_start) {
              if (!vars.empty() && vars.front().getType() == VariantType::GAP &&
                  vars.front().getStart() == current_track_super_start) {
                int old_end = vars.front().getEnd();
                vars.erase(vars.begin());
                new_vars.push_back(
                    Variant::createGap(block_offset + block_len, old_end));
              } else {
                new_vars.push_back(Variant::createGap(
                    block_offset + block_len, current_track_super_start));
              }
            }
            new_vars.insert(new_vars.end(),
                            std::make_move_iterator(vars.begin()),
                            std::make_move_iterator(vars.end()));
            vars = std::move(new_vars);
            current_track_super_start = block_offset;
          } else
            can_merge = false;
        }
      }

      if (!can_merge) {
        close_active_track();
        open_new_track(std::move(processed_seg), block_offset, block_len);
      }
    }

    if (has_active_track)
      close_active_track();

    Sequence newSeqInfo(seqID);
    for (auto &track : final_tracks) {
      newSeqInfo.addSegment(track.seg);
    }
    newBlock->addSequence(newSeqInfo);
  }
  global_timer.stop("2nd Loop");
  return newBlock;
}

BlockSet *BlockSet::createIndexedCopy(BlockManager *manager, const BlockSetID &newID) {
  BlockSet *indexedSet = manager->createBlockSet(newID);
  indexedSet->setTree(this->getTree());
  for (const auto &seq : this->getSequences()) {
    indexedSet->addSequenceName(seq);
  }

  int current_offset = 0;
  for (const auto &blkPair : this->blocks) {
    BlockID origId = blkPair.first;
    std::shared_ptr<Block> origBlk = blkPair.second;
    if (!origBlk)
      continue;

    int consLen = origBlk->getConsensus().getConsensusString().length();
    int start_idx = current_offset;
    int end_idx = current_offset + consLen;
    current_offset = end_idx;

    // 建立指向 *this (原 BlockSet) 的 index/reference 模式 Consensus (零大字串複製)
    Consensus indexedCons(this, start_idx, end_idx);

    // 強制保留原來的 BlockID
    BlockPtr newBlk = indexedSet->createBlockWithId(origId, std::move(indexedCons));

    // 複製 Sequences 與 Segments
    for (const auto &seqPair : origBlk->getSequences()) {
      newBlk->addSequence(seqPair.second);
    }
  }

  indexedSet->rebuildAllPointers();
  return indexedSet;
}

