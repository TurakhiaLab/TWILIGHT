#include "alignment.hpp"
#include "block_manager.hpp"
#include "coordinate_manager.hpp"
#include "global_alignment.hpp"
#include "timer.hpp"
#include "type.hpp"
#include "cigar_util.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <set>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>
#include <tuple>
#include <vector>

BlockSet *BlockManager::merge(BlockSet *refSet, BlockSet *qrySet,
                              AlignmentCollection &alnCollection,
                              BlockSetID newID, Tree *tree, int L_min) {
  bool DEBUG_MODE = (refSet->getSequenceCount() > 1 || qrySet->getSequenceCount() > 1);
  // bool DEBUG_MODE = false;
  auto time0 = std::chrono::high_resolution_clock::now();
  if (DEBUG_MODE)
    std::cout << "\n========================================================\n"
              << "=== TWILIGHT-MGA DYNAMIC MERGE: " << refSet->getId() << " + "
              << qrySet->getId() << " ===\n"
              << "========================================================\n";

  global_timer.start("phase1_init_coord_mgr");
  if (DEBUG_MODE) std::cout << "[Phase 1] Initialized Coordinate Manager...\n";

  BlockSet *mergedSet = createBlockSet(newID);

  refSet->setTree(tree);
  qrySet->setTree(tree);
  mergedSet->setTree(tree);

  CoordinateManager coordMgr;
  coordMgr.init(refSet, qrySet, mergedSet);

  alnCollection.ref_coverageTracker.syncFromMap(coordMgr, true); // true 代表 Ref
  alnCollection.qry_coverageTracker.syncFromMap(coordMgr, false); // false 代表 Qry

  int mergeCounter = 0;
  global_timer.stop("phase1_init_coord_mgr");

  if (DEBUG_MODE) std::cout << "[Phase 2] Iteratively merged blocks based on Minimap2 alignments...\n";

  for (auto &seq : refSet->getSequences()) mergedSet->addSequenceName(seq);
  for (auto &seq : qrySet->getSequences()) mergedSet->addSequenceName(seq);

  // 🌟 呼叫 BlockSet::createIndexedCopy 產生保留 BlockID 且零大字串複製成本的 Index 版 BlockSet
  BlockSet *refIndexSet = refSet->createIndexedCopy(this, refSet->getId() + "_indexed");
  BlockSet *qryIndexSet = qrySet->createIndexedCopy(this, qrySet->getId() + "_indexed");

  global_timer.start("phase2_iterative_merge");
  while (true) {

    global_timer.start("find_best_alignment");
    Alignments bestAlns = alnCollection.getBestAlignments(refIndexSet, qryIndexSet, coordMgr);
    global_timer.stop("find_best_alignment");

    if (bestAlns.empty())
      break;

    for (Alignment &bestAln : bestAlns) {
      if (!bestAln.valid)
        break;

      bestAln.CIGAR = compressCigar(bestAln.CIGAR);

      // 1. Detect Pure Indel
      bool is_pure_insertion = (bestAln.CIGAR.size() == 1 && bestAln.CIGAR[0].second == 'I');
      bool is_pure_deletion = (bestAln.CIGAR.size() == 1 && bestAln.CIGAR[0].second == 'D');

      int r_start = bestAln.refIdx.first, r_end = bestAln.refIdx.second;
      int q_start = std::min(bestAln.qryIdx.first, bestAln.qryIdx.second);
      int q_end = std::max(bestAln.qryIdx.first, bestAln.qryIdx.second);

      bool r_merged = false, q_merged = false;
      global_timer.start("p2_2_overlap_lookup");
      auto r_overlaps = coordMgr.getOverlapBlocks(r_start, r_end, true, r_merged);
      auto q_overlaps = coordMgr.getOverlapBlocks(q_start, q_end, false, q_merged);
      global_timer.stop("p2_2_overlap_lookup");

      if (DEBUG_MODE) {
        std::cout << "\n========================================\n";
        std::cout << "🔍 [ALIGNMENT DETAILS]\n";
        std::cout << "  ├─ CIGAR : " << cigarToStr(bestAln.CIGAR) << (bestAln.inverse ? "(-)" : "(+)") <<"\n";
        std::cout << "  ├─ Ref   : [" << r_start << ", " << r_end << ")\n";
        std::cout << "  ├─ Qry   : [" << q_start << ", " << q_end << ")\n";

        auto printOverlap = [](const auto &list) {
          std::cout << "[ ";
          bool first = true;
          for (const auto &id : list) {
            if (!first)
              std::cout << ", ";
            std::cout << id;
            first = false;
          }
          std::cout << " ]\n";
        };

        if (r_merged) {
          std::cout << "  ├─ Ref Overlap : ";
          printOverlap(r_overlaps);
        }
        if (q_merged) {
          std::cout << "  └─ Qry Overlap : ";
          printOverlap(q_overlaps);
        }
        std::cout << "========================================\n";
      }

      if (r_merged && q_merged && (*r_overlaps.begin()) == (*q_overlaps.begin())) {
        if (DEBUG_MODE)
          std::cout
              << "  -> SKIP. Both paths resolve to the exact same merged block "
              << (*r_overlaps.begin()) << ".\n";
        continue;
      }

      int rLen = r_end - r_start;
      int qLen = q_end - q_start;

      if (rLen == 0 && qLen == 0) {
        if (DEBUG_MODE) {
          std::cout << "  [Null Alignment] Skip.\n";
        }
        continue;
      }

      // =======================================================
      // Pure Indel Fast Path
      // =======================================================
      if ((is_pure_deletion && rLen > 0) || (is_pure_insertion && qLen > 0)) {
        if (DEBUG_MODE) {
          std::cout << "  -> ⚡ [PURE INDEL] Bypassing logic. Block is already "
                       "safely in mergedSet.\n";
        }
        continue;
      }

      const int min_threshold = MERGE_MIN_THRESHOLD;
      if (!bestAln.primary && !(is_pure_deletion || is_pure_insertion)) {
        int score = scoreCIGAR(bestAln.CIGAR);
        if (score < min_threshold) {
          if (DEBUG_MODE) {
            std::cout
                << "  -> SKIP. Secondary alignment score below threshold ("
                << score << " < " << min_threshold << ").\n";
          }
          continue;
        }
      }

      BlockID rCore = 0, qCore = 0;

      // =======================================================
      // 🌟 獨立能量預檢系統 (只讀取，不修改任何圖譜狀態)
      // =======================================================
      auto isAlignmentRejected =
          [&](int globalStart, int globalEnd, bool isRef,
              const std::vector<BlockID> &overlaps) -> bool {
        if (overlaps.empty())
          return false;

        const int ALPHA = MERGE_ALPHA;
        int L = globalEnd - globalStart;
        int Nc = 0;

        // 輔助工具：安全獲取 Block 深度 (包含 Ref 與 Qry 的軌道數)
        auto getBlockDepth = [&](BlockID blkId) {
          CoordinateManager::GlobalCoords coords =
              coordMgr.getGlobalCoord(blkId);
          int depth = coords.refIntervals.size() + coords.qryIntervals.size();
          if (depth == 0) {
            auto blk = mergedSet->getBlock(blkId);
            if (blk)
              depth = blk->getSequences().size();
          }
          return depth;
        };

        // 狀況 A：Alignment 只落在單一個 Block 內部
        if (overlaps.size() == 1) {
          BlockID blkId = overlaps[0];
          int pos1 = coordMgr.getLocalPos(globalStart, isRef);
          int pos2 = coordMgr.getLocalPos(globalEnd - 1, isRef);

          int localStart = std::min(pos1, pos2);
          int localEnd = std::max(pos1, pos2) + 1;
          int len = mergedSet->getBlock(blkId)->getConsensus().length();

          bool cutLeft = (localStart > SLICE_SNAP_THRESHOLD);
          bool cutRight = ((len - localEnd) > SLICE_SNAP_THRESHOLD);

          int depth = getBlockDepth(blkId);
          // if (cutLeft && cutRight) Nc += 2 * depth;      // 正中間切兩刀
          // else if (cutLeft || cutRight) Nc += 1 * depth; // 邊緣切一刀
          if (cutLeft && cutRight)
            Nc += 2; // 正中間切兩刀
          else if (cutLeft || cutRight)
            Nc += 1; // 邊緣切一刀

        }
        // 狀況 B：Alignment 跨越多個 Blocks
        else {
          // 1. 檢查左端點 (落在第一個 Block)
          BlockID firstBlk = overlaps.front();
          int p1 = coordMgr.getLocalPos(globalStart, isRef);
          int len1 = mergedSet->getBlock(firstBlk)->getConsensus().length();

          // 如果落點距離頭尾都很遠，代表這個 Block 必定會被切斷
          if (p1 > SLICE_SNAP_THRESHOLD && (len1 - p1) > SLICE_SNAP_THRESHOLD) {
            // Nc += getBlockDepth(firstBlk);
            Nc += 1;
          }

          // 2. 檢查右端點 (落在最後一個 Block)
          BlockID lastBlk = overlaps.back();
          int p2 = coordMgr.getLocalPos(globalEnd - 1, isRef);
          int len2 = mergedSet->getBlock(lastBlk)->getConsensus().length();

          if (p2 > SLICE_SNAP_THRESHOLD && (len2 - p2) > SLICE_SNAP_THRESHOLD) {
            // Nc += getBlockDepth(lastBlk);
            Nc += 1;
          }
        }

        int energy = -L + ALPHA * Nc;

        if (DEBUG_MODE) {
          std::cout << "      📊 [ENERGY-CHECK] " << (isRef ? "Ref" : "Qry")
                    << " | L: " << L << " bp, Nc: " << Nc
                    << " -> Energy: " << energy << "\n";
        }

        if (energy > 0) {
          if (DEBUG_MODE)
            std::cout << "      ❌ [ALIGNMENT REJECTED] High topology energy. "
                         "Fragmentation prevented.\n";
          return true;
        }
        return false;
      };

      // =======================================================
      // 🌟 核心 Lambda：自動處理 Concat 與 雙面 Split
      // =======================================================
      auto extractCoreBlock = [&](int globalStart, int globalEnd, bool isRef, std::vector<BlockID> &overlaps, int &out_left_snap, int &out_right_snap) -> BlockID {
        out_left_snap = 0;
        out_right_snap = 0;
        if (overlaps.empty())
          return 0;

        BlockID targetId = overlaps[0];
        if (overlaps.size() > 1) {
          if (DEBUG_MODE) {
            std::cout << "  🚨 [REJECT] extractCoreBlock: Alignment spans "
                      << overlaps.size() << " blocks. Rejecting merge.\n";
          }
          return 0;
        }

        // 2. 轉換精確 Local 邊界
        global_timer.start("ecb_1_convert_local_pos");
        int pos1 = coordMgr.getLocalPos(globalStart, isRef);
        int pos2 = coordMgr.getLocalPos(globalEnd - 1, isRef);

        int localStart = std::min(pos1, pos2);
        int localEnd = std::max(pos1, pos2) + 1;
        int currentLen = mergedSet->getBlock(targetId)->getConsensus().length();
        global_timer.stop("ecb_1_convert_local_pos");

        // 🌟 綠色通道：如果頭尾都在容忍範圍內，直接放行不切割
        if (localStart <= SLICE_SNAP_THRESHOLD &&
            (currentLen - localEnd) <= SLICE_SNAP_THRESHOLD) {
          out_left_snap = localStart;             // 記錄左邊多吃的
          out_right_snap = currentLen - localEnd; // 記錄右邊多吃的
          return targetId;
        }

        // 🌟 A. 左側防護與切割
        if (localStart > 0 && localStart <= SLICE_SNAP_THRESHOLD) {
          out_left_snap = localStart; // 記錄左邊多吃的
          localStart = 0;
        }

        if (localStart > 0) {
          global_timer.start("ecb_2_left_split");
          auto parts = mergedSet->splitSingleBlock(targetId, localStart);
          global_timer.stop("ecb_2_left_split");
          global_timer.start("ecb_2_left_update");
          coordMgr.updateAfterSplit(targetId, parts.first, parts.second, localStart);
          targetId = parts.second;
          localEnd -= localStart;
          currentLen = mergedSet->getBlock(targetId)->getConsensus().length();
          global_timer.stop("ecb_2_left_update");
        }

        // 🌟 B. 右側防護與切割
        if (localEnd < currentLen &&
            (currentLen - localEnd) <= SLICE_SNAP_THRESHOLD) {
          out_right_snap = currentLen - localEnd; // 記錄右邊多吃的
          localEnd = currentLen;
        }

        if (localEnd < currentLen) {
          global_timer.start("ecb_3_right_split");
          auto parts = mergedSet->splitSingleBlock(targetId, localEnd);
          global_timer.stop("ecb_3_right_split");
          global_timer.start("ecb_3_right_update");
          coordMgr.updateAfterSplit(targetId, parts.first, parts.second,
                                    localEnd);
          targetId = parts.first;
          global_timer.stop("ecb_3_right_update");
        }

        return targetId;
      };

      // ==========================================
      // 🌟 階段 A：獨立能量預檢 (Stateless Energy Pre-check)
      // ==========================================
      global_timer.start("p2_3_energy_precheck");
      bool energy_reject = false;

      // 🛑 核心策略：極速拒絕橫跨多個 Block (overlaps.size() > 1) 的 Alignment
      if (r_overlaps.size() > 1 || q_overlaps.size() > 1) {
        energy_reject = true;
        if (DEBUG_MODE) {
          std::cout << "  🚨 [ALIGNMENT REJECTED] Alignment spans multiple blocks (Ref: "
                    << r_overlaps.size() << ", Qry: " << q_overlaps.size()
                    << "). Bypassing merge.\n";
        }
      }

      if (!energy_reject && rLen > 0 && isAlignmentRejected(r_start, r_end, true, r_overlaps)) {
        energy_reject = true;
      }
      if (!energy_reject && qLen > 0 &&
          isAlignmentRejected(q_start, q_end, false, q_overlaps)) {
        energy_reject = true;
      }
      global_timer.stop("p2_3_energy_precheck");
      if (energy_reject) {
        if (DEBUG_MODE)
          std::cout << "  🚨 [ALIGNMENT REJECTED] Skipping this alignment to "
                       "maintain graph topology.\n";
        continue;
      }

      CigarString precomputed_cigar = bestAln.CIGAR;
      bool precomputed_inverse = bestAln.inverse;

      // 只有 B, C, D 需要重新對齊，我們在這裡先做「字串級」測試
      if (r_merged || q_merged) {
        global_timer.start("p2_4_virtual_realign");

        // 模擬 extractCoreBlock 取得字串，但不改動任何 Pointer 或區間樹
        auto getVirtualStr = [&](int gStart, int gEnd, bool isRef, const std::vector<BlockID> &ovlps) {
          if (ovlps.empty()) return std::string("");
          std::string s = "";
          for (auto id : ovlps)
            s += mergedSet->getBlock(id)->getConsensus().getConsensusString();

          int p1 = coordMgr.getLocalPos(gStart, isRef);
          int p2 = coordMgr.getLocalPos(gEnd - 1, isRef);
          int lStart = std::min(p1, p2);
          int lEnd = std::max(p1, p2) + 1;

          if (lStart <= SLICE_SNAP_THRESHOLD)
            lStart = 0;
          if (lEnd < s.length() && (s.length() - lEnd) <= SLICE_SNAP_THRESHOLD)
            lEnd = s.length();

          if (lStart >= s.length() || lStart >= lEnd)
            return std::string("");
          return s.substr(lStart, lEnd - lStart);
        };

        global_timer.start("p2_4_1_get_virtual_string");
        std::string r_virt_cons = getVirtualStr(r_start, r_end, true, r_overlaps);
        std::string q_virt_cons = getVirtualStr(q_start, q_end, false, q_overlaps);
        global_timer.stop("p2_4_1_get_virtual_string");
        bool is_r_inv = coordMgr.getIsInverse(r_start, true);
        bool is_q_inv = coordMgr.getIsInverse(q_start, false);
        precomputed_inverse = bestAln.inverse;

        if (precomputed_inverse) q_virt_cons = getReverseComplement(q_virt_cons);

        // 🚀 執行預先對齊
        global_timer.start("p2_4_2_semi_global_aln");
        precomputed_cigar = runTilingAlignment(r_virt_cons, q_virt_cons);
        global_timer.stop("p2_4_2_semi_global_aln");
        double identity = identityCIGAR(precomputed_cigar);
        double IDENTITY_THRESHOLD = 0.85;

        global_timer.stop("p2_4_virtual_realign");

        if (identity < IDENTITY_THRESHOLD) {
          if (DEBUG_MODE)
            std::cout << "  🚨 [REJECT] Alignment Identity too low ("
                      << (identity * 100) << "%). Graph untouched. Skipping.\n";
          continue; // 🎯 完美逃脫！因為還沒呼叫
                    // extractCoreBlock，圖的結構毫髮無傷！
        }
      }

      // ==========================================
      // 🌟 階段 B：正式提取與切割 (原封不動的 extractCoreBlock)
      // ==========================================
      global_timer.start("p2_5_extract_core_blocks");
      int r_left_snap = 0, r_right_snap = 0;
      int q_left_snap = 0, q_right_snap = 0;
      if (rLen > 0) {
        rCore = extractCoreBlock(r_start, r_end, true, r_overlaps, r_left_snap, r_right_snap);
        if (DEBUG_MODE) {
          std::cout << "  [REF] Extracting core block " << rCore << "...\n";
          CoordinateManager::GlobalCoords coords = coordMgr.getGlobalCoord(rCore);
          auto printIntervals = [](const std::string &label, const auto &list) {
            std::cout << "         ├─ " << label << ": ";
            if (list.empty())
              std::cout << "None";
            for (auto &iv : list) {
              std::cout << "[" << iv.start << ", " << iv.end << ") " << (iv.isInverse ? "(-)" : "(+)") << " ";
            }
            std::cout << "\n";
          };

          printIntervals("Ref", coords.refIntervals);
          printIntervals("Qry", coords.qryIntervals);
        }
      }

      // ==========================================
      // 🌟 提取 Qry Block
      // ==========================================
      if (qLen > 0) {
        qCore = extractCoreBlock(q_start, q_end, false, q_overlaps, q_left_snap, q_right_snap);
        if (DEBUG_MODE) {
          std::cout << "  [QRY] Extracting core block " << qCore << "...\n";
          CoordinateManager::GlobalCoords coords =
              coordMgr.getGlobalCoord(qCore);

          auto printIntervals = [](const std::string &label, const auto &list) {
            std::cout << "         ├─ " << label << ": ";
            if (list.empty())
              std::cout << "None";
            for (auto &iv : list) {
              std::cout << "[" << iv.start << ", " << iv.end << ") "
                        << (iv.isInverse ? "(-)" : "(+)") << " ";
            }
            std::cout << "\n";
          };

          printIntervals("Ref", coords.refIntervals);
          printIntervals("Qry", coords.qryIntervals);
        }
      }
      global_timer.stop("p2_5_extract_core_blocks");

      // 接下來你就可以拿著乾淨的 rCore 和 qCore 快樂地進行 mergeTwoBlocks 了！

      auto rBlk_ptr = mergedSet->getBlock(rCore);
      auto qBlk_ptr = mergedSet->getBlock(qCore);

      int merge_mode = 0;
      CigarString finalCigar = bestAln.CIGAR;
      bool actual_merge_inverse = bestAln.inverse;

      // ==========================================
      // 🌟 情境判定與 Realignment (B/C/D 全面套用)
      // ==========================================
      // 🌟 拓樸與同序列交叉檢查 (Universal Crossing Check)
      // ==========================================
      bool isCrossing = false;

      // 🛡️ 拓樸錨點交叉檢查 (跳過所有雜訊)
      if (!isCrossing) {
        // 取得四個方向的第一個「已合併」錨點
        auto ref_before =
            coordMgr.getNearestAnchorBlock(r_start, true, false);
        auto ref_after = coordMgr.getNearestAnchorBlock(r_end, true, true);
        auto qry_before =
            coordMgr.getNearestAnchorBlock(q_start, false, false);
        auto qry_after = coordMgr.getNearestAnchorBlock(q_end, false, true);

        if (DEBUG_MODE) {
          std::cout << "Ref: (" << ref_before.first << "," << ref_before.second<< ")\t(" << ref_after.first << "," << ref_after.second<< ")\n";
          std::cout << "Qry: (" << qry_before.first << "," << qry_before.second<< ")\t(" << qry_after.first << "," << qry_after.second<< ")\n";
        }

        if (ref_before.first != -1 && qry_before.first != -1) {
          if (ref_before.first != qry_before.first ||
              ref_before.second != qry_before.second)
            isCrossing = true;
        }
        if (ref_after.first != -1 && qry_after.first != -1) {
          if (ref_after.first != qry_after.first ||
              ref_after.second != qry_after.second)
            isCrossing = true;
        }
      }

      // ==========================================
      // 🌟 情境判定與 Realignment / Patching
      // ==========================================
      if (!r_merged && !q_merged) {
        if (DEBUG_MODE)
          std::cout
              << "  [SCENARIO A] Both ends are new. Extracting blocks...\n";

        // ---------------------------------------------------------
        // 🎯 最終結果分流
        // ---------------------------------------------------------
        if (!isCrossing) {
          if (DEBUG_MODE)
            std::cout << "      [SCENARIO A1] Collinear paths.\n";
          merge_mode = 1;
        } else {
          if (DEBUG_MODE)
            std::cout
                << "      [SCENARIO A2] Crossing or Duplication detected!\n";
          merge_mode = 2;
        }
        // =========================================================
        // 🌟 CIGAR 邊緣補綴 (Flank Patching for Snapped Regions)
        // =========================================================
        global_timer.start("p2_6_case_A_flank_patching");
        std::string rCons = rBlk_ptr->getConsensus().getConsensusString();
        std::string qCons = qBlk_ptr->getConsensus().getConsensusString();

        bool is_r_inv = coordMgr.getIsInverse(r_start, true);
        bool is_q_inv = coordMgr.getIsInverse(q_start, false);
        actual_merge_inverse = bestAln.inverse;

        // 計算 CIGAR 前端與後端各自需要補多少 bp
        // 如果 Qry 是反向的，Qry 的「實體右側」要補到 CIGAR 的「左端」
        int c_start_r = r_left_snap;
        int c_start_q = actual_merge_inverse ? q_right_snap : q_left_snap;

        int c_end_r = r_right_snap;
        int c_end_q = actual_merge_inverse ? q_left_snap : q_right_snap;

        // 輔助函數：決定要跑 Global Aln 還是直接補 Gap
        auto alignFlanks = [&](const std::string &r_seq,
                               const std::string &q_seq) -> CigarString {
          if (r_seq.empty() && q_seq.empty())
            return {};
          if (r_seq.empty())
            return {{(int)q_seq.length(), 'I'}}; // Qry 獨有 -> Insertion
          if (q_seq.empty())
            return {{(int)r_seq.length(), 'D'}}; // Ref 獨有 -> Deletion
          return runGlobalAlignment(
              r_seq, q_seq); // 兩邊都有 -> 呼叫你的 Global Alignment
        };

        // 1. 取得前端字串 (若 Qry 反向，記得取 RC)
        std::string r_start_seq = rCons.substr(0, c_start_r);
        std::string q_start_seq =
            actual_merge_inverse ? getReverseComplement(qCons.substr(
                                       qCons.length() - c_start_q, c_start_q))
                                 : qCons.substr(0, c_start_q);

        // 2. 取得後端字串 (若 Qry 反向，記得取 RC)
        std::string r_end_seq = rCons.substr(rCons.length() - c_end_r, c_end_r);
        std::string q_end_seq =
            actual_merge_inverse
                ? getReverseComplement(qCons.substr(0, c_end_q))
                : qCons.substr(qCons.length() - c_end_q, c_end_q);

        // 3. 計算外掛的 CIGAR
        CigarString start_cigar = alignFlanks(r_start_seq, q_start_seq);
        CigarString end_cigar = alignFlanks(r_end_seq, q_end_seq);

        // 4. 串接起來
        CigarString patchedCigar;
        patchedCigar.insert(patchedCigar.end(), start_cigar.begin(),
                            start_cigar.end());
        patchedCigar.insert(patchedCigar.end(), bestAln.CIGAR.begin(),
                            bestAln.CIGAR.end());
        patchedCigar.insert(patchedCigar.end(), end_cigar.begin(),
                            end_cigar.end());

        finalCigar = compressCigar(patchedCigar);
        global_timer.stop("p2_6_case_A_flank_patching");

        if (DEBUG_MODE) {
          std::cout << "      -> Patched CIGAR for Case A: ";
          printCIGAR(finalCigar);
          std::cout << "\n";
        }
      } else {
        // 情境 A3 或 B, C, D (全部觸發 Realignment)
        if (!isCrossing) {
          if (DEBUG_MODE)
            std::cout << "  [SCENARIO A3] Orthologous merge with pre-merged block (Realignment applied).\n";
          merge_mode = 1;
        } else if (r_merged && !q_merged) {
          if (DEBUG_MODE)
            std::cout << "  [SCENARIO B] Ref exists, Qry new.\n";
          merge_mode = 3;
        } else if (!r_merged && q_merged) {
          if (DEBUG_MODE)
            std::cout << "  [SCENARIO C] Qry exists, Ref new.\n";
          merge_mode = 4;
        } else {
          if (DEBUG_MODE) {
            std::cout << "  [SCENARIO D] Both exist.\n";
            std::cout << rBlk_ptr->getConsensus().getConsensusString() << '\n';
            std::cout << qBlk_ptr->getConsensus().getConsensusString() << '\n';
          }
          merge_mode = 5;
        }

        if (DEBUG_MODE)
          std::cout << "      -> Re-aligning Consensus: Ref("
                    << rBlk_ptr->getConsensus().size() << "bp) vs Qry("
                    << qBlk_ptr->getConsensus().size() << "bp)...\n";

        actual_merge_inverse = precomputed_inverse;
        finalCigar = precomputed_cigar;

        double identity = identityCIGAR(finalCigar);

        double IDENTITY_THRESHOLD = 0.85;

        if (DEBUG_MODE) {
          std::cout << "      -> Re-aligned CIGAR: ";
          printCIGAR(finalCigar);
          std::cout << "\n";
          std::cout << "      -> Alignment Identity: " << (identity * 100.0)
                    << "%\n";
        }
      }

      // 1. 🌟 [Merge 前] 抓取真正對齊的代表性 Segment (探針)
      auto getProbeSegment = [](BlockPtr blk, int start_pos) {
        for (auto& seqPair : blk->getSequences()) {
          for (auto& segPair : seqPair.second.getSegments()) {
            if (segPair.first <= start_pos && segPair.second.getEnd() > start_pos) {
              return std::make_pair(seqPair.first, segPair.second);
            }
          }
        }
        auto it = blk->getSequences().begin();
        return std::make_pair(it->first, it->second.getSegments().begin()->second);
      };

      auto r_probe = getProbeSegment(rBlk_ptr, r_start);
      std::string r_probe_seq = r_probe.first;
      int r_probe_start = r_probe.second.getStart();
      int r_probe_copy = r_probe.second.getCopyCount();

      auto q_probe = getProbeSegment(qBlk_ptr, q_start);
      std::string q_probe_seq = q_probe.first;
      int q_probe_start = q_probe.second.getStart();
      int q_probe_copy = q_probe.second.getCopyCount();

      global_timer.start("merge_blocks");
      auto mBlk = mergedSet->mergeTwoBlocks(rBlk_ptr, qBlk_ptr, finalCigar,
                                            actual_merge_inverse, merge_mode,
                                            r_probe_start, q_probe_start,
                                            r_probe_copy, q_probe_copy);
      global_timer.stop("merge_blocks");

      // 3. 🌟 [Merge 後] 在新的 mBlk 中找回那兩個探針，看它們的 Copy 變成多少
      int r_new_copy = mBlk->getSequences()[r_probe_seq].getSegment(r_probe_start).getCopyCount();
      int q_new_copy = mBlk->getSequences()[q_probe_seq].getSegment(q_probe_start).getCopyCount();

      int actual_rLen = rBlk_ptr->getConsensus().length();
      int actual_qLen = qBlk_ptr->getConsensus().length();

      // 4. 🎯 更新 CoordinateManager
      global_timer.start("p2_8_coord_tracker_update");
      bool inBoth = (merge_mode == 1 || r_new_copy == q_new_copy);
      coordMgr.updateAfterMerge(rCore, qCore, mBlk->getId(), finalCigar,
                                actual_merge_inverse, actual_rLen, actual_qLen,
                                r_new_copy, q_new_copy, inBoth,
                                r_probe_start, q_probe_start);

      BlockID mId = mBlk->getId();
      if (DEBUG_MODE)
        std::cout << "      -> Block " << mId << " created/merged (Length: "
                  << mBlk->getConsensus().length() << ")\n";
                  

      // 🌟 [Case D 特殊列印] 當為 Case D (merge_mode == 5) 且 merge 完後，整齊印出 mBlk 底下所有 Segment 資訊與 Variants
      if (merge_mode == 5 && DEBUG_MODE) {
        std::cout << "================================================================================\n";
        std::cout << "  [Case D] Merged Block ID: " << mBlk->getId()
                  << " (Consensus Length: " << mBlk->getConsensus().length() << " bp)\n";
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << "  Sequences & Segments Detail:\n";

        auto &seqs = mBlk->getSequences();
        std::vector<std::string> seq_ids;
        for (const auto &p : seqs) {
          seq_ids.push_back(p.first);
        }
        std::sort(seq_ids.begin(), seq_ids.end());

        for (size_t i = 0; i < seq_ids.size(); ++i) {
          const std::string &seq_id = seq_ids[i];
          bool is_last_seq = (i == seq_ids.size() - 1);
          std::string seq_branch = is_last_seq ? "  └─ " : "  ├─ ";
          std::string seq_indent = is_last_seq ? "     " : "  │  ";

          std::cout << seq_branch << "Sequence: " << seq_id << "\n";

          auto &segments = seqs.at(seq_id).getSegments();
          size_t seg_idx = 0;
          for (auto &seg_pair : segments) {
            bool is_last_seg = (++seg_idx == segments.size());
            std::string seg_branch = is_last_seg ? "└─ " : "├─ ";
            std::string seg_indent = is_last_seg ? "   " : "│  ";

            auto &seg = seg_pair.second;
            std::cout << seq_indent << seg_branch << "Segment: [" << seg.getStart() << " -> " << seg.getEnd() << "]"
                      << " | Strand: " << (seg.isReverse() ? "-" : "+")
                      << " | Copy: " << seg.getCopyCount() << "\n";

            auto &variants = seg.getVariants();
            if (variants.empty()) {
              std::cout << seq_indent << seg_indent << "└─ Variants: (none)\n";
            } else {
              std::cout << seq_indent << seg_indent << "├─ Variants (" << variants.size() << "):\n";
              for (size_t v_idx = 0; v_idx < variants.size(); ++v_idx) {
                auto &var = variants[v_idx];
                bool is_last_var = (v_idx == variants.size() - 1);
                std::string var_branch = is_last_var ? "└─ " : "├─ ";

                std::cout << seq_indent << seg_indent << "│  " << var_branch;
                if (var.getType() == VariantType::SNV) {
                  std::cout << "SNV pos: " << var.getStart() << ", alt: '" << var.getAlt() << "'\n";
                } else if (var.getType() == VariantType::GAP) {
                  std::cout << "GAP range: [" << var.getStart() << " -> " << var.getEnd() << "]\n";
                } else {
                  std::cout << "Variant range: [" << var.getStart() << " -> " << var.getEnd() << "]\n";
                }
              }
            }
          }
        }
        std::cout << "================================================================================\n";
      }

      
      // 🌟 將 Ref 和 Qry 的 Core 從 Set 中刪除，因為已經合體為 mId 了
      mergedSet->deleteBlock(rCore);
      mergedSet->deleteBlock(qCore);

      mergeCounter++;

      alnCollection.ref_coverageTracker.syncFromMap(coordMgr, true); // true 代表 Ref
      alnCollection.qry_coverageTracker.syncFromMap(coordMgr, false); // false 代表 Qry
      global_timer.stop("p2_8_coord_tracker_update");
      // if (DEBUG_MODE) mergedSet->debugValidateSequences(this);
    }
  }
  global_timer.stop("phase2_iterative_merge");

  if (DEBUG_MODE)
    std::cout << "\n  -> Processed " << mergeCounter << " valid alignments.\n";

  // ==========================================
  // 🌟 Phase 5: 拓撲重建與家族歸一化
  // ==========================================
  global_timer.start("phase5_rewire_edges");
  if (DEBUG_MODE)
    std::cout << "\n[Phase 5] Re-wiring Linear Pangenome Graph Edges...\n";
  for (auto &block : mergedSet->getAllBlocks()) {
    auto blk = block.lock();
    if (!blk)
      continue;
    blk->normalizeStrand();
  }

  for (auto &seq : refSet->getSequences())
    mergedSet->addSequenceName(seq);
  for (auto &seq : qrySet->getSequences())
    mergedSet->addSequenceName(seq);

  for (auto &blk : refSet->getAllBlocks()) {
    if (blk.lock()->isDistant()) {
      if (DEBUG_MODE) {
        auto b = blk.lock();
        std::cout << "  [DISTANT-REF] Adding Block " << b->getId()
                  << " (len=" << b->getConsensus().length() << ")";
        for (auto &seq : b->getSequences()) {
          for (auto &seg : seq.second.getSegments()) {
            std::cout << " | " << seq.first << ":[" << seg.second.getStart()
                      << "," << seg.second.getEnd() << ")";
          }
        }
        std::cout << "\n";
      }
      mergedSet->addBlock(blk.lock());
    }
  }
  for (auto &blk : qrySet->getAllBlocks()) {
    if (blk.lock()->isDistant()) {
      if (DEBUG_MODE) {
        auto b = blk.lock();
        std::cout << "  [DISTANT-QRY] Adding Block " << b->getId()
                  << " (len=" << b->getConsensus().length() << ")";
        for (auto &seq : b->getSequences()) {
          for (auto &seg : seq.second.getSegments()) {
            std::cout << " | " << seq.first << ":[" << seg.second.getStart()
                      << "," << seg.second.getEnd() << ")";
          }
        }
        std::cout << "\n";
      }
      mergedSet->addBlock(blk.lock());
    }
  }

  // 🔍 Diagnostic checkpoint B: after distant block additions
  if (DEBUG_MODE) std::cout << "\n  [DIAG-B] Validating AFTER distant block additions...\n";
  mergedSet->debugValidateSequences(this);

  mergedSet->rebuildLinearGraph(coordMgr);
  global_timer.stop("phase5_rewire_edges");

  // 🌟 [Phase 6] 恢復實體 Consensus 字串，並透過 Majority Voting 更新 Consensus 與 SNV
  global_timer.start("phase6_recover_consensus");
  if (DEBUG_MODE)
    std::cout << "\n[Phase 6] Recovering Consensus Strings & Majority Base Voting...\n";
  for (auto &block : mergedSet->getAllBlocks()) {
    auto blk = block.lock();
    if (!blk) continue;
    blk->refineConsensusAndVariants();
  }
  mergedSet->debugValidateSequences(this);
  global_timer.stop("phase6_recover_consensus");

  auto timeEnd = std::chrono::high_resolution_clock::now();
  if (DEBUG_MODE) {
    std::cout << "\n========================================================\n"
              << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
              << "Total Execution Time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd -
                                                                       time0)
                     .count()
              << " ms\n"
              // << "Find Best Alignment:  " << findBest << " ms\n"
              // << "Merge Blocks:         " << merge_time << " ms\n"
              << "========================================================\n\n";
  }

  // mergedSet->debugValidateLinkages(true);
  // mergedSet->debugValidateSequences(this);

  global_timer.print();

  return mergedSet;
}