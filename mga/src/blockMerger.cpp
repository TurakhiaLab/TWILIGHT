#include "alignment.hpp"
#include "block_manager.hpp"
#include "coordinate_manager.hpp"
#include "global_alignment.hpp"
#include "timer.hpp"
#include "type.hpp"
#include "cigar_util.hpp"
#include "merge_debug.hpp"

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

namespace {

/**
 * @brief Evaluates whether adding the alignment would induce excessive graph fragmentation or cross block occurrences.
 */
bool isAlignmentRejected(int globalStart, int globalEnd, bool isRef,
                         const std::vector<BlockID> &overlaps,
                         CoordinateManager &coordMgr, BlockSet *mergedSet,
                         mga::MergeDebugger &dbg) {
  if (overlaps.empty())
    return false;

  const int ALPHA = MERGE_ALPHA;
  int L = globalEnd - globalStart;
  int Nc = 0;

  // Helper: safely get Block depth (total number of Ref and Qry tracks)
  auto getBlockDepth = [&](BlockID blkId) {
    CoordinateManager::GlobalCoords coords = coordMgr.getGlobalCoord(blkId);
    int depth = coords.refIntervals.size() + coords.qryIntervals.size();
    if (depth == 0) {
      auto blk = mergedSet->getBlock(blkId);
      if (blk)
        depth = blk->getSequences().size();
    }
    return depth;
  };

  // Case A: Alignment resides entirely within a single Block
  if (overlaps.size() == 1) {
    BlockID blkId = overlaps[0];
    int pos1 = coordMgr.getLocalPos(globalStart, isRef);
    int pos2 = coordMgr.getLocalPos(globalEnd - 1, isRef);

    // Guard: overlaps.size() == 1 may still span two different occurrences of the same Block.
    // For the same occurrence, abs(localPos difference) must equal abs(globalPos difference)
    if (std::abs(pos2 - pos1) != (globalEnd - 1 - globalStart)) {
      dbg.logOccurrenceCrossingReject(false);
      return true;
    }

    int localStart = std::min(pos1, pos2);
    int localEnd = std::max(pos1, pos2) + 1;
    int len = mergedSet->getBlock(blkId)->getConsensus().length();

    bool cutLeft = (localStart > SLICE_SNAP_THRESHOLD);
    bool cutRight = ((len - localEnd) > SLICE_SNAP_THRESHOLD);

    int depth = getBlockDepth(blkId);
    if (cutLeft && cutRight)
      Nc += 2; // Split twice in the middle
    else if (cutLeft || cutRight)
      Nc += 1; // Split once at edge

  }
  // Case B: Alignment spans across multiple Blocks
  else {
    // 1. Check left endpoint (within first Block)
    BlockID firstBlk = overlaps.front();
    int p1 = coordMgr.getLocalPos(globalStart, isRef);
    int len1 = mergedSet->getBlock(firstBlk)->getConsensus().length();

    // If position is far from both ends, this Block will definitely be split
    if (p1 > SLICE_SNAP_THRESHOLD && (len1 - p1) > SLICE_SNAP_THRESHOLD) {
      // Nc += getBlockDepth(firstBlk);
      Nc += 1;
    }

    // 2. Check right endpoint (within last Block)
    BlockID lastBlk = overlaps.back();
    int p2 = coordMgr.getLocalPos(globalEnd - 1, isRef);
    int len2 = mergedSet->getBlock(lastBlk)->getConsensus().length();

    if (p2 > SLICE_SNAP_THRESHOLD && (len2 - p2) > SLICE_SNAP_THRESHOLD) {
      // Nc += getBlockDepth(lastBlk);
      Nc += 1;
    }
  }

  int energy = -L + ALPHA * Nc;
  dbg.logEnergyCalculation(isRef, L, Nc, energy);

  if (energy > 0) {
    dbg.logEnergyReject();
    return true;
  }
  return false;
}

/**
 * @brief Splits blocks if necessary and extracts clean core block for merging.
 */
BlockID extractCoreBlock(int globalStart, int globalEnd, bool isRef,
                         std::vector<BlockID> &overlaps,
                         int &out_left_snap, int &out_right_snap,
                         CoordinateManager &coordMgr, BlockSet *mergedSet,
                         mga::MergeDebugger &dbg, Timer &timer = global_timer) {
  out_left_snap = 0;
  out_right_snap = 0;
  if (overlaps.empty())
    return 0;

  BlockID targetId = overlaps[0];
  if (overlaps.size() > 1) {
    dbg.logMultiBlockSpanningReject(overlaps.size(), 0);
    return 0;
  }

  // Convert to exact local coordinates
  timer.start("ecb_1_convert_local_pos");
  int pos1 = coordMgr.getLocalPos(globalStart, isRef);
  int pos2 = coordMgr.getLocalPos(globalEnd - 1, isRef);

  // Guard: ensure coordinates belong to the same occurrence
  if (std::abs(pos2 - pos1) != (globalEnd - 1 - globalStart)) {
    dbg.logOccurrenceCrossingReject(true);
    timer.stop("ecb_1_convert_local_pos");
    return 0;
  }

  int localStart = std::min(pos1, pos2);
  int localEnd = std::max(pos1, pos2) + 1;
  int currentLen = mergedSet->getBlock(targetId)->getConsensus().length();
  timer.stop("ecb_1_convert_local_pos");

  // Fast path: if both boundaries are within snap tolerance, proceed without splitting
  if (localStart <= SLICE_SNAP_THRESHOLD &&
      (currentLen - localEnd) <= SLICE_SNAP_THRESHOLD) {
    out_left_snap = localStart;             // Record left snapped overhang
    out_right_snap = currentLen - localEnd; // Record right snapped overhang
    return targetId;
  }

  // Left and right snap guards
  if (localStart > 0 && localStart <= SLICE_SNAP_THRESHOLD) {
    out_left_snap = localStart; // Record left snapped overhang
    localStart = 0;
  }

  if (localEnd < currentLen &&
      (currentLen - localEnd) <= SLICE_SNAP_THRESHOLD) {
    out_right_snap = currentLen - localEnd; // Record right snapped overhang
    localEnd = currentLen;
  }

  // Determine single split or single-pass double split
  if (localStart > 0 && localEnd < currentLen) {
    timer.start("ecb_double_split");
    auto parts = mergedSet->splitDoubleBlock(targetId, localStart, localEnd);
    timer.stop("ecb_double_split");
    timer.start("ecb_double_update");
    coordMgr.updateAfterDoubleSplit(targetId, parts.leftID, parts.midID, parts.rightID, localStart, localEnd);
    targetId = parts.midID;
    timer.stop("ecb_double_update");
  } else if (localStart > 0) {
    timer.start("ecb_2_left_split");
    auto parts = mergedSet->splitSingleBlock(targetId, localStart);
    timer.stop("ecb_2_left_split");
    timer.start("ecb_2_left_update");
    coordMgr.updateAfterSplit(targetId, parts.first, parts.second, localStart);
    targetId = parts.second;
    timer.stop("ecb_2_left_update");
  } else if (localEnd < currentLen) {
    timer.start("ecb_3_right_split");
    auto parts = mergedSet->splitSingleBlock(targetId, localEnd);
    timer.stop("ecb_3_right_split");
    timer.start("ecb_3_right_update");
    coordMgr.updateAfterSplit(targetId, parts.first, parts.second, localEnd);
    targetId = parts.first;
    timer.stop("ecb_3_right_update");
  }

  return targetId;
}

/**
 * @brief Simulates extractCoreBlock to retrieve consensus substring without mutating the graph.
 */
std::string getVirtualConsensusString(int gStart, int gEnd, bool isRef,
                                      const std::vector<BlockID> &ovlps,
                                      CoordinateManager &coordMgr,
                                      BlockSet *mergedSet) {
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
}

/**
 * @brief Aligns overhang flanks using global alignment or pure insertion/deletion.
 */
CigarString alignFlanks(const std::string &r_seq, const std::string &q_seq) {
  if (r_seq.empty() && q_seq.empty())
    return {};
  if (r_seq.empty())
    return {{(int)q_seq.length(), 'I'}}; // Query-only -> Insertion
  if (q_seq.empty())
    return {{(int)r_seq.length(), 'D'}}; // Reference-only -> Deletion
  return runGlobalAlignment(r_seq, q_seq); // Both present -> invoke global alignment
}

} // anonymous namespace

BlockSet *BlockManager::merge(BlockSet *refSet, BlockSet *qrySet,
                              AlignmentCollection &alnCollection,
                              BlockSetID newID, Tree *tree, Option *option, int L_min) {
  // Initialize merge debugger with customizable preset or boolean flag
  bool DEBUG_MODE = false;
  mga::MergeDebugger dbg(DEBUG_MODE ? mga::MergeDebugPreset::FULL : mga::MergeDebugPreset::OFF);
  auto time0 = std::chrono::high_resolution_clock::now();
  
  dbg.logMergeHeader(refSet->getId(), qrySet->getId());

  Timer local_timer;

  local_timer.start("phase1_init_coord_mgr");
  dbg.logPhase(1, "Initialized Coordinate Manager...");

  BlockSet *mergedSet = createBlockSet(newID);

  refSet->setTree(tree);
  qrySet->setTree(tree);
  mergedSet->setTree(tree);

  CoordinateManager coordMgr;
  coordMgr.init(refSet, qrySet, mergedSet);

  alnCollection.ref_coverageTracker.syncFromMap(coordMgr, true);  // true for Ref
  alnCollection.qry_coverageTracker.syncFromMap(coordMgr, false); // false for Qry

  int mergeCounter = 0;
  local_timer.stop("phase1_init_coord_mgr");

  dbg.logPhase(2, "Iteratively merged blocks based on Minimap2 alignments...");

  for (auto &seq : refSet->getSequences()) mergedSet->addSequenceName(seq);
  for (auto &seq : qrySet->getSequences()) mergedSet->addSequenceName(seq);

  // Call BlockSet::createIndexedCopy to generate indexed BlockSet retaining BlockID with zero large-string copy cost
  BlockSetID refIndexID = refSet->getId() + "_indexed";
  BlockSetID qryIndexID = qrySet->getId() + "_indexed";
  BlockSet *refIndexSet = refSet->createIndexedCopy(this, refIndexID);
  BlockSet *qryIndexSet = qrySet->createIndexedCopy(this, qryIndexID);

  local_timer.start("phase2_iterative_merge");
  while (true) {

    local_timer.start("find_best_alignment");
    Alignments bestAlns = alnCollection.getBestAlignments(refIndexSet, qryIndexSet, coordMgr);
    local_timer.stop("find_best_alignment");

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
      local_timer.start("p2_2_overlap_lookup");
      auto r_overlaps = coordMgr.getOverlapBlocks(r_start, r_end, true, r_merged);
      auto q_overlaps = coordMgr.getOverlapBlocks(q_start, q_end, false, q_merged);
      local_timer.stop("p2_2_overlap_lookup");

      dbg.logAlignmentDetails(bestAln.CIGAR, bestAln.inverse, r_start, r_end, q_start, q_end, r_merged, q_merged, r_overlaps, q_overlaps);

      if (!r_overlaps.empty() && !q_overlaps.empty() && r_overlaps.size() == 1 && q_overlaps.size() == 1 && r_overlaps[0] == q_overlaps[0]) {
        dbg.logSkipIdenticalBlock(r_overlaps[0]);
        continue;
      }

      int rLen = r_end - r_start;
      int qLen = q_end - q_start;

      if (rLen == 0 && qLen == 0) {
        dbg.logNullAlignment();
        continue;
      }

      // =======================================================
      // Pure Indel Fast Path
      // =======================================================
      if ((is_pure_deletion && rLen > 0) || (is_pure_insertion && qLen > 0)) {
        dbg.logPureIndelBypass();
        continue;
      }

      const int min_threshold = MERGE_MIN_THRESHOLD;
      if (!bestAln.primary && !(is_pure_deletion || is_pure_insertion)) {
        int score = scoreCIGAR(bestAln.CIGAR);
        if (score < min_threshold) {
          dbg.logSecondaryScoreSkip(score, min_threshold);
          continue;
        }
      }

      BlockID rCore = 0, qCore = 0;

      // ==========================================
      // Phase A: Stateless Energy Pre-check
      // ==========================================
      local_timer.start("p2_3_energy_precheck");
      bool energy_reject = false;

      // Policy: rapidly reject alignments spanning multiple blocks (overlaps.size() > 1)
      if (r_overlaps.size() > 1 || q_overlaps.size() > 1) {
        energy_reject = true;
        dbg.logMultiBlockSpanningReject(r_overlaps.size(), q_overlaps.size());
      }

      if (!energy_reject && rLen > 0 && isAlignmentRejected(r_start, r_end, true, r_overlaps, coordMgr, mergedSet, dbg)) {
        energy_reject = true;
      }
      if (!energy_reject && qLen > 0 &&
          isAlignmentRejected(q_start, q_end, false, q_overlaps, coordMgr, mergedSet, dbg)) {
        energy_reject = true;
      }
      local_timer.stop("p2_3_energy_precheck");
      if (energy_reject) {
        dbg.logTopologySkip();
        continue;
      }

      CigarString precomputed_cigar = bestAln.CIGAR;
      bool precomputed_inverse = bestAln.inverse;

      // Only cases B, C, D require realignment; perform string-level pre-test here
      if (r_merged || q_merged) {
        local_timer.start("p2_4_virtual_realign");

        local_timer.start("p2_4_1_get_virtual_string");
        std::string r_virt_cons = getVirtualConsensusString(r_start, r_end, true, r_overlaps, coordMgr, mergedSet);
        std::string q_virt_cons = getVirtualConsensusString(q_start, q_end, false, q_overlaps, coordMgr, mergedSet);
        local_timer.stop("p2_4_1_get_virtual_string");
        bool is_r_inv = coordMgr.getIsInverse(r_start, true);
        bool is_q_inv = coordMgr.getIsInverse(q_start, false);
        precomputed_inverse = bestAln.inverse;

        if (precomputed_inverse) q_virt_cons = getReverseComplement(q_virt_cons);

        // Execute pre-alignment (Linear Gap Penalty + Tiling Alignment)
        local_timer.start("p2_4_2_semi_global_aln");
        precomputed_cigar = runTilingAlignmentLinearGap(r_virt_cons, q_virt_cons);
        local_timer.stop("p2_4_2_semi_global_aln");
        double identity = identityCIGAR(precomputed_cigar);
        double IDENTITY_THRESHOLD = 0.85;

        local_timer.stop("p2_4_virtual_realign");

        if (identity < IDENTITY_THRESHOLD) {
          dbg.logLowIdentityReject(identity);
          continue; // Early exit: extractCoreBlock has not been invoked yet, graph structure remains intact
        }
      }

      // ==========================================
      // Phase B: Formal extraction and splitting via extractCoreBlock
      // ==========================================
      local_timer.start("p2_5_extract_core_blocks");
      int r_left_snap = 0, r_right_snap = 0;
      int q_left_snap = 0, q_right_snap = 0;
      if (rLen > 0) {
        rCore = extractCoreBlock(r_start, r_end, true, r_overlaps, r_left_snap, r_right_snap, coordMgr, mergedSet, dbg, local_timer);
        dbg.logCoreExtraction("REF", rCore, coordMgr.getGlobalCoord(rCore));
      }

      // ==========================================
      // Extract Query Core Block
      // ==========================================
      if (qLen > 0) {
        qCore = extractCoreBlock(q_start, q_end, false, q_overlaps, q_left_snap, q_right_snap, coordMgr, mergedSet, dbg, local_timer);
        dbg.logCoreExtraction("QRY", qCore, coordMgr.getGlobalCoord(qCore));
      }
      local_timer.stop("p2_5_extract_core_blocks");

      // Clean rCore and qCore are now prepared for mergeTwoBlocks

      auto rBlk_ptr = mergedSet->getBlock(rCore);
      auto qBlk_ptr = mergedSet->getBlock(qCore);

      if (rCore == 0 || qCore == 0 || rCore == qCore || !rBlk_ptr || !qBlk_ptr) {
        dbg.logInvalidCoreBlocks(rCore, qCore);
        continue;
      }

      int merge_mode = 0;
      CigarString finalCigar = bestAln.CIGAR;
      bool actual_merge_inverse = bestAln.inverse;

      // ==========================================
      // Scenario Classification and Realignment (applied to B/C/D)
      // ==========================================
      // Universal Crossing Check (topology and shared sequence check)
      // ==========================================
      bool isCrossing = false;

      // Topology anchor crossing check
      if (!isCrossing) {
        // Retrieve first merged anchor block in four directions
        auto ref_before =
            coordMgr.getNearestAnchorBlock(r_start, true, false);
        auto ref_after = coordMgr.getNearestAnchorBlock(r_end, true, true);
        auto qry_before =
            coordMgr.getNearestAnchorBlock(q_start, false, false);
        auto qry_after = coordMgr.getNearestAnchorBlock(q_end, false, true);

        dbg.logAnchorNeighbors(ref_before, ref_after, qry_before, qry_after);

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

      // Check if Ref and Qry share identical sequence names (indicates intra-genome duplication/paralogy, preventing collinear merge)
      bool hasCommonSeq = false;
      if (rBlk_ptr && qBlk_ptr) {
        for (const auto &kv : rBlk_ptr->getSequences()) {
          if (qBlk_ptr->getSequences().count(kv.first)) {
            hasCommonSeq = true;
            break;
          }
        }
      }
      if (hasCommonSeq) {
        isCrossing = true;
      }

      // ==========================================
      // Scenario classification and realignment / patching
      // ==========================================
      if (!r_merged && !q_merged) {
        dbg.logScenario("SCENARIO A", "Both ends are new. Extracting blocks...");

        // ---------------------------------------------------------
        // Dispatch based on scenario
        // ---------------------------------------------------------
        if (!isCrossing) {
          dbg.logScenario("SCENARIO A1", "Collinear paths.");
          merge_mode = 1;
        } else {
          dbg.logScenario("SCENARIO A2", "Crossing or Duplication detected!");
          merge_mode = 2;
        }
        // =========================================================
        // Flank patching for snapped regions on CIGAR edges
        // =========================================================
        local_timer.start("p2_6_case_A_flank_patching");
        std::string rCons = rBlk_ptr->getConsensus().getConsensusString();
        std::string qCons = qBlk_ptr->getConsensus().getConsensusString();

        bool is_r_inv = coordMgr.getIsInverse(r_start, true);
        bool is_q_inv = coordMgr.getIsInverse(q_start, false);
        actual_merge_inverse = bestAln.inverse;

        // Compute number of bp to patch at head and tail of CIGAR
        // If Qry is inverted, Qry physical right-side patches to CIGAR left-end
        int c_start_r = r_left_snap;
        int c_start_q = actual_merge_inverse ? q_right_snap : q_left_snap;

        int c_end_r = r_right_snap;
        int c_end_q = actual_merge_inverse ? q_left_snap : q_right_snap;

        // 1. Retrieve head sequences (take reverse complement if Qry is inverted)
        std::string r_start_seq = rCons.substr(0, c_start_r);
        std::string q_start_seq =
            actual_merge_inverse ? getReverseComplement(qCons.substr(
                                       qCons.length() - c_start_q, c_start_q))
                                 : qCons.substr(0, c_start_q);

        // 2. Retrieve tail sequences (take reverse complement if Qry is inverted)
        std::string r_end_seq = rCons.substr(rCons.length() - c_end_r, c_end_r);
        std::string q_end_seq =
            actual_merge_inverse
                ? getReverseComplement(qCons.substr(0, c_end_q))
                : qCons.substr(qCons.length() - c_end_q, c_end_q);

        // 3. Compute patched CIGAR flanks
        CigarString start_cigar = alignFlanks(r_start_seq, q_start_seq);
        CigarString end_cigar = alignFlanks(r_end_seq, q_end_seq);

        // 4. Concatenate patched CIGAR
        CigarString patchedCigar;
        patchedCigar.insert(patchedCigar.end(), start_cigar.begin(),
                            start_cigar.end());
        patchedCigar.insert(patchedCigar.end(), bestAln.CIGAR.begin(),
                            bestAln.CIGAR.end());
        patchedCigar.insert(patchedCigar.end(), end_cigar.begin(),
                            end_cigar.end());

        finalCigar = compressCigar(patchedCigar);
        local_timer.stop("p2_6_case_A_flank_patching");

        dbg.logCaseAPatchedCigar(finalCigar);
      } else {
        // Scenario A3 or B, C, D (all trigger realignment)
        if (!isCrossing && !hasCommonSeq) {
          dbg.logScenario("SCENARIO A3", "Orthologous merge with pre-merged block (Realignment applied).");
          merge_mode = 1;
        } else if (r_merged && !q_merged) {
          dbg.logScenario("SCENARIO B", "Ref exists, Qry new.");
          merge_mode = 3;
        } else if (!r_merged && q_merged) {
          dbg.logScenario("SCENARIO C", "Qry exists, Ref new.");
          merge_mode = 4;
        } else {
          dbg.logScenario("SCENARIO D", "Both exist.");
          dbg.logCaseDConsensus(rBlk_ptr->getConsensus().getConsensusString(),
                                qBlk_ptr->getConsensus().getConsensusString());
          merge_mode = 5;
        }

        dbg.logRealigning(rBlk_ptr->getConsensus().size(), qBlk_ptr->getConsensus().size());

        actual_merge_inverse = precomputed_inverse;
        finalCigar = precomputed_cigar;

        double identity = identityCIGAR(finalCigar);
        dbg.logRealignedCigar(finalCigar, identity);
      }

      // 1. [Pre-merge] Retrieve actual Copy IDs corresponding to this alignment from CoordinateManager
      int r_probe_copy = coordMgr.getCopyId(r_start, true);
      int q_probe_copy = coordMgr.getCopyId(q_start, false);

      int maxRefCopy = rBlk_ptr ? rBlk_ptr->getMaxCopy() : 0;
      int maxQryCopy = qBlk_ptr ? qBlk_ptr->getMaxCopy() : 0;
      int actual_rLen = rBlk_ptr ? rBlk_ptr->getConsensus().length() : 0;
      int actual_qLen = qBlk_ptr ? qBlk_ptr->getConsensus().length() : 0;

      local_timer.start("merge_blocks");
      auto mBlk = mergedSet->mergeTwoBlocks(rBlk_ptr, qBlk_ptr, finalCigar,
                                            actual_merge_inverse, merge_mode,
                                            -1, -1,
                                            r_probe_copy, q_probe_copy);
      local_timer.stop("merge_blocks");

      // 3. [Post-merge] Compute resulting Copy IDs in the new merged block
      int r_new_copy = r_probe_copy;
      int q_new_copy = q_probe_copy;

      if (merge_mode == 1) {
        q_new_copy = r_probe_copy;
      } else if (merge_mode == 2 || merge_mode == 3) {
        q_new_copy = q_probe_copy + std::max(0, maxRefCopy + 1);
      } else if (merge_mode == 4) {
        r_new_copy = r_probe_copy + std::max(0, maxQryCopy + 1);
      } else if (merge_mode == 5) {
        if (rBlk_ptr && qBlk_ptr && rBlk_ptr->getSequences().size() >= qBlk_ptr->getSequences().size()) {
          q_new_copy = q_probe_copy + std::max(0, maxRefCopy + 1);
        } else {
          r_new_copy = r_probe_copy + std::max(0, maxQryCopy + 1);
        }
      }

      // 4. Update CoordinateManager (apply globally consistent coordinate transforms)
      local_timer.start("p2_8_coord_tracker_update");
      size_t rSegCount = 0, qSegCount = 0;
      if (rBlk_ptr) {
        for (const auto& kv : rBlk_ptr->getSequences()) rSegCount += kv.second.getSegments().size();
      }
      if (qBlk_ptr) {
        for (const auto& kv : qBlk_ptr->getSequences()) qSegCount += kv.second.getSegments().size();
      }
      coordMgr.updateAfterMerge(rCore, qCore, mBlk->getId(), finalCigar,
                                actual_merge_inverse, actual_rLen, actual_qLen,
                                merge_mode, maxRefCopy, maxQryCopy,
                                r_probe_copy, q_probe_copy,
                                rSegCount, qSegCount);
      local_timer.stop("p2_8_coord_tracker_update");

      BlockID mId = mBlk->getId();
      dbg.logMergedBlock(mId, mBlk->getConsensus().length());

      // Special logging for Case D (merge_mode == 5): print all Segments and Variants under mBlk
      if (merge_mode == 5) {
        dbg.logCaseDDetails(mBlk.get());
      }

      // Remove rCore and qCore from Set as they have merged into mId
      mergedSet->deleteBlock(rCore);
      mergedSet->deleteBlock(qCore);

      mergeCounter++;

      alnCollection.ref_coverageTracker.syncFromMap(coordMgr, true); // true for Ref
      alnCollection.qry_coverageTracker.syncFromMap(coordMgr, false); // false for Qry
      local_timer.stop("p2_8_coord_tracker_update");

      dbg.validateAfterMerge(mergedSet, this);
    }
  }
  local_timer.stop("phase2_iterative_merge");

  // Free temporary indexed BlockSets to prevent cluttering BlockManager
  this->removeBlockSet(refIndexID);
  this->removeBlockSet(qryIndexID);

  dbg.logProcessedAlignments(mergeCounter);

  // ==========================================
  // Phase 5: Topology reconstruction and strand normalization
  // ==========================================
  local_timer.start("phase5_rewire_edges");
  dbg.logPhase(5, "Re-wiring Linear Pangenome Graph Edges...");
  for (auto &block : mergedSet->getAllBlocks()) {
    auto blk = block.lock();
    if (!blk)
      continue;
    blk->normalizeStrand();
  }
  mergedSet->debugValidateCopies(true);

  for (auto &seq : refSet->getSequences())
    mergedSet->addSequenceName(seq);
  for (auto &seq : qrySet->getSequences())
    mergedSet->addSequenceName(seq);
  // Phase 5: All Ref/Qry Blocks were already placed into mergedSet custody in CoordinateManager::init,
  // so no duplicate re-insertion is required here.

  // 🔍 Diagnostic checkpoint B: after distant block additions
  dbg.validateCheckpoint("[DIAG-B] Validating AFTER distant block additions...", mergedSet, this);

  // Phase 5.5: Mark Distant Blocks and perform self-align (timed via local_timer)
  if (tree) {
    mergedSet->setDistantBlocks(*tree);
  }

  local_timer.start("phase5_5_self_align");
  mergedSet->selfAlignDistant(*option, &coordMgr);
  // mergedSet->debugValidateCopies(true);
  mergedSet->selfAlign(*option, &coordMgr);
  // mergedSet->debugValidateCopies(true);
  local_timer.stop("phase5_5_self_align");

  mergedSet->rebuildLinearGraph(coordMgr);
  local_timer.stop("phase5_rewire_edges");

  // [Phase 6] Recover concrete consensus strings and update Consensus and SNVs via majority voting
  local_timer.start("phase6_recover_consensus");
  dbg.logPhase(6, "Recovering Consensus Strings & Majority Base Voting...");
  for (auto &block : mergedSet->getAllBlocks()) {
    auto blk = block.lock();
    if (!blk) continue;
    blk->refineConsensusAndVariants();
  }
  mergedSet->debugValidateSequences(this);
  local_timer.stop("phase6_recover_consensus");

  dbg.logMergeSummary(time0);

  // Merge accumulated thread-local timer data into global_timer
  global_timer.merge(local_timer);

  return mergedSet;
}