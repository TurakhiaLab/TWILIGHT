#include "type.hpp"
#include "block_manager.hpp"
#include "global_alignment.hpp"

#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric> // for std::iota
#include <tuple>
#include <chrono>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

// ==========================================
// Helper 1: 根據 Map 動態補償 Gap (含頭尾 Padding 完美版)
// ==========================================
CigarString adjustCigarWithMap(const CigarString& origCigar, const std::vector<int>& coordMap, int superBaseLen) {
    CigarString adjCigar;
    int origRefPos = 0; 

    auto addOp = [&](int len, char op) {
        if (len <= 0) return;
        if (!adjCigar.empty() && adjCigar.back().second == op) adjCigar.back().first += len;
        else adjCigar.push_back({len, op});
    };

    // 1. 補齊前端缺失的 Gaps (Front Padding)
    if (!coordMap.empty() && coordMap[0] > 0) {
        addOp(coordMap[0], 'D');
    }

    // 2. 轉換原本的 CIGAR
    for (const auto& op : origCigar) {
        int len = op.first; char type = op.second;

        if (type == 'S' || type == 'H') continue; 
        
        if (type == 'I') {
            addOp(len, 'I'); 
        } 
        else if (type == 'M' || type == '=' || type == 'X' || type == 'D') {
            for (int i = 0; i < len; ++i) {
                if (origRefPos + 1 < (int)coordMap.size()) {
                    int p1 = coordMap[origRefPos];
                    int p2 = coordMap[origRefPos + 1];

                    if (p2 == p1) {
                        // 如果是 D，代表 New 沒有鹼基，絕對不能新增 'I'
                        if (type != 'D') addOp(1, 'I'); 
                    } else {
                        addOp(1, type);
                        int gaps = p2 - p1 - 1;
                        if (gaps > 0) addOp(gaps, 'D'); 
                    }
                    origRefPos++;
                } else {
                    // 處理 coordMap 的最後一個元素
                    if (origRefPos > 0 && coordMap[origRefPos] == coordMap[origRefPos - 1]) {
                        if (type != 'D') addOp(1, 'I');
                    } else {
                        addOp(1, type);
                    }
                    origRefPos++;
                }
            }
        }
    }

    // ========================================================
    // 🚨 終極修復：補齊尾端缺失的 Gaps (Tail Padding)
    // 當 Block 沒有被切斷，但 Mapping 提早結束時，用 'D' 把剩下的長度填滿
    // ========================================================
    int currentEnd = coordMap.empty() ? -1 : coordMap.back();
    if (superBaseLen > currentEnd + 1) {
        addOp(superBaseLen - (currentEnd + 1), 'D');
    }

    return adjCigar;
}

// ==========================================
// 輔助函數：從 SuperBlock 物理切出指定座標的積木 (基於 Segment::split)
// ==========================================
std::shared_ptr<Block> extractBlockFromSuper(BlockSet* bSet, std::shared_ptr<Block> superBlock, int start, int end) {
    int target_len = end - start;
    std::string subCons = superBlock->getConsensus().substr(start, target_len);
    auto newBlock = bSet->createBlock(subCons); 
    
    for (auto& seqPair : superBlock->getSequences()) {
        Sequence newSeqInfo(seqPair.first); 
        
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment oldSeg = segPair.second; // 拷貝出來處理
            
            // ==========================================
            // 第一刀：切掉前綴 [0, start)
            // ==========================================
            auto split1 = oldSeg.split(start);
            Segment& right_of_start = split1.second; 
            
            // 防呆：如果切完 start 後，右邊完全沒有物理序列 (純粹是 Consensus 上的 Gap)，就提早結束
            if (right_of_start.getStart() == right_of_start.getEnd()) {
                continue;
            }
            
            // ==========================================
            // 第二刀：切掉後綴 [end, total_len)
            // ==========================================
            // 注意：經過第一刀後，right_of_start 的 Variant 相對座標已經被 shift 歸零
            // 所以第二刀的相對切點必須是 target_len (即 end - start)
            auto split2 = right_of_start.split(target_len);
            Segment& target_seg = split2.first;
            
            // ==========================================
            // 驗證與裝載
            // ==========================================
            // 檢查夾在中間的這段目標區塊，是否真實擁有物理序列
            if (target_seg.getStart() != target_seg.getEnd()) {
                
                // 清除殘留的拓撲連線，因為 extract 出來的新積木，其指針會在 Phase 5 重新建立
                target_seg.setPrevBlock(std::shared_ptr<Block>(nullptr));
                target_seg.setNextBlock(std::shared_ptr<Block>(nullptr));
                
                // 此時 target_seg 的 getStart() 已經被 split() 完美映射回真實的 Whole Sequence 座標了
                newSeqInfo.getSegments()[target_seg.getStart()] = std::move(target_seg);
            }
        }
        
        if (newSeqInfo.getSegments().size() > 0) {
            newBlock->addSequence(std::move(newSeqInfo));
        }
    }
    return newBlock;
}

// ==========================================
// 主函數：Greedy Dynamic Graph Merge
// ==========================================
/*
BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, BlockBoundaries& refBounds, BlockBoundaries& qryBounds, AlignmentCollection& alnCollection, int L_min) {
    bool DEBUG_MODE = true;
    auto time0 = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) std::cout << "\n========================================================\n"
                              << "=== TWILIGHT-MGA DYNAMIC MERGE: " << refSet->getId() << " + " << qrySet->getId() << " ===\n"
                              << "========================================================\n";

    // ==========================================
    // Phase 1: 建立 Super Blocks
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 1] Concatenating Super-Blocks...\n";
    
    auto refSuperBlock = refSet->concatenateBlocks(9999991); 
    auto qrySuperBlock = qrySet->concatenateBlocks(9999992);

    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* mergedSet = createBlockSet(newID); 

    int refSuperLen = refSuperBlock->getConsensus().length();
    int qrySuperLen = qrySuperBlock->getConsensus().length();

    if (DEBUG_MODE) std::cout << "  -> Ref SuperBlock Len: " << refSuperLen << ", Qry SuperBlock Len: " << qrySuperLen << "\n";

    // ==========================================
    // Phase 2: 全域座標映射系統
    // ==========================================
    struct CoordTracker { 
        BlockID blkId = 0; 
        int localPos = -1; 
    };
    std::vector<CoordTracker> refGlobalMap(refSuperLen + 1);
    std::vector<CoordTracker> qryGlobalMap(qrySuperLen + 1);

    auto buildCoordMaps = [&](const CigarString& adjustedCigar, int qLen, int rLen, 
                              std::vector<int>& uCoordMap, std::vector<int>& stepCoordMap) {
        uCoordMap.assign(qLen + 1, 0);
        stepCoordMap.assign(rLen + 1, 0);
        
        int rPos = 0, qPos = 0, consPos = 0;
        
        for (auto op : adjustedCigar) {
            int len = op.first; char type = op.second;
            if (type == 'M' || type == '=' || type == 'X') {
                for (int i = 0; i < len; ++i) {
                    if (qPos < qLen) uCoordMap[qPos++] = consPos;
                    if (rPos < rLen) stepCoordMap[rPos++] = consPos;
                    consPos++;
                }
            } else if (type == 'D') {
                for (int i = 0; i < len; ++i) {
                    if (rPos < rLen) stepCoordMap[rPos++] = consPos;
                    consPos++;
                }
            } else if (type == 'I') {
                for (int i = 0; i < len; ++i) {
                    if (qPos < qLen) uCoordMap[qPos++] = consPos;
                    consPos++; 
                }
            } else if (type == 'S' || type == 'H') {
                for (int i = 0; i < len; ++i) qPos++;
            }
        }
        if (qPos <= qLen) uCoordMap[qPos] = consPos;
        if (rPos <= rLen) stepCoordMap[rPos] = consPos;
    };

    // ==========================================
    // Phase 3: 核心 Greedy 迴圈 (Linear Block Logic)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 3] Processing Alignments Dynamically...\n";
    int mergeCounter = 0;

    auto splitBlockSafely = [&](BlockID targetMId, int localStart, int localEnd, BlockID& outMiddleId) {
        auto updateMap = [&](std::vector<CoordTracker>& globalMap, BlockID oldID, BlockID leftID, BlockID rightID, int cutPos) {
            bool found = false;
            for (auto& tracker : globalMap) {
                if (tracker.blkId == oldID) {
                    found = true;
                    if (tracker.localPos < cutPos) tracker.blkId = leftID;
                    else { tracker.blkId = rightID; tracker.localPos -= cutPos; }
                } else if (found) {
                    break; 
                }
            }
        };

        BlockID currentId = targetMId;
        auto blk = mergedSet->getBlock(currentId);
        int currentLen = blk->getConsensus().length();

        if (localStart <= 0 && localEnd >= currentLen) {
            outMiddleId = currentId;
            return blk;
        }

        if (localStart > 0) {
            auto parts = mergedSet->splitSingleBlock(currentId, localStart);
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localStart);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localStart);
            currentId = parts.second; 
            localEnd -= localStart; 
        }

        int lenAfterFirstCut = mergedSet->getBlock(currentId)->getConsensus().length();
        if (localEnd < lenAfterFirstCut) { 
            auto parts = mergedSet->splitSingleBlock(currentId, localEnd); 
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localEnd);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localEnd);
            currentId = parts.first; 
        }

        outMiddleId = currentId;
        return mergedSet->getBlock(currentId);
    };

    auto getSafeLocalBounds = [&](const std::vector<CoordTracker>& globalMap, int start, int end, BlockID targetMId) {
        int pos1 = globalMap[start].localPos;
        int pos2 = globalMap[end - 1].localPos; 

        int local_s = std::min(pos1, pos2);
        int local_e = std::max(pos1, pos2) + 1; 

        bool isReversed = (pos1 > pos2);
        bool isStartBoundary = (start == 0 || globalMap[start].blkId != globalMap[start - 1].blkId);
        bool isEndBoundary   = (end >= globalMap.size() || globalMap[end].blkId != targetMId);

        int blockLen = mergedSet->getBlock(targetMId)->getConsensus().length();

        if (isReversed) {
            if (isStartBoundary) local_e = blockLen;
            if (isEndBoundary)   local_s = 0;
        } else {
            if (isStartBoundary) local_s = 0;
            if (isEndBoundary)   local_e = blockLen;
        }
        return std::make_pair(local_s, local_e);
    };

    auto syncTrackersFromMap = [&](std::vector<CoordTracker>& globalMap, CoverageTracker& tracker) {
        tracker.intervals.clear(); 
        BlockID currentId = 0;
        int startPos = -1;
        for (int i = 0; i < globalMap.size(); ++i) {
            if (globalMap[i].blkId != currentId) {
                if (currentId != 0 && startPos != -1) tracker.intervals[startPos] = {i, currentId}; 
                currentId = globalMap[i].blkId;
                startPos = i;
            }
        }
        if (currentId != 0 && startPos != -1) tracker.intervals[startPos] = {static_cast<int>(globalMap.size()), currentId};
    };

    auto cigarToStr = [](const CigarString& c) {
        std::string s = "";
        for (auto& op : c) s += std::to_string(op.first) + op.second;
        return s;
    };

    uint64_t findBest = 0, merge_time = 0, extract_time = 0;

    // =========================================================
    // 🌟 簡化且獨立的 Family ID 發配系統 (Ref / Qry 雙軌制)
    // =========================================================
    std::map<int, int> ref_family_map; // 專門紀錄 Ref 的舊 Family ID 對應到的新 ID
    std::map<int, int> qry_family_map; // 專門紀錄 Qry 的舊 Family ID 對應到的新 ID
    int next_family_id = 1;

    auto resolveFamily = [&](int rFam, int qFam, bool merge=false) -> int {
        if (rFam == 0 && qFam == 0) {
            if (merge) return 0;
            return next_family_id++; 
        }
        
        int r_mapped = (rFam != 0 && ref_family_map.count(rFam)) ? ref_family_map[rFam] : 0;
        int q_mapped = (qFam != 0 && qry_family_map.count(qFam)) ? qry_family_map[qFam] : 0;
        
        int final_id = 0;

        // 情境 2：兩邊都有綁定過新 ID 了
        if (r_mapped != 0 && q_mapped != 0) {
            if (r_mapped == q_mapped) {
                final_id = r_mapped; // 完美，本來就在同一個新家族
            } else {
                // ⚠️ 衝突發生：它們原本各自發展，現在透過這個 Alignment 連起來了！
                // 解法：把所有歸屬於 q_mapped 的人都過繼給 r_mapped (O(N) 掃描，因為家族數很少，速度極快)
                final_id = r_mapped;
                for (auto& kv : ref_family_map) if (kv.second == q_mapped) kv.second = final_id;
                for (auto& kv : qry_family_map) if (kv.second == q_mapped) kv.second = final_id;
            }
        } 
        // 情境 3：只有一邊有綁定過新 ID，另一邊跟著它
        else if (r_mapped != 0) {
            final_id = r_mapped;
        } 
        else if (q_mapped != 0) {
            final_id = q_mapped;
        } 
        // 情境 4：兩邊都是第一次加入家族系統
        else {
            final_id = next_family_id++;
        }
        
        // 確保這次進來的原始 ID 都有註冊到這個 final_id 上
        if (rFam != 0) ref_family_map[rFam] = final_id;
        if (qFam != 0) qry_family_map[qFam] = final_id;
        
        return final_id;
    };
    

    while (true) {
        auto best_1 = std::chrono::high_resolution_clock::now();
        Alignments bestAlns = alnCollection.getBestAlignments(refSet, qrySet, refSuperBlock, qrySuperBlock, refBounds, qryBounds);
        auto best_2 = std::chrono::high_resolution_clock::now();
        findBest += std::chrono::duration_cast<std::chrono::milliseconds>(best_2 - best_1).count();



        if (bestAlns.empty()) break;

        for (Alignment& bestAln : bestAlns) {
            if (!bestAln.valid) break;

            bestAln.CIGAR = compressCigar(bestAln.CIGAR);

            std::cout << "CIGAR: " << cigarToStr(bestAln.CIGAR) << "\n";
            int r_start = bestAln.refIdx.first, r_end = bestAln.refIdx.second;
            int q_start = std::min(bestAln.qryIdx.first, bestAln.qryIdx.second); 
            int q_end   = std::max(bestAln.qryIdx.first, bestAln.qryIdx.second);

            std::cout << "Ref: (" << r_start << ", " << r_end << "), Qry: (" << q_start << ", " << q_end << ")\n";
            

            auto r_overlaps = alnCollection.ref_coverageTracker.getOverlappingIds(r_start, r_end);
            auto q_overlaps = alnCollection.qry_coverageTracker.getOverlappingIds(q_start, q_end);
            
            bool r_merged = !r_overlaps.empty();
            bool q_merged = !q_overlaps.empty();

            if (r_merged && DEBUG_MODE) {
                std::cout << "Ref Overlap: ";
                for (auto& id : r_overlaps) std::cout << id << ", ";
                std::cout << "\n";
            }
            if (q_merged && DEBUG_MODE) { 
                std::cout << "Qry Overlap: ";
                for (auto& id : q_overlaps) std::cout << id << ", ";
                std::cout << "\n";
            }
            

            mergeCounter++;

            // ---------------------------------------------------------
            // 情境 A：兩端全新
            // ---------------------------------------------------------
            if (!r_merged && !q_merged) {
                if (DEBUG_MODE) std::cout << "  [SCENARIO A] Both ends are new. Extracting blocks...\n";

                auto getIdsBefore = [&](const CoverageTracker& ct, int pos) {
                    std::set<BlockID> res;
                    for (auto const& [s, info] : ct.intervals) if (info.end <= pos) res.insert(info.blockId);
                    return res;
                };
                auto getIdsAfter = [&](const CoverageTracker& ct, int pos) {
                    std::set<BlockID> res;
                    for (auto const& [s, info] : ct.intervals) if (s >= pos) res.insert(info.blockId);
                    return res;
                };

                std::set<BlockID> ref_before = getIdsBefore(alnCollection.ref_coverageTracker, r_start);
                std::set<BlockID> ref_after  = getIdsAfter(alnCollection.ref_coverageTracker, r_end);
                std::set<BlockID> qry_before = getIdsBefore(alnCollection.qry_coverageTracker, q_start);
                std::set<BlockID> qry_after  = getIdsAfter(alnCollection.qry_coverageTracker, q_end);

                // =======================================================
                // 🌟 2. 判定 Crossing (移除所有環狀特赦邏輯)
                // 因為序列已在預處理階段對齊 0 座標，任何拓撲交叉皆視為真實變異！
                // =======================================================
                bool isCrossing = false;
                
                // 只要 Qry 過去的積木出現在 Ref 未來，或 Qry 未來的積木出現在 Ref 過去，就是交叉！
                for(BlockID id : qry_before) {
                    if (ref_after.count(id)) { isCrossing = true; break; }
                }
                if (!isCrossing) {
                    for(BlockID id : qry_after) {
                        if (ref_before.count(id)) { isCrossing = true; break; }
                    }
                }

                // 🌟 extractBlockFromSuper 底層已自動 addBlock
                auto rBlk = extractBlockFromSuper(mergedSet, refSuperBlock, r_start, r_end);
                auto qBlk = extractBlockFromSuper(mergedSet, qrySuperBlock, q_start, q_end);
                int rLen = rBlk->getConsensus().length(), qLen = qBlk->getConsensus().length();

                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(bestAln.CIGAR, qLen, rLen, uCoordMap, stepCoordMap);
                
                // (處理反向對照表的翻轉，如我們先前討論的)
                if (bestAln.inverse) std::reverse(uCoordMap.begin(), uCoordMap.end());

                if (!isCrossing) {
                    // 🌟 A1: 完美共線性 -> Merge
                    if (DEBUG_MODE) std::cout << "      [SCENARIO A1] Collinear paths. Merging...";

                    int finalFam = resolveFamily(bestAln.refFamilyId, bestAln.qryFamilyId, true);

                    auto merge_1 = std::chrono::high_resolution_clock::now();
                    auto mBlk = mergedSet->mergeTwoBlocks(rBlk, qBlk, bestAln.CIGAR, bestAln.inverse);
                    auto merge_2 = std::chrono::high_resolution_clock::now();
                    merge_time += std::chrono::duration_cast<std::chrono::milliseconds>(merge_2 - merge_1).count();

                    mBlk->setFamilyId(finalFam); 

                    if (DEBUG_MODE) std::cout << " -> Block " << mBlk->getId() << " (Length: " << mBlk->getConsensus().length() << ")\n";

                    BlockID rootMId = mBlk->getId();
                    for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rootMId, stepCoordMap[i]};
                    for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {rootMId, uCoordMap[i]};
                } else {
                    // 🌟 A2: 真正的 Crossing -> Link
                    if (DEBUG_MODE) std::cout << "      [SCENARIO A2] Crossing detected! Tying family knot.\n";

                    int finalFam = resolveFamily(bestAln.refFamilyId, bestAln.qryFamilyId); 
                    rBlk->setFamilyId(finalFam);
                    qBlk->setFamilyId(finalFam);

                    BlockID rBlkId = rBlk->getId();
                    BlockID qBlkId = qBlk->getId();
                    for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rBlkId, stepCoordMap[i]};
                    for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {qBlkId, uCoordMap[i]};
                }
            }

            // ---------------------------------------------------------
            // 情境 B：Ref 已存在，Qry 是新的
            // ---------------------------------------------------------
            else if (r_merged && !q_merged) {
                BlockID targetMId = *r_overlaps.begin(); 
                auto bounds = getSafeLocalBounds(refGlobalMap, r_start, r_end, targetMId);
                int localStart = bounds.first, localEnd = bounds.second;

                if (DEBUG_MODE) std::cout << "  [SCENARIO B] Ref exists in Block " << targetMId 
                                          << " | Safe Bounds: [" << localStart << ", " << localEnd << ")\n";
                if (DEBUG_MODE && (localEnd-localStart) < L_min) {
                    std::cout << "    SKIP. (" << (localEnd-localStart) << " < " << L_min << ")\n";
                    continue;
                }

                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                // 🌟 發配家族牽線！
                int refFam = targetMBlk->getFamilyId();
                int finalFam = resolveFamily(refFam, bestAln.qryFamilyId);
                targetMBlk->setFamilyId(finalFam);

                if (DEBUG_MODE) {
                    std::cout << "      ✂️ [SPLIT] Original Block " << targetMId << " split at local [" << localStart << ", " << localEnd << "). Core piece -> Block " << coreMId << "\n";
                    std::cout << "      🧬 [FAMILY-LINK] RefFam(" << refFam << ") + QryFam(" << bestAln.qryFamilyId << ") -> Unified Fam: " << finalFam << "\n";
                }

                // 🌟 extractBlockFromSuper 底層已自動 addBlock
                auto qBlk = extractBlockFromSuper(mergedSet, qrySuperBlock, q_start, q_end);
                BlockID qBlkId = qBlk->getId();
                int qLen = qBlk->getConsensus().length();

                qBlk->setFamilyId(finalFam); 

                std::vector<int> targetCoordMap; 
                for (int i = r_start; i < r_end; ++i) targetCoordMap.push_back(refGlobalMap[i].localPos);

                bool isTargetReversed = (targetCoordMap.size() > 1 && targetCoordMap.front() > targetCoordMap.back());
                if (isTargetReversed) std::reverse(targetCoordMap.begin(), targetCoordMap.end());

                CigarString linkingCigar = bestAln.CIGAR;
                if (isTargetReversed) std::reverse(linkingCigar.begin(), linkingCigar.end());

                CigarString adjustedCigar = adjustCigarWithMap(linkingCigar, targetCoordMap, newSuperLen);
                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(adjustedCigar, qLen, newSuperLen, uCoordMap, stepCoordMap);

                for (auto& tracker : refGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (auto& tracker : qryGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {qBlkId, uCoordMap[i]};
            }

            // ---------------------------------------------------------
            // 情境 C：Qry 已存在，Ref 是新的
            // ---------------------------------------------------------
            else if (!r_merged && q_merged) {
                BlockID targetMId = *q_overlaps.begin(); 
                auto bounds = getSafeLocalBounds(qryGlobalMap, q_start, q_end, targetMId);
                int localStart = bounds.first, localEnd = bounds.second;

                if (DEBUG_MODE) std::cout << "  [SCENARIO C] Qry exists in Block " << targetMId 
                                          << " | Safe Bounds: [" << localStart << ", " << localEnd << ")\n";
                if (DEBUG_MODE && (localEnd-localStart) < L_min) {
                    std::cout << "    SKIP. (" << (localEnd-localStart) << " < " << L_min << ")\n";
                    continue;
                }

                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                // 🌟 發配家族牽線！
                int qryFam = targetMBlk->getFamilyId();
                int finalFam = resolveFamily(bestAln.refFamilyId, qryFam);
                targetMBlk->setFamilyId(finalFam);

                if (DEBUG_MODE) {
                    std::cout << "      ✂️ [SPLIT] Original Block " << targetMId << " split at local [" << localStart << ", " << localEnd << "). Core piece -> Block " << coreMId << "\n";
                    std::cout << "      🧬 [FAMILY-LINK] RefFam(" << bestAln.refFamilyId << ") + QryFam(" << qryFam << ") -> Unified Fam: " << finalFam << "\n";
                }

                // 🌟 extractBlockFromSuper 底層已自動 addBlock
                auto rBlk = extractBlockFromSuper(mergedSet, refSuperBlock, r_start, r_end);
                BlockID rBlkId = rBlk->getId();
                int rLen = rBlk->getConsensus().length();

                rBlk->setFamilyId(finalFam);

                std::vector<int> targetCoordMap; 
                for (int i = q_start; i < q_end; ++i) targetCoordMap.push_back(qryGlobalMap[i].localPos);

                bool isTargetReversed = (targetCoordMap.size() > 1 && targetCoordMap.front() > targetCoordMap.back());
                if (isTargetReversed) std::reverse(targetCoordMap.begin(), targetCoordMap.end());

                CigarString invertedCigar = bestAln.CIGAR;
                for (auto& op : invertedCigar) {
                    if (op.second == 'I') op.second = 'D';
                    else if (op.second == 'D') op.second = 'I';
                }
                if (isTargetReversed) std::reverse(invertedCigar.begin(), invertedCigar.end());

                CigarString adjustedCigar = adjustCigarWithMap(invertedCigar, targetCoordMap, newSuperLen);
                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(adjustedCigar, rLen, newSuperLen, uCoordMap, stepCoordMap);

                for (auto& tracker : refGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (auto& tracker : qryGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rBlkId, uCoordMap[i]};
            }

            // ---------------------------------------------------------
            // 🌟 情境 D：兩端都已存在 (純粹的家族牽線與積木切斷)
            // ---------------------------------------------------------
            else {
                BlockID mIdR = *r_overlaps.begin();
                BlockID mIdQ = *q_overlaps.begin();

                if (DEBUG_MODE) std::cout << "  [SCENARIO D] Both exist. Ref Block " << mIdR << ", Qry Block " << mIdQ << "\n";
                if (mIdR == mIdQ) {
                    if (DEBUG_MODE) std::cout << "      -> SKIP. Ref and Qry map to the exact same block.\n";
                    continue; 
                }

                auto boundsR = getSafeLocalBounds(refGlobalMap, r_start, r_end, mIdR);
                int localStartR = boundsR.first, localEndR = boundsR.second;
                bool isReversedR = (localStartR > localEndR);
                if (isReversedR) std::swap(localStartR, localEndR);

                BlockID coreMIdR;
                auto targetMBlkR = splitBlockSafely(mIdR, localStartR, localEndR, coreMIdR);

                BlockID currentMIdQ = qryGlobalMap[q_start].blkId;
                if (coreMIdR == currentMIdQ) {
                    if (DEBUG_MODE) std::cout << "      -> SKIP. After splitting Ref, Qry resolves to the same core block " << coreMIdR << ".\n";
                    continue;
                }

                auto boundsQ = getSafeLocalBounds(qryGlobalMap, q_start, q_end, currentMIdQ);
                int localStartQ = boundsQ.first, localEndQ = boundsQ.second;
                bool isReversedQ = (localStartQ > localEndQ);
                if (isReversedQ) std::swap(localStartQ, localEndQ);

                BlockID coreMIdQ;
                auto targetMBlkQ = splitBlockSafely(currentMIdQ, localStartQ, localEndQ, coreMIdQ);

                // 🌟 發配家族牽線！
                int famR = targetMBlkR->getFamilyId();
                int famQ = targetMBlkQ->getFamilyId();
                int finalFam = resolveFamily(famR, famQ);
                
                targetMBlkR->setFamilyId(finalFam);
                targetMBlkQ->setFamilyId(finalFam);

                if (DEBUG_MODE) {
                    std::cout << "      ✂️ [SPLIT-R] Original Block " << mIdR << " split at local [" << localStartR << ", " << localEndR << "). Core piece -> Block " << coreMIdR << "\n";
                    std::cout << "      ✂️ [SPLIT-Q] Original Block " << currentMIdQ << " split at local [" << localStartQ << ", " << localEndQ << "). Core piece -> Block " << coreMIdQ << "\n";
                    std::cout << "      🧬 [FAMILY-LINK] RefFam(" << famR << ") + QryFam(" << famQ << ") -> Unified Fam: " << finalFam << "\n";
                }
            }
            syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
            syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
        }   
    }

    if (DEBUG_MODE) std::cout << "\n  -> Processed " << mergeCounter << " valid alignments.\n";

    // ==========================================
    // 🌟 Phase 4: 提取未覆蓋的邊角料 (Unused Regions) - 拓撲邊界感知版
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 4] Extracting Unmapped Regions with Boundary Awareness...\n";
    int unmappedCount = 0;

    // 傳入 tracker, superBlock, 長度, 該軸的邊界 (bounds) 以及 label
    auto extractUnusedRegions = [&](const CoverageTracker& tracker, std::shared_ptr<Block> superBlock, int superLen, const std::map<int, BlockBoundary>& bounds, const std::string& label) {
        
        // 🌟 核心子函數：給定一個 Unmapped 區間，根據內部邊界把它切成多塊再萃取
        auto extractWithBounds = [&](int gapStart, int gapEnd) {
            if (gapStart >= gapEnd) return;

            int chunkStart = gapStart;
            // 找出第一個「嚴格大於」chunkStart 的邊界
            auto it = bounds.upper_bound(chunkStart);

            while (it != bounds.end() && it->first < gapEnd) {
                int bndPos = it->first;
                
                // 抽出從 chunkStart 到 bndPos 的這一段
                if (bndPos > chunkStart) {
                    auto unmappedBlk = extractBlockFromSuper(mergedSet, superBlock, chunkStart, bndPos);
                    unmappedCount++;
                    if (DEBUG_MODE) {
                        std::cout << "  -> ✂️ Extracted " << label << " Chunk [" << chunkStart << ", " << bndPos 
                                  << "] (Len: " << (bndPos - chunkStart) << ") at OLD BOUNDARY as Block ID: " << unmappedBlk->getId() << "\n";
                    }
                    chunkStart = bndPos; // 更新下一個切塊的起點
                }
                ++it;
            }

            // 抽出最後剩下的尾巴 (從最後一個邊界到 gapEnd)
            if (chunkStart < gapEnd) {
                auto unmappedBlk = extractBlockFromSuper(mergedSet, superBlock, chunkStart, gapEnd);
                unmappedCount++;
                if (DEBUG_MODE) {
                    std::cout << "  -> 🧩 Extracted " << label << " Chunk [" << chunkStart << ", " << gapEnd 
                              << "] (Len: " << (gapEnd - chunkStart) << ") as Block ID: " << unmappedBlk->getId() << "\n";
                }
            }
        };

        int currentPos = 0;
        for (auto const& [start, info] : tracker.intervals) {
            if (currentPos < start) {
                // 原本直接呼叫 extractBlockFromSuper，現在改呼叫邊界切分器
                extractWithBounds(currentPos, start);
            }
            // 推進到目前 Coverage 的終點
            currentPos = std::max(currentPos, info.end);
        }
        
        // 處理 SuperBlock 尾部尚未掃描到的區域
        if (currentPos < superLen) {
            extractWithBounds(currentPos, superLen);
        }
    };

    // 呼叫時，分別將 refBounds 和 qryBounds 傳進去
    extractUnusedRegions(alnCollection.ref_coverageTracker, refSuperBlock, refSuperLen, refBounds, "Ref");
    extractUnusedRegions(alnCollection.qry_coverageTracker, qrySuperBlock, qrySuperLen, qryBounds, "Qry");

    if (DEBUG_MODE) std::cout << "  -> Total " << unmappedCount << " unmapped boundary-aware fragments successfully integrated into Graph.\n";
    // ==========================================
    // 🌟 Phase 5: 拓撲重建與家族歸一化
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 5] Re-wiring Linear Pangenome Graph Edges...\n";

    for (auto& block: mergedSet->getAllBlocks()) {
        auto blk = block.lock();
        if (!blk) continue;
        blk->normalizeStrand();
    }
    
    
    for (auto& seq: refSet->getSequences()) mergedSet->addSequenceName(seq);
    for (auto& seq: qrySet->getSequences()) mergedSet->addSequenceName(seq);

    for (auto& blk: refSet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }
    for (auto& blk: qrySet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }

    // mergedSet->normalizeFamilyIDs();
    mergedSet->rebuildAllPointers();

    auto timeEnd = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) {
        std::cout << "\n========================================================\n"
                  << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
                  << "Total Execution Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - time0).count() << " ms\n"
                  << "Find Best Alignment:  " << findBest << " ms\n"
                  << "Merge Blocks:         " << merge_time << " ms\n"
                  << "========================================================\n\n";
    }

    return mergedSet;
}
*/


// 🌟 核心魔法：將原始序列的 CIGAR 投影到「當前積木 Consensus」的座標系上
static CigarString projectBiaxialCigar(const CigarString& orig, 
                                       const std::vector<int>& rMap, const std::vector<int>& qMap, 
                                       int rConsLen, int qConsLen, 
                                       bool inverse, int qOrigLen) 
{
    CigarString newCigar;
    int r_idx = 0, q_idx = 0;
    int curr_r = 0, curr_q = 0;
    
    for (auto op : orig) {
        int len = op.first; char type = op.second;
        for (int i = 0; i < len; ++i) {
            if (type == 'M' || type == '=' || type == 'X') {
                int next_r = rMap[r_idx++];
                
                // 處理反向對齊時的 Qry 座標映射
                int q_idx_fwd = inverse ? (qOrigLen - 1 - q_idx) : q_idx;
                int next_q = inverse ? (qConsLen - 1 - qMap[q_idx_fwd]) : qMap[q_idx_fwd];
                q_idx++;
                
                // 補齊積木 Consensus 中多出來的 Gap (其他路徑造成的 Insertion)
                if (curr_r < next_r) { newCigar.push_back({next_r - curr_r, 'D'}); curr_r = next_r; }
                if (curr_q < next_q) { newCigar.push_back({next_q - curr_q, 'I'}); curr_q = next_q; }
                
                if (!newCigar.empty() && newCigar.back().second == 'M') newCigar.back().first++;
                else newCigar.push_back({1, 'M'});
                
                curr_r++; curr_q++;
            } 
            else if (type == 'D') r_idx++;
            else if (type == 'I') q_idx++;
            else if (type == 'S' || type == 'H') q_idx++;
        }
    }
    // 收尾：補齊積木尾端的長度
    if (curr_r < rConsLen) newCigar.push_back({rConsLen - curr_r, 'D'});
    if (curr_q < qConsLen) newCigar.push_back({qConsLen - curr_q, 'I'});
    
    return compressCigar(newCigar);
}

/*
BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, BlockBoundaries& refBounds, BlockBoundaries& qryBounds, AlignmentCollection& alnCollection, int L_min) {
    bool DEBUG_MODE = true;
    auto time0 = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) std::cout << "\n========================================================\n"
                              << "=== TWILIGHT-MGA DYNAMIC MERGE: " << refSet->getId() << " + " << qrySet->getId() << " ===\n"
                              << "========================================================\n";

    // ==========================================
    // Phase 1: 建立 Super Blocks
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 1] Concatenating Super-Blocks...\n";
    
    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* mergedSet = createBlockSet(newID); 

    auto refOffsets = refSet->getAncestralBlocksOffsets();
    auto qryOffsets = qrySet->getAncestralBlocksOffsets();

    int refSuperLen = (int)(refSet->getRepresentativeConsensus()).begin()->second.size();
    int qrySuperLen = (int)(qrySet->getRepresentativeConsensus()).begin()->second.size();

    if (DEBUG_MODE) std::cout << "  -> Ref SuperBlock Len: " << refSuperLen << ", Qry SuperBlock Len: " << qrySuperLen << "\n";

    // ==========================================
    // Phase 2: 全域座標映射系統
    // ==========================================
    struct CoordTracker { 
        BlockID blkId = 0; 
        int localPos = -1; 
    };
    std::vector<CoordTracker> refGlobalMap(refSuperLen + 1);
    std::vector<CoordTracker> qryGlobalMap(qrySuperLen + 1);

    auto buildCoordMaps = [&](const CigarString& adjustedCigar, int qLen, int rLen, 
                              std::vector<int>& uCoordMap, std::vector<int>& stepCoordMap) {
        uCoordMap.assign(qLen + 1, 0);
        stepCoordMap.assign(rLen + 1, 0);
        
        int rPos = 0, qPos = 0, consPos = 0;
        
        for (auto op : adjustedCigar) {
            int len = op.first; char type = op.second;
            if (type == 'M' || type == '=' || type == 'X') {
                for (int i = 0; i < len; ++i) {
                    if (qPos < qLen) uCoordMap[qPos++] = consPos;
                    if (rPos < rLen) stepCoordMap[rPos++] = consPos;
                    consPos++;
                }
            } else if (type == 'D') {
                for (int i = 0; i < len; ++i) {
                    if (rPos < rLen) stepCoordMap[rPos++] = consPos;
                    consPos++;
                }
            } else if (type == 'I') {
                for (int i = 0; i < len; ++i) {
                    if (qPos < qLen) uCoordMap[qPos++] = consPos;
                    consPos++; 
                }
            } else if (type == 'S' || type == 'H') {
                for (int i = 0; i < len; ++i) qPos++;
            }
        }
        if (qPos <= qLen) uCoordMap[qPos] = consPos;
        if (rPos <= rLen) stepCoordMap[rPos] = consPos;
    };

    // ==========================================
    // Phase 3: 核心 Greedy 迴圈 (Linear Block Logic)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 3] Processing Alignments Dynamically...\n";
    int mergeCounter = 0;

    auto splitBlockSafely = [&](BlockID targetMId, int localStart, int localEnd, BlockID& outMiddleId) {
        
        // 🛠️ 修復 1：移除 break，確保即使有 Repeats (不連續出現) 也能更新到所有同名 Block
        auto updateMap = [&](std::vector<CoordTracker>& globalMap, BlockID oldID, BlockID leftID, BlockID rightID, int cutPos) {
            for (auto& tracker : globalMap) {
                if (tracker.blkId == oldID) {
                    if (tracker.localPos < cutPos) {
                        tracker.blkId = leftID;
                    } else {
                        tracker.blkId = rightID;
                        tracker.localPos -= cutPos;
                    }
                }
                // 🚨 已移除 else if (found) break; 
                // 這樣能強迫掃描完整個陣列，防止漏掉其他位置相同的舊 ID
            }
        };

        BlockID currentId = targetMId;
        auto blk = mergedSet->getBlock(currentId);
        int currentLen = blk->getConsensus().length();

        if (localStart <= 0 && localEnd >= currentLen) {
            outMiddleId = currentId;
            return blk;
        }

        if (localStart > 0) {
            auto parts = mergedSet->splitSingleBlock(currentId, localStart);
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localStart);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localStart);
            currentId = parts.second; 
            localEnd -= localStart; 
        }

        int lenAfterFirstCut = mergedSet->getBlock(currentId)->getConsensus().length();
        if (localEnd < lenAfterFirstCut) { 
            auto parts = mergedSet->splitSingleBlock(currentId, localEnd); 
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localEnd);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localEnd);
            currentId = parts.first; 
        }

        outMiddleId = currentId;
        return mergedSet->getBlock(currentId);
    };

    auto getSafeLocalBounds = [&](const std::vector<CoordTracker>& globalMap, int start, int end, BlockID targetMId) {
        // 🛡️ 防護 1：長度為 0 或負數的防呆
        if (start >= end) {
            if (DEBUG_MODE) std::cout << "      [DEBUG-BOUNDS] ⚠️ Warning: start >= end (" << start << " >= " << end << "). Returning {0,0}\n";
            return std::make_pair(0, 0);
        }

        // 🛡️ 防護 2：陣列全域越界防護
        if (start < 0 || end > globalMap.size()) {
            if (DEBUG_MODE) std::cout << "      [DEBUG-BOUNDS] ⚠️ Error: Out of bounds. Start: " << start << ", End: " << end << ", MapSize: " << globalMap.size() << "\n";
            return std::make_pair(0, 0);
        }

        // 🛡️ 防護 3：空指標安全防護
        auto targetBlk = mergedSet->getBlock(targetMId);
        if (!targetBlk) {
            if (DEBUG_MODE) std::cout << "      [DEBUG-BOUNDS] ⚠️ Error: targetMId " << targetMId << " not found in mergedSet!\n";
            return std::make_pair(0, 0);
        }

        // 🛡️ 防護 4：檢查邊界 ID 是否與 targetMId 一致
        if (DEBUG_MODE) {
            if (globalMap[start].blkId != targetMId || globalMap[end - 1].blkId != targetMId) {
                 std::cout << "      [DEBUG-BOUNDS] ⚠️ Warning: Boundary Block ID mismatch! Expected " << targetMId 
                           << ", but got StartBlk: " << globalMap[start].blkId << ", EndBlk: " << globalMap[end - 1].blkId << "\n";
            }
        }

        // 🌟 取出頭尾的 Local 座標
        int pos1 = globalMap[start].localPos;
        int pos2 = globalMap[end - 1].localPos; 

        // 🛠️ 修復 2：增加陣列資料毀損偵測！如果跨度大於 1，但 pos1 和 pos2 卻一樣，代表陣列沒寫好
        if (DEBUG_MODE && pos1 == pos2 && (end - start) > 1) {
            std::cout << "      [DEBUG-BOUNDS] 🚨 CRITICAL ERROR: pos1 == pos2 (" << pos1 << ") "
                      << "but global range is [" << start << ", " << end << ") (Len: " << end - start << ")!\n"
                      << "      This means your Map's localPos was corrupted or NOT initialized correctly!\n";
        }

        int local_s = std::min(pos1, pos2);
        int local_e = std::max(pos1, pos2) + 1; 

        bool isReversed = (pos1 > pos2);
        bool isStartBoundary = (start == 0 || globalMap[start].blkId != globalMap[start - 1].blkId);
        bool isEndBoundary   = (end >= globalMap.size() || globalMap[end].blkId != targetMId);

        int blockLen = targetBlk->getConsensus().length();
        int orig_s = local_s, orig_e = local_e; 

        if (isReversed) {
            if (isStartBoundary) local_e = blockLen;
            if (isEndBoundary)   local_s = 0;
        } else {
            if (isStartBoundary) local_s = 0;
            if (isEndBoundary)   local_e = blockLen;
        }

        // 🌟 強化 Debug 輸出：加上了 RawLocalPos，讓你能看見原始取出的數值
        if (DEBUG_MODE) {
            std::cout << "      [DEBUG-BOUNDS] Block " << targetMId << " (Len: " << blockLen << ") | "
                      << "Global [" << start << ", " << end << ") -> RawLocalPos {" << pos1 << ", " << pos2 << "} -> "
                      << "Local [" << orig_s << ", " << orig_e << ") | "
                      << "Dir: " << (isReversed ? "Reverse (-)" : "Forward (+)") << " | "
                      << "Snapping: Start=" << (isStartBoundary ? "Yes" : "No") << ", End=" << (isEndBoundary ? "Yes" : "No") << " | "
                      << "Final: [" << local_s << ", " << local_e << ")\n";
        }

        return std::make_pair(local_s, local_e);
    };

    auto syncTrackersFromMap = [&](std::vector<CoordTracker>& globalMap, CoverageTracker& tracker) {
        tracker.intervals.clear(); 
        BlockID currentId = 0;
        int startPos = -1;
        for (int i = 0; i < globalMap.size(); ++i) {
            if (globalMap[i].blkId != currentId) {
                if (currentId != 0 && startPos != -1) tracker.intervals[startPos] = {i, currentId}; 
                currentId = globalMap[i].blkId;
                startPos = i;
            }
        }
        if (currentId != 0 && startPos != -1) tracker.intervals[startPos] = {static_cast<int>(globalMap.size()), currentId};
    };



    uint64_t findBest = 0, merge_time = 0, extract_time = 0;

    // =========================================================
    // 🌟 簡化且獨立的 Family ID 發配系統 (Ref / Qry 雙軌制)
    // =========================================================
    std::map<int, int> ref_family_map; // 專門紀錄 Ref 的舊 Family ID 對應到的新 ID
    std::map<int, int> qry_family_map; // 專門紀錄 Qry 的舊 Family ID 對應到的新 ID
    int next_family_id = 1;

    auto resolveFamily = [&](int rFam, int qFam, bool merge=false) -> int {
        if (rFam == 0 && qFam == 0) {
            if (merge) return 0;
            return next_family_id++; 
        }
        
        int r_mapped = (rFam != 0 && ref_family_map.count(rFam)) ? ref_family_map[rFam] : 0;
        int q_mapped = (qFam != 0 && qry_family_map.count(qFam)) ? qry_family_map[qFam] : 0;
        
        int final_id = 0;

        // 情境 2：兩邊都有綁定過新 ID 了
        if (r_mapped != 0 && q_mapped != 0) {
            if (r_mapped == q_mapped) {
                final_id = r_mapped; // 完美，本來就在同一個新家族
            } else {
                // ⚠️ 衝突發生：它們原本各自發展，現在透過這個 Alignment 連起來了！
                // 解法：把所有歸屬於 q_mapped 的人都過繼給 r_mapped (O(N) 掃描，因為家族數很少，速度極快)
                final_id = r_mapped;
                for (auto& kv : ref_family_map) if (kv.second == q_mapped) kv.second = final_id;
                for (auto& kv : qry_family_map) if (kv.second == q_mapped) kv.second = final_id;
            }
        } 
        // 情境 3：只有一邊有綁定過新 ID，另一邊跟著它
        else if (r_mapped != 0) {
            final_id = r_mapped;
        } 
        else if (q_mapped != 0) {
            final_id = q_mapped;
        } 
        // 情境 4：兩邊都是第一次加入家族系統
        else {
            final_id = next_family_id++;
        }
        
        // 確保這次進來的原始 ID 都有註冊到這個 final_id 上
        if (rFam != 0) ref_family_map[rFam] = final_id;
        if (qFam != 0) qry_family_map[qFam] = final_id;
        
        return final_id;
    };
    

    while (true) {
        auto best_1 = std::chrono::high_resolution_clock::now();
        Alignments bestAlns = alnCollection.getBestAlignments(refSet, qrySet, refBounds, qryBounds);
        auto best_2 = std::chrono::high_resolution_clock::now();
        findBest += std::chrono::duration_cast<std::chrono::milliseconds>(best_2 - best_1).count();

        if (bestAlns.empty()) break;

        for (Alignment& bestAln : bestAlns) {
            if (!bestAln.valid) break;

            bestAln.CIGAR = compressCigar(bestAln.CIGAR);

            // 🌟 1. 偵測純 Indel
            bool is_pure_insertion = (bestAln.CIGAR.size() == 1 && bestAln.CIGAR[0].second == 'I');
            bool is_pure_deletion  = (bestAln.CIGAR.size() == 1 && bestAln.CIGAR[0].second == 'D');

            std::cout << "CIGAR: " << cigarToStr(bestAln.CIGAR) << "\n";
            int r_start = bestAln.refIdx.first, r_end = bestAln.refIdx.second;
            int q_start = std::min(bestAln.qryIdx.first, bestAln.qryIdx.second); 
            int q_end   = std::max(bestAln.qryIdx.first, bestAln.qryIdx.second);

            std::cout << "Ref: (" << r_start << ", " << r_end << "), Qry: (" << q_start << ", " << q_end << ")\n";

            auto r_overlaps = alnCollection.ref_coverageTracker.getOverlappingIds(r_start, r_end);
            auto q_overlaps = alnCollection.qry_coverageTracker.getOverlappingIds(q_start, q_end);
            
            bool r_merged = !r_overlaps.empty();
            bool q_merged = !q_overlaps.empty();

            if (r_merged && DEBUG_MODE) {
                std::cout << "Ref Overlap: ";
                for (auto& id : r_overlaps) std::cout << id << ", ";
                std::cout << "\n";
            }
            if (q_merged && DEBUG_MODE) { 
                std::cout << "Qry Overlap: ";
                for (auto& id : q_overlaps) std::cout << id << ", ";
                std::cout << "\n";
            }

            if (r_merged && q_merged && (*r_overlaps.begin()) == (*q_overlaps.begin()) ) {
                if (DEBUG_MODE) std::cout << "  -> SKIP. Both paths resolve to the exact same merged block " << (*r_overlaps.begin()) << ".\n";
                continue;
            }
            
            int rLen = r_end - r_start;
            int qLen = q_end - q_start;
            
            std::vector<int> rMap(rLen), qMap(qLen);
            std::iota(rMap.begin(), rMap.end(), 0);
            std::iota(qMap.begin(), qMap.end(), 0);

            // =======================================================
            // ⚡ 純 Indel 捷徑
            // =======================================================
            if (is_pure_deletion && rLen > 0) {
                if (!r_merged) {
                    if (DEBUG_MODE) std::cout << "  -> ⚡ [PURE DELETION] Bypassing logic. Extracting Ref directly.\n";
                    auto rBlk = mergedSet->addBlock(refSet->extractBlock(r_start, r_end));
                    for (int i = r_start; i < r_end; ++i) refGlobalMap[i] = {rBlk->getId(), i - r_start};
                    syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
                }
                else {
                    if (DEBUG_MODE) std::cout << "  -> SKIP. Reference already merged.\n";
                }
                continue;
            }

            if (is_pure_insertion && qLen > 0) {
                if (!q_merged) {
                    if (DEBUG_MODE) std::cout << "  -> ⚡ [PURE INSERTION] Bypassing logic. Extracting Qry directly.\n";
                    auto qBlk = mergedSet->addBlock(qrySet->extractBlock(q_start, q_end));
                    for (int i = q_start; i < q_end; ++i) qryGlobalMap[i] = {qBlk->getId(), i - q_start};
                    syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
                }
                else {
                    if (DEBUG_MODE) std::cout << "  -> SKIP. Query already merged.\n";
                }
                continue;
            }

            const int min_threshold = 100;
            if (!bestAln.primary && !is_pure_deletion && !is_pure_insertion) {
                if (scoreCIGAR(bestAln.CIGAR) < min_threshold) {
                    if (DEBUG_MODE) std::cout << "  -> SKIP. Secondary alignment score below threshold ("
                                              << scoreCIGAR(bestAln.CIGAR) << " < " << min_threshold
                                              << ").\n";
                    continue;
                }
            }

            mergeCounter++;
            BlockID rCore = 0, qCore = 0;

            // ==========================================
            // 🌟 提取 Ref Block (修正 Map 提取邏輯)
            // ==========================================
            if (rLen > 0) {
                if (!r_merged) {
                    auto rBlk = mergedSet->addBlock(refSet->extractBlock(r_start, r_end));
                    rCore = rBlk->getId();
                } else {
                    BlockID targetMId = *r_overlaps.begin();

                    std::cout << "R: Target: " << targetMId << '\t';
                    std::cout << "(" << r_start << ", " << r_end << ")\n";
                    

                    auto bounds = getSafeLocalBounds(refGlobalMap, r_start, r_end, targetMId);
                    
                    // 取得切割前的基礎偏移量
                    int baseLocalOffset = bounds.first; 
                    
                    auto blk = splitBlockSafely(targetMId, bounds.first, bounds.second, rCore);
                    rLen = blk->getConsensus().length();

                    rCore = blk->getId();

                    std::cout << " Cutout: " << rCore << '\n';
                    
                    // 🚨 修正：計算 Map 時，必須扣除切割點的 baseLocalOffset，
                    // 這樣 rMap 才會映射到切割後「新 Core Block」內部 0-indexed 的座標！
                    for(int i = 0; i < r_end - r_start; ++i) {
                        rMap[i] = std::max(0, refGlobalMap[r_start + i].localPos - baseLocalOffset);
                    }
                }
            }

            // ==========================================
            // 🌟 提取 Qry Block (修正 Map 提取邏輯)
            // ==========================================
            if (qLen > 0) {
                if (!q_merged) {
                    auto qBlk = mergedSet->addBlock(qrySet->extractBlock(q_start, q_end));
                    qCore = qBlk->getId();
                } else {
                    BlockID targetMId = *q_overlaps.begin(); 
                    auto bounds = getSafeLocalBounds(qryGlobalMap, q_start, q_end, targetMId);

                    std::cout << "Q: Target: " << targetMId << '\t';
                    std::cout << "(" << q_start << ", " << q_end << ")\n";

                    
                    int baseLocalOffset = bounds.first;

                    auto blk = splitBlockSafely(targetMId, bounds.first, bounds.second, qCore);
                    qLen = blk->getConsensus().length();

                    qCore = blk->getId();

                    std::cout << " Cutout: " << qCore << '\n';

                    // 🚨 修正：同上，扣除 baseLocalOffset
                    for(int i = 0; i < q_end - q_start; ++i) {
                        qMap[i] = std::max(0, qryGlobalMap[q_start + i].localPos - baseLocalOffset);
                    }
                }
            }

            auto rBlk_ptr = mergedSet->getBlock(rCore);
            auto qBlk_ptr = mergedSet->getBlock(qCore);

            // ==========================================
            // 🌟 Copy 軌道分發與情境判定
            // ==========================================
            int r_max = rBlk_ptr->getMaxCopy();
            int q_max = qBlk_ptr->getMaxCopy();
            int merge_mode = 0;

            CigarString finalCigar = bestAln.CIGAR;

            if (!r_merged && !q_merged) {
                if (DEBUG_MODE) std::cout << "  [SCENARIO A] Both ends are new. Extracting blocks...\n";

                auto getIdsBefore = [&](const CoverageTracker& ct, int pos) {
                    std::set<BlockID> res;
                    for (auto const& [s, info] : ct.intervals) if (info.end <= pos) res.insert(info.blockId);
                    return res;
                };
                auto getIdsAfter = [&](const CoverageTracker& ct, int pos) {
                    std::set<BlockID> res;
                    for (auto const& [s, info] : ct.intervals) if (s >= pos) res.insert(info.blockId);
                    return res;
                };

                std::set<BlockID> ref_before = getIdsBefore(alnCollection.ref_coverageTracker, r_start);
                std::set<BlockID> ref_after  = getIdsAfter(alnCollection.ref_coverageTracker, r_end);
                std::set<BlockID> qry_before = getIdsBefore(alnCollection.qry_coverageTracker, q_start);
                std::set<BlockID> qry_after  = getIdsAfter(alnCollection.qry_coverageTracker, q_end);

                bool isCrossing = false;
                for(BlockID id : qry_before) if (ref_after.count(id)) { isCrossing = true; break; }
                if (!isCrossing) {
                    for(BlockID id : qry_after) if (ref_before.count(id)) { isCrossing = true; break; }
                }

                if (!isCrossing) {
                    if (DEBUG_MODE) std::cout << "      [SCENARIO A1] Collinear paths. Sharing Copy 0.\n";
                    merge_mode = 1; // 修正: 對應你要求的 mode 1
                } else {
                    if (DEBUG_MODE) std::cout << "      [SCENARIO A2] Crossing detected! Assigning independent copies (0 and 1).\n";
                    merge_mode = 2; // 修正: 對應你要求的 mode 2
                }
            } 
            else if (r_merged && !q_merged) {
                if (DEBUG_MODE) std::cout << "  [SCENARIO B] Ref exists. Extracted core " << rCore << ", Qry new. Qry assigned to Copy " << (r_max + 1) << ".\n";
                merge_mode = 3; // 修正
            } 
            else if (!r_merged && q_merged) {
                if (DEBUG_MODE) std::cout << "  [SCENARIO C] Qry exists. Extracted core " << qCore << ", Ref new. Ref assigned to Copy " << (q_max + 1) << ".\n";
                merge_mode = 4; // 修正
            } 
            else {
                if (DEBUG_MODE) {
                    std::cout << "  [SCENARIO D] Both exist. Extracted cores " << rCore << " and " << qCore << ". Shifting Qry copies by +" << (r_max + 1) << ".\n";
                    std::cout << "      -> Re-aligning Consensus: Ref(" << rBlk_ptr->getConsensus().size() << "bp) vs Qry(" << qBlk_ptr->getConsensus().size() << "bp)...\n";
                }
                merge_mode = 5; 
                
                // =======================================================
                // 🌟 情境 D 專屬：直接對核心積木重新執行 Alignment (Re-align)
                // =======================================================
                std::string rCons = rBlk_ptr->getConsensus();
                std::string qCons = qBlk_ptr->getConsensus();
                
                // 如果是反向對齊，必須把 Qry 的字串做 Reverse Complement！
                if (bestAln.inverse) {
                    qCons = getReverseComplement(qCons); // 假設你有這個工具函數，如果沒有請看下面補充
                }

                // 呼叫你現成的 Semi-Global 函數 (允許頭尾 free gaps，適合 core block 接合)
                std::cout << rCons << '\n' << qCons << '\n';
                finalCigar = runSemiGlobalAlignment(rCons, qCons); 
                
                if (DEBUG_MODE) {
                    std::cout << "      -> Re-aligned CIGAR: "; printCIGAR(finalCigar); std::cout << "\n";
                }
            }

            // ==========================================
            // 🌟 核心合併與 CIGAR 投影 (統一使用 finalCigar)
            // ==========================================
            auto merge_1 = std::chrono::high_resolution_clock::now();
            auto mBlk = mergedSet->mergeTwoBlocks(rBlk_ptr, qBlk_ptr, finalCigar, bestAln.inverse, merge_mode);
            auto merge_2 = std::chrono::high_resolution_clock::now();
            merge_time += std::chrono::duration_cast<std::chrono::milliseconds>(merge_2 - merge_1).count();
            
            BlockID mId = mBlk->getId();
            
            if (DEBUG_MODE) {
                std::cout << "      -> Block " << mId << " created/merged (Length: " << mBlk->getConsensus().length() << ")\n";
            }

            // ==========================================
            // 🌟 結算更新 Global Tracking Maps
            // ==========================================
            std::vector<int> uCoordMap, stepCoordMap;
            buildCoordMaps(finalCigar, qLen, rLen, uCoordMap, stepCoordMap);

            for (int i = r_start; i < r_end; ++i) {
                refGlobalMap[i].blkId = mId;
                refGlobalMap[i].localPos = stepCoordMap[rMap[i - r_start]];
            }
            for (int i = q_start; i < q_end; ++i) {
                qryGlobalMap[i].blkId = mId;
                int pos_fwd = qMap[i - q_start];
                int pos_to_map = bestAln.inverse ? (qLen - 1 - pos_fwd) : pos_fwd;
                qryGlobalMap[i].localPos = uCoordMap[pos_to_map];
            }

            // 🚨 修正：全面掃描 GlobalMap，確保任何角落的 rCore/qCore 都被安全替換為 mId，消滅 Dangling Pointers！
            auto updateGlobalAfterMerge = [&](std::vector<CoordTracker>& globalMap) {
                for (auto& tracker : globalMap) {
                    if (tracker.blkId == rCore) {
                        tracker.blkId = mId;
                        if (tracker.localPos >= 0 && tracker.localPos < stepCoordMap.size()) {
                            tracker.localPos = stepCoordMap[tracker.localPos];
                        }
                    } else if (tracker.blkId == qCore) {
                        tracker.blkId = mId;
                        if (tracker.localPos >= 0 && tracker.localPos < qLen) {
                            int pos_fwd = tracker.localPos;
                            int pos_to_map = bestAln.inverse ? (qLen - 1 - pos_fwd) : pos_fwd;
                            tracker.localPos = uCoordMap[pos_to_map];
                        }
                    }
                }
            };

            updateGlobalAfterMerge(refGlobalMap);
            updateGlobalAfterMerge(qryGlobalMap);
            
            syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
            syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
        } 
    }

    if (DEBUG_MODE) std::cout << "\n  -> Processed " << mergeCounter << " valid alignments.\n";

    // ==========================================
    // 🌟 Phase 4: 提取未覆蓋的邊角料 (Unused Regions) - 拓撲邊界感知版
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 4] Extracting Unmapped Regions with Boundary Awareness...\n";
    int unmappedCount = 0;

    // 🌟 將 superBlock 改為傳入來源的圖譜指標 sourceSet
    auto extractUnusedRegions = [&](const CoverageTracker& tracker, BlockSet* sourceSet, int sourceLen, const std::map<int, BlockBoundary>& bounds, const std::string& label) {
        
        // 🌟 核心子函數：給定一個 Unmapped 區間，根據內部邊界把它切成多塊再萃取
        auto extractWithBounds = [&](int gapStart, int gapEnd) {
            if (gapStart >= gapEnd) return;

            int chunkStart = gapStart;
            // 找出第一個「嚴格大於」chunkStart 的邊界
            auto it = bounds.upper_bound(chunkStart);

            while (it != bounds.end() && it->first < gapEnd) {
                int bndPos = it->first;
                
                // 抽出從 chunkStart 到 bndPos 的這一段
                if (bndPos > chunkStart) {
                    // 🌟 取代 extractBlockFromSuper，直接從 sourceSet 提取並註冊
                    auto unmappedBlk = mergedSet->addBlock(sourceSet->extractBlock(chunkStart, bndPos));

                    unmappedCount++;
                    if (DEBUG_MODE) {
                        std::cout << "  -> ✂️ Extracted " << label << " Chunk [" << chunkStart << ", " << bndPos 
                                  << "] (Len: " << (bndPos - chunkStart) << ") at OLD BOUNDARY as Block ID: " << unmappedBlk->getId() << "\n";
                    }
                    chunkStart = bndPos; // 更新下一個切塊的起點
                }
                ++it;
            }

            // 抽出最後剩下的尾巴 (從最後一個邊界到 gapEnd)
            if (chunkStart < gapEnd) {
                // 🌟 取代 extractBlockFromSuper，直接從 sourceSet 提取並註冊
                auto unmappedBlk = mergedSet->addBlock(sourceSet->extractBlock(chunkStart, gapEnd));
                
                unmappedCount++;
                if (DEBUG_MODE) {
                    std::cout << "  -> 🧩 Extracted " << label << " Chunk [" << chunkStart << ", " << gapEnd 
                              << "] (Len: " << (gapEnd - chunkStart) << ") as Block ID: " << unmappedBlk->getId() << "\n";
                }
            }
        };

        int currentPos = 0;
        for (auto const& [start, info] : tracker.intervals) {
            if (currentPos < start) {
                extractWithBounds(currentPos, start);
            }
            // 推進到目前 Coverage 的終點
            currentPos = std::max(currentPos, info.end);
        }
        
        // 處理尾部尚未掃描到的區域
        if (currentPos < sourceLen) {
            extractWithBounds(currentPos, sourceLen);
        }
    };

    // 🌟 呼叫時，分別將原生的 refSet 和 qrySet 傳進去
    extractUnusedRegions(alnCollection.ref_coverageTracker, refSet, refSuperLen, refBounds, "Ref");
    extractUnusedRegions(alnCollection.qry_coverageTracker, qrySet, qrySuperLen, qryBounds, "Qry");

    if (DEBUG_MODE) std::cout << "  -> Total " << unmappedCount << " unmapped boundary-aware fragments successfully integrated into Graph.\n";
    // ==========================================
    // 🌟 Phase 5: 拓撲重建與家族歸一化
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 5] Re-wiring Linear Pangenome Graph Edges...\n";

    for (auto& block: mergedSet->getAllBlocks()) {
        auto blk = block.lock();
        if (!blk) continue;
        blk->normalizeStrand();
    }
    
    
    for (auto& seq: refSet->getSequences()) mergedSet->addSequenceName(seq);
    for (auto& seq: qrySet->getSequences()) mergedSet->addSequenceName(seq);

    for (auto& blk: refSet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }
    for (auto& blk: qrySet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }

    // mergedSet->normalizeFamilyIDs();
    mergedSet->rebuildAllPointers();

    auto timeEnd = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) {
        std::cout << "\n========================================================\n"
                  << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
                  << "Total Execution Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - time0).count() << " ms\n"
                  << "Find Best Alignment:  " << findBest << " ms\n"
                  << "Merge Blocks:         " << merge_time << " ms\n"
                  << "========================================================\n\n";
    }

    return mergedSet;
}
*/