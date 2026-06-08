#include "type.hpp"
#include "block.hpp"


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
// Helper 1: 根據 Map 動態補償 Gap (統一為 Target -> New 視角)
// ==========================================
// 🚨 移除了 isRefTarget！因為丟入此函數的 coordMap 永遠代表 Target，Gap 永遠是 'D'
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
                    consPos++; // Consensus 前進
                }
            } else if (type == 'D') {
                for (int i = 0; i < len; ++i) {
                    if (rPos < rLen) stepCoordMap[rPos++] = consPos;
                    consPos++; // Consensus 前進
                }
            } else if (type == 'I') {
                for (int i = 0; i < len; ++i) {
                    if (qPos < qLen) uCoordMap[qPos++] = consPos;
                    consPos++; // 🚨 核心修復：遇到 'I'，Consensus 也要前進！
                }
            } else if (type == 'S' || type == 'H') {
                for (int i = 0; i < len; ++i) {
                    qPos++; // Soft/hard clipping 消耗 Qry 但不消耗 Consensus
                }
            }
        }
        
        // 收尾：補上陣列的最後一個元素
        if (qPos <= qLen) uCoordMap[qPos] = consPos;
        if (rPos <= rLen) stepCoordMap[rPos] = consPos;
    };

    // ==========================================
    // Phase 3: 核心 Greedy 迴圈 (Linear Block Logic)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 3] Processing Alignments Dynamically...\n";
    int mergeCounter = 0;

    // =========================================================
    // 輔助工具 1：支援 [start, end) 半開區間的完美切割，告別 +- 1
    // =========================================================
    auto splitBlockSafely = [&](BlockID targetMId, int localStart, int localEnd, BlockID& outMiddleId) {
        auto updateMap = [&](std::vector<CoordTracker>& globalMap, BlockID oldID, BlockID leftID, BlockID rightID, int cutPos) {
            for (auto& tracker : globalMap) {
                if (tracker.blkId == oldID) {
                    if (tracker.localPos < cutPos) tracker.blkId = leftID;
                    else { tracker.blkId = rightID; tracker.localPos -= cutPos; }
                }
            }
        };

        BlockID currentId = targetMId;
        auto blk = mergedSet->getBlock(currentId);
        int currentLen = blk->getConsensus().length();

        // 如果範圍剛好涵蓋整個 Block [0, currentLen)，一刀都不用切！
        if (localStart <= 0 && localEnd >= currentLen) {
            outMiddleId = currentId;
            return blk;
        }

        // 第一刀：切掉前面的不要的部分 [0, localStart)
        if (localStart > 0) {
            auto parts = mergedSet->splitSingleBlock(currentId, localStart);
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localStart);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localStart);
            currentId = parts.second; 
            localEnd -= localStart; // 座標平移
        }

        // 第二刀：精準切掉後面的不要的部分
        int lenAfterFirstCut = mergedSet->getBlock(currentId)->getConsensus().length();
        if (localEnd < lenAfterFirstCut) { 
            // 🚨 這裡直接用 localEnd 切！因為 splitSingleBlock 切 localEnd，剛好代表左半邊是 [0, localEnd)
            auto parts = mergedSet->splitSingleBlock(currentId, localEnd); 
            updateMap(refGlobalMap, currentId, parts.first, parts.second, localEnd);
            updateMap(qryGlobalMap, currentId, parts.first, parts.second, localEnd);
            currentId = parts.first; 
        }

        outMiddleId = currentId;
        return mergedSet->getBlock(currentId);
    };

    // =========================================================
    // 輔助工具 2：智慧邊界探測器，轉化為半開區間 [local_s, local_e)
    // =========================================================
    auto getSafeLocalBounds = [&](const std::vector<CoordTracker>& globalMap, int start, int end, BlockID targetMId) {
        int pos1 = globalMap[start].localPos;
        int pos2 = globalMap[end - 1].localPos; // 讀取最後一個有效 index

        int local_s = std::min(pos1, pos2);
        int local_e = std::max(pos1, pos2) + 1; // 🚨 轉化為半開區間

        bool isReversed = (pos1 > pos2);
        bool isStartBoundary = (start == 0 || globalMap[start].blkId != globalMap[start - 1].blkId);
        bool isEndBoundary   = (end >= globalMap.size() || globalMap[end].blkId != targetMId);

        int blockLen = mergedSet->getBlock(targetMId)->getConsensus().length();

        // 如果對齊到邊界，強制定錨到 Block 頭尾
        if (isReversed) {
            if (isStartBoundary) local_e = blockLen;
            if (isEndBoundary)   local_s = 0;
        } else {
            if (isStartBoundary) local_s = 0;
            if (isEndBoundary)   local_e = blockLen;
        }

        return std::make_pair(local_s, local_e);
    };

    // =========================================================
    // 輔助工具 3：根據最新的 GlobalMap，重建 CoverageTracker 的雷達
    // =========================================================
    auto syncTrackersFromMap = [&](std::vector<CoordTracker>& globalMap, CoverageTracker& tracker) {
        // 🌟 1. 直接清空，保證不會有舊的殘留或交界處的髒資料
        tracker.intervals.clear(); 
        
        BlockID currentId = 0;
        int startPos = -1;
        for (int i = 0; i < globalMap.size(); ++i) {
            if (globalMap[i].blkId != currentId) {
                if (currentId != 0 && startPos != -1) {
                    // 🚨 2. 直接使用 i，因為半開區間 [start, i) 剛好完美涵蓋到 i - 1
                    tracker.intervals[startPos] = {i, currentId}; 
                }
                currentId = globalMap[i].blkId;
                startPos = i;
            }
        }
        if (currentId != 0 && startPos != -1) {
            // 🚨 3. 同理，使用 size() 而不是 size() - 1
            tracker.intervals[startPos] = {static_cast<int>(globalMap.size()), currentId};
        }
    };

    auto cigarToStr = [](const CigarString& c) {
        std::string s = "";
        for (auto& op : c) s += std::to_string(op.first) + op.second;
        return s;
    };

   while (true) {
        // 1. 改成接收 vector
        Alignments bestAlns = alnCollection.getBestAlignments(refSet, qrySet, refBounds, qryBounds, L_min);
        if (bestAlns.empty()) break;

        // 2. 把底下整坨 A, B, C, D 的邏輯包進 for 迴圈
        for (Alignment& bestAln : bestAlns) {
            if (!bestAln.valid) break;

            bestAln.CIGAR = compressCigar(bestAln.CIGAR);

            std::cout << "CIGAR: " << cigarToStr(bestAln.CIGAR) << "\n";
            // [start, end)
            int r_start = bestAln.refIdx.first, r_end = bestAln.refIdx.second;
            int q_start = std::min(bestAln.qryIdx.first, bestAln.qryIdx.second); 
            int q_end   = std::max(bestAln.qryIdx.first, bestAln.qryIdx.second);

            

            auto r_overlaps = alnCollection.ref_coverageTracker.getOverlappingIds(r_start, r_end);
            auto q_overlaps = alnCollection.qry_coverageTracker.getOverlappingIds(q_start, q_end);
            
            bool r_merged = !r_overlaps.empty();
            bool q_merged = !q_overlaps.empty();

            mergeCounter++;

            // ---------------------------------------------------------
            // 情境 A：兩端全新
            // ---------------------------------------------------------
            if (!r_merged && !q_merged) {
                if (DEBUG_MODE) std::cout << "  [SCENARIO A] Both ends are new. Extracting blocks...\n";
                bool isCrossing = false;
                bool circularPardoned = false;

                // Check Translocation
                // 1. 收集前後已經存在的 Block ID
                std::set<BlockID> ref_before, ref_after, qry_before, qry_after;
                for(int i = 0; i < r_start; ++i) if (refGlobalMap[i].blkId != 0) ref_before.insert(refGlobalMap[i].blkId);
                for(int i = r_end; i < refGlobalMap.size(); ++i) if (refGlobalMap[i].blkId != 0) ref_after.insert(refGlobalMap[i].blkId);
                for(int i = 0; i < q_start; ++i) if (qryGlobalMap[i].blkId != 0) qry_before.insert(qryGlobalMap[i].blkId);
                for(int i = q_end; i < qryGlobalMap.size(); ++i) if (qryGlobalMap[i].blkId != 0) qry_after.insert(qryGlobalMap[i].blkId);
                // 2. 傳統的易位 (Translocation) 交叉檢測
                for(BlockID id : qry_before) {
                    if (ref_after.count(id)) { isCrossing = true; break; }
                }
                if (!isCrossing) {
                    for(BlockID id : qry_after) {
                        if (ref_before.count(id)) { isCrossing = true; break; }
                    }
                }
                // =======================================================
                // 🌟 3. 環狀基因體特赦 (Circular Wrap-around Exemption)
                // =======================================================
                if (isCrossing) {
                    // 情境 1：Ref 的尾巴 (後面沒積木了) 接上 Qry 的頭 (前面沒積木了)
                    bool isRefTail_QryHead = ref_after.empty() && qry_before.empty();

                    // 情境 2：Ref 的頭 (前面沒積木了) 接上 Qry 的尾巴 (後面沒積木了)
                    bool isRefHead_QryTail = ref_before.empty() && qry_after.empty();
                    if (isRefTail_QryHead || isRefHead_QryTail) {
                        isCrossing = false; // 取消 Crossing 判定，允許 Merge！
                        circularPardoned = true;
                    }
                }
                if (DEBUG_MODE && circularPardoned) {
                    std::cout << "    -> ⭕ [CIRCULAR EXEMPTION] Wrap-around detected! Pardoning the crossing to maintain circular topology.\n";
                }
                auto rBlk = extractBlockFromSuper(mergedSet, refSuperBlock, r_start, r_end);
                auto qBlk = extractBlockFromSuper(mergedSet, qrySuperBlock, q_start, q_end);
                int rLen = rBlk->getConsensus().length(), qLen = qBlk->getConsensus().length();

                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(bestAln.CIGAR, qLen, rLen, uCoordMap, stepCoordMap);

                if (!isCrossing) {
                    // ---------------------------------------------------------
                    // 情境 A1：完美共線性 -> Merge (壓縮成同一個 Block)
                    // ---------------------------------------------------------
                    if (DEBUG_MODE) std::cout << "  [SCENARIO A1] Collinear paths. Extracting and Merging...\n";

                    auto mBlk = mergedSet->mergeTwoBlocks(rBlk, qBlk, bestAln.CIGAR, bestAln.inverse);
                    BlockID rootMId = mBlk->getId();
                    int totalConsLen = mBlk->getConsensus().length();

                    // 先行註冊全域地圖
                    for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rootMId, stepCoordMap[i]};
                    for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {rootMId, uCoordMap[i]};

                } else {
                    // ---------------------------------------------------------
                    // 情境 A2：偵測到跨越 -> Link (保留線性結構 A -> B1 -> C -> B2)
                    // ---------------------------------------------------------
                    if (DEBUG_MODE) std::cout << "  [SCENARIO A2] Crossing/Inversion detected! Linking instead of Merging to prevent loops.\n";

                    mergedSet->linkTwoBlocks(rBlk, qBlk, bestAln.CIGAR, bestAln.inverse);

                    BlockID rBlkId = rBlk->getId();
                    BlockID qBlkId = qBlk->getId();
                    int totalConsLen = rBlk->getConsensus().length();

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

                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                auto qBlk = extractBlockFromSuper(mergedSet, qrySuperBlock, q_start, q_end);
                BlockID qBlkId = qBlk->getId();
                int qLen = qBlk->getConsensus().length();

                if (DEBUG_MODE) std::cout << "    -> Split resulting Core Block ID: " << coreMId << " (Len: " << newSuperLen << ")\n"
                                          << "    -> Extracted Qry Block ID: " << qBlkId << " (Len: " << qLen << ")\n";

                std::vector<int> targetCoordMap; 
                for (int i = r_start; i < r_end; ++i) targetCoordMap.push_back(refGlobalMap[i].localPos);

                bool isTargetReversed = (targetCoordMap.size() > 1 && targetCoordMap.front() > targetCoordMap.back());
                if (isTargetReversed) std::reverse(targetCoordMap.begin(), targetCoordMap.end());

                CigarString linkingCigar = bestAln.CIGAR;
                if (isTargetReversed) std::reverse(linkingCigar.begin(), linkingCigar.end());
                bool finalInverse = bestAln.inverse ^ isTargetReversed;

                CigarString adjustedCigar = adjustCigarWithMap(linkingCigar, targetCoordMap, newSuperLen);

                if (DEBUG_MODE) {
                    std::cout << "    -> Direction: isTargetReversed=" << (isTargetReversed?"YES":"NO") 
                              << ", finalInverse=" << (finalInverse?"YES":"NO") << "\n"
                              << "    -> Orig CIGAR: " << cigarToStr(bestAln.CIGAR) << "\n"
                              << "    -> Adj  CIGAR: " << cigarToStr(adjustedCigar) << "\n";
                }

                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(adjustedCigar, qLen, newSuperLen, uCoordMap, stepCoordMap);

                mergedSet->linkTwoBlocks(targetMBlk, qBlk, adjustedCigar, finalInverse);

                for (auto& tracker : refGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (auto& tracker : qryGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {qBlkId, uCoordMap[i]};

                // syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
                // alnCollection.qry_coverageTracker.overwrite(q_start, q_end, qBlkId);
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

                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                auto rBlk = extractBlockFromSuper(mergedSet, refSuperBlock, r_start, r_end);
                BlockID rBlkId = rBlk->getId();
                int rLen = rBlk->getConsensus().length();

                if (DEBUG_MODE) std::cout << "    -> Split resulting Core Block ID: " << coreMId << " (Len: " << newSuperLen << ")\n"
                                          << "    -> Extracted Ref Block ID: " << rBlkId << " (Len: " << rLen << ")\n";

                std::vector<int> targetCoordMap; 
                for (int i = q_start; i < q_end; ++i) targetCoordMap.push_back(qryGlobalMap[i].localPos);

                bool isTargetReversed = (targetCoordMap.size() > 1 && targetCoordMap.front() > targetCoordMap.back());
                if (isTargetReversed) std::reverse(targetCoordMap.begin(), targetCoordMap.end());

                // 🚨 Scenario C 是 Qry -> Ref，必須手動倒轉 CIGAR 的 I 和 D！
                CigarString invertedCigar = bestAln.CIGAR;
                for (auto& op : invertedCigar) {
                    if (op.second == 'I') op.second = 'D';
                    else if (op.second == 'D') op.second = 'I';
                }
                if (isTargetReversed) std::reverse(invertedCigar.begin(), invertedCigar.end());
                bool finalInverse = bestAln.inverse ^ isTargetReversed;

                CigarString adjustedCigar = adjustCigarWithMap(invertedCigar, targetCoordMap, newSuperLen);

                if (DEBUG_MODE) {
                    std::cout << "    -> Direction: isTargetReversed=" << (isTargetReversed?"YES":"NO") 
                              << ", finalInverse=" << (finalInverse?"YES":"NO") << "\n"
                              << "    -> Orig CIGAR: " << cigarToStr(bestAln.CIGAR) << "\n"
                              << "    -> Inv  CIGAR: " << cigarToStr(invertedCigar) << "\n"
                              << "    -> Adj  CIGAR: " << cigarToStr(adjustedCigar) << "\n";
                }

                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(adjustedCigar, rLen, newSuperLen, uCoordMap, stepCoordMap);

                mergedSet->linkTwoBlocks(targetMBlk, rBlk, adjustedCigar, finalInverse);

                for (auto& tracker : refGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (auto& tracker : qryGlobalMap) {
                    if (tracker.blkId == coreMId && tracker.localPos < stepCoordMap.size()) tracker.localPos = stepCoordMap[tracker.localPos];
                }
                for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rBlkId, uCoordMap[i]};

                // syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
                // alnCollection.ref_coverageTracker.overwrite(r_start, r_end, rBlkId);
            }

            // ---------------------------------------------------------
            // 情境 D：兩端都已存在
            // ---------------------------------------------------------
            else {
                BlockID mIdR = *r_overlaps.begin();
                BlockID mIdQ = *q_overlaps.begin();

                if (DEBUG_MODE) std::cout << "  [SCENARIO D] Both exist. Ref Block " << mIdR << ", Qry Block " << mIdQ << "\n";

                if (mIdR == mIdQ) {
                    if (DEBUG_MODE) std::cout << "    -> [SKIP] Internal repeat within the same Block.\n";
                    continue; 
                }

                // 1. 處理 Ref 端
                auto boundsR = getSafeLocalBounds(refGlobalMap, r_start, r_end, mIdR);
                int localStartR = boundsR.first, localEndR = boundsR.second;
                bool isReversedR = (localStartR > localEndR);
                if (isReversedR) std::swap(localStartR, localEndR);

                BlockID coreMIdR;
                auto targetMBlkR = splitBlockSafely(mIdR, localStartR, localEndR, coreMIdR);
                int lenR = targetMBlkR->getConsensus().length();
                std::vector<int> targetCoordMapR;
                for (int i = r_start; i < r_end; ++i) targetCoordMapR.push_back(refGlobalMap[i].localPos);

                if (DEBUG_MODE) std::cout << "    -> Ref Split bounds: [" << localStartR << ", " << localEndR << ") "
                                          << "-> Core Block ID: " << coreMIdR << " (isReversed: " << (isReversedR?"YES":"NO") << ")\n";

                // 2. 處理 Qry 端
                BlockID currentMIdQ = qryGlobalMap[q_start].blkId;
                if (coreMIdR == currentMIdQ) {
                    if (DEBUG_MODE) std::cout << "    -> [SKIP] Ended up in the same Block after Ref split.\n";
                    continue;
                }

                auto boundsQ = getSafeLocalBounds(qryGlobalMap, q_start, q_end, currentMIdQ);
                int localStartQ = boundsQ.first, localEndQ = boundsQ.second;
                bool isReversedQ = (localStartQ > localEndQ);
                if (isReversedQ) std::swap(localStartQ, localEndQ);

                BlockID coreMIdQ;
                auto targetMBlkQ = splitBlockSafely(currentMIdQ, localStartQ, localEndQ, coreMIdQ);
                int lenQ = targetMBlkQ->getConsensus().length();
                std::vector<int> targetCoordMapQ;
                for (int i = q_start; i < q_end; ++i) targetCoordMapQ.push_back(qryGlobalMap[i].localPos);

                if (DEBUG_MODE) std::cout << "    -> Qry Split bounds: [" << localStartQ << ", " << localEndQ << ") "
                                          << "-> Core Block ID: " << coreMIdQ << " (isReversed: " << (isReversedQ?"YES":"NO") << ")\n";

                // 3. 雙向 CIGAR 映射與反轉
                if (isReversedR) std::reverse(targetCoordMapR.begin(), targetCoordMapR.end());
                if (isReversedQ) std::reverse(targetCoordMapQ.begin(), targetCoordMapQ.end());

                CigarString linkingCigar = bestAln.CIGAR;
                if (isReversedR) std::reverse(linkingCigar.begin(), linkingCigar.end());
                bool finalInverseR = bestAln.inverse ^ isReversedR;

                CigarString cigarR = adjustCigarWithMap(linkingCigar, targetCoordMapR, lenR);
                CigarString cigarR_inv;
                for (auto op : cigarR) {
                    if (op.second == 'I') cigarR_inv.push_back({op.first, 'D'});
                    else if (op.second == 'D') cigarR_inv.push_back({op.first, 'I'});
                    else cigarR_inv.push_back(op);
                }
                if (isReversedQ) std::reverse(cigarR_inv.begin(), cigarR_inv.end());
                bool finalInverseQ = finalInverseR ^ isReversedQ;

                CigarString cigarQ = adjustCigarWithMap(cigarR_inv, targetCoordMapQ, lenQ);
                CigarString finalCigar;
                for (auto op : cigarQ) {
                    if (op.second == 'I') finalCigar.push_back({op.first, 'D'});
                    else if (op.second == 'D') finalCigar.push_back({op.first, 'I'});
                    else finalCigar.push_back(op);
                }

                if (DEBUG_MODE) {
                    std::cout << "    -> Final Link Direction: finalInverseQ=" << (finalInverseQ?"YES":"NO") << "\n"
                              << "    -> Orig CIGAR:  " << cigarToStr(bestAln.CIGAR) << "\n"
                              << "    -> Final CIGAR: " << cigarToStr(finalCigar) << "\n";
                }

                std::vector<int> uCoordMap, stepCoordMap;
                buildCoordMaps(finalCigar, lenQ, lenR, uCoordMap, stepCoordMap);

                mergedSet->linkTwoBlocks(targetMBlkR, targetMBlkQ, finalCigar, finalInverseQ);

                for (auto& tracker : refGlobalMap) {
                    if (tracker.blkId == coreMIdR && tracker.localPos < stepCoordMap.size()) 
                        tracker.localPos = stepCoordMap[tracker.localPos];
                    if (tracker.blkId == coreMIdQ && tracker.localPos < uCoordMap.size()) {
                        tracker.blkId = coreMIdR;
                        tracker.localPos = uCoordMap[tracker.localPos];
                    }
                }
                for (auto& tracker : qryGlobalMap) {
                    if (tracker.blkId == coreMIdR && tracker.localPos < stepCoordMap.size()) 
                        tracker.localPos = stepCoordMap[tracker.localPos];
                    if (tracker.blkId == coreMIdQ && tracker.localPos < uCoordMap.size()) {
                        tracker.blkId = coreMIdR;
                        tracker.localPos = uCoordMap[tracker.localPos];
                    }
                }
                // syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
                // syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
            }
            syncTrackersFromMap(refGlobalMap, alnCollection.ref_coverageTracker);
            syncTrackersFromMap(qryGlobalMap, alnCollection.qry_coverageTracker);
            // mergedSet->debugValidateSegments(true);
        }   
    }
    

    if (DEBUG_MODE) std::cout << "\n  -> Processed " << mergeCounter << " valid alignments.\n";

    // ==========================================
    // Phase 4: 提取未覆蓋的邊角料 (Unused Regions)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 4] Extracting Unmapped Regions to preserve all sequences...\n";
    int unmappedCount = 0;

    auto extractUnusedRegions = [&](const CoverageTracker& tracker, std::shared_ptr<Block> superBlock, int superLen, const std::string& label) {
        int currentPos = 0;
        for (auto const& [start, info] : tracker.intervals) {
            if (currentPos < start) {
                auto unmappedBlk = extractBlockFromSuper(mergedSet, superBlock, currentPos, start);
                if (mergedSet->getBlock(unmappedBlk->getId()) == nullptr) {
                    mergedSet->addBlock(unmappedBlk);
                    unmappedCount++;
                    if (DEBUG_MODE) std::cout << "  -> Extracted " << label << " Gap [" << currentPos << ", " << start << "] as Block ID: " << unmappedBlk->getId() << "\n";
                }
            }
            currentPos = std::max(currentPos, info.end);
        }
        if (currentPos < superLen) {
            auto unmappedBlk = extractBlockFromSuper(mergedSet, superBlock, currentPos, superLen);
            if (mergedSet->getBlock(unmappedBlk->getId()) == nullptr) {
                mergedSet->addBlock(unmappedBlk);
                unmappedCount++;
                if (DEBUG_MODE) std::cout << "  -> Extracted " << label << " Tail [" << currentPos << ", " << superLen << "] as Block ID: " << unmappedBlk->getId() << "\n";
            }
        }
    };

    extractUnusedRegions(alnCollection.ref_coverageTracker, refSuperBlock, refSuperLen, "Ref");
    extractUnusedRegions(alnCollection.qry_coverageTracker, qrySuperBlock, qrySuperLen, "Qry");

    if (DEBUG_MODE) std::cout << "  -> Total " << unmappedCount << " unmapped fragments successfully integrated into Graph.\n";

    // ==========================================
    // Phase 5: 拓撲重建 (Linear Topology Reconstruction)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 5] Re-wiring Linear Pangenome Graph Edges...\n";

    for (auto& block: mergedSet->getAllBlocks()) {
        auto blk = block.lock();
        if (!blk) continue;
        blk->normalizeStrand();
    }
    
    mergedSet->rebuildAllPointers();
    for (auto& seq: refSet->getSequences()) mergedSet->addSequenceName(seq);
    for (auto& seq: qrySet->getSequences()) mergedSet->addSequenceName(seq);

    // Add distant blocks back to the merged blockset
    for (auto& blk: refSet->getAllBlocks()) {
        if (blk.lock()->isDistant()) {
            mergedSet->addBlock(blk.lock());
        }
    }
    for (auto& blk: qrySet->getAllBlocks()) {
        if (blk.lock()->isDistant()) {
            mergedSet->addBlock(blk.lock());
        }
    }

    
    auto timeEnd = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) {
        std::cout << "\n========================================================\n"
                  << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
                  << "Total Execution Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - time0).count() << " ms\n"
                  << "========================================================\n\n";
    }

    return mergedSet;
}
*/


// ==========================================
// 主函數：Greedy Dynamic Graph Merge
// ==========================================
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
    // 🌟 Phase 1.5: 家族同源並查集 (Family Union-Find Ledger)
    // ==========================================
    std::map<int, int> familyAliases;
    int nextFamId = 10000; // 從 10000 開始分配新家族

    auto getTrueFam = [&](int id) {
        if (id == 0) return 0;
        int root = id;
        while (familyAliases.count(root) && familyAliases[root] != root) root = familyAliases[root];
        int curr = id;
        while (curr != root) { int nxt = familyAliases[curr]; familyAliases[curr] = root; curr = nxt; }
        return root;
    };

    auto unifyFam = [&](int id1, int id2) {
        int r1 = getTrueFam(id1), r2 = getTrueFam(id2);
        if (r1 == 0 && r2 == 0) { r1 = ++nextFamId; familyAliases[r1] = r1; return r1; }
        if (r1 == 0) return r2;
        if (r2 == 0) return r1;
        if (r1 != r2) familyAliases[r2] = r1; 
        return r1;
    };

    auto getFamFromBlock = [&](BlockID blkId) {
        auto blk = mergedSet->getBlock(blkId);
        if (!blk) return 0;
        return getTrueFam(blk->getFamilyId());
    };

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

                    int finalFam = 0;
                    if (bestAln.refFamilyId != 0 || bestAln.qryFamilyId != 0) {
                        finalFam = unifyFam(bestAln.refFamilyId, bestAln.qryFamilyId);
                    }

                    auto merge_1 = std::chrono::high_resolution_clock::now();
                    auto mBlk = mergedSet->mergeTwoBlocks(rBlk, qBlk, bestAln.CIGAR, bestAln.inverse);
                    auto merge_2 = std::chrono::high_resolution_clock::now();
                    merge_time += std::chrono::duration_cast<std::chrono::milliseconds>(merge_2 - merge_1).count();

                    mBlk->setFamilyId(getTrueFam(finalFam)); 

                    if (DEBUG_MODE) std::cout << " -> Block " << mBlk->getId() << " (Length: " << mBlk->getConsensus().length() << ")\n";

                    BlockID rootMId = mBlk->getId();
                    for (int i = 0; i < rLen; ++i) refGlobalMap[r_start + i] = {rootMId, stepCoordMap[i]};
                    for (int i = 0; i < qLen; ++i) qryGlobalMap[q_start + i] = {rootMId, uCoordMap[i]};
                } else {
                    // 🌟 A2: 真正的 Crossing -> Link
                    if (DEBUG_MODE) std::cout << "      [SCENARIO A2] Crossing detected! Tying family knot.\n";

                    int finalFam = unifyFam(bestAln.refFamilyId, bestAln.qryFamilyId); 
                    rBlk->setFamilyId(getTrueFam(finalFam));
                    qBlk->setFamilyId(getTrueFam(finalFam));

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
                if (DEBUG_MODE) {
                    std::cout << "    SKIP for now.\n";
                    continue;
                }


                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                // 🌟 發配家族牽線！
                int refFam = getFamFromBlock(coreMId);
                int finalFam = unifyFam(refFam, bestAln.qryFamilyId);
                targetMBlk->setFamilyId(getTrueFam(finalFam));

                // 🌟 extractBlockFromSuper 底層已自動 addBlock
                auto qBlk = extractBlockFromSuper(mergedSet, qrySuperBlock, q_start, q_end);
                BlockID qBlkId = qBlk->getId();
                int qLen = qBlk->getConsensus().length();

                qBlk->setFamilyId(getTrueFam(finalFam)); 

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
                if (DEBUG_MODE) {
                    std::cout << "    SKIP for now.\n";
                    continue;
                }

                BlockID coreMId;
                auto targetMBlk = splitBlockSafely(targetMId, localStart, localEnd, coreMId);
                int newSuperLen = targetMBlk->getConsensus().length();

                // 🌟 發配家族牽線！
                int qryFam = getFamFromBlock(coreMId);
                int finalFam = unifyFam(bestAln.refFamilyId, qryFam);
                targetMBlk->setFamilyId(getTrueFam(finalFam));

                // 🌟 extractBlockFromSuper 底層已自動 addBlock
                auto rBlk = extractBlockFromSuper(mergedSet, refSuperBlock, r_start, r_end);
                BlockID rBlkId = rBlk->getId();
                int rLen = rBlk->getConsensus().length();

                rBlk->setFamilyId(getTrueFam(finalFam)); 

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
                if (DEBUG_MODE) {
                    std::cout << "    SKIP for now.\n";
                    continue;
                }
                if (mIdR == mIdQ) continue; 

                auto boundsR = getSafeLocalBounds(refGlobalMap, r_start, r_end, mIdR);
                int localStartR = boundsR.first, localEndR = boundsR.second;
                bool isReversedR = (localStartR > localEndR);
                if (isReversedR) std::swap(localStartR, localEndR);

                BlockID coreMIdR;
                auto targetMBlkR = splitBlockSafely(mIdR, localStartR, localEndR, coreMIdR);

                BlockID currentMIdQ = qryGlobalMap[q_start].blkId;
                if (coreMIdR == currentMIdQ) continue;

                auto boundsQ = getSafeLocalBounds(qryGlobalMap, q_start, q_end, currentMIdQ);
                int localStartQ = boundsQ.first, localEndQ = boundsQ.second;
                bool isReversedQ = (localStartQ > localEndQ);
                if (isReversedQ) std::swap(localStartQ, localEndQ);

                BlockID coreMIdQ;
                auto targetMBlkQ = splitBlockSafely(currentMIdQ, localStartQ, localEndQ, coreMIdQ);

                // 🌟 發配家族牽線！
                int famR = getFamFromBlock(coreMIdR);
                int famQ = getFamFromBlock(coreMIdQ);
                int finalFam = unifyFam(famR, famQ);
                
                targetMBlkR->setFamilyId(getTrueFam(finalFam));
                targetMBlkQ->setFamilyId(getTrueFam(finalFam));
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
        
        // 🚨 終極結算：確保圖譜中所有的 Family ID 都是收斂到源頭的真 ID！
        if (blk->getFamilyId() != 0) {
            blk->setFamilyId(getTrueFam(blk->getFamilyId()));
        }
    }
    
    
    for (auto& seq: refSet->getSequences()) mergedSet->addSequenceName(seq);
    for (auto& seq: qrySet->getSequences()) mergedSet->addSequenceName(seq);

    for (auto& blk: refSet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }
    for (auto& blk: qrySet->getAllBlocks()) {
        if (blk.lock()->isDistant()) mergedSet->addBlock(blk.lock());
    }

    mergedSet->normalizeFamilyIDs();
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