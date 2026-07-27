#pragma once
#include <vector>
#include <unordered_map>
#include <map>

#include "type.hpp"
#include "block_set.hpp"

enum class BlockOrigin {
    REF,      // 來自原始的 refSet
    QRY,      // 來自原始的 qrySet
    MERGED    // 在本次 merge 過程中誕生的新 Block
};

struct GlobalPosTracker {
    BlockOrigin origin;  
    BlockID blkId = 0;
    int localPos = -1;
    int copyId = 0;      
    bool isMerged = false; 
    bool isInverse = false; // 🌟 新增：記錄這段 Base 當下相對於 Global 軸是不是反向的
    bool inBoth = false;    // 🌟 新增：標註該 Base 是否來自 Mode 1 (同時存在於 Ref 與 Qry)
};

struct BlockOccurrence {
    bool isRef;          
    int globalStart;     
    int globalEnd;       
    bool isInverse = false; // 🌟 新增：記錄這個 Block 區段是不是被反轉過
};
struct BlockInterval {
    int end;
    BlockID blkId;
};

class CoordinateManager {
    private:
        std::vector<GlobalPosTracker> refMap;
        std::vector<GlobalPosTracker> qryMap;

        std::unordered_map<BlockID, std::vector<BlockOccurrence>> reverseMap;

        std::map<int, BlockInterval> refIntervals;
        std::map<int, BlockInterval> qryIntervals;

        // 🌟 記錄包含在 Mode 1 中的 VBlockIDs
        std::set<VBlockID> inBothVBlocks;

        // 🌟 優化：重複使用的緩衝區，避免頻繁 malloc/free
        std::vector<int> buffer_r_step;
        std::vector<int> buffer_q_step;

    public:
        struct GlobalCoords {
            struct IntervalInfo {
                int start;
                int end;
                bool isInverse;
            };
            std::vector<IntervalInfo> refIntervals;
            std::vector<IntervalInfo> qryIntervals;
        };

        // 🌟 快速查詢某個 Global 位置當下的正反股狀態
        bool getIsInverse(int globalPos, bool isRef) const {
            const std::vector<GlobalPosTracker>& targetMap = isRef ? refMap : qryMap;
            if (globalPos >= 0 && globalPos < targetMap.size()) {
                return targetMap[globalPos].isInverse;
            }
            return false; // 防呆預設為正向
        }

        // 🌟 快速查詢某個 Global 位置當下的 Copy ID
        int getCopyId(int globalPos, bool isRef) const {
            const std::vector<GlobalPosTracker>& targetMap = isRef ? refMap : qryMap;
            if (globalPos >= 0 && globalPos < targetMap.size()) {
                return targetMap[globalPos].copyId;
            }
            return 0;
        }

        // 🌟 快速查詢某個 Global 位置當下是否同時存在於 Ref 與 Qry (Mode 1 Merge)
        bool isInBoth(int globalPos, bool isRef) const {
            const std::vector<GlobalPosTracker>& targetMap = isRef ? refMap : qryMap;
            if (globalPos >= 0 && globalPos < targetMap.size()) {
                return targetMap[globalPos].inBoth;
            }
            return false;
        }

        // 🌟 快速查詢某個 VBlock (BlockID + copyId) 是否記錄為同時存在於 Ref 與 Qry
        bool isVBlockInBoth(BlockID blkId, int copyId) const {
            return inBothVBlocks.count({blkId, copyId}) > 0;
        }

        bool isVBlockInBoth(const VBlockID& vid) const {
            return inBothVBlocks.count(vid) > 0;
        }

        const std::set<VBlockID>& getInBothVBlocks() const {
            return inBothVBlocks;
        }

        GlobalCoords getGlobalCoord(BlockID blkId) const {
            GlobalCoords coords;
            auto it = reverseMap.find(blkId);
            if (it != reverseMap.end()) {
                for (const auto& occ : it->second) {
                    if (occ.isRef) {
                        coords.refIntervals.push_back({occ.globalStart, occ.globalEnd, occ.isInverse});
                    } else {
                        coords.qryIntervals.push_back({occ.globalStart, occ.globalEnd, occ.isInverse});
                    }
                }
            }
            return coords;
        }

        // 🌟 取得 Ref 區間樹 (唯讀參照，零拷貝開銷)
        const std::map<int, BlockInterval>& getRefIntervals() const {
            return refIntervals;
        }
    
        // 🌟 取得 Qry 區間樹 (唯讀參照，零拷貝開銷)
        const std::map<int, BlockInterval>& getQryIntervals() const {
            return qryIntervals;
        }
        void init(BlockSet* refSet, BlockSet* qrySet, BlockSet* mergedSet) {
            int refLen = (int)(refSet->getAncestralSequence()).size();
            int qryLen = (int)(qrySet->getAncestralSequence()).size();
                
            refMap.assign(refLen, GlobalPosTracker{});
            qryMap.assign(qryLen, GlobalPosTracker{});
            reverseMap.clear();
                
            refIntervals.clear();
            qryIntervals.clear();
                
            // ==========================================
            // 2. 轉移並註冊 RefSet
            // ==========================================
            int start = 0;
            std::unordered_map<BlockID, std::shared_ptr<Block>> ref_added_blocks; // 🌟 新增：紀錄已加入的 Ref 積木
                
            for (auto& vblkID : refSet->getAncestralBlocks()) {
                auto oldBlk = refSet->getBlock(vblkID.first);
                if (!oldBlk) continue;
            
                std::shared_ptr<Block> newBlk;
                // 🌟 防呆：如果這個積木已經被加過了，直接拿出來用，不要再 addBlock！
                if (ref_added_blocks.count(oldBlk->getId())) {
                    newBlk = ref_added_blocks[oldBlk->getId()];
                } else {
                    newBlk = mergedSet->addBlockReferencing(oldBlk, refSet);
                    ref_added_blocks[oldBlk->getId()] = newBlk; // 紀錄起來
                }
            
                int size = newBlk->getConsensus().size();
                int end = start + size;
            
                reverseMap[newBlk->getId()].push_back({true, start, end});
                refIntervals[start] = {end, newBlk->getId()};
            
                for (int i = start; i < end; ++i) {
                    refMap[i].blkId = newBlk->getId();
                    refMap[i].localPos = i - start;
                    refMap[i].copyId = vblkID.second; 
                    refMap[i].isMerged = false; 
                }
                start = end;
            }
        
            // ==========================================
            // 3. 轉移並註冊 QrySet
            // ==========================================
            start = 0;
            std::unordered_map<BlockID, std::shared_ptr<Block>> qry_added_blocks; // 🌟 新增：紀錄已加入的 Qry 積木
        
            for (auto& vblkID : qrySet->getAncestralBlocks()) {
                auto oldBlk = qrySet->getBlock(vblkID.first);
                if (!oldBlk) continue;
            
                std::shared_ptr<Block> newBlk;
                // 🌟 防呆：如果這個積木已經被加過了，直接拿出來用
                if (qry_added_blocks.count(oldBlk->getId())) {
                    newBlk = qry_added_blocks[oldBlk->getId()];
                } else {
                    newBlk = mergedSet->addBlockReferencing(oldBlk, qrySet);
                    qry_added_blocks[oldBlk->getId()] = newBlk; // 紀錄起來
                }
            
                int size = newBlk->getConsensus().size();
                int end = start + size;
            
                reverseMap[newBlk->getId()].push_back({false, start, end});
                qryIntervals[start] = {end, newBlk->getId()};
            
                // std::cout << oldBlk->getId() << " -> " << newBlk->getId() << ": [" << start << ", " << end << ") / " << qryLen << std::endl;
            
                for (int i = start; i < end; ++i) {
                    qryMap[i].blkId = newBlk->getId();
                    qryMap[i].localPos = i - start;
                    qryMap[i].copyId = vblkID.second; 
                    qryMap[i].isMerged = false; 
                }
                start = end;
            }
        }

        int getLocalPos(int globalPos, bool isRef) const {
            const std::vector<GlobalPosTracker>& targetMap = isRef ? refMap : qryMap;
            if (globalPos >= 0 && globalPos < targetMap.size()) {
                return targetMap[globalPos].localPos;
            }
            return 0; // 防呆
        }

        BlockIDs getOverlapBlocks(int start, int end, bool isRef, bool& alreadyMerged) const {
            const auto& intervalMap = isRef ? refIntervals : qryIntervals;
            const auto& targetMap = isRef ? refMap : qryMap;

            BlockIDs overlapIds;
            
            if (start >= end || intervalMap.empty()) {
                return overlapIds;
            }
        
            auto it = intervalMap.upper_bound(start);
        
            if (it != intervalMap.begin()) {
                auto prevIt = std::prev(it);
                if (prevIt->second.end > start) {
                    overlapIds.push_back(prevIt->second.blkId);
                    // 🌟 O(1) 檢查該區段是否已 Merge
                    if (targetMap[prevIt->first].isMerged) alreadyMerged = true;
                }
            }
        
            while (it != intervalMap.end() && it->first < end) {
                BlockID currentId = it->second.blkId;

                if (overlapIds.empty() || overlapIds.back() != currentId) {
                    if (std::find(overlapIds.begin(), overlapIds.end(), currentId) == overlapIds.end()) {
                        overlapIds.push_back(currentId);
                    }
                }
                // 🌟 O(1) 檢查該區段是否已 Merge
                if (targetMap[it->first].isMerged) alreadyMerged = true;
                it++;
            }
        
            return overlapIds;
        }

        BlockID getAdjacentBlock(int pos, bool isRef, bool lookForward) const {
            const auto& intervalMap = isRef ? refIntervals : qryIntervals;
            if (intervalMap.empty()) return 0;

            if (lookForward) {
                auto it = intervalMap.upper_bound(pos);
                if (it != intervalMap.end()) return it->second.blkId;
            } else {
                auto it = intervalMap.upper_bound(pos);
                if (it != intervalMap.begin()) {
                    --it;
                    // 確保我們找的確實是在 pos 之前的區間
                    if (it->second.end <= pos) return it->second.blkId; 
                    // 如果我們身處區間內，再往前退一個
                    if (it != intervalMap.begin()) {
                        --it;
                        return it->second.blkId;
                    }
                }
            }
            return 0;
        }   

        void updateAfterConcat(bool isRef, const std::vector<BlockID>& old_blocks, BlockID newSuperId) {
            if (old_blocks.empty()) return;

            std::vector<GlobalPosTracker>& targetMap = isRef ? refMap : qryMap;
            std::map<int, BlockInterval>& intervalMap = isRef ? refIntervals : qryIntervals;

            // 1. 取得每一個舊 Block 在 Super Block 裡面的偏移量
            std::unordered_map<BlockID, int> block_super_offsets;
            int current_offset = 0;
            for (BlockID oldId : old_blocks) {
                block_super_offsets[oldId] = current_offset;
                // 直接從 ReverseMap 中隨便找一段，就可以推算它的長度
                if (reverseMap.find(oldId) != reverseMap.end() && !reverseMap[oldId].empty()) {
                    auto& occ = reverseMap[oldId].front();
                    current_offset += (occ.globalEnd - occ.globalStart);
                }
            }
        
            std::vector<BlockOccurrence> new_occurrences;
        
            // 2. 遍歷每一個要被消滅的舊 Block
            for (BlockID oldId : old_blocks) {
                int offset = block_super_offsets[oldId];
                auto it = reverseMap.find(oldId);

                if (it != reverseMap.end()) {
                    for (auto& occ : it->second) {
                        if (occ.isRef != isRef) continue;
                    
                        // A. 更新底層一維 Global 陣列 (指標優化有助於 SIMD 向量化)
                        GlobalPosTracker* base_ptr = &targetMap[occ.globalStart];
                        int length = occ.globalEnd - occ.globalStart;
                        for (int i = 0; i < length; ++i) {
                            base_ptr[i].blkId = newSuperId;
                            base_ptr[i].localPos += offset; // 疊加新的偏移量
                            base_ptr[i].copyId = 0;         // 放棄 Copy 追蹤
                        }
                    
                        // B. 先把舊的區間從樹上拔除
                        intervalMap.erase(occ.globalStart);

                        // 收集成為新的發生紀錄
                        new_occurrences.push_back(occ);
                    }
                    // 該 ID 處理完畢，從反向字典中徹底抹除
                    reverseMap.erase(it);
                }
            }
        
            // 3. 🌟 智慧縫合連續區間 (Interval Stitching)
            if (!new_occurrences.empty()) {
                // 先依照 Global 座標排序
                std::sort(new_occurrences.begin(), new_occurrences.end(), 
                    [](const BlockOccurrence& a, const BlockOccurrence& b) {
                        return a.globalStart < b.globalStart;
                    });
                
                std::vector<BlockOccurrence> merged_occ;
                for (const auto& occ : new_occurrences) {
                    // 如果這段區間跟上一段完美接合
                    if (!merged_occ.empty() && merged_occ.back().globalEnd == occ.globalStart) {
                        // 擴充前一個區間的尾巴
                        merged_occ.back().globalEnd = occ.globalEnd;
                        // 更新區間樹上的尾巴
                        intervalMap[merged_occ.back().globalStart].end = occ.globalEnd;
                    } else {
                        // 如果沒有接合，就當作新的一段寫入
                        merged_occ.push_back(occ);
                        intervalMap[occ.globalStart] = {occ.globalEnd, newSuperId};
                    }
                }

                // 4. 註冊最乾淨、縫合好的 ReverseMap
                reverseMap[newSuperId] = merged_occ;
            }
        }
        
        void updateAfterSplit(BlockID parentID, BlockID leftID, BlockID rightID, int localCut) {
            // 1. 找出這個被切掉的母 Block 在 Global 裡所有的足跡
            auto it = reverseMap.find(parentID);
            if (it == reverseMap.end()) return;
        
            std::vector<BlockOccurrence> left_occs;
            std::vector<BlockOccurrence> right_occs;
            left_occs.reserve(it->second.size());
            right_occs.reserve(it->second.size());
        
            // 2. 遍歷每一個發生過的地方
            for (const auto& occ : it->second) {
                std::vector<GlobalPosTracker>& targetMap = occ.isRef ? refMap : qryMap;
                std::map<int, BlockInterval>& intervalMap = occ.isRef ? refIntervals : qryIntervals;
            
                int gStart = occ.globalStart;
                int gEnd = occ.globalEnd;
            
                if (gStart >= gEnd) continue;

                // 3. 🌟 O(1) 計算全域切點位置 (免除逐 Base 搜尋)
                bool leftFirst = (targetMap[gStart].localPos < localCut);
                int gCut = leftFirst ? (gStart + localCut) : (gEnd - localCut);
                gCut = std::max(gStart, std::min(gEnd, gCut));

                int lStart = leftFirst ? gStart : gCut;
                int lEnd   = leftFirst ? gCut : gEnd;

                int rStart = leftFirst ? gCut : gStart;
                int rEnd   = leftFirst ? gEnd : gCut;

                // 4. 🌟 終極紅黑樹優化：原地覆寫 (In-place) + 帶提示插入 (Hinted Insert, O(1))
                auto tree_it = intervalMap.find(gStart);

                if (lStart < lEnd && rStart < rEnd) {
                    if (leftFirst) {
                        // Left 在前 [gStart, gCut)，Right 在後 [gCut, gEnd)
                        if (tree_it != intervalMap.end()) {
                            tree_it->second = {gCut, leftID}; // 原地覆寫，0 次樹重構！
                        } else {
                            tree_it = intervalMap.insert({gStart, BlockInterval{gCut, leftID}}).first;
                        }
                        // Hinted Insert: 利用 std::next(tree_it) 當作 Hint，O(1) 插入新區段
                        intervalMap.insert(std::next(tree_it), {gCut, BlockInterval{gEnd, rightID}});
                    } else {
                        // Right 在前 [gStart, gCut)，Left 在後 [gCut, gEnd)
                        if (tree_it != intervalMap.end()) {
                            tree_it->second = {gCut, rightID}; // 原地覆寫，0 次樹重構！
                        } else {
                            tree_it = intervalMap.insert({gStart, BlockInterval{gCut, rightID}}).first;
                        }
                        // Hinted Insert: 利用 std::next(tree_it) 當作 Hint，O(1) 插入新區段
                        intervalMap.insert(std::next(tree_it), {gCut, BlockInterval{gEnd, leftID}});
                    }
                } else if (lStart < lEnd) {
                    if (tree_it != intervalMap.end()) {
                        tree_it->second = {gEnd, leftID};
                    } else {
                        intervalMap[gStart] = {gEnd, leftID};
                    }
                } else if (rStart < rEnd) {
                    if (tree_it != intervalMap.end()) {
                        tree_it->second = {gEnd, rightID};
                    } else {
                        intervalMap[gStart] = {gEnd, rightID};
                    }
                } else {
                    if (tree_it != intervalMap.end()) {
                        intervalMap.erase(tree_it);
                    }
                }

                // 5. 無分支批次更新 Left 區段 1D Tracker 與 ReverseMap
                if (lStart < lEnd) {
                    left_occs.push_back({occ.isRef, lStart, lEnd, occ.isInverse});

                    GlobalPosTracker* ptr = &targetMap[lStart];
                    int len = lEnd - lStart;
                    for (int i = 0; i < len; ++i) {
                        ptr[i].blkId = leftID;
                    }
                }

                // 6. 無分支批次更新 Right 區段 1D Tracker 與 ReverseMap
                if (rStart < rEnd) {
                    right_occs.push_back({occ.isRef, rStart, rEnd, occ.isInverse});

                    GlobalPosTracker* ptr = &targetMap[rStart];
                    int len = rEnd - rStart;
                    for (int i = 0; i < len; ++i) {
                        ptr[i].blkId = rightID;
                        ptr[i].localPos -= localCut;
                    }
                }
            }
        
            // 6. 清理並寫入新的 ReverseMap
            reverseMap.erase(it);
            if (!left_occs.empty()) reverseMap[leftID] = std::move(left_occs);
            if (!right_occs.empty()) reverseMap[rightID] = std::move(right_occs);
        }   

        void updateAfterMerge(BlockID rCore, BlockID qCore, BlockID mId,
                          const CigarString& finalCigar, bool qryInverse,
                          int actual_rLen, int actual_qLen,
                          int r_copy, int q_copy,
                          bool inBoth = false,
                          int r_probe_start = -1, int q_probe_start = -1) 
        {
            bool is_in_both = inBoth || (r_copy == q_copy);
            if (is_in_both) {
                inBothVBlocks.insert({mId, r_copy});
                inBothVBlocks.insert({mId, q_copy});
            }

            // ========================================================
            // 1. 內部直接解析 CIGAR，動態建立相對映射 (不變)
            // ========================================================
            buffer_r_step.assign(actual_rLen + 1, 0);
            buffer_q_step.assign(actual_qLen + 1, 0);
            std::vector<int>& r_stepCoordMap = buffer_r_step;
            std::vector<int>& q_uCoordMap = buffer_q_step;
            
            int rPos = 0, qPos = 0, consPos = 0;
            for (const auto& op : finalCigar) {
                int len = op.first; char type = op.second;
                if (type == 'M' || type == '=' || type == 'X') {
                    for (int i = 0; i < len; ++i) {
                        if (qPos < actual_qLen) q_uCoordMap[qPos++] = consPos;
                        if (rPos < actual_rLen) r_stepCoordMap[rPos++] = consPos;
                        consPos++;
                    }
                } else if (type == 'D') { 
                    for (int i = 0; i < len; ++i) {
                        if (rPos < actual_rLen) r_stepCoordMap[rPos++] = consPos;
                        consPos++;
                    }
                } else if (type == 'I') { 
                    for (int i = 0; i < len; ++i) {
                        if (qPos < actual_qLen) q_uCoordMap[qPos++] = consPos;
                        consPos++; 
                    }
                } else if (type == 'S' || type == 'H') {
                    for (int i = 0; i < len; ++i) qPos++;
                }
            }
            if (rPos <= actual_rLen) r_stepCoordMap[rPos] = consPos;
            if (qPos <= actual_qLen) q_uCoordMap[qPos] = consPos;
        
            // ========================================================
            // 2. 核心轉移：利用剛建好的地圖，精準打擊全域區間
            // ========================================================
            auto processCore = [&](BlockID coreId, const std::vector<int>& coordMap, bool isQryCore) {
                auto it = reverseMap.find(coreId);
                if (it == reverseMap.end()) return;

                std::vector<BlockOccurrence> occs = it->second;
                reverseMap.erase(it);

                int mapSize = static_cast<int>(coordMap.size());
                reverseMap[mId].reserve(reverseMap[mId].size() + occs.size());

                for (auto& occ : occs) {
                    std::vector<GlobalPosTracker>& targetMap = occ.isRef ? refMap : qryMap;
                    bool new_isInverse = occ.isInverse;
                    if (isQryCore && qryInverse) new_isInverse = !new_isInverse;

                    int probe_start = isQryCore ? q_probe_start : r_probe_start;
                    bool is_probe_occ = (probe_start < 0) || (occ.globalStart <= probe_start && occ.globalEnd > probe_start);

                    int base_copy = isQryCore ? q_copy : r_copy;
                    int val = is_probe_occ ? base_copy : (targetMap[occ.globalStart].copyId != base_copy ? targetMap[occ.globalStart].copyId : base_copy + 1);
                    bool current_in_both = is_probe_occ && is_in_both;
                
                    for (int i = occ.globalStart; i < occ.globalEnd; ++i) {
                        GlobalPosTracker& tracker = targetMap[i];
                        int old_local = tracker.localPos;
                        int mapped_local = occ.isInverse ? (actual_qLen - 1 - old_local) : old_local;

                        tracker.copyId = val;
                        tracker.blkId = mId;
                        tracker.isMerged = true; 
                        tracker.isInverse = new_isInverse; 
                        tracker.inBoth = current_in_both;
                        tracker.localPos = (mapped_local >= 0 && mapped_local < mapSize) ? coordMap[mapped_local] : 0;
                    }
                    
                    // 3. 更新 IntervalTree
                    std::map<int, BlockInterval>& intervalMap = occ.isRef ? refIntervals : qryIntervals;
                    intervalMap.erase(occ.globalStart);
                    intervalMap[occ.globalStart] = {occ.globalEnd, mId};
                    
                    // 4. 記錄新發生位置
                    BlockOccurrence new_occ = occ;
                    new_occ.isInverse = new_isInverse; 
                    reverseMap[mId].push_back(new_occ);
                }
            };
        
            // 分別處理 Ref Core 與 Qry Core
            processCore(rCore, r_stepCoordMap, false);
            processCore(qCore, q_uCoordMap, true);
        }
        
        // 🌟 尋找最近的「已合併」拓樸錨點 (回傳: {BlockID, CopyNumber})
        std::pair<BlockID, int> getNearestAnchorBlock(int pos, bool isRef, bool lookForward) const {
            const auto& intervalMap = isRef ? refIntervals : qryIntervals;
            const auto& targetMap = isRef ? refMap : qryMap;
        
            // 找不到時，回傳 {-1, -1} 作為空指標
            if (intervalMap.empty()) return {-1, -1};

            if (lookForward) {
                auto it = intervalMap.lower_bound(pos);
                while (it != intervalMap.end()) {
                    if (targetMap[it->first].isMerged && targetMap[it->first].inBoth) {
                        // 🌟 同時回傳 BlockID 與該位置當下的 Copy Number
                        return {it->second.blkId, targetMap[it->first].copyId};
                    }
                    it++; 
                }
            } else {
                auto it = intervalMap.upper_bound(pos);
                while (it != intervalMap.begin()) {
                    --it;
                    if (it->second.end <= pos) {
                        if (targetMap[it->first].isMerged && targetMap[it->first].inBoth) {
                            // 🌟 同時回傳 BlockID 與該位置當下的 Copy Number
                            return {it->second.blkId, targetMap[it->first].copyId};
                        }
                    }
                }
            }
            return {-1, -1}; 
        }
};