
#include "alignment.hpp"
#include "global_alignment.hpp"
#include "type.hpp"
#include "block.hpp"

// =====================================
// CoverageTracker
// =====================================

void CoverageTracker::add(int start, int end, uint64_t blockId) {
    if (start >= end) return;
    intervals[start] = {end, blockId};
}

void CoverageTracker::overwrite(int start, int end, uint64_t newBlockId) {
    if (start >= end) return;
    // 1. 找到第一個「起點大於 start」的區間
    auto it = intervals.upper_bound(start);
    
    // 2. 往前退一步，檢查前一個區間的尾巴有沒有跨越到我們的 start
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > start) {
            it = prev; // 如果有重疊，從前一個開始處理
        }
    }
    // 用來暫存被切斷後，需要保留的左右殘留區段
    std::vector<std::pair<int, SegmentInfo>> leftovers;
    // 3. 走訪所有跟 [start, end) 有重疊的舊區間
    while (it != intervals.end() && it->first < end) {
        int curr_start = it->first;
        int curr_end = it->second.end;
        BlockID curr_id = it->second.blockId;
        // (A) 如果舊區間的左邊凸出去 (Left Overhang)，把凸出去的保留下來
        if (curr_start < start) {
            leftovers.push_back({curr_start, {start, curr_id}});
        }
        
        // (B) 如果舊區間的右邊凸出去 (Right Overhang)，把凸出去的保留下來
        if (curr_end > end) {
            leftovers.push_back({end, {curr_end, curr_id}});
        }
        // (C) 刪除這個舊區間 (因為它要嘛被完全覆蓋，要嘛已經被切成左右兩塊存進 leftovers 了)
        // std::map::erase 會回傳下一個 iterator
        it = intervals.erase(it); 
    }
    // 4. 把保留下來的殘留區段塞回 Map 裡面
    for (const auto& leftover : leftovers) {
        intervals[leftover.first] = leftover.second;
    }
    // 5. 霸氣地塞入我們的新區間！
    intervals[start] = {end, newBlockId};
}

int CoverageTracker::distLeft(int pos) const {
    auto it = intervals.upper_bound(pos);
    if (it == intervals.begin()) return 999999;
    auto prev = std::prev(it);
    if (prev->second.end > pos) return 0;
    return pos - prev->second.end;
}

int CoverageTracker::distRight(int pos) const {
    auto it = intervals.upper_bound(pos);
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > pos) return 0;
    }
    if (it == intervals.end()) return 999999;
    return it->first - pos;
}

int CoverageTracker::getLeftOverlap(int start, int end) const {
    if (start >= end) return 0;
    auto it = intervals.upper_bound(start);
    if (it == intervals.begin()) return 0;
    auto prev = std::prev(it);
    // 如果前一塊積木的尾巴跨過了我們的起點，就算出它覆蓋了我們多少 bp
    if (prev->second.end > start) {
        return std::min(end, prev->second.end) - start;
    }
    return 0;
}

int CoverageTracker::getRightOverlap(int start, int end) const {
    if (start >= end) return 0;
    // 尋找涵蓋 end - 1 (最後一個鹼基) 的積木
    auto it = intervals.upper_bound(end - 1);
    if (it == intervals.begin()) return 0;
    auto prev = std::prev(it);
    // 如果積木涵蓋了我們的尾巴，算出它往回吃掉了我們多少 bp
    if (prev->second.end > end - 1) {
        return end - std::max(start, prev->first);
    }
    return 0;
}

std::set<uint64_t> CoverageTracker::getOverlappingIds(int qStart, int qEnd) const {
    std::set<uint64_t> overlappingIds;
    if (qStart >= qEnd) return overlappingIds;
    
    auto it = intervals.upper_bound(qStart);
    
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > qStart) {
            overlappingIds.insert(prev->second.blockId);
        }
    }
    
    while (it != intervals.end() && it->first < qEnd) {
        overlappingIds.insert(it->second.blockId);
        ++it;
    }
    
    return overlappingIds;
}

void CoverageTracker::getCuts(int start, int end, std::set<int>& cuts) const {
    auto it = intervals.upper_bound(start);
    
    // 1. 檢查前一個區間的尾巴是否落在我們的範圍 (start, end) 內部
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > start && prev->second.end < end) {
            cuts.insert(prev->second.end);
        }
    }
    
    // 2. 處理在範圍內的其他區間
    while (it != intervals.end() && it->first < end) {
        // 【邏輯精簡】：
        // 因為 it 是來自 upper_bound(start)，所以 it->first 絕對大於 start。
        // 因此原本的 `if (it->first > start)` 是必然成立的，可以直接拿掉，直接 insert！
        cuts.insert(it->first);
        
        // 如果這個區間的尾巴也沒有超出我們的範圍，那尾巴也是一個切點
        if (it->second.end < end) {
            cuts.insert(it->second.end);
        }
        ++it;
    }
}

int CoverageTracker::getOverlapLength(int start, int end) const {
    if (start >= end) return 0;
    int overlap = 0;
    auto it = intervals.upper_bound(start);
    
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > start) {
            overlap += std::min(end, prev->second.end) - start;
        }
    }
    
    while (it != intervals.end() && it->first < end) {
        overlap += std::min(end, it->second.end) - it->first;
        ++it;
    }
    return overlap;
}

bool CoverageTracker::isCovered(int start, int end) const {
    return !getOverlappingIds(start, end).empty();
}

// =====================================
// AlignmentCollection
// =====================================
AlignmentCollection::AlignmentCollection(Alignments& alignments, BlockSet* refBlockSet, BlockSet* qryBlockSet) {
    ref_blockset = refBlockSet;
    qry_blockset = qryBlockSet;
    for (auto& aln : alignments) addAlignment(aln);
};

void AlignmentCollection::addAlignment(Alignment aln) {
    aln.updateEnergy(ref_blockset, qry_blockset);
    queue_.push(aln);
}


/*
Alignment AlignmentCollection::getBestAlignment(BlockSet* refBlockSet, BlockSet* qryBlockSet, int L_min) {
    bool debug = false;

    while (!queue_.empty()) {
        Alignment best = queue_.top();
        queue_.pop();

        int q_len = std::abs(best.qryIdx.second - best.qryIdx.first);
        int r_len = std::abs(best.refIdx.second - best.refIdx.first);

        // 基本防呆
        if (best.alnLength < L_min || q_len < L_min || r_len < L_min) continue; 

        int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
        int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
        int check_r_start = best.refIdx.first;
        int check_r_end   = best.refIdx.second;

        // ==========================================
        // 1. 取得絕對的切割點 (只要踩線就拿切刀)
        // ==========================================
        std::set<int> q_cuts, r_cuts;
        qry_coverageTracker.getCuts(check_q_start, check_q_end, q_cuts);
        ref_coverageTracker.getCuts(check_r_start, check_r_end, r_cuts);

        // ==========================================
        // 2. 如果內部有切刀，代表跨越地盤，無情切碎！
        // ==========================================
        if (!q_cuts.empty() || !r_cuts.empty()) {
            Alignments frags = splitSingleAlignment(best, r_cuts, q_cuts);
            for (auto& frag : frags) {
                int f_q_len = std::abs(frag.qryIdx.second - frag.qryIdx.first);
                int f_r_len = std::abs(frag.refIdx.second - frag.refIdx.first);
                
                // 切出來的碎片，只要夠長，全部塞回去重新排隊！
                // 下一次彈出來的時候，它們就會是 100% 乾淨，或是 100% 覆蓋的狀態。
                if (f_q_len >= L_min && f_r_len >= L_min) {
                    frag.updateAlnLength();
                    frag.updateEnergy(refBlockSet, qryBlockSet); 
                    queue_.push(frag);
                }
            }
            continue; 
        }

        // ==========================================
        // 3. 沒有切刀 (Atomic Fragment) -> 狀態判定
        // ==========================================
        int total_q_overlap = qry_coverageTracker.getOverlapLength(check_q_start, check_q_end);
        int total_r_overlap = ref_coverageTracker.getOverlapLength(check_r_start, check_r_end);

        // 【狀態 A】：完全乾淨區段 -> 啟動邊界吸附檢查
        if (total_q_overlap == 0 && total_r_overlap == 0) {
            int qL = qry_coverageTracker.distLeft(check_q_start);
            int qR = qry_coverageTracker.distRight(check_q_end);
            int rL = ref_coverageTracker.distLeft(check_r_start);
            int rR = ref_coverageTracker.distRight(check_r_end);

            bool need_snap = false;
            int q_pad_L = 0, q_pad_R = 0, r_pad_L = 0, r_pad_R = 0;

            // 檢查引力範圍
            if (qL > 0 && qL <= L_min) { q_pad_L = qL; need_snap = true; }
            if (qR > 0 && qR <= L_min) { q_pad_R = qR; need_snap = true; }
            if (rL > 0 && rL <= L_min) { r_pad_L = rL; need_snap = true; }
            if (rR > 0 && rR <= L_min) { r_pad_R = rR; need_snap = true; }

            if (need_snap) {
                if (debug) std::cout << "[DEBUG] 🧲 SNAPPING! Ref: (+" << r_pad_L << ", +" << r_pad_R << "), Qry: (+" << q_pad_L << ", +" << q_pad_R << ")\n";
                snapAlignment(best, r_pad_L, r_pad_R, q_pad_L, q_pad_R);
                best.updateEnergy(refBlockSet, qryBlockSet);
                queue_.push(best); // 吸附完，塞回去重新競爭
                continue;
            }

            // 完全乾淨，且已經對齊邊界 (或離邊界很遠)，完美輸出！
            if (debug) std::cout << "[DEBUG] 🏆 OUTPUT (Case A): " << best.ID << " " << ((best.inverse) ? '-' : '+') << ", Energy: " << best.energy << "\n";
            return best; 
        }

        // 【狀態 B】：100% 已被覆蓋區段 (維持原樣)
        // 因為前面已經切過了，這裡只要有 overlap，就代表它是完美貼齊某個已存在 Block 的內部子片段。
        // 直接 Return 給外面的主迴圈處理！
        if (debug) std::cout << "[DEBUG] 🏆 OUTPUT (Case B/C/D): " << best.ID << " " << ((best.inverse) ? '-' : '+') << " (Merged Fragment)\n";
        return best;
    }

    if (debug) std::cout << "\n[DEBUG] ⚠️ QUEUE EMPTY!\n";
    Alignment empty_aln;
    empty_aln.ID = -1;
    empty_aln.valid = false;
    return empty_aln;
}
*/

Alignments AlignmentCollection::getBestAlignments(
    BlockSet* refBlockSet, 
    BlockSet* qryBlockSet, 
    BlockPtr refSuperBlock, 
    BlockPtr qrySuperBlock, 
    const BlockBoundaries& refBounds, 
    const BlockBoundaries& qryBounds, 
    int L_min) 
{
    bool debug = true; 
    std::vector<Alignment> final_results;
    int snapLength = 30;
        

    // =========================================================
    // ⚔️ 座標投影工具 Lambda 
    // =========================================================
    auto mapCoordinate = [](const Alignment& aln, int targetPos, bool targetIsRef) -> int {
        int rPos = aln.refIdx.first;
        int qStart = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qEnd = std::max(aln.qryIdx.first, aln.qryIdx.second);
        int qCurr = aln.inverse ? qEnd : qStart;

        // 邊界防護檢查
        if (targetIsRef && (targetPos < aln.refIdx.first || targetPos > aln.refIdx.second)) return -1;
        if (!targetIsRef && (targetPos < qStart || targetPos > qEnd)) return -1;

        for (auto op : aln.CIGAR) {
            int len = op.first; char type = op.second;
            int qNext = aln.inverse ? (qCurr - len) : (qCurr + len);

            if (type == 'M' || type == '=' || type == 'X') {
                if (targetIsRef) {
                    if (targetPos >= rPos && targetPos <= rPos + len) 
                        return aln.inverse ? (qCurr - (targetPos - rPos)) : (qCurr + (targetPos - rPos));
                } else {
                    if (targetPos >= std::min(qCurr, qNext) && targetPos <= std::max(qCurr, qNext)) 
                        return rPos + (aln.inverse ? (qCurr - targetPos) : (targetPos - qCurr));
                }
                rPos += len; 
                qCurr = qNext;
            } else if (type == 'D') {
                // Deletion 只有當 Target 是 Ref 時才有可能踩中
                if (targetIsRef && targetPos >= rPos && targetPos <= rPos + len) return qCurr;
                rPos += len;
            } else if (type == 'I' || type == 'S' || type == 'H') {
                // Insertion 只有當 Target 是 Qry 時才有可能踩中 (S/H 忽略)
                if (!targetIsRef && type == 'I' && targetPos >= std::min(qCurr, qNext) && targetPos <= std::max(qCurr, qNext)) 
                    return rPos;
                qCurr = qNext;
            }
        }
        return -1;
    };

    // 🌟 保留這兩個轉發函數，這樣你後面的 code 完全不用改！
    auto mapRefToQry = [&](const Alignment& aln, int targetR) { return mapCoordinate(aln, targetR, true); };
    auto mapQryToRef = [&](const Alignment& aln, int targetQ) { return mapCoordinate(aln, targetQ, false); };

    // =========================================================
    // 🧬 家族 ID 查表工具 (合併縮減版)
    // =========================================================
    auto getFamilyId = [&](int pos, const std::map<int, BlockBoundary>& bounds) -> int {
        if (bounds.empty()) return 0;
        auto it = bounds.upper_bound(pos);
        if (it != bounds.end()) return it->second.leftFamilyId; 
        else return bounds.rbegin()->second.rightFamilyId; 
    };

    auto finalizeAlignment = [&](Alignment& aln) {
        int check_q_start = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int check_q_end   = std::max(aln.qryIdx.first, aln.qryIdx.second);
        int check_r_start = aln.refIdx.first;
        int check_r_end   = aln.refIdx.second;

        int total_q_overlap = qry_coverageTracker.getOverlapLength(check_q_start, check_q_end);
        int total_r_overlap = ref_coverageTracker.getOverlapLength(check_r_start, check_r_end);

        if (total_q_overlap == 0 && total_r_overlap == 0) {
            int qL = qry_coverageTracker.distLeft(check_q_start);
            int qR = qry_coverageTracker.distRight(check_q_end);
            int rL = ref_coverageTracker.distLeft(check_r_start);
            int rR = ref_coverageTracker.distRight(check_r_end);

            bool need_snap = false;
            int q_pad_L = 0, q_pad_R = 0, r_pad_L = 0, r_pad_R = 0;

            if (qL > 0 && qL <= snapLength) { q_pad_L = qL; need_snap = true; }
            if (qR > 0 && qR <= snapLength) { q_pad_R = qR; need_snap = true; }
            if (rL > 0 && rL <= snapLength) { r_pad_L = rL; need_snap = true; }
            if (rR > 0 && rR <= snapLength) { r_pad_R = rR; need_snap = true; }

            if (need_snap) {
                if (debug) std::cout << "    -> [DEBUG] 🧲 SNAPPING Final Fragment! Ref: (+" << r_pad_L << ", +" << r_pad_R << "), Qry: (+" << q_pad_L << ", +" << q_pad_R << ")\n";
                snapAlignment(aln, r_pad_L, r_pad_R, q_pad_L, q_pad_R);
                aln.updateEnergy(refBlockSet, qryBlockSet);
            }
        }
    };


    auto refSeq = refBlockSet->getRepresentativeConsensus()[0].second;
    auto qrySeq = qryBlockSet->getRepresentativeConsensus()[0].second;

    while (!queue_.empty()) {
        Alignment best = queue_.top();
        queue_.pop();

        int q_len = std::abs(best.qryIdx.second - best.qryIdx.first);
        int r_len = std::abs(best.refIdx.second - best.refIdx.first);

        // 🌟 修正 1：WGA 等價模式，容許純 Insertion (r_len=0) 或純 Deletion (q_len=0)
        if (best.alnLength < L_min || std::max(q_len, r_len) < L_min) continue;


        // 確保座標邏輯：start 永遠小於 end
        int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
        int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
        int check_r_start = best.refIdx.first;
        int check_r_end   = best.refIdx.second;

        // ========================================================
        // 🚨 升級版 Debug Message：印出當前 Alignment 的雙軸範圍
        // ========================================================
        if (debug) {
            std::cout << "\n=================================================================\n";
            std::cout << "[DEBUG-BOUNDARY] 🔍 Examining Alignment ID: " << best.ID << "\n"
                      << "  ├─ Ref Bounds: [" << check_r_start << ", " << check_r_end << ") (Len: " << r_len << ")\n"
                      << "  ├─ Qry Bounds: [" << check_q_start << ", " << check_q_end << ") (Len: " << q_len << ")\n"
                      << "  └─ Strand    : " << (best.inverse ? "Reverse (-)" : "Forward (+)") << "\n"
                      << "  └─ CIGAR     : "; printCIGAR(best.CIGAR);  std::cout << '\n'; 
        }

        // ========================================================
        // ⚔️ 🛡️ 搶先防護階段：WGA 雙軸等價微幅修剪器 (Symmetric Micro-Trimmer)
        // ========================================================
        
        // 輔助 Lambda：從 CIGAR 的某一端強制修剪 N 個實體鹼基，並回傳雙軸各被剪了多少
        auto trimCIGARSide = [](CigarString& cigar, int trimLen, bool fromFront) -> std::pair<int, int> {
            int remaining = trimLen;
            int r_shaved = 0, q_shaved = 0;
            
            while (remaining > 0 && !cigar.empty()) {
                auto& op = fromFront ? cigar.front() : cigar.back();
                bool consumesRef = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D');
                bool consumesQry = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I');
                
                int op_len = op.first;
                int consume_step = std::min(remaining, op_len);
                
                if (consumesRef) r_shaved += consume_step;
                if (consumesQry) q_shaved += consume_step;
                
                if (op_len <= remaining) {
                    remaining -= op_len;
                    if (fromFront) cigar.erase(cigar.begin());
                    else           cigar.pop_back();
                } else {
                    op.first -= remaining;
                    remaining = 0;
                }
            }
            return {r_shaved, q_shaved};
        };

        // 🛠️ 1. 處理 CIGAR 起點 (Head) 的微小 Overlap
        int ref_head_ovl = ref_coverageTracker.getLeftOverlap(check_r_start, check_r_end);
        // 如果是反向，CIGAR 起點對應的是 Qry 的絕對右側 (大座標)
        int qry_head_ovl = best.inverse ? qry_coverageTracker.getRightOverlap(check_q_start, check_q_end) 
                                        : qry_coverageTracker.getLeftOverlap(check_q_start, check_q_end);
        
        int head_trim = std::max(ref_head_ovl, qry_head_ovl);

        if (head_trim > 0 && head_trim < L_min) {
            if (debug) std::cout << "  -> ✂️ [PRE-TRIM] Head Overlap (Ref:" << ref_head_ovl << ", Qry:" << qry_head_ovl << "). Shaving " << head_trim << "bp from start...\n";
            auto shaved = trimCIGARSide(best.CIGAR, head_trim, true);
            
            best.refIdx.first += shaved.first;
            
            // 🌟 修正：CIGAR 往前推
            // 正向 Qry 起點 (first) 變大
            // 反向 Qry 終點 (second) 變小
            if (!best.inverse) best.qryIdx.first += shaved.second; 
            else               best.qryIdx.second -= shaved.second; 
        }
        
        // 更新邊界供 Tail 判定使用
        check_r_start = best.refIdx.first;
        check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
        check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);

        // 🛠️ 2. 處理 CIGAR 終點 (Tail) 的微小 Overlap
        int ref_tail_ovl = ref_coverageTracker.getRightOverlap(check_r_start, check_r_end);
        // 如果是反向，CIGAR 終點對應的是 Qry 的絕對左側 (小座標)
        int qry_tail_ovl = best.inverse ? qry_coverageTracker.getLeftOverlap(check_q_start, check_q_end)
                                        : qry_coverageTracker.getRightOverlap(check_q_start, check_q_end);
                                        
        int tail_trim = std::max(ref_tail_ovl, qry_tail_ovl);

        if (tail_trim > 0 && tail_trim < L_min) {
            if (debug) std::cout << "  -> ✂️ [PRE-TRIM] Tail Overlap (Ref:" << ref_tail_ovl << ", Qry:" << qry_tail_ovl << "). Shaving " << tail_trim << "bp from end...\n";
            auto shaved = trimCIGARSide(best.CIGAR, tail_trim, false);
            
            best.refIdx.second -= shaved.first;
            
            // 🌟 修正：CIGAR 尾巴退縮
            // 正向 Qry 終點 (second) 變小
            // 反向 Qry 起點 (first) 變大
            if (!best.inverse) best.qryIdx.second -= shaved.second;
            else               best.qryIdx.first += shaved.second;
        }

        // 最終同步全局 check 邊界
        check_r_end   = best.refIdx.second;
        check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
        check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);

        best.updateAlnLength();
        if (best.alnLength < L_min) {
            if (debug) std::cout << "  -> 🗑️ [PRE-TRIM] Alignment too short after trimming. Skipping.\n";
            continue;
        }

        // ========================================================
        // 🌟 階段 0：動態對齊接回 (Snap / NW) - [覆蓋、拓撲與反向三感知版]
        // ========================================================
        
        // 🛠️ 輔助工具：計算距離最近圖譜邊界的距離
        auto getDistToBound = [&](int pos, const std::map<int, BlockBoundary>& bounds, bool lookLeft) {
            if (bounds.empty()) return -1;
            if (bounds.count(pos)) return 0; 
            if (lookLeft) {
                auto it = bounds.lower_bound(pos);
                if (it != bounds.begin()) return pos - std::prev(it)->first;
            } else {
                auto it = bounds.upper_bound(pos);
                if (it != bounds.end()) return it->first - pos;
            }
            return -1;
        };

        // 🌟 1. 取得 Coverage 距離 (反向序列時，Qry 的左右視野必須顛倒！)
        int qL_cov = best.inverse ? qry_coverageTracker.distRight(check_q_end) : qry_coverageTracker.distLeft(check_q_start);
        int rL_cov = ref_coverageTracker.distLeft(check_r_start);
        int qR_cov = best.inverse ? qry_coverageTracker.distLeft(check_q_start) : qry_coverageTracker.distRight(check_q_end);
        int rR_cov = ref_coverageTracker.distRight(check_r_end);

        // 🌟 2. 取得 Topological Bounds 距離 (同理，反向時左右視野顛倒！)
        int qL_bnd = best.inverse ? getDistToBound(check_q_end, qryBounds, false) : getDistToBound(check_q_start, qryBounds, true);
        int rL_bnd = getDistToBound(check_r_start, refBounds, true);
        int qR_bnd = best.inverse ? getDistToBound(check_q_start, qryBounds, true) : getDistToBound(check_q_end, qryBounds, false);
        int rR_bnd = getDistToBound(check_r_end, refBounds, false);

        auto getAnchorDist = [](int cov_dist, int bnd_dist) {
            if (cov_dist == 0 || bnd_dist == 0) return 0; 
            if (cov_dist > 0 && bnd_dist > 0) return std::min(cov_dist, bnd_dist);
            if (cov_dist > 0) return cov_dist;
            if (bnd_dist > 0) return bnd_dist;
            return -1; 
        };

        int qL = getAnchorDist(qL_cov, qL_bnd);
        int rL = getAnchorDist(rL_cov, rL_bnd);
        int qR = getAnchorDist(qR_cov, qR_bnd);
        int rR = getAnchorDist(rR_cov, rR_bnd);

        const int MAX_NW_LENGTH = 300;

        // --------------------------------------------------------
        // 處理左側 (Left Side / CIGAR Head) 
        // --------------------------------------------------------
        if (qL > 0 && rL > 0 && qL <= MAX_NW_LENGTH && rL <= MAX_NW_LENGTH) {
            if (debug) std::cout << "  -> [STAGE 0] 🔗 Micro-overlap resolved. Running Global Alignment to seal Left Wall (" << rL << ","<< qL << ")...";
            
            int safe_r_start = std::max(0, check_r_start - rL);
            int actual_rL = check_r_start - safe_r_start;
            std::string ref_patch = refSeq.substr(safe_r_start, actual_rL);

            int actual_qL = 0;
            std::string qry_patch = "";
            if (!best.inverse) {
                int safe_q_start = std::max(0, check_q_start - qL);
                actual_qL = check_q_start - safe_q_start;
                qry_patch = qrySeq.substr(safe_q_start, actual_qL);
            } else {
                // 反向序列的 CIGAR Head 往物理大座標延伸
                actual_qL = std::min(qL, (int)qrySeq.length() - check_q_end);
                qry_patch = qrySeq.substr(check_q_end, actual_qL);
            }
            CigarString patchCigar;
            bool use_semi = (std::max(actual_rL, actual_qL) >= 10 * std::min(actual_rL, actual_qL));

            if (debug) {
                std::cout << "  -> [STAGE 0] 🔗 Micro-overlap resolved. Sealing Left Wall (" << rL << ","<< qL << ")..." 
                          << (use_semi ? " [Using Semi-Global]" : " [Using Global]");
            }

            if (use_semi) {
                // 🌟 左側反轉策略：將字串頭尾反轉，讓 Semi-Global 將 Match 逼向右側(原 Block 邊緣)
                std::string rev_ref(ref_patch.rbegin(), ref_patch.rend());
                std::string rev_qry(qry_patch.rbegin(), qry_patch.rend());
                
                patchCigar = runSemiGlobalAlignment(rev_ref, rev_qry);
                // 🌟 把 CIGAR 陣列反轉回來，恢復正常的物理順序
                std::reverse(patchCigar.begin(), patchCigar.end());
            } else {
                patchCigar = runGlobalAlignment(ref_patch, qry_patch);
            }

            if (debug) { std::cout << " -> CIGAR: "; printCIGAR(patchCigar); std::cout << "\n"; }
            
            best.CIGAR.insert(best.CIGAR.begin(), patchCigar.begin(), patchCigar.end());
            best.refIdx.first -= actual_rL;
            check_r_start -= actual_rL;
            
            if (!best.inverse) {
                best.qryIdx.first -= actual_qL;
                check_q_start -= actual_qL;
            } else {
                // 🌟 修正：反向序列左側延伸 -> 增加最大座標 (second)
                best.qryIdx.second += actual_qL; 
                check_q_end += actual_qL;
            }
            
        } else if (rL > 0 && rL <= snapLength && (qL == 0 || qL == -1 || qL >= MAX_NW_LENGTH)) {
            if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Deletion Gap (Left) to edge...\n";
            best.CIGAR.insert(best.CIGAR.begin(), {rL, 'D'});
            best.refIdx.first -= rL; check_r_start -= rL;
        } else if (qL > 0 && qL <= snapLength && (rL == 0 || rL == -1 || rL >= MAX_NW_LENGTH)) {
            if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Insertion Gap (Left) to edge...\n";
            best.CIGAR.insert(best.CIGAR.begin(), {qL, 'I'});
            if (!best.inverse) {
                best.qryIdx.first -= qL;
                check_q_start -= qL;
            } else {
                // 🌟 修正：反向序列左側延伸 -> 增加最大座標 (second)
                best.qryIdx.second += qL; 
                check_q_end += qL;
            }
        }

        // --------------------------------------------------------
        // 處理右側 (Right Side / CIGAR Tail)
        // --------------------------------------------------------
        if (qR > 0 && rR > 0 && qR <= MAX_NW_LENGTH && rR <= MAX_NW_LENGTH) {
            if (debug) std::cout << "  -> [STAGE 0] 🔗 Micro-overlap resolved. Running Global Alignment to seal Right Wall (" << rR << ","<< qR << ")...";
            
            int actual_rR = std::min(rR, (int)refSeq.length() - check_r_end);
            int actual_qR = 0;
            std::string ref_patch = refSeq.substr(check_r_end, actual_rR);
            std::string qry_patch = "";

            if (!best.inverse) {
                actual_qR = std::min(qR, (int)qrySeq.length() - check_q_end);
                qry_patch = qrySeq.substr(check_q_end, actual_qR);
            } else {
                // 反向序列的 CIGAR Tail 往物理小座標延伸
                int safe_q_start = std::max(0, check_q_start - qR);
                actual_qR = check_q_start - safe_q_start;
                qry_patch = qrySeq.substr(safe_q_start, actual_qR);
            }

            CigarString patchCigar;
            // 🌟 判定長度差異
            bool use_semi = (std::max(actual_rR, actual_qR) >= 10 * std::min(actual_rR, actual_qR));

            if (debug) {
                std::cout << "  -> [STAGE 0] 🔗 Micro-overlap resolved. Sealing Right Wall (" << rR << ","<< qR << ")..."
                          << (use_semi ? " [Using Semi-Global]" : " [Using Global]");
            }

            if (use_semi) {
                // 🌟 右側延伸：直接使用修改後的 Semi-Global，它天然就會將 Match 逼向左側 (原 Block 邊緣)
                patchCigar = runSemiGlobalAlignment(ref_patch, qry_patch);
            } else {
                patchCigar = runGlobalAlignment(ref_patch, qry_patch);
            }

            if (debug) { std::cout << " -> CIGAR: "; printCIGAR(patchCigar); std::cout << "\n"; }

            best.CIGAR.insert(best.CIGAR.end(), patchCigar.begin(), patchCigar.end());
            best.refIdx.second += actual_rR;
            check_r_end += actual_rR;
            
            if (!best.inverse) {
                best.qryIdx.second += actual_qR;
                check_q_end += actual_qR;
            } else {
                // 🌟 修正：反向序列右側延伸 -> 減少最小座標 (first)
                best.qryIdx.first -= actual_qR; 
                check_q_start -= actual_qR;
            }
            
        } else if (rR > 0 && rR <= snapLength && (qR == 0 || qR == -1 || qR >= MAX_NW_LENGTH)) {
            if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Deletion Gap (Right) to edge...\n";
            best.CIGAR.push_back({rR, 'D'});
            best.refIdx.second += rR; check_r_end += rR;
        } else if (qR > 0 && qR <= snapLength && (rR == 0 || rR == -1 || rR >= MAX_NW_LENGTH)) {
            if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Insertion Gap (Right) to edge...\n";
            best.CIGAR.push_back({qR, 'I'});
            if (!best.inverse) {
                best.qryIdx.second += qR;
                check_q_end += qR;
            } else {
                // 🌟 修正：反向序列右側延伸 -> 減少最小座標 (first)
                best.qryIdx.first -= qR; 
                check_q_start -= qR;
            }
        }
        
        // ========================================================
        // 🌟 階段 1 & 2：一維化雙軸切點池 (1D Flattened Bi-axial Pool)
        // 核心邏輯：建立自訂排序器，將 2D 座標完美攤平為 1D 單調遞增陣列
        // ========================================================
        
        struct CutPoint {
            bool is_strict = false;
            std::string reasons = "";
        };
        
        // 🌟 修正：反向感知排序器 (Direction-Aware Comparator)
        // 確保在 1D 走訪時，完全貼合 Alignment 的物理前進方向
        auto cmp = [&best](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            if (a.first != b.first) return a.first < b.first;
            // 如果 Ref 相同 (遇到 Insertion)，依照 Strand 方向決定 Qry 的先後
            return best.inverse ? (a.second > b.second) : (a.second < b.second);
        };
        
        std::map<std::pair<int, int>, CutPoint, decltype(cmp)> unified_cuts(cmp);

        auto registerCut = [&](int r, int q, bool strict, const std::string& reason) {
            if (r == -1 || q == -1) return; 
            auto& cut = unified_cuts[{r, q}];
            cut.is_strict |= strict;
            if (!cut.reasons.empty()) cut.reasons += " + ";
            cut.reasons += reason;
        };

        // --------------------------------------------------------
        // 1-1: 收集 Coverage 絕對定錨 (STRICT)
        // --------------------------------------------------------
        std::set<int> tmp_q_cuts, tmp_r_cuts;
        qry_coverageTracker.getCuts(check_q_start, check_q_end, tmp_q_cuts);
        ref_coverageTracker.getCuts(check_r_start, check_r_end, tmp_r_cuts);
        
        if (debug && (!tmp_r_cuts.empty() || !tmp_q_cuts.empty())) {
            std::cout << "  [DEBUG-STAGE 1-1] 🚨 STRICT Coverage Edge Detected!\n";
        }

        for (int r : tmp_r_cuts) {
            int q = mapRefToQry(best, r);
            if (debug) {
                std::cout << "    ⚠️ [STRICT-REF] Hit Coverage Border @ Ref: " << r 
                          << " | Mapped Qry: " << (q != -1 ? std::to_string(q) : "Gap/Out") 
                          << " (Likely due to trailing/heading Overlap)\n";
            }
            registerCut(r, q, true, "Coverage (STRICT Ref)");
        }

        for (int q : tmp_q_cuts) {
            int r = mapQryToRef(best, q);
            if (debug) {
                std::cout << "    ⚠️ [STRICT-QRY] Hit Coverage Border @ Qry: " << q 
                          << " | Mapped Ref: " << (r != -1 ? std::to_string(r) : "Gap/Out") 
                          << " (Likely due to trailing/heading Overlap)\n";
            }
            registerCut(r, q, true, "Coverage (STRICT Qry)");
        }

        // --------------------------------------------------------
        // 1-2: 收集 Long Gaps 邊界 (WGA 等價模式)
        // --------------------------------------------------------
        if (debug) std::cout << "  [DEBUG-STAGE 1&2] Collecting Unified Cut Points...\n";
        
        int rPos = check_r_start;
        int qPos = best.inverse ? check_q_end : check_q_start;

        for (auto op : best.CIGAR) {
            int len = op.first; char type = op.second;

            if (type == 'M' || type == '=' || type == 'X') {
                rPos += len;
                qPos += best.inverse ? -len : len;
            } else if (type == 'D') { 
                if (len >= L_min) { 
                    registerCut(rPos, qPos, false, "Long Ref Deletion Start");
                    registerCut(rPos + len, qPos, false, "Long Ref Deletion End");
                }
                rPos += len;
            } else if (type == 'I') { 
                if (len >= L_min) {
                    int next_qPos = qPos + (best.inverse ? -len : len);
                    registerCut(rPos, qPos, false, "Long Qry Insertion Start");
                    registerCut(rPos, next_qPos, false, "Long Qry Insertion End");
                }
                qPos += best.inverse ? -len : len;
            } else if (type == 'S' || type == 'H') {
                qPos += best.inverse ? -len : len;
            }
        }

        // --------------------------------------------------------
        // 1-3: 收集圖譜舊邊界
        // --------------------------------------------------------
        if (debug) std::cout << "  [DEBUG-STAGE 1-3] Collecting Old Graph Boundaries...\n";

        for (auto it = refBounds.lower_bound(check_r_start); it != refBounds.end() && it->first <= check_r_end; ++it) {
            // 🌟 修正：不忽略 STRICT！根據舊圖譜的 BoundaryType 決定是否為絕對定錨
            bool is_strict = (it->second.type == BoundaryType::STRICT);
            int q = mapRefToQry(best, it->first);
            
            if (debug && is_strict) {
                std::cout << "    ⚠️ [STRICT-OLD-REF] Recovered Topological Boundary @ Ref: " << it->first 
                          << " (Reason: " << it->second.reason << ")\n";
            }
            registerCut(it->first, q, is_strict, is_strict ? "Old Ref Bound [STRICT]" : "Old Ref Bound [FLEXIBLE]");
        }

        for (auto it = qryBounds.lower_bound(check_q_start); it != qryBounds.end() && it->first <= check_q_end; ++it) {
            bool is_strict = (it->second.type == BoundaryType::STRICT);
            int r = mapQryToRef(best, it->first);
            
            if (debug && is_strict) {
                std::cout << "    ⚠️ [STRICT-OLD-QRY] Recovered Topological Boundary @ Qry: " << it->first 
                          << " (Reason: " << it->second.reason << ")\n";
            }
            registerCut(r, it->first, is_strict, is_strict ? "Old Qry Bound [STRICT]" : "Old Qry Bound [FLEXIBLE]");
        }

        // ========================================================
        // 🌟 階段 3：一維攤平掃描 (1D Flattened Array Sweep)
        // 完全捨棄複雜的迭代器操作，轉為最直觀的 1D Array 遍歷
        // ========================================================
        if (debug) std::cout << "  [DEBUG-STAGE 3] Processing 1D Flattened Cut Array...\n";

        // 萃取為 1D 陣列
        std::vector<std::pair<std::pair<int, int>, CutPoint>> flat_cuts(unified_cuts.begin(), unified_cuts.end());
        
        std::set<int> r_cuts, q_cuts;
        std::map<int, std::string> r_cut_reasons, q_cut_reasons;
        std::vector<std::pair<int, int>> active_cuts;

        // 🌟 紀錄雙切區域，利用 pair 的 (first < second) 來暗示交換方向
        std::vector<std::pair<int, int>> pending_special_swaps;

        int c_idx = 0;
        while (c_idx < flat_cuts.size()) {
            
            // ----------------------------------------------------
            // 1. 在 1D 陣列上圈出 Cluster (B 到 D)
            // ----------------------------------------------------
            int start_idx = c_idx;
            int end_idx = c_idx;
            
            while (end_idx + 1 < flat_cuts.size()) {
                int r_dist = std::abs(flat_cuts[end_idx + 1].first.first - flat_cuts[end_idx].first.first);
                int q_dist = std::abs(flat_cuts[end_idx + 1].first.second - flat_cuts[end_idx].first.second);
                
                if (std::max(r_dist, q_dist) > L_min) break;
                end_idx++;
            }

            // 🌟 核心新增：收集這個 Cluster 裡面的所有 STRICT 點
            std::vector<int> strict_in_cluster;
            for (int k = start_idx; k <= end_idx; ++k) {
                if (flat_cuts[k].second.is_strict) {
                    strict_in_cluster.push_back(k);
                }
            }

            // ====================================================
            // 🔍 [新增] Debug Message: 印出 B-D Cluster 內的所有候選點
            // ====================================================
            int cluster_size = end_idx - start_idx + 1;
            if (debug && cluster_size > 1) {
                // int cluster_size = end_idx - start_idx + 1;
                std::cout << "    -> 🎯 [CLUSTER DETECTED] B-D Range contains " << cluster_size << " candidate(s) [L_min: " << L_min << "]:\n";
                for (int k = start_idx; k <= end_idx; ++k) {
                    auto& cut = flat_cuts[k];
                    bool is_last = (k == end_idx);
                    std::cout << "       " << (is_last ? "└─ " : "├─ ") 
                              << "[" << (k - start_idx + 1) << "/" << cluster_size << "] "
                              << "Ref: " << cut.first.first << ", Qry: " << cut.first.second 
                              << " | Type: " << (cut.second.is_strict ? "STRICT" : "FLEXIBLE")
                              << " | Reason: " << cut.second.reasons << "\n";
                }
            }

            int B_R = flat_cuts[start_idx].first.first;
            int B_Q = flat_cuts[start_idx].first.second;
            int D_R = flat_cuts[end_idx].first.first;
            int D_Q = flat_cuts[end_idx].first.second;


            // ----------------------------------------------------
            // 2. 🌟 尋找外圍的安全隔離牆 (A_pair 和 E_pair) [1D 絕對空間版]
            // ----------------------------------------------------
            // std::pair<int, int> B_pair = flat_cuts[start_idx].first;
            // std::pair<int, int> D_pair = flat_cuts[end_idx].first;
            
            std::pair<int, int> A_pair = {check_r_start, best.inverse ? check_q_end : check_q_start};
            std::pair<int, int> E_pair = {check_r_end, best.inverse ? check_q_start : check_q_end};
            
            // 🧱 左牆 (A_pair)：因為 active_cuts 現在只存「過去已經切下的點」，直接拿最後一個就是最近的左牆！
            if (!active_cuts.empty()) {
                A_pair = active_cuts.back();
            }

            // 🧱 右牆 (E_pair)：從未來的 flat_cuts 中找下一個 Cluster 的起點！
            // 這樣可以保證沙盒完美地被限制在 local 區域，絕對不會跨越到下一個切點
            if (end_idx + 1 < flat_cuts.size()) {
                E_pair = flat_cuts[end_idx + 1].first;
            }

            int A_R = A_pair.first;
            int A_Q = A_pair.second;
            int E_R = E_pair.first;
            int E_Q = E_pair.second;

            // ----------------------------------------------------
            // 2.5 ⚡ 孤立點極速捷徑 (Isolated Check - 安全升級版)
            // ----------------------------------------------------
            if (start_idx == end_idx) {
                bool is_isolated = true;
                
                // 條件 1：與其他已註冊的切點 (active_cuts) 保持安全距離
                for (auto& c : active_cuts) {
                    if (std::abs(c.first - B_R) < L_min && std::abs(c.second - B_Q) < L_min) {
                        is_isolated = false; break;
                    }
                }

                // 🌟 條件 2：只與「Alignment 的絕對頭尾邊界」保持安全距離！
                // 放棄動態的 A / E，直接與這條序列最外圍的牆壁比較
                int dist_to_ref_start = std::abs(B_R - check_r_start);
                int dist_to_ref_end   = std::abs(check_r_end - B_R);
                int dist_to_qry_start = std::abs(B_Q - check_q_start);
                int dist_to_qry_end   = std::abs(check_q_end - B_Q);

                const int MIN_DIST_2_BOUNDARY = 50;

                // 如果離 Alignment 的頭或尾太近 (大於 0 且小於閾值)
                // 則剝奪 Fast-Path 資格，強制送進沙盒受審！
                if ((dist_to_ref_start > 0 && dist_to_ref_start < MIN_DIST_2_BOUNDARY) || 
                    (dist_to_ref_end > 0   && dist_to_ref_end < MIN_DIST_2_BOUNDARY) ||
                    (dist_to_qry_start > 0 && dist_to_qry_start < MIN_DIST_2_BOUNDARY) || 
                    (dist_to_qry_end > 0   && dist_to_qry_end < MIN_DIST_2_BOUNDARY)) {
                    is_isolated = false; 
                    if (debug) {
                        std::cout << "    -> 🛑 [FAST-PATH DENIED] Cut @ R:" << B_R << ", Q:" << B_Q 
                                  << " is too close to absolute alignment boundaries. Forcing Sandbox Evaluation.\n";
                    }
                }

                // 只有真正安全、遠離邊界與其他切點的點，才能走極速捷徑
                if (is_isolated) {
                    if (debug) std::cout << "    -> ⚡ [FAST-PATH] Isolated WGA Cut @ R:" << B_R << ", Q:" << B_Q << "\n";
                    r_cuts.insert(B_R);
                    q_cuts.insert(B_Q);
                    
                    std::string role = flat_cuts[start_idx].second.is_strict ? " [Forced STRICT]" : "";
                    r_cut_reasons[B_R] = flat_cuts[start_idx].second.reasons + " (Isolated)" + role;
                    q_cut_reasons[B_Q] = flat_cuts[start_idx].second.reasons + " (Isolated)" + role;
                    active_cuts.push_back({B_R, B_Q});
                    
                    c_idx = end_idx + 1;
                    continue; 
                }
            }

            // ----------------------------------------------------
            // 3. 🌟 動態 CIGAR 裁切與 A-E Micro-Sandbox 融合 (防崩潰安全版)
            // ----------------------------------------------------
            int curr_r = check_r_start;
            int curr_q = best.inverse ? check_q_end : check_q_start;
            int q_dir  = best.inverse ? -1 : 1;
            
            CigarString subCIGAR;
            int merged_consensus_len = 0;
            std::map<std::pair<int, int>, int> rq2merged; 
            
            // 嚴格遵守 Start < End 擷取原則
            int qry_box_start = std::min(A_Q, E_Q);
            int qry_box_end   = std::max(A_Q, E_Q);
            
            bool recording = false;
            // 邊界特判：如果 A_pair 剛好就是 Alignment 的最起點
            if (curr_r == A_R && curr_q == A_Q) recording = true;

            for (auto op : best.CIGAR) {
                // 如果已經走到 E 點卻還沒結束，強制跳出
                if (!recording && curr_r == E_R && curr_q == E_Q) break; 
                
                int len = op.first; char type = op.second;
                bool consumes_r = (type == 'M' || type == '=' || type == 'X' || type == 'D');
                bool consumes_q = (type == 'M' || type == '=' || type == 'X' || type == 'I');
                
                for (int i = 0; i < len; ++i) {
                    // 偵測是否抵達 A 點 (開始錄製)
                    if (!recording && curr_r == A_R && curr_q == A_Q) {
                        recording = true;
                    }
                    // 偵測是否抵達 E 點 (停止錄製並跳出)
                    if (recording && curr_r == E_R && curr_q == E_Q) {
                        recording = false;
                        break; 
                    }
                    
                    if (recording) {
                        rq2merged[{curr_r, curr_q}] = merged_consensus_len;
                        if (!subCIGAR.empty() && subCIGAR.back().second == type) {
                            subCIGAR.back().first++;
                        } else {
                            subCIGAR.push_back({1, type});
                        }
                        merged_consensus_len++;
                    }
                    
                    if (consumes_r) curr_r++;
                    if (consumes_q) curr_q += q_dir;
                }
                
                // 區段結束後，再次檢查是否剛好踩在 E_pair 邊界
                if (recording && curr_r == E_R && curr_q == E_Q) {
                    recording = false;
                    break;
                }
            }
            // 確保終點有被 Map 到，供後續安全獲取相對座標
            rq2merged[{E_R, E_Q}] = merged_consensus_len;

            if (debug) {
                std::cout << "    -> 📦 Extracting and Merging A-E Sandbox:\n"
                          << "       ├─ Ref Bounds: [" << A_R << " -> " << E_R << "] (Len: " << (E_R - A_R) << ")\n"
                          << "       ├─ Qry Bounds: [" << qry_box_start << " -> " << qry_box_end << "] (Len: " << (qry_box_end - qry_box_start) << ")\n"
                          << "       └─ Merged Consensus Length: " << merged_consensus_len << "\n";
            }
            
            BlockSet tempBlockSet("temp");
            auto localRefBlock = extractBlockFromSuper(&tempBlockSet, refSuperBlock, A_R, E_R);
            auto localQryBlock = extractBlockFromSuper(&tempBlockSet, qrySuperBlock, qry_box_start, qry_box_end);
            
            // 🌟 完美融合！CIGAR 消耗的 Ref/Qry 數量保證與抽出 Block 的實體長度 100% 吻合！
            auto mergedLocalBlock = tempBlockSet.mergeTwoBlocks(localRefBlock, localQryBlock, subCIGAR, best.inverse);

            // 輔助函式：安全獲取全域 (R,Q) 在 Merged Block 的 Local 座標
            auto getMergedCoord = [&](int r, int q) {
                auto it = rq2merged.find({r, q});
                if (it != rq2merged.end()) return it->second;
                for (auto const& [rq, mx] : rq2merged) if (rq.first >= r) return mx;
                return merged_consensus_len;
            };

            int local_B_merged = getMergedCoord(B_R, B_Q);
            int local_D_merged = getMergedCoord(D_R, D_Q);
            if (local_B_merged > local_D_merged) std::swap(local_B_merged, local_D_merged);
            
            // 呼叫評分，給定 B-D 在 Merged Consensus 上的範圍
            std::vector<Block::ColumnSplitScore> merged_scores = mergedLocalBlock->calculateSplittingScores(local_B_merged, local_D_merged);

            // ----------------------------------------------------
            // 4. 在 1D Array 範圍內執行貪婪掃描
            // ----------------------------------------------------
            while (true) {
                int best_k = -1;
                int best_cut = -1;
                int best_q_cut = -1;
                
                // 🌟 修改：起始值設為極小值，才能捕捉並比較負分
                double max_score = -1e9; 
                
                double best_perfect_bonus = 0.0;
                double best_score_L = 0.0, best_score_R = 0.0;

                // 🌟 新增：紀錄這輪存活的所有合法 Candidate，供特例判定使用
                struct Cand { int k; double score; };
                std::vector<Cand> valid_cands;

                // 遍歷當前 Cluster 內的所有 1D 節點
                for (int k = start_idx; k <= end_idx; ++k) {
                    int x_global = flat_cuts[k].first.first;
                    int qx_global = flat_cuts[k].first.second;

                    // 檢查是否已被防護網擋下
                    bool valid = true;
                    int cur_L_R = check_r_start, cur_R_R = check_r_end;
                    int cur_L_Q = check_q_start, cur_R_Q = check_q_end;

                    for (auto& c : active_cuts) {
                        int r_dist = std::abs(c.first - x_global);
                        int q_dist = std::abs(c.second - qx_global);

                        if (r_dist < L_min && q_dist < L_min) { valid = false; break; }

                        if (c.first <= x_global && c.first > cur_L_R) cur_L_R = c.first;
                        if (c.first >= x_global && c.first < cur_R_R) cur_R_R = c.first;

                        if (c.second <= qx_global && c.second > cur_L_Q) cur_L_Q = c.second;
                        if (c.second >= qx_global && c.second < cur_R_Q) cur_R_Q = c.second;
                    }
                    if (!valid) continue;

                    // 獲取精準的 1D Consensus 座標
                    int merged_x = getMergedCoord(x_global, qx_global);
                    int score_idx = merged_x - local_B_merged;
                    
                    if (score_idx < 0 || score_idx >= merged_scores.size()) continue;

                    auto s_merged = merged_scores[score_idx];
                    
                    double total_perfect = s_merged.perfect_bonus;
                    double score_L = s_merged.id_left;
                    double score_R = s_merged.id_right;

                    int L_len_R = x_global - cur_L_R;
                    int R_len_R = cur_R_R - x_global;
                    int L_len_Q = std::abs(qx_global - cur_L_Q);
                    int R_len_Q = std::abs(cur_R_Q - qx_global);

                    if (std::max(L_len_R, L_len_Q) < L_min) score_L += 10000;
                    if (std::max(R_len_R, R_len_Q) < L_min) score_R += 10000;

                    double final_score = total_perfect - ((score_L + score_R)/2);

                    // 儲存合法的候選點
                    valid_cands.push_back({k, final_score});

                    if (debug) {
                        std::cout << "      [CANDIDATE] R:" << x_global << " Q:" << qx_global 
                                  << " -> Final: " << (final_score < -5000 ? "-INF" : std::to_string(final_score)) << "\n"
                                  << "        ├─ Perfect Bonus: " << total_perfect << "\n"
                                  << "        ├─ Left Side  -> MAD Penalty Score: " << (score_L > 5000 ? "-INF (Killed by L_min)" : std::to_string(score_L)) << "\n"
                                  << "        └─ Right Side -> MAD Penalty Score: " << (score_R > 5000 ? "-INF (Killed by L_min)" : std::to_string(score_R)) << "\n";
                    }

                    if (final_score > max_score) {
                        max_score = final_score;
                        best_k = k;
                        best_cut = x_global;
                        best_q_cut = qx_global;
                    }
                }

                if (valid_cands.empty()) break; // 防護網全擋，結束掃描

                // ====================================================
                // 🌟 核心決策區：Score > 0 門檻 與 特例雙切防護
                // ====================================================
                bool is_special_dual_cut = false;
                int top1_k = -1;
                int top2_k = -1;

                // 🌟 泛用化特例判定：只要存活 >= 2 個 Candidate，且最高分 < 0
                if (valid_cands.size() >= 2 && max_score < 0) {
                    
                    // 將候選點依照分數由高到低排序，抓出前兩名
                    std::vector<Cand> sorted_cands = valid_cands;
                    std::sort(sorted_cands.begin(), sorted_cands.end(), [](const Cand& a, const Cand& b) {
                        return a.score > b.score; // 降冪排序
                    });

                    top1_k = sorted_cands[0].k;
                    top2_k = sorted_cands[1].k;

                    std::string reason1 = flat_cuts[top1_k].second.reasons;
                    std::string reason2 = flat_cuts[top2_k].second.reasons;

                    // 檢查這「最高分的兩個切點」是否剛好包含 Old 與 Long 衝突
                    bool has_old = (reason1.find("Old") != std::string::npos || reason2.find("Old") != std::string::npos);
                    bool has_long = (reason1.find("Long") != std::string::npos || reason2.find("Long") != std::string::npos);

                    if (has_old && has_long) {
                        is_special_dual_cut = true;
                    }
                }

                std::set<int> winners_to_cut;

                if (is_special_dual_cut) {
                    int cut1_R = flat_cuts[top1_k].first.first;
                    int cut1_Q = flat_cuts[top1_k].first.second;
                    int cut2_R = flat_cuts[top2_k].first.first;
                    int cut2_Q = flat_cuts[top2_k].first.second;

                    // 🌟 檢查前兩名是否包含 STRICT 邊界
                    bool has_strict = flat_cuts[top1_k].second.is_strict || 
                                      flat_cuts[top2_k].second.is_strict;

                    // 從 Reason 推導方向
                    std::string reason1 = flat_cuts[top1_k].second.reasons;
                    std::string reason2 = flat_cuts[top2_k].second.reasons;
                    bool swap_with_left = (reason1.find("End") != std::string::npos || reason2.find("End") != std::string::npos);

                    if (debug) {
                        std::cout << "      ---> 🚨 [SPECIAL DUAL-CUT] Detected Old Bound + Long Indel collision with scores < 0.\n"
                                  << "           Cutting TOP 2 scoring bounds to isolate the micro-fragment and delegating to Phase 5!\n"
                                  << "           ├─ Cut 1 (Top 1): R:" << cut1_R << " Q:" << cut1_Q << "\n"
                                  << "           ├─ Cut 2 (Top 2): R:" << cut2_R << " Q:" << cut2_Q << "\n";
                                  
                        if (has_strict) {
                            std::cout << "           └─ Action: Contains STRICT bound. Cutting BOTH and delegating directly to Phase 5 (No Swap).\n";
                        } else {
                            std::cout << "           └─ Action: Cutting BOTH bounds. Will swap with " << (swap_with_left ? "LEFT" : "RIGHT") << " fragment in Phase 4.5.\n";
                        }
                    }
                    
                    // 兩刀都必須下，以隔離微碎片
                    winners_to_cut.insert(top1_k);
                    winners_to_cut.insert(top2_k);

                    // 🌟 只有在「兩刀都不是 STRICT」的情況下，才紀錄並交給 4.5 做 Swap
                    if (!has_strict) {
                        int r_start = std::min(cut1_R, cut2_R);
                        int r_end   = std::max(cut1_R, cut2_R);
                        
                        // 利用 pair 的先後順序偷藏方向資訊
                        if (swap_with_left) pending_special_swaps.push_back({r_start, r_end});
                        else                pending_special_swaps.push_back({r_end, r_start});
                    }
                }
                else if (max_score > 0) {
                    if (debug) {
                        std::cout << "      ---> 🎉 [WINNER] Snapped to R:" << best_cut << " Q:" << best_q_cut 
                                  << " with Max Score: " << max_score << "\n";
                    }
                    winners_to_cut.insert(best_k); // 正常贏家
                } 
                else {
                    if (debug) {
                        std::cout << "      ---> 🛑 [REJECTED] Best score (" << max_score << ") <= 0. Skipping to avoid random fragmentation.\n";
                    }
                }

                // ⚠️ 圖譜最後底線：即使被 Rejected，如果是 STRICT 也必須強制下刀
                for (int s_idx : strict_in_cluster) {
                    bool is_valid = false;
                    for (auto& c : valid_cands) if (c.k == s_idx) is_valid = true;
                    
                    if (is_valid && winners_to_cut.find(s_idx) == winners_to_cut.end()) {
                        winners_to_cut.insert(s_idx);
                        if (debug) std::cout << "      ---> ⚠️ [STRICT OVERRIDE] Forced cutting STRICT bound to prevent graph topology collapse.\n";
                    }
                }

                // 將選出的 Winners 實際加入 Active Cuts
                if (!winners_to_cut.empty()) {
                    for (int idx : winners_to_cut) {
                        int cut_R = flat_cuts[idx].first.first;
                        int cut_Q = flat_cuts[idx].first.second;

                        r_cuts.insert(cut_R);
                        q_cuts.insert(cut_Q);

                        // 幫它們貼上專屬的標籤
                        std::string role = " [Score > 0 Winner]";
                        if (is_special_dual_cut) role = " [Special Dual-Cut for Phase 5]";
                        else if (flat_cuts[idx].second.is_strict && idx != best_k) role = " [Forced STRICT Override]";

                        r_cut_reasons[cut_R] = flat_cuts[idx].second.reasons + role;
                        q_cut_reasons[cut_Q] = flat_cuts[idx].second.reasons + role;
                        active_cuts.push_back({cut_R, cut_Q});
                    }
                } else {
                    break; // 這輪沒有切下任何一刀，跳出 while 貪婪迴圈
                }
            }
            c_idx = end_idx + 1;
        }

        // ========================================================
        // 後續步驟 (階段 4 初步切割 -> 階段 5 錯位補償 -> 階段 6 封裝)
        // 邏輯保持原樣，因為你的架構防護網已經建得非常穩健！
        // ========================================================
        std::vector<Alignment> frags;
        if (!r_cuts.empty() || !q_cuts.empty()) {
            frags = splitSingleAlignment(best, r_cuts, q_cuts); // q_cuts 現在由 addCut 自動處理好了
        } else {
            frags.push_back(best);
        }

        // ========================================================
        // 🌟 步驟 4.5：精確鎖定特例雙切碎片，依照解碼指令強制置換
        // ========================================================
        if (!pending_special_swaps.empty() && frags.size() >= 2) {
            for (auto& task : pending_special_swaps) {
                
                // 🌟 1. 解碼：判斷方向與真實座標
                bool swap_with_left = (task.first < task.second);
                int target_r_start = std::min(task.first, task.second);
                int target_r_end   = std::max(task.first, task.second);
                
                // 2. 尋找這個被雙切獨立出來的「微碎片」
                int micro_idx = -1;
                for (size_t i = 0; i < frags.size(); ++i) {
                    int f_r_start = std::min(frags[i].refIdx.first, frags[i].refIdx.second);
                    int f_r_end   = std::max(frags[i].refIdx.first, frags[i].refIdx.second);
                    
                    if (f_r_start == target_r_start && f_r_end == target_r_end && (f_r_end - f_r_start) > 0) {
                        micro_idx = i;
                        break;
                    }
                }

                if (micro_idx != -1) {
                    // 3. 套用解碼出的方向
                    int gap_idx = swap_with_left ? (micro_idx - 1) : (micro_idx + 1);
                    
                    if (gap_idx >= 0 && gap_idx < frags.size()) {
                        int first_idx = std::min(micro_idx, gap_idx);
                        int second_idx = std::max(micro_idx, gap_idx);
                        auto& f1 = frags[first_idx];
                        auto& f2 = frags[second_idx];

                        if (debug) {
                            std::cout << "  -> 🔀 [FRAG-SWAP] Exact Match for Dual-Cut Region R:[" << target_r_start << ", " << target_r_end << "]\n"
                                      << "       ├─ Order decoded: Swapping Frag " << micro_idx << " with Frag " << gap_idx << "\n"
                                      << "       ├─ [PRE-SWAP] Frag " << first_idx << " | CIGAR: "; printCIGAR(f1.CIGAR, false); 
                            std::cout << " | Ref: [" << f1.refIdx.first << ", " << f1.refIdx.second << "] | Qry: [" << f1.qryIdx.first << ", " << f1.qryIdx.second << "]\n"
                                      << "       └─ [PRE-SWAP] Frag " << second_idx << " | CIGAR: "; printCIGAR(f2.CIGAR, false); 
                            std::cout << " | Ref: [" << f2.refIdx.first << ", " << f2.refIdx.second << "] | Qry: [" << f2.qryIdx.first << ", " << f2.qryIdx.second << "]\n";
                        }

                        // 靈魂交換
                        std::swap(f1.CIGAR, f2.CIGAR);

                        // 重新計算座標
                        int r_shift = 0, q_shift = 0;
                        for (auto& op : f1.CIGAR) {
                            if (op.second == 'M' || op.second == '=' || op.second == 'X') { r_shift += op.first; q_shift += op.first; }
                            else if (op.second == 'D') { r_shift += op.first; }
                            else if (op.second == 'I') { q_shift += op.first; }
                        }

                        f1.refIdx.second = f1.refIdx.first + r_shift;
                        if (!f1.inverse) f1.qryIdx.second = f1.qryIdx.first + q_shift;
                        else             f1.qryIdx.second = f1.qryIdx.first - q_shift;

                        // f2 的起點無縫接住 f1 的終點
                        f2.refIdx.first = f1.refIdx.second;
                        f2.qryIdx.first = f1.qryIdx.second;

                        f1.updateAlnLength();
                        f2.updateAlnLength();

                        if (debug) {
                            std::cout << "       ├─ [POST-SWAP] Frag " << first_idx << " | CIGAR: "; printCIGAR(f1.CIGAR, false); 
                            std::cout << " | Ref: [" << f1.refIdx.first << ", " << f1.refIdx.second << "] | Qry: [" << f1.qryIdx.first << ", " << f1.qryIdx.second << "]\n"
                                      << "       └─ [POST-SWAP] Frag " << second_idx << " | CIGAR: "; printCIGAR(f2.CIGAR, false); 
                            std::cout << " | Ref: [" << f2.refIdx.first << ", " << f2.refIdx.second << "] | Qry: [" << f2.qryIdx.first << ", " << f2.qryIdx.second << "]\n";
                        }
                    }
                }
            }
        }

        // ========================================================
        // 🌟 步驟 5：非對稱單軸剝離與錯位補償制 (包含頭尾 STRICT 修剪)
        // ========================================================
        std::vector<bool> absorbed(frags.size(), false);

        auto isStrictCut = [&](int r, int q) -> bool {
            if (r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos) return true;
            if (q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos) return true;
            return false;
        };

        if (frags.size() >= 3) {
            for (size_t i = 1; i < frags.size() - 1; ++i) {
                int f_q_len = std::abs(frags[i].qryIdx.second - frags[i].qryIdx.first);
                int f_r_len = std::abs(frags[i].refIdx.second - frags[i].refIdx.first);

                // 提前取得實體座標，供 STRICT 判定使用
                int r_left_cut  = frags[i].refIdx.first;  int r_right_cut = frags[i].refIdx.second;
                int q_left_cut  = frags[i].qryIdx.first;  int q_right_cut = frags[i].qryIdx.second;

                auto isStrictR = [&](int r) { return r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos; };
                auto isStrictQ = [&](int q) { return q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos; };

                // 🌟 1. 雙 STRICT 豁免判定：只要任一軸被雙 STRICT 夾擊，絕對保留，跳過吸收！
                bool r_both_strict = isStrictR(r_left_cut) && isStrictR(r_right_cut);
                bool q_both_strict = isStrictQ(q_left_cut) && isStrictQ(q_right_cut);

                if (r_both_strict || q_both_strict) {
                    if (debug) std::cout << "      -> 🛡️ [STRICT-SHIELD] Fragment " << i << " is trapped between STRICT bounds. Absorption skipped!\n";
                    continue;
                }

                // 🌟 核心修正：分開判斷單軸是否為微小碎片
                bool is_q_micro = (f_q_len > 0 && f_q_len < 10);
                bool is_r_micro = (f_r_len > 0 && f_r_len < 10);

                // 只要有任何一軸是微碎片，就啟動剝離程序
                if (is_q_micro || is_r_micro) {
                    
                    bool is_pure_match = (frags[i].CIGAR.size() == 1 && frags[i].CIGAR[0].second == 'M');
                    if (is_pure_match && f_r_len >= 30) continue; 

                    if (!absorbed[i-1] && !absorbed[i+1]) {
                        
                        // 1. 計算左右親和力 (維持原樣)
                        int left_r_len = std::abs(frags[i-1].refIdx.second - frags[i-1].refIdx.first);
                        int left_q_len = std::abs(frags[i-1].qryIdx.second - frags[i-1].qryIdx.first);
                        int right_r_len = std::abs(frags[i+1].refIdx.second - frags[i+1].refIdx.first);
                        int right_q_len = std::abs(frags[i+1].qryIdx.second - frags[i+1].qryIdx.first);

                        auto isStrictR = [&](int r) { return r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos; };
                        auto isStrictQ = [&](int q) { return q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos; };

                        int r_left_cut  = frags[i].refIdx.first;  int r_right_cut = frags[i].refIdx.second;
                        int q_left_cut  = frags[i].qryIdx.first;  int q_right_cut = frags[i].qryIdx.second;

                        // 2. 獨立決策 (維持原樣)
                        bool ref_goes_left;
                        if (isStrictR(r_right_cut) && !isStrictR(r_left_cut)) ref_goes_left = (left_r_len > 0); 
                        else if (isStrictR(r_left_cut) && !isStrictR(r_right_cut)) ref_goes_left = !(right_r_len > 0); 
                        else ref_goes_left = (left_r_len >= right_r_len); 

                        bool qry_goes_left;
                        if (isStrictQ(q_right_cut) && !isStrictQ(q_left_cut)) qry_goes_left = (left_q_len > 0);  
                        else if (isStrictQ(q_left_cut) && !isStrictQ(q_right_cut)) qry_goes_left = !(right_q_len > 0); 
                        else qry_goes_left = (left_q_len >= right_q_len);

                        // 🌟 3. 單軸轉移量計算：只移動被判定為「微小」的那一軸！
                        int left_append_r = (is_r_micro && ref_goes_left) ? f_r_len : 0;
                        int right_prepend_r = (is_r_micro && !ref_goes_left) ? f_r_len : 0;
                        
                        int left_append_q = (is_q_micro && qry_goes_left) ? f_q_len : 0;
                        int right_prepend_q = (is_q_micro && !qry_goes_left) ? f_q_len : 0;

                        // 4. 變形並裝載至左側 (Left Block)
                        if (left_append_r > 0 && left_append_q > 0) frags[i-1].CIGAR.insert(frags[i-1].CIGAR.end(), frags[i].CIGAR.begin(), frags[i].CIGAR.end());
                        else if (left_append_r > 0) frags[i-1].CIGAR.push_back({left_append_r, 'D'});
                        else if (left_append_q > 0) frags[i-1].CIGAR.push_back({left_append_q, 'I'});

                        if (left_append_r > 0) frags[i-1].refIdx.second = frags[i].refIdx.second;
                        if (left_append_q > 0) {
                            if (!frags[i-1].inverse) frags[i-1].qryIdx.second = frags[i].qryIdx.second;
                            else                     frags[i-1].qryIdx.first  = frags[i].qryIdx.first;
                        }
                        frags[i-1].updateAlnLength();

                        // 5. 變形並裝載至右側 (Right Block)
                        if (right_prepend_r > 0 && right_prepend_q > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), frags[i].CIGAR.begin(), frags[i].CIGAR.end());
                        else if (right_prepend_r > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), {right_prepend_r, 'D'});
                        else if (right_prepend_q > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), {right_prepend_q, 'I'});

                        if (right_prepend_r > 0) frags[i+1].refIdx.first = frags[i].refIdx.first;
                        if (right_prepend_q > 0) {
                            if (!frags[i+1].inverse) frags[i+1].qryIdx.first = frags[i].qryIdx.first;
                            else                     frags[i+1].qryIdx.second = frags[i].qryIdx.second;
                        }
                        frags[i+1].updateAlnLength();

                        // 🌟 6. 決定本塊碎肉的生死 (Partial Absorption)
                        if (is_r_micro && is_q_micro) {
                            absorbed[i] = true; // 兩軸都太小，完全蒸發
                            if (debug) std::cout << "      -> 💥 [FULL-ABSORB] Fragment " << i << " fully absorbed. "
                                                 << "Ref: (" << frags[i].refIdx.first << "," << frags[i].refIdx.second << "), "
                                                 << "Qry: (" << frags[i].qryIdx.first << "," << frags[i].qryIdx.second << ").\n";
                        } else {
                            // 單軸剝離！保留宏觀 (Macro) 的那一軸，將它淨化為純 Gap Block
                            if (is_r_micro) { 
                                frags[i].CIGAR = {{f_q_len, 'I'}}; // 淨化為純 Insertion
                                if (ref_goes_left) frags[i].refIdx.first = frags[i].refIdx.second;
                                else               frags[i].refIdx.second = frags[i].refIdx.first;
                                if (debug) std::cout << "      -> ✨ [PURIFY] Fragment " << i << " stripped of Ref noise, became pure " << f_q_len << "I.\n";
                            } else if (is_q_micro) { 
                                frags[i].CIGAR = {{f_r_len, 'D'}}; // 淨化為純 Deletion
                                if (qry_goes_left) {
                                    if (!frags[i].inverse) frags[i].qryIdx.first = frags[i].qryIdx.second;
                                    else                   frags[i].qryIdx.second = frags[i].qryIdx.first;
                                } else {
                                    if (!frags[i].inverse) frags[i].qryIdx.second = frags[i].qryIdx.first;
                                    else                   frags[i].qryIdx.first = frags[i].qryIdx.second;
                                }
                                if (debug) std::cout << "      -> ✨ [PURIFY] Fragment " << i << " stripped of Qry noise, became pure " << f_r_len << "D.\n";
                            }
                            frags[i].updateAlnLength();
                        }
                    }
                }
            }
            
            // 🌟 處理「頭部碎肉 (Prefix Micro-fragment)」
            int f0_q_len = std::abs(frags[0].qryIdx.second - frags[0].qryIdx.first);
            int f0_r_len = std::abs(frags[0].refIdx.second - frags[0].refIdx.first);
            bool f0_q_micro = (f0_q_len > 0 && f0_q_len < L_min);
            bool f0_r_micro = (f0_r_len > 0 && f0_r_len < L_min);

            if (!absorbed[0] && (f0_q_micro || f0_r_micro)) {
                bool is_pure_match = (frags[0].CIGAR.size() == 1 && frags[0].CIGAR[0].second == 'M');
                if (!(is_pure_match && f0_r_len >= 30)) {
                    
                    int cut_r = frags[0].refIdx.second;
                    int cut_q = frags[0].inverse ? frags[0].qryIdx.first : frags[0].qryIdx.second;

                    if (isStrictCut(cut_r, cut_q)) {
                        if (f0_q_micro && f0_r_micro) absorbed[0] = true;
                        else {
                            if (f0_r_micro) { frags[0].CIGAR = {{f0_q_len, 'I'}}; frags[0].refIdx.second = frags[0].refIdx.first; }
                            else if (f0_q_micro) {
                                frags[0].CIGAR = {{f0_r_len, 'D'}};
                                if (!frags[0].inverse) frags[0].qryIdx.second = frags[0].qryIdx.first;
                                else                   frags[0].qryIdx.first = frags[0].qryIdx.second;
                            }
                        }
                        if (debug && absorbed[0]) std::cout << "      -> ✂️ [TERMINAL-TRIM] Prefix fully vaporized.\n";
                    } else {
                        int move_r = f0_r_micro ? f0_r_len : 0;
                        int move_q = f0_q_micro ? f0_q_len : 0;
                        
                        if (move_r > 0 && move_q > 0) frags[1].CIGAR.insert(frags[1].CIGAR.begin(), frags[0].CIGAR.begin(), frags[0].CIGAR.end());
                        else if (move_r > 0) frags[1].CIGAR.insert(frags[1].CIGAR.begin(), {move_r, 'D'});
                        else if (move_q > 0) frags[1].CIGAR.insert(frags[1].CIGAR.begin(), {move_q, 'I'});

                        if (move_r > 0) frags[1].refIdx.first = frags[0].refIdx.first;
                        if (move_q > 0) {
                            if (!frags[1].inverse) frags[1].qryIdx.first = frags[0].qryIdx.first;
                            else                   frags[1].qryIdx.second = frags[0].qryIdx.second;
                        }
                        frags[1].updateAlnLength();

                        if (f0_q_micro && f0_r_micro) absorbed[0] = true;
                        else {
                            if (f0_r_micro) { frags[0].CIGAR = {{f0_q_len, 'I'}}; frags[0].refIdx.first = frags[0].refIdx.second; }
                            else if (f0_q_micro) {
                                frags[0].CIGAR = {{f0_r_len, 'D'}};
                                if (!frags[0].inverse) frags[0].qryIdx.first = frags[0].qryIdx.second;
                                else                   frags[0].qryIdx.second = frags[0].qryIdx.first;
                            }
                        }
                        if (debug && absorbed[0]) std::cout << "      -> 💥 [PREFIX-ABSORB] Fully absorbed into Block 1.\n";
                    }
                    frags[0].updateAlnLength();
                }
            }

            // 🌟 處理「尾部碎肉 (Suffix Micro-fragment)」
            size_t last = frags.size() - 1;
            int fl_q_len = std::abs(frags[last].qryIdx.second - frags[last].qryIdx.first);
            int fl_r_len = std::abs(frags[last].refIdx.second - frags[last].refIdx.first);
            bool fl_q_micro = (fl_q_len > 0 && fl_q_len < L_min);
            bool fl_r_micro = (fl_r_len > 0 && fl_r_len < L_min);

            if (!absorbed[last] && (fl_q_micro || fl_r_micro)) {
                bool is_pure_match = (frags[last].CIGAR.size() == 1 && frags[last].CIGAR[0].second == 'M');
                if (!(is_pure_match && fl_r_len >= 30)) {

                    int cut_r = frags[last].refIdx.first;
                    int cut_q = frags[last].inverse ? frags[last].qryIdx.second : frags[last].qryIdx.first;

                    if (isStrictCut(cut_r, cut_q)) {
                        if (fl_q_micro && fl_r_micro) absorbed[last] = true;
                        else {
                            if (fl_r_micro) { frags[last].CIGAR = {{fl_q_len, 'I'}}; frags[last].refIdx.first = frags[last].refIdx.second; }
                            else if (fl_q_micro) {
                                frags[last].CIGAR = {{fl_r_len, 'D'}};
                                if (!frags[last].inverse) frags[last].qryIdx.first = frags[last].qryIdx.second;
                                else                      frags[last].qryIdx.second = frags[last].qryIdx.first;
                            }
                        }
                        if (debug && absorbed[last]) std::cout << "      -> ✂️ [TERMINAL-TRIM] Suffix fully vaporized.\n";
                    } else {
                        int move_r = fl_r_micro ? fl_r_len : 0;
                        int move_q = fl_q_micro ? fl_q_len : 0;

                        if (move_r > 0 && move_q > 0) frags[last-1].CIGAR.insert(frags[last-1].CIGAR.end(), frags[last].CIGAR.begin(), frags[last].CIGAR.end());
                        else if (move_r > 0) frags[last-1].CIGAR.push_back({move_r, 'D'});
                        else if (move_q > 0) frags[last-1].CIGAR.push_back({move_q, 'I'});

                        if (move_r > 0) frags[last-1].refIdx.second = frags[last].refIdx.second;
                        if (move_q > 0) {
                            if (!frags[last-1].inverse) frags[last-1].qryIdx.second = frags[last].qryIdx.second;
                            else                        frags[last-1].qryIdx.first = frags[last].qryIdx.first;
                        }
                        frags[last-1].updateAlnLength();

                        if (fl_q_micro && fl_r_micro) absorbed[last] = true;
                        else {
                            if (fl_r_micro) { frags[last].CIGAR = {{fl_q_len, 'I'}}; frags[last].refIdx.second = frags[last].refIdx.first; }
                            else if (fl_q_micro) {
                                frags[last].CIGAR = {{fl_r_len, 'D'}};
                                if (!frags[last].inverse) frags[last].qryIdx.second = frags[last].qryIdx.first;
                                else                      frags[last].qryIdx.first = frags[last].qryIdx.second;
                            }
                        }
                        if (debug && absorbed[last]) std::cout << "      -> 💥 [SUFFIX-ABSORB] Fully absorbed into Block " << last-1 << ".\n";
                    }
                    frags[last].updateAlnLength();
                }
            }
        }
        
        // ========================================================
        // 🌟 步驟 6：最終過濾、Family 綁定與封裝輸出
        // ========================================================
        for (size_t i = 0; i < frags.size(); ++i) {
            if (absorbed[i]) continue; // 被步驟 5 錯位補償吸收掉的，略過

            // 取得當前碎片的真實物理範圍
            int check_f_q_start = std::min(frags[i].qryIdx.first, frags[i].qryIdx.second);
            int check_f_q_end   = std::max(frags[i].qryIdx.first, frags[i].qryIdx.second);
            int check_f_r_start = frags[i].refIdx.first;
            int check_f_r_end   = frags[i].refIdx.second;

            // ==========================================
            // ❌ 已刪除危險的 Boundary-based Padding 邏輯
            // (Snapping 已經由 STAGE 0 與下方的 finalizeAlignment 基於 Coverage 完美處理)
            // ==========================================

            // 計算中心點，用來繼承舊圖譜的 Family ID
            int mid_r = (frags[i].refIdx.first + frags[i].refIdx.second) / 2;
            int mid_q = (check_f_q_start + check_f_q_end) / 2;

            frags[i].refFamilyId = getFamilyId(mid_r, refBounds);
            frags[i].qryFamilyId = getFamilyId(mid_q, qryBounds);

            if (debug) {
                if (frags[i].refFamilyId != 0 || frags[i].qryFamilyId != 0) {
                    std::cout << "      -> 🧬 [FAMILY-BOND] RefFam: " << frags[i].refFamilyId 
                              << " <=> QryFam: " << frags[i].qryFamilyId << "\n";
                }
            }

            // 呼叫輔助函數，裡面內建了最後一道安全的 Coverage Snapping 防線
            finalizeAlignment(frags[i]);
            frags[i].updateEnergy(refBlockSet, qryBlockSet);
            final_results.push_back(std::move(frags[i]));
        }

        if (!final_results.empty()) {
            if (debug) std::cout << "    ✅ [SUCCESS] Returning " << final_results.size() << " fully stabilized alignments.\n";
            return final_results;
        }
    }

    if (debug) std::cout << "\n[DEBUG] ⚠️ QUEUE EMPTY!\n";
    return final_results; 
}

// =====================================
// Alignment
// =====================================

void Alignment::updateAlnLength() {
    alnLength = 0;
    for (const auto& op : CIGAR) {
        alnLength += op.first; // 取出 pair 中的 int (長度) 並累加
    }
}
