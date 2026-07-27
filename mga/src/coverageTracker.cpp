#include "alignment.hpp"
#include "coordinate_manager.hpp"
#include "type.hpp"
#include <map>
#include <vector>
#include <set>
#include <algorithm>

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

// 🌟 升級版：直接從 CoordinateManager 的快速區間樹同步 Coverage
void CoverageTracker::syncFromMap(const CoordinateManager& coordMgr, bool isRefAxis) {
    this->intervals.clear(); 
    
    // 1. 根據目前追蹤的是 Ref 還是 Qry 軸，直接拉取對應的紅黑樹區間
    // 💡 註：請確保你在 CoordinateManager 類別中，有為 refIntervals 與 qryIntervals 提供 getter 
    //    或者將 CoverageTracker 設為 CoordinateManager 的 friend class
    const std::map<int, BlockInterval>& sourceIntervals = isRefAxis ? coordMgr.getRefIntervals() : coordMgr.getQryIntervals();

    // 2. 直接進行樹對樹的無縫平移投影，完全免去 O(N) 一維掃描的開銷
    for (const auto& [startPos, intervalInfo] : sourceIntervals) {
        // 排除未分配的空區塊 (ID = 0)
        if (intervalInfo.blkId == 0) continue;

        // 轉移至 Coverage 系統：intervals[start] = {end, blockId}
        this->intervals[startPos] = {intervalInfo.end, intervalInfo.blkId};
    }
}
