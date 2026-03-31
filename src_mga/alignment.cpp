#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <list>
#include <fstream>

struct PrimaryTracker {
    struct SegmentInfo {
        int end;
        int alnId; // 記錄這塊地盤屬於哪一個 Primary Alignment
    };
    
    // key: start_coordinate, value: {end_coordinate, alnId}
    std::map<int, SegmentInfo> intervals; 

    // 新增 Primary 區間 (保證不重疊，直接插入)
    void add(int start, int end, int alnId) {
        if (start >= end) return;
        intervals[start] = {end, alnId};
    }

    // 取得所有涵蓋 [start, end) 範圍內的 Primary IDs
    std::set<int> getOverlappingIds(int qStart, int qEnd) const {
        std::set<int> overlappingIds;
        if (qStart >= qEnd) return overlappingIds;

        // 找到第一個 start > qStart 的區間
        auto it = intervals.upper_bound(qStart);
        
        // 往前退一步，檢查前一個區間是否跨越了 qStart
        if (it != intervals.begin()) {
            auto prev = std::prev(it);
            if (prev->second.end > qStart) {
                overlappingIds.insert(prev->second.alnId);
            }
        }

        // 繼續往後找，直到區間的起點超出 qEnd
        while (it != intervals.end() && it->first < qEnd) {
            overlappingIds.insert(it->second.alnId);
            ++it;
        }

        return overlappingIds;
    }

    // 取得切點 (供 splitSingleAlignment 使用)
    void getCuts(int start, int end, std::set<int>& cuts) const {
        auto it = intervals.upper_bound(start);
        if (it != intervals.begin()) {
            auto prev = std::prev(it);
            if (prev->second.end > start && prev->second.end < end) {
                cuts.insert(prev->second.end);
            }
        }
        while (it != intervals.end() && it->first < end) {
            if (it->first > start) cuts.insert(it->first);
            if (it->second.end < end) cuts.insert(it->second.end);
            ++it;
        }
    }
    
    // 檢查是否被覆蓋 (保留原功能)
    bool isCovered(int start, int end) const {
        return !getOverlappingIds(start, end).empty();
    }
};


void printCoordinate(int st, int en) {
    std::cerr << "(" << st << "," << en << "]";
}

void printCoordinate(std::pair<int,int> r) {
    std::cerr << "(" << r.first << "," << r.second << "]";
}

// Helper to check if a CIGAR op consumes Reference
bool consumesRef(char op) {
    return (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X');
}

// Helper to check if a CIGAR op consumes Query
bool consumesQry(char op) {
    return (op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X');
}

mga::AlnChain::AlnChain(int id, float score, IntVec& chain) {
    this->identifier = id;
    this->score = score;
    this->chainedAln = chain;
}

void mga::Alignment::show() {
    std::string cigarString = "";
    for (auto op: this->CIGAR) cigarString += (std::to_string(op.first) + op.second);
    std::cerr << "Alignment{id=" << this->identifier
              << ", ref=[" << this->refIdx.first << "," << this->refIdx.second << ")"
              << ", qry=[" << this->qryIdx.first << "," << this->qryIdx.second << ")"
              << ", inverse=" << this->inverse
              << ", CIGAR=" << cigarString
              << ", alnScore=" << this->alnScore << "\n";
    return;
}


mga::chainVec mga::getAlignmentChains(alnVec& alignments) {
    // Minimap2 chaining algorithm
    chainVec chains;

    if (alignments.empty()) return chains;

    // Sort alignments by Strand and Reference Start 
    std::sort(alignments.begin(), alignments.end(), [](const Alignment& a, const Alignment& b) {
        if (a.inverse == b.inverse) return a.refIdx.second < b.refIdx.second;
        return (b.inverse);
    });

    int n = alignments.size();
    std::vector<int> parents(n, -1);
    float w_avg = 0;
    // Constants for Chaining
    const int MAX_DIST = 10000;  // Max distance to look for a chain
    const int MAX_SKIP = 50;     // Max number of predecessors to check
    
    // Gap Penalties (Minimap2)
    // Minimap2 chan gap penalty
    auto alpha = [](const Alignment& i, const Alignment& j) -> float {
        float a = std::min((i.qryIdx.second - j.qryIdx.second), (i.refIdx.second - j.refIdx.second));
        float b = i.alnScore;
        return std::min(a, b);
    };
    auto beta = [](const Alignment& i, const Alignment& j, float w_avg) -> float {
        int l = (i.qryIdx.second - j.qryIdx.second) - (i.refIdx.second - j.refIdx.second);
        if (l == 0) return 0.0f;
        return (0.01f * w_avg * std::abs(l)) + (0.5f * std::log2(std::abs(l))); 
    };
    
    // Dynamic Programming
    for (int i = 0; i < n; ++i) {
        // Initialize
        alignments[i].chainScore = alignments[i].alnScore;

        int skip_count = 0;

        // Look backwards for the best predecessor i
        for (int j = i - 1; j >= 0; --j) {
            
            int dr = alignments[i].refIdx.second - alignments[j].refIdx.second;
            int dq = std::abs(alignments[i].qryIdx.second - alignments[j].qryIdx.second);

            // If distance is too large
            if (dr > MAX_DIST || dq > MAX_DIST) break;

            // If not on the same strand
            if (alignments[j].inverse != alignments[i].inverse) continue;
            
            // If overlapping
            int d_ref = alignments[i].refIdx.first - alignments[j].refIdx.second;
            int d_qry = (!alignments[i].inverse) ? alignments[i].qryIdx.first - alignments[j].qryIdx.second : 
                                                   alignments[j].qryIdx.first - alignments[i].qryIdx.second;
            if (d_ref < 0 || d_qry < 0) continue; 
            
            // Calculate Chain Score
            int gap_len = std::min(d_ref, d_qry);
            // If gap too large
            if (gap_len > MAX_DIST) continue; 

            float current_w_avg = static_cast<float>(dr + dq) / 2.0f;

            float chain_score = alignments[j].chainScore + alpha(alignments[i], alignments[j]) - beta(alignments[i], alignments[j], current_w_avg);
            
            if (chain_score > alignments[i].chainScore) {
                alignments[i].chainScore = chain_score;
                parents[i] = j;
            }
            skip_count++;
            if (skip_count >= MAX_SKIP) break;
        }
    }

    

    while (true) {
        // Find Global Maximum Score
        int id = chains.size();
        int best_idx = -1;
        float max_score = -1e9;
        for (int k = 0; k < n; ++k) {
            if (!alignments[k].used) {
                if (alignments[k].chainScore > max_score) {
                    max_score = alignments[k].chainScore;
                    best_idx = k;
                }
            }
        }
        if (best_idx == -1) break;
        IntVec aln;
        int curr = best_idx;
        while (curr != -1) {
            alignments[curr].used = true;
            aln.push_back(curr);
            curr = parents[curr];
            if (alignments[curr].used) break;
        }
        AlnChain chain (id, max_score, aln);
        chains.push_back(chain);
    }

    // Sort chains from high score to low
    std::sort(chains.begin(), chains.end(), [](const AlnChain& a, const AlnChain& b) {
        return a.score > b.score;
    });
    return chains;
}


inline int calcDrop(int L, char type, int rNeeded, int qNeeded) {
    if (rNeeded <= 0 && qNeeded <= 0) return 0;

    bool rC = consumesRef(type);
    bool qC = consumesQry(type);
    
    int cut = 0;
    if (rC && qC) {
        cut = std::max(rNeeded, qNeeded); 
    } else if (rC && !qC) {
        cut = (qNeeded > 0) ? L : rNeeded; 
    } else if (!rC && qC) {
        cut = (rNeeded > 0) ? L : qNeeded; 
    } else {
        cut = L; 
    }
    return std::min(cut, L);
}

bool trimAlignmentExact(mga::Alignment& aln, int reqRefHead, int reqQryHead, int reqRefTail, int reqQryTail) {
    
    int headIdx = 0, headOffset = 0;
    int rDroppedHead = 0, qDroppedHead = 0;
    
    // From Head
    for (int i = 0; i < aln.CIGAR.size(); ++i) {
        int L = aln.CIGAR[i].first;
        char type = aln.CIGAR[i].second;
        
        int rNeeded = reqRefHead - rDroppedHead;
        int qNeeded = reqQryHead - qDroppedHead;

        int drop = calcDrop(L, type, rNeeded, qNeeded);
        
        if (drop == 0) { headIdx = i; headOffset = 0; break; }

        rDroppedHead += consumesRef(type) ? drop : 0;
        qDroppedHead += consumesQry(type) ? drop : 0;

        if (drop < L) { headIdx = i; headOffset = drop; break; }
        if (i == aln.CIGAR.size() - 1 && drop == L) return false;
    }

    // From Tail
    int tailIdx = aln.CIGAR.size() - 1, tailOffset = 0; 
    int rDroppedTail = 0, qDroppedTail = 0;
    
    for (int i = aln.CIGAR.size() - 1; i >= headIdx; --i) {
        int L = (i == headIdx) ? (aln.CIGAR[i].first - headOffset) : aln.CIGAR[i].first;
        char type = aln.CIGAR[i].second;
        
        int rNeeded = reqRefTail - rDroppedTail;
        int qNeeded = reqQryTail - qDroppedTail;

        int drop = calcDrop(L, type, rNeeded, qNeeded);
        
        if (drop == 0) { tailIdx = i; tailOffset = drop; break; }

        rDroppedTail += consumesRef(type) ? drop : 0;
        qDroppedTail += consumesQry(type) ? drop : 0;

        if (drop < L) { tailIdx = i; tailOffset = drop; break; }
        if (i == headIdx && drop == L) return false;
    }

    // Reconstruct CIGAR
    std::vector<std::pair<int, char>> newCIGAR;
    for (int i = headIdx; i <= tailIdx; ++i) {
        int finalL = aln.CIGAR[i].first;
        if (i == headIdx) finalL -= headOffset;
        if (i == tailIdx) finalL -= tailOffset;
        if (finalL > 0) newCIGAR.push_back({finalL, aln.CIGAR[i].second});
    }

    if (newCIGAR.empty()) return false;

    aln.CIGAR = std::move(newCIGAR);
    aln.refIdx.first += rDroppedHead;
    aln.refIdx.second -= rDroppedTail;

    if (!aln.inverse) {
        aln.qryIdx.first += qDroppedHead;
        aln.qryIdx.second -= qDroppedTail;
    } else {
        aln.qryIdx.second -= qDroppedHead;
        aln.qryIdx.first += qDroppedTail;
    }

    return true;
}


/*
void mga::identifyPrimaryAlignments(alnVec& alignments, chainVec& chains) {
    const int MAX_OVERLAP = 200;
    IntVec primaryIdx;

    for (const auto& chain : chains) {
        for (int alnId : chain.chainedAln) {
            Alignment& aln = alignments[alnId];
            bool discard = false;
            
            int rHeadOv = 0, rTailOv = 0, qHeadOv = 0, qTailOv = 0;

            for (int pIdx : primaryIdx) {
                const auto& pAln = alignments[pIdx];

                int rStart = std::max(aln.refIdx.first, pAln.refIdx.first);
                int rEnd = std::min(aln.refIdx.second, pAln.refIdx.second);
                int rOv = std::max(0, rEnd - rStart);

                int qStart = std::max(aln.qryIdx.first, pAln.qryIdx.first);
                int qEnd = std::min(aln.qryIdx.second, pAln.qryIdx.second);
                int qOv = std::max(0, qEnd - qStart);

                if (rOv > MAX_OVERLAP || qOv > MAX_OVERLAP) {
                    discard = true;
                    break;
                }

                if (rOv > 0) {
                    // pAln inside aln
                    if (pAln.refIdx.first > aln.refIdx.first && pAln.refIdx.second < aln.refIdx.second) {
                        discard = true; 
                        break; 
                    }

                    if (pAln.refIdx.first <= aln.refIdx.first) {
                        rHeadOv = std::max(rHeadOv, rOv);
                    } 
                    else {
                        rTailOv = std::max(rTailOv, rOv);
                    }
                }
                if (qOv > 0) {
                    // pAln inside aln
                    if (pAln.qryIdx.first > aln.qryIdx.first && pAln.qryIdx.second < aln.qryIdx.second) {
                        discard = true; 
                        break; 
                    }

                    int qHeadOv = 0, qTailOv = 0;
                    if (pAln.qryIdx.first <= aln.qryIdx.first) {
                        qHeadOv = std::max(qHeadOv, qOv);
                    } else {
                        qTailOv = std::max(qTailOv, qOv);
                    }
                }
            }

            if (!discard) {
                if (aln.inverse) {
                    auto temp = qHeadOv; qHeadOv = qTailOv; qTailOv = temp;
                }
                if (rHeadOv > 0 || rTailOv > 0 || qHeadOv > 0 || qTailOv > 0) {
                    if (trimAlignmentExact(aln, rHeadOv, rTailOv, qHeadOv, qTailOv)) {
                        aln.type = PRIMARY;
                        primaryIdx.push_back(alnId);
                    }
                }
            }
        }
    }

    int pCount = primaryIdx.size();
    std::cout << pCount << '\\' << (alignments.size() - pCount) << '\n';
}
*/

void mga::resolveGlobalOverlaps(alnVec& alignments) {
    
    // Collect all end points of alignments
    std::set<int> refEndpointsSet;
    std::set<int> qryEndpointsSet;
    for (const auto& aln : alignments) {
        if (!aln.valid) continue;
        refEndpointsSet.insert(aln.refIdx.first);
        refEndpointsSet.insert(aln.refIdx.second);
        qryEndpointsSet.insert(aln.qryIdx.first);
        qryEndpointsSet.insert(aln.qryIdx.second);
    }

    std::vector<int> allRefCuts(refEndpointsSet.begin(), refEndpointsSet.end());
    std::vector<int> allQryCuts(qryEndpointsSet.begin(), qryEndpointsSet.end());

    std::vector<Alignment> slicedAlignments;

    int maxAlnID = 0;
    for (const auto& aln : alignments) {
        if (aln.identifier > maxAlnID) {
            maxAlnID = aln.identifier;
        }
    }
    int nextAlnID = maxAlnID + 1; // 新 Alignment 的 ID 起始值

    // Split alignments based on end points
    for (const auto& aln : alignments) {
        if (!aln.valid) continue;

        // Find cut points in the alignment
        std::deque<int> rCuts;
        for (int c : allRefCuts) {
            if (c > aln.refIdx.first && c < aln.refIdx.second) rCuts.push_back(c);
        }

        std::deque<int> qCuts;
        if (!aln.inverse) {
            for (int c : allQryCuts) {
                if (c > aln.qryIdx.first && c < aln.qryIdx.second) qCuts.push_back(c);
            }
        } else {
            for (auto it = allQryCuts.rbegin(); it != allQryCuts.rend(); ++it) {
                if (*it > aln.qryIdx.first && *it < aln.qryIdx.second) qCuts.push_back(*it);
            }
        }
 
        Alignment currChunk = aln; 
        currChunk.CIGAR.clear();
        
        int cRef = aln.refIdx.first;
        int cQry = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first;

        currChunk.refIdx.first = cRef;
        if (!aln.inverse) currChunk.qryIdx.first = cQry;
        else currChunk.qryIdx.second = cQry; 

        for (const auto& op : aln.CIGAR) {
            int L = op.first;
            char type = op.second;
            bool rC = consumesRef(type);
            bool qC = consumesQry(type);

            while (L > 0) {
                int step = L;
                if (rC && !rCuts.empty()) {
                    step = std::min(step, rCuts.front() - cRef);
                }
                if (qC && !qCuts.empty()) {
                    int qDist = aln.inverse ? (cQry - qCuts.front()) : (qCuts.front() - cQry);
                    step = std::min(step, qDist);
                }

                if (step > 0) {
                    currChunk.CIGAR.push_back({step, type});
                    cRef += rC ? step : 0;
                    cQry += qC ? (aln.inverse ? -step : step) : 0;
                    L -= step;
                }

                bool hitCut = false;
                if (rC && !rCuts.empty() && cRef == rCuts.front()) {
                    rCuts.pop_front();
                    hitCut = true;
                }
                if (qC && !qCuts.empty() && cQry == qCuts.front()) {
                    qCuts.pop_front();
                    hitCut = true;
                }

                if (hitCut) {
                    currChunk.refIdx.second = cRef;
                    if (!aln.inverse) currChunk.qryIdx.second = cQry;
                    else currChunk.qryIdx.first = cQry; 

                    if (!currChunk.CIGAR.empty()) {
                        currChunk.identifier = nextAlnID++;
                        slicedAlignments.push_back(currChunk);
                    }

                    currChunk.CIGAR.clear();
                    currChunk.refIdx.first = cRef;
                    if (!aln.inverse) currChunk.qryIdx.first = cQry;
                    else currChunk.qryIdx.second = cQry;
                }
            }
        }
        if (!currChunk.CIGAR.empty()) {
            currChunk.refIdx.second = cRef;
            if (!aln.inverse) currChunk.qryIdx.second = cQry;
            else currChunk.qryIdx.first = cQry;
            currChunk.identifier = nextAlnID++;
            slicedAlignments.push_back(currChunk);
        }
    }

    // ==========================================
    // 步驟 3：長度過濾與 Secondary (Duplication) 標記
    // ==========================================
    std::vector<Alignment> finalAlignments;
    std::vector<int> primaryIndices; // 紀錄在 finalAlignments 中的 Index

    const int MIN_LEN = 200;

    // 先把 PRIMARY 放進去，建立基準
    for (auto& aln : slicedAlignments) {
        int rLen = aln.refIdx.second - aln.refIdx.first;
        int qLen = aln.qryIdx.second - aln.qryIdx.first;
        if (std::max(rLen, qLen) < MIN_LEN) continue; 

        if (aln.type == PRIMARY) {
            finalAlignments.push_back(aln);
            primaryIndices.push_back(finalAlignments.size() - 1);
        }
    }

    // 再處理 NON-PRIMARY
    for (auto& aln : slicedAlignments) {
        if (aln.type == PRIMARY) continue; 
        
        int rLen = aln.refIdx.second - aln.refIdx.first;
        int qLen = aln.qryIdx.second - aln.qryIdx.first;
        if (std::max(rLen, qLen) < MIN_LEN) continue;

        int refOverlapPrimaryIdx = -1;
        int qryOverlapPrimaryIdx = -1;

        // 檢查與 Primary 的重疊情形 (因為 Primary 互不重疊，所以一旦找到就可以 break)
        for (int pIdx : primaryIndices) {
            const auto& pAln = finalAlignments[pIdx];
            
            if (aln.refIdx.first == pAln.refIdx.first && aln.refIdx.second == pAln.refIdx.second) {
                refOverlapPrimaryIdx = pIdx;
            }
            if (aln.qryIdx.first == pAln.qryIdx.first && aln.qryIdx.second == pAln.qryIdx.second) {
                qryOverlapPrimaryIdx = pIdx;
            }

            // 如果兩軸都找到了對應的 Primary，就可以提早結束搜尋
            if (refOverlapPrimaryIdx != -1 && qryOverlapPrimaryIdx != -1) break; 
        }

        bool hasRefOverlap = (refOverlapPrimaryIdx != -1);
        bool hasQryOverlap = (qryOverlapPrimaryIdx != -1);

        if (hasRefOverlap && hasQryOverlap) {
            // 【狀況 A：兩軸都有重疊 (Bridge Alignment)】
            // 推斷 Paralogs：如果 Ref 撞到的 Primary 和 Qry 撞到的不是同一個
            if (refOverlapPrimaryIdx != qryOverlapPrimaryIdx) {
                int p1_id = finalAlignments[refOverlapPrimaryIdx].identifier;
                int p2_id = finalAlignments[qryOverlapPrimaryIdx].identifier;
                
                finalAlignments[refOverlapPrimaryIdx].paralogs.push_back(p2_id);
                finalAlignments[qryOverlapPrimaryIdx].paralogs.push_back(p1_id);
            }
            // 完成橋接任務，直接捨棄這條第三方 Alignment
            continue; 

        } else if (hasRefOverlap || hasQryOverlap) {
            // 【狀況 B：單一軸有重疊 (Duplication)】
            aln.type = SECONDARY;
            int targetIdx = hasRefOverlap ? refOverlapPrimaryIdx : qryOverlapPrimaryIdx;
            finalAlignments[targetIdx].duplications.push_back(aln.identifier);
            finalAlignments.push_back(aln);

        } else {
            // 【狀況 C：兩軸都沒有重疊】
            aln.type = PRIMARY;
            primaryIndices.push_back(finalAlignments.size());
            finalAlignments.push_back(aln);
        }
    }

    // 更新原陣列
    alignments = std::move(finalAlignments);
}

void mga::fillUnalignedRegions(alnVec& alignments, int refTotalLen, int qryTotalLen) {
    
    std::vector<std::pair<int, int>> refIntervals;
    std::vector<std::pair<int, int>> qryIntervals;

    int maxID = 0;
    for (const auto& aln : alignments) {
        if (aln.identifier > maxID) maxID = aln.identifier;

        if ((aln.type == mga::PRIMARY || aln.type == mga::SECONDARY) && aln.valid) {
            int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
            int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
            refIntervals.push_back({rMin, rMax});

            int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
            int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
            qryIntervals.push_back({qMin, qMax});
        }
    }

    auto findGaps = [&](std::vector<std::pair<int, int>>& intervals, int totalLen) -> std::vector<std::pair<int, int>> {
        std::vector<std::pair<int, int>> gaps;
        if (totalLen <= 0) return gaps; 
        if (intervals.empty()) {
            gaps.push_back({0, totalLen});
            return gaps;
        }
        std::sort(intervals.begin(), intervals.end());
        int currEnd = 0;
        for (const auto& interval : intervals) {
            int start = std::max(0, interval.first);
            int end = std::min(totalLen, interval.second);
            if (start > currEnd) gaps.push_back({currEnd, start});
            currEnd = std::max(currEnd, end);
        }
        if (currEnd < totalLen) gaps.push_back({currEnd, totalLen});
        return gaps;
    };

    auto rawRefGaps = findGaps(refIntervals, refTotalLen);
    auto rawQryGaps = findGaps(qryIntervals, qryTotalLen);

    // ==========================================
    // Phase: Merge Small Gaps into Preceding Alignments
    // ==========================================
    
    // 輔助函數：尋找在 targetEndPos 結束的 Alignment (Primary 優先)
    auto findPrecedingAlignment = [&](int targetEndPos, bool isRef) -> int {
        int bestAlnIndex = -1;
        int bestPriority = 999; // 1 for PRIMARY, 2 for SECONDARY

        for (size_t i = 0; i < alignments.size(); ++i) {
            const auto& aln = alignments[i];
            if (!aln.valid) continue;
            if (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY) continue;

            int alnEnd = isRef ? std::max(aln.refIdx.first, aln.refIdx.second)
                               : std::max(aln.qryIdx.first, aln.qryIdx.second);

            if (alnEnd == targetEndPos) {
                int priority = (aln.type == mga::PRIMARY) ? 1 : 2;
                if (priority < bestPriority) {
                    bestPriority = priority;
                    bestAlnIndex = i;
                }
            }
        }
        return bestAlnIndex;
    };

    std::vector<std::pair<int, int>> finalRefGaps;
    std::vector<std::pair<int, int>> finalQryGaps;
    int mergedRefCount = 0, mergedQryCount = 0;

    // 處理 Ref Gaps (Deletions)
    for (const auto& gap : rawRefGaps) {
        int len = gap.second - gap.first;
        if (len < 100 && gap.first > 0) {
            int mergeIdx = findPrecedingAlignment(gap.first, true);
            if (mergeIdx != -1) {
                auto& aln = alignments[mergeIdx];
                // 擴展座標
                if (aln.refIdx.first < aln.refIdx.second) aln.refIdx.second = gap.second;
                else aln.refIdx.first = gap.second; 
                
                // 更新 CIGAR (若最後一個字元剛好也是 'D'，則直接數字相加)
                if (!aln.CIGAR.empty() && aln.CIGAR.back().second == 'D') {
                    aln.CIGAR.back().first += len;
                } else {
                    aln.CIGAR.push_back({len, 'D'});
                }
                mergedRefCount++;
                continue; // 成功合併，跳過加入 finalRefGaps
            }
        }
        finalRefGaps.push_back(gap);
    }

    // 處理 Qry Gaps (Insertions)
    for (const auto& gap : rawQryGaps) {
        int len = gap.second - gap.first;
        if (len < 100 && gap.first > 0) {
            int mergeIdx = findPrecedingAlignment(gap.first, false);
            if (mergeIdx != -1) {
                auto& aln = alignments[mergeIdx];
                // 擴展座標
                if (aln.qryIdx.first < aln.qryIdx.second) aln.qryIdx.second = gap.second;
                else aln.qryIdx.first = gap.second;
                
                // 更新 CIGAR
                if (!aln.CIGAR.empty() && aln.CIGAR.back().second == 'I') {
                    aln.CIGAR.back().first += len;
                } else {
                    aln.CIGAR.push_back({len, 'I'});
                }
                mergedQryCount++;
                continue; 
            }
        }
        finalQryGaps.push_back(gap);
    }

    // ==========================================
    // Phase: Print & Create Remaining UNALIGNED Blocks
    // ==========================================
    
    std::cout << "\n[Merger] --- UNALIGNED Regions Detail ---\n";
    std::cout << "Merged " << mergedRefCount << " small Ref gaps and " 
              << mergedQryCount << " small Qry gaps into preceding alignments.\n\n";

    auto printGapList = [&](const std::string& label, const std::vector<std::pair<int, int>>& gaps) {
        long long totalGapLen = 0;
        std::cout << ">> " << label << " Gaps:\n";
        if (gaps.empty()) {
            std::cout << "   (None)\n";
            return 0LL;
        }
        for (size_t i = 0; i < gaps.size(); ++i) {
            int len = gaps[i].second - gaps[i].first;
            totalGapLen += len;
            std::cout << "   [" << i + 1 << "] Range: [" << gaps[i].first << ", " << gaps[i].second 
                      << "] | Length: " << len << " bp";
            if (len > 1000) std::cout << " [LARGE]";
            std::cout << "\n";
        }
        return totalGapLen;
    };

    long long refUnalignedLen = printGapList("Reference (Missing in Query)", finalRefGaps);
    long long qryUnalignedLen = printGapList("Query (Insertion/Novel)", finalQryGaps);

    // 計算 Coverage (基於合併後的 finalGaps)
    long long refCoveredLen = std::max(0LL, (long long)refTotalLen - refUnalignedLen);
    double refCovPct = refTotalLen > 0 ? (refCoveredLen * 100.0 / refTotalLen) : 0.0;

    long long qryCoveredLen = std::max(0LL, (long long)qryTotalLen - qryUnalignedLen);
    double qryCovPct = qryTotalLen > 0 ? (qryCoveredLen * 100.0 / qryTotalLen) : 0.0;

    // --- 將剩下沒有被合併的大塊 Gap 轉為 UNALIGNED ---
    for (const auto& gap : finalRefGaps) {
        int len = gap.second - gap.first;
        if (len < 1) continue; 
        mga::Alignment refSeg;
        refSeg.identifier = ++maxID;
        refSeg.refIdx = gap;
        refSeg.qryIdx = {-1, -1};
        refSeg.valid = true;
        refSeg.type = mga::UNALIGNED; 
        refSeg.CIGAR.push_back({len, 'D'}); 
        alignments.push_back(refSeg);
    }

    for (const auto& gap : finalQryGaps) {
        int len = gap.second - gap.first;
        if (len < 1) continue;
        mga::Alignment qrySeg;
        qrySeg.identifier = ++maxID;
        qrySeg.refIdx = {-1, -1};
        qrySeg.qryIdx = gap;
        qrySeg.valid = true;
        qrySeg.type = mga::UNALIGNED;
        qrySeg.CIGAR.push_back({len, 'I'}); 
        alignments.push_back(qrySeg);
    }
    
    std::cout << "\n--- Final Coverage Summary ---\n";
    std::cout << "Ref Coverage: " << std::fixed << std::setprecision(2) << refCovPct << "% (" << refCoveredLen << "/" << refTotalLen << " bp)\n";
    std::cout << "Qry Coverage: " << std::fixed << std::setprecision(2) << qryCovPct << "% (" << qryCoveredLen << "/" << qryTotalLen << " bp)\n";
    std::cout << "------------------------------\n\n";
}

bool mga::validateCoverage(const alnVec& alignments, int refTotalLen, int qryTotalLen) {
    std::cerr << "\n[DEBUG] Validating Sequence Coverage...\n";
    
    auto checkSingleSide = [&](const std::string& label, int totalLen, bool isRef) -> bool {
        if (totalLen <= 0) return true; // 防呆

        std::vector<std::pair<int, int>> intervals;
        
        for (const auto& aln : alignments) {
            if (!aln.valid) continue; 
            
            int start, end;
            if (isRef) {
                if (aln.refIdx.first == -1) continue; 
                // 修正：對 Reference 加上 min/max
                start = std::min(aln.refIdx.first, aln.refIdx.second);
                end = std::max(aln.refIdx.first, aln.refIdx.second);
            } else {
                if (aln.qryIdx.first == -1) continue; 
                start = std::min(aln.qryIdx.first, aln.qryIdx.second);
                end = std::max(aln.qryIdx.first, aln.qryIdx.second);
            }
            
            if (start < 0) start = 0;
            if (end > totalLen) end = totalLen;
            if (end > start) intervals.push_back({start, end});
        }

        std::sort(intervals.begin(), intervals.end());

        std::vector<std::pair<int, int>> gaps;
        int currCovered = 0;

        for (const auto& interval : intervals) {
            if (interval.first > currCovered) {
                gaps.push_back({currCovered, interval.first});
            }
            currCovered = std::max(currCovered, interval.second);
        }

        if (currCovered < totalLen) {
            gaps.push_back({currCovered, totalLen});
        }

        if (gaps.empty()) {
            std::cerr << "  > " << label << ": OK (Fully Covered " << totalLen << " bp)\n";
            return true;
        } else {
            std::cerr << "  > " << label << ": FAILED! Found " << gaps.size() << " gaps:\n";
            for (const auto& gap : gaps) {
                std::cerr << "    - Gap: [" << gap.first << ", " << gap.second << ") Len: " << (gap.second - gap.first) << "\n";
            }
            return false;
        }
    };

    bool refOk = checkSingleSide("Reference", refTotalLen, true);
    bool qryOk = checkSingleSide("Query    ", qryTotalLen, false);

    if (refOk && qryOk) {
        std::cerr << "[DEBUG] Validation PASSED. All sequences are fully covered.\n";
        return true;
    } else {
        std::cerr << "[DEBUG] Validation FAILED.\n";
        return false;
    }
}

// Helper to check if a variation overlaps with the first split part
bool isVarInPart1(const Variation& v, int offset) {
    return v.start < offset;
}

// Helper to check if a variation overlaps with the second split part
bool isVarInPart2(const Variation& v, int offset) {
    return v.end > offset;
}

std::pair<std::shared_ptr<Block>, std::shared_ptr<Block>> Block::split(int offset, ID new_id_1, ID new_id_2) {
    // 0. Boundary Check
    if (offset <= 0 || offset >= consensus_sequence_.length()) {
        throw std::out_of_range("Split offset out of bounds");
    }

    // 1. Split Consensus Sequence
    std::string seq_str_1 = consensus_sequence_.substr(0, offset);
    std::string seq_str_2 = consensus_sequence_.substr(offset);

    auto b1 = std::make_shared<Block>(new_id_1, seq_str_1);
    auto b2 = std::make_shared<Block>(new_id_2, seq_str_2);

    // 2. Split Sequences and Variations
    for (const auto& seq : sequences_) {
        // --- Calculate split point in RAW sequence ---
        // 我們需要知道 Consensus 的 'offset' 對應到 raw_sequence 的第幾個 base
        // Logic: 
        //   Raw Length = Offset - (Total length of GAPs before offset)
        //   SNV 不影響長度，所以不用扣除
        
        int raw_split_len = offset;
        for (const auto& var : seq.variations) {
            if (var.type == Variation::GAP) {
                // 計算這個 GAP 在 offset 之前佔了多少長度
                int gap_start = var.start;
                int gap_end = var.end;
                
                // 找出 GAP 與 [0, offset) 的交集長度
                int overlap_start = std::max(gap_start, 0);
                int overlap_end = std::min(gap_end, offset);
                
                if (overlap_end > overlap_start) {
                    raw_split_len -= (overlap_end - overlap_start);
                }
            }
        }

        // --- Create SequenceInfo for Block 1 ---
        SequenceInfo s1;
        s1.sequence_id = seq.sequence_id;
        s1.start_coordinate = seq.start_coordinate;
        s1.end_coordinate = seq.start_coordinate + raw_split_len; // Update coordinate
        
        // Filter Variations for Block 1
        for (const auto& var : seq.variations) {
            if (var.start < offset) {
                if (var.end <= offset) {
                    // Case 1: Variation fully inside Block 1
                    s1.variations.push_back(var);
                } else {
                    // Case 2: Variation spans across the cut (Must be a GAP)
                    // Cut the GAP: [var.start, offset)
                    if (var.type == Variation::GAP) {
                        s1.variations.emplace_back(var.start, offset);
                    }
                }
            }
        }
        b1->addSequence(s1);

        // --- Create SequenceInfo for Block 2 ---
        SequenceInfo s2;
        s2.sequence_id = seq.sequence_id;
        s2.start_coordinate = s1.end_coordinate; // Start where s1 ended
        s2.end_coordinate = seq.end_coordinate;

        // Filter and Shift Variations for Block 2
        for (const auto& var : seq.variations) {
            if (var.end > offset) {
                if (var.start >= offset) {
                    // Case 3: Variation fully inside Block 2
                    // Shift coordinates by -offset
                    if (var.type == Variation::SNV) {
                         s2.variations.emplace_back(var.start - offset, var.alt);
                    } else { // GAP
                         s2.variations.emplace_back(var.start - offset, var.end - offset);
                    }
                } else {
                    // Case 4: Variation spans across the cut (Must be a GAP)
                    // Cut the GAP: [0, var.end - offset) relative to Block 2
                    if (var.type == Variation::GAP) {
                        s2.variations.emplace_back(0, var.end - offset);
                    }
                }
            }
        }
        b2->addSequence(s2);
    }

    // 3. Topology Update (Graph Rewiring)
    // 這裡我們需要使用 shared_from_this() 來取得目前 Block 的 pointer

    // A. 連接 b1 與 b2 (b1 -> b2)
    b1->addNextBlock(b2);
    b2->addPrevBlock(b1);

    // B. 處理 Prev Blocks (原本指向這顆 Block 的，現在要指向 b1)
    for (auto& weak_prev : prev_blocks_) {
        if (auto prev = weak_prev.lock()) {
            // 讓 b1 指向 prev
            b1->addPrevBlock(prev);
            
            // 讓 prev 指向 b1，並移除舊的 (this)
            prev->addNextBlock(b1);
            prev->removeNextBlock(shared_from_this()); 
        }
    }

    // C. 處理 Next Blocks (原本被這顆 Block 指向的，現在要被 b2 指向)
    for (auto& weak_next : next_blocks_) {
        if (auto next = weak_next.lock()) {
            // 讓 b2 指向 next
            b2->addNextBlock(next);
            
            // 讓 next 指向 b2，並移除舊的 (this)
            next->addPrevBlock(b2);
            next->removePrevBlock(shared_from_this());
        }
    }

    return {b1, b2};
}

void mga::collectCutPoints(const std::vector<mga::Alignment>& alignments, std::set<int>& refCuts, std::set<int>& qryCuts) {
    for (const auto& aln : alignments) {
        if (!aln.valid) continue;
        if (aln.refIdx.first != -1) {
            refCuts.insert(aln.refIdx.first);
            refCuts.insert(aln.refIdx.second);
        }
        if (aln.qryIdx.first != -1) {
            qryCuts.insert(std::min(aln.qryIdx.first, aln.qryIdx.second));
            qryCuts.insert(std::max(aln.qryIdx.first, aln.qryIdx.second));
        }
    }
    std::cerr << "Total cut points: \n"
              << "Ref: " << refCuts.size() << "\n"
              << "Qry: " << qryCuts.size() << "\n";
}

void mga::identifyPrimaryAlignments(alnVec& alignments, chainVec& chains) {
    PrimaryTracker refMask, qryMask;
    std::vector<Alignment> validAlns;
    std::vector<Alignment> remainingAlns;

    // 1. 攤平所有的 Alignment 並收集
    for (const auto& chain : chains) {
        for (int oldAlnId : chain.chainedAln) {
            auto& aln = alignments[oldAlnId];
            if (aln.valid) {
                // 檢查是否兩端都是 main sequence (透過你在 consensus 產生的 "_main" 後綴)
                // 註：請確保 aln.refName 和 aln.qryName 符合你 struct 裡的變數名稱
                bool isRefMain = (aln.refName.find("_main") != std::string::npos);
                bool isQryMain = (aln.qryName.find("_main") != std::string::npos);

                if (isRefMain && isQryMain) {
                    validAlns.push_back(aln);
                } else {
                    // 如果牽涉到 remaining block，賦予新 Type 並暫存，先不參與 Greedy 切割
                    aln.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(aln);
                }
            }
        }
    }

    // 2. 依照 Alignment Score 由大到小排序 (Greedy 核心)
    std::sort(validAlns.begin(), validAlns.end(), [](const Alignment& a, const Alignment& b) {
        return a.alnScore > b.alnScore;
    });

    std::vector<Alignment> primaryList;
    std::vector<Alignment> secondaryList;
    const int MIN_ALN_LEN = 100; 
    int nextGlobalId = 1;

    // ==========================================
    // Phase 1: Greedy 填滿 Ref 與 Qry 軸
    // ==========================================
    for (const auto& aln : validAlns) {
        std::set<int> rCuts, qCuts;

        // 1-1. 找出與「已經存在的 Primary」重疊的切點
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        refMask.getCuts(rMin, rMax, rCuts);

        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
        qryMask.getCuts(qMin, qMax, qCuts);

        // 1-2 & 1-3. 依據投影點切碎 Alignment
        std::vector<Alignment> fragments;
        if (!rCuts.empty() || !qCuts.empty()) {
            fragments = splitSingleAlignment(aln, rCuts, qCuts);
        } else {
            fragments.push_back(aln);
        }

        // 檢驗每一個碎片，決定是 Primary 還是 Secondary
        for (auto& frag : fragments) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            std::set<int> refOverlapIds = refMask.getOverlappingIds(fRMin, fRMax);
            std::set<int> qryOverlapIds = qryMask.getOverlappingIds(fQMin, fQMax);

            bool refCov = !refOverlapIds.empty();
            bool qryCov = !qryOverlapIds.empty();

            frag.identifier = nextGlobalId++;
            frag.duplications.clear();
            frag.paralogs.clear();

            if (!refCov && !qryCov) {
                // 完全沒有 Overlap -> 檢查長度
                // if ((fRMax - fRMin) >= MIN_ALN_LEN) {
                    frag.type = mga::PRIMARY;
                    refMask.add(fRMin, fRMax, frag.identifier);
                    qryMask.add(fQMin, fQMax, frag.identifier);
                    primaryList.push_back(frag);
                // }
            } else {
                // 有 Overlap -> 設為 Secondary
                frag.type = mga::SECONDARY;
                secondaryList.push_back(frag);
            }
        }
    }

    // ==========================================
    // Phase 2: 動態切割 (Boundary Tolerance > 50)
    // ==========================================

    std::cerr << "Primary: " << primaryList.size() << "\t"
              << "Secondary (Raw): " << secondaryList.size() << "\n";
    
    // 【修改點】：在收集切點之前，先做 Redundant 與長度過濾！
    std::vector<Alignment> validSecondaries;
    for (auto& sec : secondaryList) {
        if (std::abs(sec.refIdx.second - sec.refIdx.first) < MIN_ALN_LEN) continue;

        int sRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int sRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int sQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int sQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        // 透過原本 Phase 1 建立的 Mask 來查表
        std::set<int> rIds = refMask.getOverlappingIds(sRMin, sRMax);
        std::set<int> qIds = qryMask.getOverlappingIds(sQMin, sQMax);

        bool isRedundant = false;
        if (!rIds.empty() && !qIds.empty()) {
            isRedundant = true;
            for (int id : rIds) {
                if (qIds.count(id) == 0) { isRedundant = false; break; }
            }
            for (int id : qIds) {
                if (rIds.count(id) == 0) { isRedundant = false; break; }
            }
        }

        // 只有非冗餘的 Paralog / Duplication 才有資格留下來參與切割
        if (!isRedundant) {
            validSecondaries.push_back(sec);
        }
    }
    secondaryList = std::move(validSecondaries);

    std::cerr << "Secondary (Filtered): " << secondaryList.size() << "\n";

    // 將留下來的 Secondary 依長度由大到小排序
    std::sort(secondaryList.begin(), secondaryList.end(), [](const Alignment& a, const Alignment& b) {
        return std::abs(a.refIdx.second - a.refIdx.first) > std::abs(b.refIdx.second - b.refIdx.first);
    });

    const int TOLERANCE = 50;
    std::set<int> dynamicRefCuts;
    std::set<int> dynamicQryCuts;

    // [Pass 2-A]: 蒐集所有需要切開的座標 (此時 secondaryList 已經是完全乾淨的了)
    for (const auto& sec : secondaryList) {
        int sRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int sRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int sQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int sQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> rIds = refMask.getOverlappingIds(sRMin, sRMax);
        for (int pId : rIds) {
            auto it = std::find_if(primaryList.begin(), primaryList.end(), [pId](const Alignment& a){ return a.identifier == pId; });
            if (it != primaryList.end()) {
                int pRMin = std::min(it->refIdx.first, it->refIdx.second);
                int pRMax = std::max(it->refIdx.first, it->refIdx.second);
                if (sRMin - pRMin > TOLERANCE) dynamicRefCuts.insert(sRMin);
                if (pRMax - sRMax > TOLERANCE) dynamicRefCuts.insert(sRMax);
            }
        }

        std::set<int> qIds = qryMask.getOverlappingIds(sQMin, sQMax);
        for (int pId : qIds) {
            auto it = std::find_if(primaryList.begin(), primaryList.end(), [pId](const Alignment& a){ return a.identifier == pId; });
            if (it != primaryList.end()) {
                int pQMin = std::min(it->qryIdx.first, it->qryIdx.second);
                int pQMax = std::max(it->qryIdx.first, it->qryIdx.second);
                if (sQMin - pQMin > TOLERANCE) dynamicQryCuts.insert(sQMin);
                if (pQMax - sQMax > TOLERANCE) dynamicQryCuts.insert(sQMax);
            }
        }
    }

    // [Pass 2-B]: 執行全局切割，並用新的 Tracker 註冊
    std::vector<Alignment> finalPrimaries;
    std::vector<Alignment> finalSecondaries;
    PrimaryTracker finalRefMask, finalQryMask; // 建立全新的地圖，因為 ID 會重洗

    for (const auto& p : primaryList) {
        auto frags = splitSingleAlignment(p, dynamicRefCuts, dynamicQryCuts);
        for (auto& frag : frags) {
            frag.identifier = nextGlobalId++; // 賦予新生 Primary 全新 ID
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);
            
            finalRefMask.add(fRMin, fRMax, frag.identifier);
            finalQryMask.add(fQMin, fQMax, frag.identifier);
            finalPrimaries.push_back(frag);
        }
    }

    for (const auto& sec : secondaryList) {
        if (sec.refIdx.second - sec.refIdx.first < TOLERANCE) continue;
        if (std::abs(sec.qryIdx.second - sec.qryIdx.first) < TOLERANCE) continue;
        // 同樣用這些切點切 Secondary，保證邊界完美對齊！
        auto frags = splitSingleAlignment(sec, dynamicRefCuts, dynamicQryCuts);
        for (auto& frag : frags) {
            frag.identifier = nextGlobalId++;
            finalSecondaries.push_back(frag);
        }
    }

    // [Pass 2-C]: 根據全新的 ID 與邊界，互相綁定關係
    std::map<int, std::vector<int>> primaryToDuplications;
    std::map<int, std::vector<int>> primaryToParalogs;

    for (auto& sec : finalSecondaries) {
        int fRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int fRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int fQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int fQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> refOverlapIds = finalRefMask.getOverlappingIds(fRMin, fRMax);
        std::set<int> qryOverlapIds = finalQryMask.getOverlappingIds(fQMin, fQMax);

        bool refCov = !refOverlapIds.empty();
        bool qryCov = !qryOverlapIds.empty();

        if (refCov && qryCov) {
            // 2-2: 兩端都 Overlap -> Paralog
            // 核心邏輯修改：讓對應到的兩塊 Primary "互設為 Paralog"
            for (int rPId : refOverlapIds) {
                for (int qPId : qryOverlapIds) {
                    if (rPId != qPId) { // 避免自己設自己為 Paralog
                        primaryToParalogs[rPId].push_back(qPId);
                        primaryToParalogs[qPId].push_back(rPId);
                    }
                }
            }
            // 依然可以把 Secondary ID 存入 Primary 的 duplication 以利後續追蹤
            for (int pId : refOverlapIds) primaryToDuplications[pId].push_back(sec.identifier);
            for (int pId : qryOverlapIds) primaryToDuplications[pId].push_back(sec.identifier);
            
        } else if (refCov && !qryCov) {
            // 2-1: 單邊 Overlap -> Duplication (Ref端)
            for (int pId : refOverlapIds) {
                sec.duplications.push_back(pId); // Secondary 記住自己攀附誰
                primaryToDuplications[pId].push_back(sec.identifier); // Primary 記住自己身上有誰
            }
        } else if (!refCov && qryCov) {
            // 2-1: 單邊 Overlap -> Duplication (Qry端)
            for (int pId : qryOverlapIds) {
                sec.duplications.push_back(pId);
                primaryToDuplications[pId].push_back(sec.identifier);
            }
        }
    }

    // [Pass 2-D]: 將收集好的關係寫回最終的 Primary 物件中
    for (auto& p : finalPrimaries) {
        int pId = p.identifier;
        if (primaryToDuplications.count(pId)) {
            std::set<int> uniqueDups(primaryToDuplications[pId].begin(), primaryToDuplications[pId].end());
            p.duplications.assign(uniqueDups.begin(), uniqueDups.end());
        }
        if (primaryToParalogs.count(pId)) {
            std::set<int> uniqueParas(primaryToParalogs[pId].begin(), primaryToParalogs[pId].end());
            p.paralogs.assign(uniqueParas.begin(), uniqueParas.end());
        }
    }

    // 3. 收尾：覆寫回原來的陣列
    alignments.clear();
    alignments.insert(alignments.end(), finalPrimaries.begin(), finalPrimaries.end());
    alignments.insert(alignments.end(), finalSecondaries.begin(), finalSecondaries.end());

    alignments.insert(alignments.end(), remainingAlns.begin(), remainingAlns.end());

    // ==========================================
    // Phase 3: Debugging Message (依照 X 軸座標排序並印出詳細資訊)
    // ==========================================
    std::vector<const Alignment*> primaryAlns;
    std::map<int, const Alignment*> idToAln; 
    int pCount = 0, sCount = 0;
    
    for (const auto& aln : alignments) {
        idToAln[aln.identifier] = &aln;
        if (aln.type == mga::PRIMARY) {
            primaryAlns.push_back(&aln);
            pCount++;
        } else {
            sCount++;
        }
    }

    std::sort(primaryAlns.begin(), primaryAlns.end(), [](const Alignment* a, const Alignment* b) {
        return std::min(a->refIdx.first, a->refIdx.second) < std::min(b->refIdx.first, b->refIdx.second);
    });

    std::cout << "\n[Merger] --- Primary Alignment Tracking (Sorted by Ref X-axis) ---\n";
    for (const auto* pAln : primaryAlns) {
        int pId = pAln->identifier;
        std::cout << "Primary ID: " << pId 
                  << " | Ref: [" << pAln->refIdx.first << ", " << pAln->refIdx.second << "]"
                  << " | Qry: [" << pAln->qryIdx.first << ", " << pAln->qryIdx.second << "]\t" << "Strand: " << pAln->inverse << '\n';
        
        // 印出所有攀附在這個 Primary 身上的 Secondary 碎塊
        std::cout << "  ├─ Duplications / Attached Fragments (Secondary):\n";
        if (!pAln->duplications.empty()) {
            for (int sId : pAln->duplications) {
                const auto* sAln = idToAln[sId]; 
                std::cout << "  │    - ID: " << sId 
                          << " | Ref: [" << sAln->refIdx.first << ", " << sAln->refIdx.second << "]"
                          << " | Qry: [" << sAln->qryIdx.first << ", " << sAln->qryIdx.second << "]\t" << "Strand: " << sAln->inverse << '\n';
            }
        } else {
            std::cout << "  │    None\n";
        }

        // 直接印出它 Paralog 連接到的其他 Primary ID
        std::cout << "  └─ Paralogs (Linked Primaries):\n";
        if (!pAln->paralogs.empty()) {
            std::cout << "       - Linked directly to Primary ID: ";
            for (size_t i = 0; i < pAln->paralogs.size(); ++i) {
                std::cout << pAln->paralogs[i] << (i + 1 == pAln->paralogs.size() ? "" : ", ");
            }
            std::cout << "\n";
        } else {
            std::cout << "       None\n";
        }
        std::cout << "\n";
    }
    
    std::cout << "[Merger] Identified " << pCount << " Primary alignments and " 
              << sCount << " Secondary exact-fragments.\n";

    std::cout << "[Greedy Merger] Output " << finalPrimaries.size() << " Primary and " 
              << finalSecondaries.size() << " Secondary alignments after dynamic cutting.\n";
}