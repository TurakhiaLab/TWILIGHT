#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <list>
#include <fstream>
#include "global_alignment.hpp"

// =====================================
// CoverageTracker
// =====================================

void mga::CoverageTracker::add(int start, int end, int alnId) {
    if (start >= end) return;
    intervals[start] = {end, alnId};
}

std::set<int> mga::CoverageTracker::getOverlappingIds(int qStart, int qEnd) const {
    std::set<int> overlappingIds;
    if (qStart >= qEnd) return overlappingIds;
    auto it = intervals.upper_bound(qStart);
    
    if (it != intervals.begin()) {
        auto prev = std::prev(it);
        if (prev->second.end > qStart) {
            overlappingIds.insert(prev->second.alnId);
        }
    }
    
    while (it != intervals.end() && it->first < qEnd) {
        overlappingIds.insert(it->second.alnId);
        ++it;
    }
    return overlappingIds;
}

void mga::CoverageTracker::getCuts(int start, int end, std::set<int>& cuts) const {
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

bool mga::CoverageTracker::isCovered(int start, int end) const {
    return !getOverlappingIds(start, end).empty();
}


// =====================================
// Alignment
// =====================================


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



void forceTrimLeft(mga::Cigar& cigar, int trim_r, int trim_q) {
    int r_discarded = 0, q_discarded = 0;
    mga::Cigar new_cigar;
    bool trimming = true;
    
    for (auto& op : cigar) {
        if (!trimming) { new_cigar.push_back(op); continue; }
        
        int r_len = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') ? op.first : 0;
        int q_len = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I') ? op.first : 0;
        
        // 整個 op 都能丟掉
        if (r_discarded + r_len <= trim_r && q_discarded + q_len <= trim_q) {
            r_discarded += r_len; q_discarded += q_len; continue;
        }
        
        // 必須精細拆解 op
        for (int i = 0; i < op.first; ++i) {
            bool consumes_r = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D');
            bool consumes_q = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I');
            
            if (trimming) {
                if (consumes_r && r_discarded < trim_r && consumes_q && q_discarded < trim_q) {
                    r_discarded++; q_discarded++;
                } else if (consumes_r && r_discarded < trim_r && !consumes_q) {
                    r_discarded++;
                } else if (!consumes_r && consumes_q && q_discarded < trim_q) {
                    q_discarded++;
                } else if (consumes_r && r_discarded < trim_r && consumes_q && q_discarded >= trim_q) {
                    r_discarded++; new_cigar.push_back({1, 'I'}); // 扣 R 不扣 Q -> Q 變為 Insertion
                } else if (consumes_q && q_discarded < trim_q && consumes_r && r_discarded >= trim_r) {
                    q_discarded++; new_cigar.push_back({1, 'D'}); // 扣 Q 不扣 R -> R 變為 Deletion
                } else {
                    trimming = false; new_cigar.push_back({1, op.second});
                }
            } else {
                new_cigar.push_back({1, op.second});
            }
        }
    }
    
    cigar.clear();
    for (auto& op : new_cigar) {
        if (!cigar.empty() && cigar.back().second == op.second) cigar.back().first += op.first;
        else cigar.push_back(op);
    }
}

void forceTrimRight(mga::Cigar& cigar, int trim_r, int trim_q) {
    int r_discarded = 0, q_discarded = 0;
    mga::Cigar new_cigar;
    bool trimming = true;
    
    for (auto it = cigar.rbegin(); it != cigar.rend(); ++it) {
        auto op = *it;
        if (!trimming) { new_cigar.push_back(op); continue; }
        
        int r_len = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') ? op.first : 0;
        int q_len = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I') ? op.first : 0;
        
        if (r_discarded + r_len <= trim_r && q_discarded + q_len <= trim_q) {
            r_discarded += r_len; q_discarded += q_len; continue;
        }
        
        for (int i = 0; i < op.first; ++i) {
            bool consumes_r = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D');
            bool consumes_q = (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I');
            
            if (trimming) {
                if (consumes_r && r_discarded < trim_r && consumes_q && q_discarded < trim_q) {
                    r_discarded++; q_discarded++;
                } else if (consumes_r && r_discarded < trim_r && !consumes_q) {
                    r_discarded++;
                } else if (!consumes_r && consumes_q && q_discarded < trim_q) {
                    q_discarded++;
                } else if (consumes_r && r_discarded < trim_r && consumes_q && q_discarded >= trim_q) {
                    r_discarded++; new_cigar.push_back({1, 'I'});
                } else if (consumes_q && q_discarded < trim_q && consumes_r && r_discarded >= trim_r) {
                    q_discarded++; new_cigar.push_back({1, 'D'});
                } else {
                    trimming = false; new_cigar.push_back({1, op.second});
                }
            } else {
                new_cigar.push_back({1, op.second});
            }
        }
    }
    
    std::reverse(new_cigar.begin(), new_cigar.end());
    cigar.clear();
    for (auto& op : new_cigar) {
        if (!cigar.empty() && cigar.back().second == op.second) cigar.back().first += op.first;
        else cigar.push_back(op);
    }
}

mga::Cigar mga::compressCigar(const mga::Cigar& cigar) {
    mga::Cigar compressed;
    for (auto& op : cigar) {
        if (!compressed.empty() && compressed.back().second == op.second) {
            compressed.back().first += op.first;
        } else if (op.first > 0) {
            compressed.push_back(op);
        }
    }
    return compressed;
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
bool isVarInPart1(Variation& v, int offset) {
    return v.getStart() < offset;
}

// Helper to check if a variation overlaps with the second split part
bool isVarInPart2(Variation& v, int offset) {
    return v.getEnd() > offset;
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

mga::alnVec phase_0_preprocessing_alignments(
    mga::alnVec& alignments, 
    BlockSet* ref_blockset, 
    BlockSet* qry_blockset, 
    std::unordered_set<int>& ref_breakpoints, 
    std::unordered_set<int>& qry_breakpoints,
    std::unordered_map<std::string, int>& ref_offsets,
    std::unordered_map<std::string, int>& qry_offsets
) {
    auto buildOffsets = [](BlockSet* bs, std::unordered_map<std::string, int>& offsets, std::unordered_set<int>& breakpoints) {
        int current_offset = 0;
        int contig_idx = 0;
        const int MIN_BACKBONE_LEN = 200;

        const auto& cache = bs->getRepresentativeBlocks();
        const auto& starts = bs->getChunkStarts(); 

        for (size_t i = 0; i < starts.size() - 1; ++i) {
            int start_idx = starts[i];
            int end_idx = starts[i+1];
            
            int chunk_len = 0;
            for (int j = start_idx; j < end_idx; ++j) {
                chunk_len += bs->getBlock(cache[j])->getConsensus().length();
            }
            
            std::string name;
            if (chunk_len >= MIN_BACKBONE_LEN) {
                name = bs->getId() + "_main_contig_" + std::to_string(contig_idx++);
            } else {
                name = bs->getId() + "_" + std::to_string(cache[start_idx]);
            }

            offsets[name] = current_offset;
            breakpoints.insert(current_offset); 
            
            current_offset += chunk_len;
        }
    };

    buildOffsets(ref_blockset, ref_offsets, ref_breakpoints);
    buildOffsets(qry_blockset, qry_offsets, qry_breakpoints);

    for (auto& aln : alignments) {
        if (!aln.valid) continue;
            
        auto ref_it = ref_offsets.find(aln.refName);
        if (ref_it != ref_offsets.end()) {
            int offset = ref_it->second;
            aln.refIdx.first += offset;
            aln.refIdx.second += offset;
        }

        auto qry_it = qry_offsets.find(aln.qryName);
        if (qry_it != qry_offsets.end()) {
            int offset = qry_it->second;
            aln.qryIdx.first += offset;
            aln.qryIdx.second += offset;
        }
    }

    std::vector<mga::Alignment> validAlns;
    for (int alnId = 0; alnId < alignments.size(); ++alnId) {
        auto& aln = alignments[alnId];
        if (aln.valid && aln.alnScore > 0) {
            validAlns.push_back(aln);
        }
    }

    std::sort(validAlns.begin(), validAlns.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        return a.alnScore > b.alnScore;
    });

    return validAlns;
}

mga::alnVec phase_1_greedy_select_primary_aln(
    mga::alnVec& validAlns, 
    mga::alnVec& primaryList, 
    mga::alnVec& secondaryList, 
    PrimaryTracker& refMask, 
    PrimaryTracker& qryMask, 
    int& nextGlobalId
) {
    const int MIN_ALN_LEN = 50; // 定義骨幹與保留區段的最低門檻

    mga::alnVec remainingAlns;

    for (const auto& aln : validAlns) {
        std::set<int> rCuts, qCuts;

        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        refMask.getCuts(rMin, rMax, rCuts);

        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
        qryMask.getCuts(qMin, qMax, qCuts);

        mga::alnVec fragments;
        if (!rCuts.empty() || !qCuts.empty()) {
            fragments = splitSingleAlignment(aln, rCuts, qCuts);
        } else {
            fragments.push_back(aln);
        }

        for (auto& frag : fragments) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            frag.identifier = nextGlobalId++;
            frag.duplications.clear();
            frag.paralogs.clear();

            if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
                frag.type = mga::REMAINING_ALN;
                remainingAlns.push_back(frag);
                continue; 
            }

            std::set<int> refOverlapIds = refMask.getOverlappingIds(fRMin, fRMax);
            std::set<int> qryOverlapIds = qryMask.getOverlappingIds(fQMin, fQMax);

            bool refCov = !refOverlapIds.empty();
            bool qryCov = !qryOverlapIds.empty();

            if (!refCov && !qryCov) {
                frag.type = mga::PRIMARY;
                refMask.add(fRMin, fRMax, frag.identifier);
                qryMask.add(fQMin, fQMax, frag.identifier);
                primaryList.push_back(frag);
            } else {
                frag.type = mga::SECONDARY;
                secondaryList.push_back(frag);
            }
        }
    }

    return remainingAlns;
}

// 注意：函數簽名增加了 ref_offsets 與 qry_offsets
mga::alnVec phase2_chain_and_fill_gaps(
    mga::alnVec& primaryList,
    const std::unordered_set<int>& ref_breakpoints_set,
    const std::unordered_set<int>& qry_breakpoints_set,
    mga::stringMap& ref_seqs,
    mga::stringMap& qry_seqs,
    const std::unordered_map<std::string, int>& ref_offsets,
    const std::unordered_map<std::string, int>& qry_offsets,
    int& nextGlobalId
) {
    bool DEBUG_MODE = false; // 開啟除錯訊息

    // ==========================================
    // 0. 輔助工具 (Lambdas)
    // ==========================================
    
    // 取得 Offset 小工具 (如果是 _main 則沒有 offset)
    auto getRefOffset = [&](const std::string& name) -> int {
        if (name.find("_main") != std::string::npos) return 0;
        auto it = ref_offsets.find(name);
        return (it != ref_offsets.end()) ? it->second : 0;
    };
    auto getQryOffset = [&](const std::string& name) -> int {
        if (name.find("_main") != std::string::npos) return 0;
        auto it = qry_offsets.find(name);
        return (it != qry_offsets.end()) ? it->second : 0;
    };

    // 安全字串擷取 (自動攔截越界，並印出防護 Log)
    auto safeSubstr = [](const std::string& seq, int local_pos, int len, const std::string& context, const std::string& seqName) -> std::string {
        if (seq.empty()) return "";
        if (local_pos < 0) {
            std::cerr << "[WARNING][" << context << "] " << seqName << " local_pos (" << local_pos << ") < 0. 修正為 0.\n";
            len += local_pos; 
            local_pos = 0;
        }
        if (len <= 0) return "";
        if (local_pos >= static_cast<int>(seq.size())) {
            std::cerr << "[WARNING][" << context << "] " << seqName << " local_pos (" << local_pos << ") >= size (" << seq.size() << "). 回傳空字串.\n";
            return "";
        }
        int actualLen = std::min(len, static_cast<int>(seq.size()) - local_pos);
        return seq.substr(local_pos, actualLen);
    };

    auto trimCigarEnds = [](mga::Cigar& cig, bool trimFront, bool trimBack) {
        if (trimBack) {
            while (!cig.empty() && (cig.back().second == 'I' || cig.back().second == 'D')) cig.pop_back();
        }
        if (trimFront) {
            int dropIdx = 0;
            while (dropIdx < cig.size() && (cig[dropIdx].second == 'I' || cig[dropIdx].second == 'D')) dropIdx++;
            if (dropIdx > 0) cig.erase(cig.begin(), cig.begin() + dropIdx);
        }
    };

    auto getConsumed = [](const mga::Cigar& cig, int& consR, int& consQ) {
        consR = 0; consQ = 0;
        for (const auto& op : cig) {
            if (op.second == 'M' || op.second == '=') { consR += op.first; consQ += op.first; }
            else if (op.second == 'D') consR += op.first;
            else if (op.second == 'I') consQ += op.first;
        }
    };

    // Breakpoint 搜尋工具
    std::vector<int> ref_bps(ref_breakpoints_set.begin(), ref_breakpoints_set.end());
    std::sort(ref_bps.begin(), ref_bps.end());
    std::vector<int> qry_bps(qry_breakpoints_set.begin(), qry_breakpoints_set.end());
    std::sort(qry_bps.begin(), qry_bps.end());

    auto has_breakpoint = [](int start, int end, const std::vector<int>& bps) {
        if (start >= end) return false;
        auto it = std::upper_bound(bps.begin(), bps.end(), start);
        return it != bps.end() && *it < end;
    };
    auto get_prev_bp = [&](int pos, const std::vector<int>& bps) -> int {
        auto it = std::upper_bound(bps.begin(), bps.end(), pos);
        return (it == bps.begin()) ? 0 : *(--it);
    };
    auto get_next_bp = [&](int pos, const std::vector<int>& bps, int maxLen) -> int {
        auto it = std::lower_bound(bps.begin(), bps.end(), pos);
        return (it == bps.end()) ? maxLen : *it;
    };

    // ==========================================
    // 1-3. 排序與建立 Backbone Chains
    // ==========================================
    std::vector<mga::Alignment> qrySorted = primaryList;
    std::sort(qrySorted.begin(), qrySorted.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        int minA = std::min(a.qryIdx.first, a.qryIdx.second);
        int minB = std::min(b.qryIdx.first, b.qryIdx.second);
        return (minA == minB) ? a.identifier < b.identifier : minA < minB;
    });

    std::unordered_map<int, int> qryRank;
    for (size_t i = 0; i < qrySorted.size(); ++i) qryRank[qrySorted[i].identifier] = i;

    std::sort(primaryList.begin(), primaryList.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        int minA = std::min(a.refIdx.first, a.refIdx.second);
        int minB = std::min(b.refIdx.first, b.refIdx.second);
        return (minA == minB) ? a.identifier < b.identifier : minA < minB;
    });

    std::vector<std::vector<mga::Alignment>> backboneChains;
    std::vector<mga::Alignment> currentChain;

    for (const auto& aln : primaryList) {
        if (currentChain.empty()) {
            currentChain.push_back(aln);
        } else {
            const auto& prev = currentChain.back();
            bool sameStrand = (aln.inverse == prev.inverse);
            bool sameSeq = (aln.refName == prev.refName && aln.qryName == prev.qryName);
            bool validNext = false;
            
            if (sameStrand && sameSeq) {
                int prevREnd = std::max(prev.refIdx.first, prev.refIdx.second);
                int currRStart = std::min(aln.refIdx.first, aln.refIdx.second);
                
                if (currRStart >= prevREnd && !has_breakpoint(prevREnd, currRStart, ref_bps)) {
                    int currQRank = qryRank[aln.identifier];
                    int prevQRank = qryRank[prev.identifier];

                    if (!aln.inverse) { 
                        int prevQEnd = std::max(prev.qryIdx.first, prev.qryIdx.second);
                        int currQStart = std::min(aln.qryIdx.first, aln.qryIdx.second);
                        if (currQStart >= prevQEnd && currQRank == prevQRank + 1 && !has_breakpoint(prevQEnd, currQStart, qry_bps)) validNext = true;
                    } else { 
                        int prevQStart = std::min(prev.qryIdx.first, prev.qryIdx.second);
                        int currQEnd = std::max(aln.qryIdx.first, aln.qryIdx.second);
                        if (currQEnd <= prevQStart && currQRank == prevQRank - 1 && !has_breakpoint(currQEnd, prevQStart, qry_bps)) validNext = true;
                    }
                }
            }
            if (validNext) currentChain.push_back(aln);
            else { backboneChains.push_back(currentChain); currentChain.clear(); currentChain.push_back(aln); }
        }
    }
    if (!currentChain.empty()) backboneChains.push_back(currentChain);

    // ==========================================
    // 4. Gap Filling 與 雙向視窗延伸
    // ==========================================
    const int MAX_GAP_SIZE = 2000;
    const float MIN_IDENTITY = 0.85f; 
    std::vector<std::vector<mga::Alignment>> filledChains;

    for (auto& chain : backboneChains) {
        if (chain.empty()) continue;

        std::vector<mga::Alignment> currentFilledChain;
        currentFilledChain.push_back(chain[0]);

        for (size_t i = 1; i < chain.size(); ++i) {
            mga::Alignment& prev = currentFilledChain.back();
            mga::Alignment curr = chain[i]; 

            int prevREnd = std::max(prev.refIdx.first, prev.refIdx.second);
            int currRStart = std::min(curr.refIdx.first, curr.refIdx.second);
            int gapR = currRStart - prevREnd;

            int gapQ = 0, prevQEnd, currQStart;
            if (!curr.inverse) {
                prevQEnd = std::max(prev.qryIdx.first, prev.qryIdx.second);
                currQStart = std::min(curr.qryIdx.first, curr.qryIdx.second);
                gapQ = currQStart - prevQEnd;
            } else {
                prevQEnd = std::min(prev.qryIdx.first, prev.qryIdx.second);
                currQStart = std::max(curr.qryIdx.first, curr.qryIdx.second);
                gapQ = prevQEnd - currQStart;
            }

            if (gapR > MAX_GAP_SIZE || gapQ > MAX_GAP_SIZE || gapR < 0 || gapQ < 0) {
                currentFilledChain.push_back(curr);
                continue;
            }

            int minGap = std::min(gapR, gapQ);
            int maxGap = std::max(gapR, gapQ);
            
            // 取得全域座標轉區域座標的 Offset
            int rOffset = getRefOffset(curr.refName);
            int qOffset = getQryOffset(curr.qryName);

            if (DEBUG_MODE) {
                std::cout << "\n------------------------------------------------------\n";
                std::cout << "[INFO] Processing Gap: Ref=" << gapR << ", Qry=" << gapQ 
                          << " (Offsets: R=" << rOffset << ", Q=" << qOffset << ")\n";
            }

            // =========================================
            // 策略 A: 一方為 0 的處理 (純 InDel)
            // =========================================
            if (minGap == 0) {
                if (maxGap < 10) {
                    if (DEBUG_MODE) std::cout << "  -> [Case A] 極小純 InDel (<10bp)。直接轉為 CIGAR 吸收至 Prev Alignment。\n";
                    if (gapR > 0) prev.CIGAR.push_back({gapR, 'D'});
                    if (gapQ > 0) prev.CIGAR.push_back({gapQ, 'I'});
                    prev.refIdx.second = currRStart;
                    if (!prev.inverse) prev.qryIdx.second = currQStart;
                    else prev.qryIdx.first = currQStart;
                    prev.CIGAR = mga::compressCigar(prev.CIGAR);
                } else {
                    if (DEBUG_MODE) std::cout << "  -> [Case A] 大型純 InDel/SV (長度 " << maxGap << "bp)。不處理，直接保留 Gap。\n";
                }
                currentFilledChain.push_back(curr);
                continue; 
            }

            // =========================================
            // 策略 B: 雙向視窗延伸 (Dual-directional Window Extension)
            // =========================================
            bool useWindow = (maxGap > minGap * 1.5) || (maxGap <= 50);

            if (useWindow && ref_seqs.count(curr.refName) && qry_seqs.count(curr.qryName)) {
                if (DEBUG_MODE) std::cout << "  -> [Case B] 差距過大或極小區塊，觸發【雙向視窗延伸】。\n";
                
                int winR = std::max(15, std::min(gapR, static_cast<int>(minGap * 2.0)));
                int winQ = std::max(15, std::min(gapQ, static_cast<int>(minGap * 2.0)));

                AlnResult resLeft, resRight;
                std::string rL, qL, rR, qR;

                // [Left Window] 
                rL = safeSubstr(ref_seqs[curr.refName], prevREnd - rOffset, winR, "LeftWin_Ref", curr.refName);
                if (!curr.inverse) {
                    qL = safeSubstr(qry_seqs[curr.qryName], prevQEnd - qOffset, winQ, "LeftWin_Qry_Fwd", curr.qryName);
                } else {
                    qL = safeSubstr(qry_seqs[curr.qryName], prevQEnd - winQ - qOffset, winQ, "LeftWin_Qry_Rev", curr.qryName);
                    qL = mga::getReverseComplement(qL);
                }
                resLeft = runSemiGlobalAlignment(rL, qL);

                // [Right Window]
                rR = safeSubstr(ref_seqs[curr.refName], currRStart - winR - rOffset, winR, "RightWin_Ref", curr.refName);
                if (!curr.inverse) {
                    qR = safeSubstr(qry_seqs[curr.qryName], currQStart - winQ - qOffset, winQ, "RightWin_Qry_Fwd", curr.qryName);
                } else {
                    qR = safeSubstr(qry_seqs[curr.qryName], currQStart - qOffset, winQ, "RightWin_Qry_Rev", curr.qryName);
                    qR = mga::getReverseComplement(qR);
                }
                resRight = runSemiGlobalAlignment(rR, qR);

                bool okLeft = resLeft.success && resLeft.identity >= MIN_IDENTITY;
                bool okRight = resRight.success && resRight.identity >= MIN_IDENTITY;

                if (okLeft && okRight) {
                    if (resLeft.score >= resRight.score) okRight = false;
                    else okLeft = false;
                }

                if (okLeft) {
                    if (DEBUG_MODE) std::cout << "    => [Result] Left Window 勝出 (Score: " << resLeft.score << ", Id: " << resLeft.identity << ")! 成功向右延伸 Prev Alignment。\n";
                    auto cig = mga::parser::parseCigar(resLeft.cigar);
                    trimCigarEnds(cig, false, true); 
                    if (!cig.empty()) {
                        int cR, cQ; getConsumed(cig, cR, cQ);
                        prev.CIGAR.insert(prev.CIGAR.end(), cig.begin(), cig.end());
                        prev.CIGAR = mga::compressCigar(prev.CIGAR);
                        prev.refIdx.second += cR;
                        if (!prev.inverse) prev.qryIdx.second += cQ;
                        else prev.qryIdx.first -= cQ;
                    }
                } else if (okRight) {
                    if (DEBUG_MODE) std::cout << "    => [Result] Right Window 勝出 (Score: " << resRight.score << ", Id: " << resRight.identity << ")! 成功向左延伸 Curr Alignment。\n";
                    auto cig = mga::parser::parseCigar(resRight.cigar);
                    trimCigarEnds(cig, true, false); 
                    if (!cig.empty()) {
                        int cR, cQ; getConsumed(cig, cR, cQ);
                        curr.CIGAR.insert(curr.CIGAR.begin(), cig.begin(), cig.end());
                        curr.CIGAR = mga::compressCigar(curr.CIGAR);
                        curr.refIdx.first -= cR;
                        if (!curr.inverse) curr.qryIdx.first -= cQ;
                        else curr.qryIdx.second += cQ;
                    }
                } else {
                    if (DEBUG_MODE) std::cout << "    => [Result] 雙向延伸皆失敗 (Identity 或 Score 不足)。放生此 Gap。\n";
                }
                currentFilledChain.push_back(curr);
                continue; 
            }

            // =========================================
            // 策略 C: 正常全 Gap 填補
            // =========================================
            int diff = std::abs(gapR - gapQ);
            int threshold = std::max(20, static_cast<int>(minGap * 0.15));
            bool isSemiGlobal = (diff > threshold);

            if (DEBUG_MODE) std::cout << "  -> [Case C] 觸發【全區塊 Gap 填補】 (使用 " << (isSemiGlobal ? "Semi-Global" : "Global") << " Alignment)。\n";

            std::string refSeq = "", qrySeq = "";
            if (ref_seqs.count(curr.refName) && qry_seqs.count(curr.qryName)) {
                refSeq = safeSubstr(ref_seqs[curr.refName], prevREnd - rOffset, gapR, "FullGap_Ref", curr.refName);
                if (!curr.inverse) {
                    qrySeq = safeSubstr(qry_seqs[curr.qryName], prevQEnd - qOffset, gapQ, "FullGap_Qry_Fwd", curr.qryName);
                } else {
                    qrySeq = safeSubstr(qry_seqs[curr.qryName], currQStart - qOffset, gapQ, "FullGap_Qry_Rev", curr.qryName);
                    qrySeq = mga::getReverseComplement(qrySeq);
                }
            }

            AlnResult result;
            if (!refSeq.empty() && !qrySeq.empty()) {
                if (!isSemiGlobal) result = runGlobalAlignment(refSeq, qrySeq);
                else result = runSemiGlobalAlignment(refSeq, qrySeq);
            } else { result.success = false; }

            if (result.success && result.identity >= MIN_IDENTITY && result.score >= 100) {
                auto cigar = mga::parser::parseCigar(result.cigar);
                bool acceptGapAln = true;
                int trimRHead = 0, trimQHead = 0, trimRTail = 0, trimQTail = 0;
                float qCoverage = 1.0f;

                if (isSemiGlobal && !cigar.empty()) {
                    int matchBases = 0;
                    for (auto& op : cigar) if (op.second == 'M' || op.second == '=') matchBases += op.first;
                    
                    qCoverage = (gapQ > 0) ? (float)matchBases / gapQ : 0;
                    if (qCoverage < 0.5f && gapQ >= 50) acceptGapAln = false; 

                    if (cigar.front().second == 'I') trimQHead = cigar.front().first;
                    if (cigar.front().second == 'D') trimRHead = cigar.front().first;
                    if (cigar.back().second == 'I') trimQTail = cigar.back().first;
                    if (cigar.back().second == 'D') trimRTail = cigar.back().first;
                }

                if (acceptGapAln && !cigar.empty()) {
                    if (DEBUG_MODE) std::cout << "    => [Result] 填補成功! (Score: " << result.score << ", Id: " << result.identity << ") 建立新的 Primary GapBlock。\n";
                    mga::Alignment gapAln;
                    gapAln.identifier = nextGlobalId++;
                    gapAln.type = mga::PRIMARY;
                    gapAln.valid = true;
                    gapAln.used = false;
                    gapAln.chainScore = 0;
                    
                    gapAln.refName = curr.refName;
                    gapAln.qryName = curr.qryName;
                    gapAln.inverse = curr.inverse;
                    gapAln.CIGAR = cigar;
                    gapAln.alnScore = result.score;
                    
                    gapAln.refIdx.first = prevREnd + trimRHead;
                    gapAln.refIdx.second = currRStart - trimRTail;
                    if (!curr.inverse) {
                        gapAln.qryIdx.first = prevQEnd + trimQHead;
                        gapAln.qryIdx.second = currQStart - trimQTail;
                    } else {
                        gapAln.qryIdx.first = currQStart + trimQTail;
                        gapAln.qryIdx.second = prevQEnd - trimQHead;
                    }
                    currentFilledChain.push_back(gapAln);
                } else {
                    if (DEBUG_MODE) std::cout << "    => [Result] 填補失敗：Coverage 過低 (" << qCoverage << ")。丟棄此 Alignment。\n";
                }
            } else {
                if (DEBUG_MODE) {
                    if (!result.success) std::cout << "    => [Result] 填補失敗：序列擷取異常。\n";
                    else std::cout << "    => [Result] 填補失敗：分數 (" << result.score << ") 或 Identity (" << result.identity << ") 未達標。\n";
                }
            }
            currentFilledChain.push_back(curr);
        }
        filledChains.push_back(currentFilledChain);
    }

    // ==========================================
    // 5. 合併連續片段
    // ==========================================
    std::vector<mga::Alignment> newPrimaryList;
    for (auto& chain : filledChains) {
        if (chain.empty()) continue;
        std::vector<mga::Alignment> currentMergedChain;
        mga::Alignment mergedAln = chain[0];
        
        for (size_t i = 1; i < chain.size(); ++i) {
            const mga::Alignment& curr = chain[i];
            bool contiguousRef = (mergedAln.refIdx.second == curr.refIdx.first);
            bool contiguousQry = (!mergedAln.inverse) ? (mergedAln.qryIdx.second == curr.qryIdx.first) : (mergedAln.qryIdx.first == curr.qryIdx.second);
            
            if (contiguousRef && contiguousQry && (mergedAln.inverse == curr.inverse)) {
                mergedAln.refIdx.second = curr.refIdx.second;
                if (!mergedAln.inverse) mergedAln.qryIdx.second = curr.qryIdx.second;
                else mergedAln.qryIdx.first = curr.qryIdx.first; 
                mergedAln.CIGAR.insert(mergedAln.CIGAR.end(), curr.CIGAR.begin(), curr.CIGAR.end());
                mergedAln.CIGAR = mga::compressCigar(mergedAln.CIGAR);
                mergedAln.alnScore += curr.alnScore;
            } else {
                currentMergedChain.push_back(mergedAln);
                mergedAln = curr;
            }
        }
        currentMergedChain.push_back(mergedAln);
        newPrimaryList.insert(newPrimaryList.end(), currentMergedChain.begin(), currentMergedChain.end());
    }

    // ==========================================
    // 5.5 斷點空間填滿 (Breakpoint Flank Filling)
    // ==========================================
    const int BP_EXTEND_LIMIT = 300; 

    for (auto& aln : newPrimaryList) {
        if (!ref_seqs.count(aln.refName) || !qry_seqs.count(aln.qryName)) continue;
        
        int refLen = ref_seqs[aln.refName].size(); 
        int qryLen = qry_seqs[aln.qryName].size();

        // 取得 Offset
        int rOffset = getRefOffset(aln.refName);
        int qOffset = getQryOffset(aln.qryName);

        // --- 尾部延伸 (Tail Flank) ---
        // 注意：這裡的 refLen 只是 Local 長度，我們加上 Offset 變成該 block 的終點全域座標
        int globalRefEnd = rOffset + refLen;
        int nextR_BP = get_next_bp(aln.refIdx.second, ref_bps, globalRefEnd);
        int distR_tail = nextR_BP - aln.refIdx.second;

        if (distR_tail > 0 && distR_tail <= BP_EXTEND_LIMIT) {
            int nextQ_BP, distQ_tail;
            int globalQryEnd = qOffset + qryLen;
            if (!aln.inverse) {
                nextQ_BP = get_next_bp(aln.qryIdx.second, qry_bps, globalQryEnd);
                distQ_tail = nextQ_BP - aln.qryIdx.second;
            } else {
                nextQ_BP = get_prev_bp(aln.qryIdx.first, qry_bps);
                distQ_tail = aln.qryIdx.first - nextQ_BP;
            }

            int winQ = std::max(15, std::min(distQ_tail, static_cast<int>(distR_tail * 2.0)));
            if (winQ > 0) {
                std::string rSeq = safeSubstr(ref_seqs[aln.refName], aln.refIdx.second - rOffset, distR_tail, "BP_Tail_Ref", aln.refName);
                std::string qSeq;
                if (!aln.inverse) {
                    qSeq = safeSubstr(qry_seqs[aln.qryName], aln.qryIdx.second - qOffset, winQ, "BP_Tail_Qry_Fwd", aln.qryName);
                } else {
                    qSeq = safeSubstr(qry_seqs[aln.qryName], aln.qryIdx.first - winQ - qOffset, winQ, "BP_Tail_Qry_Rev", aln.qryName);
                    qSeq = mga::getReverseComplement(qSeq); 
                }

                if (!rSeq.empty() && !qSeq.empty()) {
                    AlnResult res = runSemiGlobalAlignment(rSeq, qSeq);
                    if (res.success && res.identity >= MIN_IDENTITY) {
                        auto cig = mga::parser::parseCigar(res.cigar);
                        trimCigarEnds(cig, false, true); 
                        if (!cig.empty()) {
                            int cR, cQ; getConsumed(cig, cR, cQ);
                            aln.CIGAR.insert(aln.CIGAR.end(), cig.begin(), cig.end());
                            aln.refIdx.second += cR;
                            if (!aln.inverse) aln.qryIdx.second += cQ;
                            else aln.qryIdx.first -= cQ;
                        }
                    }
                }
            }
        }

        // --- 頭部延伸 (Head Flank) ---
        int prevR_BP = get_prev_bp(aln.refIdx.first, ref_bps);
        int distR_head = aln.refIdx.first - prevR_BP;

        if (distR_head > 0 && distR_head <= BP_EXTEND_LIMIT) {
            int prevQ_BP, distQ_head;
            int globalQryEnd = qOffset + qryLen;
            if (!aln.inverse) {
                prevQ_BP = get_prev_bp(aln.qryIdx.first, qry_bps);
                distQ_head = aln.qryIdx.first - prevQ_BP;
            } else {
                prevQ_BP = get_next_bp(aln.qryIdx.second, qry_bps, globalQryEnd);
                distQ_head = prevQ_BP - aln.qryIdx.second;
            }

            int winQ = std::max(15, std::min(distQ_head, static_cast<int>(distR_head * 2.0)));
            if (winQ > 0) {
                std::string rSeq = safeSubstr(ref_seqs[aln.refName], prevR_BP - rOffset, distR_head, "BP_Head_Ref", aln.refName);
                std::string qSeq;
                if (!aln.inverse) {
                    qSeq = safeSubstr(qry_seqs[aln.qryName], aln.qryIdx.first - winQ - qOffset, winQ, "BP_Head_Qry_Fwd", aln.qryName);
                } else {
                    qSeq = safeSubstr(qry_seqs[aln.qryName], aln.qryIdx.second - qOffset, winQ, "BP_Head_Qry_Rev", aln.qryName);
                    qSeq = mga::getReverseComplement(qSeq); 
                }

                if (!rSeq.empty() && !qSeq.empty()) {
                    AlnResult res = runSemiGlobalAlignment(rSeq, qSeq);
                    if (res.success && res.identity >= MIN_IDENTITY) {
                        auto cig = mga::parser::parseCigar(res.cigar);
                        trimCigarEnds(cig, true, false); 
                        if (!cig.empty()) {
                            int cR, cQ; getConsumed(cig, cR, cQ);
                            aln.CIGAR.insert(aln.CIGAR.begin(), cig.begin(), cig.end());
                            aln.refIdx.first -= cR;
                            if (!aln.inverse) aln.qryIdx.first -= cQ;
                            else aln.qryIdx.second += cQ;
                        }
                    }
                }
            }
        }
        
        aln.CIGAR = mga::compressCigar(aln.CIGAR);
    }

    return newPrimaryList;
}

/*
void phase3_split_and_snap(
    mga::alnVec& alignments,
    mga::alnVec& remainingAlns,
    mga::stringMap& ref_seqs,
    mga::stringMap& qry_seqs,
    const std::unordered_map<std::string, int>& ref_offsets,
    const std::unordered_map<std::string, int>& qry_offsets,
    std::vector<mga::Alignment>& finalPrimaries,
    std::vector<mga::Alignment>& finalSecondaries,
    int& nextGlobalId
) {
    const int MIN_ALN_LEN = 100;
    const int TOLERANCE = 50;

    auto mergeContiguousFrags = [](std::vector<mga::Alignment>& frags, int minLen) {
        if (frags.empty()) return frags;
        std::vector<mga::Alignment> cleaned;
        cleaned.push_back(frags[0]);

        for (size_t i = 1; i < frags.size(); ++i) {
            auto& prev = cleaned.back();
            auto& curr = frags[i];
            
            int prevRLen = std::abs(prev.refIdx.second - prev.refIdx.first);
            int prevQLen = std::abs(prev.qryIdx.second - prev.qryIdx.first);
            int currRLen = std::abs(curr.refIdx.second - curr.refIdx.first);
            int currQLen = std::abs(curr.qryIdx.second - curr.qryIdx.first);

            // 只要兩者其中一個是碎片，就把它們融合！
            if (prevRLen < minLen || prevQLen < minLen || currRLen < minLen || currQLen < minLen) {
                // 安全地擴展 Bounding Box，保留原本的 first/second 方向性
                if (prev.refIdx.first < prev.refIdx.second) {
                    prev.refIdx.first = std::min({prev.refIdx.first, prev.refIdx.second, curr.refIdx.first, curr.refIdx.second});
                    prev.refIdx.second = std::max({prev.refIdx.first, prev.refIdx.second, curr.refIdx.first, curr.refIdx.second});
                } else {
                    prev.refIdx.first = std::max({prev.refIdx.first, prev.refIdx.second, curr.refIdx.first, curr.refIdx.second});
                    prev.refIdx.second = std::min({prev.refIdx.first, prev.refIdx.second, curr.refIdx.first, curr.refIdx.second});
                }
                
                if (prev.qryIdx.first < prev.qryIdx.second) {
                    prev.qryIdx.first = std::min({prev.qryIdx.first, prev.qryIdx.second, curr.qryIdx.first, curr.qryIdx.second});
                    prev.qryIdx.second = std::max({prev.qryIdx.first, prev.qryIdx.second, curr.qryIdx.first, curr.qryIdx.second});
                } else {
                    prev.qryIdx.first = std::max({prev.qryIdx.first, prev.qryIdx.second, curr.qryIdx.first, curr.qryIdx.second});
                    prev.qryIdx.second = std::min({prev.qryIdx.first, prev.qryIdx.second, curr.qryIdx.first, curr.qryIdx.second});
                }

                // CIGAR 與分數合併
                prev.CIGAR.insert(prev.CIGAR.end(), curr.CIGAR.begin(), curr.CIGAR.end());
                prev.CIGAR = mga::compressCigar(prev.CIGAR);
                prev.alnScore += curr.alnScore;
            } else {
                cleaned.push_back(curr);
            }
        }
        return cleaned;
    };

    // 1. 從 alignments 中分離出 Primary 和 Secondary (並初步過濾太短的雜訊)
    mga::alnVec primaryList;
    mga::alnVec secondaryList;
    for (auto& aln : alignments) {
        if (!aln.valid) continue;
        int rLen = std::abs(aln.refIdx.second - aln.refIdx.first);
        int qLen = std::abs(aln.qryIdx.second - aln.qryIdx.first);
        
        // 過濾掉一開始就小於 100bp 的異常 Alignment
        if (rLen < MIN_ALN_LEN || qLen < MIN_ALN_LEN) {
            aln.type = mga::REMAINING_ALN;
            remainingAlns.push_back(aln);
            continue;
        }

        if (aln.type == mga::PRIMARY) primaryList.push_back(aln);
        else if (aln.type == mga::SECONDARY) secondaryList.push_back(aln);
        else remainingAlns.push_back(aln);
    }

    // 3. 收集全域斷點
    std::set<int> globalRefCuts, globalQryCuts;
    for (const auto& p : primaryList) {
        globalRefCuts.insert(std::min(p.refIdx.first, p.refIdx.second));
        globalRefCuts.insert(std::max(p.refIdx.first, p.refIdx.second));
        globalQryCuts.insert(std::min(p.qryIdx.first, p.qryIdx.second));
        globalQryCuts.insert(std::max(p.qryIdx.first, p.qryIdx.second));
    }
    for (const auto& sec : secondaryList) {
        globalRefCuts.insert(std::min(sec.refIdx.first, sec.refIdx.second));
        globalRefCuts.insert(std::max(sec.refIdx.first, sec.refIdx.second));
        globalQryCuts.insert(std::min(sec.qryIdx.first, sec.qryIdx.second));
        globalQryCuts.insert(std::max(sec.qryIdx.first, sec.qryIdx.second));
    }

    auto deduplicateCuts = [](std::set<int>& cuts, int tol) {
        if (cuts.empty()) return;
        std::set<int> cleanCuts;
        int lastCut = -1e9; 
        for (int c : cuts) {
            if (c - lastCut >= tol) { 
                cleanCuts.insert(c);
                lastCut = c;
            }
        }
        cuts = cleanCuts;
    };
    deduplicateCuts(globalRefCuts, TOLERANCE);
    deduplicateCuts(globalQryCuts, TOLERANCE);

    auto filterCutsForAlignment = [](const std::set<int>& globalCuts, int start, int end, int minLen) {
        std::set<int> validCuts;
        int lastCut = start;
        for (int c : globalCuts) {
            if (c > start && c < end) {
                if (c - lastCut >= minLen && end - c >= minLen) {
                    validCuts.insert(c);
                    lastCut = c;
                }
            }
        }
        return validCuts;
    };

    finalPrimaries.clear();
    finalSecondaries.clear();
    PrimaryTracker finalRefMask, finalQryMask;
    std::unordered_map<int, mga::Alignment> finalPrimaryMap;

    // 5. 切割 Primary Alignments (含碎片吸收)
    for (const auto& p : primaryList) {
        int pRMin = std::min(p.refIdx.first, p.refIdx.second);
        int pRMax = std::max(p.refIdx.first, p.refIdx.second);
        int pQMin = std::min(p.qryIdx.first, p.qryIdx.second);
        int pQMax = std::max(p.qryIdx.first, p.qryIdx.second);

        std::set<int> rCuts = filterCutsForAlignment(globalRefCuts, pRMin, pRMax, MIN_ALN_LEN);
        std::set<int> qCuts = filterCutsForAlignment(globalQryCuts, pQMin, pQMax, MIN_ALN_LEN);

        auto frags = splitSingleAlignment(p, rCuts, qCuts); 
        
        // 【執行吸收】：把因雙重切點產生的微小碎片無縫融合！
        frags = mergeContiguousFrags(frags, MIN_ALN_LEN);

        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            // 雙重保險：如果整條合併後還是太短，就丟棄
            if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
                frag.type = mga::REMAINING_ALN;
                remainingAlns.push_back(frag);
                continue; 
            }

            frag.identifier = nextGlobalId++;
            finalRefMask.add(fRMin, fRMax, frag.identifier);
            finalQryMask.add(fQMin, fQMax, frag.identifier);
            finalPrimaries.push_back(frag);
            finalPrimaryMap[frag.identifier] = frag; 
        }
    }

    // 6. 切割 Secondary Alignments (含碎片吸收)
    for (const auto& sec : secondaryList) {
        int sRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int sRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int sQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int sQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> rCuts = filterCutsForAlignment(globalRefCuts, sRMin, sRMax, MIN_ALN_LEN);
        std::set<int> qCuts = filterCutsForAlignment(globalQryCuts, sQMin, sQMax, MIN_ALN_LEN);

        auto frags = splitSingleAlignment(sec, rCuts, qCuts);
        
        // 【執行吸收】
        frags = mergeContiguousFrags(frags, MIN_ALN_LEN);

        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
                frag.type = mga::REMAINING_ALN;
                remainingAlns.push_back(frag);
                continue; 
            }
            
            frag.identifier = nextGlobalId++;
            std::set<int> rIds = finalRefMask.getOverlappingIds(fRMin, fRMax);
            std::set<int> qIds = finalQryMask.getOverlappingIds(fQMin, fQMax);

            if (rIds.empty() && qIds.empty()) {
                frag.type = mga::PRIMARY;
                finalRefMask.add(fRMin, fRMax, frag.identifier);
                finalQryMask.add(fQMin, fQMax, frag.identifier);
                finalPrimaries.push_back(frag);
                finalPrimaryMap[frag.identifier] = frag;
            } else {
                finalSecondaries.push_back(frag);
            }
        }
    }

    // ==========================================
    // 7. 完美一對一關係綁定 & 序列延伸對齊 (Sequence Snapping)
    // ==========================================
    std::map<int, std::vector<int>> primaryToDuplications;
    std::map<int, std::vector<int>> primaryToParalogs;

    auto getStrictOverlaps = [&](int sMin, int sMax, const std::set<int>& rawIds, bool isRef) {
        std::set<int> strictIds;
        for (int pId : rawIds) {
            auto it = finalPrimaryMap.find(pId);
            if (it != finalPrimaryMap.end()) {
                int pMin = isRef ? std::min(it->second.refIdx.first, it->second.refIdx.second) : std::min(it->second.qryIdx.first, it->second.qryIdx.second);
                int pMax = isRef ? std::max(it->second.refIdx.first, it->second.refIdx.second) : std::max(it->second.qryIdx.first, it->second.qryIdx.second);
                int oMin = std::max(sMin, pMin);
                int oMax = std::min(sMax, pMax);
                if (oMax - oMin >= MIN_ALN_LEN / 2) strictIds.insert(pId);
            }
        }
        return strictIds;
    };

    auto getRefOffset = [&](const std::string& name) -> int {
        if (name.find("_main") != std::string::npos) return 0;
        auto it = ref_offsets.find(name);
        return (it != ref_offsets.end()) ? it->second : 0;
    };
    auto getQryOffset = [&](const std::string& name) -> int {
        if (name.find("_main") != std::string::npos) return 0;
        auto it = qry_offsets.find(name);
        return (it != qry_offsets.end()) ? it->second : 0;
    };

    auto safeSubstr = [](const std::string& seq, int local_pos, int len, const std::string& context, const std::string& seqName) -> std::string {
        if (seq.empty()) return "";
        if (local_pos < 0) { len += local_pos; local_pos = 0; }
        if (len <= 0) return "";
        if (local_pos >= static_cast<int>(seq.size())) return "";
        int actualLen = std::min(len, static_cast<int>(seq.size()) - local_pos);
        return seq.substr(local_pos, actualLen);
    };

    auto padAlignmentToBoundary = [&](mga::Alignment& aln, int tgtRMin, int tgtRMax, int tgtQMin, int tgtQMax) {
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        if (rMin == tgtRMin && rMax == tgtRMax && qMin == tgtQMin && qMax == tgtQMax) return;

        int skipR = std::max(0, tgtRMin - rMin);
        int skipQ = aln.inverse ? std::max(0, qMax - tgtQMax) : std::max(0, tgtQMin - qMin);
        int targetR = tgtRMax - tgtRMin;
        int targetQ = tgtQMax - tgtQMin;

        mga::Cigar finalCigar;
        int padRFront = std::max(0, rMin - tgtRMin);
        int padQFront = aln.inverse ? std::max(0, tgtQMax - qMax) : std::max(0, qMin - tgtQMin);

        int rOffset = getRefOffset(aln.refName);
        int qOffset = getQryOffset(aln.qryName);

        // 1. Head Sequence Padding
        if (padRFront > 0 && padQFront > 0 && ref_seqs.count(aln.refName) && qry_seqs.count(aln.qryName)) {
            std::string rSeq = safeSubstr(ref_seqs[aln.refName], tgtRMin - rOffset, padRFront, "SnapHead_Ref", aln.refName);
            std::string qSeq;
            if (!aln.inverse) {
                qSeq = safeSubstr(qry_seqs[aln.qryName], tgtQMin - qOffset, padQFront, "SnapHead_Qry_Fwd", aln.qryName);
            } else {
                qSeq = safeSubstr(qry_seqs[aln.qryName], qMax - qOffset, padQFront, "SnapHead_Qry_Rev", aln.qryName);
                qSeq = mga::getReverseComplement(qSeq);
            }
            AlnResult res = runSemiGlobalAlignment(rSeq, qSeq);
            if (res.success && res.identity >= 0.65f) { 
                auto cig = mga::parser::parseCigar(res.cigar);
                finalCigar.insert(finalCigar.end(), cig.begin(), cig.end());
            } else {
                finalCigar.push_back({padRFront, 'D'}); finalCigar.push_back({padQFront, 'I'});
            }
        } else {
            if (padRFront > 0) finalCigar.push_back({padRFront, 'D'});
            if (padQFront > 0) finalCigar.push_back({padQFront, 'I'});
        }

        // 2. 1-bp Stepping Trimming
        int curR = 0, curQ = 0, totalR = padRFront, totalQ = padQFront; 
        auto addOp = [&](char t) {
            if (!finalCigar.empty() && finalCigar.back().second == t) finalCigar.back().first++;
            else finalCigar.push_back({1, t});
        };

        for (auto op : aln.CIGAR) {
            int l = op.first; char t = op.second;
            if (t == 'S' || t == 'H') t = 'I';
            bool rCons = (t == 'M' || t == '=' || t == 'X' || t == 'D');
            bool qCons = (t == 'M' || t == '=' || t == 'X' || t == 'I');

            for (int i = 0; i < l; ++i) {
                bool useR = false, useQ = false;
                if (rCons) { if (curR < skipR) curR++; else if (totalR < targetR) { useR = true; totalR++; } }
                if (qCons) { if (curQ < skipQ) curQ++; else if (totalQ < targetQ) { useQ = true; totalQ++; } }

                if (useR && useQ) addOp((t == 'M' || t == '=' || t == 'X') ? t : 'M'); 
                else if (useR) addOp('D'); 
                else if (useQ) addOp('I'); 
            }
        }

        // 3. Tail Sequence Padding
        int padRTail = targetR - totalR;
        int padQTail = targetQ - totalQ;

        if (padRTail > 0 && padQTail > 0 && ref_seqs.count(aln.refName) && qry_seqs.count(aln.qryName)) {
            std::string rSeq = safeSubstr(ref_seqs[aln.refName], rMax - rOffset, padRTail, "SnapTail_Ref", aln.refName);
            std::string qSeq;
            if (!aln.inverse) {
                qSeq = safeSubstr(qry_seqs[aln.qryName], qMax - qOffset, padQTail, "SnapTail_Qry_Fwd", aln.qryName);
            } else {
                qSeq = safeSubstr(qry_seqs[aln.qryName], tgtQMin - qOffset, padQTail, "SnapTail_Qry_Rev", aln.qryName);
                qSeq = mga::getReverseComplement(qSeq);
            }
            AlnResult res = runSemiGlobalAlignment(rSeq, qSeq);
            if (res.success && res.identity >= 0.65f) {
                auto cig = mga::parser::parseCigar(res.cigar);
                finalCigar.insert(finalCigar.end(), cig.begin(), cig.end());
                totalR += padRTail; totalQ += padQTail; 
            }
        }
        
        // 4. Fallback Tail Gaps
        if (totalR < targetR) finalCigar.push_back({targetR - totalR, 'D'});
        if (totalQ < targetQ) finalCigar.push_back({targetQ - totalQ, 'I'});

        aln.CIGAR = mga::compressCigar(finalCigar);
        aln.refIdx.first = aln.refIdx.first < aln.refIdx.second ? tgtRMin : tgtRMax;
        aln.refIdx.second = aln.refIdx.first < aln.refIdx.second ? tgtRMax : tgtRMin;
        aln.qryIdx.first = aln.qryIdx.first < aln.qryIdx.second ? tgtQMin : tgtQMax;
        aln.qryIdx.second = aln.qryIdx.first < aln.qryIdx.second ? tgtQMax : tgtQMin;
    };

    // 8. 處理 Final Secondaries
    for (auto& sec : finalSecondaries) {
        int fRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int fRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int fQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int fQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> rRawIds = finalRefMask.getOverlappingIds(fRMin, fRMax);
        std::set<int> qRawIds = finalQryMask.getOverlappingIds(fQMin, fQMax);

        std::set<int> rIds = getStrictOverlaps(fRMin, fRMax, rRawIds, true);
        std::set<int> qIds = getStrictOverlaps(fQMin, fQMax, qRawIds, false);

        int targetRMin = fRMin, targetRMax = fRMax;
        int targetQMin = fQMin, targetQMax = fQMax;

        for (int pId : rIds) {
            auto& pAln = finalPrimaryMap[pId];
            int pRMin = std::min(pAln.refIdx.first, pAln.refIdx.second);
            int pRMax = std::max(pAln.refIdx.first, pAln.refIdx.second);
            if (std::abs(fRMin - pRMin) <= TOLERANCE) targetRMin = pRMin;
            if (std::abs(fRMax - pRMax) <= TOLERANCE) targetRMax = pRMax;
        }

        for (int pId : qIds) {
            auto& pAln = finalPrimaryMap[pId];
            int pQMin = std::min(pAln.qryIdx.first, pAln.qryIdx.second);
            int pQMax = std::max(pAln.qryIdx.first, pAln.qryIdx.second);
            if (std::abs(fQMin - pQMin) <= TOLERANCE) targetQMin = pQMin;
            if (std::abs(fQMax - pQMax) <= TOLERANCE) targetQMax = pQMax;
        }

        padAlignmentToBoundary(sec, targetRMin, targetRMax, targetQMin, targetQMax);

        bool refCov = !rIds.empty();
        bool qryCov = !qIds.empty();

        if (refCov && qryCov) {
            for (int rPId : rIds) {
                for (int qPId : qIds) {
                    if (rPId != qPId) { 
                        primaryToParalogs[rPId].push_back(qPId);
                        primaryToParalogs[qPId].push_back(rPId);
                    }
                    sec.paralogs.push_back(rPId);
                    sec.paralogs.push_back(qPId);
                }
            }
            for (int pId : rIds) primaryToDuplications[pId].push_back(sec.identifier);
            for (int pId : qIds) primaryToDuplications[pId].push_back(sec.identifier);
            
        } else if (refCov && !qryCov) {
            for (int pId : rIds) {
                sec.duplications.push_back(pId); 
                primaryToDuplications[pId].push_back(sec.identifier);
            }
        } else if (!refCov && qryCov) {
            for (int pId : qIds) {
                sec.duplications.push_back(pId);
                primaryToDuplications[pId].push_back(sec.identifier);
            }
        }
    }

    // 9. 關係寫回 Primary
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
}
*/

void phase3_split_and_snap(
    mga::alnVec& alignments,
    mga::alnVec& remainingAlns,
    mga::stringMap& ref_seqs,
    mga::stringMap& qry_seqs,
    const std::unordered_map<std::string, int>& ref_offsets,
    const std::unordered_map<std::string, int>& qry_offsets,
    std::vector<mga::Alignment>& finalPrimaries,
    std::vector<mga::Alignment>& finalSecondaries,
    int& nextGlobalId
) {
    const int MIN_INITIAL_LEN = 100; 
    const int MIN_FRAG_LEN = 15;     

    mga::alnVec primaryList, secondaryList;
    std::set<int> uniqueRefCuts, uniqueQryCuts; // 暫存用以去重

    // 1. 分離與初步收集切點
    for (auto& aln : alignments) {
        if (!aln.valid) continue;
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
        
        if ((rMax - rMin) < MIN_INITIAL_LEN || (qMax - qMin) < MIN_INITIAL_LEN) {
            aln.type = mga::REMAINING_ALN;
            remainingAlns.push_back(aln);
            continue;
        }

        if (aln.type == mga::PRIMARY) primaryList.push_back(aln);
        else if (aln.type == mga::SECONDARY) secondaryList.push_back(aln);
        else remainingAlns.push_back(aln);

        // 收集所有絕對端點
        uniqueRefCuts.insert(rMin); uniqueRefCuts.insert(rMax);
        uniqueQryCuts.insert(qMin); uniqueQryCuts.insert(qMax);
    }

    // 【優化核心 1】：將 Set 轉為連續記憶體的 Vector，準備進行極速 Binary Search
    std::vector<int> globalRefCuts(uniqueRefCuts.begin(), uniqueRefCuts.end());
    std::vector<int> globalQryCuts(uniqueQryCuts.begin(), uniqueQryCuts.end());

    // 擷取局部切點的 Lambda 函數 (O(log N) 搜尋 + O(K) 複製)
    auto getLocalCuts = [](int minPos, int maxPos, const std::vector<int>& globalCuts) {
        std::set<int> localCuts;
        auto it_start = std::lower_bound(globalCuts.begin(), globalCuts.end(), minPos + 1);
        auto it_end = std::lower_bound(globalCuts.begin(), globalCuts.end(), maxPos);
        for (auto it = it_start; it != it_end; ++it) {
            localCuts.insert(*it);
        }
        return localCuts;
    };

    finalPrimaries.clear();
    finalSecondaries.clear();
    PrimaryTracker finalRefMask, finalQryMask;
    std::unordered_map<int, mga::Alignment> finalPrimaryMap;

    // ==========================================
    // 3. 嚴格切割 Primary Alignments
    // ==========================================
    for (const auto& p : primaryList) {
        int rMin = std::min(p.refIdx.first, p.refIdx.second);
        int rMax = std::max(p.refIdx.first, p.refIdx.second);
        int qMin = std::min(p.qryIdx.first, p.qryIdx.second);
        int qMax = std::max(p.qryIdx.first, p.qryIdx.second);

        // 【優化核心 2】：只傳遞落在該區間內的切點，不把整個世界的切點傳進去！
        std::set<int> localRCuts = getLocalCuts(rMin, rMax, globalRefCuts);
        std::set<int> localQCuts = getLocalCuts(qMin, qMax, globalQryCuts);

        auto frags = splitSingleAlignment(p, localRCuts, localQCuts); 

        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            if ((fRMax - fRMin) < MIN_FRAG_LEN || (fQMax - fQMin) < MIN_FRAG_LEN) {
                frag.type = mga::REMAINING_ALN;
                remainingAlns.push_back(frag);
                continue; 
            }

            frag.identifier = nextGlobalId++;
            finalRefMask.add(fRMin, fRMax, frag.identifier);
            finalQryMask.add(fQMin, fQMax, frag.identifier);
            finalPrimaries.push_back(frag);
            finalPrimaryMap[frag.identifier] = frag; 
        }
    }

    // ==========================================
    // 4. 嚴格切割 Secondary Alignments
    // ==========================================
    for (const auto& sec : secondaryList) {
        int rMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int rMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int qMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int qMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> localRCuts = getLocalCuts(rMin, rMax, globalRefCuts);
        std::set<int> localQCuts = getLocalCuts(qMin, qMax, globalQryCuts);

        auto frags = splitSingleAlignment(sec, localRCuts, localQCuts);

        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            if ((fRMax - fRMin) < MIN_FRAG_LEN || (fQMax - fQMin) < MIN_FRAG_LEN) {
                frag.type = mga::REMAINING_ALN;
                remainingAlns.push_back(frag);
                continue; 
            }
            
            frag.identifier = nextGlobalId++;
            std::set<int> rIds = finalRefMask.getOverlappingIds(fRMin, fRMax);
            std::set<int> qIds = finalQryMask.getOverlappingIds(fQMin, fQMax);

            if (rIds.empty() && qIds.empty()) {
                frag.type = mga::PRIMARY;
                finalRefMask.add(fRMin, fRMax, frag.identifier);
                finalQryMask.add(fQMin, fQMax, frag.identifier);
                finalPrimaries.push_back(frag);
                finalPrimaryMap[frag.identifier] = frag;
            } else {
                finalSecondaries.push_back(frag);
            }
        }
    }
    // ==========================================
    // 5. 綁定 Paralog / Duplication 關係
    // ==========================================
    std::map<int, std::vector<int>> primaryToDuplications;
    std::map<int, std::vector<int>> primaryToParalogs;

    // 尋找高重疊率 (嚴格覆蓋) 的關聯 ID
    auto getStrictOverlaps = [&](int sMin, int sMax, const std::set<int>& rawIds, bool isRef) {
        std::set<int> strictIds;
        for (int pId : rawIds) {
            auto it = finalPrimaryMap.find(pId);
            if (it != finalPrimaryMap.end()) {
                int pMin = isRef ? std::min(it->second.refIdx.first, it->second.refIdx.second) : std::min(it->second.qryIdx.first, it->second.qryIdx.second);
                int pMax = isRef ? std::max(it->second.refIdx.first, it->second.refIdx.second) : std::max(it->second.qryIdx.first, it->second.qryIdx.second);
                
                // 因為我們是用嚴格全域切點切割的，所以如果它們有交集，邊界通常是完全一致的。
                // 這裡要求至少有 80% 的重疊率才算綁定，防止微小交錯。
                int oMin = std::max(sMin, pMin);
                int oMax = std::min(sMax, pMax);
                int overlapLen = oMax - oMin;
                int fragLen = sMax - sMin;
                if (overlapLen > 0 && (double)overlapLen / fragLen >= 0.8) {
                    strictIds.insert(pId);
                }
            }
        }
        return strictIds;
    };

    for (auto& sec : finalSecondaries) {
        int fRMin = std::min(sec.refIdx.first, sec.refIdx.second);
        int fRMax = std::max(sec.refIdx.first, sec.refIdx.second);
        int fQMin = std::min(sec.qryIdx.first, sec.qryIdx.second);
        int fQMax = std::max(sec.qryIdx.first, sec.qryIdx.second);

        std::set<int> rRawIds = finalRefMask.getOverlappingIds(fRMin, fRMax);
        std::set<int> qRawIds = finalQryMask.getOverlappingIds(fQMin, fQMax);

        std::set<int> rIds = getStrictOverlaps(fRMin, fRMax, rRawIds, true);
        std::set<int> qIds = getStrictOverlaps(fQMin, fQMax, qRawIds, false);

        // 注意：這裡完全不需要呼叫 padAlignmentToBoundary！
        // CIGAR 會保持 minimap2 最初的真實面貌，不會無端產生 Gap。

        bool refCov = !rIds.empty();
        bool qryCov = !qIds.empty();

        if (refCov && qryCov) {
            for (int rPId : rIds) {
                for (int qPId : qIds) {
                    if (rPId != qPId) { 
                        primaryToParalogs[rPId].push_back(qPId);
                        primaryToParalogs[qPId].push_back(rPId);
                    }
                    sec.paralogs.push_back(rPId);
                    sec.paralogs.push_back(qPId);
                }
            }
            for (int pId : rIds) primaryToDuplications[pId].push_back(sec.identifier);
            for (int pId : qIds) primaryToDuplications[pId].push_back(sec.identifier);
            
        } else if (refCov && !qryCov) {
            for (int pId : rIds) {
                sec.duplications.push_back(pId); 
                primaryToDuplications[pId].push_back(sec.identifier);
            }
        } else if (!refCov && qryCov) {
            for (int pId : qIds) {
                sec.duplications.push_back(pId);
                primaryToDuplications[pId].push_back(sec.identifier);
            }
        }
    }

    // 6. 關係寫回 Primary
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

    // 7. 更新傳入的 alignments 容器
    alignments.clear();
    alignments.insert(alignments.end(), finalPrimaries.begin(), finalPrimaries.end());
    alignments.insert(alignments.end(), finalSecondaries.begin(), finalSecondaries.end());
}

void debug_check_fragmentation(
    const mga::alnVec& alignments,
    const std::unordered_set<int>& ref_breakpoints_set,
    const std::unordered_set<int>& qry_breakpoints_set,
    int min_len_threshold = 100
) {
    std::cout << "\n=======================================================\n";
    std::cout << "=== [DEBUG] Alignment Fragmentation Sanity Check ===\n";
    std::cout << "=======================================================\n";

    std::unordered_map<std::string, std::set<int>> seqCuts;

    std::unordered_set<std::string> seenRefs;
    std::unordered_set<std::string> seenQrys;
    
    for (const auto& aln : alignments) {
        if (!aln.valid) continue;
        seenRefs.insert(aln.refName);
        seenQrys.insert(aln.qryName);

        // 收集 Alignment 的端點
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        seqCuts[aln.refName].insert(rMin);
        seqCuts[aln.refName].insert(rMax);
        seqCuts[aln.qryName].insert(qMin);
        seqCuts[aln.qryName].insert(qMax);
        
        // 同時檢查單一 Alignment 自身的長度
        int rLen = rMax - rMin;
        int qLen = qMax - qMin;
        if (rLen > 0 && rLen < min_len_threshold) {
            std::cout << "[WARN] 單體過短 (Ref): Alignment ID " << aln.identifier 
                      << " on " << aln.refName << " 長度為 " << rLen << " (" << rMin << "-" << rMax << ")\n";
        }
        if (qLen > 0 && qLen < min_len_threshold) {
            std::cout << "[WARN] 單體過短 (Qry): Alignment ID " << aln.identifier 
                      << " on " << aln.qryName << " 長度為 " << qLen << " (" << qMin << "-" << qMax << ")\n";
        }
    }

    // 將 Breakpoints 注入對應的 Sequence
    for (const std::string& refName : seenRefs) {
        for (int bp : ref_breakpoints_set) seqCuts[refName].insert(bp);
    }
    for (const std::string& qryName : seenQrys) {
        for (int bp : qry_breakpoints_set) seqCuts[qryName].insert(bp);
    }

    // 3. 掃描每個 Sequence 上的端點間距
    int totalIssues = 0;

    for (const auto& kv : seqCuts) {
        const std::string& seqName = kv.first;
        const std::set<int>& cuts = kv.second;

        if (cuts.size() < 2) continue; // 只有一個端點，無法構成區間

        auto it = cuts.begin();
        int prevCut = *it;
        ++it;

        for (; it != cuts.end(); ++it) {
            int currCut = *it;
            int dist = currCut - prevCut;

            // 如果距離大於 0 (不是同一個點) 且小於閾值
            if (dist > 0 && dist < min_len_threshold) {
                std::cout << "[ERROR] 碎片化區間發現! Sequence: " << seqName 
                          << " | 區間: [" << prevCut << ", " << currCut << "] | 距離: " << dist << " bp\n";
                totalIssues++;
            }
            prevCut = currCut;
        }
    }

    // 4. 總結報告
    std::cout << "-------------------------------------------------------\n";
    if (totalIssues == 0) {
        std::cout << "[PASS] 完美! 沒有發現任何小於 " << min_len_threshold << " bp 的碎片區間。\n";
    } else {
        std::cout << "[FAIL] 警告: 共發現 " << totalIssues << " 個過短的區間。\n";
        std::cout << "請檢查 phase3 的 padding 或 trimming 邏輯是否有漏網之魚。\n";
    }
    std::cout << "=======================================================\n\n";
}



void mga::identifyPrimaryAlignments(alnVec& alignments, BlockSet* ref_blockset, BlockSet* qry_blockset, stringMap& ref_seqs, stringMap& qry_seqs) {
    bool DEBUG_MODE = false; 
    
    PrimaryTracker refMask, qryMask;
    std::vector<Alignment> validAlns;
    std::unordered_set<int> ref_breakpoints, qry_breakpoints;
    std::unordered_map<std::string, int> ref_offsets, qry_offsets;

    // ==========================================
    // Phase 0: Pre-processing alignments
    // ==========================================
    validAlns = phase_0_preprocessing_alignments(
        alignments, 
        ref_blockset, 
        qry_blockset, 
        ref_breakpoints, 
        qry_breakpoints,
        ref_offsets,
        qry_offsets
    );

    // ==========================================
    // Phase 1: Greedy Select Primary Alignments
    // ==========================================
    int nextGlobalId = 1;
    alnVec primaryList, secondaryList;
    alnVec remainingAlns = phase_1_greedy_select_primary_aln(
        validAlns, 
        primaryList, 
        secondaryList, 
        refMask, 
        qryMask,
        nextGlobalId
    );
    
    // ==========================================
    // Phase 2: Chain and fill gaps
    // ==========================================
    // primaryList = phase2_chain_and_fill_gaps(
    //     primaryList,
    //     ref_breakpoints,
    //     qry_breakpoints,
    //     ref_seqs,
    //     qry_seqs,
    //     ref_offsets,
    //     qry_offsets,
    //     nextGlobalId
    // );

    alnVec finalPrimaries, finalSecondaries;

    

    // ==========================================
    // 【修復核心 1】：把 Phase 1 & 2 的結果包裝回 alignments，供 Phase 3 讀取
    // ==========================================
    alignments.clear();
    alignments.insert(alignments.end(), primaryList.begin(), primaryList.end());
    alignments.insert(alignments.end(), secondaryList.begin(), secondaryList.end());

    // ==========================================
    // Phase 3: Collect Global break points
    // ==========================================
    phase3_split_and_snap(
        alignments,
        remainingAlns,
        ref_seqs,
        qry_seqs,
        ref_offsets,
        qry_offsets,
        finalPrimaries,
        finalSecondaries,
        nextGlobalId
    );

    

    // ==========================================
    // 【修復核心 2】：將 Phase 3 切好的終極結果，寫回 alignments
    // 確保接下來的 Debug Print 能抓到最新的狀態
    // ==========================================
    alignments.clear();
    alignments.insert(alignments.end(), finalPrimaries.begin(), finalPrimaries.end());
    alignments.insert(alignments.end(), finalSecondaries.begin(), finalSecondaries.end());
    alignments.insert(alignments.end(), remainingAlns.begin(), remainingAlns.end());


    // std::unordered_set<int> ref_bp, qry_bp;
    // int ro = 0, qo = 0;
    // for (auto id: ref_blockset->getRepresentativeBlocks()) {
    //     ref_bp.insert(ro);
    //     ro += ref_blockset->getBlock(id)->getConsensus().size();
    // }
    // for (auto id: qry_blockset->getRepresentativeBlocks()) {
    //     qry_bp.insert(qo);
    //     qo += qry_blockset->getBlock(id)->getConsensus().size();
    // }
    // debug_check_fragmentation(
    //     alignments,
    //     ref_bp,
    //     qry_bp,
    //     100
    // );
    


    // ==========================================
    // Phase 3 & 4: Coverage Calculation & Debugging
    // ==========================================
    std::vector<std::pair<int, int>> refPrimInts, qryPrimInts;
    std::vector<std::pair<int, int>> refTotalInts, qryTotalInts;
    int maxRefCoord = 0, maxQryCoord = 0;

    for (const auto& aln : finalPrimaries) {
        auto rPair = std::make_pair(std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second));
        auto qPair = std::make_pair(std::min(aln.qryIdx.first, aln.qryIdx.second), std::max(aln.qryIdx.first, aln.qryIdx.second));
        refPrimInts.push_back(rPair);
        qryPrimInts.push_back(qPair);
        refTotalInts.push_back(rPair);
        qryTotalInts.push_back(qPair);
    }
    for (const auto& aln : finalSecondaries) {
        refTotalInts.push_back({std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second)});
        qryTotalInts.push_back({std::min(aln.qryIdx.first, aln.qryIdx.second), std::max(aln.qryIdx.first, aln.qryIdx.second)});
    }

    auto calculateMergedCoverage = [](std::vector<std::pair<int, int>>& intervals, int& maxCoord) -> int {
        if (intervals.empty()) return 0;
        std::sort(intervals.begin(), intervals.end());
        int totalCovered = 0;
        int currentStart = intervals[0].first;
        int currentEnd = intervals[0].second;
        maxCoord = std::max(maxCoord, currentEnd);

        for (size_t i = 1; i < intervals.size(); ++i) {
            maxCoord = std::max(maxCoord, intervals[i].second);
            if (intervals[i].first <= currentEnd) {
                currentEnd = std::max(currentEnd, intervals[i].second); 
            } else {
                totalCovered += (currentEnd - currentStart);            
                currentStart = intervals[i].first;
                currentEnd = intervals[i].second;
            }
        }
        totalCovered += (currentEnd - currentStart);
        return totalCovered;
    };

    int dummyMax = 0;
    int primRefCov = calculateMergedCoverage(refPrimInts, dummyMax);
    int primQryCov = calculateMergedCoverage(qryPrimInts, dummyMax);
    int totalRefCov = calculateMergedCoverage(refTotalInts, maxRefCoord);
    int totalQryCov = calculateMergedCoverage(qryTotalInts, maxQryCoord);
    
    double primRefRatio = (maxRefCoord > 0) ? (primRefCov * 100.0 / maxRefCoord) : 0.0;
    double primQryRatio = (maxQryCoord > 0) ? (primQryCov * 100.0 / maxQryCoord) : 0.0;
    double totalRefRatio = (maxRefCoord > 0) ? (totalRefCov * 100.0 / maxRefCoord) : 0.0;
    double totalQryRatio = (maxQryCoord > 0) ? (totalQryCov * 100.0 / maxQryCoord) : 0.0;

    if (true) {
        std::vector<const Alignment*> primaryAlns;
        std::vector<const Alignment*> secondaryAlns;
        
        for (const auto& aln : alignments) {
            if (aln.type == mga::PRIMARY) primaryAlns.push_back(&aln);
            else if (aln.type == mga::SECONDARY) secondaryAlns.push_back(&aln);
        }

        auto sortByRef = [](const Alignment* a, const Alignment* b) {
            return std::min(a->refIdx.first, a->refIdx.second) < std::min(b->refIdx.first, b->refIdx.second);
        };
        std::sort(primaryAlns.begin(), primaryAlns.end(), sortByRef);
        std::sort(secondaryAlns.begin(), secondaryAlns.end(), sortByRef);

        if (DEBUG_MODE) {
            std::cout << "\n[Merger] --- Final Alignment Tracking ---\n";
            std::cout << "\n>>> PRIMARY ALIGNMENTS <<<\n";
            for (const auto* pAln : primaryAlns) {
                int pId = pAln->identifier;
                int refLen = std::abs(pAln->refIdx.second - pAln->refIdx.first);
                std::cout << "Primary ID: " << pId 
                          << " | Ref: [" << pAln->refIdx.first << ", " << pAln->refIdx.second << "] (Len: " << refLen << " bp)"
                          << " | Qry: [" << pAln->qryIdx.first << ", " << pAln->qryIdx.second << "]\tStrand: " << pAln->inverse << '\n';
            }

            std::cout << "\n>>> SECONDARY ALIGNMENTS <<<\n";
            for (const auto* sAln : secondaryAlns) {
                int sId = sAln->identifier;
                int refLen = std::abs(sAln->refIdx.second - sAln->refIdx.first);
                std::cout << "Secondary ID: " << sId 
                          << " | Ref: [" << sAln->refIdx.first << ", " << sAln->refIdx.second << "] (Len: " << refLen << " bp)"
                          << " | Qry: [" << sAln->qryIdx.first << ", " << sAln->qryIdx.second << "]\tStrand: " << sAln->inverse << '\n';

                if (!sAln->paralogs.empty()) {
                    std::set<int> uniqueLinks(sAln->paralogs.begin(), sAln->paralogs.end());
                    std::cout << "    └─ Type: Paralog (Links Primary IDs: ";
                    for (auto it = uniqueLinks.begin(); it != uniqueLinks.end(); ++it) std::cout << *it << (std::next(it) == uniqueLinks.end() ? "" : ", ");
                    std::cout << ")\n";
                } else if (!sAln->duplications.empty()) {
                    std::set<int> uniqueDups(sAln->duplications.begin(), sAln->duplications.end());
                    std::cout << "    └─ Type: Duplication (Attached to Primary IDs: ";
                    for (auto it = uniqueDups.begin(); it != uniqueDups.end(); ++it) std::cout << *it << (std::next(it) == uniqueDups.end() ? "" : ", ");
                    std::cout << ")\n";
                } else {
                    std::cout << "    └─ Type: Unknown / Orphaned\n";
                }
            }
        }
        std::cout << "\n[Merger] Identified " << primaryAlns.size() << " Primary and " 
                  << secondaryAlns.size() << " Secondary alignments.\n";

        std::cout << "\n[Coverage Info]\n"
                  << "--- Primary Coverage (Only Primary Regions) ---\n"
                  << "  - Ref: " << primRefCov << " bp (Est. Ratio: " << primRefRatio << "%)\n"
                  << "  - Qry: " << primQryCov << " bp (Est. Ratio: " << primQryRatio << "%)\n"
                  << "--- Total Coverage (Primary + Uncovered Secondary) ---\n"
                  << "  - Ref: " << totalRefCov << " bp (Est. Ratio: " << totalRefRatio << "%)\n"
                  << "  - Qry: " << totalQryCov << " bp (Est. Ratio: " << totalQryRatio << "%)\n";
    }
}


std::vector<std::pair<int, char>> mga::miniGlobalAlignment(const std::string& ref, const std::string& qry) {
    int n = ref.length();
    int m = qry.length();
    if (n == 0 && m == 0) return {};
    
    if (n == 0) return {{m, 'I'}};
    if (m == 0) return {{n, 'D'}};

    const int MATCH = 2, MISMATCH = -4, GAP = -4;
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    
    for (int i = 1; i <= n; i++) dp[i][0] = i * GAP;
    for (int j = 1; j <= m; j++) dp[0][j] = j * GAP;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score = (ref[i-1] == qry[j-1]) ? MATCH : MISMATCH;
            dp[i][j] = std::max({dp[i-1][j-1] + score, dp[i-1][j] + GAP, dp[i][j-1] + GAP});
        }
    }

    std::vector<std::pair<int, char>> cigar;
    int i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + ((ref[i-1] == qry[j-1]) ? MATCH : MISMATCH)) {
            cigar.push_back({1, (ref[i-1] == qry[j-1]) ? '=' : 'X'});
            i--; j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + GAP) {
            cigar.push_back({1, 'D'});
            i--;
        } else {
            cigar.push_back({1, 'I'});
            j--;
        }
    }
    std::reverse(cigar.begin(), cigar.end());
    
    std::vector<std::pair<int, char>> compressed_cigar;
    for (auto& op : cigar) {
        if (!compressed_cigar.empty() && compressed_cigar.back().second == op.second) {
            compressed_cigar.back().first += op.first;
        } else {
            compressed_cigar.push_back(op);
        }
    }
    return compressed_cigar;
}

int mga::getNearestBoundary(int pos, const std::set<int>& bnds, int threshold) {
    auto it = bnds.lower_bound(pos);
    int best = -1;
    int min_dist = threshold + 1;
    
    if (it != bnds.end() && *it - pos <= threshold) {
        best = *it;
        min_dist = *it - pos;
    }
    if (it != bnds.begin()) {
        --it;
        if (pos - *it < min_dist) {
            best = *it;
        }
    }
    return best;
}


/*
void mga::snapAlignmentsToBlockBoundaries(alnVec& alignments, BlockSet* bs1, BlockSet* bs2, const std::string& refSeq, const std::string& qrySeq, int threshold) {
    
    bool DEBUG_MODE = false; // 隨時可以關閉

    std::set<int> refBnds, qryBnds;
    int cur = 0;
    for (auto& b : bs1->getRepresentativeBlocks()) { refBnds.insert(cur); cur += bs1->getBlock(b)->getConsensus().length(); }
    refBnds.insert(cur);
    cur = 0;
    for (auto& b : bs2->getRepresentativeBlocks()) { qryBnds.insert(cur); cur += bs2->getBlock(b)->getConsensus().length(); }
    qryBnds.insert(cur);

    auto printCigarSummary = [](const mga::Cigar& c) {
        if (c.empty()) return std::string("[]");
        std::string s = "";
        if (c.size() <= 10) {
            for (auto& op : c) s += std::to_string(op.first) + op.second;
        } else {
            for (size_t i = 0; i < 5; ++i) s += std::to_string(c[i].first) + c[i].second;
            s += "...";
            for (size_t i = c.size() - 5; i < c.size(); ++i) s += std::to_string(c[i].first) + c[i].second;
        }
        return s;
    };

    auto getBndContext = [](int pos, const std::set<int>& bnds) {
        if (bnds.empty()) return std::string("[]");
        auto it = bnds.lower_bound(pos);
        std::string res = "[ ";
        if (it != bnds.begin()) { auto prev = it; --prev; res += std::to_string(*prev) + " "; }
        if (it != bnds.end()) {
            res += std::to_string(*it) + " ";
            auto next = it; ++next;
            if (next != bnds.end()) res += std::to_string(*next) + " ";
        }
        res += "]";
        return res;
    };

    for (auto& aln : alignments) {
        if (!aln.valid || aln.type != mga::PRIMARY) continue;

        if (DEBUG_MODE) {
            std::cout << "\n============================================================\n"
                      << "[DEBUG SNAP] Processing Alignment (Strand: " << (aln.inverse ? "-" : "+") << ")\n"
                      << "  Original Ref: [" << aln.refIdx.first << ", " << aln.refIdx.second << ")\n"
                      << "  Original Qry: [" << aln.qryIdx.first << ", " << aln.qryIdx.second << ")\n"
                      << "  Original CIGAR: " << printCigarSummary(aln.CIGAR) << "\n";
        }

        // ========================================================
        // 處理左端 (Left Boundary)
        // ========================================================
        int ref_start = aln.refIdx.first;
        int qry_start_cigar = (!aln.inverse) ? aln.qryIdx.first : aln.qryIdx.second;
        
        int nearest_ref = getNearestBoundary(ref_start, refBnds, threshold);
        int nearest_qry = getNearestBoundary(qry_start_cigar, qryBnds, threshold);

        if (DEBUG_MODE) {
            std::cout << "------------------------------------------------------------\n"
                      << " [Left Boundary Check]\n"
                      << "  Ref Start: " << ref_start << " | Nearest Bnds: " << getBndContext(ref_start, refBnds) << " -> Target: " << nearest_ref << "\n"
                      << "  Qry Start: " << qry_start_cigar << " | Nearest Bnds: " << getBndContext(qry_start_cigar, qryBnds) << " -> Target: " << nearest_qry << "\n";
        }

        if (nearest_ref != -1 && nearest_qry != -1) {
            int dist_r = ref_start - nearest_ref; 
            int dist_q = (!aln.inverse) ? (qry_start_cigar - nearest_qry) : (nearest_qry - qry_start_cigar);
            
            if (dist_r > 0 && dist_q > 0) {
                // Case 1-2 & 1-3: 兩邊都在外，使用 Global Alignment Extension
                int max_q_ext = (!aln.inverse) ? qry_start_cigar : (qrySeq.length() - qry_start_cigar);
                int ref_ext_len = std::min({dist_r, 500, ref_start});
                int qry_ext_len = std::min({dist_q, 500, max_q_ext});
                
                std::string ref_ext_seq = refSeq.substr(ref_start - ref_ext_len, ref_ext_len);
                std::string qry_ext_seq = (!aln.inverse) ? qrySeq.substr(qry_start_cigar - qry_ext_len, qry_ext_len) : getReverseComplement(qrySeq.substr(qry_start_cigar, qry_ext_len));

                auto ext_cigar = miniGlobalAlignment(ref_ext_seq, qry_ext_seq);
                aln.CIGAR.insert(aln.CIGAR.begin(), ext_cigar.begin(), ext_cigar.end());
                aln.CIGAR = mga::compressCigar(aln.CIGAR);
                
                aln.refIdx.first -= ref_ext_len;
                if (!aln.inverse) aln.qryIdx.first -= qry_ext_len;
                else aln.qryIdx.second += qry_ext_len;
                
                if (DEBUG_MODE) {
                    std::cout << "  ✅ Action: [EXTENSION] (Case 1-2 & 1-3)\n"
                              << "     Extending Ref by " << ref_ext_len << " | Seq: " << ref_ext_seq << "\n"
                              << "     Extending Qry by " << qry_ext_len << " | Seq: " << qry_ext_seq << "\n"
                              << "     CIGAR Updated: " << printCigarSummary(aln.CIGAR) << "\n";
                }
            } else {
                // Case 1-1 & 1-5: 退到內並 Pad I/D (Trim-and-Pad 邏輯)
                int r = aln.refIdx.first;
                int q = (!aln.inverse) ? aln.qryIdx.first : aln.qryIdx.second;
                mga::Cigar new_cigar;
                bool trimming = true;
                
                for (auto& op : aln.CIGAR) {
                    if (!trimming) { new_cigar.push_back(op); continue; }
                    int op_len_remaining = op.first;
                    for (int i = 0; i < op.first; ++i) {
                        bool ref_ok = (r >= nearest_ref);
                        bool qry_ok = (!aln.inverse) ? (q >= nearest_qry) : (q <= nearest_qry);
                        
                        if (ref_ok && qry_ok) {
                            trimming = false;
                            if (op_len_remaining > 0) new_cigar.push_back({op_len_remaining, op.second});
                            break;
                        }
                        if (op.second == 'M' || op.second == '=' || op.second == 'X') { r++; (!aln.inverse) ? q++ : q--; }
                        else if (op.second == 'D') { r++; }
                        else if (op.second == 'I') { (!aln.inverse) ? q++ : q--; }
                        op_len_remaining--;
                    }
                }
                
                int pad_r = r - nearest_ref;
                int pad_q = (!aln.inverse) ? (q - nearest_qry) : (nearest_qry - q);
                
                if (pad_q > 0) new_cigar.insert(new_cigar.begin(), {pad_q, 'I'});
                if (pad_r > 0) new_cigar.insert(new_cigar.begin(), {pad_r, 'D'});
                
                aln.CIGAR = mga::compressCigar(new_cigar);
                aln.refIdx.first = nearest_ref;
                if (!aln.inverse) aln.qryIdx.first = nearest_qry;
                else aln.qryIdx.second = nearest_qry;

                if (DEBUG_MODE) {
                    std::cout << "  ✂️ Action: [TRIM & PAD] (Case 1-1 & 1-5)\n"
                              << "     Trimmed until Ref=" << r << ", Qry=" << q << "\n"
                              << "     Padded with " << pad_r << "D and " << pad_q << "I\n"
                              << "     New Ref Start: " << aln.refIdx.first << " | New Qry Bound: " 
                              << ((!aln.inverse) ? std::to_string(aln.qryIdx.first) : std::to_string(aln.qryIdx.second)) << "\n"
                              << "     CIGAR Updated: " << printCigarSummary(aln.CIGAR) << "\n";
                }
            }
        } else {
            if (DEBUG_MODE) std::cout << "  ⚠️ Action: [SKIP] No valid boundary found within threshold.\n";
        }

        // ========================================================
        // 處理右端 (Right Boundary)
        // ========================================================
        int ref_end = aln.refIdx.second;
        int qry_end_cigar = (!aln.inverse) ? aln.qryIdx.second : aln.qryIdx.first;

        nearest_ref = getNearestBoundary(ref_end, refBnds, threshold);
        nearest_qry = getNearestBoundary(qry_end_cigar, qryBnds, threshold);

        if (DEBUG_MODE) {
            std::cout << "------------------------------------------------------------\n"
                      << " [Right Boundary Check]\n"
                      << "  Ref End: " << ref_end << " | Nearest Bnds: " << getBndContext(ref_end, refBnds) << " -> Target: " << nearest_ref << "\n"
                      << "  Qry End: " << qry_end_cigar << " | Nearest Bnds: " << getBndContext(qry_end_cigar, qryBnds) << " -> Target: " << nearest_qry << "\n";
        }

        if (nearest_ref != -1 && nearest_qry != -1) {
            int dist_r = nearest_ref - ref_end; 
            int dist_q = (!aln.inverse) ? (nearest_qry - qry_end_cigar) : (qry_end_cigar - nearest_qry);
            
            if (dist_r > 0 && dist_q > 0) {
                // Case 1-2 & 1-3: Extension
                int max_q_ext = (!aln.inverse) ? (qrySeq.length() - qry_end_cigar) : qry_end_cigar;
                int ref_ext_len = std::min({dist_r, 500, (int)refSeq.length() - ref_end});
                int qry_ext_len = std::min({dist_q, 500, max_q_ext});
                
                std::string ref_ext_seq = refSeq.substr(ref_end, ref_ext_len);
                std::string qry_ext_seq = (!aln.inverse) ? qrySeq.substr(qry_end_cigar, qry_ext_len) : getReverseComplement(qrySeq.substr(qry_end_cigar - qry_ext_len, qry_ext_len));

                auto ext_cigar = miniGlobalAlignment(ref_ext_seq, qry_ext_seq);
                aln.CIGAR.insert(aln.CIGAR.end(), ext_cigar.begin(), ext_cigar.end());
                aln.CIGAR = mga::compressCigar(aln.CIGAR);
                
                aln.refIdx.second += ref_ext_len;
                if (!aln.inverse) aln.qryIdx.second += qry_ext_len;
                else aln.qryIdx.first -= qry_ext_len;
                
                if (DEBUG_MODE) {
                    std::cout << "  ✅ Action: [EXTENSION] (Case 1-2 & 1-3)\n"
                              << "     Extending Ref by " << ref_ext_len << " | Seq: " << ref_ext_seq << "\n"
                              << "     Extending Qry by " << qry_ext_len << " | Seq: " << qry_ext_seq << "\n"
                              << "     CIGAR Updated: " << printCigarSummary(aln.CIGAR) << "\n";
                }
            } else {
                // Case 1-1 & 1-5: 退到內並 Pad I/D (Trim-and-Pad 邏輯 - 逆向走訪)
                int r = aln.refIdx.second;
                int q = (!aln.inverse) ? aln.qryIdx.second : aln.qryIdx.first;
                mga::Cigar new_cigar;
                bool trimming = true;
                
                for (auto it = aln.CIGAR.rbegin(); it != aln.CIGAR.rend(); ++it) {
                    auto op = *it;
                    if (!trimming) { new_cigar.push_back(op); continue; }
                    int op_len_remaining = op.first;
                    for (int i = 0; i < op.first; ++i) {
                        bool ref_ok = (r <= nearest_ref);
                        bool qry_ok = (!aln.inverse) ? (q <= nearest_qry) : (q >= nearest_qry);
                        
                        if (ref_ok && qry_ok) {
                            trimming = false;
                            if (op_len_remaining > 0) new_cigar.push_back({op_len_remaining, op.second});
                            break;
                        }
                        if (op.second == 'M' || op.second == '=' || op.second == 'X') { r--; (!aln.inverse) ? q-- : q++; }
                        else if (op.second == 'D') { r--; }
                        else if (op.second == 'I') { (!aln.inverse) ? q-- : q++; }
                        op_len_remaining--;
                    }
                }
                std::reverse(new_cigar.begin(), new_cigar.end());
                
                int pad_r = nearest_ref - r;
                int pad_q = (!aln.inverse) ? (nearest_qry - q) : (q - nearest_qry);
                
                if (pad_q > 0) new_cigar.push_back({pad_q, 'I'});
                if (pad_r > 0) new_cigar.push_back({pad_r, 'D'});
                
                aln.CIGAR = mga::compressCigar(new_cigar);
                aln.refIdx.second = nearest_ref;
                if (!aln.inverse) aln.qryIdx.second = nearest_qry;
                else aln.qryIdx.first = nearest_qry;

                if (DEBUG_MODE) {
                    std::cout << "  ✂️ Action: [TRIM & PAD] (Case 1-1 & 1-5)\n"
                              << "     Trimmed until Ref=" << r << ", Qry=" << q << "\n"
                              << "     Padded with " << pad_r << "D and " << pad_q << "I\n"
                              << "     New Ref End: " << aln.refIdx.second << " | New Qry Bound: " 
                              << ((!aln.inverse) ? std::to_string(aln.qryIdx.second) : std::to_string(aln.qryIdx.first)) << "\n"
                              << "     CIGAR Updated: " << printCigarSummary(aln.CIGAR) << "\n";
                }
            }
        } else {
            if (DEBUG_MODE) std::cout << "  ⚠️ Action: [SKIP] No valid boundary found within threshold.\n";
        }
    }
}
*/