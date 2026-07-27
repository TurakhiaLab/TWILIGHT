
#include "alignment.hpp"
#include "global_alignment.hpp"

#include <list>
#include <fstream>

#include <string>
#include <vector>



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

/*
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
*/

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

/*
void mga::identifyPrimaryAlignments(alnVec& alignments, chainVec& chains) {
    bool DEBUG_MODE = false; 
    
    PrimaryTracker refMask, qryMask;
    std::vector<Alignment> validAlns;
    std::vector<Alignment> remainingAlns;

    // 1. 攤平所有的 Alignment 並收集
    for (const auto& chain : chains) {
        for (int oldAlnId : chain.chainedAln) {
            auto& aln = alignments[oldAlnId];
            if (aln.valid && aln.alnScore > 0) {
                bool isRefMain = (aln.refName.find("_main") != std::string::npos);
                bool isQryMain = (aln.qryName.find("_main") != std::string::npos);

                if (isRefMain && isQryMain) {
                    validAlns.push_back(aln);
                } else {
                    aln.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(aln);
                }
            }
        }
    }

    std::cout << "Main Alignments Count: " << validAlns.size() << "\n";
    std::cout << "Remaining Alignments Count: " << remainingAlns.size() << "\n";
    if (remainingAlns.size() > 0) exit(1);

    // 2. 依照 Alignment Score 由大到小排序
    std::sort(validAlns.begin(), validAlns.end(), [](const Alignment& a, const Alignment& b) {
        return a.alnScore > b.alnScore;
    });

    std::vector<Alignment> primaryList;
    std::vector<Alignment> secondaryList;
    std::unordered_map<int, Alignment> debugPrimaryMap; 
    
    const int MIN_ALN_LEN = 100; 
    const int TOLERANCE = 50;  
    int nextGlobalId = 1;

    if (DEBUG_MODE) std::cout << "\n============================================\n"
                              << "=== PHASE 1: GREEDY PRIMARY ASSIGNMENT ===\n"
                              << "============================================\n";

    // ==========================================
    // Phase 1: Greedy 初步分類
    // ==========================================
    for (const auto& aln : validAlns) {
        std::set<int> rCuts, qCuts;

        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        refMask.getCuts(rMin, rMax, rCuts);

        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
        qryMask.getCuts(qMin, qMax, qCuts);

        std::vector<Alignment> fragments;
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

            std::set<int> refOverlapIds = refMask.getOverlappingIds(fRMin, fRMax);
            std::set<int> qryOverlapIds = qryMask.getOverlappingIds(fQMin, fQMax);

            bool refCov = !refOverlapIds.empty();
            bool qryCov = !qryOverlapIds.empty();

            frag.identifier = nextGlobalId++;
            frag.duplications.clear();
            frag.paralogs.clear();

            if (!refCov && !qryCov) {
                if (std::abs(fRMax - fRMin) >= MIN_ALN_LEN && std::abs(fQMax - fQMin) >= MIN_ALN_LEN) {
                    frag.type = mga::PRIMARY;
                    refMask.add(fRMin, fRMax, frag.identifier);
                    qryMask.add(fQMin, fQMax, frag.identifier);
                    primaryList.push_back(frag);
                    debugPrimaryMap[frag.identifier] = frag;
                } else {
                    frag.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(frag);
                }
            } else {
                if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
                    frag.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(frag);
                } else {
                    frag.type = mga::SECONDARY;
                    secondaryList.push_back(frag);
                }
            }
        }
    }

    if (DEBUG_MODE) std::cout << "\n============================================\n"
                              << "=== PHASE 2: GLOBAL ATOMIC SYNCHRONIZATION ===\n"
                              << "============================================\n";

    // ==========================================
    // Phase 2: 全局斷點收集與原子切割 
    // ==========================================
    std::vector<Alignment> validSecondaries;
    for (auto& sec : secondaryList) {
        if (std::abs(sec.refIdx.second - sec.refIdx.first) >= MIN_ALN_LEN) {
            validSecondaries.push_back(sec);
        }
    }
    secondaryList = std::move(validSecondaries);

    std::set<int> globalRefCuts;
    std::set<int> globalQryCuts;

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
    // [修改重點]: 移除 Secondary 的斷點收集，Secondary 無權決定 Graph 的全局切線

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

    std::vector<Alignment> finalPrimaries;
    std::vector<Alignment> finalSecondaries;
    PrimaryTracker finalRefMask, finalQryMask;
    std::unordered_map<int, Alignment> finalPrimaryMap;

    for (const auto& p : primaryList) {
        auto frags = splitSingleAlignment(p, globalRefCuts, globalQryCuts);
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
            finalRefMask.add(fRMin, fRMax, frag.identifier);
            finalQryMask.add(fQMin, fQMax, frag.identifier);
            finalPrimaries.push_back(frag);
            finalPrimaryMap[frag.identifier] = frag; 
        }
    }

    for (const auto& sec : secondaryList) {
        auto frags = splitSingleAlignment(sec, globalRefCuts, globalQryCuts);
        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
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
    // Phase 2-C: 完美一對一關係綁定 & 邊界強制對齊 (Boundary Snapping)
    // ==========================================
    std::map<int, std::vector<int>> primaryToDuplications;
    std::map<int, std::vector<int>> primaryToParalogs;

    auto getStrictOverlaps = [&](int sMin, int sMax, const std::set<int>& rawIds, bool isRef) {
        std::set<int> strictIds;
        for (int pId : rawIds) {
            auto it = finalPrimaryMap.find(pId);
            if (it != finalPrimaryMap.end()) {
                int pMin = isRef ? std::min(it->second.refIdx.first, it->second.refIdx.second)
                                 : std::min(it->second.qryIdx.first, it->second.qryIdx.second);
                int pMax = isRef ? std::max(it->second.refIdx.first, it->second.refIdx.second)
                                 : std::max(it->second.qryIdx.first, it->second.qryIdx.second);
                
                int oMin = std::max(sMin, pMin);
                int oMax = std::min(sMax, pMax);
                int overlapLen = oMax - oMin;

                if (overlapLen >= MIN_ALN_LEN / 2) strictIds.insert(pId);
            }
        }
        return strictIds;
    };

    // 【核心升級】：1-bp Stepping Engine，支援同時 Trimming (超出裁切) 與 Padding (不足補齊)
    auto padAlignmentToBoundary = [](Alignment& aln, int tgtRMin, int tgtRMax, int tgtQMin, int tgtQMax) {
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        if (rMin == tgtRMin && rMax == tgtRMax && qMin == tgtQMin && qMax == tgtQMax) return;

        // 計算要從 CIGAR 頭部砍掉 (Skip) 的量：當 Secondary 超出 Primary 的時候
        int skipR = std::max(0, tgtRMin - rMin);
        int skipQ = aln.inverse ? std::max(0, qMax - tgtQMax) : std::max(0, tgtQMin - qMin);

        // 目標的精準長度
        int targetR = tgtRMax - tgtRMin;
        int targetQ = tgtQMax - tgtQMin;

        mga::Cigar finalCigar;

        // 1. 如果 Secondary 不足，要在前面補 Gap
        int padRFront = std::max(0, rMin - tgtRMin);
        int padQFront = aln.inverse ? std::max(0, tgtQMax - qMax) : std::max(0, qMin - tgtQMin);
        
        if (padRFront > 0) finalCigar.push_back({padRFront, 'D'});
        if (padQFront > 0) finalCigar.push_back({padQFront, 'I'});

        // 2. 利用 1-bp 微步進，精準截取中間符合 Target 範圍的 CIGAR
        int curR = 0, curQ = 0; 
        int totalR = padRFront, totalQ = padQFront; 

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
                
                // 走訪 Ref
                if (rCons) {
                    if (curR < skipR) curR++; // Skip 掉超出的部分
                    else if (totalR < targetR) { useR = true; totalR++; } // 進入保留區
                }
                // 走訪 Qry
                if (qCons) {
                    if (curQ < skipQ) curQ++; // Skip 掉超出的部分
                    else if (totalQ < targetQ) { useQ = true; totalQ++; } // 進入保留區
                }

                // 根據雙方的走訪狀況組合 CIGAR
                if (useR && useQ) addOp((t == 'M' || t == '=' || t == 'X') ? t : 'M'); 
                else if (useR) addOp('D'); 
                else if (useQ) addOp('I'); 
            }
        }

        // 3. 如果到了尾巴 Secondary 還是不足，補齊 Gap
        if (totalR < targetR) {
            addOp('D'); finalCigar.back().first += (targetR - totalR - 1);
        }
        if (totalQ < targetQ) {
            addOp('I'); finalCigar.back().first += (targetQ - totalQ - 1);
        }

        // 4. 清理與合併相鄰同類的 CIGAR 碎片
        mga::Cigar cleanCigar;
        for (auto op : finalCigar) {
            if (!cleanCigar.empty() && cleanCigar.back().second == op.second) {
                cleanCigar.back().first += op.first;
            } else {
                cleanCigar.push_back(op);
            }
        }
        aln.CIGAR = cleanCigar;

        // 5. 套用精準的新座標
        if (aln.refIdx.first < aln.refIdx.second) {
            aln.refIdx.first = tgtRMin; aln.refIdx.second = tgtRMax;
        } else {
            aln.refIdx.first = tgtRMax; aln.refIdx.second = tgtRMin;
        }

        if (aln.qryIdx.first < aln.qryIdx.second) {
            aln.qryIdx.first = tgtQMin; aln.qryIdx.second = tgtQMax;
        } else {
            aln.qryIdx.first = tgtQMax; aln.qryIdx.second = tgtQMin;
        }
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

        // 執行完美補償與裁切
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

    // [Pass 2-D]: 關係寫回 Primary
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

    alignments.clear();
    alignments.insert(alignments.end(), finalPrimaries.begin(), finalPrimaries.end());
    alignments.insert(alignments.end(), finalSecondaries.begin(), finalSecondaries.end());
    alignments.insert(alignments.end(), remainingAlns.begin(), remainingAlns.end());

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

    // if (DEBUG_MODE) {
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
*/

struct CutPoint {
    int r;
    int q;
    bool operator<(const CutPoint& other) const {
        return r < other.r; 
    }
};

// 輔助函式：根據 Reference 座標，沿著 CIGAR 找出精準的 Query 2D 對應點
CutPoint getCutFromRef(const mga::Alignment& aln, int targetR) {
    int rPos = aln.refIdx.first;
    int qPos = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first;
    int qDir = aln.inverse ? -1 : 1;
    for (const auto& op : aln.CIGAR) {
        int len = op.first; char t = op.second;
        bool cR = (t == 'M' || t == '=' || t == 'X' || t == 'D');
        bool cQ = (t == 'M' || t == '=' || t == 'X' || t == 'I');
        if (cR && targetR >= std::min(rPos, rPos + len) && targetR <= std::max(rPos, rPos + len)) {
            int diff = std::abs(targetR - rPos);
            return {targetR, qPos + (cQ ? diff * qDir : 0)};
        }
        if (cR) rPos += len;
        if (cQ) qPos += len * qDir;
    }
    return {targetR, qPos}; 
};

// 【核心引擎 A】：完美 2D 變形與補齊引擎 (Snapping & Padding)
void padAlignmentToBoundary(mga::Alignment& aln, int tgtRMin, int tgtRMax, int tgtQMin, int tgtQMax) {
    int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
    int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
    int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
    int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

    if (rMin == tgtRMin && rMax == tgtRMax && qMin == tgtQMin && qMax == tgtQMax) return;
    
    int skipRFront = std::max(0, tgtRMin - rMin);
    int skipRBack = std::max(0, rMax - tgtRMax);
    int skipQFront = aln.inverse ? std::max(0, qMax - tgtQMax) : std::max(0, tgtQMin - qMin);
    int skipQBack = aln.inverse ? std::max(0, qMin - tgtQMin) : std::max(0, qMax - tgtQMax);

    int padRFront = std::max(0, rMin - tgtRMin);
    int padRBack = std::max(0, tgtRMax - rMax);
    int padQFront = aln.inverse ? std::max(0, tgtQMax - qMax) : std::max(0, qMin - tgtQMin);
    int padQBack = aln.inverse ? std::max(0, qMin - tgtQMin) : std::max(0, tgtQMax - qMax);

    int targetMiddleR = std::max(0, std::min(rMax, tgtRMax) - std::max(rMin, tgtRMin));
    int targetMiddleQ = std::max(0, std::min(qMax, tgtQMax) - std::max(qMin, tgtQMin));

    auto alignGap = [&](int gapR, int gapQ, bool isFront) -> mga::Cigar {
        mga::Cigar res;
        if (gapR == 0 && gapQ == 0) return res;
        if (gapR == 0) { res.push_back({gapQ, 'I'}); return res; }
        if (gapQ == 0) { res.push_back({gapR, 'D'}); return res; }
        int alignLenR = gapR;
        int alignLenQ = gapQ;
        if (gapR >= 100 || gapQ >= 100) {
            int minLen = std::min(gapR, gapQ);
            alignLenR = std::min(gapR, minLen * 2);
            alignLenQ = std::min(gapQ, minLen * 2);
        }
        int leftoverR = gapR - alignLenR;
        int leftoverQ = gapQ - alignLenQ;
        // [TODO] 這裡可以抽換成真實的序列抓取與 NW/KSW2
        mga::Cigar alignedCigar;
        int dummyMatch = std::min(alignLenR, alignLenQ);
        if (dummyMatch > 0) alignedCigar.push_back({dummyMatch, 'M'});
        if (alignLenR > dummyMatch) alignedCigar.push_back({alignLenR - dummyMatch, 'D'});
        if (alignLenQ > dummyMatch) alignedCigar.push_back({alignLenQ - dummyMatch, 'I'});
        if (isFront) {
            if (leftoverR > 0) res.push_back({leftoverR, 'D'});
            if (leftoverQ > 0) res.push_back({leftoverQ, 'I'});
            res.insert(res.end(), alignedCigar.begin(), alignedCigar.end());
        } else {
            res.insert(res.end(), alignedCigar.begin(), alignedCigar.end());
            if (leftoverR > 0) res.push_back({leftoverR, 'D'});
            if (leftoverQ > 0) res.push_back({leftoverQ, 'I'});
        }
        return res;
    };
    
    mga::Cigar finalCigar;
    mga::Cigar frontCigar = alignGap(padRFront, padQFront, true);
    finalCigar.insert(finalCigar.end(), frontCigar.begin(), frontCigar.end());

    int curR = 0, curQ = 0, savedR = 0, savedQ = 0;
    auto addMiddleOp = [&](char t) {
        if (!finalCigar.empty() && finalCigar.back().second == t) finalCigar.back().first++;
        else finalCigar.push_back({1, t});
    };

    for (auto op : aln.CIGAR) {
        int l = op.first; char t = op.second;
        if (t == 'S' || t == 'H') t = 'I';
        bool rCons = (t == 'M' || t == '=' || t == 'X' || t == 'D');
        bool qCons = (t == 'M' || t == '=' || t == 'X' || t == 'I');
        for (int i = 0; i < l; ++i) {
            bool useR = rCons && (curR >= skipRFront) && (savedR < targetMiddleR);
            bool useQ = qCons && (curQ >= skipQFront) && (savedQ < targetMiddleQ);
            if (rCons) curR++;
            if (qCons) curQ++;
            if (useR && useQ) { addMiddleOp((t == 'M' || t == '=' || t == 'X') ? t : 'M'); savedR++; savedQ++; } 
            else if (useR) { addMiddleOp('D'); savedR++; } 
            else if (useQ) { addMiddleOp('I'); savedQ++; }
        }
    }

    mga::Cigar backCigar = alignGap(padRBack, padQBack, false);
    finalCigar.insert(finalCigar.end(), backCigar.begin(), backCigar.end());

    mga::Cigar cleanCigar;
    for (auto op : finalCigar) {
        if (!cleanCigar.empty() && cleanCigar.back().second == op.second) cleanCigar.back().first += op.first;
        else if (op.first > 0) cleanCigar.push_back(op);
    }
    aln.CIGAR = cleanCigar;

    if (aln.refIdx.first < aln.refIdx.second) { aln.refIdx.first = tgtRMin; aln.refIdx.second = tgtRMax; } 
    else { aln.refIdx.first = tgtRMax; aln.refIdx.second = tgtRMin; }

    if (aln.qryIdx.first < aln.qryIdx.second) { aln.qryIdx.first = tgtQMin; aln.qryIdx.second = tgtQMax; } 
    else { aln.qryIdx.first = tgtQMax; aln.qryIdx.second = tgtQMin; }
};

// 【核心引擎 B】：完美 2D 切割引擎 (2D Perfect Splitting)
std::vector<mga::Alignment> splitSingleAlignment2D (const mga::Alignment& aln, const std::vector<CutPoint>& syncedCuts, int tolerance) {
    std::vector<mga::Alignment> frags;
    if (syncedCuts.empty()) {
        frags.push_back(aln);
        return frags;
    }

    std::vector<CutPoint> waypoints;
    waypoints.push_back({std::min(aln.refIdx.first, aln.refIdx.second), 
                         aln.inverse ? std::max(aln.qryIdx.first, aln.qryIdx.second) : std::min(aln.qryIdx.first, aln.qryIdx.second)});
    
    for (const auto& cut : syncedCuts) {
        if (cut.r > waypoints.front().r && cut.r < std::max(aln.refIdx.first, aln.refIdx.second)) {
            waypoints.push_back(cut);
        }
    }
    std::sort(waypoints.begin(), waypoints.end());
    waypoints.push_back({std::max(aln.refIdx.first, aln.refIdx.second), 
                         aln.inverse ? std::min(aln.qryIdx.first, aln.qryIdx.second) : std::max(aln.qryIdx.first, aln.qryIdx.second)});
    int currentCutIdx = 1;
    int rPos = waypoints[0].r, qPos = waypoints[0].q;
    int qDir = aln.inverse ? -1 : 1;
    mga::Alignment currAln = aln;
    currAln.CIGAR.clear();
    int currRStart = rPos, currQStart = qPos;
    for (const auto& op : aln.CIGAR) {
        if (currentCutIdx >= waypoints.size()) break;
        int len = op.first; char type = op.second;
        while (len > 0 && currentCutIdx < waypoints.size()) {
            int targetR = waypoints[currentCutIdx].r;
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');
            int step = len;
            if (consumesRef) {
                int distR = targetR - rPos;
                if (distR > 0 && distR < step) step = distR;
            }
            if (!currAln.CIGAR.empty() && currAln.CIGAR.back().second == type) currAln.CIGAR.back().first += step; 
            else currAln.CIGAR.push_back({step, type});
            if (consumesRef) rPos += step;
            if (consumesQry) qPos += step * qDir;
            len -= step;
            // 強制在 Target R 切斷並對齊 Target Q
            if (consumesRef && rPos == targetR) {
                currAln.refIdx.first = currRStart; currAln.refIdx.second = rPos;
                currAln.qryIdx.first = currQStart; currAln.qryIdx.second = qPos;
                int targetQ = waypoints[currentCutIdx].q;
                int tgtRMin = std::min(currRStart, targetR);
                int tgtRMax = std::max(currRStart, targetR);
                int tgtQMin = std::min(currQStart, targetQ);
                int tgtQMax = std::max(currQStart, targetQ);
                
                padAlignmentToBoundary(currAln, tgtRMin, tgtRMax, tgtQMin, tgtQMax);
                frags.push_back(currAln);
                currRStart = targetR; currQStart = targetQ;
                rPos = targetR; qPos = targetQ;
                currAln = aln; currAln.CIGAR.clear();
                currentCutIdx++;
            }
        }
    }
    // 收尾防呆
    if (std::abs(waypoints.back().r - currRStart) > tolerance || std::abs(waypoints.back().q - currQStart) > tolerance) {
         if (!currAln.CIGAR.empty() || currRStart != waypoints.back().r || currQStart != waypoints.back().q) {
            currAln.refIdx.first = currRStart; currAln.refIdx.second = waypoints.back().r;
            currAln.qryIdx.first = currQStart; currAln.qryIdx.second = waypoints.back().q;
            
            int tgtRMin = std::min(currRStart, waypoints.back().r);
            int tgtRMax = std::max(currRStart, waypoints.back().r);
            int tgtQMin = std::min(currQStart, waypoints.back().q);
            int tgtQMax = std::max(currQStart, waypoints.back().q);
            
            padAlignmentToBoundary(currAln, tgtRMin, tgtRMax, tgtQMin, tgtQMax);
            frags.push_back(currAln);
         }
    }
    return frags;
};

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

std::unordered_map<std::string, std::vector<int>> collect_and_classified_alignments(mga::alnVec& alignments, mga::alnVec& mainAlns, mga::alnVec& remaining2remainingAlns, mga::alnVec& remainingAlns) {
    // ==========================================
    // Phase 0. Collect and Classified Alignments
    // ==========================================
    std::unordered_map<std::string, int> bestRemAlnIds;
    for (int alnId = 0; alnId < alignments.size(); ++alnId) {
        auto& aln = alignments[alnId];
        if (aln.valid && aln.alnScore > 0) {
            bool isRefMain = (aln.refName.find("_main") != std::string::npos);
            bool isQryMain = (aln.qryName.find("_main") != std::string::npos);
            if (!isRefMain) {
                if (bestRemAlnIds.find(aln.refName) == bestRemAlnIds.end() || aln.alnScore > alignments[bestRemAlnIds[aln.refName]].alnScore) bestRemAlnIds[aln.refName] = alnId;
            }
            if (!isQryMain) {
                if (bestRemAlnIds.find(aln.qryName) == bestRemAlnIds.end() || aln.alnScore > alignments[bestRemAlnIds[aln.qryName]].alnScore) bestRemAlnIds[aln.qryName] = alnId;
            }
        }
    }

    std::unordered_set<int> topRemIds;
    for (const auto& pair : bestRemAlnIds) topRemIds.insert(pair.second);
    
    std::unordered_map<std::string, std::vector<int>> remainingSeqToAlns;

    for (int alnId = 0; alnId < alignments.size(); ++alnId) {
        auto& aln = alignments[alnId];
        if (aln.valid && aln.alnScore > 0) {
            bool isRefMain = (aln.refName.find("_main") != std::string::npos);
            bool isQryMain = (aln.qryName.find("_main") != std::string::npos);
            if (isRefMain && isQryMain) mainAlns.push_back(aln);
            else if (topRemIds.count(alnId)) {
                if (isRefMain || isQryMain) mainAlns.push_back(aln);
                else remaining2remainingAlns.push_back(aln);
            }
            else { aln.type = mga::REMAINING_ALN; remainingAlns.push_back(aln); }
        }
    }

    std::sort(mainAlns.begin(), mainAlns.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        return a.alnScore > b.alnScore;
    });

    return remainingSeqToAlns;
}

void mga::identifyPrimaryAlignments(alnVec& alignments, chainVec& chains, stringMap& ref_seqs, stringMap& qry_seqs) {
    bool DEBUG_MODE = true; 

    std::cout << "Total Alignments: " << alignments.size() << '\n';
    std::vector<Alignment> mainAlns;
    std::vector<Alignment> remaining2remainingAlns;
    std::vector<Alignment> remainingAlns;
    auto remainingSeqToAlns = collect_and_classified_alignments(alignments, mainAlns, remaining2remainingAlns, remainingAlns);


    // ==========================================
    // Phase 1: Build Main Backbone (Alignments without any overlaps)
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n=== PHASE 1: BUILDING BACKBONE ===\n";
    
    std::unordered_map<std::string, PrimaryTracker> trackers;
    std::vector<Alignment> backbonePrimaries;
    std::vector<Alignment> finalSecondaries;
    std::unordered_map<int, Alignment> backbonePrimaryMap;
    
    const int MIN_ALN_LEN = 100; 
    const int MAX_OVERLAP = 50;
    const int TOLERANCE = 50;  
    int nextGlobalId = 1;

    for (auto& aln : mainAlns) {
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        std::set<int> rCuts, qCuts;
        trackers[aln.refName].getCuts(rMin, rMax, rCuts);
        trackers[aln.qryName].getCuts(qMin, qMax, qCuts);

        std::vector<Alignment> fragments;
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

            std::set<int> refOverlapIds = trackers[frag.refName].getOverlappingIds(fRMin, fRMax);
            std::set<int> qryOverlapIds = trackers[frag.qryName].getOverlappingIds(fQMin, fQMax);

            if (refOverlapIds.empty() && qryOverlapIds.empty()) {
                if (std::abs(fRMax - fRMin) >= MIN_ALN_LEN && std::abs(fQMax - fQMin) >= MIN_ALN_LEN) {
                    mga::Alignment p = frag;
                    p.identifier = nextGlobalId++;
                    p.type = mga::PRIMARY;
                    trackers[p.refName].add(fRMin, fRMax, p.identifier);
                    trackers[p.qryName].add(fQMin, fQMax, p.identifier);
                    backbonePrimaries.push_back(p);
                    backbonePrimaryMap[p.identifier] = p;

                    bool isRefMain = (p.refName.find("_main") != std::string::npos);
                    bool isQryMain = (p.qryName.find("_main") != std::string::npos);
                    if (!isRefMain) remainingSeqToAlns[p.refName].push_back(p.identifier);
                    if (!isQryMain) remainingSeqToAlns[p.qryName].push_back(p.identifier);
                } else {
                    mga::Alignment rem = frag;
                    rem.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(rem);
                }
            } else {
                if (std::abs(fRMax - fRMin) >= MIN_ALN_LEN && std::abs(fQMax - fQMin) >= MIN_ALN_LEN) {
                    mga::Alignment sec = frag;
                    sec.type = mga::SECONDARY;
                    finalSecondaries.push_back(sec);

                    bool isRefMain = (sec.refName.find("_main") != std::string::npos);
                    bool isQryMain = (sec.qryName.find("_main") != std::string::npos);
                    if (!isRefMain) remainingSeqToAlns[sec.refName].push_back(sec.identifier);
                    if (!isQryMain) remainingSeqToAlns[sec.qryName].push_back(sec.identifier);
                } else {
                    mga::Alignment rem = frag;
                    rem.type = mga::REMAINING_ALN;
                    remainingAlns.push_back(rem);
                }
            }
        }
    }

    // ==========================================
    // Phase 1.5: Sort Backbone Primaries into Chains
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n=== PHASE 1.5: SORT BACKBONE INTO CHAINS ===\n";

    std::vector<Alignment> qrySorted = backbonePrimaries;
    std::sort(qrySorted.begin(), qrySorted.end(), [](const Alignment& a, const Alignment& b) {
        int minA = std::min(a.qryIdx.first, a.qryIdx.second);
        int minB = std::min(b.qryIdx.first, b.qryIdx.second);
        if (minA == minB) return a.identifier < b.identifier;
        return minA < minB;
    });

    std::unordered_map<int, int> qryRank;
    for (size_t i = 0; i < qrySorted.size(); ++i) {
        qryRank[qrySorted[i].identifier] = i;
    }

    std::sort(backbonePrimaries.begin(), backbonePrimaries.end(), [](const Alignment& a, const Alignment& b) {
        int minA = std::min(a.refIdx.first, a.refIdx.second);
        int minB = std::min(b.refIdx.first, b.refIdx.second);
        if (minA == minB) return a.identifier < b.identifier;
        return minA < minB;
    });

    std::vector<std::vector<Alignment>> backboneChains;
    std::vector<Alignment> currentChain;

    for (const auto& aln : backbonePrimaries) {
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
                
                if (currRStart >= prevREnd) {
                    int currQRank = qryRank[aln.identifier];
                    int prevQRank = qryRank[prev.identifier];

                    if (!aln.inverse) { // Forward strand
                        int prevQEnd = std::max(prev.qryIdx.first, prev.qryIdx.second);
                        int currQStart = std::min(aln.qryIdx.first, aln.qryIdx.second);
                        if (currQStart >= prevQEnd && currQRank == prevQRank + 1) {
                            validNext = true;
                        }
                    } else { // Reverse strand: Ref increases but Qry decreases
                        int prevQStart = std::min(prev.qryIdx.first, prev.qryIdx.second);
                        int currQEnd = std::max(aln.qryIdx.first, aln.qryIdx.second);
                        if (currQEnd <= prevQStart && currQRank == prevQRank - 1) {
                            validNext = true;
                        }
                    }
                }
            }
            
            if (validNext) {
                currentChain.push_back(aln);
            } else {
                backboneChains.push_back(currentChain);
                currentChain.clear();
                currentChain.push_back(aln);
            }
        }
    }
    if (!currentChain.empty()) {
        backboneChains.push_back(currentChain);
    }

    if (DEBUG_MODE) {
        int totalBackboneAlignments = 0;
        std::cout << "Identified " << backboneChains.size() << " backbone chains.\n";
        for (size_t i = 0; i < backboneChains.size(); ++i) {
            std::cout << "  Chain " << i + 1 << ": " << backboneChains[i].size() << " alignments.\n";
            for (const auto& aln : backboneChains[i]) {
                ++totalBackboneAlignments;
                std::string cigarString = "";
                for (auto op: aln.CIGAR) cigarString += (std::to_string(op.first) + op.second);
                std::cout << "    Ref: [" << aln.refIdx.first << ", " << aln.refIdx.second << ") | "
                          << "Qry: [" << aln.qryIdx.first << ", " << aln.qryIdx.second << ") | "
                          << "Strand: " << (aln.inverse ? "-" : "+") << " | "
                          << "CIGAR: " << cigarString << "\n";
            }
        }
        std::cout << "Total Backbone Alignments: " << totalBackboneAlignments << "\n";
        std::cout << "Total Secondary Alignments: " << finalSecondaries.size() << "\n";
        for (const auto& aln : finalSecondaries) {
            std::string cigarString = "";
            for (auto op: aln.CIGAR) cigarString += (std::to_string(op.first) + op.second);
            std::cout << "    Ref: [" << aln.refIdx.first << ", " << aln.refIdx.second << ") | "
                      << "Qry: [" << aln.qryIdx.first << ", " << aln.qryIdx.second << ") | "
                      << "Strand: " << (aln.inverse ? "-" : "+") << " | "
                      << "CIGAR: " << cigarString << "\n";
        }
    }

    // ==========================================
    // Phase 2: Gap Filling within Chains
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n=== PHASE 2: GAP FILLING ===\n";

    const int MAX_GAP_SIZE = 2000;
    const float MIN_IDENTITY = 0.85f; // Only accept alignments with identity higher than MIN_IDENTITY

    std::vector<std::vector<Alignment>> filledChains;

    for (auto& chain : backboneChains) {
        if (chain.empty()) continue;

        std::vector<Alignment> currentFilledChain;
        currentFilledChain.push_back(chain[0]);

        for (size_t i = 1; i < chain.size(); ++i) {
            const auto& prev = chain[i - 1];
            const auto& curr = chain[i];

            int prevREnd = std::max(prev.refIdx.first, prev.refIdx.second);
            int currRStart = std::min(curr.refIdx.first, curr.refIdx.second);
            int gapR = currRStart - prevREnd;

            int gapQ = 0;
            int prevQEnd, currQStart;
            
            if (!curr.inverse) {
                prevQEnd = std::max(prev.qryIdx.first, prev.qryIdx.second);
                currQStart = std::min(curr.qryIdx.first, curr.qryIdx.second);
                gapQ = currQStart - prevQEnd;
            } else {
                // Reverse strand: prev 雖然在 Ref 上較前，但在 Qry 上座標較大
                prevQEnd = std::min(prev.qryIdx.first, prev.qryIdx.second);
                currQStart = std::max(curr.qryIdx.first, curr.qryIdx.second);
                gapQ = prevQEnd - currQStart; // 反向計算距離
            }

            // 1. 檢查是否超過最大 Gap 限制
            if (gapR > MAX_GAP_SIZE || gapQ > MAX_GAP_SIZE || gapR < 0 || gapQ < 0) {
                // 不填補，直接加入目前的 Alignment
                currentFilledChain.push_back(curr);
                continue;
            }

            // 2. 判斷長度一致性，決定對齊策略
            int diff = std::abs(gapR - gapQ);
            int minGap = std::min(gapR, gapQ);
            
            // Threshold: 基礎容忍度 20bp，或是較小 Gap 的 15% (確保 block merging 的長度合理性)
            int threshold = std::max(20, static_cast<int>(minGap * 0.15));

            // 提取序列
            std::string refSeq = "";
            std::string qrySeq = "";
            if (ref_seqs.count(curr.refName) && qry_seqs.count(curr.qryName)) {
                refSeq = ref_seqs[curr.refName].substr(prevREnd, gapR);
                if (!curr.inverse) {
                    qrySeq = qry_seqs[curr.qryName].substr(prevQEnd, gapQ);
                } else {
                    // 對於反股，要用 currQStart 當起點往右抓 gapQ 長度，然後作 Reverse Complement
                    qrySeq = qry_seqs[curr.qryName].substr(currQStart, gapQ);
                    qrySeq = getReverseComplement(qrySeq);
                }
            }

            AlnResult result;
            bool isSemiGlobal = (diff > threshold);

            if (refSeq.size() > 0 && qrySeq.size() > 0) {
                if (diff <= threshold) {
                    // 長度接近 -> Global Alignment (Needleman-Wunsch)
                    if (DEBUG_MODE) std::cout << "Filling Gap (Global): Ref=" << gapR << ", Qry=" << gapQ << "\n";
                    result = runGlobalAlignment(refSeq, qrySeq);
                } else {
                    // 長度懸殊 -> Semi-Global Alignment
                    if (DEBUG_MODE) std::cout << "Filling Gap (Semi-Global): Ref=" << gapR << ", Qry=" << gapQ << "\n";
                    result = runSemiGlobalAlignment(refSeq, qrySeq);
                }
                std::cout << result.success << '\t' << result.identity << '\t' << result.score << '\t' << result.cigar << '\n';
            }

            // 3. 檢查 Identity 並決定是否 Accept
            if (result.success && result.identity >= MIN_IDENTITY) {
                auto cigar = parser::parseCigar(result.cigar);
                
                int trimRHead = 0, trimQHead = 0;
                int trimRTail = 0, trimQTail = 0;

                if (isSemiGlobal && !cigar.empty()) {
                    auto& frontOp = cigar.front();
                    if (frontOp.second == 'I' || frontOp.second == 'D') {
                        if (frontOp.first >= MIN_ALN_LEN) {
                            if (frontOp.second == 'I') trimQHead = frontOp.first;
                            if (frontOp.second == 'D') trimRHead = frontOp.first;
                        }
                    }
                    
                    auto& backOp = cigar.back();
                    if (cigar.size() > 1 && (backOp.second == 'I' || backOp.second == 'D')) {
                        if (backOp.first >= MIN_ALN_LEN) {
                            if (backOp.second == 'I') trimQTail = backOp.first;
                            if (backOp.second == 'D') trimRTail = backOp.first;
                        }
                    }
                    
                    if (trimRHead > 0 || trimQHead > 0) forceTrimLeft(cigar, trimRHead, trimQHead);
                    if (trimRTail > 0 || trimQTail > 0) forceTrimRight(cigar, trimRTail, trimQTail);
                }

                if (!cigar.empty()) {
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
                }
            }

            currentFilledChain.push_back(curr);
        }
        filledChains.push_back(currentFilledChain);
    }

    // ==========================================
    // Phase 2.5: Merge Contiguous Alignments
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n=== PHASE 2.5: MERGE CONTIGUOUS ALIGNMENTS ===\n";

    std::vector<std::vector<Alignment>> mergedChains;
    std::vector<Alignment> newBackbonePrimaries;

    for (auto& chain : filledChains) {
        if (chain.empty()) continue;
        
        std::vector<Alignment> currentMergedChain;
        Alignment mergedAln = chain[0];
        
        for (size_t i = 1; i < chain.size(); ++i) {
            const Alignment& curr = chain[i];
            
            bool contiguousRef = (mergedAln.refIdx.second == curr.refIdx.first);
            bool contiguousQry = false;
            
            if (!mergedAln.inverse) {
                contiguousQry = (mergedAln.qryIdx.second == curr.qryIdx.first);
            } else {
                contiguousQry = (mergedAln.qryIdx.first == curr.qryIdx.second);
            }
            
            if (contiguousRef && contiguousQry && (mergedAln.inverse == curr.inverse)) {
                // Merge curr into mergedAln
                mergedAln.refIdx.second = curr.refIdx.second;
                
                if (!mergedAln.inverse) {
                    mergedAln.qryIdx.second = curr.qryIdx.second;
                } else {
                    mergedAln.qryIdx.first = curr.qryIdx.first; // 反股，範圍是往左延伸 (數值越小)
                }
                
                // Append CIGAR & Compress
                mergedAln.CIGAR.insert(mergedAln.CIGAR.end(), curr.CIGAR.begin(), curr.CIGAR.end());
                mergedAln.CIGAR = compressCigar(mergedAln.CIGAR);
                
                mergedAln.alnScore += curr.alnScore;
            } else {
                // Not contiguous, push the mergedAln and start a new one
                currentMergedChain.push_back(mergedAln);
                mergedAln = curr;
            }
        }
        currentMergedChain.push_back(mergedAln);
        mergedChains.push_back(currentMergedChain);
        
        newBackbonePrimaries.insert(newBackbonePrimaries.end(), currentMergedChain.begin(), currentMergedChain.end());
    }

    // 更新回原本的 backbone variables
    backbonePrimaries = newBackbonePrimaries;
    backbonePrimaryMap.clear();
    trackers.clear();

    for (auto& aln : backbonePrimaries) {
        backbonePrimaryMap[aln.identifier] = aln;
        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);
        trackers[aln.refName].add(rMin, rMax, aln.identifier);
        trackers[aln.qryName].add(qMin, qMax, aln.identifier);
    }
    backboneChains = mergedChains;

    if (DEBUG_MODE) {
        std::cout << "Identified " << backboneChains.size() << " backbone chains AFTER merging.\n";
        for (size_t i = 0; i < backboneChains.size(); ++i) {
            std::cout << "  Chain " << i + 1 << ": " << backboneChains[i].size() << " alignments.\n";
            for (const auto& aln : backboneChains[i]) {
                std::string cigarString = "";
                for (auto op: aln.CIGAR) cigarString += (std::to_string(op.first) + op.second);
                std::cout << "    Ref: [" << aln.refIdx.first << ", " << aln.refIdx.second << ") | "
                          << "Qry: [" << aln.qryIdx.first << ", " << aln.qryIdx.second << ") | "
                          << "Strand: " << (aln.inverse ? "-" : "+") << " | "
                          << "CIGAR: " << cigarString << "\n";
            }
        }
    }

    // ==========================================
    // Phase 3: Integrate Secondary Alignments
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n=== PHASE 3: INTEGRATE SECONDARY ALIGNMENTS ===\n";

    std::set<int> globalRefCuts;
    std::set<int> globalQryCuts;

    for (const auto& p : backbonePrimaries) {
        globalRefCuts.insert(std::min(p.refIdx.first, p.refIdx.second));
        globalRefCuts.insert(std::max(p.refIdx.first, p.refIdx.second));
        globalQryCuts.insert(std::min(p.qryIdx.first, p.qryIdx.second));
        globalQryCuts.insert(std::max(p.qryIdx.first, p.qryIdx.second));
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

    std::vector<Alignment> processedSecondaries;
    std::map<int, std::vector<int>> primaryToDuplications;
    std::map<int, std::vector<int>> primaryToParalogs;

    for (const auto& sec : finalSecondaries) {
        auto frags = splitSingleAlignment(sec, globalRefCuts, globalQryCuts);
        for (auto& frag : frags) {
            int fRMin = std::min(frag.refIdx.first, frag.refIdx.second);
            int fRMax = std::max(frag.refIdx.first, frag.refIdx.second);
            int fQMin = std::min(frag.qryIdx.first, frag.qryIdx.second);
            int fQMax = std::max(frag.qryIdx.first, frag.qryIdx.second);

            if (std::abs(fRMax - fRMin) < MIN_ALN_LEN || std::abs(fQMax - fQMin) < MIN_ALN_LEN) {
                continue; 
            }

            std::set<int> rRawIds = trackers[frag.refName].getOverlappingIds(fRMin, fRMax);
            std::set<int> qRawIds = trackers[frag.qryName].getOverlappingIds(fQMin, fQMax);

            auto getStrictOverlaps = [&](int sMin, int sMax, const std::set<int>& rawIds, bool isRef) {
                std::set<int> strictIds;
                for (int pId : rawIds) {
                    auto it = backbonePrimaryMap.find(pId);
                    if (it != backbonePrimaryMap.end()) {
                        int pMin = isRef ? std::min(it->second.refIdx.first, it->second.refIdx.second)
                                         : std::min(it->second.qryIdx.first, it->second.qryIdx.second);
                        int pMax = isRef ? std::max(it->second.refIdx.first, it->second.refIdx.second)
                                         : std::max(it->second.qryIdx.first, it->second.qryIdx.second);
                        
                        int oMin = std::max(sMin, pMin);
                        int oMax = std::min(sMax, pMax);
                        if (oMax - oMin >= MIN_ALN_LEN / 2) strictIds.insert(pId);
                    }
                }
                return strictIds;
            };

            std::set<int> rIds = getStrictOverlaps(fRMin, fRMax, rRawIds, true);
            std::set<int> qIds = getStrictOverlaps(fQMin, fQMax, qRawIds, false);

            int targetRMin = fRMin, targetRMax = fRMax;
            int targetQMin = fQMin, targetQMax = fQMax;

            for (int pId : rIds) {
                auto& pAln = backbonePrimaryMap[pId];
                int pRMin = std::min(pAln.refIdx.first, pAln.refIdx.second);
                int pRMax = std::max(pAln.refIdx.first, pAln.refIdx.second);
                if (std::abs(fRMin - pRMin) <= TOLERANCE) targetRMin = pRMin;
                if (std::abs(fRMax - pRMax) <= TOLERANCE) targetRMax = pRMax;
            }

            for (int pId : qIds) {
                auto& pAln = backbonePrimaryMap[pId];
                int pQMin = std::min(pAln.qryIdx.first, pAln.qryIdx.second);
                int pQMax = std::max(pAln.qryIdx.first, pAln.qryIdx.second);
                if (std::abs(fQMin - pQMin) <= TOLERANCE) targetQMin = pQMin;
                if (std::abs(fQMax - pQMax) <= TOLERANCE) targetQMax = pQMax;
            }

            // === Boundary Sequence Alignment Padding Logic ===
            int padRFront = std::max(0, fRMin - targetRMin);
            int padRBack = std::max(0, targetRMax - fRMax);
            int padQFront = frag.inverse ? std::max(0, targetQMax - fQMax) : std::max(0, fQMin - targetQMin);
            int padQBack = frag.inverse ? std::max(0, fQMin - targetQMin) : std::max(0, targetQMax - fQMax);

            int skipR = std::max(0, targetRMin - fRMin);
            int skipQ = frag.inverse ? std::max(0, fQMax - targetQMax) : std::max(0, targetQMin - fQMin);

            auto alignGap = [&](int gapR, int gapQ, bool isFront) -> mga::Cigar {
                mga::Cigar res;
                if (gapR == 0 && gapQ == 0) return res;
                if (gapR == 0) { res.push_back({gapQ, 'I'}); return res; }
                if (gapQ == 0) { res.push_back({gapR, 'D'}); return res; }
                
                std::string refGapSeq = "";
                std::string qryGapSeq = "";
                if (ref_seqs.count(frag.refName) && qry_seqs.count(frag.qryName)) {
                    if (isFront) {
                        refGapSeq = ref_seqs[frag.refName].substr(targetRMin, gapR);
                        qryGapSeq = qry_seqs[frag.qryName].substr(frag.inverse ? targetQMax - gapQ : targetQMin, gapQ);
                    } else {
                        refGapSeq = ref_seqs[frag.refName].substr(fRMax, gapR);
                        qryGapSeq = qry_seqs[frag.qryName].substr(frag.inverse ? targetQMin : fQMax, gapQ);
                    }
                    if (frag.inverse) qryGapSeq = getReverseComplement(qryGapSeq);
                }
                
                if (!refGapSeq.empty() && !qryGapSeq.empty()) {
                    AlnResult alnRes = runGlobalAlignment(refGapSeq, qryGapSeq);
                    if (alnRes.success) return parser::parseCigar(alnRes.cigar);
                }
                
                // Fallback: 如果抓不到序列 (極端例外狀況)，退回基礎 I/D Padding
                int dummyMatch = std::min(gapR, gapQ);
                res.push_back({dummyMatch, 'M'});
                if (gapR > dummyMatch) res.push_back({gapR - dummyMatch, 'D'});
                if (gapQ > dummyMatch) res.push_back({gapQ - dummyMatch, 'I'});
                return res;
            };

            mga::Cigar frontCigar = alignGap(padRFront, padQFront, true);
            mga::Cigar backCigar = alignGap(padRBack, padQBack, false);

            // 核心片段提取 (剔除超出邊界的部分)
            mga::Cigar midCigar;
            int curR = 0, curQ = 0;
            int targetMidR = fRMax - fRMin - skipR;
            int targetMidQ = fQMax - fQMin - skipQ;
            int savedR = 0, savedQ = 0;

            auto addMidOp = [&](char t) {
                if (!midCigar.empty() && midCigar.back().second == t) midCigar.back().first++;
                else midCigar.push_back({1, t});
            };

            for (auto op : frag.CIGAR) {
                int l = op.first; char t = op.second;
                if (t == 'S' || t == 'H') t = 'I';
                bool rCons = (t == 'M' || t == '=' || t == 'X' || t == 'D');
                bool qCons = (t == 'M' || t == '=' || t == 'X' || t == 'I');

                for (int i = 0; i < l; ++i) {
                    bool useR = false, useQ = false;
                    if (rCons) { if (curR >= skipR && savedR < targetMidR) useR = true; curR++; }
                    if (qCons) { if (curQ >= skipQ && savedQ < targetMidQ) useQ = true; curQ++; }

                    if (useR && useQ) addMidOp((t == 'M' || t == '=' || t == 'X') ? t : 'M');
                    else if (useR) addMidOp('D');
                    else if (useQ) addMidOp('I');
                    
                    if (useR) savedR++;
                    if (useQ) savedQ++;
                }
            }

            mga::Cigar finalCigar;
            finalCigar.insert(finalCigar.end(), frontCigar.begin(), frontCigar.end());
            finalCigar.insert(finalCigar.end(), midCigar.begin(), midCigar.end());
            finalCigar.insert(finalCigar.end(), backCigar.begin(), backCigar.end());
            
            frag.CIGAR = compressCigar(finalCigar);
            frag.refIdx.first = (frag.refIdx.first < frag.refIdx.second) ? targetRMin : targetRMax;
            frag.refIdx.second = (frag.refIdx.first < frag.refIdx.second) ? targetRMax : targetRMin;
            frag.qryIdx.first = (frag.qryIdx.first < frag.qryIdx.second) ? targetQMin : targetQMax;
            frag.qryIdx.second = (frag.qryIdx.first < frag.qryIdx.second) ? targetQMax : targetQMin;
            frag.identifier = nextGlobalId++;

            // === Relation Binding ===
            bool refCov = !rIds.empty();
            bool qryCov = !qIds.empty();

            if (refCov && qryCov) {
                for (int rPId : rIds) {
                    for (int qPId : qIds) {
                        if (rPId != qPId) { 
                            primaryToParalogs[rPId].push_back(qPId);
                            primaryToParalogs[qPId].push_back(rPId);
                        }
                        frag.paralogs.push_back(rPId);
                        frag.paralogs.push_back(qPId);
                    }
                }
                for (int pId : rIds) primaryToDuplications[pId].push_back(frag.identifier);
                for (int pId : qIds) primaryToDuplications[pId].push_back(frag.identifier);
                
            } else if (refCov && !qryCov) {
                for (int pId : rIds) {
                    frag.duplications.push_back(pId); 
                    primaryToDuplications[pId].push_back(frag.identifier);
                }
            } else if (!refCov && qryCov) {
                for (int pId : qIds) {
                    frag.duplications.push_back(pId);
                    primaryToDuplications[pId].push_back(frag.identifier);
                }
            }

            processedSecondaries.push_back(frag);
        }
    }

    finalSecondaries = processedSecondaries;

    // Write relations back to Primary
    for (auto& p : backbonePrimaries) {
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

    alignments.clear();
    alignments.insert(alignments.end(), backbonePrimaries.begin(), backbonePrimaries.end());
    alignments.insert(alignments.end(), finalSecondaries.begin(), finalSecondaries.end());
    alignments.insert(alignments.end(), remainingAlns.begin(), remainingAlns.end());
    
    
    // ==========================================
    // Phase 4: Coverage Calculation & Debugging
    // ==========================================
    std::vector<std::pair<int, int>> refPrimInts, qryPrimInts;
    std::vector<std::pair<int, int>> refTotalInts, qryTotalInts;
    int maxRefCoord = 0, maxQryCoord = 0;

    for (const auto& aln : backbonePrimaries) {
        auto rPair = std::make_pair(std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second));
        auto qPair = std::make_pair(std::min(aln.qryIdx.first, aln.qryIdx.second), std::max(aln.qryIdx.first, aln.qryIdx.second));
        refPrimInts.push_back(rPair); qryPrimInts.push_back(qPair);
        refTotalInts.push_back(rPair); qryTotalInts.push_back(qPair);
        if (aln.refName.find("_main") != std::string::npos) {
            auto rPair = std::make_pair(std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second));
            refPrimInts.push_back(rPair); refTotalInts.push_back(rPair);
        }
        if (aln.qryName.find("_main") != std::string::npos) {
            auto qPair = std::make_pair(std::min(aln.qryIdx.first, aln.qryIdx.second), std::max(aln.qryIdx.first, aln.qryIdx.second));
            qryPrimInts.push_back(qPair); qryTotalInts.push_back(qPair);
        }
    }
    for (const auto& aln : finalSecondaries) {
        refTotalInts.push_back({std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second)});
        qryTotalInts.push_back({std::min(aln.qryIdx.first, aln.qryIdx.second), std::max(aln.qryIdx.first, aln.qryIdx.second)});
        if (aln.refName.find("_main") != std::string::npos)
            refTotalInts.push_back({std::min(aln.refIdx.first, aln.refIdx.second), std::max(aln.refIdx.first, aln.refIdx.second)});
        if (aln.qryName.find("_main") != std::string::npos)
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
            if (intervals[i].first <= currentEnd) { currentEnd = std::max(currentEnd, intervals[i].second); } 
            else { totalCovered += (currentEnd - currentStart); currentStart = intervals[i].first; currentEnd = intervals[i].second; }
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
                std::cout << "Primary ID: " << pAln->identifier 
                          << " | Ref: [" << pAln->refIdx.first << ", " << pAln->refIdx.second << "]"
                          << " | Qry: [" << pAln->qryIdx.first << ", " << pAln->qryIdx.second << "]\tStrand: " << pAln->inverse << '\n';
            }

            std::cout << "\n>>> SECONDARY ALIGNMENTS <<<\n";
            for (const auto* sAln : secondaryAlns) {
                std::cout << "Secondary ID: " << sAln->identifier 
                          << " | Ref: [" << sAln->refIdx.first << ", " << sAln->refIdx.second << "]"
                          << " | Qry: [" << sAln->qryIdx.first << ", " << sAln->qryIdx.second << "]\tStrand: " << sAln->inverse << '\n';
            }
        }
        std::cout << "\n[Merger] Identified " << primaryAlns.size() << " Primary and " << secondaryAlns.size() << " Secondary alignments.\n";
        std::cout << "\n[Coverage Info]\n"
                  << "--- Primary Coverage ---\n"
                  << "  - Ref: " << primRefCov << " bp (" << primRefRatio << "%)\n"
                  << "  - Qry: " << primQryCov << " bp (" << primQryRatio << "%)\n"
                  << "--- Total Coverage ---\n"
                  << "  - Ref: " << totalRefCov << " bp (" << totalRefRatio << "%)\n"
                  << "  - Qry: " << totalQryCov << " bp (" << totalQryRatio << "%)\n";
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

std::string mga::getReverseComplement(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break; case 'T': c = 'A'; break;
            case 'C': c = 'G'; break; case 'G': c = 'C'; break;
            case 'a': c = 't'; break; case 't': c = 'a'; break;
            case 'c': c = 'g'; break; case 'g': c = 'c'; break;
        }
    }
    return seq;
}



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
                aln.CIGAR = compressCigar(aln.CIGAR);
                
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
                
                aln.CIGAR = compressCigar(new_cigar);
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
                aln.CIGAR = compressCigar(aln.CIGAR);
                
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
                
                aln.CIGAR = compressCigar(new_cigar);
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