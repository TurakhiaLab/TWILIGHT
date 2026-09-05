
#include "alignment.hpp"
#include "global_alignment.hpp"
#include "type.hpp"
#include "block_set.hpp"
#include "coordinate_manager.hpp"
#include "timer.hpp"
#include "cigar_util.hpp"

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

int Alignment::mapCoordinate(int targetPos, bool targetIsRef) {
    auto aln = *this;
    int rPos = aln.refIdx.first;
    int qStart = std::min(aln.qryIdx.first, aln.qryIdx.second);
    int qEnd = std::max(aln.qryIdx.first, aln.qryIdx.second);
    int qCurr = aln.inverse ? qEnd : qStart;
        
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
            if (targetIsRef && targetPos >= rPos && targetPos <= rPos + len) return qCurr;
            rPos += len;
        } else if (type == 'I' || type == 'S' || type == 'H') {
            if (!targetIsRef && type == 'I' && targetPos >= std::min(qCurr, qNext) && targetPos <= std::max(qCurr, qNext)) 
                return rPos;
            qCurr = qNext;
        }
    }
    return -1;
};


void Alignment::finalizeAlignment(
    CoverageTracker& ref_coverageTracker, 
    CoverageTracker& qry_coverageTracker, 
    const std::map<int, std::string>& r_cut_reasons, 
    const std::map<int, std::string>& q_cut_reasons, 
    int snapLength) 
{
    bool debug = false;
    
    int check_q_start = std::min(this->qryIdx.first, this->qryIdx.second);
    int check_q_end   = std::max(this->qryIdx.first, this->qryIdx.second);
    int check_r_start = this->refIdx.first;
    int check_r_end   = this->refIdx.second;
    int total_q_overlap = qry_coverageTracker.getOverlapLength(check_q_start, check_q_end);
    int total_r_overlap = ref_coverageTracker.getOverlapLength(check_r_start, check_r_end);

    if (total_q_overlap == 0 && total_r_overlap == 0) {
        int qL = qry_coverageTracker.distLeft(check_q_start);
        int qR = qry_coverageTracker.distRight(check_q_end);
        int rL = ref_coverageTracker.distLeft(check_r_start);
        int rR = ref_coverageTracker.distRight(check_r_end);

        auto isStrictR = [&](int r) { return r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos; };
        auto isStrictQ = [&](int q) { return q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos; };

        bool need_snap = false;
        int q_pad_L = 0, q_pad_R = 0, r_pad_L = 0, r_pad_R = 0;

        if (!isStrictQ(check_q_start) && qL > 0 && qL <= snapLength) { q_pad_L = qL; need_snap = true; }
        if (!isStrictQ(check_q_end)   && qR > 0 && qR <= snapLength) { q_pad_R = qR; need_snap = true; }
        if (!isStrictR(check_r_start) && rL > 0 && rL <= snapLength) { r_pad_L = rL; need_snap = true; }
        if (!isStrictR(check_r_end)   && rR > 0 && rR <= snapLength) { r_pad_R = rR; need_snap = true; }

        if (need_snap) {
            if (debug) std::cout << "    -> [DEBUG] 🧲 SNAPPING Final Fragment! Ref: (+" << r_pad_L << ", +" << r_pad_R << "), Qry: (+" << q_pad_L << ", +" << q_pad_R << ")\n";
            snapAlignment(*this, r_pad_L, r_pad_R, q_pad_L, q_pad_R);
        }
    }
}

// =========================================================
// Helper 2: CIGAR 修剪工具 (從頭或尾強制修剪 N 個實體鹼基)
// =========================================================
static std::pair<int, int> trimCIGARSide(CigarString& cigar, int trimLen, bool fromFront) {
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
}

// =========================================================
// Helper 3: 錯位補償與微小碎片吸收 (原 Stage 5 核心)
// =========================================================
static void absorbMicroFragments(
    std::vector<Alignment>& frags, 
    const std::map<int, std::string>& r_cut_reasons, 
    const std::map<int, std::string>& q_cut_reasons, 
    int L_min, bool debug) 
{
    if (frags.size() < 2) return;
    std::vector<bool> absorbed(frags.size(), false);

    auto isStrictCut = [&](int r, int q) -> bool {
        if (r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos) return true;
        if (q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos) return true;
        return false;
    };

    // 處理中間的 Fragment (i = 1 to size-2)
    if (frags.size() >= 3) {
        for (size_t i = 1; i < frags.size() - 1; ++i) {
            int f_q_len = std::abs(frags[i].qryIdx.second - frags[i].qryIdx.first);
            int f_r_len = std::abs(frags[i].refIdx.second - frags[i].refIdx.first);
            int r_left_cut = frags[i].refIdx.first;  int r_right_cut = frags[i].refIdx.second;
            int q_left_cut = frags[i].qryIdx.first;  int q_right_cut = frags[i].qryIdx.second;

            auto isStrictR = [&](int r) { return r_cut_reasons.count(r) && r_cut_reasons.at(r).find("STRICT") != std::string::npos; };
            auto isStrictQ = [&](int q) { return q_cut_reasons.count(q) && q_cut_reasons.at(q).find("STRICT") != std::string::npos; };

            if ((isStrictR(r_left_cut) && isStrictR(r_right_cut)) || (isStrictQ(q_left_cut) && isStrictQ(q_right_cut))) continue;

            bool is_q_micro = (f_q_len > 0 && f_q_len < 10);
            bool is_r_micro = (f_r_len > 0 && f_r_len < 10);

            if ((is_q_micro || is_r_micro) && !absorbed[i-1] && !absorbed[i+1]) {
                bool is_pure_match = (frags[i].CIGAR.size() == 1 && frags[i].CIGAR[0].second == 'M');
                if (is_pure_match && f_r_len >= 30) continue; 

                int left_r_len = std::abs(frags[i-1].refIdx.second - frags[i-1].refIdx.first);
                int left_q_len = std::abs(frags[i-1].qryIdx.second - frags[i-1].qryIdx.first);
                int right_r_len = std::abs(frags[i+1].refIdx.second - frags[i+1].refIdx.first);
                int right_q_len = std::abs(frags[i+1].qryIdx.second - frags[i+1].qryIdx.first);

                bool ref_goes_left = (isStrictR(r_right_cut) && !isStrictR(r_left_cut)) ? true : 
                                     (isStrictR(r_left_cut) && !isStrictR(r_right_cut)) ? false : (left_r_len >= right_r_len); 

                bool qry_goes_left = (isStrictQ(q_right_cut) && !isStrictQ(q_left_cut)) ? true : 
                                     (isStrictQ(q_left_cut) && !isStrictQ(q_right_cut)) ? false : (left_q_len >= right_q_len);

                int left_append_r = (is_r_micro && ref_goes_left) ? f_r_len : 0;
                int right_prepend_r = (is_r_micro && !ref_goes_left) ? f_r_len : 0;
                int left_append_q = (is_q_micro && qry_goes_left) ? f_q_len : 0;
                int right_prepend_q = (is_q_micro && !qry_goes_left) ? f_q_len : 0;

                // 修改左側
                if (left_append_r > 0 && left_append_q > 0) frags[i-1].CIGAR.insert(frags[i-1].CIGAR.end(), frags[i].CIGAR.begin(), frags[i].CIGAR.end());
                else if (left_append_r > 0) frags[i-1].CIGAR.push_back({left_append_r, 'D'});
                else if (left_append_q > 0) frags[i-1].CIGAR.push_back({left_append_q, 'I'});

                if (left_append_r > 0) frags[i-1].refIdx.second = frags[i].refIdx.second;
                if (left_append_q > 0) {
                    if (!frags[i-1].inverse) frags[i-1].qryIdx.second = frags[i].qryIdx.second;
                    else                     frags[i-1].qryIdx.first  = frags[i].qryIdx.first;
                }
                frags[i-1].updateAlnLength();

                // 修改右側
                if (right_prepend_r > 0 && right_prepend_q > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), frags[i].CIGAR.begin(), frags[i].CIGAR.end());
                else if (right_prepend_r > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), {right_prepend_r, 'D'});
                else if (right_prepend_q > 0) frags[i+1].CIGAR.insert(frags[i+1].CIGAR.begin(), {right_prepend_q, 'I'});

                if (right_prepend_r > 0) frags[i+1].refIdx.first = frags[i].refIdx.first;
                if (right_prepend_q > 0) {
                    if (!frags[i+1].inverse) frags[i+1].qryIdx.first = frags[i].qryIdx.first;
                    else                     frags[i+1].qryIdx.second = frags[i].qryIdx.second;
                }
                frags[i+1].updateAlnLength();

                if (is_r_micro && is_q_micro) absorbed[i] = true;
                else {
                    if (is_r_micro) { 
                        frags[i].CIGAR = {{f_q_len, 'I'}};
                        frags[i].refIdx.first = ref_goes_left ? frags[i].refIdx.second : frags[i].refIdx.first;
                        frags[i].refIdx.second = frags[i].refIdx.first;
                    } else if (is_q_micro) { 
                        frags[i].CIGAR = {{f_r_len, 'D'}};
                        if (qry_goes_left) frags[i].qryIdx.first = frags[i].qryIdx.second;
                        else               frags[i].qryIdx.second = frags[i].qryIdx.first;
                    }
                    frags[i].updateAlnLength();
                }
            }
        }
    }

    // 處理頭部 (Prefix)
    int f0_q_len = std::abs(frags[0].qryIdx.second - frags[0].qryIdx.first);
    int f0_r_len = std::abs(frags[0].refIdx.second - frags[0].refIdx.first);
    if (!absorbed[0] && (f0_q_len < L_min || f0_r_len < L_min)) {
        absorbed[0] = true; // 簡化處理：頭部雜訊太小直接丟棄 (可依需求展開為原先詳細邏輯)
    }

    // 處理尾部 (Suffix)
    size_t last = frags.size() - 1;
    int fl_q_len = std::abs(frags[last].qryIdx.second - frags[last].qryIdx.first);
    int fl_r_len = std::abs(frags[last].refIdx.second - frags[last].refIdx.first);
    if (!absorbed[last] && (fl_q_len < L_min || fl_r_len < L_min)) {
        absorbed[last] = true;
    }

    // 🌟 在函數內直接過濾掉被吸收的碎塊，還給主程式一個乾淨的 Vector！
    std::vector<Alignment> cleaned_frags;
    for (size_t i = 0; i < frags.size(); ++i) {
        if (!absorbed[i]) cleaned_frags.push_back(std::move(frags[i]));
    }
    frags = std::move(cleaned_frags);
}

// =========================================================
// Helper Structures for getBestAlignments
// =========================================================
struct CutPoint { 
    bool is_strict = false; 
    std::string reasons = ""; 
};

struct CutPointCmp {
    bool inverse;
    CutPointCmp(bool inv) : inverse(inv) {}
    bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
        if (a.first != b.first) return a.first < b.first;
        return inverse ? (a.second > b.second) : (a.second < b.second);
    }
};

// =========================================================
// Static Helper Functions for getBestAlignments
// =========================================================

static int getDynamicDistToBound(const CoordinateManager& coordMgr, int pos, bool isRef, bool lookLeft) {
    const auto& intervalMap = isRef ? coordMgr.getRefIntervals() : coordMgr.getQryIntervals();
    if (intervalMap.empty()) return -1;

    auto it = intervalMap.upper_bound(pos);

    if (lookLeft) {
        int max_left_bound = -1;
        if (it != intervalMap.begin()) {
            auto prevIt = std::prev(it);
            if (pos >= prevIt->first && pos <= prevIt->second.end) {
                if (pos == prevIt->first || pos == prevIt->second.end) return 0;
                max_left_bound = std::max(max_left_bound, prevIt->first);
            } else if (prevIt->second.end <= pos) {
                max_left_bound = std::max(max_left_bound, prevIt->second.end);
            }
        }
        return (max_left_bound != -1) ? (pos - max_left_bound) : -1;
    } else {
        int min_right_bound = -1;
        if (it != intervalMap.begin()) {
            auto prevIt = std::prev(it);
            if (pos >= prevIt->first && pos <= prevIt->second.end) {
                if (pos == prevIt->first || pos == prevIt->second.end) return 0;
                min_right_bound = prevIt->second.end;
            }
        }
        if (it != intervalMap.end()) {
            if (min_right_bound == -1 || it->first < min_right_bound) {
                min_right_bound = it->first;
            }
        }
        return (min_right_bound != -1) ? (min_right_bound - pos) : -1;
    }
}

/**
 * @brief Trims the alignment edges if they overlap with already-covered genomic regions.
 * 
 * This function checks if the head or tail of the current best alignment overlaps with 
 * existing alignments tracked by the CoverageTracker. If the overlap length is smaller 
 * than L_min, it shaves off the CIGAR operations and updates the coordinate boundaries 
 * to prevent double-counting of aligned regions.
 * 
 * @param best The alignment to be trimmed (modified in place).
 * @param ref_tracker Coverage tracker for the reference axis.
 * @param qry_tracker Coverage tracker for the query axis.
 * @param L_min Minimum alignment length threshold.
 */
static void preTrimAlignment(Alignment& best, const CoverageTracker& ref_tracker, const CoverageTracker& qry_tracker, int L_min) {
    
    int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
    int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
    int check_r_start = best.refIdx.first;
    int check_r_end   = best.refIdx.second;

    int ref_head_ovl = ref_tracker.getLeftOverlap(check_r_start, check_r_end);
    int qry_head_ovl = best.inverse ? qry_tracker.getRightOverlap(check_q_start, check_q_end) 
                                    : qry_tracker.getLeftOverlap(check_q_start, check_q_end);
    int head_trim = std::max(ref_head_ovl, qry_head_ovl);

    if (head_trim > 0 && head_trim < L_min) {
        auto shaved = trimCIGARSide(best.CIGAR, head_trim, true);
        best.refIdx.first += shaved.first;
        if (!best.inverse) best.qryIdx.first += shaved.second; 
        else               best.qryIdx.second -= shaved.second; 
    }
    
    check_r_start = best.refIdx.first;
    check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
    check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);

    int ref_tail_ovl = ref_tracker.getRightOverlap(check_r_start, check_r_end);
    int qry_tail_ovl = best.inverse ? qry_tracker.getLeftOverlap(check_q_start, check_q_end)
                                    : qry_tracker.getRightOverlap(check_q_start, check_q_end);
    int tail_trim = std::max(ref_tail_ovl, qry_tail_ovl);

    if (tail_trim > 0 && tail_trim < L_min) {
        auto shaved = trimCIGARSide(best.CIGAR, tail_trim, false);
        best.refIdx.second -= shaved.first;
        if (!best.inverse) best.qryIdx.second -= shaved.second;
        else               best.qryIdx.first += shaved.second;
    }

    best.updateAlnLength();
}

static void snapToBoundaries(
    Alignment& best, 
    const CoordinateManager& coordMgr, 
    const std::string& refSeq, 
    const std::string& qrySeq, 
    int L_min, 
    bool debug) 
{
    int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
    int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
    int check_r_start = best.refIdx.first;
    int check_r_end   = best.refIdx.second;

    int qL = best.inverse ? getDynamicDistToBound(coordMgr, check_q_end, false, false) 
                          : getDynamicDistToBound(coordMgr, check_q_start, false, true);
    int rL = getDynamicDistToBound(coordMgr, check_r_start, true, true);
    int qR = best.inverse ? getDynamicDistToBound(coordMgr, check_q_start, false, true) 
                          : getDynamicDistToBound(coordMgr, check_q_end, false, false);
    int rR = getDynamicDistToBound(coordMgr, check_r_end, true, false);

    // --------------------------------------------------------
    // 處理左側 (Left Side / CIGAR Head) 
    // --------------------------------------------------------
    if (qL > 0 && rL > 0 && qL <= MAX_NW_LENGTH && rL <= MAX_NW_LENGTH) {
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
            std::string rev_ref(ref_patch.rbegin(), ref_patch.rend());
            std::string rev_qry(qry_patch.rbegin(), qry_patch.rend());
            patchCigar = runSemiGlobalAlignment(rev_ref, rev_qry);
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
            best.qryIdx.second += actual_qL; 
            check_q_end += actual_qL;
        }
        
    } else if (rL > 0 && rL <= SNAP_LENGTH && (qL == 0 || qL == -1 || qL >= MAX_NW_LENGTH)) {
        if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Deletion Gap (Left) to edge...\n";
        best.CIGAR.insert(best.CIGAR.begin(), {rL, 'D'});
        best.refIdx.first -= rL; check_r_start -= rL;
    } else if (qL > 0 && qL <= SNAP_LENGTH && (rL == 0 || rL == -1 || rL >= MAX_NW_LENGTH)) {
        if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Insertion Gap (Left) to edge...\n";
        best.CIGAR.insert(best.CIGAR.begin(), {qL, 'I'});
        if (!best.inverse) {
            best.qryIdx.first -= qL;
            check_q_start -= qL;
        } else {
            best.qryIdx.second += qL; 
            check_q_end += qL;
        }
    }

    // --------------------------------------------------------
    // 處理右側 (Right Side / CIGAR Tail)
    // --------------------------------------------------------
    if (qR > 0 && rR > 0 && qR <= MAX_NW_LENGTH && rR <= MAX_NW_LENGTH) {
        int actual_rR = std::min(rR, (int)refSeq.length() - check_r_end);
        int actual_qR = 0;
        std::string ref_patch = refSeq.substr(check_r_end, actual_rR);
        std::string qry_patch = "";

        if (!best.inverse) {
            actual_qR = std::min(qR, (int)qrySeq.length() - check_q_end);
            qry_patch = qrySeq.substr(check_q_end, actual_qR);
        } else {
            int safe_q_start = std::max(0, check_q_start - qR);
            actual_qR = check_q_start - safe_q_start;
            qry_patch = qrySeq.substr(safe_q_start, actual_qR);
        }

        CigarString patchCigar;
        bool use_semi = (std::max(actual_rR, actual_qR) >= 10 * std::min(actual_rR, actual_qR));

        if (debug) {
            std::cout << "  -> [STAGE 0] 🔗 Micro-overlap resolved. Sealing Right Wall (" << rR << ","<< qR << ")..."
                      << (use_semi ? " [Using Semi-Global]" : " [Using Global]");
        }

        if (use_semi) {
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
            best.qryIdx.first -= actual_qR; 
            check_q_start -= actual_qR;
        }
        
    } else if (rR > 0 && rR <= SNAP_LENGTH && (qR == 0 || qR == -1 || qR >= MAX_NW_LENGTH)) {
        if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Deletion Gap (Right) to edge...\n";
        best.CIGAR.push_back({rR, 'D'});
        best.refIdx.second += rR; check_r_end += rR;
    } else if (qR > 0 && qR <= SNAP_LENGTH && (rR == 0 || rR == -1 || rR >= MAX_NW_LENGTH)) {
        if (debug) std::cout << "  -> [STAGE 0] 🩹 Appending Insertion Gap (Right) to edge...\n";
        best.CIGAR.push_back({qR, 'I'});
        if (!best.inverse) {
            best.qryIdx.second += qR;
            check_q_end += qR;
        } else {
            best.qryIdx.first -= qR; 
            check_q_start -= qR;
        }
    }
}

static std::map<std::pair<int, int>, CutPoint, CutPointCmp> collectCutPoints(
    const Alignment& best, 
    const CoordinateManager& coordMgr, 
    const CoverageTracker& ref_tracker, 
    const CoverageTracker& qry_tracker, 
    int L_min, 
    bool debug) 
{
    int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
    int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
    int check_r_start = best.refIdx.first;
    int check_r_end   = best.refIdx.second;

    std::map<std::pair<int, int>, CutPoint, CutPointCmp> unified_cuts(CutPointCmp(best.inverse));

    auto registerCut = [&](int r, int q, bool strict, const std::string& reason) {
        if (r == -1 || q == -1) return; 
        auto& cut = unified_cuts[{r, q}];
        cut.is_strict |= strict;
        if (!cut.reasons.empty()) cut.reasons += " + ";
        cut.reasons += reason;
    };

    auto mapRefToQry = [&](const Alignment& aln, int targetR) { 
        Alignment temp = aln; 
        return temp.mapCoordinate(targetR, true); 
    };
    auto mapQryToRef = [&](const Alignment& aln, int targetQ) { 
        Alignment temp = aln; 
        return temp.mapCoordinate(targetQ, false); 
    };

    // --------------------------------------------------------
    // 1-1: Collect Strict Cuts From Previous BlockSet
    // --------------------------------------------------------
    std::set<int> tmp_q_cuts, tmp_r_cuts;
    qry_tracker.getCuts(check_q_start, check_q_end, tmp_q_cuts);
    ref_tracker.getCuts(check_r_start, check_r_end, tmp_r_cuts);
    
    if (debug && (!tmp_r_cuts.empty() || !tmp_q_cuts.empty())) {
        std::cout << "  [DEBUG-STAGE 1-1] 💎 STRICT Coverage Edge Detected!\n";
    }

    for (int r : tmp_r_cuts) {
        int q = mapRefToQry(best, r);
        if (debug) {
            std::cout << "     [STRICT-REF] Hit Coverage Border @ Ref: " << r 
                      << " | Mapped Qry: " << (q != -1 ? std::to_string(q) : "Gap/Out") 
                      << " (Likely due to trailing/heading Overlap)\n";
        }
        registerCut(r, q, true, "Coverage (STRICT Ref)");
    }

    for (int q : tmp_q_cuts) {
        int r = mapQryToRef(best, q);
        if (debug) {
            std::cout << "     [STRICT-QRY] Hit Coverage Border @ Qry: " << q 
                      << " | Mapped Ref: " << (r != -1 ? std::to_string(r) : "Gap/Out") 
                      << " (Likely due to trailing/heading Overlap)\n";
        }
        registerCut(r, q, true, "Coverage (STRICT Qry)");
    }

    // --------------------------------------------------------
    // 1-2: Collect Long Gaps From Alignments
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
    // 1-3: Dynamic Graph Boundaries (Powered by CoordinateManager)
    // --------------------------------------------------------
    if (debug) std::cout << "  [DEBUG-STAGE 1-3] Collecting Dynamic Graph Boundaries from CoordinateManager...\n";

    auto collectDynamicBoundaries = [&](bool isRef, int start_pos, int end_pos) {
        std::set<int> boundaries;
        const auto& intervalMap = isRef ? coordMgr.getRefIntervals() : coordMgr.getQryIntervals();
        if (intervalMap.empty()) return boundaries;

        auto it = intervalMap.upper_bound(start_pos);
        if (it != intervalMap.begin()) {
            --it;
        }

        while (it != intervalMap.end() && it->first <= end_pos) {
            if (it->first >= start_pos && it->first <= end_pos) {
                boundaries.insert(it->first);
            }
            if (it->second.end >= start_pos && it->second.end <= end_pos) {
                boundaries.insert(it->second.end);
            }
            ++it;
        }
        return boundaries;
    };

    std::set<int> dynamic_ref_bounds = collectDynamicBoundaries(true, check_r_start, check_r_end);
    std::set<int> dynamic_qry_bounds = collectDynamicBoundaries(false, check_q_start, check_q_end);

    for (int r_bound : dynamic_ref_bounds) {
        int q = mapRefToQry(best, r_bound);
        if (debug) {
            std::cout << "    ⚠️ [STRICT-DYNAMIC-REF] Recovered Topological Boundary @ Ref: " << r_bound 
                      << " | Mapped Qry: " << (q != -1 ? std::to_string(q) : "Gap/Out") << "\n";
        }
        registerCut(r_bound, q, true, "Dynamic Ref Bound [STRICT]");
    }

    for (int q_bound : dynamic_qry_bounds) {
        int r = mapQryToRef(best, q_bound);
        if (debug) {
            std::cout << "    ⚠️ [STRICT-DYNAMIC-QRY] Recovered Topological Boundary @ Qry: " << q_bound 
                      << " | Mapped Ref: " << (r != -1 ? std::to_string(r) : "Gap/Out") << "\n";
        }
        registerCut(r, q_bound, true, "Dynamic Qry Bound [STRICT]");
    }

    return unified_cuts;
}

static void evaluateCutClusters(
    const Alignment& best,
    const std::vector<std::pair<std::pair<int, int>, CutPoint>>& flat_cuts,
    BlockSet* refBlockSet,
    BlockSet* qryBlockSet,
    const CoordinateManager& coordMgr,
    int L_min,
    bool debug,
    std::set<int>& r_cuts,
    std::set<int>& q_cuts,
    std::map<int, std::string>& r_cut_reasons,
    std::map<int, std::string>& q_cut_reasons)
{
    int check_q_start = std::min(best.qryIdx.first, best.qryIdx.second);
    int check_q_end   = std::max(best.qryIdx.first, best.qryIdx.second);
    int check_r_start = best.refIdx.first;
    int check_r_end   = best.refIdx.second;

    std::vector<std::pair<int, int>> active_cuts;
    
    int c_idx = 0;
    while (c_idx < flat_cuts.size()) {
        global_timer.start("ecc_1_cluster_grouping");
        int start_idx = c_idx;
        int end_idx = c_idx;
        
        while (end_idx + 1 < flat_cuts.size()) {
            int r_dist = std::abs(flat_cuts[end_idx + 1].first.first - flat_cuts[end_idx].first.first);
            int q_dist = std::abs(flat_cuts[end_idx + 1].first.second - flat_cuts[end_idx].first.second);
            
            if (std::max(r_dist, q_dist) > L_min) break;
            end_idx++;
        }

        std::vector<int> strict_in_cluster;
        for (int k = start_idx; k <= end_idx; ++k) {
            if (flat_cuts[k].second.is_strict) {
                strict_in_cluster.push_back(k);
            }
        }
        global_timer.stop("ecc_1_cluster_grouping");

        int cluster_size = end_idx - start_idx + 1;
        if (debug && cluster_size > 1) {
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

        std::pair<int, int> A_pair = {check_r_start, best.inverse ? check_q_end : check_q_start};
        std::pair<int, int> E_pair = {check_r_end, best.inverse ? check_q_start : check_q_end};
        
        if (!active_cuts.empty()) {
            A_pair = active_cuts.back();
        }

        if (end_idx + 1 < flat_cuts.size()) {
            E_pair = flat_cuts[end_idx + 1].first;
        }

        int A_R = A_pair.first;
        int A_Q = A_pair.second;
        int E_R = E_pair.first;
        int E_Q = E_pair.second;

        if (start_idx == end_idx) {
            bool is_isolated = true;
            
            for (auto& c : active_cuts) {
                if (std::abs(c.first - B_R) < L_min && std::abs(c.second - B_Q) < L_min) {
                    is_isolated = false; break;
                }
            }

            int dist_to_ref_start = std::abs(B_R - check_r_start);
            int dist_to_ref_end   = std::abs(check_r_end - B_R);
            int dist_to_qry_start = std::abs(B_Q - check_q_start);
            int dist_to_qry_end   = std::abs(check_q_end - B_Q);

            bool is_strict = flat_cuts[start_idx].second.is_strict;

            bool is_near_boundary = ((dist_to_ref_start > 0 && dist_to_ref_start < MIN_DIST_2_BOUNDARY) || 
                                     (dist_to_ref_end > 0   && dist_to_ref_end < MIN_DIST_2_BOUNDARY) ||
                                     (dist_to_qry_start > 0 && dist_to_qry_start < MIN_DIST_2_BOUNDARY) || 
                                     (dist_to_qry_end > 0   && dist_to_qry_end < MIN_DIST_2_BOUNDARY));

            if (is_near_boundary) {
                if (is_strict) {
                    if (debug) std::cout << "    -> ⚠️ [BOUNDARY-STRICT] Cut @ R:" << B_R << ", Q:" << B_Q << " is near boundary but STRICT. Fast-Path granted.\n";
                } else {
                    if (debug) std::cout << "    -> 🗑️ [BOUNDARY-TRASH] Cut @ R:" << B_R << ", Q:" << B_Q  << " is fuzzy edge noise (< 50bp). Discarding completely!\n";
                    c_idx = end_idx + 1;
                    continue;
                }
            }

            if (is_isolated) {
                if (debug) std::cout << "    -> ⚡ [FAST-PATH] Isolated WGA Cut @ R:" << B_R << ", Q:" << B_Q << "\n";
                r_cuts.insert(B_R);
                q_cuts.insert(B_Q);
                
                std::string role = is_strict ? " [Forced STRICT]" : "";
                r_cut_reasons[B_R] = flat_cuts[start_idx].second.reasons + " (Isolated)" + role;
                q_cut_reasons[B_Q] = flat_cuts[start_idx].second.reasons + " (Isolated)" + role;
                active_cuts.push_back({B_R, B_Q});
                
                c_idx = end_idx + 1;
                continue; 
            }
        }

        int curr_r = check_r_start;
        int curr_q = best.inverse ? check_q_end : check_q_start;
        int q_dir  = best.inverse ? -1 : 1;
        
        int merged_consensus_len = 0;
        
        int qry_box_start = std::min(A_Q, E_Q);
        int qry_box_end   = std::max(A_Q, E_Q);
        
        // 🌟 優化：使用 Vector 取代極度緩慢的 std::map
        int r_box_len = E_R - A_R;
        int q_box_len = qry_box_end - qry_box_start;
        std::vector<int> r_map(std::max(1, r_box_len + 2), -1);
        std::vector<int> q_map(std::max(1, q_box_len + 2), -1);

        bool recording = false;
        if (curr_r == A_R && curr_q == A_Q) recording = true;

        for (auto op : best.CIGAR) {
            if (!recording && curr_r == E_R && curr_q == E_Q) break; 
            
            int len = op.first; char type = op.second;
            bool consumes_r = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumes_q = (type == 'M' || type == '=' || type == 'X' || type == 'I');
            
            for (int i = 0; i < len; ++i) {
                if (!recording && curr_r == A_R && curr_q == A_Q) {
                    recording = true;
                }
                if (recording && curr_r == E_R && curr_q == E_Q) {
                    recording = false;
                    break; 
                }
                
                if (recording) {
                    int r_idx = curr_r - A_R;
                    if (r_idx >= 0 && r_idx < r_map.size() && r_map[r_idx] == -1) r_map[r_idx] = merged_consensus_len;
                    
                    int q_idx = std::abs(curr_q - A_Q);
                    if (q_idx >= 0 && q_idx < q_map.size() && q_map[q_idx] == -1) q_map[q_idx] = merged_consensus_len;

                    merged_consensus_len++;
                }
                
                if (consumes_r) curr_r++;
                if (consumes_q) curr_q += q_dir;
            }
            
            if (recording && curr_r == E_R && curr_q == E_Q) {
                recording = false;
                break;
            }
        }
        
        if (E_R - A_R >= 0 && E_R - A_R < r_map.size()) r_map[E_R - A_R] = merged_consensus_len;
        if (std::abs(E_Q - A_Q) < q_map.size()) q_map[std::abs(E_Q - A_Q)] = merged_consensus_len;

        auto getMergedCoord = [&](int r, int q) {
            int r_idx = r - A_R;
            int q_idx = std::abs(q - A_Q);
            
            int r_val = (r_idx >= 0 && r_idx < r_map.size()) ? r_map[r_idx] : -1;
            int q_val = (q_idx >= 0 && q_idx < q_map.size()) ? q_map[q_idx] : -1;

            if (r_val != -1 && q_val != -1) return std::max(r_val, q_val); 
            if (r_val != -1) return r_val;
            if (q_val != -1) return q_val;
            
            for (int i = std::max(0, r_idx); i < r_map.size(); ++i) {
                if (r_map[i] != -1) return r_map[i];
            }
            return merged_consensus_len;
        };

        if (debug) {
            std::cout << "    -> 📦 Extracting and Merging A-E Sandbox:\n"
                      << "       ├─ Ref Bounds: [" << A_R << " -> " << E_R << "] (Len: " << (E_R - A_R) << ")\n"
                      << "       ├─ Qry Bounds: [" << qry_box_start << " -> " << qry_box_end << "] (Len: " << (qry_box_end - qry_box_start) << ")\n"
                      << "       └─ Merged Consensus Length: " << merged_consensus_len << "\n";
        }
        int local_B_merged = getMergedCoord(B_R, B_Q);
        int local_D_merged = getMergedCoord(D_R, D_Q);
        if (local_B_merged > local_D_merged) std::swap(local_B_merged, local_D_merged);

        int left_flank_len = local_B_merged;
        int right_flank_len = merged_consensus_len - local_D_merged;

        // 🌟 核心優化 (方案 A)：僅針對精確的 B-D 區域進行小沙盒 Extract 與 Merge (縮小 90% 以上計算規模)
        int bd_ref_start = std::min(B_R, D_R);
        int bd_ref_end   = std::max(B_R, D_R);
        int bd_qry_start = std::min(B_Q, D_Q);
        int bd_qry_end   = std::max(B_Q, D_Q);

        // 截取 B-D 區間內的 subCIGAR
        CigarString bd_subCIGAR;
        int bd_merged_len = 0;
        curr_r = check_r_start;
        curr_q = best.inverse ? check_q_end : check_q_start;
        bool bd_recording = false;

        for (auto op : best.CIGAR) {
            if (!bd_recording && curr_r == bd_ref_end && curr_q == (best.inverse ? bd_qry_start : bd_qry_end)) break;

            int len = op.first; char type = op.second;
            bool consumes_r = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumes_q = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            for (int i = 0; i < len; ++i) {
                if (!bd_recording && curr_r == bd_ref_start && curr_q == (best.inverse ? bd_qry_end : bd_qry_start)) {
                    bd_recording = true;
                }
                if (bd_recording && curr_r == bd_ref_end && curr_q == (best.inverse ? bd_qry_start : bd_qry_end)) {
                    bd_recording = false;
                    break;
                }

                if (bd_recording) {
                    if (!bd_subCIGAR.empty() && bd_subCIGAR.back().second == type) {
                        bd_subCIGAR.back().first++;
                    } else {
                        bd_subCIGAR.push_back({1, type});
                    }
                    bd_merged_len++;
                }

                if (consumes_r) curr_r++;
                if (consumes_q) curr_q += q_dir;
            }

            if (bd_recording && curr_r == bd_ref_end && curr_q == (best.inverse ? bd_qry_start : bd_qry_end)) {
                bd_recording = false;
                break;
            }
        }

        global_timer.start("ecc_2_sandbox_extract");
        // 🌟 輕量化沙盒 (Option 1): 直接由 bd_subCIGAR 與邊界構造極速評分 Block，省去全圖解構與重構開銷
        std::string dummy_cons(std::max(1, bd_merged_len), 'N');
        auto mergedLocalBlock = std::make_shared<Block>(999990, std::move(dummy_cons));

        int ref_total_len = bd_ref_end - bd_ref_start;
        int qry_total_len = bd_qry_end - bd_qry_start;

        Sequence refSeq("Ref");
        Sequence qrySeq("Qry");

        Segment refSeg(0, ref_total_len);
        Segment qrySeg(0, qry_total_len);

        int curr_pos = 0;
        for (auto op : bd_subCIGAR) {
            int len = op.first;
            char type = op.second;

            if (type == 'I') {
                refSeg.getVariants().push_back(Variant::createGap(curr_pos, curr_pos + len));
            } else if (type == 'D') {
                qrySeg.getVariants().push_back(Variant::createGap(curr_pos, curr_pos + len));
            }

            curr_pos += len;
        }

        refSeq.getSegments()[0] = std::move(refSeg);
        qrySeq.getSegments()[0] = std::move(qrySeg);

        mergedLocalBlock->getSequences()["Ref"] = std::move(refSeq);
        mergedLocalBlock->getSequences()["Qry"] = std::move(qrySeq);
        global_timer.stop("ecc_2_sandbox_extract");

        global_timer.start("ecc_2_sandbox_merge");
        // 輕量沙盒模式下合流已在 Extract 中完成
        global_timer.stop("ecc_2_sandbox_merge");

        global_timer.start("ecc_3_calc_split_scores");
        // 傳入左/右側 Flank 補償長度，評分邏輯 100% 完全相同
        std::vector<Block::ColumnSplitScore> merged_scores = mergedLocalBlock->calculateSplittingScores(0, mergedLocalBlock->getConsensus().length() - 1, left_flank_len, right_flank_len);
        global_timer.stop("ecc_3_calc_split_scores");

        global_timer.start("ecc_4_select_best_candidate");
        int best_k = -1;
        int best_cut = -1;
        int best_q_cut = -1;
        double max_score = -1e9; 
        
        for (int k = start_idx; k <= end_idx; ++k) {
            int x_global = flat_cuts[k].first.first;
            int qx_global = flat_cuts[k].first.second;

            int merged_x = getMergedCoord(x_global, qx_global);
            int score_idx = merged_x - local_B_merged;
            
            if (score_idx < 0 || score_idx >= merged_scores.size()) continue;

            auto s_merged = merged_scores[score_idx];
            double total_perfect = s_merged.perfect_bonus;
            double score_L = s_merged.id_left;
            double score_R = s_merged.id_right;

            int L_len_R = x_global - check_r_start;
            int R_len_R = check_r_end - x_global;
            int L_len_Q = std::abs(qx_global - check_q_start);
            int R_len_Q = std::abs(check_q_end - qx_global);

            if (std::max(L_len_R, L_len_Q) < L_min) score_L += 10000;
            if (std::max(R_len_R, R_len_Q) < L_min) score_R += 10000;

            double final_score = total_perfect - ((score_L + score_R)/2);

            if (debug) {
                std::cout << "      [CANDIDATE] R:" << x_global << " Q:" << qx_global 
                          << " -> Final: " << (final_score < -5000 ? "-INF" : std::to_string(final_score)) << "\n";
            }

            if (final_score > max_score) {
                max_score = final_score;
                best_k = k;
                best_cut = x_global;
                best_q_cut = qx_global;
            }
        }

        if (!strict_in_cluster.empty()) {
            if (debug) std::cout << "      ---> ⚠️ [STRICT OVERRIDE] STRICT points found in cluster. Forcing cuts to prevent topology collapse.\n";
            for (int s_idx : strict_in_cluster) {
                int cut_R = flat_cuts[s_idx].first.first;
                int cut_Q = flat_cuts[s_idx].first.second;

                r_cuts.insert(cut_R);
                q_cuts.insert(cut_Q);

                r_cut_reasons[cut_R] = flat_cuts[s_idx].second.reasons + " [Forced STRICT]";
                q_cut_reasons[cut_Q] = flat_cuts[s_idx].second.reasons + " [Forced STRICT]";
                active_cuts.push_back({cut_R, cut_Q});
            }
        } else if (best_k != -1) {
            if (debug) {
                std::cout << "      ---> 🎉 [WINNER] Snapped to R:" << best_cut << " Q:" << best_q_cut 
                          << " with Max Score: " << max_score << " (Best in Cluster)\n";
            }
            r_cuts.insert(best_cut);
            q_cuts.insert(best_q_cut);
            
            r_cut_reasons[best_cut] = flat_cuts[best_k].second.reasons + " [Best Score Winner]";
            q_cut_reasons[best_q_cut] = flat_cuts[best_k].second.reasons + " [Best Score Winner]";
            active_cuts.push_back({best_cut, best_q_cut});
        }
        global_timer.stop("ecc_4_select_best_candidate");
        c_idx = end_idx + 1;
    }
}

Alignments AlignmentCollection::getBestAlignments(BlockSet* refBlockSet, BlockSet* qryBlockSet, CoordinateManager& coordMgr, int L_min) {
    // bool debug = true; 
    bool debug = false; 
    
    std::vector<Alignment> final_results;

    // ⚠️ 必須用值拷貝而非引用！index-mode Consensus 在 buildCaches 過程中會
    //    重建原始 BlockSet 的 ancestral_seq_cache，導致引用懸空。
    const std::string& refSeq = refBlockSet->getAncestralSequence();
    const std::string& qrySeq = qryBlockSet->getAncestralSequence();

    // bool debug = (refBlockSet->getSequenceCount() >= 20 || qryBlockSet->getSequenceCount() >= 20);
  
    while (!queue_.empty()) {
        Alignment best = queue_.top();
        queue_.pop();

        int q_len = std::abs(best.qryIdx.second - best.qryIdx.first);
        int r_len = std::abs(best.refIdx.second - best.refIdx.first);

        if (best.alnLength < L_min || std::max(q_len, r_len) < L_min) continue;

        if (debug) {
            std::cout << "\n=================================================================\n";
            std::cout << "🎯 [NEW ALIGNMENT DEQUEUED] Score/Energy: " << best.energy 
                      << " | ID: " << best.ID << " | TYPE: " << (best.primary ? "PRIMARY" : "SECONDARY") << "\n"
                      << "   ├─ Ref Range: [" << best.refIdx.first << ", " << best.refIdx.second << ") (Len: " << r_len << ")\n"
                      << "   ├─ Qry Range: [" << std::min(best.qryIdx.first, best.qryIdx.second) << ", " 
                                             << std::max(best.qryIdx.first, best.qryIdx.second) << ") (Len: " << q_len << ")\n"
                      << "   ├─ Strand   : " << (best.inverse ? "Reverse (-)" : "Forward (+)") << "\n"
                      << "   ├─ Aln Length: " << best.alnLength << "\n"
                      << "   └─ Raw CIGAR: " << best.CIGAR << "\n";
        }

        // ========================================================
        // ✂️ STAGE 0 & PRE-TRIM: 處理 CIGAR 頭尾 Overlap 與邊界貼合 (NW)
        // ========================================================
        global_timer.start("gba_1_pre_trim_snap");
        preTrimAlignment(best, ref_coverageTracker, qry_coverageTracker, L_min);
        if (best.alnLength < L_min) {
            global_timer.stop("gba_1_pre_trim_snap");
            continue;
        }

        snapToBoundaries(best, coordMgr, refSeq, qrySeq, L_min, debug);
        best.updateAlnLength();
        global_timer.stop("gba_1_pre_trim_snap");

        if (best.alnLength < L_min) continue;

        // ========================================================
        // 🌟 STAGE 1 & 2：一維化雙軸切點池 (1D Flattened Bi-axial Pool)
        // ========================================================
        global_timer.start("gba_2_collect_cut_points");
        auto unified_cuts = collectCutPoints(best, coordMgr, ref_coverageTracker, qry_coverageTracker, L_min, debug);
        global_timer.stop("gba_2_collect_cut_points");
        
        // ========================================================
        // 🌟 STAGE 3：一維攤平掃描與微沙盒融合評分 (1D Flattened Array Sweep)
        // ========================================================
        if (debug) std::cout << "  [DEBUG-STAGE 3] Processing 1D Flattened Cut Array...\n";

        std::set<int> r_cuts, q_cuts;
        std::map<int, std::string> r_cut_reasons, q_cut_reasons;
        std::vector<std::pair<std::pair<int, int>, CutPoint>> flat_cuts(unified_cuts.begin(), unified_cuts.end());
        
        global_timer.start("gba_3_eval_cut_clusters");
        evaluateCutClusters(best, flat_cuts, refBlockSet, qryBlockSet, coordMgr, L_min, debug, r_cuts, q_cuts, r_cut_reasons, q_cut_reasons);
        global_timer.stop("gba_3_eval_cut_clusters");

        // ========================================================
        // 🌟 STAGE 4：執行切割
        // ========================================================
        global_timer.start("gba_4_split_alignment");
        std::vector<Alignment> frags;
        if (!r_cuts.empty() || !q_cuts.empty()) {
            frags = splitSingleAlignment(best, r_cuts, q_cuts); 
        } else {
            frags.push_back(best);
        }
        global_timer.stop("gba_4_split_alignment");

        // ========================================================
        // 🌟 STAGE 5：非對稱單軸剝離與微小碎片吸收
        // ========================================================
        global_timer.start("gba_5_absorb_fragments");
        absorbMicroFragments(frags, r_cut_reasons, q_cut_reasons, L_min, debug);
        global_timer.stop("gba_5_absorb_fragments");

        // ========================================================
        // 🌟 STAGE 6：最終過濾、Family 綁定與封裝輸出
        // ========================================================
        global_timer.start("gba_6_finalize_output");
        for (auto& frag : frags) {
            frag.finalizeAlignment(ref_coverageTracker, qry_coverageTracker, r_cut_reasons, q_cut_reasons);
            final_results.push_back(std::move(frag));
        }
        global_timer.stop("gba_6_finalize_output");

        if (!final_results.empty()) return final_results;
    }

    return final_results; 
}

// =====================================
// Alignment
// =====================================

void Alignment::updateAlnLength() {
    alnLength = getCigarLength(CIGAR);
}
