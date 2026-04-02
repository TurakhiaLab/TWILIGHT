#include "mga.hpp"


SequenceInfo::ID SequenceInfo::getID() const {
    return sequence_id;
}

Segment& SequenceInfo::getSegment(int start_coordinate) {
    return info[start_coordinate];
}

std::map<int, Segment>& SequenceInfo::getSegments() {
    return info;
}


bool SequenceInfo::addSegments(int start, int end, variationList& var) {
    if (info.find(start) != info.end()) return false;
    Segment seg(start, end, var);
    info[start] = seg;
    return true;
}

bool SequenceInfo::addSegments(int start, int end) {
    if (info.find(start) != info.end()) return false;
    Segment seg(start, end);
    info[start] = seg;
    return true;
}

bool SequenceInfo::addSegments(Segment& seg) {
    if (info.find(seg.getStart()) != info.end()) return false;
    info[seg.getStart()] = seg;
    return true;
}


void Segment::reverseComplement(int consLen) {
    is_reverse = !is_reverse;
    for (auto& var : variations) {
        var.reverseComplement(consLen);
    }
    std::reverse(variations.begin(), variations.end());
}

std::pair<Segment, Segment> Segment::split(int localCut) {
    Segment left(*this);
    Segment right(*this);
    
    left.variations.clear();
    right.variations.clear();
    
    int gap_bases_before_cut = 0;
    
    for (auto& var : variations) {
        if (var.getEnd() <= localCut) {
            left.variations.push_back(var);
            if (var.getType() == Variation::GAP) {
                gap_bases_before_cut += (var.getEnd() - var.getStart());
            }
        } else if (var.getStart() >= localCut) {
            Variation right_var = var;
            right_var.shift(-localCut); 
            right.variations.push_back(right_var);
        } else {
            // Gap 橫跨切點
            Variation left_gap = Variation::createGap(var.getStart(), localCut);
            left.variations.push_back(left_gap);
            gap_bases_before_cut += (localCut - var.getStart());
            
            Variation right_gap = Variation::createGap(localCut, var.getEnd());
            right_gap.shift(-localCut);
            right.variations.push_back(right_gap);
        }
    }
    
    int consumed_seq_len = localCut - gap_bases_before_cut;
    
    if (!is_reverse) {
        left.end_coordinate = start_coordinate + consumed_seq_len;
        right.start_coordinate = left.end_coordinate;
    } else {
        left.start_coordinate = end_coordinate - consumed_seq_len;
        right.end_coordinate = left.start_coordinate;
    }
    
    return {left, right};
}

mga::Cigar adjustCigarWithVariations(mga::Cigar& origCigar, Segment& refSeg, Segment& qrySeg, bool qryInverse, int qryConsLen) {
    mga::Cigar newCigar;
    
    std::map<int, int> refGaps;
    for (auto& v : refSeg.getVariations()) {
        if (v.getType() == Variation::GAP) {
            refGaps[v.getStart()] = v.getEnd() - v.getStart();
        }
    }

    std::map<int, int> qryGaps;
    for (auto& v : qrySeg.getVariations()) {
        if (v.getType() == Variation::GAP) {
            int s = v.getStart();
            int e = v.getEnd();
            // 【核心修復】：在這裡精準計算反轉座標，絕不依賴外部函式
            if (qryInverse) {
                int old_s = s;
                s = qryConsLen - e;
                e = qryConsLen - old_s;
            }
            qryGaps[s] = e - s;
        }
    }

    int rConsPos = 0; 
    int qConsPos = 0; 

    auto appendOp = [&](int len, char op) {
        if (len > 0) {
            if (!newCigar.empty() && newCigar.back().second == op) {
                newCigar.back().first += len;
            } else {
                newCigar.push_back({len, op});
            }
        }
    };

    for (const auto& op : origCigar) {
        int remain = op.first;
        char type = op.second;

        while (remain > 0) {
            int nextRefGapDist = remain + 1;
            int nextQryGapDist = remain + 1;
            
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            if (consumesRef) {
                auto it = refGaps.lower_bound(rConsPos);
                if (it != refGaps.end() && it->first < rConsPos + remain) {
                    nextRefGapDist = it->first - rConsPos;
                }
            }
            if (consumesQry) {
                auto it = qryGaps.lower_bound(qConsPos);
                if (it != qryGaps.end() && it->first < qConsPos + remain) {
                    nextQryGapDist = it->first - qConsPos;
                }
            }

            int step = std::min({remain, nextRefGapDist, nextQryGapDist});

            if (step > 0) {
                appendOp(step, type);
                if (consumesRef) rConsPos += step;
                if (consumesQry) qConsPos += step;
                remain -= step;
            }

            if (consumesRef && refGaps.count(rConsPos)) {
                int gLen = refGaps[rConsPos];
                appendOp(gLen, 'D'); 
                rConsPos += gLen;
            }
            if (consumesQry && qryGaps.count(qConsPos)) {
                int gLen = qryGaps[qConsPos];
                appendOp(gLen, 'I'); 
                qConsPos += gLen;
            }
        }
    }

    while (refGaps.count(rConsPos)) {
        int gLen = refGaps[rConsPos];
        appendOp(gLen, 'D');
        rConsPos += gLen;
    }
    while (qryGaps.count(qConsPos)) {
        int gLen = qryGaps[qConsPos];
        appendOp(gLen, 'I');
        qConsPos += gLen;
    }

    return newCigar;
}

mga::Cigar extractSubCigar(const mga::Cigar& origCigar, int refOffset, int refLen) {
    mga::Cigar subCigar;
    int currentRef = 0;     // 目前走到了原本 Reference 的哪個位置
    int extractedRef = 0;   // 已經提取了多少 Reference 的長度

    // 輔助 Lambda：安全地將操作加入 subCigar，並自動合併連續的相同操作
    auto appendOp = [&](int len, char type) {
        if (len <= 0) return;
        if (!subCigar.empty() && subCigar.back().second == type) {
            subCigar.back().first += len;
        } else {
            subCigar.push_back({len, type});
        }
    };

    for (const auto& op : origCigar) {
        if (extractedRef >= refLen) break; // 已經提取足夠長度，提早結束

        int opLen = op.first;
        char type = op.second;

        bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
        int refOpLen = consumesRef ? opLen : 0;

        // 【過濾階段】：跳過我們不需要的區段
        if (consumesRef) {
            // 如果這個操作完全落在 refOffset 之前 (或剛好切齊)
            if (currentRef + refOpLen <= refOffset) {
                currentRef += refOpLen;
                continue;
            }
        } else {
            // 如果是 Insertion (不消耗 Ref)，且發生在我們開始提取的點「之前」，就跳過。
            // (注意：用嚴格小於 '<'，代表如果 Insertion 剛好發生在 refOffset 切點上，我們將其保留給右邊的 fragment！)
            if (currentRef < refOffset) {
                continue;
            }
        }

        // 【提取階段】：這個操作涵蓋到了我們要的範圍
        int useLen = opLen;

        // 1. 如果操作橫跨了起點，修剪掉前面多餘的部分
        if (consumesRef && currentRef < refOffset) {
            int trim = refOffset - currentRef;
            useLen -= trim;
            currentRef += trim;
        }

        // 2. 如果操作橫跨了終點，修剪掉後面多餘的部分
        if (consumesRef && (extractedRef + useLen > refLen)) {
            useLen = refLen - extractedRef;
        }

        // 加入 CIGAR
        appendOp(useLen, type);

        // 推進座標
        if (consumesRef) {
            currentRef += useLen;
            extractedRef += useLen;
        }
    }
    
    return subCigar;
}

mga::Cigar adjustCigarWithVariations(const mga::Cigar& origCigar, Segment& refSeg, Segment& qrySeg) {
    mga::Cigar newCigar;
    
    // 1. 將 Variations 整理成 Map，方便快速查找 (Key: Consensus上的起點, Value: Gap長度)
    std::map<int, int> refGaps;
    for (auto& v : refSeg.getVariations()) {
        if (v.getType() == Variation::GAP) {
            refGaps[v.getStart()] = v.getEnd() - v.getStart();
        }
    }
    std::map<int, int> qryGaps;
    for (auto& v : qrySeg.getVariations()) {
        if (v.getType() == Variation::GAP) {
            qryGaps[v.getStart()] = v.getEnd() - v.getStart();
        }
    }

    int rConsPos = 0; // 追蹤目前走到 Ref Consensus 的哪個位置
    int qConsPos = 0; // 追蹤目前走到 Qry Consensus 的哪個位置

    auto appendOp = [&](int len, char op) {
        if (len > 0) {
            if (!newCigar.empty() && newCigar.back().second == op) {
                newCigar.back().first += len;
            } else {
                newCigar.push_back({len, op});
            }
        }
    };

    // 2. 逐一消耗原始 CIGAR
    for (const auto& op : origCigar) {
        int remainingLen = op.first;
        char type = op.second;

        // 如果這個操作很長，我們可能需要因為中途的 GAP 把它分段處理
        while (remainingLen > 0) {
            
            // A. 先處理堆積在當前 Consensus 座標上的歷史 GAP
            // Ref 的 GAP 代表 RefConsensus 變長了，所以在這裡補上 Deletion
            while (refGaps.count(rConsPos)) {
                int gLen = refGaps[rConsPos];
                appendOp(gLen, 'D');
                rConsPos += gLen;
            }
            // Qry 的 GAP 代表 QryConsensus 變長了，所以在這裡補上 Insertion
            while (qryGaps.count(qConsPos)) {
                int gLen = qryGaps[qConsPos];
                appendOp(gLen, 'I');
                qConsPos += gLen;
            }

            // B. 計算我們這一步最多可以走多遠，才不會撞到下一個 GAP
            int nextRefGapDist = remainingLen;
            int nextQryGapDist = remainingLen;
            
            if (type == 'M' || type == '=' || type == 'X') {
                auto rIt = refGaps.upper_bound(rConsPos);
                if (rIt != refGaps.end() && rIt->first < rConsPos + remainingLen) {
                    nextRefGapDist = rIt->first - rConsPos;
                }
                auto qIt = qryGaps.upper_bound(qConsPos);
                if (qIt != qryGaps.end() && qIt->first < qConsPos + remainingLen) {
                    nextQryGapDist = qIt->first - qConsPos;
                }
            } else if (type == 'D') {
                auto rIt = refGaps.upper_bound(rConsPos);
                if (rIt != refGaps.end() && rIt->first < rConsPos + remainingLen) {
                    nextRefGapDist = rIt->first - rConsPos;
                }
            } else if (type == 'I') {
                auto qIt = qryGaps.upper_bound(qConsPos);
                if (qIt != qryGaps.end() && qIt->first < qConsPos + remainingLen) {
                    nextQryGapDist = qIt->first - qConsPos;
                }
            }

            // C. 取安全距離前進
            int stepLen = std::min({remainingLen, nextRefGapDist, nextQryGapDist});
            
            appendOp(stepLen, type);
            
            // D. 推進 Consensus 座標
            if (type == 'M' || type == '=' || type == 'X') {
                rConsPos += stepLen;
                qConsPos += stepLen;
            } else if (type == 'D') {
                rConsPos += stepLen;
            } else if (type == 'I') {
                qConsPos += stepLen;
            }
            
            remainingLen -= stepLen;
        }
    }

    // 3. 收尾：如果整個原始序列都走完了，但尾巴還有殘留的 GAP，全部補齊
    while (refGaps.count(rConsPos)) {
        int gLen = refGaps[rConsPos];
        appendOp(gLen, 'D');
        rConsPos += gLen;
    }
    while (qryGaps.count(qConsPos)) {
        int gLen = qryGaps[qConsPos];
        appendOp(gLen, 'I');
        qConsPos += gLen;
    }

    return newCigar;
}