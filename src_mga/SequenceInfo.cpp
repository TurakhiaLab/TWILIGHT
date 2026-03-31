#include "mga.hpp"


SequenceInfo::ID SequenceInfo::getID() const {
    return sequence_id;
}

SequenceInfo::Segment& SequenceInfo::getSegment(int start_coordinate) {
    return info[start_coordinate];
}

std::unordered_map<int, SequenceInfo::Segment>& SequenceInfo::getSegments() {
    return info;
}


bool SequenceInfo::addSegments(int start, int end, variationList& var) {
    if (info.find(start) != info.end()) return false;
    Segment seg(start, end, var);
    info[start] = seg;
    return true;
}


void SequenceInfo::Segment::reverseComplement(int consLen) {
    is_reverse = !is_reverse;
    for (auto& var : variations) {
        var.reverseComplement(consLen);
    }
    std::reverse(variations.begin(), variations.end());
}

mga::Cigar adjustCigarWithVariations(mga::Cigar& origCigar, SequenceInfo::Segment& refSeg, SequenceInfo::Segment& qrySeg) {
    mga::Cigar newCigar;
    
    // 將 Variations 整理成 Map，方便快速查找 (Key: consensus 座標, Value: Gap 長度)
    std::map<int, int> refGaps;
    for (auto& v : refSeg.getVariations()) {
        if (v.getType() == SequenceInfo::Variation::GAP) refGaps[v.getStart()] = v.getEnd() - v.getStart();
    }

    std::map<int, int> qryGaps;
    for (auto& v : qrySeg.getVariations()) {
        if (v.getType() == SequenceInfo::Variation::GAP) qryGaps[v.getStart()] = v.getEnd() - v.getStart();
    }

    int rConsPos = 0; // 目前在 Ref Consensus 上的座標
    int qConsPos = 0; // 目前在 Qry Consensus 上的座標

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
        int len = op.first;
        char type = op.second;

        // 在消耗原本的 CIGAR 之前，先檢查「目前的 Consensus 座標」是否有之前的 Merge 留下的 Gap
        // Ref Consensus 如果有 Gap，代表它變長了，對新的 Alignment 來說相當於一塊 Deletion ('D')
        while (refGaps.count(rConsPos)) {
            int gLen = refGaps[rConsPos];
            appendOp(gLen, 'D');
            rConsPos += gLen;
        }
        // Qry Consensus 如果有 Gap，代表它變長了，相當於一塊 Insertion ('I')
        while (qryGaps.count(qConsPos)) {
            int gLen = qryGaps[qConsPos];
            appendOp(gLen, 'I');
            qConsPos += gLen;
        }

        // 加上原本的 CIGAR 操作，並推進 Consensus 座標
        appendOp(len, type);
        if (type == 'M' || type == '=' || type == 'X') {
            rConsPos += len;
            qConsPos += len;
        } else if (type == 'D') {
            rConsPos += len;
        } else if (type == 'I') {
            qConsPos += len;
        }
    }

    // 處理尾巴可能殘留的 Gap
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