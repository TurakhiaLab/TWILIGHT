
#include "sequence.hpp"

// --- Segment ---

std::pair<Segment, Segment> Segment::split(int localCut) {
    Segment left(*this);
    Segment right(*this);
    
    left.variations.clear();
    right.variations.clear();
    
    int gap_bases_before_cut = 0;
    
    for (auto& var : variations) {
        if (var.getEnd() <= localCut) {
            left.variations.push_back(var);
            if (var.getType() == VariantType::GAP) {
                gap_bases_before_cut += (var.getEnd() - var.getStart());
            }
        } else if (var.getStart() >= localCut) {
            Variant right_var = var;
            int newStart = right_var.getStart() - localCut;
            int newEnd = right_var.getEnd() - localCut;
            if (newStart >= 0) { // 🌟 ROOT CAUSE FIX: 確保右半部 Variant 起點絕不為負數
                if (right_var.getType() == VariantType::GAP) {
                    right.variations.push_back(Variant::createGap(newStart, newEnd));
                } else {
                    right.variations.push_back(Variant(newStart, right_var.getAlt()));
                }
            }
        } else {
            // Straddling the cut point
            if (var.getType() == VariantType::GAP) {
                if (localCut > var.getStart()) {
                    Variant left_gap = Variant::createGap(var.getStart(), localCut);
                    left.variations.push_back(left_gap);
                    gap_bases_before_cut += (localCut - var.getStart());
                }
                if (var.getEnd() > localCut) {
                    int rStart = 0; // 右半部從 offset 0 開始
                    int rEnd = var.getEnd() - localCut;
                    if (rEnd > 0) {
                        right.variations.push_back(Variant::createGap(rStart, rEnd));
                    }
                }
            } else {
                // SNV at or straddling cut point
                if (var.getStart() < localCut) {
                    left.variations.push_back(var);
                } else {
                    int newStart = var.getStart() - localCut;
                    if (newStart >= 0) {
                        right.variations.push_back(Variant(newStart, var.getAlt()));
                    }
                }
            }
        }
    }
    
    int consumed_seq_len = localCut - gap_bases_before_cut;
    auto start_coordinate = this->getStart();
    auto end_coordinate = this->getEnd();
    
    if (!is_reverse) {
        left.setEnd(start_coordinate + consumed_seq_len);
        right.setStart(left.getEnd());
    } else {
        left.setStart(end_coordinate - consumed_seq_len);
        right.setEnd(left.getStart());
    }
    
    return {left, right};
}