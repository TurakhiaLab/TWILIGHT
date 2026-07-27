
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
            right_var.shift(-localCut); 
            right.variations.push_back(right_var);
        } else {
            // Gap across the cut point
            Variant left_gap = Variant::createGap(var.getStart(), localCut);
            left.variations.push_back(left_gap);
            gap_bases_before_cut += (localCut - var.getStart());
            
            Variant right_gap = Variant::createGap(localCut, var.getEnd());
            right_gap.shift(-localCut);
            right.variations.push_back(right_gap);
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