#pragma once

#include "type.hpp"
#include "alignment.hpp"
#include "cigar_util.hpp"

class Variant {
private:
  VariantType type;
  Range range;
  char alt;

public:
  Variant(int pos_, char alt_)
      : type(VariantType::SNV), range({pos_, pos_ + 1}), alt(alt_) {};
  Variant(int start_, int end_)
      : type(VariantType::GAP), range({start_, end_}), alt('-') {};
  ~Variant() = default;

  // --- Getter ---
  VariantType getType() const { return type; }
  Range getRange() const { return range; }
  int getStart() const { return range.first; }
  int getEnd() const { return range.second; }
  char getAlt() const { return alt; }

  // --- Setter ---
  void setStart(int st_) { range.first = st_; }
  void setEnd(int en_) { range.second = en_; }
  void setRange(int start_, int end_) { range = {start_, end_}; }

  // --- Helper ---
  static Variant createGap(int start, int end) {
    Variant v(start, end);
    return v;
  }
  void reverseComplement(int consLen) {
    int oldStart = range.first;
    range.first = consLen - range.second;
    range.second = consLen - oldStart;
    if (type == VariantType::SNV)
      alt = complement(alt);
  }
  void shift(int offset) {
    range.first += offset;
    range.second += offset;
  }
};
