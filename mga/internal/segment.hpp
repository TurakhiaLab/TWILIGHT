#pragma once

#include <algorithm>
#include <utility>

#include "type.hpp"
#include "variant.hpp"

class Segment {
private:
  Range range;
  bool is_reverse;
  Variants variations;
  BlockWeakPtr prev_block;
  BlockWeakPtr next_block;
  int copy = 0;

public:
  Segment(int start, int end)
      : range({start, end}), is_reverse(false), variations() {};
  Segment(int start, int end, Variants &var)
      : range({start, end}), is_reverse(false), variations(var) {};
  Segment(int start, int end, bool reverse)
      : range({start, end}), is_reverse(reverse), variations() {};
  Segment(int start, int end, Variants &var, bool reverse)
      : range({start, end}), is_reverse(reverse), variations(var) {};
  Segment() = default;
  ~Segment() = default;
  Segment(const Segment &) = default;
  Segment &operator=(const Segment &) = default;

  // --- Getter ---
  Variants &getVariants() { return variations; }
  const Variants &getVariants() const { return variations; }
  Range getRange() const { return range; }
  int getStart() const { return range.first; }
  int getEnd() const { return range.second; }
  BlockWeakPtr getPrevBlock() const { return prev_block; }
  BlockWeakPtr getNextBlock() const { return next_block; }
  bool isReverse() const { return is_reverse; }
  int getCopyCount() const { return copy; }

  // --- Setter ---
  void setStart(int st) { range.first = st; }
  void setEnd(int en) { range.second = en; }
  void setReverse(bool reverse) { is_reverse = reverse; }
  void resetReverse() { is_reverse = false; }
  void setPrevBlock(BlockPtr block) { prev_block = block; }
  void setNextBlock(BlockPtr block) { next_block = block; }
  void setPrevBlock(BlockWeakPtr block) { prev_block = block; }
  void setNextBlock(BlockWeakPtr block) { next_block = block; }
  void setCopyCount(int copy_) { copy = copy_; }

  // --- Helper ---
  void reverseVariants(int consLen) {
    is_reverse = !is_reverse;
    for (auto &var : variations) {
      var.reverseComplement(consLen);
    }
    std::reverse(variations.begin(), variations.end());
  }
  std::pair<Segment, Segment> split(int localCut);
};
