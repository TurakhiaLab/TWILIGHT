#pragma once

#include <map>
#include "type.hpp"
#include "segment.hpp"

class Sequence {
private:
  SequenceID ID;
  Segments segments;

public:
  Sequence(SequenceID id_) : ID(id_) {};
  Sequence() = default;
  ~Sequence() = default;

  // --- Getter ---
  SequenceID getId() const { return ID; };
  Segment &getSegment(int start_coordinate) {
    return segments[start_coordinate];
  };
  std::map<int, Segment> &getSegments() { return segments; };
  const std::map<int, Segment> &getSegments() const { return segments; };

  // --- Helper ---
  bool addSegment(int start, int end, Variants &var) {
    if (segments.find(start) != segments.end())
      return false;
    Segment seg(start, end, var);
    segments[start] = seg;
    return true;
  }
  bool addSegment(int start, int end) {
    if (segments.find(start) != segments.end())
      return false;
    Segment seg(start, end);
    segments[start] = seg;
    return true;
  }
  bool addSegment(Segment &seg) {
    if (segments.find(seg.getStart()) != segments.end())
      return false;
    segments[seg.getStart()] = seg;
    return true;
  }
  void reverseSegments(int len) {
    for (auto &seg : segments) {
      seg.second.reverseVariants(len);
    }
  }
};
