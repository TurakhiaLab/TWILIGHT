#pragma once

#include <iostream>
#include <string>

#include "type.hpp"

struct BlockBoundary {
  BlockID leftBlockId;
  BlockID rightBlockId;

  int leftConsensusEndPos;

  BoundaryType type = BoundaryType::FLEXIBLE;

  std::string reason; // For debugging

  void print() const {
    std::cout << "[Boundary] " << leftBlockId;
    std::cout << " -> " << rightBlockId;
    std::cout << " | AbsPos: " << leftConsensusEndPos << " | Type: ";
    switch (type) {
    case BoundaryType::FLEXIBLE:
      std::cout << "FLEXIBLE";
      break;
    case BoundaryType::STRICT:
      std::cout << "STRICT";
      break;
    case BoundaryType::PUSH_RIGHT_ONLY:
      std::cout << "PUSH_RIGHT_ONLY";
      break;
    case BoundaryType::PUSH_LEFT_ONLY:
      std::cout << "PUSH_LEFT_ONLY";
      break;
    }
    std::cout << " | Reason: " << reason << "\n";
  }
};
