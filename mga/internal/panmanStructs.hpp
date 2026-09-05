#pragma once

#include <string>
#include <vector>
#include <map>
#include <cstdint>

namespace mga {
namespace panman {

struct NucMutData {
    int32_t primaryBlockId = 0;
    int32_t nucPosition = 0;
    int32_t nucGapPosition = -1;
    bool nucGapExist = false;
    uint32_t mutInfo = 0; // Packed wire format
};

struct BlockMutData {
    int32_t primaryBlockId = 0;
    int32_t secondaryBlockId = -1;
    int blockMutInfo = 0; // 0: Insert, 1: Delete, 2: Present/No-change
    bool inversion = false;
    int64_t chrIdx = 0;
};

struct NodeData {
    std::string identifier;
    std::vector<NucMutData> nucMutations;
    std::vector<BlockMutData> blockMutations;
    std::vector<std::string> annotations;
    std::vector<NodeData*> children;
};

struct GapListData {
    int32_t primaryBlockId = 0;
    int32_t secondaryBlockId = -1;
    std::vector<int32_t> nucPosition;
    std::vector<int32_t> nucGapLength;
};

struct ConsensusBlockData {
    std::vector<uint32_t> consensusSeq; // Packed 4-bit DNA symbols
    std::vector<std::pair<int64_t, bool>> blockIds; // (blockId, blockGapExist)
};

} // namespace panman
} // namespace mga
