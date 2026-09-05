#ifndef TYPE_HPP
#define TYPE_HPP

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <memory>

// Forward declaration
class Variant;
class Segment;
class Sequence;
class Block;
class BlockSet;
class BlockManager;
class CoordinateManager;

struct Alignment;
struct BlockBoundary;
struct CoverageTracker;

namespace phylogeny {
class Tree;
class Node;
} // namespace phylogeny

// General
using String = std::string;
using Strings = std::vector<std::string>;
using StringPair = std::pair<std::string, std::string>;
using StringPairs = std::vector<StringPair>;

struct SequenceRef {
    std::string name;
    const std::string& seq;
};
using SequenceRefs = std::vector<SequenceRef>;
using StringMap = std::unordered_map<std::string, std::string>;
using Booleans = std::vector<bool>;
using Integers = std::vector<int>;

// Alignment related
using AlignmentId = int32_t;
using CigarOp = std::pair<int, char>;
using CigarString = std::vector<CigarOp>;
using Alignments = std::vector<Alignment>;
using Range = std::pair<int, int>;

// Variation
using Variants = std::vector<Variant>;
enum class VariantType {
  SNV,
  GAP,
};

// Segment
using Segments = std::map<int, Segment>;

// SequenceInfo
using SequenceID = std::string;
using Sequences = std::unordered_map<SequenceID, Sequence>;

// Block
using BlockID = uint32_t;
using BlockIDs = std::vector<BlockID>;
using BlockPtr = std::shared_ptr<Block>;
using BlockPtrPair = std::pair<BlockPtr, BlockPtr>;
using BlockWeakPtr = std::weak_ptr<Block>;
using BlockWeakPtrs = std::vector<std::weak_ptr<Block>>;
using Blocks = std::unordered_map<BlockID, BlockPtr>;
enum class BoundaryType { FLEXIBLE, STRICT, PUSH_RIGHT_ONLY, PUSH_LEFT_ONLY };
using BlockBoundaries = std::map<int, BlockBoundary>;
using VBlockID = std::pair<BlockID, int>;
using VBlockIDs = std::vector<VBlockID>;

// BlockSet
using BlockSetID = std::string;
using BlockSetPtr = std::shared_ptr<BlockSet>;
using BlockSets = std::unordered_map<BlockSetID, BlockSetPtr>;
using BlockSetPtrs = std::vector<BlockSet *>;

// BlockSet
using BlockManagerPtr = std::unique_ptr<BlockManager>;

// Phylogeny
using Tree = phylogeny::Tree;
using Node = phylogeny::Node;
using NodePair = std::pair<Node *, Node *>;
using NodePairs = std::vector<NodePair>;

// ============================================================================
// Centralized Constants for Alignment, Merging, and Graph Processing
// ============================================================================

// Snapping window size used during alignment processing (bp)
constexpr int SNAP_LENGTH = 100;

// Maximum sequence length for running dynamic programming (Needleman-Wunsch) alignment
constexpr int MAX_NW_LENGTH = 300;

// Minimum distance to the boundary to classify a cut point as non-noise (bp)
constexpr int MIN_DIST_2_BOUNDARY = 100;

// Default tree lookdown depth for finding distant/ancestral blocks
constexpr int DEFAULT_LOOKDOWN_DEPTH = 5;

// Default minimum alignment length threshold for merging (bp)
constexpr int DEFAULT_L_MIN = 100;

// Minimum consensus length for block realignment/filtering (bp)
constexpr int BLOCK_MIN_LENGTH = 100;

// Maximum consensus length for block realignment/filtering (bp)
constexpr int BLOCK_MAX_LENGTH = 50000;

// Minimum sequence identity threshold for block realignment/filtering (percentage)
constexpr float BLOCK_MIN_IDENTITY = 99.95f;

// Minimum sequence depth (coverage) threshold for block realignment/filtering
constexpr int BLOCK_MIN_DEPTH = 8;

// Maximum allowed gap length inside a segment during block splitting (bp)
constexpr int MAX_GAP_LEN = 50;

// Gap size threshold to classify a gap as a non-linear long gap (bp)
constexpr int LONG_GAP_THRESHOLD = 100;

// Minimum score threshold to accept secondary alignments during block merging
constexpr int MERGE_MIN_THRESHOLD = 100;

// Snapping threshold for core block boundary extraction (bp)
constexpr int SLICE_SNAP_THRESHOLD = 50;

// Cost weight multiplier for alignment penalty during block merging decision
constexpr int MERGE_ALPHA = 100;

#endif // TYPE_HPP