#ifndef TYPE_HPP
#define TYPE_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include <memory>

// Forward declaration
class Variant;
class Segment;
class Sequence;
class Block;
class BlockSet;
class BlockManager;
struct Alignment;
struct BlockBoundary;

namespace phylogeny {
    class Tree;
    class Node;
}
    



// General
using StringPair = std::pair<std::string, std::string>;
using StringPairs = std::vector<StringPair>;
using StringMap = std::unordered_map<std::string, std::string>;


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
using FamilyID = uint32_t;
using BlockIDs = std::vector<BlockID>;
using BlockPtr = std::shared_ptr<Block>;
using BlockWeakPtr = std::weak_ptr<Block>;
using BlockWeakPtrs = std::vector<std::weak_ptr<Block>>;
using Blocks = std::unordered_map<BlockID, BlockPtr>;
enum class BoundaryType {
    FLEXIBLE,        
    STRICT,          
    PUSH_RIGHT_ONLY, 
    PUSH_LEFT_ONLY   
};
using BlockBoundaries = std::map<int, BlockBoundary>;



// BlockSet
using BlockSetID = std::string;
using BlockSetPtr = std::shared_ptr<BlockSet>;
using BlockSets = std::unordered_map<BlockSetID, BlockSetPtr>;
using BlockSetPtrs = std::vector<BlockSet*>;

// BlockSet
using BlockManagerPtr = std::unique_ptr<BlockManager>;

// Phylogeny
using Tree = phylogeny::Tree;
using Node = phylogeny::Node;
using NodePair = std::pair<Node*, Node*>;
using NodePairs = std::vector<NodePair>;



#endif // TYPE_HPP