#pragma once

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

#include "type.hpp"
#include "block.hpp"

// Forward declaration
class BlockManager;
class CoordinateManager;
struct Option;
struct Alignment;
using Alignments = std::vector<Alignment>;

class BlockSet {

public:
  struct SegNode; // Forward declaration

private:
  friend class BlockManager; // Allow BlockManager to modify ID

  BlockSetID ID;
  Blocks blocks;

  std::unordered_set<std::string> sequence_names;

  bool is_cached;
  VBlockIDs linear_block_cache;
  VBlockIDs ancestral_block_cache;
  std::string ancestral_seq_cache;
  std::map<int, VBlockID> ancestral_offset_cache;

  std::atomic<BlockID> next_block_id_{1};

  Tree *tree_ptr;

  // --- Merger Helper ---
  struct MergeMappingData {
    Consensus mergedCons;
    std::vector<int> refOldToNew;
    std::vector<int> qryOldToNew;
    Variants newRefGaps;
    Variants newQryGaps;
    CigarString refConsensusChanges;
    CigarString qryConsensusChanges;
  };
  MergeMappingData calculateMergedConsensusAndMappings(
      std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock,
      const CigarString &cigar, bool inverse,
      const std::vector<Segment *> &refSegsFlat,
      const std::vector<Segment *> &qrySegsFlat);
  void updateSegmentVariations(std::vector<Segment *> &refSegsFlat,
                               std::vector<Segment *> &qrySegsFlat,
                               const MergeMappingData &mapData, bool inverse,
                               int refConsLen, int qryConsLen);

  void debugValidateBlockMerge(std::shared_ptr<Block> refBlock,
                               std::shared_ptr<Block> qryBlock);

public:
  struct SegNode {
    int start;
    int end;
    BlockID blkId;
  };

  BlockSet(BlockSetID id) : ID(id), is_cached(false), tree_ptr(nullptr) {};
  ~BlockSet() = default;

  // --- Getter ---
  BlockSetID getId() const { return ID; }
  BlockPtr getBlock(BlockID id) {
    return (blocks.find(id) == blocks.end()) ? nullptr : blocks[id];
  };
  std::unordered_set<std::string> getSequences() { return sequence_names; }
  const VBlockIDs &getLinearizeBlocks();
  const VBlockIDs &getAncestralBlocks();
  const std::string &getAncestralSequence();
  const std::map<int, VBlockID> &getAncestralBlocksOffsets();

  BlockWeakPtrs getAllBlocks() {
    BlockWeakPtrs all_blocks;
    all_blocks.reserve(blocks.size());
    for (const auto &pair : blocks)
      all_blocks.push_back(pair.second);
    return all_blocks;
  };
  size_t getSequenceCount() { return sequence_names.size(); }

  // --- Setter ---
  void invalidateRepCache() { is_cached = false; }
  void setTree(Tree *tree) { tree_ptr = tree; }
  Tree *getTree() const { return tree_ptr; }
  void setDistantBlocks(Tree &tree, int lookdownDepth = DEFAULT_LOOKDOWN_DEPTH);

  // --- Helper ---
  void print(std::ostream &os) const;
  void buildCaches();
  void rebuildAllPointers();
  bool detectVBlockCycle(const CoordinateManager &coordMgr, bool verbose = false) const;
  void rebuildLinearGraph(const CoordinateManager &coordMgr);
  BlockPtr createBlockWithId(BlockID id, Consensus consensus) {
    auto new_block = std::make_shared<Block>(id, std::move(consensus));
    blocks[id] = new_block;
    if (id >= next_block_id_) {
      next_block_id_ = id + 1;
    }
    invalidateRepCache();
    return new_block;
  }
  BlockPtr createBlock(Consensus consensus) {
    BlockID new_id = next_block_id_++;
    return createBlockWithId(new_id, std::move(consensus));
  }
  BlockPtr createBlock(const std::string &consensus) {
    return createBlock(Consensus(consensus));
  }
  BlockSet* createIndexedCopy(BlockManager* manager, const BlockSetID& newID);
  bool deleteBlock(BlockID id) {
    invalidateRepCache();
    return blocks.erase(id) > 0;
  }
  void clearBlocks() {
    blocks.clear();
    invalidateRepCache();
    next_block_id_ = 1;
  }
  void addSequenceName(std::string seqName) { sequence_names.insert(seqName); }
  BlockPtr addBlock(BlockPtr oldBlock);
  BlockPtr addBlockReferencing(BlockPtr oldBlock, BlockSet* sourceSet);
  BlockPtr extractBlock(int extract_start, int extract_end);

  std::string reconstructSequence(const std::string &seqName);
  BlockPtr concatAllBlocks(BlockID superId);
  BlockPtr concatBlocks(const BlockIDs &target_blocks, bool dryRun = false);

  std::pair<BlockID, BlockID> splitSingleBlock(int parentID, int localCut);

  // --- Debug/Validate ---
  void debugValidateSegments(bool verbose = false);
  void debugValidateLinkages(bool verbose = false);
  void debugValidateQuality(bool verbose = false);
  void debugValidateSequences(BlockManager *manager, bool verbose = false);
  void debugValidateLinearizedBlocks(bool verbose = false);
  void debugValidateBlocks(bool verbose = false);

  // --- Merger ---
  BlockPtr mergeTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock,
                          const CigarString &cigar, bool inverse, int mode,
                          int r_probe_start = 0, int q_probe_start = 0,
                          int r_probe_copy = 0, int q_probe_copy = 0);
  void linkTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock,
                     const CigarString &cigar, bool inverse);

  // --- refine ---
  void refineGraph();
  void splitBlocksByLongGaps();
  void extractMicroSegments();
  BlockIDs absorbMicroBlocksIndividual();
  void cleanGappedColumns();
  BlockIDs splitMultiBlocks(BlockID parentID, const std::vector<int> &cuts);
  void reconnectBlocks();

  void realignBlocks(std::string tempDir);
  bool realignBlock(BlockID blkId, std::string tempDir, int iterations = 1);

  // --- Self alignment ---
  Alignments selfAlign(Option &option, CoordinateManager *coordMgr = nullptr, double min_ratio = 0.10);
  Alignments selfAlignDistant(Option &option, CoordinateManager *coordMgr = nullptr, double min_ratio = 0.10);
};
