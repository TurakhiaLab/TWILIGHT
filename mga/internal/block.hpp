#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <utility>

#include "type.hpp"
#include "sequence.hpp"
#include "alignment.hpp"
#include "consensus.hpp"

// Forward declaration
class BlockSet;

std::shared_ptr<Block> extractBlockFromSuper(BlockSet *bSet,
                                             std::shared_ptr<Block> superBlock,
                                             int start, int end);

class Block {
private:
  BlockID ID;
  Consensus consensus;
  Sequences sequences;
  // BlockWeakPtr prev_block;
  // BlockWeakPtr next_block;
  BlockWeakPtrs prev_blocks;
  BlockWeakPtrs next_blocks;

  Booleans distant;

public:
  Block(BlockID id, Consensus consensus)
      : ID(id), consensus(std::move(consensus)),
        distant(Booleans(1, false)) {};

  Block(BlockID id, std::string consensus_str)
      : ID(id), consensus(Consensus(std::move(consensus_str))),
        distant(Booleans(1, false)) {};

  ~Block() = default;
  Block(const Block &) = delete;
  Block &operator=(const Block &) = delete;

  // --- Getter ---
  BlockID getId() const { return ID; }
  const Consensus &getConsensus() const { return consensus; }
  Consensus &getConsensus() { return consensus; }
  std::string getConsensusString() const { return consensus.getConsensusString(); }
  std::string getConsensusAsFasta(const std::string &name) const {
    return (">" + name + "\n" + consensus.getConsensusString() + "\n");
  }
  Sequences &getSequences() { return sequences; }
  BlockWeakPtr getPrevBlock(int copy) const {
    if (copy >= 0 && copy < static_cast<int>(prev_blocks.size()))
      return prev_blocks[copy];
    return BlockWeakPtr();
  };
  BlockWeakPtr getNextBlock(int copy) const {
    if (copy >= 0 && copy < static_cast<int>(next_blocks.size()))
      return next_blocks[copy];
    return BlockWeakPtr();
  };
  bool isCoreBlock() const { return distant.size() == 1; }
  bool isAllDistant() const {
    if (distant.empty()) return false;
    for (auto d : distant) {
      if (!d) return false;
    }
    return true;
  }
  bool isDistant(int copy = -1) const {
    if (copy != -1) {
      if (copy >= 0 && copy < static_cast<int>(distant.size()))
        return distant[copy];
      return false;
    }
    return isAllDistant();
  }
  int getMaxCopy() const;

  // --- Setter ---
  void setId(BlockID id_) { ID = id_; }
  void setConsensus(const Consensus &consensus_) { consensus = consensus_; }
  void setConsensus(Consensus &&consensus_) { consensus = std::move(consensus_); }
  void setConsensus(const std::string &consensus_str) { consensus = Consensus(consensus_str); }
  void setSequences(Sequences &sequences_) {
    sequences.clear();
    sequences = sequences_;
  }
  void setPrevBlock(int copy, BlockPtr block_ptr) {
    if (copy >= (int)prev_blocks.size()) {
      prev_blocks.resize(copy + 1);
    }
    prev_blocks[copy] = block_ptr;
  }
  void setNextBlock(int copy, BlockPtr block_ptr) {
    if (copy >= (int)next_blocks.size()) {
      next_blocks.resize(copy + 1);
    }
    next_blocks[copy] = block_ptr;
  }
  void setDistant(bool dis, int copy) {
    if (copy >= (int)distant.size()) {
      distant.resize(copy + 1, false);
    }
    distant[copy] = dis;
  }
  void clearPrevBlocksMap() { prev_blocks.clear(); }
  void clearNextBlocksMap() { next_blocks.clear(); }

  // --- Helper ---
  void addSequence(const Sequence &seq) { sequences[seq.getId()] = seq; }
  void reverse() {
    int len = consensus.length();
    consensus.setStoredString(getReverseComplement(consensus.getConsensusString()));
    consensus.setUseStoredString(true);
    for (auto &seq : sequences) {
      seq.second.reverseSegments(len);
    }
  }
  BlockPtrPair split(int cut, int copy = -1) const;
  void applyCopyAssignment(int target_copy, int shift_amount);

  struct ColumnSplitScore {
    double perfect_bonus;
    double id_left;
    double id_right;
  };

  std::vector<ColumnSplitScore> calculateSplittingScores(int local_start,
                                                         int local_end,
                                                         int left_flank_len = 0,
                                                         int right_flank_len = 0,
                                                         bool debug = false);

  void print(std::ostream &os = std::cout, int copy = -1) const;
  bool normalizeStrand();
  void refine();
  void refineConsensusAndVariants();

  // std::pair<std::shared_ptr<Block>, std::shared_ptr<Block>> split(int offset,
  // ID new_id_1, ID new_id_2);
};