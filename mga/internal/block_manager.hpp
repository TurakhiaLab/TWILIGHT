#pragma once

#include <unordered_map>
#include <string>
#include <memory>
#include <iostream>
#include <vector>

#include "type.hpp"
#include "block_set.hpp"
#include "alignment.hpp"
#include "option.hpp"

class BlockManager {
private:
  BlockSets blocksets;
  std::unordered_map<std::string, int> sequence_lengths;
  std::unordered_map<std::string, std::string> sequences;

public:
  BlockManager() = default;
  ~BlockManager() = default;

  BlockManager(const BlockManager &) = delete;
  BlockManager &operator=(const BlockManager &) = delete;

  // --- Getter ---
  BlockSet *getBlockSet(BlockSetID id) const {
    return (blocksets.find(id) == blocksets.end()) ? nullptr
                                                   : blocksets.at(id).get();
  };
  BlockSets getBlockSets() const { return blocksets; };
  BlockSetPtrs getAllBlockSetPtrs() const {
    std::vector<BlockSet *> all_sets;
    all_sets.reserve(blocksets.size());
    for (const auto &pair : blocksets)
      all_sets.push_back(pair.second.get());
    return all_sets;
  }

  BlockSet *createBlockSet(BlockSetID id) {
    auto new_set = std::make_unique<BlockSet>(id);
    auto ptr = new_set.get();
    blocksets[id] = std::move(new_set);
    return ptr;
  }
  bool changeBlockSetId(BlockSetID old_id, BlockSetID new_id) {
    if (blocksets.count(new_id) || old_id == new_id)
      return false;
    auto node_handle = blocksets.extract(old_id);
    if (node_handle.empty())
      return false;
    node_handle.key() = new_id;
    node_handle.mapped()->ID = new_id;
    blocksets.insert(std::move(node_handle));
    return true;
  }
  bool removeBlockSet(BlockSetID id) {
    auto it = blocksets.find(id);
    if (it != blocksets.end()) {
      blocksets.erase(it);
      return true;
    }
    return false;
  }

  void print(std::ostream &os = std::cout) const;
  BlockSet *merge(BlockSet *refSet, BlockSet *qrySet,
                  AlignmentCollection &alnCollection, BlockSetID newID,
                  Tree *tree = nullptr, int L_min = DEFAULT_L_MIN);
  void orientCircularGenomes(Option &option, std::string refSequenceName = "");

  void updateLongestSequences();
  void addSequenceLength(std::string seqName, int seqLen) {
    sequence_lengths[seqName] = seqLen;
  }
  int getSequenceLength(std::string seqName) {
    return sequence_lengths[seqName];
  }
  void addSequence(std::string &seqName, std::string &seq) {
    if (sequences.find(seqName) != sequences.end()) {
      std::cerr << "ERROR: Sequence " << seqName << " already exists.\n";
      return;
    }
    sequences[seqName] = seq;
  };
  std::string &getSequence(std::string seqName) { return sequences[seqName]; }
};
