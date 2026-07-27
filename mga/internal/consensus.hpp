#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <iostream>

#include "type.hpp"

// Forward declaration to break circular dependency with block_set.hpp
class BlockSet;

class Consensus {
private:
  mutable bool use_stored_string;
  mutable std::string stored_string;

  mutable BlockSet* block_set_ptr;
  mutable int start_idx;
  mutable int end_idx;
  mutable std::map<int, std::string> insertions;

  std::string buildStringFromIndex() const;

public:
  Consensus() 
      : use_stored_string(true), 
        stored_string(""), 
        block_set_ptr(nullptr), 
        start_idx(-1), 
        end_idx(-1) {}

  Consensus(const std::string& str) 
      : use_stored_string(true), 
        stored_string(str), 
        block_set_ptr(nullptr), 
        start_idx(-1), 
        end_idx(-1) {}

  Consensus(std::string&& str) noexcept
      : use_stored_string(true), 
        stored_string(std::move(str)), 
        block_set_ptr(nullptr), 
        start_idx(-1), 
        end_idx(-1) {}

  Consensus(BlockSet* block_set, int start, int end, const std::map<int, std::string>& insts = {}) 
      : use_stored_string(false), 
        stored_string(""), 
        block_set_ptr(block_set), 
        start_idx(start), 
        end_idx(end), 
        insertions(insts) {}

  Consensus(BlockSet* block_set, int start, int end, std::map<int, std::string>&& insts) noexcept
      : use_stored_string(false), 
        stored_string(""), 
        block_set_ptr(block_set), 
        start_idx(start), 
        end_idx(end), 
        insertions(std::move(insts)) {}

  ~Consensus() = default;

  // --- Getters ---
  bool isUsingStoredString() const { return use_stored_string; }
  const std::string& getStoredString() const { return stored_string; }
  BlockSet* getBlockSetPtr() const { return block_set_ptr; }
  int getStartIdx() const { return start_idx; }
  int getEndIdx() const { return end_idx; }
  const std::map<int, std::string>& getInsertions() const { return insertions; }

  // --- Setters ---
  void setUseStoredString(bool use_stored) { use_stored_string = use_stored; }
  void setStoredString(const std::string& str) { stored_string = str; }
  void setStoredString(std::string&& str) { stored_string = std::move(str); }
  void setBlockSetPtr(BlockSet* block_set) { block_set_ptr = block_set; }
  void setStartIdx(int start) { start_idx = start; }
  void setEndIdx(int end) { end_idx = end; }
  void setInsertions(const std::map<int, std::string>& insts) { insertions = insts; }
  void setInsertions(std::map<int, std::string>&& insts) { insertions = std::move(insts); }

  // --- Sequence Recovery ---
  void recoverStoredString() const;
  std::string getConsensusString() const;
  char charAt(size_t pos) const;

  size_t getLength() const;
  void merge(const Consensus& B, const CigarString& cigar, bool inverse = false);

  size_t length() const { return getLength(); }
  size_t size() const { return getLength(); }
  bool empty() const { return getLength() == 0; }

  // --- Substring ---
  Consensus substr(size_t pos, size_t count = std::string::npos) const;

  // --- Debug Print ---
  void print(std::ostream& os = std::cout) const;
};
