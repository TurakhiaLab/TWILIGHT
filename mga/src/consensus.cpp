#include "consensus.hpp"
#include "block_set.hpp"
#include "cigar_util.hpp"
#include <algorithm>
#include <limits>
#include <map>
#include <vector>
#include <string>
#include <utility>

std::string Consensus::buildStringFromIndex() const {
  if (!block_set_ptr) {
    return "";
  }
  const std::string& rep_consensus = block_set_ptr->getAncestralSequence();
  if (rep_consensus.empty()) {
    return "";
  }
  
  std::string result = "";
  for (int i = start_idx; i <= end_idx; ++i) {
    auto it = insertions.find(i);
    if (it != insertions.end()) {
      result += it->second;
    }
    if (i < end_idx && i >= 0 && i < (int)rep_consensus.size()) {
      result += rep_consensus[i];
    }
  }
  return result;
}

void Consensus::recoverStoredString() const {
  if (use_stored_string) {
    return;
  }
  stored_string = buildStringFromIndex();
  use_stored_string = true;

  // 清理舊的 index / reference 結構，轉換為純 string 儲存模式
  insertions.clear();
  start_idx = -1;
  end_idx = -1;
  block_set_ptr = nullptr;
}

std::string Consensus::getConsensusString() const {
  if (use_stored_string) {
    return stored_string;
  }
  return buildStringFromIndex();
}

char Consensus::charAt(size_t pos) const {
  if (use_stored_string) {
    if (pos < stored_string.length()) {
      return stored_string[pos];
    }
    return '\0';
  }

  if (!block_set_ptr) {
    return '\0';
  }

  const std::string& rep_consensus = block_set_ptr->getAncestralSequence();
  size_t cur_local = 0;

  for (int i = start_idx; i <= end_idx; ++i) {
    auto it = insertions.find(i);
    if (it != insertions.end()) {
      const std::string& ins = it->second;
      size_t ins_len = ins.length();
      if (pos >= cur_local && pos < cur_local + ins_len) {
        return ins[pos - cur_local];
      }
      cur_local += ins_len;
    }

    if (i < end_idx && i >= 0 && i < (int)rep_consensus.size()) {
      if (pos == cur_local) {
        return rep_consensus[i];
      }
      cur_local += 1;
    }

    if (cur_local > pos) {
      break;
    }
  }

  return '\0';
}

size_t Consensus::getLength() const {
  if (use_stored_string) {
    return stored_string.length();
  }
  size_t len = 0;
  if (start_idx != -1 && end_idx != -1) {
    len += (end_idx - start_idx);
  }
  for (const auto& pair : insertions) {
    len += pair.second.length();
  }
  return len;
}

void Consensus::merge(const Consensus& B, const CigarString& cigar, bool inverse) {
  if (use_stored_string || !block_set_ptr) {
    std::string seqA = getConsensusString();
    std::string seqB = B.getConsensusString();
    size_t lenA = seqA.length();
    size_t lenB = seqB.length();

    std::string new_seq = "";
    int aPos = 0;
    int bPos = 0;

    for (const auto& op : cigar) {
      int len = op.first;
      char type = op.second;

      if (type == 'M' || type == '=' || type == 'X') {
        if (aPos + len <= (int)lenA) {
          new_seq += seqA.substr(aPos, len);
        } else if (aPos < (int)lenA) {
          new_seq += seqA.substr(aPos);
        }
        aPos += len;
        bPos += len;
      } else if (type == 'D') {
        if (aPos + len <= (int)lenA) {
          new_seq += seqA.substr(aPos, len);
        } else if (aPos < (int)lenA) {
          new_seq += seqA.substr(aPos);
        }
        aPos += len;
      } else if (type == 'I') {
        std::string B_patch;
        if (!inverse) {
          if (bPos + len <= (int)lenB) {
            B_patch = seqB.substr(bPos, len);
          } else if (bPos < (int)lenB) {
            B_patch = seqB.substr(bPos);
          }
        } else {
          size_t real_b_pos = (lenB >= (size_t)bPos + len) ? (lenB - bPos - len) : 0;
          size_t real_len = std::min((size_t)len, lenB - real_b_pos);
          B_patch = getReverseComplement(seqB.substr(real_b_pos, real_len));
        }
        new_seq += B_patch;
        bPos += len;
      }
    }
    if (aPos < (int)lenA) {
      new_seq += seqA.substr(aPos);
    }

    stored_string = std::move(new_seq);
    use_stored_string = true;
    return;
  }
  
  const std::string& rep_consensus = block_set_ptr->getAncestralSequence();
  if (rep_consensus.empty()) {
    return;
  }
  
  size_t lenA = getLength();
  
  struct AMapNode { 
    int idx; 
    int offset; 
  };
  std::vector<AMapNode> A_map(lenA + 1);
  
  int cur_a = 0;
  for (int i = start_idx; i <= end_idx; ++i) {
    auto it = insertions.find(i);
    if (it != insertions.end()) {
      const std::string& ins = it->second;
      for (size_t k = 0; k < ins.length(); ++k) {
        if (cur_a < (int)A_map.size()) {
          A_map[cur_a++] = {i, (int)k};
        }
      }
    }
    if (i < end_idx && i >= 0 && i < (int)rep_consensus.size()) {
      if (cur_a < (int)A_map.size()) {
        A_map[cur_a++] = {i, -1};
      }
    }
  }
  while (cur_a < (int)A_map.size()) {
    A_map[cur_a++] = {end_idx, -1};
  }
  
  int aPos = 0;
  int bPos = 0;
  
  std::map<std::pair<int, int>, std::string> accumulated_inserts;
  
  for (const auto& op : cigar) {
    int len = op.first;
    char type = op.second;
    
    if (type == 'M' || type == '=' || type == 'X') {
      aPos += len;
      bPos += len;
    } else if (type == 'D') {
      aPos += len;
    } else if (type == 'I') {
      std::string B_patch;
      if (!inverse) {
        B_patch = B.substr(bPos, len).getConsensusString();
      } else {
        size_t lenB = B.getLength();
        size_t real_b_pos = (lenB >= (size_t)bPos + len) ? (lenB - bPos - len) : 0;
        B_patch = getReverseComplement(B.substr(real_b_pos, len).getConsensusString());
      }
      
      if (aPos < (int)A_map.size()) {
        AMapNode node = A_map[aPos];
        accumulated_inserts[{node.idx, node.offset}] += B_patch;
      }
      bPos += len;
    }
  }
  
  std::map<int, std::vector<std::pair<int, std::string>>> insertions_to_apply;
  for (const auto& pair : accumulated_inserts) {
    insertions_to_apply[pair.first.first].push_back({pair.first.second, pair.second});
  }
  
  for (auto& pair : insertions_to_apply) {
    int i = pair.first;
    auto& list = pair.second;
    
    std::sort(list.begin(), list.end(), [](const auto& a, const auto& b) {
      return a.first > b.first;
    });
    
    for (const auto& item : list) {
      int offset = item.first;
      const std::string& patch = item.second;
      if (offset >= 0) {
        insertions[i].insert(offset, patch);
      } else {
        insertions[i] += patch;
      }
    }
  }
}

Consensus Consensus::substr(size_t pos, size_t count) const {
  if (use_stored_string) {
    if (pos > stored_string.length()) {
      return Consensus("");
    }
    return Consensus(stored_string.substr(pos, count));
  }

  size_t target_start = pos;
  size_t target_end = (count == std::string::npos || pos + count < pos)
                          ? std::numeric_limits<size_t>::max()
                          : pos + count;

  // ⚡ Fast Path 1: 當無任何 insertions 時，直接 O(1) 算術定位，零迴圈耗時
  if (insertions.empty()) {
    if (start_idx == -1 || end_idx == -1 || start_idx >= end_idx) {
      return Consensus("");
    }
    size_t cons_len = (size_t)(end_idx - start_idx);
    if (target_start >= cons_len) {
      return Consensus("");
    }
    size_t actual_count = std::min(count, cons_len - target_start);
    int new_start_idx = start_idx + (int)target_start;
    int new_end_idx = new_start_idx + (int)actual_count;
    return Consensus(block_set_ptr, new_start_idx, new_end_idx);
  }

  // ⚡ General Path: 區間跳躍演算法 (Jump Search)
  int new_start_idx = -1;
  int new_end_idx = -1;
  std::map<int, std::string> new_insertions;

  int curr_i = start_idx;
  size_t cur_local = 0;
  bool found_any = false;

  auto ins_it = insertions.begin();

  while (curr_i <= end_idx && cur_local < target_end) {
    int next_ins_i = end_idx;
    while (ins_it != insertions.end() && ins_it->first < curr_i) {
      ++ins_it;
    }
    if (ins_it != insertions.end() && ins_it->first <= end_idx) {
      next_ins_i = ins_it->first;
    }

    if (curr_i < next_ins_i) {
      size_t gap_len = (size_t)(next_ins_i - curr_i);
      size_t gap_start_local = cur_local;
      size_t gap_end_local = cur_local + gap_len;

      size_t ov_start = std::max(gap_start_local, target_start);
      size_t ov_end = std::min(gap_end_local, target_end);

      if (ov_start < ov_end) {
        int sub_start_i = curr_i + (int)(ov_start - gap_start_local);
        int sub_end_i = curr_i + (int)(ov_end - gap_start_local);
        if (!found_any) {
          new_start_idx = sub_start_i;
          found_any = true;
        }
        new_end_idx = sub_end_i;
      }

      cur_local += gap_len;
      curr_i = next_ins_i;
      if (cur_local >= target_end) break;
    }

    if (ins_it != insertions.end() && ins_it->first == curr_i && curr_i <= end_idx) {
      const std::string& ins = ins_it->second;
      size_t ins_len = ins.length();
      size_t ins_start = cur_local;
      size_t ins_end = cur_local + ins_len;

      size_t ov_start = std::max(ins_start, target_start);
      size_t ov_end = std::min(ins_end, target_end);

      if (ov_start < ov_end) {
        size_t sub_start = ov_start - ins_start;
        size_t sub_len = ov_end - ov_start;
        new_insertions[curr_i] = ins.substr(sub_start, sub_len);

        if (!found_any) {
          new_start_idx = curr_i;
          found_any = true;
        }
        new_end_idx = curr_i;
      }
      cur_local += ins_len;
      ++ins_it;

      if (curr_i < end_idx) {
        size_t rep_pos = cur_local;
        if (rep_pos >= target_start && rep_pos < target_end) {
          if (!found_any) {
            new_start_idx = curr_i;
            found_any = true;
          }
          new_end_idx = curr_i + 1;
        }
        cur_local += 1;
      }
      curr_i += 1;
    } else {
      break;
    }
  }

  if (!found_any) {
    return Consensus("");
  }

  return Consensus(block_set_ptr, new_start_idx, new_end_idx, std::move(new_insertions));
}

void Consensus::print(std::ostream& os) const {
  os << "Consensus Info:\n";
  if (block_set_ptr) {
    os << "  Reference BlockSet ID: " << block_set_ptr->getId() << "\n";
  } else {
    os << "  Reference BlockSet ID: None (Stored String)\n";
  }
  os << "  Index: [" << start_idx << ", " << end_idx << "]\n";
  os << "  Insertions: ";
  if (insertions.empty()) {
    os << "None\n";
  } else {
    os << "\n";
    for (const auto& pair : insertions) {
      os << "    Position " << pair.first << ": \"" << pair.second << "\"\n";
    }
  }
  os << "  Consensus String Length: " << getLength() << "\n";
}
