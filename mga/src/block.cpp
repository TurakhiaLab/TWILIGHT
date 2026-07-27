
#include "block.hpp"
#include "timer.hpp"
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>


// =======================
// Block Implementation
// =======================
void Block::print(std::ostream& os, int copy) const {
    
    
    os << "┌────────────────────────────────────────────────────────────────────────┐\n";
    os << "│ 📦 BLOCK ID: " << std::left << std::setw(8) << ID;

    // 🌟 根據傳入的 copy 決定印出 Fam ID 還是特定 Copy ID
    if (copy == -1) {
        os << " │ All Homologous ";
    } else {
        os << " │ Copy: " << std::setw(8) << copy;
    }

    std::string distStr;
    if (copy != -1) distStr = this->isDistant(copy) ? "YES" : "NO ";

    os << " │ Len: " << std::setw(6) << consensus.length() << " bp";
    if (copy != -1) os << " │ Distant: " << distStr << " │\n";
    os << "├────────────────────────────────────────────────────────────────────────┤\n";
    
    // ==========================================
    // 1. 讀取並印出相鄰的 Block ID (支援 unordered_map 展開)
    // ==========================================
    os << "  🔗 Linkage:\n";
    auto print_link = [&](int c) {
        std::shared_ptr<Block> p, n;
        p = prev_blocks[c].lock();
        n = next_blocks[c].lock();
        
        os << "     Copy " << std::setw(2) << c << ": ";
        if (p) os << "[Block " << p->getId() << "]"; else os << "[None]";
        os << " <─── (Current) ───> ";
        if (n) os << "[Block " << n->getId() << "]"; else os << "[None]";
        os << "\n";
    };

    if (copy != -1) {
        print_link(copy); // 只印出指定的 copy 路線
    } else {
        for (int c = 0; c < prev_blocks.size(); ++c) print_link(c); // 印出所有分岔路口
    }
    
    // ==========================================
    // 2. 迭代印出底下的所有 Sequences 與 Segments 座標
    // ==========================================
    os << "  🧬 Mapping Sequences:\n";
    
    // 🌟 先預先過濾出合法的 Sequence 與其對應的 Segments，以維持樹狀結構排版的正確性
    auto& mutableSeqs = const_cast<Sequences&>(sequences);
    std::vector<std::pair<std::string, std::vector<Segment*>>> valid_data;
    
    for (const auto& seqPair : sequences) {
        std::vector<Segment*> valid_segs;
        auto& segments = mutableSeqs[seqPair.first].getSegments();
        for (auto& segPairInner : segments) {
            int current_copy = segPairInner.second.getCopyCount(); // 依據你的實作名稱調整
            if (copy == -1 || current_copy == copy) {
                valid_segs.push_back(&segPairInner.second);
            }
        }
        if (!valid_segs.empty()) {
            valid_data.push_back({seqPair.first, valid_segs});
        }
    }

    if (valid_data.empty()) {
        os << "     └─ (No mapped sequences" << (copy != -1 ? " for this copy" : "") << ")\n";
    } else {
        for (size_t i = 0; i < valid_data.size(); ++i) {
            bool isLastSeq = (i == valid_data.size() - 1);
            os << "     " << (isLastSeq ? "└─" : "├─") << " Seq: " << valid_data[i].first << "\n";
            
            auto& segs = valid_data[i].second;
            for (size_t j = 0; j < segs.size(); ++j) {
                bool isLastSeg = (j == segs.size() - 1);
                
                os << "     " << (isLastSeq ? " " : "│") << "  " << (isLastSeg ? "└─" : "├─") 
                   << " Seg: [" << segs[j]->getStart() << " -> " << segs[j]->getEnd() << "]";
                
                // 🌟 如果是 -1 模式，在最後面補上 Copy Number 標籤
                if (copy == -1) {
                    os << " (Copy: " << segs[j]->getCopyCount() << ")";
                }
                os << "\n";
            }
        }
    }
    os << "└────────────────────────────────────────────────────────────────────────┘\n";
}

bool Block::normalizeStrand() {
    bool DEBUG_MODE = false;


    int forward_count = 0;
    int reverse_count = 0;

    // 統計這個 Block 內所有 Segment 的走向
    for (auto& seqPair : sequences) {
        for (auto& segPair : seqPair.second.getSegments()) {
            if (segPair.second.isReverse()) {
                reverse_count++;
            } else {
                forward_count++;
            }
        }
    }

    // ==========================================
    // 🌟 輸出統計與決策結果
    // ==========================================
    if (DEBUG_MODE) {
    std::cout << "  [STRAND-NORM] Block ID: " << this->getId() 
              << " | Forward (+): " << forward_count 
              << " | Reverse (-): " << reverse_count 
              << "  => ";
    }

    // 如果反股 (Inverse) 佔多數，就執行全局翻轉
    if (reverse_count > forward_count) {
        if (DEBUG_MODE) std::cout << "ACTION: FLIP (Reverse Majority)\n";
        this->reverse(); // 呼叫你原本寫好的 reverse 函數
        return true;     // 回傳 true 代表發生了翻轉
    }
    
    if (DEBUG_MODE) std::cout << "ACTION: KEEP (Forward Majority or Tie)\n";
    return false; 
}


// =========================================================
// 🌟 Block 內部：純粹提供幾何特徵，不牽涉 L_min 防護
// =========================================================
std::vector<Block::ColumnSplitScore> Block::calculateSplittingScores(int local_start, int local_end, int left_flank_len, int right_flank_len, bool debug) {
    const double ALPHA = 2.0;
    const double BETA = 1.0;
    const int MIN_BLOCK_LENGTH = 30;

    int scan_len = local_end - local_start + 1;
    std::vector<ColumnSplitScore> scores(scan_len, {0.0, 0.0, 0.0});
    if (scan_len <= 0) return scores;

    int cons_len = this->getConsensus().length();

    // =======================================================
    // 🌟 1. 提取序列資訊與初始狀態 (Initial State Setup)
    // =======================================================
    int N = 0;
    std::vector<std::string> seq_identifiers; // 用於 Debug 印出序列名稱
    for (auto& [seqName, seqData] : this->getSequences()) {
        for (auto& [segId, segNode] : seqData.getSegments()) {
            N++;
            if (debug) {
                seq_identifiers.push_back(seqName + "_seg" + std::to_string(segId));
            }
        }
    }
    if (N == 0) return scores;

    std::vector<int> init_left_lens(N, 0);
    std::vector<int> init_right_lens(N, 0);

    // 🌟 記憶體扁平化：使用單一 1D vector 避免 $2N$ 次 2D Vector 動態記憶體分配
    std::vector<uint8_t> is_gap_at(N * scan_len, 0);
    std::vector<uint8_t> is_boundary_at(N * scan_len, 0);

    int s = 0;
    for (auto& [seqName, seqData] : this->getSequences()) {
        for (auto& [segId, segNode] : seqData.getSegments()) {
            int total_len = std::abs(segNode.getEnd() - segNode.getStart());
            int gaps_before = 0;

            for (auto& var : segNode.getVariants()) {
                if (var.getType() == VariantType::GAP) { 
                    int v_start = var.getStart();
                    int v_end = var.getEnd();

                    // 🧮 計算 local_start 之前總共有多長的 GAP
                    if (v_start < local_start) {
                        gaps_before += std::min(local_start, v_end) - v_start;
                    }

                    // 🚧 標記掃描區間內的 GAP 狀態
                    int overlap_start = std::max(local_start, v_start);
                    int overlap_end = std::min(local_end + 1, v_end);
                    for (int p = overlap_start; p < overlap_end; ++p) {
                        is_gap_at[s * scan_len + (p - local_start)] = 1;
                    }

                    // 🎯 標記 GAP 的起點與終點邊界
                    if (v_start >= local_start && v_start <= local_end) {
                        is_boundary_at[s * scan_len + (v_start - local_start)] = 1;
                    }
                    if (v_end >= local_start && v_end <= local_end) {
                        is_boundary_at[s * scan_len + (v_end - local_start)] = 1;
                    }
                }
            }

            // 直接算出起點的 Left / Right 長度（包含 Flank 補償）
            init_left_lens[s] = left_flank_len + (local_start - gaps_before);
            init_right_lens[s] = right_flank_len + (total_len - (local_start - gaps_before));
            
            s++;
        }
    }

    // 🌟 核心演算法解耦：計算 Prefix Sum 陣列
    // gaps_cnt[seq_idx * (scan_len + 1) + i] 代表 sequence s 在區間 [0, i) 內累計的 GAP 數量
    std::vector<int> gaps_cnt(N * (scan_len + 1), 0);
    for (int seq_idx = 0; seq_idx < N; ++seq_idx) {
        int accum = 0;
        int seq_offset = seq_idx * scan_len;
        int cnt_offset = seq_idx * (scan_len + 1);
        gaps_cnt[cnt_offset] = 0;
        for (int i = 0; i < scan_len; ++i) {
            if (is_gap_at[seq_offset + i]) accum++;
            gaps_cnt[cnt_offset + i + 1] = accum;
        }
    }

    if (debug) {
        std::cout << "  [BLOCK-SCORE-DEBUG] 🏁 Initial State Matrix for Sandbox Range [" << local_start << " -> " << local_end << "]\n";
        for (int s_idx = 0; s_idx < N; ++s_idx) {
            std::cout << "    Seq [" << seq_identifiers[s_idx] << "] -> Init Left Len: " << init_left_lens[s_idx] << ", Init Right Len: " << init_right_lens[s_idx] << "\n";
        }
    }

    // =======================================================
    // 🌟 2. 獨立 Column 分數計算 (支援 TBB 平行化)
    // =======================================================
    auto compute_column = [&](int i) {
        int current_global_pos = local_start + i;
        double bonus = 0.0;
        double left_div = 0.0;
        double right_div = 0.0;
        
        int left_N = 0; 
        int right_N = 0;

        int leftBlockLen = left_flank_len + local_start + i;
        int rightBlockLen = right_flank_len + cons_len - (local_start + i);

        if (debug) {
            std::cout << "    --------------------------------------------------------\n"
                      << "    📊 [COLUMN " << i << "] Global Pos: " << current_global_pos << "\n"
                      << "      ├─ Block Geometry: LeftBlockLen = " << leftBlockLen << ", RightBlockLen = " << rightBlockLen << "\n";
        }

        for (int seq_idx = 0; seq_idx < N; ++seq_idx) {
            int seq_offset = seq_idx * scan_len;
            int cnt_offset = seq_idx * (scan_len + 1);

            // 🌟 O(1) 數學公式獨立計算當前 Column 的 Left / Right 長度
            int num_gaps_before_i = gaps_cnt[cnt_offset + i];
            int non_gaps_before_i = i - num_gaps_before_i;

            int cur_left_len  = init_left_lens[seq_idx] + non_gaps_before_i;
            int cur_right_len = init_right_lens[seq_idx] - non_gaps_before_i;

            if (is_boundary_at[seq_offset + i]) {
                bonus += 1.0; 
                if (debug) {
                    std::cout << "      ├─ 🎉 [BONUS HIT] Boundary detected in " << seq_identifiers[seq_idx] << "\n";
                }
            }
            
            if (debug) {
                std::cout << "      ├─ [" << seq_identifiers[seq_idx] << "] Profile: "
                          << "Gap=" << (is_gap_at[seq_offset + i] ? "YES" : "NO ") << " | "
                          << "Seq_Left=" << cur_left_len << ", Seq_Right=" << cur_right_len;
                if (cur_left_len > 0) std::cout << " (Ratio_L: " << static_cast<double>(leftBlockLen) / cur_left_len << ")";
                if (cur_right_len > 0) std::cout << " (Ratio_R: " << static_cast<double>(rightBlockLen) / cur_right_len << ")";
                std::cout << "\n";
            }

            if (cur_left_len > 0) {
                left_N++;
                left_div += static_cast<double>(leftBlockLen) / static_cast<double>(cur_left_len);
            }
            if (cur_right_len > 0) {
                right_N++;
                right_div += static_cast<double>(rightBlockLen) / static_cast<double>(cur_right_len);
            }
        }

        double mad_left = (left_N == 0) ? 0.0 : (left_div / left_N);
        double mad_right = (right_N == 0) ? 0.0 : (right_div / right_N);
        
        scores[i].perfect_bonus = ALPHA * bonus;
        scores[i].id_left  = (leftBlockLen <= MIN_BLOCK_LENGTH)  ? 10000 : BETA * (mad_left - 1.0);  
        scores[i].id_right = (rightBlockLen <= MIN_BLOCK_LENGTH) ? 10000 : BETA * (mad_right - 1.0);

        if (debug) {
            std::cout << "      └─ 🧮 Result -> Bonus: " << scores[i].perfect_bonus 
                      << " | Score_L: " << (scores[i].id_left == 10000 ? "10000 [PENALTY]" : std::to_string(scores[i].id_left))
                      << " | Score_R: " << (scores[i].id_right == 10000 ? "10000 [PENALTY]" : std::to_string(scores[i].id_right)) << "\n";
        }
    };

    if (debug) {
        for (int i = 0; i < scan_len; ++i) compute_column(i);
    } else {
        tbb::parallel_for(tbb::blocked_range<int>(0, scan_len), [&](const tbb::blocked_range<int>& r) {
            for (int i = r.begin(); i < r.end(); ++i) {
                compute_column(i);
            }
        });
    }

    return scores;
}

BlockPtrPair Block::split(int cut, int copy) const {
    if (cut <= 0 || cut >= (int)consensus.length()) {
        return {nullptr, nullptr};
    }

    // global_timer.start("Block::consensus");
    // consensus.print();
    auto leftBlock =  std::make_shared<Block>(999991, std::move(consensus.substr(0, cut)));
    auto rightBlock = std::make_shared<Block>(999992, std::move(consensus.substr(cut)));
    global_timer.stop("Block::consensus");
    // leftBlock->getConsensus().print();
    // rightBlock->getConsensus().print();
    // global_timer.print();
    std::set<int> active_copies;
    global_timer.start("Block::split");
    for (auto [seqID, seqInfo] : this->sequences) {
        Sequence leftSeqInfo(seqID);
        Sequence rightSeqInfo(seqID);
        
        for (auto& [segStart, seg] : seqInfo.getSegments()) {
            Segment oldSeg = seg;
            if (copy != -1 && oldSeg.getCopyCount() != copy) {
                continue;
            }
            active_copies.insert(oldSeg.getCopyCount());

            auto splitSegs = oldSeg.split(cut);
            Segment leftSeg = splitSegs.first;
            Segment rightSeg = splitSegs.second;

            bool validLeft = (leftSeg.getStart() != leftSeg.getEnd());
            bool validRight = (rightSeg.getStart() != rightSeg.getEnd());

            if (validLeft && validRight) {
                if (!oldSeg.isReverse()) {
                    leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                    leftSeg.setNextBlock(rightBlock);
                    rightSeg.setPrevBlock(leftBlock);
                    rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
                } else {
                    rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                    rightSeg.setNextBlock(leftBlock);
                    leftSeg.setPrevBlock(rightBlock);
                    leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
                }
            } else if (validLeft && !validRight) {
                leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
            } else if (!validLeft && validRight) {
                rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
            }

            if (validLeft)  leftSeqInfo.addSegment(leftSeg);
            if (validRight) rightSeqInfo.addSegment(rightSeg);
        }

        if (!leftSeqInfo.getSegments().empty())  leftBlock->addSequence(std::move(leftSeqInfo));
        if (!rightSeqInfo.getSegments().empty()) rightBlock->addSequence(std::move(rightSeqInfo));
    }
    global_timer.stop("Block::split");

    global_timer.start("Block::SetPtr");

    for (int c = 0; c < (int)this->prev_blocks.size(); ++c) {
        if (copy != -1 && c != copy) continue;
        if (auto p = this->prev_blocks[c].lock()) {
            leftBlock->setPrevBlock(c, p);
        }
    }
    for (int c = 0; c < (int)this->next_blocks.size(); ++c) {
        if (copy != -1 && c != copy) continue;
        if (auto n = this->next_blocks[c].lock()) {
            rightBlock->setNextBlock(c, n);
        }
    }

    for (int c : active_copies) {
        leftBlock->setNextBlock(c, rightBlock);
        rightBlock->setPrevBlock(c, leftBlock);
    }

    global_timer.stop("Block::SetPtr");

    return {leftBlock, rightBlock};
}

int Block::getMaxCopy() const {
    int max_c = -1;
    for (auto& [seq, info] : this->sequences) {
        for (auto& [s, seg] : const_cast<Sequence&>(info).getSegments()) {
            max_c = std::max(max_c, seg.getCopyCount());
        }
    }
    return std::max(0, max_c);
}

void Block::applyCopyAssignment(int target_copy, int shift_amount) {
    for (auto& [seq, info] : this->sequences) {
        for (auto& [s, seg] : info.getSegments()) {
            if (target_copy != -1) seg.setCopyCount(target_copy);
            else if (shift_amount > 0) seg.setCopyCount(seg.getCopyCount() + shift_amount);
        }
    }
}

/*
void Block::print(std::ostream& os) const {
    if (sequences_.size() <= 1) return;

    const std::string border = "==================================================";
    const std::string subBorder = "--------------------------------------------------";

    os << "\n" << border << "\n";
    os << " BLOCK ID      : " << id_ << "\n";
    os << " CONSENSUS     : " << consensus_sequence_.size() << "\n";

    auto print_links = [&](const std::string& label, const auto& links) {
        os << label;
        if (links.empty()) {
            os << "[None]\n";
        } else {
            os << "{ ";
            for (const auto& weak_link : links) {
                if (auto link = weak_link.lock()) {
                    os << link->getId() << " ";
                } else {
                    os << "[exp] ";
                }
            }
            os << "}\n";
        }
    };
    
    print_links(" PREV BLOCKS   : ", prev_blocks_);
    print_links(" NEXT BLOCKS   : ", next_blocks_);

    os << subBorder << "\n";
    os << " CONTAINED SEQUENCES (" << sequences_.size() << ")\n";
    os << subBorder << "\n";

    if (sequences_.empty()) {
        os << "  (No sequences mapped)\n";
    }

    // for (const auto& seq : sequences_) {
    //     os << "  > ID: " << std::left << std::setw(15) << seq.sequence_id 
    //        << " Range: [" << seq.start_coordinate << "-" << seq.end_coordinate << "]\n";
    // }
    os << border << "\n";
}

void Block::refine() {
    int consLen = consensus_sequence_.length();
    if (consLen == 0) return;

    // ==========================================
    // 內部神器 1：提取純粹的 Native Forward 序列
    // ==========================================
    auto getNativeForwardSequence = [&](Segment& seg) {
        std::string mutated = consensus_sequence_;
        for (auto& v : seg.getVariants()) {
            if (v.getType() == Variation::SNV && v.getStart() < mutated.length()) {
                mutated[v.getStart()] = v.getAlt();
            }
        }
        std::string raw = "";
        int cur = 0;
        std::vector<Variation> sorted_vars = seg.getVariants();
        std::sort(sorted_vars.begin(), sorted_vars.end(), [](Variation& a, Variation& b){
            return a.getStart() < b.getStart();
        });
        for (auto& v : sorted_vars) {
            if (v.getType() == Variation::GAP) {
                if (v.getStart() > cur) raw += mutated.substr(cur, v.getStart() - cur);
                cur = v.getEnd();
            }
        }
        if (cur < mutated.length()) raw += mutated.substr(cur);
        
        if (seg.isReverse()) {
            std::string rc = raw;
            std::reverse(rc.begin(), rc.end());
            for (char& c : rc) {
                switch (c) {
                    case 'A': c = 'T'; break; case 'T': c = 'A'; break;
                    case 'C': c = 'G'; break; case 'G': c = 'C'; break;
                    case 'a': c = 't'; break; case 't': c = 'a'; break;
                    case 'c': c = 'g'; break; case 'g': c = 'c'; break;
                }
            }
            raw = rc;
        }
        return raw;
    };

    // ==========================================
    // 內部神器 2：取得 Segment 去頭去尾 Gap 後的核心區間
    // ==========================================
    auto extractCoreVars = [](Segment& seg, int cLen) {
        auto& vars = seg.getVariants();
        if (vars.empty()) return std::make_pair(0, cLen); 
        int c_start = (vars.front().getType() == Variation::GAP && vars.front().getStart() == 0) ? vars.front().getEnd() : 0;
        int c_end = (vars.back().getType() == Variation::GAP && vars.back().getEnd() == cLen) ? vars.back().getStart() : cLen;
        if (vars.size() == 1 && c_start > 0 && c_end < cLen) return std::make_pair(cLen, 0); 
        return std::make_pair(c_start, c_end);
    };

    // ==========================================
    // Phase 0: 孤島區塊大一統 (Singleton Block Collapse)
    // ==========================================
    if (sequences_.size() == 1) {
        auto& seqPair = *sequences_.begin();
        auto& segments_map = seqPair.second.getSegments();
        
        if (!segments_map.empty()) {
            bool all_contiguous = true;
            auto it = segments_map.begin();
            auto prev_it = it;
            ++it;
            while (it != segments_map.end()) {
                if (prev_it->second.getEnd() != it->second.getStart()) {
                    all_contiguous = false;
                    break;
                }
                prev_it = it; ++it;
            }

            if (all_contiguous) {
                std::string new_consensus = "";
                for (auto& pair : segments_map) {
                    new_consensus += getNativeForwardSequence(pair.second);
                }

                Segment first_seg = segments_map.begin()->second;
                Segment last_seg = segments_map.rbegin()->second;
                Segment merged_seg(first_seg.getStart(), last_seg.getEnd());
                merged_seg.setReverse(false); 

                // 🌟 1. 先萃取出絕對的「左邊」與「右邊」鄰居
                auto left_conn = first_seg.isReverse() ? first_seg.getNextBlock() : first_seg.getPrevBlock();
                auto right_conn = last_seg.isReverse() ? last_seg.getPrevBlock() : last_seg.getNextBlock();

                // 🌟 2. 因為 merged_seg 強制被設為正股 (Forward)，所以左邊接 Prev，右邊接 Next
                merged_seg.setPrevBlock(left_conn.lock());
                merged_seg.setNextBlock(right_conn.lock());
                
                merged_seg.getVariants().clear(); 

                consensus_sequence_ = new_consensus;
                segments_map.clear();
                segments_map[merged_seg.getStart()] = merged_seg;
                return;
            }
        }
    }

    // ==========================================
    // Phase 1.1: Self-Healing Segment Terminal Gaps 
    // ==========================================
    for (auto& seqPair : sequences_) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            int seq_len = seg.getEnd() - seg.getStart();
            int total_gap = 0;
            for (auto& v : seg.getVariants()) {
                if (v.getType() == Variation::GAP) total_gap += (v.getEnd() - v.getStart());
            }
            int claimed_cons = seq_len + total_gap;
            
            if (claimed_cons < consLen) {
                seg.getVariants().push_back(Variation::createGap(claimed_cons, consLen));
            }
            
            auto& vars = seg.getVariants();
            std::sort(vars.begin(), vars.end(), [](Variation& a, Variation& b){
                if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                return a.getType() > b.getType(); 
            });
            std::vector<Variation> clean_vars;
            for (auto& v : vars) {
                if (!clean_vars.empty() && clean_vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && clean_vars.back().getEnd() >= v.getStart()) {
                    int mStart = clean_vars.back().getStart();
                    int mEnd = std::max(clean_vars.back().getEnd(), v.getEnd());
                    clean_vars.pop_back();
                    clean_vars.push_back(Variation::createGap(mStart, mEnd));
                } else {
                    clean_vars.push_back(v);
                }
            }
            vars = std::move(clean_vars);
        }
    }

    // ==========================================
    // Phase 1.2: 計算 Coverage 並精準剔除 All-Gap Columns
    // ==========================================
    std::vector<int> coverage(consLen, 0);
    for (auto& seqPair : sequences_) {
        for (auto& segPair : seqPair.second.getSegments()) {
            for (int i = 0; i < consLen; ++i) coverage[i]++;
            for (auto& var : segPair.second.getVariants()) {
                if (var.getType() == Variation::GAP) {
                    int vs = std::max(0, std::min(consLen, var.getStart()));
                    int ve = std::max(0, std::min(consLen, var.getEnd()));
                    for (int i = vs; i < ve; ++i) coverage[i]--;
                }
            }
        }
    }

    std::string new_consensus = "";
    new_consensus.reserve(consLen);
    std::vector<int> offset(consLen + 1, 0); 
    int current_offset = 0;
    for (int i = 0; i < consLen; ++i) {
        if (coverage[i] <= 0) current_offset++; 
        else new_consensus += consensus_sequence_[i];
        offset[i + 1] = current_offset; 
    }
    
    consensus_sequence_ = new_consensus; 
    consLen = consensus_sequence_.length();

    for (auto& seqPair : sequences_) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second; 
            std::vector<Variation> new_vars;
            for (auto& var : seg.getVariants()) {
                int safe_v_start = std::max(0, std::min((int)offset.size()-1, var.getStart()));
                int safe_v_end   = std::max(0, std::min((int)offset.size()-1, var.getEnd()));
                int new_v_start = safe_v_start - offset[safe_v_start];
                int new_v_end   = safe_v_end - offset[safe_v_end];

                if (new_v_start < new_v_end) {
                    if (var.getType() == Variation::GAP) new_vars.push_back(Variation::createGap(new_v_start, new_v_end));
                    else new_vars.push_back(Variation(new_v_start, var.getAlt()));
                }
            }
            seg.getVariants() = std::move(new_vars); 
        }
    }

    // ==========================================
    // Phase 2: Anchor-based Chunk Merging 
    // ==========================================
    for (auto& seqPair : sequences_) {
        auto& segments_map = seqPair.second.getSegments();
        if (segments_map.empty()) continue;

        std::vector<Segment> seg_list;
        for (auto& p : segments_map) seg_list.push_back(p.second);

        std::vector<std::vector<Segment>> chunks;
        chunks.push_back({seg_list[0]});
        for (size_t i = 1; i < seg_list.size(); ++i) {
            if (chunks.back().back().getEnd() == seg_list[i].getStart()) {
                chunks.back().push_back(seg_list[i]);
            } else {
                chunks.push_back({seg_list[i]});
            }
        }

        std::map<int, Segment> new_segments_map;
        
        for (size_t chunk_idx = 0; chunk_idx < chunks.size(); ++chunk_idx) {
            auto& chunk = chunks[chunk_idx];
            
            if (chunk.size() == 1) {
                new_segments_map[chunk[0].getStart()] = chunk[0];
                continue;
            }

            int max_len = -1;
            int anchor_idx = -1;
            for (size_t i = 0; i < chunk.size(); ++i) {
                int len = chunk[i].getEnd() - chunk[i].getStart();
                if (len > max_len) { max_len = len; anchor_idx = i; }
            }
            Segment anchor = chunk[anchor_idx];
            bool target_strand = anchor.isReverse();

            std::vector<std::string> native_fwds(chunk.size());
            for (size_t i = 0; i < chunk.size(); ++i) {
                native_fwds[i] = getNativeForwardSequence(chunk[i]);
            }

            std::string full_seq = "";
            int L_add = 0, R_add = 0; 

            auto rc = [](const std::string& s) {
                std::string rev = s; std::reverse(rev.begin(), rev.end());
                for(char& c: rev){ if(c=='A') c='T'; else if(c=='T') c='A'; else if(c=='C') c='G'; else if(c=='G') c='C'; }
                return rev;
            };

            if (!target_strand) { 
                for (size_t i = 0; i < chunk.size(); ++i) {
                    full_seq += native_fwds[i];
                    if (i < anchor_idx) L_add += native_fwds[i].length();
                    if (i > anchor_idx) R_add += native_fwds[i].length();
                }
            } else { 
                for (int i = chunk.size() - 1; i >= 0; --i) {
                    std::string rev_seq = rc(native_fwds[i]);
                    full_seq += rev_seq;
                    if (i > anchor_idx) L_add += rev_seq.length(); 
                    if (i < anchor_idx) R_add += rev_seq.length(); 
                }
            }

            auto [a_start, a_end] = extractCoreVars(anchor, consLen);
            int G_left = a_start;
            int G_right = consLen - a_end;

            int over_L = std::max(0, L_add - G_left);
            int over_R = std::max(0, R_add - G_right);

            if (over_L > 0 || over_R > 0) {
                std::string prepend_str = (over_L > 0) ? full_seq.substr(0, over_L) : "";
                std::string append_str = (over_R > 0) ? full_seq.substr(full_seq.length() - over_R, over_R) : "";

                consensus_sequence_ = prepend_str + consensus_sequence_ + append_str;
                int old_consLen = consLen;
                consLen = consensus_sequence_.length();

                // 🌟 【修正 2】：全域 Gap 延長函數 (確保不遺漏任何備份)
                auto applyConsensusExpansion = [&](Segment& sg) {
                    auto& vars = sg.getVariants();
                    for (auto& v : vars) v.shift(over_L); // 全部往右推

                    // 處理左邊界：延長原有的端點 Gap，或補上新的 Gap
                    if (!vars.empty() && vars.front().getType() == Variation::GAP && vars.front().getStart() == over_L) {
                        vars.front().setStart(0);
                    } else if (over_L > 0) {
                        vars.insert(vars.begin(), Variation::createGap(0, over_L));
                    }

                    // 處理右邊界：延長原有的端點 Gap，或補上新的 Gap
                    if (!vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == old_consLen + over_L) {
                        vars.back().setEnd(consLen);
                    } else if (over_R > 0) {
                        vars.push_back(Variation::createGap(old_consLen + over_L, consLen));
                    }
                };

                // A. 更新 Block 中其他的 Sequences
                for (auto& sq : sequences_) {
                    for (auto& sg : sq.second.getSegments()) applyConsensusExpansion(sg.second);
                }
                
                // B. 同步更新本 Sequence 中「已經處理完」的 Chunk
                for (auto& p : new_segments_map) applyConsensusExpansion(p.second);
                
                // C. 【關鍵修復】同步更新本 Sequence 中「尚未處理」的備份 Chunk！
                for (size_t future_c = chunk_idx + 1; future_c < chunks.size(); ++future_c) {
                    for (auto& seg : chunks[future_c]) applyConsensusExpansion(seg);
                }
            }

            int new_core_start = (a_start + over_L) - L_add;
            int new_core_end = new_core_start + full_seq.length();

            Segment merged_seg(chunk.front().getStart(), chunk.back().getEnd());
            merged_seg.setReverse(target_strand);

            Segment& first_seg = chunk.front();
            Segment& last_seg = chunk.back();
                    
            // 🌟 1. 一樣先萃取出絕對的「左邊」與「右邊」鄰居
            auto left_conn = first_seg.isReverse() ? first_seg.getNextBlock() : first_seg.getPrevBlock();
            auto right_conn = last_seg.isReverse() ? last_seg.getPrevBlock() : last_seg.getNextBlock();
                    
            // 🌟 2. 根據 merged_seg 最終決定要當正股還是反股，來分配接頭
            if (!target_strand) { 
                // merged_seg 是正股：左進右出
                merged_seg.setPrevBlock(left_conn.lock());
                merged_seg.setNextBlock(right_conn.lock());
            } else { 
                // merged_seg 是反股：右進左出
                merged_seg.setPrevBlock(right_conn.lock());
                merged_seg.setNextBlock(left_conn.lock());
            }

            std::vector<Variation> new_vars;
            if (new_core_start > 0) new_vars.push_back(Variation::createGap(0, new_core_start));
            
            for (int i = 0; i < full_seq.length(); ++i) {
                int pos = new_core_start + i;
                if (full_seq[i] != consensus_sequence_[pos]) {
                    new_vars.push_back(Variation(pos, full_seq[i]));
                }
            }
            
            if (new_core_end < consLen) new_vars.push_back(Variation::createGap(new_core_end, consLen));

            std::vector<Variation> clean_vars;
            for (auto& v : new_vars) {
                if (!clean_vars.empty() && clean_vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && clean_vars.back().getEnd() >= v.getStart()) {
                    int mStart = clean_vars.back().getStart();
                    int mEnd = std::max(clean_vars.back().getEnd(), v.getEnd());
                    clean_vars.pop_back();
                    clean_vars.push_back(Variation::createGap(mStart, mEnd));
                } else clean_vars.push_back(v);
            }

            merged_seg.getVariants() = clean_vars;
            new_segments_map[merged_seg.getStart()] = merged_seg;
        }
        
        segments_map = std::move(new_segments_map);
    }
}

*/

void Block::refineConsensusAndVariants() {
    // 1. 先確保原生的 Consensus index/reference 恢復並儲存為實體 stored_string
    consensus.recoverStoredString();
    std::string cons_str = consensus.getConsensusString();
    int consLen = cons_str.length();
    if (consLen == 0) return;

return;
    // 2. 統計每個位點各鹼基的覆蓋數量
    std::vector<std::unordered_map<char, int>> base_counts(consLen);
    
    for (auto& seqPair : sequences) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            
            std::vector<bool> is_gap(consLen, false);
            std::unordered_map<int, char> snv_at_pos;
            
            for (auto& var : seg.getVariants()) {
                if (var.getType() == VariantType::GAP) {
                    for (int pos = var.getStart(); pos < var.getEnd() && pos < consLen; ++pos) {
                        is_gap[pos] = true;
                    }
                } else if (var.getType() == VariantType::SNV) {
                    if (var.getStart() < consLen) {
                        snv_at_pos[var.getStart()] = var.getAlt();
                    }
                }
            }
            
            for (int i = 0; i < consLen; ++i) {
                if (is_gap[i]) continue;
                
                auto it = snv_at_pos.find(i);
                if (it != snv_at_pos.end()) {
                    base_counts[i][it->second]++;
                } else {
                    base_counts[i][cons_str[i]]++;
                }
            }
        }
    }

    // 3. 多數決 (Majority Voting) 找出新的 Consensus 鹼基
    std::string new_cons = cons_str;
    bool consensus_changed = false;

    for (int i = 0; i < consLen; ++i) {
        if (base_counts[i].empty()) continue;
        
        char current_base = cons_str[i];
        int current_count = base_counts[i].count(current_base) ? base_counts[i][current_base] : 0;
        
        char max_base = current_base;
        int max_count = current_count;
        
        for (const auto& pair : base_counts[i]) {
            if (pair.second > max_count) {
                max_count = pair.second;
                max_base = pair.first;
            }
        }
        
        if (max_base != current_base) {
            new_cons[i] = max_base;
            consensus_changed = true;
        }
    }

    // 4. 若 Consensus 發生更新，精確重構各 Segment 的 SNV Variants
    if (consensus_changed) {
        for (auto& seqPair : sequences) {
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;
                auto& vars = seg.getVariants();
                
                std::vector<bool> is_gap(consLen, false);
                std::unordered_map<int, char> existing_snv;
                
                for (auto& var : vars) {
                    if (var.getType() == VariantType::GAP) {
                        for (int pos = var.getStart(); pos < var.getEnd() && pos < consLen; ++pos) {
                            is_gap[pos] = true;
                        }
                    } else if (var.getType() == VariantType::SNV) {
                        existing_snv[var.getStart()] = var.getAlt();
                    }
                }
                
                std::vector<Variant> updated_vars;
                
                // 保留所有 GAP
                for (auto& var : vars) {
                    if (var.getType() == VariantType::GAP) {
                        updated_vars.push_back(var);
                    }
                }
                
                // 處理所有位點的有效鹼基
                for (int i = 0; i < consLen; ++i) {
                    if (is_gap[i]) continue;
                    
                    char seg_effective_base = cons_str[i];
                    auto it = existing_snv.find(i);
                    if (it != existing_snv.end()) {
                        seg_effective_base = it->second;
                    }
                    
                    if (seg_effective_base != new_cons[i]) {
                        updated_vars.push_back(Variant(i, seg_effective_base));
                    }
                }
                
                std::sort(updated_vars.begin(), updated_vars.end(), [](Variant& a, Variant& b) {
                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                    return a.getType() > b.getType();
                });
                
                seg.getVariants() = std::move(updated_vars);
            }
        }
    }

    // 5. 將修正後的多數 Consensus 實體字串寫回並切換模式
    consensus.setStoredString(std::move(new_cons));
    consensus.setUseStoredString(true);
}