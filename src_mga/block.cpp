
#include "block.hpp"
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <functional>
#include <unordered_set>


// =======================
// Block Implementation
// =======================
void Block::print(std::ostream& os) const {
    
    std::string distStr = distant ? "YES" : "NO "; 
    os << "┌────────────────────────────────────────────────────────────────────────┐\n";
    os << "│ 📦 BLOCK ID: " << std::left << std::setw(8) << ID 
       << " │ Fam ID: " << std::setw(6) << family_ID 
       << " │ Len: " << std::setw(6) << consensus.length() << " bp"
       << " │ Distant: " << distStr << " │\n";
    os << "├────────────────────────────────────────────────────────────────────────┤\n";
    
    // 1. 讀取並印出相鄰的 Block ID
    auto prev = prev_block.lock();
    auto next = next_block.lock();
    os << "  🔗 Linkage: ";
    if (prev) os << "[Block " << prev->getId() << "]"; else os << "[None]";
    os << " <─── (Current) ───> ";
    if (next) os << "[Block " << next->getId() << "]"; else os << "[None]";
    os << "\n";
    
    // 2. 迭代印出底下的所有 Sequences 與 Segments 座標
    os << "  🧬 Mapping Sequences:\n";
    if (sequences.empty()) {
        os << "     └─ (No mapped sequences)\n";
    } else {
        size_t seqIdx = 0;
        for (const auto& seqPair : sequences) {
            seqIdx++;
            // 判斷是否為最後一條 sequence，用來決定樹狀圖線條
            bool isLastSeq = (seqIdx == sequences.size());
            os << "     " << (isLastSeq ? "└─" : "├─") << " Seq: " << seqPair.first << "\n";
            
            // 為了在 const 函數中讀取 segments，我們可以使用 const 走訪
            // 如果你的 getSegments() 沒有 const 多載，可以使用 const_cast 繞過
            auto& mutableSeqs = const_cast<Sequences&>(sequences);
            auto& segments = mutableSeqs[seqPair.first].getSegments();
            
            size_t segIdx = 0;
            for (auto& segPairInner : segments) {
                segIdx++;
                bool isLastSeg = (segIdx == segments.size());
                
                // 結構線條排版
                os << "     " << (isLastSeq ? " " : "│") << "  " << (isLastSeg ? "└─" : "├─") 
                   << " Seg: [" << segPairInner.second.getStart() << " -> " 
                   << segPairInner.second.getEnd() << "]\n";
            }
        }
    }
    os << "└────────────────────────────────────────────────────────────────────────┘\n";
}


bool Block::normalizeStrand() {
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
    // 如果反股 (Inverse) 佔多數，就執行全局翻轉
    if (reverse_count > forward_count) {
        this->reverse(); // 呼叫你原本寫好的 reverse 函數
        return true;     // 回傳 true 代表發生了翻轉
    }
    return false; 
}


// =========================================================
// 🌟 Block 內部：純粹提供幾何特徵，不牽涉 L_min 防護
// =========================================================
std::vector<Block::ColumnSplitScore> Block::calculateSplittingScores(int local_start, int local_end) {
    // 🌟 控制此函數內部的 Debug 訊息開關
    bool debug = false; 

    const double ALPHA = 2.0;
    const double BETA = 1.0;
    const int MIN_BLOCK_LENGTH = 30;

    int scan_len = local_end - local_start + 1;
    std::vector<ColumnSplitScore> scores(scan_len, {0.0, 0.0, 0.0});
    if (scan_len <= 0) return scores;

    int cons_len = this->getConsensus().length();

    // =======================================================
    // 🌟 1. 提取序列資訊並計算初始狀態 (Initial State Setup)
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

    std::vector<int> left_lens(N, 0);
    std::vector<int> right_lens(N, 0);
    
    std::vector<std::vector<bool>> is_gap_at(N, std::vector<bool>(scan_len, false));
    std::vector<std::vector<bool>> is_boundary_at(N, std::vector<bool>(scan_len, false));

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

                    // 🚧 標記掃描區間內的 GAP 凍結狀態
                    int overlap_start = std::max(local_start, v_start);
                    int overlap_end = std::min(local_end + 1, v_end);
                    for (int p = overlap_start; p < overlap_end; ++p) {
                        is_gap_at[s][p - local_start] = true;
                    }

                    // 🎯 精準標記 GAP 的起點與終點邊界
                    if (v_start >= local_start && v_start <= local_end) {
                        is_boundary_at[s][v_start - local_start] = true;
                    }
                    if (v_end >= local_start && v_end <= local_end) {
                        is_boundary_at[s][v_end - local_start] = true;
                    }
                }
            }

            // 直接算出起點的 Left / Right 長度
            left_lens[s] = local_start - gaps_before;
            right_lens[s] = total_len - left_lens[s];
            
            s++;
        }
    }

    if (debug) {
        std::cout << "  [BLOCK-SCORE-DEBUG] 🏁 Initial State Matrix for Sandbox Range [" << local_start << " -> " << local_end << "]\n";
        for (int s_idx = 0; s_idx < N; ++s_idx) {
            std::cout << "    Seq [" << seq_identifiers[s_idx] << "] -> Init Left Len: " << left_lens[s_idx] << ", Init Right Len: " << right_lens[s_idx] << "\n";
        }
    }

    // =======================================================
    // 🌟 2. 狀態機推進：掃描 B -> D 區間 (State Progression)
    // =======================================================
    for (int i = 0; i < scan_len; ++i) {
        int current_global_pos = local_start + i;
        double bonus = 0.0;
        double left_div = 0.0;
        double right_div = 0.0;
        
        int left_N = 0; 
        int right_N = 0;

        int leftBlockLen = local_start + i;
        int rightBlockLen = cons_len - leftBlockLen;

        if (debug) {
            std::cout << "    --------------------------------------------------------\n"
                      << "    📊 [COLUMN " << i << "] Global Pos: " << current_global_pos << "\n"
                      << "      ├─ Block Geometry: LeftBlockLen = " << leftBlockLen << ", RightBlockLen = " << rightBlockLen << "\n";
        }

        for (int s = 0; s < N; ++s) {
            // 如果踩到 GAP 邊界，直接加上 Bonus
            if (is_boundary_at[s][i]) {
                bonus += 1.0; 
                if (debug) {
                    std::cout << "      ├─ 🎉 [BONUS HIT] Boundary detected in " << seq_identifiers[s] << "\n";
                }
            }
            
            if (debug) {
                std::cout << "      ├─ [" << seq_identifiers[s] << "] Profile: "
                          << "Gap=" << (is_gap_at[s][i] ? "YES" : "NO ") << " | "
                          << "Seq_Left=" << left_lens[s] << ", Seq_Right=" << right_lens[s];
                if (left_lens[s] > 0) std::cout << " (Ratio_L: " << static_cast<double>(leftBlockLen) / left_lens[s] << ")";
                if (right_lens[s] > 0) std::cout << " (Ratio_R: " << static_cast<double>(rightBlockLen) / right_lens[s] << ")";
                std::cout << "\n";
            }

            if (left_lens[s] > 0) {
                left_N++;
                left_div += static_cast<double>(leftBlockLen) / static_cast<double>(left_lens[s]);
            }
            if (right_lens[s] > 0) {
                right_N++;
                right_div += static_cast<double>(rightBlockLen) / static_cast<double>(right_lens[s]);
            }
        }

        double mad_left = (left_N == 0) ? 0.0 : (left_div / left_N);
        double mad_right = (right_N == 0) ? 0.0 : (right_div / right_N);
        
        // 📝 裝填當前 Column 的回傳數據
        scores[i].perfect_bonus = ALPHA * bonus;
        scores[i].id_left  = (leftBlockLen <= MIN_BLOCK_LENGTH)  ? 10000 : BETA * (mad_left - 1.0);  
        scores[i].id_right = (rightBlockLen <= MIN_BLOCK_LENGTH) ? 10000 : BETA * (mad_right - 1.0);

        if (debug) {
            std::cout << "      └─ 🧮 Result -> Bonus: " << scores[i].perfect_bonus 
                      << " | Score_L: " << (scores[i].id_left == 10000 ? "10000 [PENALTY]" : std::to_string(scores[i].id_left))
                      << " | Score_R: " << (scores[i].id_right == 10000 ? "10000 [PENALTY]" : std::to_string(scores[i].id_right)) << "\n";
        }

        // 🚀 狀態推進
        if (debug && i < scan_len - 1) std::cout << "      🚀 Advancing State to Next Column...\n";
        for (int s = 0; s < N; ++s) {
            if (!is_gap_at[s][i]) {
                left_lens[s] += 1;
                right_lens[s] -= 1;
            }
        }
    }

    return scores;
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