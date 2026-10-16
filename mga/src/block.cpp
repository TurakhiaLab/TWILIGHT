
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

    // Determine whether to print family ID or specific copy ID based on input copy
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
    // 1. Read and print adjacent Block IDs (supports unordered_map traversal)
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
        print_link(copy); // Print only the specified copy route
    } else {
        for (int c = 0; c < prev_blocks.size(); ++c) print_link(c); // Print all junction routes
    }
    
    // ==========================================
    // 2. Iterate and print all underlying Sequences and Segment coordinates
    // ==========================================
    os << "  🧬 Mapping Sequences:\n";
    
    // Pre-filter valid Sequences and their corresponding Segments to maintain tree layout correctness
    auto& mutableSeqs = const_cast<Sequences&>(sequences);
    std::vector<std::pair<std::string, std::vector<Segment*>>> valid_data;
    
    for (const auto& seqPair : sequences) {
        std::vector<Segment*> valid_segs;
        auto& segments = mutableSeqs[seqPair.first].getSegments();
        for (auto& segPairInner : segments) {
            int current_copy = segPairInner.second.getCopyCount(); // Adjust according to implementation
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
                
                // In -1 mode, append copy number label at the end
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

    // Count strand orientation of all Segments within this Block
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
    // Output statistics and decision results
    // ==========================================
    if (DEBUG_MODE) {
    std::cout << "  [STRAND-NORM] Block ID: " << this->getId() 
              << " | Forward (+): " << forward_count 
              << " | Reverse (-): " << reverse_count 
              << "  => ";
    }

    // Flip globally if reverse strand is the majority
    if (reverse_count > forward_count) {
        if (DEBUG_MODE) std::cout << "ACTION: FLIP (Reverse Majority)\n";
        this->reverse(); // Call reverse function
        return true;     // Return true indicating flip occurred
    }
    
    if (DEBUG_MODE) std::cout << "ACTION: KEEP (Forward Majority or Tie)\n";
    return false; 
}


// =========================================================
// Block internal: provides geometric features only, without L_min protection
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
    // 1. Extract sequence information and initial state (Initial State Setup)
    // =======================================================
    int N = 0;
    std::vector<std::string> seq_identifiers; // Sequence identifiers for debug logging
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

    // Flatten memory: use a single 1D vector to avoid 2N dynamic 2D vector allocations
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

                    // Calculate total GAP length preceding local_start
                    if (v_start < local_start) {
                        gaps_before += std::min(local_start, v_end) - v_start;
                    }

                    // Mark GAP status within scan interval
                    int overlap_start = std::max(local_start, v_start);
                    int overlap_end = std::min(local_end + 1, v_end);
                    for (int p = overlap_start; p < overlap_end; ++p) {
                        is_gap_at[s * scan_len + (p - local_start)] = 1;
                    }

                    // Mark start and end boundaries of GAP
                    if (v_start >= local_start && v_start <= local_end) {
                        is_boundary_at[s * scan_len + (v_start - local_start)] = 1;
                    }
                    if (v_end >= local_start && v_end <= local_end) {
                        is_boundary_at[s * scan_len + (v_end - local_start)] = 1;
                    }
                }
            }

            // Calculate initial Left / Right lengths directly (including flank compensation)
            init_left_lens[s] = left_flank_len + (local_start - gaps_before);
            init_right_lens[s] = right_flank_len + (total_len - (local_start - gaps_before));
            
            s++;
        }
    }

    // Core algorithm decoupling: compute prefix sum array
    // gaps_cnt[seq_idx * (scan_len + 1) + i] represents cumulative GAP count in [0, i) for sequence s
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
    // 2. Independent column score computation (supports TBB parallelization)
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

            // O(1) formula to compute Left / Right lengths independently for current column
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


void Block::refineConsensusAndVariants() {
    // 1. Recover native consensus index/reference and store as concrete stored_string
    consensus.recoverStoredString();
    std::string cons_str = consensus.getConsensusString();
    int consLen = cons_str.length();
    if (consLen == 0) return;

return;
    // 2. Count base coverage at each position
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

    // 3. Majority voting to identify updated consensus bases
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

    // 4. If consensus changed, reconstruct SNV variants for each segment
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
                
                // Retain all GAP variants
                for (auto& var : vars) {
                    if (var.getType() == VariantType::GAP) {
                        updated_vars.push_back(var);
                    }
                }
                
                // Process effective bases across all positions
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

    // 5. Write back corrected majority consensus string and switch mode
    consensus.setStoredString(std::move(new_cons));
    consensus.setUseStoredString(true);
}