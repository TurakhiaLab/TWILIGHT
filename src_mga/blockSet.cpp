
#include "block.hpp"
#include "phylogeny.hpp"


#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for_each.h>
#include <queue>
#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;


// =======================
// BlockSet Implementation
// =======================

int Alignment::countVariationsInRange(BlockSet* blockSet, int aln_start, int aln_end) {
    if (!blockSet) return 0;

    int total_vars = 0;

    std::vector<BlockID> linearized = blockSet->getLinearizeBlocks();
        
    int current_offset = 0;
    for (BlockID id : linearized) {
        auto block = blockSet->getBlock(id);
        if (!block) continue;
        int block_len = block->getConsensus().length();
        int b_start = current_offset;
        int b_end = current_offset + block_len;
        current_offset = b_end; // 更新下一個 Block 的起點
        // 1. 檢查這個 Block 是否與 Alignment 的區間有重疊
        if (b_end <= aln_start || b_start >= aln_end) {
            continue; // 完全沒有重疊，跳過
        }
        // 2. 計算重疊區域，並轉換為 Block 的「內部相對座標 (Local Coordinates)」
        int overlap_start = std::max(b_start, aln_start);
        int overlap_end   = std::min(b_end, aln_end);
        
        int local_start = overlap_start - b_start;
        int local_end   = overlap_end - b_start;
        // 3. 走訪這個 Block 內所有序列的 Segment，計算落在 local 區間內的 Variation
        for (auto& seq_pair : block->getSequences()) {
            for (auto& seg_pair : seq_pair.second.getSegments()) {
                auto& segment = seg_pair.second;
                
                for (auto& var : segment.getVariants()) {
                    int v_start = var.getStart();
                    int v_end = var.getEnd();
                    // 如果 Variation 的位置落在我們 Alignment 的重疊範圍內
                    if (v_end > local_start && v_start < local_end) {
                        total_vars++;
                    }
                }
            }
        }
    }
    return total_vars;
}

void Alignment::updateEnergy(BlockSet* refBlockSet, BlockSet* qryBlockSet, double beta) {
    int q_vars = countVariationsInRange(qryBlockSet, refIdx.first, refIdx.second);
    int r_vars = countVariationsInRange(refBlockSet, qryIdx.first, qryIdx.second);
    int total_variation_count = q_vars + r_vars;
        
    double aln_len = static_cast<double>(alnLength); 
        
    energy = -aln_len + beta * total_variation_count;
}

void BlockSet::rebuildAllPointers() {
    auto allBlocks = this->getAllBlocks(); // 假設回傳 std::vector<std::shared_ptr<Block>>
    if (allBlocks.empty()) return;

    // ==========================================
    // 1. 清理舊指標並收集所有 Segments (循序執行確保安全)
    // ==========================================
    // 這裡我們直接使用 shared_ptr，避免 weak_ptr 在轉換時造成的生命週期遺失
    struct SegRef { 
        Segment* seg; 
        std::shared_ptr<Block> blk; 
    };
    std::map<std::string, std::vector<SegRef>> seqTracks;

    for (auto& blk : allBlocks) {
        // [任務 A] 清除 Block 本身的舊指標
        blk.lock()->setPrevBlock(nullptr);
        blk.lock()->setNextBlock(nullptr);

        for (auto& seqPair : blk.lock()->getSequences()) {
            for (auto& segPairInner : seqPair.second.getSegments()) {
                // [任務 B] 清除 Segment 本身的舊指標
                segPairInner.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                segPairInner.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                
                // 將 Segment 的實體記憶體位置與其隸屬的 Block 綁定收集
                seqTracks[seqPair.first].push_back({ &segPairInner.second, blk.lock() });
            }
        }
    }

    // ==========================================
    // 2. 重新接線 Segment 層級的指標 (TBB 平行加速)
    // ==========================================
    std::vector<std::vector<SegRef>*> trackPtrs;
    trackPtrs.reserve(seqTracks.size());
    for (auto& trackPair : seqTracks) {
        trackPtrs.push_back(&trackPair.second);
    }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, trackPtrs.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto& track = *(trackPtrs[i]);
                
                // 依照絕對物理座標排序 (保證 0 -> N 的順序)
                std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                    return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
                });

                // 進行雙向綁定
                for (size_t j = 1; j < track.size(); ++j) {
                    auto& prevRef = track[j-1]; // 物理位置在左側
                    auto& currRef = track[j];   // 物理位置在右側

                    // 根據 Strand 決定進入與離開的方向
                    if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                    else prevRef.seg->setPrevBlock(currRef.blk); 

                    if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                    else currRef.seg->setNextBlock(prevRef.blk); 
                }
            }
        }
    );

    // ==========================================
    // 3. 建立 Block 層級的「嚴格線性」指標 (DFS-based Sequence Topological Sort)
    // ==========================================
    // 核心思想：只允許真實序列走過的路徑建立連線。利用 DFS 自然忽略環狀基因體的 back-edge。

    std::unordered_map<BlockID, std::vector<BlockID>> adj;
    std::unordered_map<BlockID, int> inDegree;
    std::unordered_map<BlockID, int> minCoord;
    std::unordered_map<BlockID, std::shared_ptr<Block>> idToBlock;

    for (auto& blk : allBlocks) {
        auto sharedBlk = blk.lock();
        if (sharedBlk) {
            BlockID id = sharedBlk->getId();
            idToBlock[id] = sharedBlk;
            inDegree[id] = 0;
            
            // 收集最小座標，僅用於在遇到「完全沒有頭」的純環狀基因體時，決定從哪裡剪開
            int mC = std::numeric_limits<int>::max();
            for (auto& seqPair : sharedBlk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    mC = std::min(mC, std::min(segPair.second.getStart(), segPair.second.getEnd()));
                }
            }
            minCoord[id] = mC;
        }
    }

    // 嚴格依據真實 sequence 軌跡建立有向圖
    for (auto& trackPair : seqTracks) {
        auto& track = trackPair.second;
        for (size_t j = 1; j < track.size(); ++j) {
            BlockID u = track[j-1].blk->getId(); // 離開
            BlockID v = track[j].blk->getId();   // 進入
            if (u != v) {
                // 確保不重複加入相同的邊
                if (std::find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
                    adj[u].push_back(v);
                    inDegree[v]++;
                }
            }
        }
    }

    // 準備所有的節點，並決定 DFS 的「起跑順序」
    std::vector<BlockID> nodes;
    nodes.reserve(idToBlock.size());
    for (auto& p : idToBlock) nodes.push_back(p.first);

    // 🚨 核心防護 1：排序起點！永遠從「入度最少 (最像源頭)」且「座標最小」的 Block 開始走
    std::sort(nodes.begin(), nodes.end(), [&](BlockID a, BlockID b) {
        if (inDegree[a] != inDegree[b]) return inDegree[a] < inDegree[b];
        if (minCoord[a] != minCoord[b]) return minCoord[a] < minCoord[b];
        return a < b;
    });

    // 🚨 核心防護 2：讓 DFS 在分岔路口時，永遠優先走向座標較小的鄰居
    for (auto& p : adj) {
        std::sort(p.second.begin(), p.second.end(), [&](BlockID a, BlockID b) {
            if (minCoord[a] != minCoord[b]) return minCoord[a] < minCoord[b];
            return a < b;
        });
    }

    std::unordered_map<BlockID, int> state; // 0=未走訪, 1=正在走訪(偵測環用), 2=已走完
    std::vector<BlockID> postOrder;

    // DFS 遞迴引擎
    std::function<void(BlockID)> dfs = [&](BlockID u) {
        state[u] = 1; // 標記為正在走訪
        for (BlockID v : adj[u]) {
            if (state[v] == 0) {
                dfs(v); // 繼續順著 Sequence 往下走
            }
            // 如果 state[v] == 1，代表這是一條「尾接頭」的繞回環狀邊 (Back-edge)，DFS 會優雅地忽略它！
        }
        state[u] = 2; // 走投無路，標記為走完
        postOrder.push_back(u); // 存入 Post-order
    };

    // 啟動 DFS
    for (BlockID u : nodes) {
        if (state[u] == 0) {
            dfs(u);
        }
    }

    // DFS 的 Post-order 反轉過來，就是完美的拓撲順序！
    std::reverse(postOrder.begin(), postOrder.end());

    // 根據最終完美的 DFS 順序，重新串接 1D Backbone
    for (size_t i = 0; i < postOrder.size(); ++i) {
        auto currBlk = idToBlock[postOrder[i]];
        if (i > 0) {
            currBlk->setPrevBlock(idToBlock[postOrder[i-1]]);
        } else {
            currBlk->setPrevBlock(nullptr); 
        }
        
        if (i < postOrder.size() - 1) {
            currBlk->setNextBlock(idToBlock[postOrder[i+1]]);
        } else {
            currBlk->setNextBlock(nullptr); 
        }
    }
}

void BlockSet::normalizeFamilyIDs() {
    bool debug = true;
    if (debug) std::cout << "\n[BlockSet] 🧬 Normalizing Family IDs...\n";

    std::unordered_map<FamilyID, FamilyID> oldToNewMap;
    FamilyID next_new_id = 1;
    
    std::unordered_map<FamilyID, BlockIDs> new_family_index;

    int updated_blocks = 0;

    for (auto& [blkId, blkPtr] : this->blocks) {
        if (!blkPtr) continue;

        FamilyID old_fam_id = blkPtr->getFamilyId(); 
        
        if (old_fam_id == 0) {
            new_family_index[0].push_back(blkId);
            continue;
        }

        if (oldToNewMap.find(old_fam_id) == oldToNewMap.end()) {
            oldToNewMap[old_fam_id] = next_new_id++;
        }

        FamilyID new_fam_id = oldToNewMap[old_fam_id];
        
        blkPtr->setFamilyId(new_fam_id); 
        updated_blocks++;

        new_family_index[new_fam_id].push_back(blkId);
    }

    this->family_index = std::move(new_family_index);

    if (debug) {
        std::cout << "  -> ✅ Normalized " << oldToNewMap.size() << " unique families.\n"
                  << "  -> 📊 Total blocks updated: " << updated_blocks << "\n"
                  << "  -> 🆔 Max Family ID is now: " << (next_new_id - 1) << "\n";
    }
}

BlockIDs BlockSet::getLinearizeBlocks() {

    if (is_cached) return linear_block_cache;

    BlockIDs linearized_order;

    BlockWeakPtr random_start = (this->getAllBlocks())[0];

    // Get starting block
    while (random_start.lock()->getPrevBlock().lock() != nullptr) {
        random_start = random_start.lock()->getPrevBlock().lock();
    }

    linearized_order.push_back(random_start.lock()->getId());
    
    while (random_start.lock()->getNextBlock().lock() != nullptr) {
        random_start = random_start.lock()->getNextBlock().lock();
        linearized_order.push_back(random_start.lock()->getId());
    }

    is_cached = true;
    this->linear_block_cache = linearized_order;
    return linearized_order;
}

BlockIDs BlockSet::getAncestralBlocks() {

    BlockIDs ancestral_blocks;

    for (auto& blkID: this->getLinearizeBlocks()) {
        auto blk = this->getBlock(blkID);
        if (blk->isDistant()) continue;
        ancestral_blocks.push_back(blk->getId());
    }

    return ancestral_blocks;
}


StringPairs BlockSet::getRemainingBlockConsensus() {
    StringPairs remaining;
    auto linearized_blocks = this->getLinearizeBlocks();
    for (const auto& id : linearized_blocks) {
        auto blk = this->getBlock(id);
        if (blk->isDistant()) {
            std::string block_id = this->ID + "_" + std::to_string(blk->getId());
            remaining.push_back({block_id, blk->getConsensus()});
        }
    }
    return remaining;
}

StringPairs BlockSet::getRepresentativeConsensus() {
    
    auto ancestral_blocks = this->getAncestralBlocks();
    StringPairs representative;

    std::string consensus = "";

    for (auto block_id: ancestral_blocks) {
        auto blk = this->getBlock(block_id);
        consensus += this->getBlock(block_id)->getConsensus();
    }
    representative.push_back({this->ID, consensus});
    return representative;
}

std::string BlockSet::reconstructSequence(const std::string& seqName) {
    struct SegNode {
        int start;
        int end;
        bool isRev;
        Variants vars;
        BlockPtr blk;
    };
    std::vector<SegNode> ordered_segments;

    // 1. 收集該序列散落在全圖的所有 Segments
    for (auto& blkPair : blocks) {
        std::shared_ptr<Block> blk = blkPair.second;
        auto& seqs = blk->getSequences();
        
        auto it = seqs.find(seqName);
        if (it != seqs.end()) {
            for (auto& segPair : it->second.getSegments()) {
                ordered_segments.push_back({
                    segPair.second.getStart(),
                    segPair.second.getEnd(),
                    segPair.second.isReverse(),
                    segPair.second.getVariants(),
                    blk
                });
            }
        }
    }

    if (ordered_segments.empty()) {
        std::cerr << "[Warning] Sequence '" << seqName << "' not found in BlockSet " << ID << ".\n";
        return "";
    }

    // 2. 依照基因體真實座標排序 (保證從 5' 端一路拼到 3' 端)
    std::sort(ordered_segments.begin(), ordered_segments.end(), [](const SegNode& a, const SegNode& b) {
        return std::min(a.start, a.end) < std::min(b.start, b.end);
    });

    // 3. 依序重建序列
    std::string reconstructedSeq = "";
    for (auto& node : ordered_segments) {
        std::string block_seq = node.blk->getConsensus();

        // 步驟 A: 應用 SNV
        for (auto& v : node.vars) {
            if (v.getType() == VariantType::SNV) {
                // 安全檢查，避免越界
                if (v.getStart() < block_seq.length()) {
                    block_seq[v.getStart()] = v.getAlt();
                }
            }
        }

        // 步驟 B: 應用 GAP
        std::string seg_seq = "";
        int cur = 0;
        
        Variants sorted_vars = node.vars;
        std::sort(sorted_vars.begin(), sorted_vars.end(), [](Variant& a, Variant& b) {
            return a.getStart() < b.getStart();
        });

        for (auto& v : sorted_vars) {
            if (v.getType() == VariantType::GAP) {
                if (v.getStart() > cur) {
                    seg_seq += block_seq.substr(cur, v.getStart() - cur);
                }
                cur = v.getEnd(); // 直接跳過 GAP
            }
        }
        if (cur < block_seq.length()) {
            seg_seq += block_seq.substr(cur); 
        }

        // 步驟 C: 如果是反股 (-)，進行 Reverse Complement
        if (node.isRev) {
            std::string rc = seg_seq;
            std::reverse(rc.begin(), rc.end());
            for (char& c : rc) {
                switch (c) {
                    case 'A': c = 'T'; break; case 'T': c = 'A'; break;
                    case 'C': c = 'G'; break; case 'G': c = 'C'; break;
                    case 'a': c = 't'; break; case 't': c = 'a'; break;
                    case 'c': c = 'g'; break; case 'g': c = 'c'; break;
                }
            }
            seg_seq = rc;
        }

        reconstructedSeq += seg_seq;
    }

    return reconstructedSeq;
}

BlockPtr BlockSet::concatenateBlocks(BlockID superId) {

    bool DEBUG_MODE = false;

    auto time0 = std::chrono::high_resolution_clock::now();
    // auto consensusBlocks = this->getRepresentativeBlocks(); 
    BlockIDs consensusBlocks;

    for (auto& blkID: this->getAncestralBlocks()) {
        auto blk = this->getBlock(blkID);
        consensusBlocks.push_back(blk->getId());
    }
    if (consensusBlocks.empty()) return nullptr;
    
    // ==========================================
    // 步驟 1: 預先建立 Consensus 全局尺與基礎資料
    // ==========================================
    std::string super_consensus = "";
    std::unordered_map<BlockID, int> block_super_offsets; // 記錄每個 Block 在 SuperBlock 的起點
    std::unordered_map<BlockID, int> block_lengths;       // 預存長度，避免在 TBB 內並行讀取 Block
    
    int current_offset = 0;
    for (auto blkID : consensusBlocks) {
        auto blk = this->getBlock(blkID);
        block_super_offsets[blkID] = current_offset;
        
        int blk_len = blk->getConsensus().length();
        block_lengths[blkID] = blk_len;
        
        super_consensus += blk->getConsensus();
        current_offset += blk_len;
    }
    int total_super_len = super_consensus.length();

    auto time1 = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) std::cout << "Step 1 (Offsets & Consensus) took: " 
                              << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << " ms\n";

    // ==========================================
    // 步驟 2: 建立 Inverted Index (SeqID -> Segments)
    // ==========================================
    struct SegEntry {
        BlockID blkID;
        Segment seg;
    };
    std::unordered_map<std::string, std::vector<SegEntry>> seq_to_segments;

    for (auto blkID : consensusBlocks) {
        auto blk = this->getBlock(blkID);
        for (auto& seqPair : blk->getSequences()) {
            const std::string& seqID = seqPair.first;
            for (auto& segPairInner : seqPair.second.getSegments()) {
                seq_to_segments[seqID].push_back({blkID, segPairInner.second});
            }
        }
    }

    // 將資料轉移到 vector 以便讓 TBB 進行一維平行切分
    std::vector<std::string> unique_seqs;
    std::vector<std::vector<SegEntry>*> seq_entries_ptrs;
    unique_seqs.reserve(seq_to_segments.size());
    seq_entries_ptrs.reserve(seq_to_segments.size());

    for (auto& kv : seq_to_segments) {
        unique_seqs.push_back(kv.first);
        seq_entries_ptrs.push_back(&kv.second);
    }

    struct Track {
        std::string seqID;
        Segment seg;
    };
    // TBB 平行寫入的絕對安全陣列
    std::vector<std::vector<Track>> final_tracks_results(unique_seqs.size());

    // ==========================================
    // 步驟 3: TBB 平行處理 (Sequence-Centric 核心)
    // ==========================================
    tbb::parallel_for(tbb::blocked_range<size_t>(0, unique_seqs.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const std::string& seqID = unique_seqs[i];
                std::vector<SegEntry>& entries = *(seq_entries_ptrs[i]);
                std::vector<Track>& final_tracks = final_tracks_results[i];

                if (entries.empty()) continue;

                // [核心防線]：根據生物真實座標(start)遞增排序，保證順著基因體方向(5'->3')組裝
                std::sort(entries.begin(), entries.end(), [](SegEntry& a, SegEntry& b) {
                    return a.seg.getStart() < b.seg.getStart();
                });

                Track current_track;
                current_track.seqID = seqID;
                bool has_active_track = false;
                
                // 正股使用 end 追蹤向右延伸；反股使用 start 追蹤向左延伸
                int current_track_super_end = 0;   
                int current_track_super_start = 0; 

                // --- Helper 1: 收尾並儲存目前的 Track ---
                auto close_active_track = [&]() {
                    auto& vars = current_track.seg.getVariants();
                    if (!current_track.seg.isReverse()) {
                        // 正股 Suffix Gap (補齊右邊)
                        if (current_track_super_end < total_super_len) {
                            if (!vars.empty() && vars.back().getType() == VariantType::GAP && vars.back().getEnd() == current_track_super_end) {
                                int old_start = vars.back().getStart();
                                vars.pop_back();
                                vars.push_back(Variant::createGap(old_start, total_super_len));
                            } else {
                                vars.push_back(Variant::createGap(current_track_super_end, total_super_len));
                            }
                        }
                    } else {
                        // 反股 Suffix Gap (生物學尾端 = 物理最左邊，補齊 0 到 start)
                        if (current_track_super_start > 0) {
                            if (!vars.empty() && vars.front().getType() == VariantType::GAP && vars.front().getStart() == current_track_super_start) {
                                int old_end = vars.front().getEnd();
                                vars.erase(vars.begin());
                                Variants new_vars;
                                new_vars.push_back(Variant::createGap(0, old_end));
                                new_vars.insert(new_vars.end(), std::make_move_iterator(vars.begin()), std::make_move_iterator(vars.end()));
                                vars = std::move(new_vars);
                            } else {
                                Variants new_vars;
                                new_vars.push_back(Variant::createGap(0, current_track_super_start));
                                new_vars.insert(new_vars.end(), std::make_move_iterator(vars.begin()), std::make_move_iterator(vars.end()));
                                vars = std::move(new_vars);
                            }
                        }
                    }
                    final_tracks.push_back(std::move(current_track));
                    has_active_track = false;
                };

                // --- Helper 2: 開啟全新的 Track ---
                auto open_new_track = [&](Segment seg, int offset, int len) {
                    current_track.seg = std::move(seg);
                    if (!current_track.seg.isReverse()) {
                        // 正股 Prefix Gap (補齊左邊)
                        if (offset > 0) {
                            auto& vars = current_track.seg.getVariants();
                            Variants new_vars;
                            new_vars.reserve(vars.size() + 1);
                            new_vars.push_back(Variant::createGap(0, offset));
                            new_vars.insert(new_vars.end(), std::make_move_iterator(vars.begin()), std::make_move_iterator(vars.end()));
                            vars = std::move(new_vars);
                        }
                        current_track_super_end = offset + len;
                    } else {
                        // 反股 Prefix Gap (生物學開頭 = 物理最右邊，補齊 end 到 total_len)
                        // 注意：反股座標是遞減的，所以 Prefix Gap 是在陣列的最尾端
                        if (offset + len < total_super_len) {
                            current_track.seg.getVariants().push_back(Variant::createGap(offset + len, total_super_len));
                        }
                        current_track_super_start = offset;
                    }
                    has_active_track = true;
                };

                // 開始遍歷 Segments
                for (size_t j = 0; j < entries.size(); ++j) {
                    SegEntry& entry = entries[j];
                    int block_offset = block_super_offsets.at(entry.blkID);
                    int block_len = block_lengths.at(entry.blkID);

                    Segment processed_seg = entry.seg; 
                    for (auto& var : processed_seg.getVariants()) {
                        var.shift(block_offset);
                    }

                    if (!has_active_track) {
                        open_new_track(std::move(processed_seg), block_offset, block_len);
                        continue;
                    } 

                    bool is_same_strand = (current_track.seg.isReverse() == processed_seg.isReverse());
                    bool can_merge = false;
                    
                    // 1. 生物學相鄰檢查 (無論正反股，既然已按 Start 排序，必定是前者的 End 接後者的 Start)
                    if (is_same_strand && current_track.seg.getEnd() == processed_seg.getStart()) {
                        can_merge = true;
                    }

                    // 2. 依照 Strand 進行物理座標防呆與合併
                    if (can_merge) {
                        if (!current_track.seg.isReverse()) {
                            // 【正股 Forward 合併邏輯】
                            if (block_offset >= current_track_super_end) {
                                current_track.seg.setEnd(processed_seg.getEnd());
                                auto& vars = current_track.seg.getVariants();
                                
                                if (block_offset > current_track_super_end) { // 處理中間 Gap
                                    if (!vars.empty() && vars.back().getType() == VariantType::GAP && vars.back().getEnd() == current_track_super_end) {
                                        int old_start = vars.back().getStart();
                                        vars.pop_back();
                                        vars.push_back(Variant::createGap(old_start, block_offset));
                                    } else {
                                        vars.push_back(Variant::createGap(current_track_super_end, block_offset));
                                    }
                                }
                                // 正股向右延伸：把新 Variations 加在尾巴
                                vars.insert(vars.end(), std::make_move_iterator(processed_seg.getVariants().begin()), std::make_move_iterator(processed_seg.getVariants().end()));
                                current_track_super_end = block_offset + block_len;
                            } else {
                                can_merge = false; // Duplication 往回跳，觸發斷點
                            }
                        } else {
                            // 【反股 Reverse 合併邏輯】(修復核心)
                            // 物理座標檢查：B 的右邊界 (offset+len) 必須 <= A 的左邊界 (current_super_start)
                            if (block_offset + block_len <= current_track_super_start) {
                                current_track.seg.setEnd(processed_seg.getEnd()); 
                                
                                auto& vars = current_track.seg.getVariants();
                                Variants new_vars;
                                new_vars.reserve(processed_seg.getVariants().size() + 1 + vars.size());
                                
                                // 反股向左延伸：B 的物理座標較小，Variations 必須插在最前面
                                new_vars.insert(new_vars.end(), std::make_move_iterator(processed_seg.getVariants().begin()), std::make_move_iterator(processed_seg.getVariants().end()));
                                
                                // 處理中間 Gap (如果 B 和 A 中間有跳過 Block)
                                if (block_offset + block_len < current_track_super_start) {
                                    if (!vars.empty() && vars.front().getType() == VariantType::GAP && vars.front().getStart() == current_track_super_start) {
                                        int old_end = vars.front().getEnd();
                                        vars.erase(vars.begin());
                                        new_vars.push_back(Variant::createGap(block_offset + block_len, old_end));
                                    } else {
                                        new_vars.push_back(Variant::createGap(block_offset + block_len, current_track_super_start));
                                    }
                                }
                                
                                // 把 A 原本的 Variations 接在後面
                                new_vars.insert(new_vars.end(), std::make_move_iterator(vars.begin()), std::make_move_iterator(vars.end()));
                                
                                vars = std::move(new_vars);
                                current_track_super_start = block_offset; // 更新最左邊界
                            } else {
                                can_merge = false; // Duplication 往右跳 (在反股等於往回跳)，觸發斷點
                            }
                        }
                    }

                    // 【情境 C：遇到 Breakpoint，無法合併】
                    if (!can_merge) {
                        close_active_track();
                        open_new_track(std::move(processed_seg), block_offset, block_len);
                    }
                }

                // 迴圈結束，收尾最後一個 Track
                if (has_active_track) {
                    close_active_track();
                }
            }
        }
    );

    // ==========================================
    // 步驟 4: 快速組裝 Final Super Block
    // ==========================================
    auto super_block = std::make_shared<Block>(superId, super_consensus);
    
    for (size_t i = 0; i < unique_seqs.size(); ++i) {
        if (final_tracks_results[i].empty()) continue;
        
        Sequence seq_info(unique_seqs[i]);
        for (auto& track : final_tracks_results[i]) {
            seq_info.addSegment(track.seg); // Segment 會自動以 start 為 Key 放入 map
        }
        super_block->addSequence(seq_info);
    }
    
    auto time2 = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) std::cout << "Step 3+4 (TBB Sequence Assembly) took: " 
                              << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << " ms\n";

    return super_block;
}

// ==========================================
// 將外部的 Block 加入此 BlockSet 並賦予新 ID
// ==========================================
BlockPtr BlockSet::addBlock(BlockPtr oldBlock, FamilyID familyId) {
    if (!oldBlock) return nullptr;

    BlockID newId = next_block_id_++;

    FamilyID assignedFamilyId = (familyId == 0) ? newId : familyId;

    // 3. 利用舊 Block 的 Consensus 建立全新的 Block
    auto newBlock = std::make_shared<Block>(newId, oldBlock->getConsensus());

    // 【新增】將決定的 Family ID 賦予給這個新 Block
    newBlock->setFamilyId(assignedFamilyId);

    // 4. 深拷貝：將所有的 SequenceInfo 複製過去
    for (const auto& seqPair : oldBlock->getSequences()) {
        newBlock->addSequence(seqPair.second);
    }

    // 5. 註冊進這個 BlockSet 的 Dictionary 中
    blocks[newId] = newBlock;

    // 6. 【新增】更新 Family Index 字典
    // 把這個新 Block 的 ID 加到它所屬家族的清單中
    family_index[assignedFamilyId].push_back(newId);

    invalidateRepCache();
    
    return newBlock; // 回傳新建立的 Block 智慧指標
}

std::pair<BlockID, BlockID> BlockSet::splitSingleBlock(int parentID, int localCut) {
    auto parent = this->getBlock(parentID);
    if (!parent) return {(uint64_t)-1, (uint64_t)-1};

    bool debug = false;

    // 1. 建立新 Block (左右半部)
    auto left = this->createBlock(parent->getConsensus().substr(0, localCut));
    auto right = this->createBlock(parent->getConsensus().substr(localCut));

    // ==========================================
    // 2. TBB 平行化：切割 Sequence 與 Segment
    // ==========================================
    std::vector<std::string> seqIDs;
    seqIDs.reserve(parent->getSequences().size());
    for (auto& kv : parent->getSequences()) {
        seqIDs.push_back(kv.first);
    }

    // 用來儲存平行切割結果的暫存結構，避免 Thread Contention
    struct SplitResult {
        Sequence leftSeq;
        Sequence rightSeq;
        bool hasLeft = false;
        bool hasRight = false;
    };
    std::vector<SplitResult> splitResults(seqIDs.size());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, seqIDs.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const std::string& seqID = seqIDs[i];
                auto& parentSeqInfo = parent->getSequences().at(seqID);

                Sequence leftSeqInfo(seqID);
                Sequence rightSeqInfo(seqID);

                for (auto& segPair : parentSeqInfo.getSegments()) {
                    Segment oldSeg = segPair.second; // 拷貝出來處理
                    auto splitSegs = oldSeg.split(localCut);
                    Segment& leftSeg = splitSegs.first;
                    Segment& rightSeg = splitSegs.second;

                    // 1. 先判斷這個 Segment 切出來後，是否真實擁有物理序列 (非純 Gap)
                    bool validLeft = (leftSeg.getStart() != leftSeg.getEnd());
                    bool validRight = (rightSeg.getStart() != rightSeg.getEnd());

                    // 2. 內部接線：嚴格遵守正反股走向，且「只對真正存在的 Segment 接線」
                    if (validLeft && validRight) {
                        // 兩邊都有肉：互相連接，並對外連接
                        if (!oldSeg.isReverse()) {
                            leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                            leftSeg.setNextBlock(right);
                            rightSeg.setPrevBlock(left);
                            rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
                        } else {
                            rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                            rightSeg.setNextBlock(left);
                            leftSeg.setPrevBlock(right);
                            leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
                        }
                    } 
                    else if (validLeft && !validRight) {
                        // 只有左邊有肉：左邊直接繼承原 Segment 的所有對外連接
                        leftSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                        leftSeg.setNextBlock(oldSeg.getNextBlock().lock());
                    } 
                    else if (!validLeft && validRight) {
                        // 只有右邊有肉：右邊直接繼承原 Segment 的所有對外連接
                        rightSeg.setPrevBlock(oldSeg.getPrevBlock().lock());
                        rightSeg.setNextBlock(oldSeg.getNextBlock().lock());
                    }
                    // 如果兩邊都沒肉 (!validLeft && !validRight)，那就什麼都不用接，直接丟棄

                    // 3. 將真實存在的 Segment 放入 Map 裡
                    if (validLeft) {
                        leftSeqInfo.getSegments()[leftSeg.getStart()] = leftSeg;
                    }
                    if (validRight) {
                        rightSeqInfo.getSegments()[rightSeg.getStart()] = rightSeg;
                    }
                }

                // 將結果存入專屬的 index，確保 Thread Safe
                splitResults[i].leftSeq = std::move(leftSeqInfo);
                splitResults[i].rightSeq = std::move(rightSeqInfo);
                splitResults[i].hasLeft = !splitResults[i].leftSeq.getSegments().empty();
                splitResults[i].hasRight = !splitResults[i].rightSeq.getSegments().empty();
            }
        }
    );

    // 主執行緒快速合併結果 (將 Map 搬進 left/right block)
    for (size_t i = 0; i < seqIDs.size(); ++i) {
        if (splitResults[i].hasLeft) left->addSequence(std::move(splitResults[i].leftSeq));
        if (splitResults[i].hasRight) right->addSequence(std::move(splitResults[i].rightSeq));
    }

    // ==========================================
    // 3. 【修復核心】：消滅全圖掃描，改為「鄰居局部掃描」+ TBB 平行接線
    // ==========================================
    std::unordered_set<std::shared_ptr<Block>> neighbors;
    
    // 必須加入 left 和 right 來解開 Self-loop
    neighbors.insert(left);
    neighbors.insert(right);

    // 收集真正有牽連的鄰居 (只看 parent 原本的連線)
    for (auto& seqPair : parent->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            if (auto p = segPair.second.getPrevBlock().lock()) neighbors.insert(p);
            if (auto n = segPair.second.getNextBlock().lock()) neighbors.insert(n);
        }
    }
    neighbors.erase(parent); // parent 即將被刪除，不用幫它接線

    // 將鄰居轉為 Vector 以供 TBB 平行處理
    std::vector<std::shared_ptr<Block>> neighbor_vec(neighbors.begin(), neighbors.end());
    std::atomic<int> rewiredCount{0};

    // TBB 平行接線：因為每條 Thread 處理不同的 Neighbor Block，
    // 其內部的 Segment 也是獨立的，因此絕對 Thread Safe！
    tbb::parallel_for(tbb::blocked_range<size_t>(0, neighbor_vec.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            int local_rewired = 0;
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto currentBlock = neighbor_vec[i];

                for (auto& seqPair : currentBlock->getSequences()) {
                    std::string seqID = seqPair.first;
                    for (auto& segPair : seqPair.second.getSegments()) {
                        Segment& seg = segPair.second;

                        // 檢查 Prev：如果這段序列是從 parent 來的
                        if (seg.getPrevBlock().lock() == parent) {
                            bool found = false;
                            if (left->getSequences().count(seqID)) {
                                for (auto& lSeg : left->getSequences().at(seqID).getSegments()) {
                                    if (lSeg.second.getEnd() == seg.getStart()) {
                                        seg.setPrevBlock(left);
                                        found = true; local_rewired++; break;
                                    }
                                }
                            }
                            if (!found && right->getSequences().count(seqID)) {
                                for (auto& rSeg : right->getSequences().at(seqID).getSegments()) {
                                    if (rSeg.second.getEnd() == seg.getStart()) {
                                        seg.setPrevBlock(right);
                                        local_rewired++; break;
                                    }
                                }
                            }
                        }

                        // 檢查 Next：如果這段序列下一步要走到 parent
                        if (seg.getNextBlock().lock() == parent) {
                            bool found = false;
                            if (left->getSequences().count(seqID)) {
                                for (auto& lSeg : left->getSequences().at(seqID).getSegments()) {
                                    if (lSeg.second.getStart() == seg.getEnd()) {
                                        seg.setNextBlock(left);
                                        found = true; local_rewired++; break;
                                    }
                                }
                            }
                            if (!found && right->getSequences().count(seqID)) {
                                for (auto& rSeg : right->getSequences().at(seqID).getSegments()) {
                                    if (rSeg.second.getStart() == seg.getEnd()) {
                                        seg.setNextBlock(right);
                                        local_rewired++; break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            rewiredCount += local_rewired; // Atomic 累加
        }
    );

    if (debug) std::cout << "[DEBUG-SPLIT] Block " << parentID << " (Len: " << parent->getConsensus().size() << ") cut at " << localCut 
                         << " -> L: " << left->getId() << " (Len: " << left->getConsensus().size() << "), R: " << right->getId() << " (Len: " << right->getConsensus().size() 
                         << " | Rewired pointers: " << rewiredCount.load() << "\n";

    // 4. 安全刪除舊 Block
    this->deleteBlock(parent->getId());
    
    return {left->getId(), right->getId()};
}

void BlockSet::print(std::ostream& os) const {
    os << "\n==========================================================================\n";
    os << " 🌐 BLOCK SET ID: " << ID << " | Total Blocks In Map: " << blocks.size() << "\n";
    os << "==========================================================================\n";
    
    // 因為 getLinearizeBlocks() 在宣告中是非 const，這裡我們透過 const_cast 來安全調用
    auto& mutableSet = const_cast<BlockSet&>(*this);
    BlockIDs linearBlocks = mutableSet.getLinearizeBlocks();
    
    if (linearBlocks.empty()) {
        os << "  ⚠️  [Warning] Graph is empty or contains no linearized backbone blocks.\n";
        os << "==========================================================================\n\n";
        return;
    }
    
    // 依照 Linearized 順序逐一印出 Block
    for (size_t i = 0; i < linearBlocks.size(); ++i) {
        BlockID bid = linearBlocks[i];
        auto blk = mutableSet.getBlock(bid);
        
        if (blk) {
            blk->print(os);
            // 如果後面還有 Block，印出一個漂亮的拓撲流向箭頭
            if (i + 1 < linearBlocks.size()) {
                os << "                                   │\n";
                os << "                                   ▼\n";
            }
        } else {
            os << "  ❌ [ERROR] Block ID " << bid << " listed in linear backbone but missing from map!\n";
        }
    }
    os << "==========================================================================\n\n";
}

void BlockSet::setDistantBlocks(Tree& tree, int lookdownDepth) {
    bool debug = true;
    std::string targetNodeId = this->getId(); 

    auto it = tree.allNodes.find(targetNodeId);
    if (it == tree.allNodes.end()) return;
    Node* currentNode = it->second;

    if (currentNode->is_leaf() || currentNode->children.size() < 2) return; 

    std::vector<std::unordered_set<std::string>> sequenceSets;
    tree.getSubLineages(currentNode, lookdownDepth, 0, sequenceSets);

    uint64_t minSequenceSets = 1ULL << lookdownDepth;
    if (sequenceSets.size() < minSequenceSets) return;

    auto all_blocks = this->getAllBlocks();
    std::atomic<int> distantCount{0};
    std::atomic<int> coreCount{0};

    tbb::parallel_for_each(all_blocks.begin(), all_blocks.end(), [&](const auto& weak_blk) {
        auto blk = weak_blk.lock();
        if (!blk) return; // 🚨 注意：在 Lambda 裡，continue 要改成 return

        std::unordered_set<std::string> blockSeqNames;
        for (const auto& seqPair : blk->getSequences()) {
            blockSeqNames.insert(seqPair.first);
        }

        int supportedSets = 0;
        for (const auto& leafSet : sequenceSets) {
            for (const auto& seqName : blockSeqNames) {
                if (leafSet.find(seqName) != leafSet.end()) {
                    supportedSets++;
                    break; // 這個 Set 投下支持票，跳出檢查下一個 Set
                }
            }
        }

        if (supportedSets >= 2) {
            blk->setDistant(false); 
            coreCount++;    // std::atomic 支援直接安全的 ++
        } else {
            blk->setDistant(true);  
            distantCount++; // std::atomic 支援直接安全的 ++
        }
    });

    if (debug) { // DEBUG_MODE
        std::cout << "[INFO] Node " << targetNodeId 
                  << " (Lookdown " << lookdownDepth << " levels -> " << sequenceSets.size() << " sets)"
                  << " | Core=" << coreCount << ", Distant=" << distantCount << "\n";
    }
}

BlockBoundaries BlockSet::extractBlockBoundaries() {
    BlockBoundaries boundaries;
    bool debug = false; 

    BlockIDs consensusBlocks;

    for (auto& blkID: this->getAncestralBlocks()) {
        auto blk = this->getBlock(blkID);
        consensusBlocks.push_back(blk->getId());
    }
    auto block_id = std::move(consensusBlocks);
    if (block_id.size() < 2) return boundaries;

    int currentAbsolutePos = 0; 

    for (size_t i = 0; i < block_id.size() - 1; ++i) {
        int leftID = block_id[i];
        int rightID = block_id[i + 1];
        
        auto leftBlk = blocks[leftID];
        auto rightBlk = blocks[rightID];
        
        if (!leftBlk || !rightBlk) {
            if (leftBlk) currentAbsolutePos += leftBlk->getConsensus().length();
            continue;
        }

        currentAbsolutePos += leftBlk->getConsensus().length();

        BlockBoundary bnd;
        bnd.leftBlockId = leftBlk->getId();
        bnd.rightBlockId = rightBlk->getId();
        bnd.leftConsensusEndPos = currentAbsolutePos; 
        bnd.type = BoundaryType::FLEXIBLE; 
        bnd.reason = "FLEXIBLE Match"; 

        // 🌟 寫入左右積木的 Family ID
        bnd.leftFamilyId = leftBlk->getFamilyId();
        bnd.rightFamilyId = rightBlk->getFamilyId();

        bool isStrict = false;
        bool cannot_push_left = false;  // Sequence 存在於 left，不存在於 right
        bool cannot_push_right = false; // Sequence 存在於 right，不存在於 left
        bool has_intersection = false;

        auto& leftSeqs = leftBlk->getSequences();
        auto& rightSeqs = rightBlk->getSequences();

        for (auto& leftSeqPair : leftSeqs) {
            const std::string& seqName = leftSeqPair.first;
            auto rightSeqIt = rightSeqs.find(seqName);
            
            if (rightSeqIt != rightSeqs.end()) {
                has_intersection = true;
                
                auto& leftSegs = leftSeqPair.second.getSegments();
                auto& rightSegs = rightSeqIt->second.getSegments();
                
                if (!leftSegs.empty() && !rightSegs.empty()) {
                    auto leftSeg = leftSegs.rbegin()->second;
                    auto rightSeg = rightSegs.begin()->second;
                    
                    bool leftStrand = leftSeg.isReverse();
                    bool rightStrand = rightSeg.isReverse();
                    
                    // ==========================================
                    // 規則 1: Strict (Strand 相反)
                    // ==========================================
                    if (leftStrand != rightStrand) {
                        isStrict = true;
                        bnd.reason = "Strand Inversion on " + seqName;
                        break;
                    }
                    
                    // ==========================================
                    // 規則 2b: Strict (Order 沒有對起來)
                    // 根據絕對座標邏輯：start 恆小於 end
                    // ==========================================
                    bool isContiguous = false;
                    if (!leftStrand) {
                        // 正向: 左積木的尾巴 必須接上 右積木的頭
                        if (leftSeg.getEnd() == rightSeg.getStart()) isContiguous = true;
                    } else {
                        // 反向: 序列從右向左讀，所以左積木的頭 必須接上 右積木的尾巴
                        if (leftSeg.getStart() == rightSeg.getEnd()) isContiguous = true;
                    }
                    
                    if (!isContiguous) {
                        isStrict = true;
                        bnd.reason = "Coordinate Discontinuity on " + seqName;
                        break;
                    }
                }
            } else {
                // ==========================================
                // 規則 4: 在 Left 不在 Right -> 不能 push left (會吃掉特有序列)
                // ==========================================
                cannot_push_left = true;
            }
        }

        // 如果在前面的檢查中已經觸發 STRICT，直接寫入並跳到下一個迴圈
        if (isStrict) {
            bnd.type = BoundaryType::STRICT;
            boundaries[currentAbsolutePos] = bnd; 
            continue;
        }

        // 走訪 RightSeqs，檢查是否存在於 Right 但不存在於 Left
        for (const auto& rightSeqPair : rightSeqs) {
            if (leftSeqs.find(rightSeqPair.first) == leftSeqs.end()) {
                // ==========================================
                // 規則 3: 在 Right 不在 Left -> 不能 push right
                // ==========================================
                cannot_push_right = true;
                break; // 只要找到一條符合的就可以提早結束
            }
        }

        // ==========================================
        // 綜合判定最終的邊界屬性
        // ==========================================
        if (!has_intersection && (!leftSeqs.empty() || !rightSeqs.empty())) {
            // 規則 2a: 完全沒有交集
            bnd.type = BoundaryType::STRICT;
            bnd.reason = "No Sequence Intersection";
        } 
        else if (cannot_push_left && cannot_push_right) {
            // 規則 5: 同時不能 push right 也不能 push left
            bnd.type = BoundaryType::STRICT;
            bnd.reason = "Bidirectional Exclusivity (Cannot Push L/R)";
        } 
        else if (cannot_push_left) {
            // 不能 push left，表示只能向右推
            bnd.type = BoundaryType::PUSH_RIGHT_ONLY;
            bnd.reason = "Left exclusive seq (Cannot Push Left)";
        } 
        else if (cannot_push_right) {
            // 不能 push right，表示只能向左推
            bnd.type = BoundaryType::PUSH_LEFT_ONLY;
            bnd.reason = "Right exclusive seq (Cannot Push Right)";
        }

        // 如果依然是 FLEXIBLE 且帶有 Family 特徵，給它客製化的 Reason 方便 Debug
        if (bnd.isFamilySeam() && bnd.type == BoundaryType::FLEXIBLE) {
            bnd.reason = "FLEXIBLE Family Boundary";
            if (bnd.isHomoLeft())  bnd.reason += " [Homo_Left Fam:" + std::to_string(bnd.leftFamilyId) + "]";
            if (bnd.isHomoRight()) bnd.reason += " [Homo_Right Fam:" + std::to_string(bnd.rightFamilyId) + "]";
        }

        boundaries[currentAbsolutePos] = bnd; 
    }

    if (debug) {
        std::cout << "\n========================================================\n";
        std::cout << "=== EXTRACTED BLOCK BOUNDARIES (" << boundaries.size() << " total) ===\n";
        std::cout << "========================================================\n";
        for (auto& pair : boundaries) {
            std::cout << "[AbsPos: " << pair.first << "] ";
            
            std::string leftStr = std::to_string(pair.second.leftBlockId) + 
                                  (pair.second.leftFamilyId != 0 ? "(Fam:" + std::to_string(pair.second.leftFamilyId) + ")" : "");
            std::string rightStr = std::to_string(pair.second.rightBlockId) + 
                                   (pair.second.rightFamilyId != 0 ? "(Fam:" + std::to_string(pair.second.rightFamilyId) + ")" : "");
            
            std::string typeStr = (pair.second.type == BoundaryType::STRICT) ? "STRICT" : 
                                  (pair.second.type == BoundaryType::FLEXIBLE) ? "FLEXIBLE" : 
                                  (pair.second.type == BoundaryType::PUSH_LEFT_ONLY) ? "PUSH_LEFT_ONLY" : "PUSH_RIGHT_ONLY";
                                  
            std::cout << leftStr << " -> " << rightStr 
                      << " | Type: " << typeStr << " | Reason: " << pair.second.reason << "\n";
        }
        std::cout << "========================================================\n\n";
    }

    return boundaries;
}

bool BlockSet::realignBlock(BlockID blkId, std::string tempDir, int iterations) {

    bool DEBUG_MODE = true; // 🌟 幫你預設開啟，可以隨時關掉
    bool showImproved = true;
    bool alignAll = false;
    const int minLength = 100;
    const int maxLength = 50000;
    const float minIdentity = 99.95;
    const int minDepth = 8;


    struct SegMeta {
        std::string orig_seq_id;
        int start_coord;
        int end_coord;
        bool is_reverse;
    };

    std::string dipper_path = "/home/y3tseng@AD.UCSD.EDU/DIPPER/bin/dipper";
    std::string twilight_path_long = "/home/y3tseng@AD.UCSD.EDU/TWILIGHT/bin/twilight";
    std::string twilight_path_short = "/home/y3tseng@AD.UCSD.EDU/TWILIGHT_consistency/TWILIGHT/bin/twilight";
    

    auto blk = this->getBlock(blkId);
    if (!blk) {
        std::cerr << "[Error] Block " << blkId << " not found in BlockSet " << this->getId() << "\n";
        return false;
    }

    std::string oldConsensus = blk->getConsensus();
    int oldLen = oldConsensus.length();

    std::string twilight_path = (oldLen > 5000) ? twilight_path_long : twilight_path_short;
    
    // ==========================================
    // 步驟 1: 計算當前的 Length 與 Alignment Identity
    // ==========================================
    int total_bases = 0;
    int mismatch_bases = 0;
    int totalSegs = 0;

    for (auto& seqPair : blk->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            total_bases += oldLen; 
            for (auto& var : seg.getVariants()) {
                if (var.getType() == VariantType::SNV) {
                    mismatch_bases += 1;
                } else if (var.getType() == VariantType::GAP) {
                    mismatch_bases += (var.getEnd() - var.getStart());
                }
            }
            ++totalSegs;
        }
    }
    
    double oldIdentity = total_bases > 0 ? 100.0 * (1.0 - (double)mismatch_bases / total_bases) : 0.0;

    // 🌟 需求 1：Identity >= 99% 就直接跳過，不浪費算力！
    if (!alignAll) {
        if (oldIdentity >= minIdentity || (oldLen < minLength || oldLen > maxLength) || totalSegs < minDepth) {
            if (DEBUG_MODE) {
                std::ostringstream oss;
                oss << "[Realign Skip] Block " << blkId << " Identity is already " << oldIdentity << "%\n";
                std::cout << oss.str();
            }
            return false;
        }
    }

    // ==========================================
    // 步驟 2: Reconstruct Raw Sequences & Lookup Table
    // ==========================================
    std::map<std::string, SegMeta> lookupTable;
    std::string blockDir = tempDir + "/" + this->getId() + "_blk" + std::to_string(blkId);
    fs::create_directories(blockDir);

    std::string initFastaName = blockDir + "/init.fa";
    std::ofstream outFasta(initFastaName);

    if (!outFasta.is_open()) {
        std::cerr << "[Error] Cannot write to temporary FASTA file: " << initFastaName << "\n";
        return false;
    }

    for (auto& seqPair : blk->getSequences()) {
        std::string seqName = seqPair.first;
        int fragmentIndex = 1;
        
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            std::string rawSeq = "";
            int cons_pos = 0;
            
            Variants sortedVars = seg.getVariants();
            std::sort(sortedVars.begin(), sortedVars.end(), [](Variant& a, Variant& b) {
                return a.getStart() < b.getStart();
            });

            for (auto& var : sortedVars) {
                if (var.getStart() > cons_pos) {
                    rawSeq += oldConsensus.substr(cons_pos, var.getStart() - cons_pos);
                }
                if (var.getType() == VariantType::SNV) {
                    rawSeq += var.getAlt();
                    cons_pos = var.getStart() + 1;
                } else if (var.getType() == VariantType::GAP) {
                    cons_pos = var.getEnd(); 
                }
            }
            if (cons_pos < oldLen) {
                rawSeq += oldConsensus.substr(cons_pos);
            }

            std::string headerName = seqName + "." + std::to_string(fragmentIndex++);
            lookupTable[headerName] = {seqName, seg.getStart(), seg.getEnd(), seg.isReverse()};

            outFasta << ">" << headerName << "\n" << rawSeq << "\n";
        }
    }
    outFasta.close();

    // ==========================================
    // 步驟 3: 執行 DIPPER + TWILIGHT
    // ==========================================
    std::string currentFasta = initFastaName;
    for (int i = 1; i <= iterations; ++i) {
        std::string msaOut  = blockDir + "/msa_iter" + std::to_string(i) + ".fa";
        std::string treeOut = blockDir + "/tree_iter" + std::to_string(i) + ".nwk";
        
        std::string dipperCmd = (i == 1) ?    
                    dipper_path + " -i r -o t -I " + currentFasta + " -O " + treeOut + " > /dev/null 2>&1":
                    dipper_path + " -i m -o t -I " + currentFasta + " -O " + treeOut + " > /dev/null 2>&1";
        std::string twilightCmd = twilight_path + " --cpu-only -C 8 -i " + initFastaName + " -t " + treeOut + " -o " + msaOut + " > /dev/null 2>&1";
        // std::string twilightCmd = twilight_path + " -C 8 --cpu-only -i " + initFastaName + " -t " + treeOut + " -o " + msaOut;
        

        int ret1 = system(dipperCmd.c_str());
        int ret2 = system(twilightCmd.c_str());
        
        if (ret1 != 0 || ret2 != 0) {
            std::cerr << "[Warning] Block " << blkId << ": MSA tools returned non-zero at iteration " << i << ".\n";
            exit(1);
        }
        
        currentFasta = msaOut; 
    }

    // ==========================================
    // 步驟 4: 解析最終 MSA、重新計算 Consensus 與 Variations
    // ==========================================
    std::ifstream inMsa(currentFasta);
    if (!inMsa.is_open()) {
        std::cerr << "[Error] Cannot open final MSA file: " << currentFasta << " for Block " << blkId << "\n";
        fs::remove_all(blockDir);
        return false;
    }

    std::map<std::string, std::string> alignedSeqs;
    std::string line, curHeader = "";
    while (std::getline(inMsa, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (line[0] == '>') {
            curHeader = line.substr(1);
        } else {
            alignedSeqs[curHeader] += line;
        }
    }
    inMsa.close();

    if (alignedSeqs.empty()) {
        fs::remove_all(blockDir);
        return false;
    }

    int msaLen = alignedSeqs.begin()->second.length();
    std::string newConsensus = "";

    for (int i = 0; i < msaLen; ++i) {
        int counts[256] = {0};
        for (const auto& pair : alignedSeqs) {
            char base = pair.second[i];
            if (base != '-') counts[(unsigned char)base]++;
        }
        
        char bestBase = 'A'; 
        int maxCount = -1;
        for (int b = 0; b < 256; ++b) {
            if (b != '-' && counts[b] > maxCount) {
                maxCount = counts[b];
                bestBase = (char)b;
            }
        }
        newConsensus += bestBase;
    }

    std::unordered_map<std::string, Sequence> newSeqMap;
    int new_mismatch_bases = 0;
    int new_total_bases = 0;

    auto oldSeqs = blk->getSequences();

    for (const auto& pair : alignedSeqs) {
        std::string header = pair.first;
        std::string aln = pair.second;
        SegMeta meta = lookupTable[header];

        if (newSeqMap.find(meta.orig_seq_id) == newSeqMap.end()) {
            newSeqMap[meta.orig_seq_id] = Sequence(meta.orig_seq_id);
        }

        Segment newSeg = oldSeqs.at(meta.orig_seq_id).getSegments().at(meta.start_coord);
        
        std::vector<Variant> newVars;
        int gapStart = -1;

        for (int i = 0; i < msaLen; ++i) {
            new_total_bases++;
            if (aln[i] == '-') {
                if (gapStart == -1) gapStart = i;
                new_mismatch_bases++;
            } else {
                if (gapStart != -1) { 
                    newVars.push_back(Variant::createGap(gapStart, i));
                    gapStart = -1;
                }
                if (aln[i] != newConsensus[i]) { 
                    newVars.push_back(Variant(i, aln[i]));
                    new_mismatch_bases++;
                }
            }
        }
        if (gapStart != -1) {
            newVars.push_back(Variant::createGap(gapStart, msaLen));
        }

        newSeg.getVariants() = std::move(newVars);
        newSeqMap[meta.orig_seq_id].getSegments()[meta.start_coord] = newSeg;
    }

    double newIdentity = new_total_bases > 0 ? 100.0 * (1.0 - (double)new_mismatch_bases / new_total_bases) : 0.0;
    
    // 🌟 需求 2：有進步才覆蓋，沒進步就直接丟掉，保持原本的狀態
    if (newIdentity > oldIdentity) {
        blk->setConsensus(newConsensus);
        blk->setSequences(newSeqMap);

        // if (DEBUG_MODE) {
        if (showImproved) {
            std::ostringstream oss;
            oss << "[Realign ✅] Block " << blkId << " Improved! " 
                << oldLen << "bp -> " << msaLen << "bp | Id: " 
                << oldIdentity << "% -> " << newIdentity << "%\n";
            std::cout << oss.str();
        }
    } else {
        if (DEBUG_MODE || alignAll) {
            std::ostringstream oss;
            oss << "[Realign ❌] Block " << blkId << " No Improvement. "
                << "Id: " << oldIdentity << "% -> " << newIdentity << "%. Reverting changes.\n";
            std::cout << oss.str();
        }
        fs::remove_all(blockDir);
        return false;
    }
    
    // 確保暫存檔案被清理
    fs::remove_all(blockDir);
    return true;
}

void BlockSet::realignBlocks(std::string tempDir) {

    bool DEBUG_MODE = true;
    
    int currentBlock = 0, totalBlock = getLinearizeBlocks().size();

    // for (size_t i = 0; i < all_blocks.size(); ++i) {
    for (size_t i = 0; i < totalBlock; ++i) {
        bool improved = this->realignBlock(this->getBlock(linear_block_cache[i])->getId(), tempDir);
        ++currentBlock;
        if (currentBlock % 100 == 0) std::cout << '[' << currentBlock << '/' << totalBlock << "]\n";
    }
}


/*
std::vector<Block::ID> BlockSet::splitMultiBlocks(Block::ID parentID, const std::vector<int>& cuts) {
    auto parent = this->getBlock(parentID);
    if (!parent || cuts.empty()) return {parentID};

    auto time0 = std::chrono::high_resolution_clock::now();
    bool debug = false;


    // 1. 過濾與排序切點
    std::vector<int> validCuts;
    int consLen = parent->getConsensus().length();
    
    std::vector<int> sortedCuts = cuts;
    std::sort(sortedCuts.begin(), sortedCuts.end());
    
    for (int c : sortedCuts) {
        if (c > 0 && c < consLen && (validCuts.empty() || c != validCuts.back())) {
            validCuts.push_back(c);
        }
    }
    if (validCuts.empty()) return {parentID};

    auto time1 = std::chrono::high_resolution_clock::now();

    // 2. 建立新 Blocks (K 刀產生 K+1 塊)
    std::vector<std::shared_ptr<Block>> newBlocks;
    int prevCut = 0;
    // std::cout << "Length: " <<consLen << " Cuts: ";
    for (int cut : validCuts) {
        // std::cout << cut << ",";
        newBlocks.push_back(this->createBlock(parent->getConsensus().substr(prevCut, cut - prevCut)));
        prevCut = cut;
    }
    // std::cout << '\n';
    newBlocks.push_back(this->createBlock(parent->getConsensus().substr(prevCut)));
    
    auto time2 = std::chrono::high_resolution_clock::now();
    // ==========================================
    // 3. TBB 平行切割 Sequence 與 Segment
    // ==========================================
    std::vector<std::string> seqIDs;
    seqIDs.reserve(parent->getSequences().size());
    for (auto& kv : parent->getSequences()) seqIDs.push_back(kv.first);

    struct SplitResult {
        std::vector<SequenceInfo> chunkSeqs;
    };
    std::vector<SplitResult> splitResults(seqIDs.size());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, seqIDs.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const std::string& seqID = seqIDs[i];
                auto& parentSeqInfo = parent->getSequences().at(seqID);

                // 預先分配記憶體 (避免 chunkSeqs 內部動態擴容)
                std::vector<SequenceInfo> chunkSeqs(newBlocks.size(), SequenceInfo(seqID));

                for (auto& segPair : parentSeqInfo.getSegments()) {
                    
                    // 🚀 最佳化 1：直接把舊 Segment 的靈魂抽出來，不拷貝！
                    Segment currentSeg = std::move(segPair.second); 
                    
                    std::vector<Segment> subSegs;
                    subSegs.reserve(newBlocks.size()); // 預先分配，避免 vector 擴容鎖

                    int currentOffset = 0;
                    for (int cut : validCuts) {
                        int relativeCut = cut - currentOffset; 
                        auto splitSegs = currentSeg.split(relativeCut);
                        
                        // 🚀 最佳化 2：強制使用 std::move 轉移所有權
                        subSegs.push_back(std::move(splitSegs.first));  
                        currentSeg = std::move(splitSegs.second);       
                        currentOffset = cut;
                    }
                    subSegs.push_back(std::move(currentSeg)); 

                    // 內部接線
                    for (size_t j = 0; j < subSegs.size(); ++j) {
                        Segment& seg = subSegs[j]; // 這裡用 Reference 就好，不要拷貝

                        std::shared_ptr<Block> prevPtr;
                        std::shared_ptr<Block> nextPtr;

                        if (!seg.isReverse()) { // 注意：舊的 segPair.second 已經被掏空，這裡改看 seg
                            prevPtr = (j == 0) ? seg.getPrevBlock().lock() : newBlocks[j - 1];
                            nextPtr = (j == subSegs.size() - 1) ? seg.getNextBlock().lock() : newBlocks[j + 1];
                        } else {
                            prevPtr = (j == subSegs.size() - 1) ? seg.getPrevBlock().lock() : newBlocks[j + 1];
                            nextPtr = (j == 0) ? seg.getNextBlock().lock() : newBlocks[j - 1];
                        }

                        seg.setPrevBlock(prevPtr);
                        seg.setNextBlock(nextPtr);

                        if (seg.getStart() != seg.getEnd()) {
                            // 🚀 最佳化 3：使用 emplace 直接在 Map 內部建構，結合 std::move 達成完美零拷貝！
                            chunkSeqs[j].getSegments().emplace(seg.getStart(), std::move(seg));
                        }
                    }
                }
                splitResults[i].chunkSeqs = std::move(chunkSeqs);
            }
        }
    );
    auto time3 = std::chrono::high_resolution_clock::now();

    // std::cout << "seqs: " << seqIDs.size() << " Splits: " << validCuts.size() << "\t";
    // std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << "ms\n";

    
    // 主執行緒快速將暫存結果推入新的 Blocks
    for (size_t i = 0; i < seqIDs.size(); ++i) {
        for (size_t bIdx = 0; bIdx < newBlocks.size(); ++bIdx) {
            if (!splitResults[i].chunkSeqs[bIdx].getSegments().empty()) {
                newBlocks[bIdx]->addSequence(std::move(splitResults[i].chunkSeqs[bIdx]));
            }
        }
    }

    auto time4 = std::chrono::high_resolution_clock::now();

    // ==========================================
    // 4. TBB 平行外部鄰居接線
    // ==========================================
    
    std::unordered_set<std::shared_ptr<Block>> neighbors;
    for (auto& nb : newBlocks) neighbors.insert(nb); // 解開 Self-loop

    for (auto& seqPair : parent->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            if (auto p = segPair.second.getPrevBlock().lock()) neighbors.insert(p);
            if (auto n = segPair.second.getNextBlock().lock()) neighbors.insert(n);
        }
    }
    neighbors.erase(parent);

    std::vector<std::shared_ptr<Block>> neighbor_vec(neighbors.begin(), neighbors.end());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, neighbor_vec.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto currentBlock = neighbor_vec[i];

                for (auto& seqPair : currentBlock->getSequences()) {
                    std::string seqID = seqPair.first;
                    for (auto& segPair : seqPair.second.getSegments()) {
                        Segment& seg = segPair.second;

                        // 檢查 Prev
                        if (seg.getPrevBlock().lock() == parent) {
                            for (auto& nb : newBlocks) {
                                if (nb->getSequences().count(seqID)) {
                                    bool found = false;
                                    for (auto& nSeg : nb->getSequences().at(seqID).getSegments()) {
                                        if (nSeg.second.getEnd() == seg.getStart()) {
                                            seg.setPrevBlock(nb); found = true; break;
                                        }
                                    }
                                    if (found) break; // 提早 Break 節省時間
                                }
                            }
                        }
                        
                        // 檢查 Next
                        if (seg.getNextBlock().lock() == parent) {
                            for (auto& nb : newBlocks) {
                                if (nb->getSequences().count(seqID)) {
                                    bool found = false;
                                    for (auto& nSeg : nb->getSequences().at(seqID).getSegments()) {
                                        if (nSeg.second.getStart() == seg.getEnd()) {
                                            seg.setNextBlock(nb); found = true; break;
                                        }
                                    }
                                    if (found) break; // 提早 Break 節省時間
                                }
                            }
                        }
                    }
                }
            }
        }
    );
    

    // 5. 安全刪除舊 Block
    this->deleteBlock(parent->getId());

    auto time5 = std::chrono::high_resolution_clock::now();

    std::vector<Block::ID> resultIDs;
    for (auto& nb : newBlocks) resultIDs.push_back(nb->getId());

    auto time6 = std::chrono::high_resolution_clock::now();

    // std::cout << "Time 0: " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << "ms\n";
    // std::cout << "Time 1: " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "ms\n";
    // std::cout << "Time 2: " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << "ms\n";
    // std::cout << "Time 3: " << std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count() << "ms\n";
    // std::cout << "Time 4: " << std::chrono::duration_cast<std::chrono::milliseconds>(time5 - time4).count() << "ms\n";
    // std::cout << "Time 5: " << std::chrono::duration_cast<std::chrono::milliseconds>(time6 - time5).count() << "ms\n";
    

    return resultIDs;
}

void BlockSet::updateSegmentLinks(std::shared_ptr<Block> oldBlk, std::shared_ptr<Block> newBlk) {
    if (!oldBlk || !newBlk) return;

    bool debug = false;
    int updateCount = 0;
    
    // 這裡維持全域掃描是安全的，因為 oldBlk 是被「完全取代 (Merged)」，不牽涉到切分左右的座標問題
    for (auto& blockPair : blocks_) {
        auto currentBlock = blockPair.second;
        
        for (auto& seqPair : currentBlock->getSequences()) {
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;

                if (auto pBlk = seg.getPrevBlock().lock()) {
                    if (pBlk == oldBlk) {
                        seg.setPrevBlock(newBlk);
                        updateCount++;
                    }
                }

                if (auto nBlk = seg.getNextBlock().lock()) {
                    if (nBlk == oldBlk) {
                        seg.setNextBlock(newBlk);
                        updateCount++;
                    }
                }
            }
        }
    }
    
    if (updateCount > 0) {
        if (debug) std::cout << "[DEBUG-LINK] Merged Block " << oldBlk->getId() << " into " << newBlk->getId() 
                             << " | Re-wired " << updateCount << " external segment pointers.\n";
    }
}

std::map<int, Block::ID> BlockSet::splitBlocksByCuts(const std::set<int>& cuts) {
    std::map<int, Block::ID> blocksMap;
    
    // 取得最初始的 Consensus Blocks
    std::vector<Block::ID> involvedBlocks = this->getRepresentativeBlocks();
    if (involvedBlocks.empty()) return blocksMap;
    std::cerr << "Representative Blocks: " << involvedBlocks.size() << "\n";


    int currentStart = 0;
    std::shared_ptr<Block> currentBlk = this->getBlock(involvedBlocks[0]); 
    
    // 輔助函數：沿著 longest_sequence_ 找出這個 Block 的下一個 Block
    auto getNextBlockViaLongestSeq = [&](std::shared_ptr<Block> blk) -> std::shared_ptr<Block> {
        if (!blk) return nullptr;
        auto& seqs = blk->getSequences();
        auto seqIt = seqs.find(longest_sequence_);
        if (seqIt == seqs.end()) return nullptr;

        // 因為這條 sequence 可能在這個 block 被切成好幾個 segment
        // 我們要找這個 block 裡面 "最後一個" segment 的下一個指標
        int maxEnd = -1;
        std::weak_ptr<Block> nextWeak;
        for (auto& segPair : seqIt->second.getSegments()) {
            Segment& seg = segPair.second;
            int segEnd = std::max(seg.getStart(), seg.getEnd());
            if (segEnd > maxEnd) {
                maxEnd = segEnd;
                nextWeak = seg.getNextBlock();
            }
        }
        return nextWeak.lock();
    };

    auto cutsIt = cuts.begin();

    // 只要還有 Block 還沒走完，就繼續走
    while (currentBlk) {
        
        int blkLen = currentBlk->getConsensus().length();
        int blkEnd = currentStart + blkLen;

        // 檢查是否有切點落在這個 Block 內部 (不包含邊界)
        bool cutOccurred = false;
        while (cutsIt != cuts.end()) {
            int cut = *cutsIt;
            
            if (cut <= currentStart) {
                // 已經過去的切點，跳過
                cutsIt++;
                continue;
            }
            
            if (cut < blkEnd) {
                // 切點精準落在當前 Block 內部！執行切割
                int localCut = cut - currentStart;
                
                auto parts = this->splitSingleBlock(currentBlk->getId(), localCut);
                
                // 記錄前半段
                blocksMap[currentStart] = parts.first;
                
                // 將 currentBlk 替換成切出來的後半段，等待下一個迴圈檢查
                // (注意：此時 currentStart 會前進到 cut 的位置)
                currentBlk = this->getBlock(parts.second);
                currentStart = cut;
                cutOccurred = true;
                break; // 跳出內層迴圈，讓外層 while 重新評估這個新的 currentBlk
            }
            
            // 如果 cut >= blkEnd，代表這把刀在後面的 Block，保留 cutsIt 等待前進
            break; 
        }

        // 如果這個 Block 內部沒有發生任何切割，我們就把他完整記錄下來，並前進到下一個
        if (!cutOccurred) {
            blocksMap[currentStart] = currentBlk->getId();
            currentStart += currentBlk->getConsensus().length();
            
            // 沿著主幹前進
            currentBlk = getNextBlockViaLongestSeq(currentBlk);
        }
    }
    
    return blocksMap;
}
*/
/*
// Newer and Faster Version (functionality hasn'y been tested)
std::shared_ptr<Block> BlockSet::mergeTwoBlocks(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, const mga::Cigar& cigar, bool inverse) 
{
    bool debug = false;
    // 1. 取得舊 Consensus 與 SequenceInfo Map
    std::string refSeq = refBlock->getConsensus();
    std::string qrySeq = qryBlock->getConsensus();

    if (debug) std::cout << "\n[DEBUG-PRE-VALIDATION] Validating \n";
              
    for (auto& seqPair : refBlock->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariants()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    if (debug) {
                        gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                    }
                }
            }
            
            if (debug) {
                int calculatedConsensusLen = origLen + totalGapLen;
                std::cerr << "  [REF] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]\n";
                if (!gapDetails.empty()) {
                    std::cout << "      -> Contains " << gapDetails.size() << " GAPs: ";
                    for (const auto& detail : gapDetails) { std::cout << detail << " "; }
                    std::cout << "\n";
                } else {
                    std::cout << "      -> No GAPs.\n";
                }
            }
        }
    }
    
    for (auto& seqPair : qryBlock->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariants()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    if (debug) {
                        gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                    }
                }
            }
            
            if (debug) {
                int calculatedConsensusLen = origLen + totalGapLen;
                std::cerr << "  [QRY] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]\n";
                if (!gapDetails.empty()) {
                    std::cout << "      -> Contains " << gapDetails.size() << " GAPs: ";
                    for (const auto& detail : gapDetails) { std::cout << detail << " "; }
                    std::cout << "\n";
                } else {
                    std::cout << "      -> No GAPs.\n";
                }
            }
        }
    }
    if (debug) std::cout << "------------------------------------------------------------\n";
    
    auto refSeqs = refBlock->getSequences();
    auto qrySeqs = qryBlock->getSequences();

    // 2. 處理反股 (Reverse Complement) - 只反轉 Consensus 字串
    if (inverse) {
        auto rcString = [](const std::string& s) {
            std::string rc = s;
            std::reverse(rc.begin(), rc.end());
            for (char& c : rc) {
                switch (c) {
                    case 'A': c = 'T'; break; case 'T': c = 'A'; break;
                    case 'C': c = 'G'; break; case 'G': c = 'C'; break;
                    case 'a': c = 't'; break; case 't': c = 'a'; break;
                    case 'c': c = 'g'; break; case 'g': c = 'c'; break;
                }
            }
            return rc;
        };
        qrySeq = rcString(qrySeq);
    }

    // 3. CIGAR 長度嚴格校驗
    int cRefLen = 0, cQryLen = 0;
    for (const auto& op : cigar) {
        if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') cRefLen += op.first;
        if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I') cQryLen += op.first;
    }

    if (cRefLen != refSeq.length() || cQryLen != qrySeq.length()) {
        std::cerr << "\n[CRITICAL ERROR] CIGAR length mismatch!\n"
                  << "  Ref Block Len: " << refSeq.length() << " vs CIGAR Ref: " << cRefLen << "\n"
                  << "  Qry Block Len: " << qrySeq.length() << " vs CIGAR Qry: " << cQryLen << "\n";
        return refBlock; 
    }

    // ==========================================
    // [TBB 前置準備]: 將 Map 中的 Segment 取出成為連續的 Pointer Array
    // ==========================================
    std::vector<Segment*> refSegsFlat;
    for (auto& seqPair : refSeqs) {
        for (auto& segPairInner : seqPair.second.getSegments()) {
            refSegsFlat.push_back(&segPairInner.second);
        }
    }
    std::vector<Segment*> qrySegsFlat;
    for (auto& seqPair : qrySeqs) {
        for (auto& segPairInner : seqPair.second.getSegments()) {
            qrySegsFlat.push_back(&segPairInner.second);
        }
    }
    size_t numRefSegs = refSegsFlat.size();
    size_t totalSegs = numRefSegs + qrySegsFlat.size();
    int refConsLen = refSeq.length();
    int qryConsLen = qryBlock->getConsensus().length(); // 確保拿到原始長度

    // 4. 走訪 CIGAR，建構新的 Consensus 與座標對應表
    std::string mergedConsensus = "";
    mergedConsensus.reserve(refSeq.length() + qrySeq.length()); 

    std::vector<int> refOldToNew(refSeq.length() + 1, 0);
    std::vector<int> qryOldToNew(qrySeq.length() + 1, 0);

    std::vector<Variation> newRefGaps;
    std::vector<Variation> newQryGaps;

    // 【新增追蹤】：記錄舊共識跟新共識發生差異的座標
    std::vector<std::pair<int, char>> refConsensusChanges;
    std::vector<std::pair<int, char>> qryConsensusChanges;

    int rPos = 0, qPos = 0, mPos = 0; 

    auto getBaseFromSeg = [](Segment& seg, int pos, char defaultBase, bool needRc, int consLen) -> char {
        int lookupPos = pos;
        if (needRc) lookupPos = consLen - 1 - pos; 
        auto& vars = seg.getVariants();
        auto it = std::lower_bound(vars.begin(), vars.end(), lookupPos, 
            [](Variation& v, int p) { return v.getStart() < p; });
        
        if (it != vars.end() && it->getStart() == lookupPos && it->getType() == Variation::SNV) {
            char alt = it->getAlt();
            if (needRc) {
                switch (alt) {
                    case 'A': return 'T'; case 'T': return 'A';
                    case 'C': return 'G'; case 'G': return 'C';
                    case 'a': return 't'; case 't': return 'a';
                    case 'c': return 'g'; case 'g': return 'c';
                }
            }
            return alt;
        }
        return defaultBase;
    };

    // 用 std::array 讓 TBB Reduce 可以配發在 Stack 上，極大化效能
    using FreqArray = std::array<int, 256>;

    for (const auto& op : cigar) {
        int len = op.first;
        char type = op.second;

        if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                char rBase = refSeq[rPos];
                char qBase = qrySeq[qPos];
            
                if (rBase == qBase) {
                    mergedConsensus += rBase;
                } else {
                    FreqArray finalCounts = tbb::parallel_reduce(
                        tbb::blocked_range<size_t>(0, totalSegs),
                        FreqArray{}, // 初始化為全 0
                        [&](const tbb::blocked_range<size_t>& r, FreqArray localCounts) -> FreqArray {
                            for (size_t idx = r.begin(); idx != r.end(); ++idx) {
                                if (idx < numRefSegs) {
                                    char b = getBaseFromSeg(*(refSegsFlat[idx]), rPos, rBase, false, refConsLen);
                                    localCounts[(unsigned char)b]++;
                                } else {
                                    char b = getBaseFromSeg(*(qrySegsFlat[idx - numRefSegs]), qPos, qBase, inverse, qryConsLen);
                                    localCounts[(unsigned char)b]++;
                                }
                            }
                            return localCounts;
                        },
                        [](FreqArray a, const FreqArray& b) -> FreqArray {
                            for(int k=0; k<256; ++k) a[k] += b[k];
                            return a;
                        }
                    );
                
                    char bestBase = rBase; 
                    int maxFreq = -1; 
                    for (int k = 0; k < 256; ++k) {
                        if (finalCounts[k] > maxFreq) {
                            maxFreq = finalCounts[k];
                            bestBase = (char)k;
                        }
                    }
                    mergedConsensus += bestBase;
                }

                // 【核心修復 1】：如果舊的共識鹼基輸給了投票，記錄下來！
                if (rBase != mergedConsensus.back()) refConsensusChanges.push_back({rPos, rBase});
                if (qBase != mergedConsensus.back()) qryConsensusChanges.push_back({qPos, qBase});
                
                refOldToNew[rPos] = mPos;
                qryOldToNew[qPos] = mPos;
                rPos++; qPos++; mPos++;
            }
        }
        else if (type == 'I') { 
            mergedConsensus += qrySeq.substr(qPos, len);
            for(int i = 0; i < len; ++i) qryOldToNew[qPos + i] = mPos + i; 
            newRefGaps.push_back(Variation::createGap(mPos, mPos + len));
            qPos += len; mPos += len;
        } 
        else if (type == 'D') { 
            mergedConsensus += refSeq.substr(rPos, len);
            for(int i = 0; i < len; ++i) refOldToNew[rPos + i] = mPos + i; 
            newQryGaps.push_back(Variation::createGap(mPos, mPos + len));
            rPos += len; mPos += len;
        }
    }
    refOldToNew[rPos] = mPos;
    qryOldToNew[qPos] = mPos;

    // 5. 建立新的 Merged Block
    auto mergedBlock = this->createBlock(mergedConsensus);
    if (debug) std::cout << "[DEBUG] Mer: " << refBlock->getId() << " & " << qryBlock->getId() << " -> " << mergedBlock->getId() << '\n';

    // ==========================================
    // [TBB 加速 2] 平行處理所有 Segment 的變異轉換
    // 建立任務列，避免在平行迴圈中對 Map 執行 Insert 操作
    // ==========================================
    struct SegUpdateTask {
        Segment* seg;
        bool isQrySide;
        int originalConsLen;
        const std::vector<int>* oldToNew;
        const std::vector<Variation>* inducedGaps;
        const std::vector<std::pair<int, char>>* consensusChanges;
    };

    std::vector<SegUpdateTask> updateTasks;
    updateTasks.reserve(totalSegs);
    for (auto* seg : refSegsFlat) updateTasks.push_back({seg, false, refConsLen, &refOldToNew, &newRefGaps, &refConsensusChanges});
    for (auto* seg : qrySegsFlat) updateTasks.push_back({seg, true, qryConsLen, &qryOldToNew, &newQryGaps, &qryConsensusChanges});

    // 平行執行每個 Segment 內部的耗時操作 (如建立、排序 Variation、過濾疊加 Gaps)
    tbb::parallel_for(tbb::blocked_range<size_t>(0, updateTasks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto& task = updateTasks[i];
                Segment& seg = *(task.seg); // 直接使用指標就地修改，完全 Thread-Safe
                
                if (task.isQrySide && inverse) seg.reverseComplement(task.originalConsLen);
                
                std::vector<Variation> segmentGaps;
                std::vector<Variation> candidateSnvs;
                
                for (auto& var : seg.getVariants()) {
                    if (var.getType() == Variation::GAP) {
                        if (var.getStart() >= task.oldToNew->size() || var.getEnd() >= task.oldToNew->size()) continue;
                        int newStart = (*task.oldToNew)[var.getStart()];
                        int newEnd = (*task.oldToNew)[var.getEnd()];
                        segmentGaps.push_back(Variation::createGap(newStart, newEnd));
                    } else {
                        if (var.getStart() >= task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[var.getStart()];
                        if (var.getAlt() != mergedConsensus[newPos]) {
                            candidateSnvs.push_back(Variation(newPos, var.getAlt()));
                        }
                    }
                }
                
                segmentGaps.insert(segmentGaps.end(), task.inducedGaps->begin(), task.inducedGaps->end());

                std::sort(segmentGaps.begin(), segmentGaps.end(), [](Variation& a, Variation& b) {
                    return a.getStart() < b.getStart();
                });

                std::vector<Variation> mergedGaps;
                for (auto& gap : segmentGaps) {
                    if (mergedGaps.empty()) {
                        mergedGaps.push_back(gap);
                    } else {
                        auto& lastGap = mergedGaps.back();
                        if (lastGap.getEnd() >= gap.getStart()) { 
                            int mStart = lastGap.getStart();
                            int mEnd = std::max(lastGap.getEnd(), gap.getEnd());
                            mergedGaps.pop_back();
                            mergedGaps.push_back(Variation::createGap(mStart, mEnd));
                        } else {
                            mergedGaps.push_back(gap);
                        }
                    }
                }

                for (auto& change : *(task.consensusChanges)) {
                    int oldPos = change.first;
                    char oldConsBase = change.second;
                    
                    bool hasOldSnv = false;
                    for (auto& v : seg.getVariants()) {
                        if (v.getType() == Variation::SNV && v.getStart() == oldPos) {
                            hasOldSnv = true; break;
                        }
                    }
                    
                    if (!hasOldSnv) {
                        if (oldPos >= task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[oldPos];
                        if (oldConsBase != mergedConsensus[newPos]) { 
                            candidateSnvs.push_back(Variation(newPos, oldConsBase));
                        }
                    }
                }
                
                std::vector<Variation> finalSnvs;
                for (auto& snv : candidateSnvs) {
                    bool inGap = false;
                    for (auto& gap : mergedGaps) {
                        if (snv.getStart() >= gap.getStart() && snv.getStart() < gap.getEnd()) {
                            inGap = true; break;
                        }
                    }
                    if (!inGap) finalSnvs.push_back(snv);
                }
                
                std::vector<Variation> finalVars;
                finalVars.reserve(mergedGaps.size() + finalSnvs.size());
                finalVars.insert(finalVars.end(), mergedGaps.begin(), mergedGaps.end());
                finalVars.insert(finalVars.end(), finalSnvs.begin(), finalSnvs.end());

                std::sort(finalVars.begin(), finalVars.end(), [](Variation& a, Variation& b) {
                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                    return a.getType() > b.getType(); 
                });

                std::vector<Variation> cleanedVars;
                for (auto& var : finalVars) {
                    if (cleanedVars.empty()) {
                        cleanedVars.push_back(var);
                    } else {
                        auto& last = cleanedVars.back();
                        if (last.getStart() == var.getStart() && last.getType() == Variation::SNV && var.getType() == Variation::SNV) continue;
                        cleanedVars.push_back(var);
                    }
                }

                // 就地覆寫這個 Segment 專屬的 Variation
                seg.getVariants() = std::move(cleanedVars);
            }
        }
    );

    // ==========================================
    // 6. 循序且安全地將更新完畢的資料掛載回 Merged Block
    // (將 CPU 繁重的任務抽離給 TBB，掛載則使用安全的循序寫入)
    // ==========================================
    auto& mergedSeqs = mergedBlock->getSequences();
    
    for (auto& seqPair : refSeqs) {
        if (mergedSeqs.find(seqPair.first) == mergedSeqs.end()) {
            mergedSeqs[seqPair.first] = SequenceInfo(seqPair.first);
        }
        for (auto& segPairInner : seqPair.second.getSegments()) {
            mergedSeqs[seqPair.first].getSegments()[segPairInner.second.getStart()] = std::move(segPairInner.second);
        }
    }
    
    for (auto& seqPair : qrySeqs) {
        if (mergedSeqs.find(seqPair.first) == mergedSeqs.end()) {
            mergedSeqs[seqPair.first] = SequenceInfo(seqPair.first);
        }
        for (auto& segPairInner : seqPair.second.getSegments()) {
            mergedSeqs[seqPair.first].getSegments()[segPairInner.second.getStart()] = std::move(segPairInner.second);
        }
    }
    
    // 7. Validation 驗證合併後的 Segment 長度
    int expectedConsensusLen = mergedBlock->getConsensus().length();
    if (debug) std::cout << "\n[DEBUG-VALIDATION] Validating Merged Block ID: " << mergedBlock->getId() 
                         << " | Expected Consensus Length: " << expectedConsensusLen << "\n";
              
    for (auto& seqPair : mergedBlock->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariants()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    if (debug) {
                        gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                    }
                }
            }
            
            int calculatedConsensusLen = origLen + totalGapLen;
            if (calculatedConsensusLen != expectedConsensusLen) {
                std::cerr << "  [WARNING] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]"
                          << " | OrigLen: " << origLen << " + Gaps: " << totalGapLen 
                          << " = " << calculatedConsensusLen 
                          << " (Mismatch with " << expectedConsensusLen << ")\n";
            }
        }
    }
    
    // ==========================================
    // 8. 專屬 Debug 訊息：印出 < 500bp 區塊的合併細節
    // ==========================================
    // if (debug && expectedConsensusLen < 500) {
    if (debug) {
        size_t max_char = 500;
        std::cout << "\n============================================================\n"
                  << "=== [DEBUG-MERGE-DETAILS] Short Block Merge Info (< 500bp) ===\n"
                  << "============================================================\n";
        
        std::cout << "[Input Blocks]\n"
                  << "  Ref Block ID : " << refBlock->getId() << " (Len: " << refBlock->getConsensus().length() << ")\n"
                  << "  Ref Consensus: " << refSeq.substr(0, std::min(refBlock->getConsensus().length(), max_char)) << "\n"
                  << "  Qry Block ID : " << qryBlock->getId() << " (Len: " << qryBlock->getConsensus().length() << ")\n"
                  << "  Qry Consensus: " << qrySeq.substr(0, std::min(qryBlock->getConsensus().length(), max_char)) << (inverse ? " (Reversed)" : "") << "\n";
                  
        std::cout << "\n[Alignment CIGAR]\n  ";
        for (const auto& op : cigar) {
            std::cout << op.first << op.second;
        }
        std::cout << "\n";

        std::cout << "\n[Merged Result]\n"
                  << "  Merged Block ID : " << mergedBlock->getId() << " (Len: " << expectedConsensusLen << ")\n"
                  << "  Merged Consensus: " << mergedConsensus.substr(0, std::min(expectedConsensusLen, (int)max_char)) << "\n";

        std::cout << "\n[Segment & Variation Tracking]\n";
        for (auto& seqPair : mergedBlock->getSequences()) {
            std::cout << "  Sequence: " << seqPair.first << "\n";
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;
                std::cout << "    ├─ Seg [" << seg.getStart() << " -> " << seg.getEnd() << "]" 
                          << (seg.isReverse() ? " (-)" : " (+)") << "\n";
                
                auto& vars = seg.getVariants();
                if (vars.empty()) {
                    std::cout << "    │    └─ (No Variations)\n";
                } else {
                    for (size_t i = 0; i < vars.size(); ++i) {
                        auto& v = vars[i];
                        std::string branch = (i == vars.size() - 1) ? "    │    └─ " : "    │    ├─ ";
                        if (v.getType() == Variation::GAP) {
                            std::cout << branch << "GAP [" << v.getStart() << " -> " << v.getEnd() << "] (Len: " << (v.getEnd() - v.getStart()) << ")\n";
                        } else {
                            std::cout << branch << "SNV at " << v.getStart() << " (Alt: " << v.getAlt() << ")\n";
                        }
                    }
                }
            }
        }
        std::cout << "============================================================\n";
    }

    return mergedBlock;
}

void BlockSet::rebuildDictionary(std::map<int, SegNode>& dict, const std::string& targetSeqName) {
    dict.clear(); // 清空舊字典
    
    // 走訪目前所有的 Block
    for (const auto& blockPair : blocks_) {
        Block::ID blkId = blockPair.first;
        std::shared_ptr<Block> blk = blockPair.second;
        
        // 找到我們要的這條 Sequence (例如 Self-Mapping 的那條)
        auto seqIt = blk->getSequences().find(targetSeqName);
        if (seqIt != blk->getSequences().end()) {
            
            // 將這個 Block 裡面的所有 Segment 加進字典
            for (auto& segPair : seqIt->second.getSegments()) {
                Segment& seg = segPair.second;
                
                // 字典的 Key 必須是小到大，確保正反股都能正確放入字典
                int s = std::min(seg.getStart(), seg.getEnd());
                int e = std::max(seg.getStart(), seg.getEnd());
                
                dict[s] = {s, e, blkId};
            }
        }
    }
}


void BlockSet::refineBlocks() {
    // 把所有的 Block 取出來轉成 Vector，方便 TBB 切分任務
    std::vector<std::shared_ptr<Block>> all_blocks;
    all_blocks.reserve(blocks_.size());
    for (auto& kv : blocks_) {
        all_blocks.push_back(kv.second);
    }

    // 使用 TBB 平行處理所有 Block
    tbb::parallel_for(tbb::blocked_range<size_t>(0, all_blocks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                // 每個 Block 獨立清理自己的空殼 Column 跟縫合 Segment
                all_blocks[i]->refine();
            }
        }
    );

    // [選擇性]：如果 refine 之後有 Block 裡面的 Sequence 全空了，可以在這裡把它刪除
    std::vector<Block::ID> emptyBlocks;
    for (auto blk : all_blocks) {
        if (blk->getSequences().empty() || blk->getConsensus().empty()) {
            emptyBlocks.push_back(blk->getId());
        }
    }
    for (auto id : emptyBlocks) {
        this->deleteBlock(id);
    }
}
*/

/*
void BlockSet::refineFast() {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet RefineFast: SuperBlock Splitting + Fast Absorb ===\n"
                              << "============================================================\n";

    auto time0 = std::chrono::high_resolution_clock::now();
    // ==========================================
    // 內部神器 1：基於真實基因體座標的全局指標重建
    // ==========================================
    auto rebuildAllPointers = [&]() {
        auto allBlocks = this->getAllBlocks();

        tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto blk = allBlocks[i];
                    blk->clearLinkages();
                    for (auto& seqPair : blk->getSequences()) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            segPair.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                            segPair.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                        }
                    }
                }
            }
        );
        
        struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
        std::unordered_map<std::string, std::vector<SegRef>> seqTracks;
        for (auto blk : allBlocks) {
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    seqTracks[seqPair.first].push_back({ &segPairInner.second, blk });
                }
            }
        }

        std::vector<std::vector<SegRef>*> trackPtrs;
        trackPtrs.reserve(seqTracks.size());
        for (auto& trackPair : seqTracks) {
            trackPtrs.push_back(&trackPair.second);
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, trackPtrs.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto& track = *(trackPtrs[i]);
                    std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                        return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
                    });
                    for (size_t j = 0; j < track.size(); ++j) {
                        if (j > 0) {
                            auto& prevRef = track[j-1];
                            auto& currRef = track[j];
                            if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                            else prevRef.seg->setPrevBlock(currRef.blk); 
                            if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                            else currRef.seg->setNextBlock(prevRef.blk); 
                        }
                    }
                }
            }
        );

        for (auto blk : allBlocks) {
            std::set<std::shared_ptr<Block>> prevBlocksSet;
            std::set<std::shared_ptr<Block>> nextBlocksSet;
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    if (auto p = segPair.second.getPrevBlock().lock()) prevBlocksSet.insert(p);
                    if (auto n = segPair.second.getNextBlock().lock()) nextBlocksSet.insert(n);
                }
            }
            for (auto p : prevBlocksSet) blk->addPrevBlock(p);
            for (auto n : nextBlocksSet) blk->addNextBlock(n);
        }
    };

    // ==========================================
    // 內部神器 2：合併區塊 (為 Phase 4 準備)
    // ==========================================
    auto concatBlocks = [&](std::shared_ptr<Block> left, std::shared_ptr<Block> right) {
        int leftLen = left->getConsensus().length();
        int rightLen = right->getConsensus().length();
        auto newBlock = this->createBlock(left->getConsensus() + right->getConsensus());
        
        std::set<std::string> allSeqs;
        for (auto& kv : left->getSequences()) allSeqs.insert(kv.first);
        for (auto& kv : right->getSequences()) allSeqs.insert(kv.first);
        
        for (const auto& seqID : allSeqs) {
            SequenceInfo newSeqInfo(seqID);
            std::vector<Segment> leftSegs, rightSegs;
            if (left->getSequences().count(seqID)) {
                for (auto& kv : left->getSequences().at(seqID).getSegments()) leftSegs.push_back(kv.second);
            }
            if (right->getSequences().count(seqID)) {
                for (auto& kv : right->getSequences().at(seqID).getSegments()) rightSegs.push_back(kv.second);
            }
            
            std::vector<bool> rUsed(rightSegs.size(), false);
            
            for (auto& lSeg : leftSegs) {
                bool matched = false;
                for (size_t i = 0; i < rightSegs.size(); ++i) {
                    if (rUsed[i]) continue;
                    auto& rSeg = rightSegs[i];
                    
                    bool sameStrand = (lSeg.isReverse() == rSeg.isReverse());
                    bool isContiguousFwd = sameStrand && !lSeg.isReverse() && (lSeg.getEnd() == rSeg.getStart());
                    bool isContiguousRev = sameStrand && lSeg.isReverse() && (lSeg.getStart() == rSeg.getEnd());
                    
                    if (isContiguousFwd || isContiguousRev) {
                        Segment newSeg = lSeg;
                        newSeg.setStart(std::min(lSeg.getStart(), rSeg.getStart())); 
                        newSeg.setEnd(std::max(lSeg.getEnd(), rSeg.getEnd())); 

                        auto& vars = newSeg.getVariants();
                        for (auto v : rSeg.getVariants()) { 
                            v.shift(leftLen); 
                            if (!vars.empty() && vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && vars.back().getEnd() == v.getStart()) {
                                int oldStart = vars.back().getStart();
                                vars.pop_back();
                                vars.push_back(Variation::createGap(oldStart, v.getEnd()));
                            } else {
                                vars.push_back(v); 
                            }
                        }
                        
                        newSeg.setNextBlock(rSeg.getNextBlock().lock());
                        newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                        rUsed[i] = true;
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    Segment newSeg = lSeg;
                    auto& vars = newSeg.getVariants();
                    if (!vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == leftLen) {
                        int oldStart = vars.back().getStart();
                        vars.pop_back();
                        vars.push_back(Variation::createGap(oldStart, leftLen + rightLen));
                    } else {
                        vars.push_back(Variation::createGap(leftLen, leftLen + rightLen));
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            for (size_t i = 0; i < rightSegs.size(); ++i) {
                if (!rUsed[i]) {
                    Segment newSeg = rightSegs[i];
                    newSeg.getVariants().clear(); 
                    newSeg.getVariants().push_back(Variation::createGap(0, leftLen));
                    
                    auto& vars = newSeg.getVariants();
                    for (auto v : rightSegs[i].getVariants()) {
                        v.shift(leftLen);
                        if (!vars.empty() && vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && vars.back().getEnd() == v.getStart()) {
                            int oldStart = vars.back().getStart();
                            vars.pop_back();
                            vars.push_back(Variation::createGap(oldStart, v.getEnd()));
                        } else {
                            vars.push_back(v);
                        }
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            if (!newSeqInfo.getSegments().empty()) newBlock->addSequence(newSeqInfo);
        }
        
        std::set<std::shared_ptr<Block>> neighbors;
        auto addNeighbors = [&](std::shared_ptr<Block> b) {
            for (auto& seq: b->getSequences()) for (auto& seg: seq.second.getSegments()) if (auto sp = seg.second.getPrevBlock().lock()) neighbors.insert(sp);
            for (auto& seq: b->getSequences()) for (auto& seg: seq.second.getSegments()) if (auto sp = seg.second.getNextBlock().lock()) neighbors.insert(sp);
        };
        addNeighbors(left); addNeighbors(right);
        neighbors.erase(left); neighbors.erase(right);

        for (auto& neighbor : neighbors) {
            for (auto& seqPair : neighbor->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    Segment& seg = segPair.second;
                    if (seg.getPrevBlock().lock() == left || seg.getPrevBlock().lock() == right) seg.setPrevBlock(newBlock);
                    if (seg.getNextBlock().lock() == left || seg.getNextBlock().lock() == right) seg.setNextBlock(newBlock);
                }
            }
        }
        
        this->deleteBlock(left->getId());
        this->deleteBlock(right->getId());
        return newBlock;
    };

    // 初始化指標
    rebuildAllPointers();
    auto time1 = std::chrono::high_resolution_clock::now();
    // ==========================================
    // 階段 1: 找出長 Gap 並獨立 (Batching)
    // ==========================================
    bool splitOccurred = true;
    while (splitOccurred) {
        splitOccurred = false;
        auto currentBlocks = this->getAllBlocks(); 
        for (auto blk : currentBlocks) {
            if (!this->getBlock(blk->getId())) continue;

            int consLen = blk->getConsensus().length();
            int margin = 50; 

            std::set<int> boundary_set;
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    for (auto& v : segPair.second.getVariants()) {
                        if (v.getType() == Variation::GAP && (v.getEnd() - v.getStart() > 100)) {
                            if (v.getStart() > 0 && v.getStart() < consLen) boundary_set.insert(v.getStart());
                            if (v.getEnd() > 0 && v.getEnd() < consLen) boundary_set.insert(v.getEnd());
                        }
                    }
                }
            }

            if (boundary_set.empty()) continue;

            std::vector<int> bnds;
            bnds.push_back(0);
            for (int b : boundary_set) bnds.push_back(b);
            bnds.push_back(consLen);

            std::vector<int> valid_cuts;
            for (size_t i = 1; i < bnds.size() - 1; ++i) {
                int b = bnds[i];
                if ((b - bnds[i-1]) >= margin || (bnds[i+1] - b) >= margin) {
                    valid_cuts.push_back(b);
                }
            }

            if (!valid_cuts.empty()) {
                this->splitMultiBlocks(blk->getId(), valid_cuts);
                splitOccurred = true; 
            }
        }
    }

    auto time2 = std::chrono::high_resolution_clock::now();
    // ==========================================
    // 階段 2: 移除 Pure Gap Segments
    // ==========================================
    for (auto blk : this->getAllBlocks()) {
        std::vector<std::string> seqsToRemove;
        for (auto& seqPair : blk->getSequences()) {
            std::string seqID = seqPair.first;
            std::vector<int> segsToRemove;
            for (auto& segPairInner : seqPair.second.getSegments()) {
                if (segPairInner.second.getStart() == segPairInner.second.getEnd()) {
                    segsToRemove.push_back(segPairInner.first);
                }
            }
            for (int sCoord : segsToRemove) seqPair.second.getSegments().erase(sCoord);
            if (seqPair.second.getSegments().empty()) seqsToRemove.push_back(seqID);
        }
        for (const auto& s : seqsToRemove) blk->getSequences().erase(s);
    }
    rebuildAllPointers();
    auto time3 = std::chrono::high_resolution_clock::now();

    // ==========================================
    // 階段 3 & 4: Unzip (拆解嵌合體) 與 Local Merge
    // ==========================================
    bool topologyChanged = true;
    while (topologyChanged) {
        topologyChanged = false;

        bool unzippedOccurred = true;
        while (unzippedOccurred) {
            unzippedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            for (auto blk : currentBlocks) {
                if (!this->getBlock(blk->getId())) continue;
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

                struct Bucket { Block* pBlk; Block* nBlk; std::map<std::string, Segment> segs; };
                std::vector<Bucket> buckets;

                for (auto& seqPair : blk->getSequences()) {
                    std::string seqID = seqPair.first;
                    for (auto& segPairInner : seqPair.second.getSegments()) {
                        Segment& seg = segPairInner.second;
                        Block* pBlk = seg.getPrevBlock().lock().get();
                        Block* nBlk = seg.getNextBlock().lock().get();

                        bool placed = false;
                        for (auto& bucket : buckets) {
                            if (bucket.pBlk == pBlk && bucket.nBlk == nBlk && bucket.segs.count(seqID) == 0) {
                                bucket.segs[seqID] = seg; placed = true; break;
                            }
                        }
                        if (!placed) buckets.push_back({pBlk, nBlk, {{seqID, seg}}});
                    }
                }

                if (buckets.size() > 1) {
                    for (auto& bucket : buckets) {
                        auto newBlock = this->createBlock(blk->getConsensus());
                        for (auto& kv : bucket.segs) {
                            SequenceInfo newSeqInfo(kv.first);
                            newSeqInfo.getSegments()[kv.second.getStart()] = kv.second; 
                            newBlock->addSequence(newSeqInfo);
                            
                            // 局部指標重接
                            if (auto pBlk = kv.second.getPrevBlock().lock()) {
                                if (pBlk->getSequences().count(kv.first)) {
                                    for (auto& pSegPair : pBlk->getSequences().at(kv.first).getSegments()) {
                                        if (pSegPair.second.getNextBlock().lock() == blk) {
                                            pSegPair.second.setNextBlock(newBlock);
                                        } else if (pSegPair.second.getPrevBlock().lock() == blk) {
                                            pSegPair.second.setPrevBlock(newBlock);
                                        }
                                    }
                                }
                            }
                            if (auto nBlk = kv.second.getNextBlock().lock()) {
                                if (nBlk->getSequences().count(kv.first)) {
                                    for (auto& nSegPair : nBlk->getSequences().at(kv.first).getSegments()) {
                                        if (nSegPair.second.getPrevBlock().lock() == blk) {
                                            nSegPair.second.setPrevBlock(newBlock);
                                        } else if (nSegPair.second.getNextBlock().lock() == blk) {
                                            nSegPair.second.setNextBlock(newBlock);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    this->deleteBlock(blk->getId());
                    unzippedOccurred = true;
                    topologyChanged = true; 
                }
            }
        }

        bool mergedOccurred = true;
        while (mergedOccurred) {
            mergedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            std::unordered_set<Block::ID> processed_this_round; 

            for (auto blk : currentBlocks) {
                if (processed_this_round.count(blk->getId())) continue;
                if (!this->getBlock(blk->getId())) continue; 
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

                // --- 檢查 Prev 合併 ---
                std::shared_ptr<Block> sharedPrev = nullptr;
                bool allSamePrev = true;
                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPair : seqPair.second.getSegments()) {
                        auto pBlk = segPair.second.getPrevBlock().lock();
                        if (!pBlk) { allSamePrev = false; break; }
                        if (!sharedPrev) sharedPrev = pBlk;
                        else if (sharedPrev != pBlk) { allSamePrev = false; break; }
                    }
                    if (!allSamePrev) break;
                }

                bool sameCompositionPrev = true;
                if (allSamePrev && sharedPrev) {
                    if (sharedPrev->getSequences().size() != blk->getSequences().size()) sameCompositionPrev = false;
                    else {
                        for (auto& kv : blk->getSequences()) {
                            if (sharedPrev->getSequences().find(kv.first) == sharedPrev->getSequences().end()) { sameCompositionPrev = false; break; }
                        }
                    }
                }

                if (allSamePrev && sharedPrev && sharedPrev != blk && sameCompositionPrev && !processed_this_round.count(sharedPrev->getId())) {
                    auto newBlk = concatBlocks(sharedPrev, blk);
                    processed_this_round.insert(blk->getId());
                    processed_this_round.insert(sharedPrev->getId());
                    processed_this_round.insert(newBlk->getId());
                    mergedOccurred = true;
                    topologyChanged = true;
                    continue; 
                }

                // --- 檢查 Next 合併 ---
                std::shared_ptr<Block> sharedNext = nullptr;
                bool allSameNext = true;
                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPair : seqPair.second.getSegments()) {
                        auto nBlk = segPair.second.getNextBlock().lock();
                        if (!nBlk) { allSameNext = false; break; }
                        if (!sharedNext) sharedNext = nBlk;
                        else if (sharedNext != nBlk) { allSameNext = false; break; }
                    }
                    if (!allSameNext) break;
                }

                bool sameCompositionNext = true;
                if (allSameNext && sharedNext) {
                    if (sharedNext->getSequences().size() != blk->getSequences().size()) sameCompositionNext = false;
                    else {
                        for (auto& kv : blk->getSequences()) {
                            if (sharedNext->getSequences().find(kv.first) == sharedNext->getSequences().end()) { sameCompositionNext = false; break; }
                        }
                    }
                }

                if (allSameNext && sharedNext && sharedNext != blk && sameCompositionNext && !processed_this_round.count(sharedNext->getId())) {
                    auto newBlk = concatBlocks(blk, sharedNext);
                    processed_this_round.insert(blk->getId());
                    processed_this_round.insert(sharedNext->getId());
                    processed_this_round.insert(newBlk->getId());
                    mergedOccurred = true;
                    topologyChanged = true;
                }
            }
        }
    } 
    
    if (topologyChanged) rebuildAllPointers();
    auto time4 = std::chrono::high_resolution_clock::now();
    // ==========================================
    // 階段 5: 全新 Micro-block 無痛吸收 (極速版)
    // ==========================================
    bool absorptionOccurred = true;
    while (absorptionOccurred) {
        absorptionOccurred = false;
        auto currentBlocks = this->getAllBlocks();

        std::vector<std::shared_ptr<Block>> microBlocks;
        for (auto blk : currentBlocks) {
            if (!this->getBlock(blk->getId())) continue;
            int len = blk->getConsensus().length();
            if (len > 0 && len < 50) microBlocks.push_back(blk);
        }

        // 這兩行非常耗時，現在我們讓它每次 while 迴圈只建立一次！
        std::unordered_map<std::string, std::unordered_map<int, std::pair<std::shared_ptr<Block>, Segment>>> hostStartMap;
        std::unordered_map<std::string, std::unordered_map<int, std::pair<std::shared_ptr<Block>, Segment>>> hostEndMap;

        for (auto b : currentBlocks) {
            if (!this->getBlock(b->getId())) continue;
            for (auto& seqPair : b->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    hostStartMap[seqPair.first][segPair.second.getStart()] = {b, segPair.second};
                    hostEndMap[seqPair.first][segPair.second.getEnd()] = {b, segPair.second};
                }
            }
        }

        for (auto mBlk : microBlocks) {
            if (!this->getBlock(mBlk->getId())) continue;
            if (mBlk->getSequences().empty()) continue;

            std::vector<std::string> seqIDs;
            for (auto& kv : mBlk->getSequences()) seqIDs.push_back(kv.first);

            for (const auto& mSeqID : seqIDs) {
                if (mBlk->getSequences().count(mSeqID) == 0) continue;
                auto& mSegMap = mBlk->getSequences().at(mSeqID).getSegments();
                if (mSegMap.empty()) continue;

                Segment mSeg = mSegMap.begin()->second;
                std::shared_ptr<Block> hostToMatch = nullptr;
                Segment hSeg;
                bool hostIsBeforeMicro = false;

                if (hostEndMap[mSeqID].count(mSeg.getStart())) {
                    hostToMatch = hostEndMap[mSeqID][mSeg.getStart()].first;
                    hSeg = hostEndMap[mSeqID][mSeg.getStart()].second;
                    hostIsBeforeMicro = true;
                } else if (hostStartMap[mSeqID].count(mSeg.getEnd())) {
                    hostToMatch = hostStartMap[mSeqID][mSeg.getEnd()].first;
                    hSeg = hostStartMap[mSeqID][mSeg.getEnd()].second;
                    hostIsBeforeMicro = false;
                }

                if (hostToMatch && this->getBlock(hostToMatch->getId()) && hostToMatch != mBlk) {

                    // 💡【核心加速】：在動刀前，先把牽涉到的舊座標從字典中拔除！
                    hostStartMap[mSeqID].erase(hSeg.getStart());
                    hostEndMap[mSeqID].erase(hSeg.getEnd());
                    hostStartMap[mSeqID].erase(mSeg.getStart());
                    hostEndMap[mSeqID].erase(mSeg.getEnd());

                    std::string mSeq = "";
                    std::string mCons = mBlk->getConsensus();
                    for (int i = 0; i < mCons.length(); ++i) {
                        bool isGap = false;
                        char c = mCons[i];
                        for (auto& v : mSeg.getVariants()) {
                            if (v.getType() == Variation::GAP && i >= v.getStart() && i < v.getEnd()) {
                                isGap = true; break;
                            } else if (v.getType() == Variation::SNV && i == v.getStart()) {
                                c = v.getAlt();
                            }
                        }
                        if (!isGap) mSeq += c;
                    }

                    if (mSeg.isReverse() != hSeg.isReverse()) {
                        std::reverse(mSeq.begin(), mSeq.end());
                        for (char& c : mSeq) {
                            if (c == 'A') c = 'T'; else if (c == 'T') c = 'A';
                            else if (c == 'C') c = 'G'; else if (c == 'G') c = 'C';
                        }
                    }

                    int mLen = mSeq.length();
                    if (mLen == 0) {
                        mBlk->getSequences().at(mSeqID).getSegments().erase(mSeg.getStart());
                        if (mBlk->getSequences().at(mSeqID).getSegments().empty()) mBlk->getSequences().erase(mSeqID);
                        continue;
                    }

                    bool isFront = (!hSeg.isReverse()) ? !hostIsBeforeMicro : hostIsBeforeMicro;
                    std::string hostCons = hostToMatch->getConsensus();
                    auto& realHSeg = hostToMatch->getSequences().at(mSeqID).getSegments().at(hSeg.getStart());
                    auto& vars = realHSeg.getVariants();

                    // === [原封不動保留你的 Gap Padding 與 Consensus 擴充邏輯] ===
                    if (isFront) {
                        int gap_end = 0;
                        if (!vars.empty() && vars.front().getType() == Variation::GAP && vars.front().getStart() == 0) gap_end = vars.front().getEnd();
                        int available_gap = gap_end;
                        int offset = 0;

                        if (mLen <= available_gap) {
                            offset = available_gap - mLen;
                            if (mLen == available_gap) vars.erase(vars.begin());
                            else vars.front().setEnd(available_gap - mLen);
                        } else {
                            int overflow = mLen - available_gap;
                            hostCons = mSeq.substr(0, overflow) + hostCons;
                            hostToMatch->setConsensus(hostCons);

                            for (auto& seqPair : hostToMatch->getSequences()) {
                                for (auto& segPair : seqPair.second.getSegments()) {
                                    Segment& s = segPair.second;
                                    auto& s_vars = s.getVariants();
                                    for (auto& v : s_vars) v.shift(overflow);
                                    if (&s != &realHSeg) {
                                        if (!s_vars.empty() && s_vars.front().getType() == Variation::GAP && s_vars.front().getStart() == overflow) {
                                            s_vars.front().setStart(0);
                                        } else s_vars.insert(s_vars.begin(), Variation::createGap(0, overflow));
                                    }
                                }
                            }
                            if (available_gap > 0 && !vars.empty() && vars.front().getType() == Variation::GAP && vars.front().getStart() == overflow) {
                                vars.erase(vars.begin());
                            }
                        }

                        for (int i = 0; i < mLen; ++i) {
                            int c_pos = offset + i;
                            if (mSeq[i] != hostCons[c_pos]) vars.push_back(Variation(c_pos, mSeq[i]));
                        }

                        if (!realHSeg.isReverse()) {
                            realHSeg.setStart(mSeg.getStart());
                            realHSeg.setPrevBlock(mSeg.getPrevBlock().lock());
                        } else {
                            realHSeg.setEnd(mSeg.getEnd());
                            realHSeg.setPrevBlock(mSeg.getPrevBlock().lock());
                        }
                    } else { 
                        int hLen = hostCons.length();
                        int gap_start = hLen;
                        if (!vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == hLen) gap_start = vars.back().getStart();
                        int available_gap = hLen - gap_start;
                        int offset = gap_start;

                        if (mLen <= available_gap) {
                            if (mLen == available_gap) vars.pop_back();
                            else vars.back().setStart(gap_start + mLen);
                        } else {
                            int overflow = mLen - available_gap;
                            hostCons += mSeq.substr(available_gap, overflow);
                            hostToMatch->setConsensus(hostCons);

                            for (auto& seqPair : hostToMatch->getSequences()) {
                                for (auto& segPair : seqPair.second.getSegments()) {
                                    Segment& s = segPair.second;
                                    auto& s_vars = s.getVariants();
                                    if (&s != &realHSeg) {
                                        if (!s_vars.empty() && s_vars.back().getType() == Variation::GAP && s_vars.back().getEnd() == hLen) {
                                            s_vars.back().setEnd(hLen + overflow);
                                        } else s_vars.push_back(Variation::createGap(hLen, hLen + overflow));
                                    }
                                }
                            }
                            if (available_gap > 0 && !vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == hLen) {
                                vars.pop_back();
                            }
                        }

                        for (int i = 0; i < mLen; ++i) {
                            int c_pos = offset + i;
                            if (mSeq[i] != hostCons[c_pos]) vars.push_back(Variation(c_pos, mSeq[i]));
                        }

                        if (!realHSeg.isReverse()) {
                            realHSeg.setEnd(mSeg.getEnd());
                            realHSeg.setNextBlock(mSeg.getNextBlock().lock());
                        } else {
                            realHSeg.setStart(mSeg.getStart());
                            realHSeg.setNextBlock(mSeg.getNextBlock().lock());
                        }
                    }

                    std::sort(vars.begin(), vars.end(), [](Variation& a, Variation& b){
                        if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                        return a.getType() > b.getType();
                    });
                    std::vector<Variation> mergedVars;
                    for (auto& v : vars) {
                        if (mergedVars.empty()) mergedVars.push_back(v);
                        else {
                            auto& last = mergedVars.back();
                            if (last.getType() == Variation::GAP && v.getType() == Variation::GAP && last.getEnd() >= v.getStart()) {
                                int mStart = last.getStart();
                                int mEnd = std::max(last.getEnd(), v.getEnd());
                                mergedVars.pop_back();
                                mergedVars.push_back(Variation::createGap(mStart, mEnd));
                            } else mergedVars.push_back(v);
                        }
                    }
                    vars = mergedVars;
                    // === [邏輯結束] ===

                    // 💡【核心加速】：把合併後的新宿主座標，重新寫回 Map 供後續 Micro Blocks 使用！
                    int newStart = realHSeg.getStart();
                    if (newStart != hSeg.getStart()) {
                        auto& hostSegsMap = hostToMatch->getSequences().at(mSeqID).getSegments();
                        Segment movedSeg = std::move(realHSeg);
                        hostSegsMap.erase(hSeg.getStart());
                        hostSegsMap[newStart] = std::move(movedSeg);
                        
                        // 將位移後的 Segment 註冊回 Map
                        hostStartMap[mSeqID][newStart] = {hostToMatch, hostSegsMap[newStart]};
                        hostEndMap[mSeqID][hostSegsMap[newStart].getEnd()] = {hostToMatch, hostSegsMap[newStart]};
                    } else {
                        // 座標沒變，直接把更新好的 Segment 寫回 Map
                        hostStartMap[mSeqID][realHSeg.getStart()] = {hostToMatch, realHSeg};
                        hostEndMap[mSeqID][realHSeg.getEnd()] = {hostToMatch, realHSeg};
                    }

                    if (auto pB = mSeg.getPrevBlock().lock()) {
                        if (pB->getSequences().count(mSeqID)) {
                            for (auto& pSegPair : pB->getSequences().at(mSeqID).getSegments()) {
                                if (pSegPair.second.getNextBlock().lock() == mBlk) pSegPair.second.setNextBlock(hostToMatch);
                                else if (pSegPair.second.getPrevBlock().lock() == mBlk) pSegPair.second.setPrevBlock(hostToMatch);
                            }
                        }
                    }
                    if (auto nB = mSeg.getNextBlock().lock()) {
                        if (nB->getSequences().count(mSeqID)) {
                            for (auto& nSegPair : nB->getSequences().at(mSeqID).getSegments()) {
                                if (nSegPair.second.getPrevBlock().lock() == mBlk) nSegPair.second.setPrevBlock(hostToMatch);
                                else if (nSegPair.second.getNextBlock().lock() == mBlk) nSegPair.second.setNextBlock(hostToMatch);
                            }
                        }
                    }

                    mBlk->getSequences().at(mSeqID).getSegments().erase(mSeg.getStart());
                    if (mBlk->getSequences().at(mSeqID).getSegments().empty()) {
                        mBlk->getSequences().erase(mSeqID);
                    }

                    absorptionOccurred = true;
                    // 💡【移除 Break】：不再提早中斷，繼續處理同一個 mBlk 裡的其他 Sequence！
                }
            } 

            if (mBlk->getSequences().empty()) {
                this->deleteBlock(mBlk->getId());
            }

            // 💡【移除 Break】：不再提早跳出微區塊迴圈，讓所有 microBlocks 在同一個 pass 被徹底掃蕩！
        }
    }
    auto time5 = std::chrono::high_resolution_clock::now();

    // ==========================================
    // 最終清理收尾
    // ==========================================
    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    rebuildAllPointers();


    std::cout << "Time 1  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << " ms\n";
    std::cout << "Time 2  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << " ms\n";
    std::cout << "Time 3  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << " ms\n";
    std::cout << "Time 4  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count() << " ms\n";
    std::cout << "Time 5  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time5 - time4).count() << " ms\n";
              


}




void BlockSet::refine() {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                               << "=== BlockSet Refine: Optimizing Graph Topology ===\n"
                               << "============================================================\n";


    const int MAX_GAP_LEN = 50; //Longest allowed gap in a segment
    // ==========================================
    // 內部神器：基於真實基因體座標的全局指標重建！
    // ==========================================
    auto rebuildAllPointers = [&]() {
        auto allBlocks = this->getAllBlocks();

        tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto blk = allBlocks[i];
                    blk->clearLinkages();
                    for (auto& seqPair : blk->getSequences()) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            segPair.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                            segPair.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                        }
                    }
                }
            }
        );
        struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
        std::unordered_map<std::string, std::vector<SegRef>> seqTracks;
        for (auto blk : allBlocks) {
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    seqTracks[seqPair.first].push_back({ &segPairInner.second, blk });
                }
            }
        }

        std::vector<std::vector<SegRef>*> trackPtrs;
        trackPtrs.reserve(seqTracks.size());
        for (auto& trackPair : seqTracks) {
            trackPtrs.push_back(&trackPair.second);
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, trackPtrs.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto& track = *(trackPtrs[i]);
                    std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                        return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
                    });
                    for (size_t j = 0; j < track.size(); ++j) {
                        if (j > 0) {
                            auto& prevRef = track[j-1];
                            auto& currRef = track[j];
                            if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                            else prevRef.seg->setPrevBlock(currRef.blk);
                            if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                            else currRef.seg->setNextBlock(prevRef.blk);
                        }
                    }
                }
            }
        );

        for (auto blk : allBlocks) {
            std::set<std::shared_ptr<Block>> prevBlocksSet;
            std::set<std::shared_ptr<Block>> nextBlocksSet;
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    if (auto p = segPair.second.getPrevBlock().lock()) prevBlocksSet.insert(p);
                    if (auto n = segPair.second.getNextBlock().lock()) nextBlocksSet.insert(n);
                }
            }
            for (auto p : prevBlocksSet) blk->addPrevBlock(p);
            for (auto n : nextBlocksSet) blk->addNextBlock(n);
        }
    };
    auto getRawSeq = [](const std::string& cons, Segment& seg) {
        std::string raw = "";
        int pos = 0;
        for (auto& var : seg.getVariants()) {
            if (var.getStart() > pos) raw += cons.substr(pos, var.getStart() - pos);
            if (var.getType() == Variation::SNV) {
                raw += var.getAlt();
                pos = var.getStart() + 1;
            } else if (var.getType() == Variation::GAP) {
                pos = var.getEnd();
            }
        }
        if (pos < cons.length()) raw += cons.substr(pos);
        return raw;
    };

    // ==========================================
    // 階段 0: 抽出寄生短 Segment
    // ==========================================
    for (auto blk : this->getAllBlocks()) {
        int consLen = blk->getConsensus().length();
        if (consLen < 200) continue;

        std::vector<std::pair<std::string, int>> segsToRemove;
        for (auto& seqPair : blk->getSequences()) {
            std::string seqID = seqPair.first;
            for (auto& segPairInner : seqPair.second.getSegments()) {
                Segment& seg = segPairInner.second;
                int actualLen = std::abs(seg.getEnd() - seg.getStart());
                if (actualLen > 0 && actualLen < 50 && actualLen < (consLen * 0.1)) {
                    segsToRemove.push_back({seqID, seg.getStart()});
                    std::string rawSeq = getRawSeq(blk->getConsensus(), seg);
                    auto isolatedBlk = this->createBlock(rawSeq);
                    Segment cleanSeg = seg;
                    cleanSeg.getVariants().clear();
                    SequenceInfo newSeqInfo(seqID);
                    newSeqInfo.getSegments()[cleanSeg.getStart()] = cleanSeg;
                    isolatedBlk->addSequence(newSeqInfo);
                }
            }
        }
        for (auto& rmPair : segsToRemove) {
            blk->getSequences()[rmPair.first].getSegments().erase(rmPair.second);
            if (blk->getSequences()[rmPair.first].getSegments().empty()) blk->getSequences().erase(rmPair.first);
        }
    }
    rebuildAllPointers();

    // ==========================================
    // 階段 1: 找出長 Gap 並獨立 (Batching)
    // ==========================================
    bool splitOccurred = true;
    while (splitOccurred) {
        splitOccurred = false;
        auto currentBlocks = this->getAllBlocks();
        for (auto blk : currentBlocks) {
            if (!this->getBlock(blk->getId())) continue; // 防呆，避免處理到已經被切掉的

            int cutPos = -1;
            int consLen = blk->getConsensus().length();
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    for (auto& v : segPair.second.getVariants()) {
                        if (v.getType() == Variation::GAP && (v.getEnd() - v.getStart() > MAX_GAP_LEN)) {
                            if (v.getStart() > 0 && v.getStart() < consLen) { cutPos = v.getStart(); break; }
                            if (v.getEnd() > 0 && v.getEnd() < consLen) { cutPos = v.getEnd(); break; }
                        }
                    }
                    if (cutPos != -1) break;
                }
                if (cutPos != -1) break;
            }
            if (cutPos != -1) {
                this->splitSingleBlock(blk->getId(), cutPos);
                splitOccurred = true;
                // 【優化】：取消 break，讓它在同一次掃描中把其他有長 Gap 的 Block 也切完
            }
        }
    }

    // ==========================================
    // 階段 2: 移除 Pure Gap Segments
    // ==========================================
    for (auto blk : this->getAllBlocks()) {
        std::vector<std::string> seqsToRemove;
        for (auto& seqPair : blk->getSequences()) {
            std::string seqID = seqPair.first;
            std::vector<int> segsToRemove;
            for (auto& segPairInner : seqPair.second.getSegments()) {
                if (segPairInner.second.getStart() == segPairInner.second.getEnd()) {
                    segsToRemove.push_back(segPairInner.first);
                }
            }
            for (int sCoord : segsToRemove) seqPair.second.getSegments().erase(sCoord);
            if (seqPair.second.getSegments().empty()) seqsToRemove.push_back(seqID);
        }
        for (const auto& s : seqsToRemove) blk->getSequences().erase(s);
    }
    rebuildAllPointers();

    // ==========================================
    // 太極迴圈核心工具：合併區塊 (Concat)
    // ==========================================
    auto concatBlocks = [&](std::shared_ptr<Block> left, std::shared_ptr<Block> right) {
        int leftLen = left->getConsensus().length();
        int rightLen = right->getConsensus().length();
        auto newBlock = this->createBlock(left->getConsensus() + right->getConsensus());
        std::set<std::string> allSeqs;
        for (auto& kv : left->getSequences()) allSeqs.insert(kv.first);
        for (auto& kv : right->getSequences()) allSeqs.insert(kv.first);
        for (const auto& seqID : allSeqs) {
            SequenceInfo newSeqInfo(seqID);
            std::vector<Segment> leftSegs, rightSegs;
            if (left->getSequences().count(seqID)) {
                for (auto& kv : left->getSequences().at(seqID).getSegments()) leftSegs.push_back(kv.second);
            }
            if (right->getSequences().count(seqID)) {
                for (auto& kv : right->getSequences().at(seqID).getSegments()) rightSegs.push_back(kv.second);
            }
            std::vector<bool> rUsed(rightSegs.size(), false);
            for (auto& lSeg : leftSegs) {
                bool matched = false;
                for (size_t i = 0; i < rightSegs.size(); ++i) {
                    if (rUsed[i]) continue;
                    auto& rSeg = rightSegs[i];
                    bool isContiguousFwd = (!lSeg.isReverse() && !rSeg.isReverse() && lSeg.getEnd() == rSeg.getStart());
                    bool isContiguousRev = (lSeg.isReverse() && rSeg.isReverse() && lSeg.getStart() == rSeg.getEnd());
                    if (isContiguousFwd || isContiguousRev) {
                        Segment newSeg = lSeg;
                        if (isContiguousFwd) newSeg.setEnd(rSeg.getEnd());
                        else newSeg.setStart(rSeg.getStart());

                        auto& vars = newSeg.getVariants();
                        for (auto v : rSeg.getVariants()) {
                            v.shift(leftLen);
                            if (!vars.empty() && vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && vars.back().getEnd() == v.getStart()) {
                                int oldStart = vars.back().getStart();
                                vars.pop_back();
                                vars.push_back(Variation::createGap(oldStart, v.getEnd()));
                            } else {
                                vars.push_back(v);
                            }
                        }
                        if (!lSeg.isReverse()) newSeg.setNextBlock(rSeg.getNextBlock().lock());
                        else newSeg.setPrevBlock(rSeg.getPrevBlock().lock());
                        newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                        rUsed[i] = true;
                        matched = true;
                        break;
                    }
                }
                if (!matched) {
                    Segment newSeg = lSeg;
                    auto& vars = newSeg.getVariants();
                    if (!vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == leftLen) {
                        int oldStart = vars.back().getStart();
                        vars.pop_back();
                        vars.push_back(Variation::createGap(oldStart, leftLen + rightLen));
                    } else {
                        vars.push_back(Variation::createGap(leftLen, leftLen + rightLen));
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            for (size_t i = 0; i < rightSegs.size(); ++i) {
                if (!rUsed[i]) {
                    Segment newSeg = rightSegs[i];
                    newSeg.getVariants().clear();
                    newSeg.getVariants().push_back(Variation::createGap(0, leftLen));
                    auto& vars = newSeg.getVariants();
                    for (auto v : rightSegs[i].getVariants()) {
                        v.shift(leftLen);
                        if (!vars.empty() && vars.back().getType() == Variation::GAP && v.getType() == Variation::GAP && vars.back().getEnd() == v.getStart()) {
                            int oldStart = vars.back().getStart();
                            vars.pop_back();
                            vars.push_back(Variation::createGap(oldStart, v.getEnd()));
                        } else {
                            vars.push_back(v);
                        }
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            if (!newSeqInfo.getSegments().empty()) newBlock->addSequence(newSeqInfo);
        }
        // ==========================================
        // 【終極效能修復】：局部鄰居更新 (Local Update)
        // 拋棄掃描全部 Graph，只找出真正有牽連的鄰居進行修改
        // ==========================================
        std::set<std::shared_ptr<Block>> neighbors;
        auto addNeighbors = [&](std::shared_ptr<Block> b) {
            for (auto& seq: b->getSequences()) for (auto& seg: seq.second.getSegments()) if (auto sp = seg.second.getPrevBlock().lock()) neighbors.insert(sp);
            for (auto& seq: b->getSequences()) for (auto& seg: seq.second.getSegments()) if (auto sp = seg.second.getNextBlock().lock()) neighbors.insert(sp);
        };
        addNeighbors(left);
        addNeighbors(right);
        neighbors.erase(left);
        neighbors.erase(right);

        for (auto& neighbor : neighbors) {
            for (auto& seqPair : neighbor->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    Segment& seg = segPair.second;
                    if (seg.getPrevBlock().lock() == left || seg.getPrevBlock().lock() == right) seg.setPrevBlock(newBlock);
                    if (seg.getNextBlock().lock() == left || seg.getPrevBlock().lock() == right) seg.setNextBlock(newBlock);
                }
            }
        }
        this->deleteBlock(left->getId());
        this->deleteBlock(right->getId());
        return newBlock;
    };

    // ==========================================
    // 太極迴圈：Phase 3 (Unzip) 與 Phase 4 (Concat)
    // ==========================================
    bool topologyChanged = true;
    while (topologyChanged) {
        topologyChanged = false;

        bool unzippedOccurred = true;
        while (unzippedOccurred) {
            unzippedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            for (auto blk : currentBlocks) {
                if (!this->getBlock(blk->getId())) continue;
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

                struct Bucket {
                    Block* pBlk; Block* nBlk; std::map<std::string, Segment> segs;
                };
                std::vector<Bucket> buckets;

                for (auto& seqPair : blk->getSequences()) {
                    std::string seqID = seqPair.first;
                    for (auto& segPairInner : seqPair.second.getSegments()) {
                        Segment& seg = segPairInner.second;
                        Block* pBlk = (!seg.isReverse()) ? seg.getPrevBlock().lock().get() : seg.getNextBlock().lock().get();
                        Block* nBlk = (!seg.isReverse()) ? seg.getNextBlock().lock().get() : seg.getPrevBlock().lock().get();

                        bool placed = false;
                        for (auto& bucket : buckets) {
                            if (bucket.pBlk == pBlk && bucket.nBlk == nBlk && bucket.segs.count(seqID) == 0) {
                                bucket.segs[seqID] = seg; placed = true; break;
                            }
                        }
                        if (!placed) buckets.push_back({pBlk, nBlk, {{seqID, seg}}});
                    }
                }

                if (buckets.size() > 1) {
                    for (auto& bucket : buckets) {
                        auto newBlock = this->createBlock(blk->getConsensus());
                        for (auto& kv : bucket.segs) {
                            SequenceInfo newSeqInfo(kv.first);
                            newSeqInfo.getSegments()[kv.second.getStart()] = kv.second;
                            newBlock->addSequence(newSeqInfo);
                        }
                    }
                    this->deleteBlock(blk->getId());
                    unzippedOccurred = true;
                    topologyChanged = true;
                    // 【優化】：取消 break，一次性 Unzip 到底
                }
            }
            // 只有在這個完整 Pass 中有發生 Unzip 時，才做一次指標重建
            if (unzippedOccurred) rebuildAllPointers();
        }

        bool mergedOccurred = true;
        while (mergedOccurred) {
            mergedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            std::unordered_set<Block::ID> processed_this_round; // 【優化】：防止在同一個 Pass 合併重複的 Block

            for (auto blk : currentBlocks) {
                if (processed_this_round.count(blk->getId())) continue;
                if (!this->getBlock(blk->getId())) continue;
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

                // --- 檢查 Prev 合併 ---
                std::shared_ptr<Block> sharedPrev = nullptr;
                bool allSamePrev = true;
                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPair : seqPair.second.getSegments()) {
                        auto pBlk = (!segPair.second.isReverse()) ? segPair.second.getPrevBlock().lock() : segPair.second.getNextBlock().lock();
                        if (!pBlk) { allSamePrev = false; break; }
                        if (!sharedPrev) sharedPrev = pBlk;
                        else if (sharedPrev != pBlk) { allSamePrev = false; break; }
                    }
                    if (!allSamePrev) break;
                }

                bool sameCompositionPrev = true;
                if (allSamePrev && sharedPrev) {
                    if (sharedPrev->getSequences().size() != blk->getSequences().size()) sameCompositionPrev = false;
                    else {
                        for (auto& kv : blk->getSequences()) {
                            if (sharedPrev->getSequences().find(kv.first) == sharedPrev->getSequences().end()) { sameCompositionPrev = false; break; }
                        }
                    }
                }

                if (allSamePrev && sharedPrev && sharedPrev != blk && sameCompositionPrev && !processed_this_round.count(sharedPrev->getId())) {
                    auto newBlk = concatBlocks(sharedPrev, blk);
                    processed_this_round.insert(blk->getId());
                    processed_this_round.insert(sharedPrev->getId());
                    processed_this_round.insert(newBlk->getId());
                    mergedOccurred = true;
                    topologyChanged = true;
                    continue; // 成功合併 Prev，跳過這個小區塊的 Next 檢查
                }

                // --- 檢查 Next 合併 ---
                std::shared_ptr<Block> sharedNext = nullptr;
                bool allSameNext = true;
                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPair : seqPair.second.getSegments()) {
                        auto nBlk = (!segPair.second.isReverse()) ? segPair.second.getNextBlock().lock() : segPair.second.getPrevBlock().lock();
                        if (!nBlk) { allSameNext = false; break; }
                        if (!sharedNext) sharedNext = nBlk;
                        else if (sharedNext != nBlk) { allSameNext = false; break; }
                    }
                    if (!allSameNext) break;
                }

                bool sameCompositionNext = true;
                if (allSameNext && sharedNext) {
                    if (sharedNext->getSequences().size() != blk->getSequences().size()) sameCompositionNext = false;
                    else {
                        for (auto& kv : blk->getSequences()) {
                            if (sharedNext->getSequences().find(kv.first) == sharedNext->getSequences().end()) { sameCompositionNext = false; break; }
                        }
                    }
                }

                if (allSameNext && sharedNext && sharedNext != blk && sameCompositionNext && !processed_this_round.count(sharedNext->getId())) {
                    auto newBlk = concatBlocks(blk, sharedNext);
                    processed_this_round.insert(blk->getId());
                    processed_this_round.insert(sharedNext->getId());
                    processed_this_round.insert(newBlk->getId());
                    mergedOccurred = true;
                    topologyChanged = true;
                }
            }
        }
    }

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    rebuildAllPointers();
}


void BlockSet::refine_new(BlockSet* refSet, BlockSet* qrySet) {
    bool DEBUG_MODE = true;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet Refine New: Sequence-Boundary Driven Cut ===\n"
                              << "============================================================\n";

    // ==========================================
    // 內部工具：重建拓撲指標 (原封不動保留你的神兵利器)
    // ==========================================
    auto rebuildAllPointers = [&]() {
        auto allBlocks = this->getAllBlocks();

        tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto blk = allBlocks[i];
                    blk->clearLinkages();
                    for (auto& seqPair : blk->getSequences()) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            segPair.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                            segPair.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                        }
                    }
                }
            }
        );
        struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
        std::unordered_map<std::string, std::vector<SegRef>> seqTracks;
        for (auto blk : allBlocks) {
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    seqTracks[seqPair.first].push_back({ &segPairInner.second, blk });
                }
            }
        }

        std::vector<std::vector<SegRef>*> trackPtrs;
        trackPtrs.reserve(seqTracks.size());
        for (auto& trackPair : seqTracks) {
            trackPtrs.push_back(&trackPair.second);
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, trackPtrs.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto& track = *(trackPtrs[i]);
                    std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                        return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
                    });
                    for (size_t j = 0; j < track.size(); ++j) {
                        if (j > 0) {
                            auto& prevRef = track[j-1];
                            auto& currRef = track[j];
                            if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                            else prevRef.seg->setPrevBlock(currRef.blk);
                            if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                            else currRef.seg->setNextBlock(prevRef.blk);
                        }
                    }
                }
            }
        );

        for (auto blk : allBlocks) {
            std::set<std::shared_ptr<Block>> prevBlocksSet;
            std::set<std::shared_ptr<Block>> nextBlocksSet;
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    if (auto p = segPair.second.getPrevBlock().lock()) prevBlocksSet.insert(p);
                    if (auto n = segPair.second.getNextBlock().lock()) nextBlocksSet.insert(n);
                }
            }
            for (auto p : prevBlocksSet) blk->addPrevBlock(p);
            for (auto n : nextBlocksSet) blk->addNextBlock(n);
        }
    };

    // ==========================================
    // 階段 1: 靠 Segment 自身端點與長 Gap 找出切割點
    // ==========================================
    std::map<Block::ID, std::set<int>> blockCuts;
    const int LONG_GAP_THRESHOLD = 100; // 認定為非線性的長 Gap 門檻

    for (auto blk : this->getAllBlocks()) {
        int consLen = blk->getConsensus().length();

        for (auto& seqPair : blk->getSequences()) {
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;
                
                // 【核心 1】：拓撲自然斷點 (某條 Sequence 在這塊 Block 的中途加入或離開)
                if (seg.getStart() > 0 && seg.getStart() < consLen) {
                    blockCuts[blk->getId()].insert(seg.getStart());
                }
                if (seg.getEnd() > 0 && seg.getEnd() < consLen) {
                    blockCuts[blk->getId()].insert(seg.getEnd());
                }

                // 【核心 2】：找尋序列不線性的地方 (如你提到的長 Gap)
                for (auto& var : seg.getVariants()) {
                    if (var.getType() == Variation::GAP) {
                        int gapLen = var.getEnd() - var.getStart();
                        if (gapLen >= LONG_GAP_THRESHOLD) {
                            if (var.getStart() > 0 && var.getStart() < consLen) {
                                blockCuts[blk->getId()].insert(var.getStart());
                            }
                            if (var.getEnd() > 0 && var.getEnd() < consLen) {
                                blockCuts[blk->getId()].insert(var.getEnd());
                            }
                        }
                    }
                }
            }
        }
    }

    // ==========================================
    // 階段 2: 執行切割 (附帶小 Block 防呆提醒)
    // ==========================================
    for (auto& kv : blockCuts) {
        Block::ID blkId = kv.first;
        std::vector<int> cuts(kv.second.begin(), kv.second.end());
        std::sort(cuts.rbegin(), cuts.rend()); // ★ 降序排序：由後往前切，保證前面座標不偏移

        auto blk = this->getBlock(blkId);
        if (!blk) continue;
        int prevConsEnd = blk->getConsensus().length();

        for (int cutPos : cuts) {
            int newBlkSize = prevConsEnd - cutPos;
            
            // 【除錯提醒】：如果有小於 100 bp 的碎塊被切出來，警告你！
            if (DEBUG_MODE && newBlkSize > 0 && newBlkSize < 100) {
                std::cout << "[DEBUG] 警告: 正在 Block " << blkId << " 的位置 " << cutPos 
                          << " 下刀，這會切出一個長度僅 " << newBlkSize << " bp 的微小 Block。\n";
            }
            
            this->splitSingleBlock(blkId, cutPos);
            prevConsEnd = cutPos;
        }

        if (DEBUG_MODE && prevConsEnd > 0 && prevConsEnd < 100) {
            std::cout << "[DEBUG] 警告: Block " << blkId << " 被切完後，最左邊殘留的區塊長度僅 " << prevConsEnd << " bp。\n";
        }
    }

    // ==========================================
    // 階段 3: 清理 Pure Gap Segments 與空 Block
    // ==========================================
    if (DEBUG_MODE) std::cout << "  -> Cleaning up pure gap segments and empty blocks...\n";
    for (auto blk : this->getAllBlocks()) {
        if (!blk) continue;
        std::vector<std::string> seqsToRemove;
        
        for (auto& seqPair : blk->getSequences()) {
            std::string seqID = seqPair.first;
            std::vector<int> segsToRemove;
            
            for (auto& segPairInner : seqPair.second.getSegments()) {
                // 如果 Start == End，代表這是一段切完後剩下的 Pure Gap，必須移除
                if (segPairInner.second.getStart() == segPairInner.second.getEnd()) {
                    segsToRemove.push_back(segPairInner.first); 
                }
            }
            
            for (int sCoord : segsToRemove) {
                seqPair.second.getSegments().erase(sCoord);
            }
            
            if (seqPair.second.getSegments().empty()) {
                seqsToRemove.push_back(seqID);
            }
        }
        
        for (const auto& s : seqsToRemove) {
            blk->getSequences().erase(s);
        }
    }

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty() || blk->getConsensus().empty()) {
            emptyBlocks.push_back(blk->getId());
        }
    }
    for (auto id : emptyBlocks) {
        this->deleteBlock(id);
    }

    // ==========================================
    // 階段 4: 重新連接 Graph Edge
    // ==========================================
    rebuildAllPointers(); 

    if (DEBUG_MODE) std::cout << "=== BlockSet Refine New: Completed ===\n\n";
}

std::vector<Block::ID> BlockSet::absorbMicroBlocks() {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << ">>> Executing Fast Micro-Block Absorption...\n";

    // 🌟 新增：用來記錄哪些大 Block 被「污染」了
    std::unordered_set<Block::ID> modified_hosts;

    bool absorptionOccurred = true;
    while (absorptionOccurred) {
        absorptionOccurred = false;
        auto currentBlocks = this->getAllBlocks();

        std::vector<std::shared_ptr<Block>> microBlocks;
        for (auto blk : currentBlocks) {
            if (!this->getBlock(blk->getId())) continue;
            int len = blk->getConsensus().length();
            if (len > 0 && len < 50) microBlocks.push_back(blk);
        }

        std::unordered_map<std::string, std::unordered_map<int, std::pair<std::shared_ptr<Block>, Segment>>> hostStartMap;
        std::unordered_map<std::string, std::unordered_map<int, std::pair<std::shared_ptr<Block>, Segment>>> hostEndMap;

        for (auto b : currentBlocks) {
            if (!this->getBlock(b->getId())) continue;
            for (auto& seqPair : b->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    hostStartMap[seqPair.first][segPair.second.getStart()] = {b, segPair.second};
                    hostEndMap[seqPair.first][segPair.second.getEnd()] = {b, segPair.second};
                }
            }
        }

        for (auto mBlk : microBlocks) {
            if (!this->getBlock(mBlk->getId())) continue;
            if (mBlk->getSequences().empty()) continue;

            std::vector<std::string> seqIDs;
            for (auto& kv : mBlk->getSequences()) seqIDs.push_back(kv.first);

            for (const auto& mSeqID : seqIDs) {
                if (mBlk->getSequences().count(mSeqID) == 0) continue;
                auto& mSegMap = mBlk->getSequences().at(mSeqID).getSegments();
                if (mSegMap.empty()) continue;

                Segment mSeg = mSegMap.begin()->second;
                std::shared_ptr<Block> hostToMatch = nullptr;
                Segment hSeg;
                bool hostIsBeforeMicro = false;

                if (hostEndMap[mSeqID].count(mSeg.getStart())) {
                    hostToMatch = hostEndMap[mSeqID][mSeg.getStart()].first;
                    hSeg = hostEndMap[mSeqID][mSeg.getStart()].second;
                    hostIsBeforeMicro = true;
                } else if (hostStartMap[mSeqID].count(mSeg.getEnd())) {
                    hostToMatch = hostStartMap[mSeqID][mSeg.getEnd()].first;
                    hSeg = hostStartMap[mSeqID][mSeg.getEnd()].second;
                    hostIsBeforeMicro = false;
                }

                if (hostToMatch && this->getBlock(hostToMatch->getId()) && hostToMatch != mBlk) {

                    // 💡 [O(1) 更新]：拔除舊座標索引
                    hostStartMap[mSeqID].erase(hSeg.getStart());
                    hostEndMap[mSeqID].erase(hSeg.getEnd());
                    hostStartMap[mSeqID].erase(mSeg.getStart());
                    hostEndMap[mSeqID].erase(mSeg.getEnd());

                    std::string mSeq = "";
                    std::string mCons = mBlk->getConsensus();
                    for (int i = 0; i < mCons.length(); ++i) {
                        bool isGap = false;
                        char c = mCons[i];
                        for (auto& v : mSeg.getVariants()) {
                            if (v.getType() == Variation::GAP && i >= v.getStart() && i < v.getEnd()) {
                                isGap = true; break;
                            } else if (v.getType() == Variation::SNV && i == v.getStart()) {
                                c = v.getAlt();
                            }
                        }
                        if (!isGap) mSeq += c;
                    }

                    if (mSeg.isReverse() != hSeg.isReverse()) {
                        std::reverse(mSeq.begin(), mSeq.end());
                        for (char& c : mSeq) {
                            if (c == 'A') c = 'T'; else if (c == 'T') c = 'A';
                            else if (c == 'C') c = 'G'; else if (c == 'G') c = 'C';
                        }
                    }

                    int mLen = mSeq.length();
                    if (mLen == 0) {
                        mBlk->getSequences().at(mSeqID).getSegments().erase(mSeg.getStart());
                        if (mBlk->getSequences().at(mSeqID).getSegments().empty()) mBlk->getSequences().erase(mSeqID);
                        continue;
                    }

                    bool isFront = (!hSeg.isReverse()) ? !hostIsBeforeMicro : hostIsBeforeMicro;
                    std::string hostCons = hostToMatch->getConsensus();
                    auto& realHSeg = hostToMatch->getSequences().at(mSeqID).getSegments().at(hSeg.getStart());
                    auto& vars = realHSeg.getVariants();

                    if (isFront) {
                        int gap_end = 0;
                        if (!vars.empty() && vars.front().getType() == Variation::GAP && vars.front().getStart() == 0) gap_end = vars.front().getEnd();
                        int available_gap = gap_end;
                        int offset = 0;

                        if (mLen <= available_gap) {
                            offset = available_gap - mLen;
                            if (mLen == available_gap) vars.erase(vars.begin());
                            else vars.front().setEnd(available_gap - mLen);
                        } else {
                            int overflow = mLen - available_gap;
                            hostCons = mSeq.substr(0, overflow) + hostCons;
                            hostToMatch->setConsensus(hostCons);

                            for (auto& seqPair : hostToMatch->getSequences()) {
                                for (auto& segPair : seqPair.second.getSegments()) {
                                    Segment& s = segPair.second;
                                    auto& s_vars = s.getVariants();
                                    for (auto& v : s_vars) v.shift(overflow);
                                    if (&s != &realHSeg) {
                                        if (!s_vars.empty() && s_vars.front().getType() == Variation::GAP && s_vars.front().getStart() == overflow) {
                                            s_vars.front().setStart(0);
                                        } else s_vars.insert(s_vars.begin(), Variation::createGap(0, overflow));
                                    }
                                }
                            }
                            if (available_gap > 0 && !vars.empty() && vars.front().getType() == Variation::GAP && vars.front().getStart() == overflow) {
                                vars.erase(vars.begin());
                            }
                        }

                        for (int i = 0; i < mLen; ++i) {
                            int c_pos = offset + i;
                            if (mSeq[i] != hostCons[c_pos]) vars.push_back(Variation(c_pos, mSeq[i]));
                        }

                    } else { 
                        int hLen = hostCons.length();
                        int gap_start = hLen;
                        if (!vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == hLen) gap_start = vars.back().getStart();
                        int available_gap = hLen - gap_start;
                        int offset = gap_start;

                        if (mLen <= available_gap) {
                            if (mLen == available_gap) vars.pop_back();
                            else vars.back().setStart(gap_start + mLen);
                        } else {
                            int overflow = mLen - available_gap;
                            hostCons += mSeq.substr(available_gap, overflow);
                            hostToMatch->setConsensus(hostCons);

                            for (auto& seqPair : hostToMatch->getSequences()) {
                                for (auto& segPair : seqPair.second.getSegments()) {
                                    Segment& s = segPair.second;
                                    auto& s_vars = s.getVariants();
                                    if (&s != &realHSeg) {
                                        if (!s_vars.empty() && s_vars.back().getType() == Variation::GAP && s_vars.back().getEnd() == hLen) {
                                            s_vars.back().setEnd(hLen + overflow);
                                        } else s_vars.push_back(Variation::createGap(hLen, hLen + overflow));
                                    }
                                }
                            }
                            if (available_gap > 0 && !vars.empty() && vars.back().getType() == Variation::GAP && vars.back().getEnd() == hLen) {
                                vars.pop_back();
                            }
                        }

                        for (int i = 0; i < mLen; ++i) {
                            int c_pos = offset + i;
                            if (mSeq[i] != hostCons[c_pos]) vars.push_back(Variation(c_pos, mSeq[i]));
                        }

                    }

                    std::shared_ptr<Block> microLeftBlock  = (!mSeg.isReverse()) ? mSeg.getPrevBlock().lock() : mSeg.getNextBlock().lock();
                    std::shared_ptr<Block> microRightBlock = (!mSeg.isReverse()) ? mSeg.getNextBlock().lock() : mSeg.getPrevBlock().lock();

                    if (hostIsBeforeMicro) {
                        // 情況 A: Host 在左，Micro 在右。Host 的實體右端 (End) 要延伸。
                        realHSeg.setEnd(mSeg.getEnd());
                        
                        // Host 需要接管 Micro 右邊的積木，接在自己的「實體右端」
                        if (!realHSeg.isReverse()) {
                            realHSeg.setNextBlock(microRightBlock); // 正向的右端是 Next
                        } else {
                            realHSeg.setPrevBlock(microRightBlock); // 反向的右端是 Prev
                        }
                    } else {
                        // 情況 B: Micro 在左，Host 在右。Host 的實體左端 (Start) 要延伸。
                        realHSeg.setStart(mSeg.getStart());
                        
                        // Host 需要接管 Micro 左邊的積木，接在自己的「實體左端」
                        if (!realHSeg.isReverse()) {
                            realHSeg.setPrevBlock(microLeftBlock);  // 正向的左端是 Prev
                        } else {
                            realHSeg.setNextBlock(microLeftBlock);  // 反向的左端是 Next
                        }
                    }
                    std::sort(vars.begin(), vars.end(), [](Variation& a, Variation& b){
                        if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                        return a.getType() > b.getType();
                    });
                    std::vector<Variation> mergedVars;
                    for (auto& v : vars) {
                        if (mergedVars.empty()) mergedVars.push_back(v);
                        else {
                            auto& last = mergedVars.back();
                            if (last.getType() == Variation::GAP && v.getType() == Variation::GAP && last.getEnd() >= v.getStart()) {
                                int mStart = last.getStart();
                                int mEnd = std::max(last.getEnd(), v.getEnd());
                                mergedVars.pop_back();
                                mergedVars.push_back(Variation::createGap(mStart, mEnd));
                            } else mergedVars.push_back(v);
                        }
                    }
                    vars = mergedVars;

                    int newStart = realHSeg.getStart();
                    if (newStart != hSeg.getStart()) {
                        auto& hostSegsMap = hostToMatch->getSequences().at(mSeqID).getSegments();
                        Segment movedSeg = std::move(realHSeg);
                        hostSegsMap.erase(hSeg.getStart());
                        hostSegsMap[newStart] = std::move(movedSeg);
                        
                        hostStartMap[mSeqID][newStart] = {hostToMatch, hostSegsMap[newStart]};
                        hostEndMap[mSeqID][hostSegsMap[newStart].getEnd()] = {hostToMatch, hostSegsMap[newStart]};
                    } else {
                        hostStartMap[mSeqID][realHSeg.getStart()] = {hostToMatch, realHSeg};
                        hostEndMap[mSeqID][realHSeg.getEnd()] = {hostToMatch, realHSeg};
                    }

                    if (auto pB = mSeg.getPrevBlock().lock()) {
                        if (pB->getSequences().count(mSeqID)) {
                            for (auto& pSegPair : pB->getSequences().at(mSeqID).getSegments()) {
                                if (pSegPair.second.getNextBlock().lock() == mBlk) pSegPair.second.setNextBlock(hostToMatch);
                                else if (pSegPair.second.getPrevBlock().lock() == mBlk) pSegPair.second.setPrevBlock(hostToMatch);
                            }
                        }
                    }
                    if (auto nB = mSeg.getNextBlock().lock()) {
                        if (nB->getSequences().count(mSeqID)) {
                            for (auto& nSegPair : nB->getSequences().at(mSeqID).getSegments()) {
                                if (nSegPair.second.getPrevBlock().lock() == mBlk) nSegPair.second.setPrevBlock(hostToMatch);
                                else if (nSegPair.second.getNextBlock().lock() == mBlk) nSegPair.second.setNextBlock(hostToMatch);
                            }
                        }
                    }

                    mBlk->getSequences().at(mSeqID).getSegments().erase(mSeg.getStart());
                    if (mBlk->getSequences().at(mSeqID).getSegments().empty()) {
                        mBlk->getSequences().erase(mSeqID);
                    }

                    // 🌟 核心記錄：將這個成功吸收別人的大 Block 標記為 Dirty
                    modified_hosts.insert(hostToMatch->getId());

                    absorptionOccurred = true;
                }
            } 

            if (mBlk->getSequences().empty()) {
                this->deleteBlock(mBlk->getId());
            }
        }
    }

    // 🌟 結束前：將 Set 轉換成 Vector 回傳
    return std::vector<Block::ID>(modified_hosts.begin(), modified_hosts.end());
}

void BlockSet::absorbSingletons() {
    bool DEBUG_MODE = true;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet Phase 4: Singleton Absorber (Verbose Mode) ===\n"
                              << "============================================================\n";

    // ==========================================
    // 內部工具 1：還原未對齊的真實序列 (考慮 Strand)
    // ==========================================
    auto getRawSeq = [](const std::string& cons, Segment& seg) {
        std::string raw = "";
        int pos = std::min(seg.getStart(), seg.getEnd());
        int end = std::max(seg.getStart(), seg.getEnd());
        
        for (auto var : seg.getVariants()) {
            if (var.getStart() >= end || var.getEnd() <= pos) continue;
            if (var.getStart() > pos) raw += cons.substr(pos, var.getStart() - pos);
            if (var.getType() == Variation::SNV) {
                raw += var.getAlt();
                pos = var.getStart() + 1;
            } else if (var.getType() == Variation::GAP) {
                pos = var.getEnd();
            }
        }
        if (pos < end) raw += cons.substr(pos, end - pos);
        
        if (seg.isReverse()) {
            std::reverse(raw.begin(), raw.end());
            for(char& c : raw) {
                if (c=='A') c='T'; else if(c=='T') c='A'; else if(c=='C') c='G'; else if(c=='G') c='C';
            }
        }
        return raw;
    };

    // ==========================================
    // 內部工具 2：擴充 Consensus 並同步推移所有變異
    // ==========================================
    auto expandConsensus = [&](std::shared_ptr<Block> blk, int pos, const std::string& insSeq, const std::string& targetSeqID) {
        int insLen = insSeq.length();
        std::string cons = blk->getConsensus();
        cons.insert(pos, insSeq);
        blk->setConsensus(cons);

        for (auto& seqPair : blk->getSequences()) {
            std::string sID = seqPair.first;
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;
                auto& vars = seg.getVariants();
                
                for (auto& v : vars) {
                    if (v.getStart() >= pos) {
                        v.shift(insLen);
                    } else if (v.getType() == Variation::GAP && v.getStart() < pos && v.getEnd() > pos) {
                        v.setEnd(v.getEnd() + insLen);
                    }
                }

                if (sID != targetSeqID) {
                    bool coveredByGap = false;
                    for (auto& v : vars) {
                        if (v.getType() == Variation::GAP && v.getStart() <= pos && v.getEnd() >= (pos + insLen)) {
                            coveredByGap = true; break;
                        }
                    }
                    if (!coveredByGap) {
                        vars.push_back(Variation::createGap(pos, pos + insLen));
                    }
                }
                
                if (seg.getEnd() > pos) seg.setEnd(seg.getEnd() + insLen);
                if (seg.getStart() > pos) seg.setStart(seg.getStart() + insLen);
            }
        }
    };

    // ==========================================
    // 主迴圈：不斷尋找並吸收 Singleton
    // ==========================================
    bool absorptionOccurred = true;
    while (absorptionOccurred) {
        absorptionOccurred = false;
        auto currentBlocks = this->getAllBlocks();

        for (auto mBlk : currentBlocks) {
            if (!this->getBlock(mBlk->getId())) continue;
            if (mBlk->getSequences().size() != 1) continue; // 目標鎖定：孤兒 Singleton

            std::string seqID = mBlk->getSequences().begin()->first;
            Segment mSeg = mBlk->getSequences().begin()->second.getSegments().begin()->second;
            std::string mSeqStr = getRawSeq(mBlk->getConsensus(), mSeg);
            if (mSeqStr.empty()) continue;

            if (DEBUG_MODE) {
                std::cout << "--------------------------------------------------------\n";
                std::cout << "[Target] 發現 Singleton! Block ID=" << mBlk->getId() 
                          << " | Seq=" << seqID << " | 真實長度=" << mSeqStr.length() << " bp\n";
            }

            // 嘗試吸收的內部函式
            auto attemptAbsorb = [&](std::shared_ptr<Block> hBlk, bool isHostPrev) -> bool {
                std::string direction = isHostPrev ? "Prev" : "Next";
                
                if (!hBlk || !this->getBlock(hBlk->getId())) {
                    if (DEBUG_MODE) std::cout << "  -> 嘗試 " << direction << " Host: 失敗 (Host Block 缺失或無效)\n";
                    return false;
                }
                
                if (DEBUG_MODE) std::cout << "  -> 嘗試 " << direction << " Host (Block ID=" << hBlk->getId() << ")...\n";

                if (hBlk->getSequences().count(seqID) == 0) {
                    if (DEBUG_MODE) std::cout << "    x 失敗: Host 中找不到相同的 Sequence (" << seqID << ")\n";
                    return false;
                }

                Segment* hSegPtr = nullptr;
                for (auto& segPair : hBlk->getSequences().at(seqID).getSegments()) {
                    Segment& s = segPair.second;
                    if (isHostPrev) {
                        if ((!s.isReverse() && s.getNextBlock().lock() == mBlk) || (s.isReverse() && s.getPrevBlock().lock() == mBlk)) { hSegPtr = &s; break; }
                    } else {
                        if ((!s.isReverse() && s.getPrevBlock().lock() == mBlk) || (s.isReverse() && s.getNextBlock().lock() == mBlk)) { hSegPtr = &s; break; }
                    }
                }
                
                if (!hSegPtr) {
                    if (DEBUG_MODE) std::cout << "    x 失敗: Host 中有該序列，但拓撲上並未與 Singleton 直接相連\n";
                    return false;
                }

                int consLen = hBlk->getConsensus().length();
                int targetGapStart = -1, targetGapEnd = -1;
                int gapVarIdx = -1;

                bool matchAtEnd = (!hSegPtr->isReverse() && isHostPrev) || (hSegPtr->isReverse() && !isHostPrev);
                auto& vars = hSegPtr->getVariants();
                
                for (size_t i = 0; i < vars.size(); ++i) {
                    if (vars[i].getType() == Variation::GAP) {
                        if (matchAtEnd && vars[i].getEnd() == consLen) {
                            targetGapStart = vars[i].getStart(); targetGapEnd = vars[i].getEnd(); gapVarIdx = i; break;
                        } else if (!matchAtEnd && vars[i].getStart() == 0) {
                            targetGapStart = vars[i].getStart(); targetGapEnd = vars[i].getEnd(); gapVarIdx = i; break;
                        }
                    }
                }

                if (gapVarIdx == -1) {
                    if (DEBUG_MODE) std::cout << "    x 失敗: 成功接壤，但邊界處沒有找到 Variation::GAP 可供吸收\n";
                    return false;
                }

                int gapLen = targetGapEnd - targetGapStart;
                if (DEBUG_MODE) std::cout << "    v 找到邊界 GAP: [" << targetGapStart << " - " << targetGapEnd << "] (長度: " << gapLen << " bp)\n";

                std::string gapCons = hBlk->getConsensus().substr(targetGapStart, gapLen);
                std::string alignedQry = mSeqStr;
                
                if (hSegPtr->isReverse() != mSeg.isReverse()) {
                    std::reverse(alignedQry.begin(), alignedQry.end());
                    for(char& c : alignedQry) {
                        if (c=='A') c='T'; else if(c=='T') c='A'; else if(c=='C') c='G'; else if(c=='G') c='C';
                    }
                }

                if (DEBUG_MODE) std::cout << "    - 執行 Alignment (GAP vs Singleton)... ";
                
                AlnResult res = runSemiGlobalAlignment(gapCons, alignedQry);
                
                if (!res.success) {
                    if (DEBUG_MODE) std::cout << "失敗 (對齊演算法回傳錯誤)\n";
                    return false;
                }
                if (res.identity < 0.6f) {
                    if (DEBUG_MODE) std::cout << "失敗 (Identity: " << res.identity << " < 0.60 門檻)\n";
                    return false;
                }

                if (DEBUG_MODE) std::cout << "成功! (Identity: " << res.identity << ")\n";

                auto cigar = mga::parser::parseCigar(res.cigar);
                
                int consR = 0;
                for (auto op : cigar) if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') consR += op.first;
                if (consR < gapCons.length()) cigar.push_back({gapCons.length() - consR, 'D'});

                vars.erase(vars.begin() + gapVarIdx);

                int curRefPos = targetGapStart;
                int curQryPos = 0;
                std::vector<Variation> newVars;

                for (auto op : cigar) {
                    int len = op.first; char t = op.second;
                    if (t == 'S' || t == 'H') { curQryPos += len; continue; }
                    
                    if (t == 'M' || t == '=' || t == 'X') {
                        for (int i = 0; i < len; ++i) {
                            if (gapCons[curRefPos - targetGapStart] != alignedQry[curQryPos]) {
                                newVars.push_back(Variation(curRefPos, alignedQry[curQryPos]));
                            }
                            curRefPos++; curQryPos++;
                        }
                    } else if (t == 'D') {
                        newVars.push_back(Variation::createGap(curRefPos, curRefPos + len));
                        curRefPos += len;
                    } else if (t == 'I') {
                        if (DEBUG_MODE) std::cout << "    ! 觸發 Consensus 擴展 (長度: " << len << " bp)\n";
                        std::string insSeq = alignedQry.substr(curQryPos, len);
                        expandConsensus(hBlk, curRefPos, insSeq, seqID);
                        curRefPos += len;
                        curQryPos += len;
                    }
                }

                vars.insert(vars.end(), newVars.begin(), newVars.end());
                
                std::shared_ptr<Block> microLeftBlock  = (!mSeg.isReverse()) ? mSeg.getPrevBlock().lock() : mSeg.getNextBlock().lock();
                std::shared_ptr<Block> microRightBlock = (!mSeg.isReverse()) ? mSeg.getNextBlock().lock() : mSeg.getPrevBlock().lock();

                if (isHostPrev) {
                    if (!hSegPtr->isReverse()) hSegPtr->setNextBlock(microRightBlock);
                    else hSegPtr->setPrevBlock(microRightBlock);
                } else {
                    if (!hSegPtr->isReverse()) hSegPtr->setPrevBlock(microLeftBlock);
                    else hSegPtr->setNextBlock(microLeftBlock);
                }
                
                mBlk->getSequences().erase(seqID);
                
                if (DEBUG_MODE) std::cout << "    >>> ★ 吸收完成! 成功將 Singleton 融入 Host ★\n";
                return true;
            };

            std::shared_ptr<Block> pBlk = (!mSeg.isReverse()) ? mSeg.getPrevBlock().lock() : mSeg.getNextBlock().lock();
            std::shared_ptr<Block> nBlk = (!mSeg.isReverse()) ? mSeg.getNextBlock().lock() : mSeg.getPrevBlock().lock();
            
            if (pBlk && attemptAbsorb(pBlk, true)) { absorptionOccurred = true; continue; }
            if (nBlk && attemptAbsorb(nBlk, false)) { absorptionOccurred = true; continue; }
            
            if (DEBUG_MODE) std::cout << "  -> [結論] 兩側 Host 皆無法吸收，保留此 Singleton。\n";
        }
    }

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    
    if (DEBUG_MODE && !emptyBlocks.empty()) {
        std::cout << "--------------------------------------------------------\n";
        std::cout << "-> 正在清除 " << emptyBlocks.size() << " 個已被掏空的 Singleton Blocks...\n";
    }
    
    for (auto id : emptyBlocks) this->deleteBlock(id);

    auto rebuildAllPointers = [&]() {
        auto allBlocks = this->getAllBlocks();

        tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto blk = allBlocks[i];
                    blk->clearLinkages();
                    for (auto& seqPair : blk->getSequences()) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            segPair.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                            segPair.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                        }
                    }
                }
            }
        );
        
        struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
        std::unordered_map<std::string, std::vector<SegRef>> seqTracks;
        for (auto blk : allBlocks) {
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    seqTracks[seqPair.first].push_back({ &segPairInner.second, blk });
                }
            }
        }

        std::vector<std::vector<SegRef>*> trackPtrs;
        trackPtrs.reserve(seqTracks.size());
        for (auto& trackPair : seqTracks) {
            trackPtrs.push_back(&trackPair.second);
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, trackPtrs.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    auto& track = *(trackPtrs[i]);
                    std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                        return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
                    });
                    for (size_t j = 0; j < track.size(); ++j) {
                        if (j > 0) {
                            auto& prevRef = track[j-1];
                            auto& currRef = track[j];
                            if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                            else prevRef.seg->setPrevBlock(currRef.blk); 
                            if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                            else currRef.seg->setNextBlock(prevRef.blk); 
                        }
                    }
                }
            }
        );

        for (auto blk : allBlocks) {
            std::set<std::shared_ptr<Block>> prevBlocksSet;
            std::set<std::shared_ptr<Block>> nextBlocksSet;
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    if (auto p = segPair.second.getPrevBlock().lock()) prevBlocksSet.insert(p);
                    if (auto n = segPair.second.getNextBlock().lock()) nextBlocksSet.insert(n);
                }
            }
            for (auto p : prevBlocksSet) blk->addPrevBlock(p);
            for (auto n : nextBlocksSet) blk->addNextBlock(n);
        }
    };
    
    rebuildAllPointers();
    if (DEBUG_MODE) std::cout << "=== Singleton Absorber Completed ===\n\n";
}
*/