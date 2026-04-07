#include "mga.hpp"

// =======================
// BlockSet Implementation
// =======================

BlockSet::BlockSet(SetId id) : id_(id) {}

BlockSet::SetId BlockSet::getId() const { return id_; }

std::shared_ptr<Block> BlockSet::createBlock(const std::string& consensus) {
    Block::ID new_id = next_block_id_++;
    auto new_block = std::make_shared<Block>(new_id, consensus);
    blocks_[new_id] = new_block;
    invalidateRepCache();
    return new_block;
}

std::shared_ptr<Block> BlockSet::getBlock(Block::ID id) {
    auto it = blocks_.find(id);
    if (it != blocks_.end()) {
        return it->second;
    }
    return nullptr;
}

std::vector<std::shared_ptr<Block>> BlockSet::getAllBlocks() {
    std::vector<std::shared_ptr<Block>> all_blocks;
    all_blocks.reserve(blocks_.size());
    for (const auto& pair : blocks_) {
        all_blocks.push_back(pair.second);
    }
    return all_blocks;
}

void BlockSet::updateLongestSequence(std::unordered_map<std::string, int>& sequence_lengths) {
    int longestLen = -1;
    for (auto& seq: seqs) {
        if (sequence_lengths[seq] > longestLen) {
            longest_sequence_ = seq;
        }
    }
}

/*
std::vector<Block::ID> BlockSet::getRepresentativeBlocks() {
    std::vector<Block::ID> representative_blocks;
    std::unordered_set<Block::ID> visited;

    bool debug = false;

    std::cout << "\n[DEBUG-REP] === Extracting Representative Blocks for [" << longest_sequence_ << "] ===\n";

    // ==========================================
    // 步驟 1: 找出 longest_sequence_ 真正的起點 (最小座標)
    // ==========================================
    std::shared_ptr<Block> start_block = nullptr;
    int min_start_coord = std::numeric_limits<int>::max();

    for (const auto& pair : blocks_) {
        auto blk = pair.second;
        auto& sequences = blk->getSequences();
        
        auto seq_it = sequences.find(longest_sequence_);
        if (seq_it != sequences.end()) {
            for (const auto& seg_pair : seq_it->second.getSegments()) {
                if (seg_pair.first < min_start_coord) {
                    min_start_coord = seg_pair.first;
                    start_block = blk;
                }
            }
        }
    }

    if (!start_block) {
        std::cerr << "[DEBUG-REP] ERROR: Sequence [" << longest_sequence_ << "] not found in any block.\n"; 
        return representative_blocks;
    }

    if (debug) std::cout << "[DEBUG-REP] Global Start Coordinate: " << min_start_coord << " in Block ID: " << start_block->getId() << "\n";

    // ==========================================
    // 步驟 2: 從起點開始，順流而下 (Forward Traversal)
    // ==========================================
    std::shared_ptr<Block> current_block = start_block;
    int current_coord = min_start_coord;
    int total_traversed_length = 0;

    if (debug) std::cout << "[DEBUG-REP] Traversal Path:\n";

    while (current_block) {
        auto& sequences = current_block->getSequences();
        auto seq_it = sequences.find(longest_sequence_);
        
        if (seq_it == sequences.end()) {
            std::cerr << "[DEBUG-REP] WARNING: Graph broken. Sequence missing in Block ID: " << current_block->getId() << "\n";
            break;
        }

        auto& segments = seq_it->second.getSegments();
        auto seg_it = segments.find(current_coord);
        
        if (seg_it == segments.end()) {
            std::cerr << "[DEBUG-REP] WARNING: Sequence gap or coordinate mismatch. Looking for coord " 
                      << current_coord << " in Block ID: " << current_block->getId() << "\n";
            break; 
        }

        Segment& current_segment = seg_it->second;
        int seg_length = std::abs(current_segment.getEnd() - current_segment.getStart());

        if (visited.find(current_block->getId()) == visited.end()) {
            visited.insert(current_block->getId());
            representative_blocks.push_back(current_block->getId());
            total_traversed_length += seg_length;
            
            if (debug) std::cout << "  -> Block ID: " << current_block->getId() 
                                 << " | Seg [" << current_segment.getStart() << "->" << current_segment.getEnd() 
                                 << "] | SegLen: " << seg_length << "\n";
        } else {
            if (debug) std::cout << "  -> Block ID: " << current_block->getId() << " (Already visited, skipping loop/paralog)\n";
        }

        current_coord = current_segment.getEnd();
        current_block = current_segment.getNextBlock().lock(); 
    }

    std::cout << "[DEBUG-REP] === Traversal Summary ===\n";
    std::cout << "  - Total Representative Blocks: " << representative_blocks.size() << "\n";
    std::cout << "  - Total Sequence Length Traversed: " << total_traversed_length << " bp\n";
    std::cout << "------------------------------------------------------------\n";

    return representative_blocks;
}
*/

std::vector<Block::ID> BlockSet::getRepresentativeBlocks() {
    // 1. 如果已經計算過，直接一秒回傳快取結果！
    if (is_rep_cached_) {
        return representative_cache_;
    }

    bool debug = true;
    if (debug) std::cout << "\n[DEBUG-REP] === Extracting Representative Blocks (Fast Sort) ===\n";

    std::vector<std::pair<int, Block::ID>> backbone_candidates;
    std::vector<Block::ID> remaining_blocks;

    // 2. 走訪全圖 (線性掃描，不走複雜指標)
    for (auto& pair : blocks_) {
        Block::ID id = pair.first;
        auto blk = pair.second;
        
        auto seq_it = blk->getSequences().find(longest_sequence_);
        if (seq_it != blk->getSequences().end()) {
            // 如果是 Backbone，找出它在這條序列上最左邊的座標
            int min_coord = std::numeric_limits<int>::max();
            for (auto& seg_pair : seq_it->second.getSegments()) {
                int start = std::min(seg_pair.second.getStart(), seg_pair.second.getEnd());
                if (start < min_coord) min_coord = start;
            }
            backbone_candidates.push_back({min_coord, id});
        } else {
            // 如果沒有包含最長序列，歸類為 Remaining
            remaining_blocks.push_back(id);
        }
    }

    // 3. 極速排序
    // Backbone 依照基因體座標排序 (還原真實物理順序)
    std::sort(backbone_candidates.begin(), backbone_candidates.end());
    // Remaining 依照 ID 排序 (確保每次執行結果固定不變)
    std::sort(remaining_blocks.begin(), remaining_blocks.end());

    // 4. 寫入快取
    representative_cache_.clear();
    representative_cache_.reserve(backbone_candidates.size() + remaining_blocks.size());

    for (const auto& item : backbone_candidates) {
        representative_cache_.push_back(item.second);
    }
    for (Block::ID id : remaining_blocks) {
        representative_cache_.push_back(id);
    }

    // 標記快取為有效
    is_rep_cached_ = true;

    if (debug) {
        std::cout << "[DEBUG-REP] === Fast Traversal Summary ===\n";
        std::cout << "  - Backbone Blocks: " << backbone_candidates.size() << "\n";
        std::cout << "  - Appended Remaining Blocks: " << remaining_blocks.size() << "\n";
        std::cout << "  - Total Blocks: " << representative_cache_.size() << "\n";
        std::cout << "------------------------------------------------------------\n";
    }

    return representative_cache_;
}

std::vector<std::shared_ptr<Block>> BlockSet::getRemainingBlocks() {
    std::vector<std::shared_ptr<Block>> remaining_blocks;
    std::set<Block::ID> consensus_ids;
    auto representative_blocks = this->getRepresentativeBlocks();
    for (const auto& id : representative_blocks) {
        consensus_ids.insert(id);
    }
    for (const auto& pair : blocks_) {
        if (consensus_ids.find(pair.first) == consensus_ids.end()) {
            remaining_blocks.push_back(pair.second);
        }
    }
    return remaining_blocks;
}

bool BlockSet::deleteBlock(Block::ID id) {
    invalidateRepCache();
    return blocks_.erase(id) > 0;
}

void BlockSet::getRepresentativeAndRemaining(std::vector<std::pair<std::string, std::string>>& representative, std::vector<std::pair<std::string, std::string>>& remaining) {
    representative.clear();
    remaining.clear();
    std::string main_name = id_ + "_main";
    std::string main_sequence = "";
    std::vector<Block::ID> remaining_id;

    auto representative_blocks = this->getRepresentativeBlocks();
    std::set<Block::ID> repBlocks(representative_blocks.begin(), representative_blocks.end());
    for (auto& blocks: this->blocks_) {
        if (repBlocks.find(blocks.first) != repBlocks.end()) {
            main_sequence += blocks.second->getConsensus();
        } else {
            remaining_id.push_back(blocks.first);
        }
    }

    representative.push_back({main_name, main_sequence});
    std::sort(remaining_id.begin(), remaining_id.end());

    for (auto& blockID : remaining_id) {
        std::string name = id_ + "_" + std::to_string(blockID);
        remaining.push_back({name, blocks_[blockID]->getConsensus()});
    }
    return;
}

// ==========================================
// 清空所有 Blocks 並安全釋放記憶體
// ==========================================
void BlockSet::clearBlocks() {
    blocks_.clear();
    invalidateRepCache();
    next_block_id_ = 1; 
}

// ==========================================
// 將外部的 Block 加入此 BlockSet 並賦予新 ID
// ==========================================
std::shared_ptr<Block> BlockSet::addBlock(std::shared_ptr<Block> oldBlock) {
    if (!oldBlock) return 0; // 防呆機制

    // 1. 取得專屬於這個 BlockSet 的新 ID
    Block::ID newId = next_block_id_++;

    // 2. 利用舊 Block 的 Consensus 建立全新的 Block
    auto newBlock = std::make_shared<Block>(newId, oldBlock->getConsensus());

    // 3. 深拷貝：將所有的 SequenceInfo 複製過去
    for (const auto& seqPair : oldBlock->getSequences()) {
        newBlock->addSequence(seqPair.second);
    }

    // 4. 拷貝屬性狀態 (例如是否來自 Primary)
    newBlock->setFromPrimary(oldBlock->isFromPrimary());

    // 註：如果你的老 Block 裡面有 duplications_ 或 paralogs_ 需要拷貝，
    // 需要在 Block class 補上 getter 才能在這裡一起 copy 過去。
    // 否則新 Block 預設這兩個 set 是空的，通常在重新 merge 的情境下空的是合理的。

    // 5. 註冊進這個 BlockSet 的 Dictionary 中
    blocks_[newId] = newBlock;
    invalidateRepCache();
    return blocks_[newId]; // 回傳被賦予的全新 ID
}


void BlockSet::print(std::ostream& os) const {
    os << "\n";
    os << "------------------------------------------------------------\n";
    os << ">>> BlockSet ID: " << id_ << " (Contains " << blocks_.size() << " blocks)\n";
    os << "------------------------------------------------------------\n";
    
    for (const auto& pair : blocks_) {
        if (pair.second) {
            pair.second->print(os);
        }
    }
}

void BlockSet::selfMapping(Option& option) {

    bool debug = false;

    std::string temp_dir = option.tempDir;
    std::string seqName = blocks_.begin()->second->getSequences().begin()->first;
    std::string seq = blocks_.begin()->second->getConsensus();
    std::string seqFile = temp_dir + "/" +  seqName + ".fa";
    std::string pafFile = temp_dir + "/" +  seqName + "-self.paf"; // 建議副檔名改為 paf
    mga::stringPairVec seqPair (1, {seqName, seq});
    
    // 【新增統計】：記錄原始長度
    int total_original_length = seq.length();
    
    std::cout << "[SelfMapping] 0. Writing and Self-aligning [" << seqName << "] (Len: " << total_original_length << " bp)...\n";
    mga::io::writeAlignment(seqFile, seqPair, false, false);

    std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
    std::string command;
    int system_ret; 

    command = minimap2_path + " -cx asm5 -g 500 -r 500,500 -N 20 -X " + seqFile + " " + seqFile + " > " + pafFile;
    system_ret = system(command.c_str());
    if (system_ret != 0) {
        std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
    }

    std::cout << "[SelfMapping] 1. Parsing PAF and Deduplicating...\n";
    std::vector<mga::Alignment> rawAlignments = mga::parser::parseMinimap2PAF(pafFile);
    std::vector<mga::Alignment> uniqueAlignments;

    for (size_t i = 0; i < rawAlignments.size(); ++i) {
        if (i % 2 == 0) { 
            if (rawAlignments[i].alnScore > 0) uniqueAlignments.push_back(rawAlignments[i]);
        }
    }

    std::sort(uniqueAlignments.begin(), uniqueAlignments.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        return a.alnScore > b.alnScore; 
    });

    std::cout << "[SelfMapping] 2. Initializing Coordinate Dictionary...\n";

    // 字典：Key 是 start_coordinate，Value 是區間資訊
    std::map<int, SegNode> dict;

    // 一開始，整條 Sequence 就是一個完整的 Segment
    auto initialBlock = blocks_.begin()->second;
    int seqLen = total_original_length; 
    dict[0] = {0, seqLen, initialBlock->getId()};

    auto splitDictSegment = [&](int cutPos) {
        auto it = dict.upper_bound(cutPos); it--;
        if (it->first == cutPos) return; // 已經是邊界，不用切

        int start = it->second.start;
        Block::ID targetBlkId = it->second.blkId;
        std::shared_ptr<Block> targetBlk = this->getBlock(targetBlkId);
        if (!targetBlk) return;

        // 1. 找出這個 Block 對應的 Segment，判斷正反股
        Segment* targetSeg = nullptr;
        for (auto& seqPair : targetBlk->getSequences()) {
            if (seqPair.second.getSegments().count(start)) {
                targetSeg = &(seqPair.second.getSegments()[start]);
                break;
            }
        }
        if (!targetSeg) return;

        // 2. 算好真實的 Consensus Cut
        int origLocalCut = targetSeg->isReverse() ? (targetSeg->getEnd() - cutPos) : (cutPos - targetSeg->getStart());
        int consensusCut = origLocalCut;
        for (auto& v : targetSeg->getVariations()) {
            if (v.getType() == Variation::GAP && v.getStart() <= consensusCut) {
                consensusCut += (v.getEnd() - v.getStart());
            }
        }

        if (consensusCut <= 0 || consensusCut >= targetBlk->getConsensus().length()) return; 

        // 3. 物理切割 Block
        auto parts = this->splitSingleBlock(targetBlkId, consensusCut);
        if (parts.first == -1) return; 

        // 4. 直接呼叫重建！一行解決原本幾十行的字典同步惡夢！
        this->rebuildDictionary(dict, seqName); 
    };

    std::cout << "[SelfMapping] 3. Greedy Interval Merging...\n";

    int total_collapsed_bases = 0;
    int merge_operations_count = 0;

    int alnCount = 0;

    for (const auto& aln : uniqueAlignments) {

        alnCount++;

        int subAlnCount = 0;

        
        int rPos = aln.refIdx.first;
        int qPos = aln.qryIdx.first;
        int alnRefEnd = aln.refIdx.second;
        int alnQryEnd = aln.qryIdx.second;

        // 手動走訪 CIGAR，確保 rPos 和 qPos 絕對同步
        int cigarIdx = 0;
        if (aln.CIGAR.empty()) continue;
        int opRemain = aln.CIGAR[0].first;
        char opType = aln.CIGAR[0].second;

        // 分塊處理這條 Alignment
        while (rPos < alnRefEnd && qPos < alnQryEnd && cigarIdx < aln.CIGAR.size()) {

            subAlnCount++;
            if (subAlnCount > 20) exit(1);
            
            // 1. 查字典，找出目前的 rPos 和 qPos 屬於哪兩個 Segment 區間
            auto rIt = dict.upper_bound(rPos); rIt--;
            auto qIt = dict.upper_bound(qPos); qIt--;
            SegNode rSeg = rIt->second;
            SegNode qSeg = qIt->second;

            // 2. 決定這一個 Chunk 最多能走多遠 (不能超過當前 Segment 的邊界)
            int rRemain = rSeg.end - rPos;
            int qRemain = qSeg.end - qPos;
            
            mga::Cigar fragCigar;
            int chunkRefLen = 0;
            int chunkQryLen = 0;

            while (cigarIdx < aln.CIGAR.size()) {
                bool consumesRef = (opType == 'M' || opType == '=' || opType == 'X' || opType == 'D');
                bool consumesQry = (opType == 'M' || opType == '=' || opType == 'X' || opType == 'I');

                int maxStep = opRemain;
                if (consumesRef && maxStep > rRemain - chunkRefLen) maxStep = rRemain - chunkRefLen;
                if (consumesQry && maxStep > qRemain - chunkQryLen) maxStep = qRemain - chunkQryLen;

                if (maxStep == 0) break; // 撞到某一個 Segment 邊界了，結算這個 Chunk

                if (!fragCigar.empty() && fragCigar.back().second == opType) {
                    fragCigar.back().first += maxStep;
                } else {
                    fragCigar.push_back({maxStep, opType});
                }
                
                if (consumesRef) chunkRefLen += maxStep;
                if (consumesQry) chunkQryLen += maxStep;

                opRemain -= maxStep;
                if (opRemain == 0) {
                    cigarIdx++;
                    if (cigarIdx < aln.CIGAR.size()) {
                        opRemain = aln.CIGAR[cigarIdx].first;
                        opType = aln.CIGAR[cigarIdx].second;
                    }
                }
            } // 結束 CIGAR 讀取

            // ==========================================
            // === 新增的 Debug 輸出：顯示 Chunk 資訊 ===
            // ==========================================
            if (debug) {
                std::cout << "\n------------------------------------------------------------\n";
                std::cout << "[DEBUG] Alignment [" << alnCount << "-" << subAlnCount << "/" << uniqueAlignments.size() << "]\n";
                std::cout << "[DEBUG] Direction: " << (aln.inverse ? "Reverse Complement (-)" : "Forward (+)") << "\n";
                std::cout << "[DEBUG] Ref Raw Chunk: [" << rPos << " -> " << (rPos + chunkRefLen) << ") | Length: " << chunkRefLen << " | Curr Blk ID: " << rSeg.blkId << "\n";
                std::cout << "[DEBUG] Qry Raw Chunk: [" << qPos << " -> " << (qPos + chunkQryLen) << ") | Length: " << chunkQryLen << " | Curr Blk ID: " << qSeg.blkId << "\n";
                std::cout << "[DEBUG] Fragment CIGAR: ";
                for (auto op : fragCigar) std::cout << op.first << op.second;
                std::cout << "\n";
            }

            // 如果已經在同一個 Block 且不同 Segment 裡 (之前 Merge 過了)，跳過
            if (rSeg.blkId == qSeg.blkId) {
                bool skip = false;
                auto r_st = rPos, r_en = rPos + chunkRefLen;
                auto q_st = qPos, q_en = qPos + chunkQryLen;
                if (r_st <= q_en && q_st <= r_en) {
                    if (debug) std::cout << "[DEBUG] Action: Skip (Overlapped) \t" 
                                         << "[" << r_st << "," << r_en << ")/[" << q_st << "," << q_en << ")\n";
                    skip = true;
                }
                else {
                    int dist = (r_st < q_st) ? (q_st - r_en) : (r_st - q_en);
                    bool tooClose = (dist < 100);
                    if (tooClose) {
                        if (debug) std::cout << "[DEBUG] Action: Skip (Too closed) \t" 
                                             << "[" << r_st << "," << r_en << ")/[" << q_st << "," << q_en << ")\n";
                        skip = true;
                    }
                }
                if (!skip && rSeg.start != qSeg.start) {
                    if (debug) std::cout << "[DEBUG] Action: Skip (Already merged in Block " << rSeg.blkId << ")\n";
                    skip = true;
                }
                if (skip) {
                    rPos += chunkRefLen;
                    qPos += chunkQryLen;
                    continue;
                }
            }
            

            // 2. 【核心修改：先過濾短片段】
            // 如果這段實際對齊的片段太短 (< 100bp)，不值得 Merge，直接跳過並推進座標
            if (chunkRefLen < 100 || chunkQryLen < 100) {
                if (debug) std::cout << "[DEBUG] Action: Skip (Too short: Ref " << chunkRefLen << "bp, Qry " << chunkQryLen << "bp)\n";
                rPos += chunkRefLen;
                qPos += chunkQryLen;
                continue;
            }

            // 3. 【確定要 Merge 後，才進行 Overhang 切割與字典更新】
            // (註：splitDictSegment 函式內部已經是「先切 Block，再 Update 字典」的順序了)
            
            if (rPos - rSeg.start > 100) {
                splitDictSegment(rPos);
                rIt = dict.upper_bound(rPos); rIt--; // 重新獲取切完後包含 rPos 的右半部
                rSeg = rIt->second;
            }
            if (rSeg.end - (rPos + chunkRefLen) > 100) {
                splitDictSegment(rPos + chunkRefLen);
                rIt = dict.upper_bound(rPos); rIt--; // 重新獲取切完後包含 rPos 的左半部
                rSeg = rIt->second;
            }

            if (qPos - qSeg.start > 100) {
                splitDictSegment(qPos);
                qIt = dict.upper_bound(qPos); qIt--;
                qSeg = qIt->second;
            }
            if (qSeg.end - (qPos + chunkQryLen) > 100) {
                splitDictSegment(qPos + chunkQryLen);
                qIt = dict.upper_bound(qPos); qIt--;
                qSeg = qIt->second;
            }

            // 4. 取得原始 Overhang (維持原本的寫法)
            int leftRefOverhang = rPos - rSeg.start;
            int rightRefOverhang = rSeg.end - (rPos + chunkRefLen);
            int leftQryOverhang = qPos - qSeg.start;
            int rightQryOverhang = qSeg.end - (qPos + chunkQryLen);

            // 取出精確的 Block 與 Segment
            std::shared_ptr<Block> refBlk = this->getBlock(rSeg.blkId);
            std::shared_ptr<Block> qryBlk = this->getBlock(qSeg.blkId);

            Segment* targetRefSeg = nullptr;
            for (auto& seqPair : refBlk->getSequences()) {
                if (seqPair.second.getSegments().count(rSeg.start)) {
                    targetRefSeg = &(seqPair.second.getSegments()[rSeg.start]);
                    break;
                }
            }
            Segment* targetQrySeg = nullptr;
            for (auto& seqPair : qryBlk->getSequences()) {
                if (seqPair.second.getSegments().count(qSeg.start)) {
                    targetQrySeg = &(seqPair.second.getSegments()[qSeg.start]);
                    break;
                }
            }

            if (!targetRefSeg || !targetQrySeg) {
                rPos += chunkRefLen; qPos += chunkQryLen;
                continue;
            }

            // ==========================================
            // 5. 【終極修復】：反股時，CIGAR 的左右必須對調！
            // ==========================================
            int padLeftQry = aln.inverse ? rightQryOverhang : leftQryOverhang;
            int padRightQry = aln.inverse ? leftQryOverhang : rightQryOverhang;

            auto addOp = [&](mga::Cigar& c, int len, char op) {
                if (len == 0) return;
                if (!c.empty() && c.back().second == op) c.back().first += len;
                else c.push_back({len, op});
            };

            mga::Cigar paddedCigar; // 這裡裝的是純原始長度
            addOp(paddedCigar, leftRefOverhang, 'D');
            addOp(paddedCigar, padLeftQry, 'I');
            for (auto op : fragCigar) addOp(paddedCigar, op.first, op.second);
            addOp(paddedCigar, rightRefOverhang, 'D');
            addOp(paddedCigar, padRightQry, 'I');

            // ==========================================
            // 6. 【終極修復】：反股時，必須把 Qry Segment 也反轉，
            // 讓 GAP 的座標能夠與反向的 CIGAR 完美對應！
            // ==========================================
            mga::Cigar adjustedCigar = adjustCigarWithVariations(
                paddedCigar, 
                *targetRefSeg, 
                *targetQrySeg, 
                aln.inverse, 
                qryBlk->getConsensus().length()
            );

            // 交給 adjustCigar 自動膨脹，它會自動把 GAP 補在正確的 Overhang 裡！
            // mga::Cigar adjustedCigar = adjustCigarWithVariations(paddedCigar, *targetRefSeg, tempQrySeg);
            
            // ==========================================
            // === 新增的 Debug 輸出：顯示 Adjusted CIGAR ===
            // ==========================================
            if (debug) {
                std::cout << "[DEBUG] Adjusted CIGAR (Consensus space): ";
                for (auto op : adjustedCigar) std::cout << op.first << op.second;
                std::cout << "\n";
            }

            auto mergedBlock = this->mergeTwoBlocks(refBlk, qryBlk, adjustedCigar, aln.inverse);

            // ==========================================
            // === 新增的 Debug 輸出：顯示結果 ===
            // ==========================================
            if (debug) std::cout << "[DEBUG] -> Merged into New Block ID: " << mergedBlock->getId() << " (Consensus Len: " << mergedBlock->getConsensus().length() << ")\n";

            this->updateSegmentLinks(refBlk, mergedBlock);
            this->updateSegmentLinks(qryBlk, mergedBlock);
            this->deleteBlock(refBlk->getId());
            this->deleteBlock(qryBlk->getId());

            this->rebuildDictionary(dict, seqName);

            total_collapsed_bases += chunkRefLen; // 用實際配對長度統計
            merge_operations_count++;

            // 推進下一輪座標
            rPos += chunkRefLen;
            qPos += chunkQryLen;
        }
    }

    // ==========================================
    // Phase 4: 列印統計數據與最終 Block 內容
    // ==========================================
    double collapsed_ratio = (total_original_length > 0) ? 
        ((double)total_collapsed_bases / total_original_length) * 100.0 : 0.0;

    std::cout << "\n[SelfMapping] === Self-Mapping Summary ===\n";
    std::cout << "  - Total Original Length: " << total_original_length << " bp\n";
    std::cout << "  - Duplications Merged:   " << total_collapsed_bases << " bp (" 
              << std::fixed << std::setprecision(2) << collapsed_ratio << "% collapsed)\n";
    std::cout << "  - Total Merge Operations: " << merge_operations_count << "\n";
    std::cout << "==========================================\n\n";

    if (debug) {
        std::cout << "[SelfMapping] === Final Block Topology & Segments ===\n";
        for (const auto& blockPair : this->blocks_) {
            auto blk = blockPair.second;

            std::cout << "Block ID: " << blk->getId() << "\t| Consensus Len: " << blk->getConsensus().length() << "\n";

            for (auto& seqPair : blk->getSequences()) {
                std::cout << "  ├─ Sequence: " << seqPair.first << "\n";
                for (auto& segPair : seqPair.second.getSegments()) {
                    Segment& seg = segPair.second;
                    std::cout << "  │    └─ Segment Range: [" << seg.getStart() << ", " << seg.getEnd() << "]"
                              << "\t Strand: " << (seg.isReverse() ? "(-)" : "(+)");
                    if (seg.getNextBlock().lock()) std::cout << "    Next: " << seg.getNextBlock().lock()->getId();              
                    std::cout << "\n";
                }
            }
            std::cout << "  └--------------------------------------------------\n";
        }
        std::cout << "\n";
    }
}

std::pair<uint64_t, uint64_t> BlockSet::splitSingleBlock(int parentID, int localCut) {
    auto parent = this->getBlock(parentID);
    if (!parent) return {-1, -1};

    bool debug = false;

    // 1. 建立新 Block
    auto left = this->createBlock(parent->getConsensus().substr(0, localCut));
    auto right = this->createBlock(parent->getConsensus().substr(localCut));

    // 2. 切割所有 Sequence 與 Segment
    for (auto& seqPair : parent->getSequences()) {
        std::string seqID = seqPair.first;
        SequenceInfo leftSeqInfo(seqID);
        SequenceInfo rightSeqInfo(seqID);

        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& oldSeg = segPair.second;
            auto splitSegs = oldSeg.split(localCut);
            Segment& leftSeg = splitSegs.first;
            Segment& rightSeg = splitSegs.second;

            // 內部接線：嚴格遵守正反股的生物學走向
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

            if (leftSeg.getStart() != leftSeg.getEnd()) {
                leftSeqInfo.getSegments()[leftSeg.getStart()] = leftSeg;
            }
            if (rightSeg.getStart() != rightSeg.getEnd()) {
                rightSeqInfo.getSegments()[rightSeg.getStart()] = rightSeg;
            }
        }

        if (!leftSeqInfo.getSegments().empty()) left->addSequence(leftSeqInfo); 
        if (!rightSeqInfo.getSegments().empty()) right->addSequence(rightSeqInfo);
    }

    // ==========================================
    // 3. 【修復核心】：利用生物學座標進行精準的外部接線
    // ==========================================
    int rewiredCount = 0;
    
    // 掃描全圖 (包含剛建好的 left 和 right，這樣才能解開 Self-loop)
    for (auto& blockPair : blocks_) {
        auto currentBlock = blockPair.second;
        if (currentBlock == parent) continue; // 略過即將刪除的 parent

        for (auto& seqPair : currentBlock->getSequences()) {
            std::string seqID = seqPair.first;
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment& seg = segPair.second;

                // 檢查 Prev：如果這段序列是從 parent 來的
                if (seg.getPrevBlock().lock() == parent) {
                    bool found = false;
                    // 生物學上，上一段的「終點 (End)」一定等於這段的「起點 (Start)」
                    if (left->getSequences().count(seqID)) {
                        for (auto& lSeg : left->getSequences().at(seqID).getSegments()) {
                            if (lSeg.second.getEnd() == seg.getStart()) {
                                seg.setPrevBlock(left);
                                found = true; rewiredCount++; break;
                            }
                        }
                    }
                    if (!found && right->getSequences().count(seqID)) {
                        for (auto& rSeg : right->getSequences().at(seqID).getSegments()) {
                            if (rSeg.second.getEnd() == seg.getStart()) {
                                seg.setPrevBlock(right);
                                rewiredCount++; break;
                            }
                        }
                    }
                }

                // 檢查 Next：如果這段序列下一步要走到 parent
                if (seg.getNextBlock().lock() == parent) {
                    bool found = false;
                    // 生物學上，下一段的「起點 (Start)」一定等於這段的「終點 (End)」
                    if (left->getSequences().count(seqID)) {
                        for (auto& lSeg : left->getSequences().at(seqID).getSegments()) {
                            if (lSeg.second.getStart() == seg.getEnd()) {
                                seg.setNextBlock(left);
                                found = true; rewiredCount++; break;
                            }
                        }
                    }
                    if (!found && right->getSequences().count(seqID)) {
                        for (auto& rSeg : right->getSequences().at(seqID).getSegments()) {
                            if (rSeg.second.getStart() == seg.getEnd()) {
                                seg.setNextBlock(right);
                                rewiredCount++; break;
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug) std::cout << "[DEBUG-SPLIT] Block " << parentID << " cut at " << localCut 
                         << " -> L: " << left->getId() << ", R: " << right->getId() 
                         << " | Rewired pointers: " << rewiredCount << "\n";

    // 4. 安全刪除舊 Block
    this->deleteBlock(parent->getId());
    
    return {left->getId(), right->getId()};
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


/* Older Version (slower)
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
            
            // 原始座標消耗長度 (取絕對值以防反股)
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariations()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                }
            }
            
            int calculatedConsensusLen = origLen + totalGapLen;
            if (debug) {
                std::cerr << "  [REF] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]\n";
                
                if (!gapDetails.empty()) {
                    std::cout << "      -> Contains " << gapDetails.size() << " GAPs: ";
                    for (const auto& detail : gapDetails) {
                        std::cout << detail << " ";
                    }
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
            
            // 原始座標消耗長度 (取絕對值以防反股)
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariations()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                }
            }
            
            int calculatedConsensusLen = origLen + totalGapLen;
            if (debug) {
                std::cerr << "  [QRY] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]\n";

                if (!gapDetails.empty()) {
                    std::cout << "      -> Contains " << gapDetails.size() << " GAPs: ";
                    for (const auto& detail : gapDetails) {
                        std::cout << detail << " ";
                    }
                    std::cout << "\n";
                } else {
                    std::cout << "      -> No GAPs.\n";
                }
            }
        }
    }
    if (debug) std::cout << "------------------------------------------------------------\n";
    
    // 取出拷貝，避免改動原始 Block 的資料
    auto refSeqs = refBlock->getSequences();
    auto qrySeqs = qryBlock->getSequences();

    // 2. 處理反股 (Reverse Complement) - 只反轉 Consensus 字串，不動 Segment！
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

    // 3. CIGAR 長度嚴格校驗 (Heap 守護者)
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

    // 4. 走訪 CIGAR，建構新的 Consensus 與座標對應表
    std::string mergedConsensus = "";
    mergedConsensus.reserve(refSeq.length() + qrySeq.length()); 

    std::vector<int> refOldToNew(refSeq.length() + 1, 0);
    std::vector<int> qryOldToNew(qrySeq.length() + 1, 0);

    std::vector<Variation> newRefGaps;
    std::vector<Variation> newQryGaps;

    int rPos = 0, qPos = 0, mPos = 0; 

    // 【核心修正】：新增 Helper 函式，安全地從 Segment 中提取真實鹼基
    auto getBaseFromSeg = [](Segment seg, int pos, char defaultBase, bool needRc, int consLen) -> char {
        if (needRc) {
            seg.reverseComplement(consLen); // 這裡只翻轉拷貝，絕對不影響外部狀態
        }
        auto& vars = seg.getVariations();
        auto it = std::lower_bound(vars.begin(), vars.end(), pos, 
            [](Variation& v, int p) { return v.getStart() < p; });
        
        if (it != vars.end() && it->getStart() == pos && it->getType() == Variation::SNV) {
            return it->getAlt();
        }
        return defaultBase;
    };

    for (const auto& op : cigar) {
        int len = op.first;
        char type = op.second;

        if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                char rBase = refSeq.at(rPos);
                char qBase = qrySeq.at(qPos);
            
                if (rBase == qBase) {
                    mergedConsensus += rBase;
                } else {
                    // Mismatch 投票處理
                    std::map<char, int> baseFreq;
                    baseFreq[rBase] = 0; 
                    
                    for (auto& seqPair : refSeqs) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            baseFreq[getBaseFromSeg(segPair.second, rPos, rBase, false, refBlock->getConsensus().length())]++;
                        }
                    }
                    for (auto& seqPair : qrySeqs) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            // qry 端的 Segment 需要根據 inverse 決定是否即時翻轉後再取鹼基
                            baseFreq[getBaseFromSeg(segPair.second, qPos, qBase, inverse, qryBlock->getConsensus().length())]++;
                        }
                    }
                
                    char bestBase = rBase; 
                    int maxFreq = baseFreq[rBase]; 
                    for (const auto& kv : baseFreq) {
                        if (kv.second > maxFreq) {
                            maxFreq = kv.second;
                            bestBase = kv.first;
                        }
                    }
                    mergedConsensus += bestBase;
                }
                
                refOldToNew.at(rPos) = mPos;
                qryOldToNew.at(qPos) = mPos;
                rPos++; qPos++; mPos++;
            }
        }
        else if (type == 'I') { 
            mergedConsensus += qrySeq.substr(qPos, len);
            for(int i = 0; i < len; ++i) qryOldToNew.at(qPos + i) = mPos + i; 
            newRefGaps.push_back(Variation::createGap(mPos, mPos + len));
            qPos += len; mPos += len;
        } 
        else if (type == 'D') { 
            mergedConsensus += refSeq.substr(rPos, len);
            for(int i = 0; i < len; ++i) refOldToNew.at(rPos + i) = mPos + i; 
            newQryGaps.push_back(Variation::createGap(mPos, mPos + len));
            rPos += len; mPos += len;
        }
    }
    refOldToNew.at(rPos) = mPos;
    qryOldToNew.at(qPos) = mPos;

    // 5. 建立新的 Merged Block
    auto mergedBlock = this->createBlock(mergedConsensus);
    if (debug) std::cout << "[DEBUG] Merged: " << refBlock->getId() << " & " << qryBlock->getId() << " -> " << mergedBlock->getId() << '\n';

    // 6. 更新 SequenceInfo/Segment 並寫入新 Block
    auto updateAndAddSeqs = [&](std::unordered_map<std::string, SequenceInfo>& seqs, 
                                const std::vector<int>& oldToNew, 
                                const std::vector<Variation>& inducedGaps,
                                bool isQrySide, int originalConsLen) {
        
        auto& mergedSeqs = mergedBlock->getSequences();

        for (auto& seqPair : seqs) {
            std::string seqID = seqPair.first;
            
            if (mergedSeqs.find(seqID) == mergedSeqs.end()) {
                mergedSeqs[seqID] = SequenceInfo(seqID);
            }
            
            auto& targetSeq = mergedSeqs[seqID];
            
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment seg = segPair.second; // 複製原始未翻轉的 Segment
                
                // ==========================================
                // === [DEBUG-VAR] 開始追蹤 Variation 變化 ===
                // ==========================================
                if (debug) {
                    std::cout << "\n[DEBUG-VAR] Processing Sequence: " << seqID 
                              << " | Original Seg [" << seg.getStart() << " -> " << seg.getEnd() << "]"
                              << " | Source: " << (isQrySide ? "Qry Block" : "Ref Block") << "\n";
                
                    if (!seg.getVariations().empty()) {
                        std::cout << "[DEBUG-VAR]   -> Initial Variations:\n";
                        for (auto& v : seg.getVariations()) {
                            if (v.getType() == Variation::GAP) std::cout << "       - GAP: [" << v.getStart() << "->" << v.getEnd() << "]\n";
                            else std::cout << "       - SNV: " << v.getStart() << " (Alt: " << v.getAlt() << ")\n";
                        }
                    }
                }

                // 【核心修復】：在這裡才對 Qry 端的 Segment 進行座標翻轉！
                if (isQrySide && inverse) {
                    seg.reverseComplement(originalConsLen);
                    if (debug) {
                        std::cout << "[DEBUG-VAR]   -> Applied ReverseComplement (originalConsLen=" << originalConsLen << ")\n";
                        std::cout << "[DEBUG-VAR]   -> Variations after ReverseComplement:\n";
                        for (auto& v : seg.getVariations()) {
                            if (v.getType() == Variation::GAP) std::cout << "       - GAP: [" << v.getStart() << "->" << v.getEnd() << "]\n";
                            else std::cout << "       - SNV: " << v.getStart() << " (Alt: " << v.getAlt() << ")\n";
                        }
                    }
                }
                
                std::vector<Variation> newVars;
                
                if (debug) std::cout << "[DEBUG-VAR]   -> Mapping to new coordinates (oldToNew lookup):\n";
                for (auto& var : seg.getVariations()) {
                    try {
                        int newStart = oldToNew.at(var.getStart());
                        
                        if (var.getType() == Variation::SNV) {
                            newVars.push_back(Variation(newStart, var.getAlt()));
                            if (debug) std::cout << "       - SNV: " << var.getStart() << " -> " << newStart << "\n";
                        } else {
                            int newEnd = oldToNew.at(var.getEnd());
                            newVars.push_back(Variation::createGap(newStart, newEnd));
                            if (debug) std::cout << "       - GAP: [" << var.getStart() << "->" << var.getEnd() 
                                                 << "] -> [" << newStart << "->" << newEnd << "]\n";
                        }
                    } catch (const std::out_of_range& e) {
                        std::cerr << "       - [ERROR] Variation coordinate " << var.getStart() << " out of bounds!\n";
                        continue; 
                    }
                }
                
                if (!inducedGaps.empty()) {
                    if (debug) std::cout << "[DEBUG-VAR]   -> Adding " << inducedGaps.size() << " induced GAPs from CIGAR alignment.\n";
                }
                newVars.insert(newVars.end(), inducedGaps.begin(), inducedGaps.end());

                std::sort(newVars.begin(), newVars.end(), [](Variation& a, Variation& b) {
                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                    return a.getType() > b.getType(); // GAP 優先
                });
                
                std::vector<Variation> cleanedVars;
                int mergedGapCount = 0;
                for (auto& var : newVars) {
                    if (cleanedVars.empty()) {
                        cleanedVars.push_back(var);
                    } else {
                        auto& last = cleanedVars.back();
                        if (last.getType() == Variation::GAP && var.getType() == Variation::GAP && last.getEnd() >= var.getStart()) {
                            int mStart = last.getStart();
                            int mEnd = std::max(last.getEnd(), var.getEnd());
                            cleanedVars.pop_back();
                            cleanedVars.push_back(Variation::createGap(mStart, mEnd));
                            mergedGapCount++;
                        } else {
                            cleanedVars.push_back(var);
                        }
                    }
                }
                
                if (mergedGapCount > 0) {
                    if (debug) std::cout << "[DEBUG-VAR]   -> Merged " << mergedGapCount << " overlapping/adjacent GAPs.\n";
                }

                if (debug) {
                    std::cout << "[DEBUG-VAR]   -> Final Variation List for this Segment:\n";
                    for (auto& v : cleanedVars) {
                        if (v.getType() == Variation::GAP) {
                            std::cout << "       - GAP: [" << v.getStart() << "->" << v.getEnd() << "] (Len: " << (v.getEnd() - v.getStart()) << ")\n";
                        } else {
                            std::cout << "       - SNV: " << v.getStart() << " (Alt: " << v.getAlt() << ")\n";
                        }
                    }
                }

                seg.getVariations() = std::move(cleanedVars);
                targetSeq.getSegments()[seg.getStart()] = seg;
            }
        }
    };

    // 分別處理 Ref 和 Qry 的 Sequence 更新
    updateAndAddSeqs(refSeqs, refOldToNew, newRefGaps, false, refBlock->getConsensus().length());
    updateAndAddSeqs(qrySeqs, qryOldToNew, newQryGaps, true, qryBlock->getConsensus().length());

    // ==========================================
    // 7. Validation 驗證合併後的 Segment 長度
    // ==========================================
    int expectedConsensusLen = mergedBlock->getConsensus().length();
    if (debug) std::cout << "\n[DEBUG-VALIDATION] Validating Merged Block ID: " << mergedBlock->getId() 
                         << " | Expected Consensus Length: " << expectedConsensusLen << "\n";
              
    for (auto& seqPair : mergedBlock->getSequences()) {
        for (auto& segPair : seqPair.second.getSegments()) {
            Segment& seg = segPair.second;
            
            // 原始座標消耗長度 (取絕對值以防反股)
            int origLen = std::abs(seg.getEnd() - seg.getStart());
            
            int totalGapLen = 0;
            std::vector<std::string> gapDetails;
            
            for (auto& var : seg.getVariations()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    gapDetails.push_back("[" + std::to_string(var.getStart()) + "->" + std::to_string(var.getEnd()) + ", L:" + std::to_string(gapLen) + "]");
                }
            }
            
            int calculatedConsensusLen = origLen + totalGapLen;
            
            
            if (calculatedConsensusLen != expectedConsensusLen) {
                std::cerr << "  [WARNING] Sequence: " << seqPair.first 
                          << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]"
                          << " | OrigLen: " << origLen << " + Gaps: " << totalGapLen 
                          << " = " << calculatedConsensusLen 
                          << " (Mismatch with " << expectedConsensusLen << ")\n";
            } else {
                if (debug) std::cout << "  [OK] Sequence: " << seqPair.first 
                                     << " | Seg [" << seg.getStart() << ", " << seg.getEnd() << "]"
                                     << " exactly matches consensus length (" << calculatedConsensusLen << ").\n";
            }

            if (debug) {
                if (!gapDetails.empty()) {
                    std::cout << "      -> Contains " << gapDetails.size() << " GAPs: ";
                    for (const auto& detail : gapDetails) {
                        std::cout << detail << " ";
                    }
                    std::cout << "\n";
                } else {
                    std::cout << "      -> No GAPs.\n";
                }
            }
        }
    }
    if (debug) std::cout << "------------------------------------------------------------\n";

    return mergedBlock;
}
*/

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
            
            for (auto& var : seg.getVariations()) {
                if (var.getType() == Variation::GAP) {
                    int gapLen = var.getEnd() - var.getStart();
                    totalGapLen += gapLen;
                    // 【優化】：只在 debug 模式下才配發並組裝字串
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
            
            for (auto& var : seg.getVariations()) {
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
    
    // 取出拷貝，避免改動原始 Block 的資料
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

    // 4. 走訪 CIGAR，建構新的 Consensus 與座標對應表
    std::string mergedConsensus = "";
    mergedConsensus.reserve(refSeq.length() + qrySeq.length()); 

    std::vector<int> refOldToNew(refSeq.length() + 1, 0);
    std::vector<int> qryOldToNew(qrySeq.length() + 1, 0);

    std::vector<Variation> newRefGaps;
    std::vector<Variation> newQryGaps;

    int rPos = 0, qPos = 0, mPos = 0; 

    // 【核心優化】：Pass-by-reference 且使用虛擬座標映射，不配發任何記憶體！
    auto getBaseFromSeg = [](Segment& seg, int pos, char defaultBase, bool needRc, int consLen) -> char {
        int lookupPos = pos;
        if (needRc) {
            lookupPos = consLen - 1 - pos; // 虛擬翻轉，不改變陣列本身
        }
        
        auto& vars = seg.getVariations();
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
                    // 【優化】：使用原生 C-Array 取代 std::map，完全零動態配置
                    int baseFreq[256] = {0}; 
                    
                    for (auto& seqPair : refSeqs) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            char b = getBaseFromSeg(segPair.second, rPos, rBase, false, refBlock->getConsensus().length());
                            baseFreq[(unsigned char)b]++;
                        }
                    }
                    for (auto& seqPair : qrySeqs) {
                        for (auto& segPair : seqPair.second.getSegments()) {
                            char b = getBaseFromSeg(segPair.second, qPos, qBase, inverse, qryBlock->getConsensus().length());
                            baseFreq[(unsigned char)b]++;
                        }
                    }
                
                    char bestBase = rBase; 
                    int maxFreq = -1; 
                    for (int k = 0; k < 256; ++k) {
                        if (baseFreq[k] > maxFreq) {
                            maxFreq = baseFreq[k];
                            bestBase = (char)k;
                        }
                    }
                    mergedConsensus += bestBase;
                }
                
                // 【優化】：利用 [] 取代 .at() 省去邊界檢查
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
    if (debug) std::cout << "[DEBUG] Merged: " << refBlock->getId() << " & " << qryBlock->getId() << " -> " << mergedBlock->getId() << '\n';

    // 6. 更新 SequenceInfo/Segment 並寫入新 Block
    auto updateAndAddSeqs = [&](std::unordered_map<std::string, SequenceInfo>& seqs, 
                                const std::vector<int>& oldToNew, 
                                const std::vector<Variation>& inducedGaps,
                                bool isQrySide, int originalConsLen) {
        
        auto& mergedSeqs = mergedBlock->getSequences();

        for (auto& seqPair : seqs) {
            std::string seqID = seqPair.first;
            
            if (mergedSeqs.find(seqID) == mergedSeqs.end()) {
                mergedSeqs[seqID] = SequenceInfo(seqID);
            }
            
            auto& targetSeq = mergedSeqs[seqID];
            
            for (auto& segPair : seqPair.second.getSegments()) {
                Segment seg = segPair.second; 
                
                if (debug) {
                    std::cout << "\n[DEBUG-VAR] Processing Sequence: " << seqID 
                              << " | Original Seg [" << seg.getStart() << " -> " << seg.getEnd() << "]"
                              << " | Source: " << (isQrySide ? "Qry Block" : "Ref Block") << "\n";
                }

                if (isQrySide && inverse) {
                    seg.reverseComplement(originalConsLen);
                }
                
                std::vector<Variation> newVars;
                // 【優化】：提前配發足夠的容量，避免多次 reallocation
                newVars.reserve(seg.getVariations().size() + inducedGaps.size());
                
                for (auto& var : seg.getVariations()) {
                    // 利用 [] 取代 .at() 提速
                    if (var.getStart() >= oldToNew.size()) continue; // 安全防呆
                    int newStart = oldToNew[var.getStart()];
                    
                    if (var.getType() == Variation::SNV) {
                        newVars.push_back(Variation(newStart, var.getAlt()));
                    } else {
                        if (var.getEnd() >= oldToNew.size()) continue; // 安全防呆
                        int newEnd = oldToNew[var.getEnd()];
                        newVars.push_back(Variation::createGap(newStart, newEnd));
                    }
                }
                
                newVars.insert(newVars.end(), inducedGaps.begin(), inducedGaps.end());

                std::sort(newVars.begin(), newVars.end(), [](Variation& a, Variation& b) {
                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                    return a.getType() > b.getType(); // GAP 優先
                });
                
                std::vector<Variation> cleanedVars;
                cleanedVars.reserve(newVars.size());
                
                for (auto& var : newVars) {
                    if (cleanedVars.empty()) {
                        cleanedVars.push_back(var);
                    } else {
                        auto& last = cleanedVars.back();
                        if (last.getType() == Variation::GAP && var.getType() == Variation::GAP && last.getEnd() >= var.getStart()) {
                            int mStart = last.getStart();
                            int mEnd = std::max(last.getEnd(), var.getEnd());
                            cleanedVars.pop_back();
                            cleanedVars.push_back(Variation::createGap(mStart, mEnd));
                        } else {
                            cleanedVars.push_back(var);
                        }
                    }
                }

                seg.getVariations() = std::move(cleanedVars);
                targetSeq.getSegments()[seg.getStart()] = seg;
            }
        }
    };

    updateAndAddSeqs(refSeqs, refOldToNew, newRefGaps, false, refBlock->getConsensus().length());
    updateAndAddSeqs(qrySeqs, qryOldToNew, newQryGaps, true, qryBlock->getConsensus().length());

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
            
            for (auto& var : seg.getVariations()) {
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
    if (debug) std::cout << "------------------------------------------------------------\n";

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

void BlockSet::debugValidateSegments(bool verbose) {
    std::cout << "\n============================================================\n"
              << "=== BlockSet Debug Validation: " << id_ << " ===\n"
              << "============================================================\n";

    int totalBlocks = 0;
    int totalSegments = 0;
    int errorCount = 0;

    for (const auto& blockPair : blocks_) {
        std::shared_ptr<Block> blk = blockPair.second;
        int consLen = blk->getConsensus().length();
        totalBlocks++;

        // 計算這個 Block 內共有多少個 Segment
        int segmentsInBlock = 0;
        for (auto& seqPair : blk->getSequences()) {
            segmentsInBlock += seqPair.second.getSegments().size();
        }

        if (verbose) std::cout << "[Block ID: " << blk->getId() << "] "
                               << "Consensus Len: " << consLen 
                               << " | Sequences: " << blk->getSequences().size() 
                               << " | Segments: " << segmentsInBlock << "\n";

        for (auto& seqPair : blk->getSequences()) {
            std::string seqName = seqPair.first;
            auto& segments = seqPair.second.getSegments();
            
            for (auto& segPairInner : segments) {
                Segment& seg = segPairInner.second;
                totalSegments++;

                int start = seg.getStart();
                int end = seg.getEnd();
                int coordDiff = std::abs(end - start);
                
                int gapLen = 0;
                for (auto& var : seg.getVariations()) {
                    if (var.getType() == Variation::GAP) {
                        gapLen += (var.getEnd() - var.getStart());
                    }
                }

                int calculatedLen = coordDiff + gapLen;

                if (calculatedLen != consLen) {
                    std::cerr << "  ❌ [ERROR] Seq: " << seqName 
                              << " | Seg: [" << start << ", " << end << "] "
                              << (seg.isReverse() ? "(-)" : "(+)")
                              << "\n      => CoordDiff (" << coordDiff << ") + Gaps (" << gapLen 
                              << ") = " << calculatedLen << " != Consensus (" << consLen << ")\n";
                    errorCount++;
                }
            }
        }
    }
    
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Validation Complete. Checked " << totalBlocks << " blocks and " << totalSegments << " segments.\n";
    if (errorCount == 0) {
        std::cout << "🎉 PERFECT! All segment lengths match their block consensus length perfectly.\n";
    } else {
        std::cout << "🚨 CRITICAL: FOUND " << errorCount << " LENGTH MISMATCH ERROR(S)!\n";
    }
    std::cout << "============================================================\n\n";
}

void BlockSet::debugValidateLinkages(bool verbose) {
    std::cout << "\n============================================================\n"
              << "=== Topology & Linkage Debug Validation: " << id_ << " ===\n"
              << "============================================================\n";

    struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
    std::map<std::string, std::vector<SegRef>> seqTracks;
    
    // 1. 收集所有 Segment
    for (const auto& blockPair : blocks_) {
        for (auto& seqPair : blockPair.second->getSequences()) {
            for (auto& segPairInner : seqPair.second.getSegments()) {
                seqTracks[seqPair.first].push_back({ &segPairInner.second, blockPair.second });
            }
        }
    }

    int pointerErrorCount = 0;
    int coordinateGapCount = 0;

    // 2. 針對每一條 Sequence 進行排序與檢查
    for (auto& trackPair : seqTracks) {
        auto& track = trackPair.second;
        
        // 依照基因體座標嚴格排序
        std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
            return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
        });

        if (verbose) std::cout << "Checking Sequence Track: " << trackPair.first << " (" << track.size() << " segments)\n";

        for (size_t i = 0; i < track.size(); ++i) {
            auto& currRef = track[i];

            if (i > 0) {
                auto& prevRef = track[i-1];

                // --- A. 檢查座標是否連續 ---
                int prevEnd = std::max(prevRef.seg->getStart(), prevRef.seg->getEnd());
                int currStart = std::min(currRef.seg->getStart(), currRef.seg->getEnd());
                
                if (prevEnd != currStart) {
                    std::cerr << "  ⚠️ [WARNING] Coordinate Gap: Block " << prevRef.blk->getId() 
                              << " end(" << prevEnd << ") != Block " << currRef.blk->getId() 
                              << " start(" << currStart << ")\n";
                    coordinateGapCount++;
                }

                // --- B. 檢查拓撲指標 (考慮正反股的邏輯) ---
                std::shared_ptr<Block> expectedNextForPrev = currRef.blk;
                std::shared_ptr<Block> expectedPrevForCurr = prevRef.blk;

                // 根據正反股，找出實際儲存的指標方向
                std::shared_ptr<Block> actualNextForPrev = (!prevRef.seg->isReverse()) ? prevRef.seg->getNextBlock().lock() : prevRef.seg->getPrevBlock().lock();
                std::shared_ptr<Block> actualPrevForCurr = (!currRef.seg->isReverse()) ? currRef.seg->getPrevBlock().lock() : currRef.seg->getNextBlock().lock();

                if (actualNextForPrev != expectedNextForPrev) {
                    std::cerr << "  ❌ [ERROR] Broken Pointer (Forward): Block " << prevRef.blk->getId() 
                              << " does NOT point to Block " << currRef.blk->getId() << "\n";
                    pointerErrorCount++;
                }

                if (actualPrevForCurr != expectedPrevForCurr) {
                    std::cerr << "  ❌ [ERROR] Broken Pointer (Backward): Block " << currRef.blk->getId() 
                              << " does NOT point back to Block " << prevRef.blk->getId() << "\n";
                    pointerErrorCount++;
                }
            }
        }
    }

    std::cout << "------------------------------------------------------------\n";
    if (pointerErrorCount == 0 && coordinateGapCount == 0) {
        std::cout << "🎉 PERFECT! All Graph linkages and coordinates are contiguous and sound.\n";
    } else {
        std::cout << "🚨 SUMMARY: Found " << pointerErrorCount << " Pointer Error(s) and " 
                  << coordinateGapCount << " Coordinate Gap(s)!\n";
    }
    std::cout << "============================================================\n\n";
}

void BlockSet::debugValidateQuality(bool verbose) {
    std::cout << "\n============================================================\n"
              << "=== Pangenome Graph Quality Report: " << id_ << " ===\n"
              << "============================================================\n";

    if (blocks_.empty()) {
        std::cout << "  [Warning] Graph is empty. No metrics to calculate.\n";
        return;
    }

    // ---------------------------------------------------------
    // 1. 動態計算每條 Sequence 的真實基因體長度 (找出最長者)
    // ---------------------------------------------------------
    std::map<std::string, int> seqRealLengths;
    int totalSequenceCount = 0;
    
    for (auto& blkPair : blocks_) {
        for (auto& seqPair : blkPair.second->getSequences()) {
            for (auto& segPair : seqPair.second.getSegments()) {
                seqRealLengths[seqPair.first] += std::abs(segPair.second.getEnd() - segPair.second.getStart());
            }
        }
    }
    totalSequenceCount = seqRealLengths.size();
    
    int maxSeqLen = 0;
    std::string longestSeqName = "";
    for (const auto& kv : seqRealLengths) {
        if (kv.second > maxSeqLen) {
            maxSeqLen = kv.second;
            longestSeqName = kv.first;
        }
    }

    // 計算閾值
    int threshold95 = std::ceil(totalSequenceCount * 0.95);
    int threshold90 = std::ceil(totalSequenceCount * 0.90);

    // ---------------------------------------------------------
    // 2. 核心數據收集迴圈
    // ---------------------------------------------------------
    uint64_t totalConsensusLen = 0;
    uint64_t singletonLenSum = 0;
    int singletonCount = 0;
    
    // 用於計算 Global Identity (排除 Singleton)
    uint64_t globalVarLen = 0;
    uint64_t globalDenominator = 0;

    // 用於計算 N50
    std::vector<int> allBlockLengths;
    
    // 用於計算 Core vs Accessory
    int coreBlocksCount = 0;
    int softCore95BlocksCount = 0;
    int softCore90BlocksCount = 0;
    int accessoryBlocksCount = 0;
    
    uint64_t coreLenSum = 0;
    uint64_t softCore95LenSum = 0;
    uint64_t softCore90LenSum = 0;
    uint64_t accessoryLenSum = 0;

    if (verbose) std::cout << "[Block-level Identity Info]\n";

    for (const auto& blkPair : blocks_) {
        std::shared_ptr<Block> blk = blkPair.second;
        int consLen = blk->getConsensus().length();
        totalConsensusLen += consLen;
        allBlockLengths.push_back(consLen);

        int blockSegCount = 0;
        uint64_t blockVarLen = 0;
        std::set<std::string> uniqueSeqsInBlock;

        for (auto& seqPair : blk->getSequences()) {
            uniqueSeqsInBlock.insert(seqPair.first);
            for (auto& segPair : seqPair.second.getSegments()) {
                blockSegCount++;
                for (auto& var : segPair.second.getVariations()) {
                    if (var.getType() == Variation::SNV) {
                        blockVarLen += 1;
                    } else if (var.getType() == Variation::GAP) {
                        blockVarLen += (var.getEnd() - var.getStart());
                    }
                }
            }
        }

        // 統計 Core, Soft-core (95%, 90%), Accessory
        int uniqueSeqCount = uniqueSeqsInBlock.size();
        
        if (uniqueSeqCount == totalSequenceCount && totalSequenceCount > 1) {
            coreBlocksCount++;
            coreLenSum += consLen;
        } 
        
        if (uniqueSeqCount >= threshold95 && totalSequenceCount > 1) {
            softCore95BlocksCount++;
            softCore95LenSum += consLen;
        }
        
        if (uniqueSeqCount >= threshold90 && totalSequenceCount > 1) {
            softCore90BlocksCount++;
            softCore90LenSum += consLen;
        }
        
        if (uniqueSeqCount < totalSequenceCount || totalSequenceCount <= 1) {
             accessoryBlocksCount++;
             accessoryLenSum += consLen;
        }


        // 判斷 Singleton
        if (blockSegCount == 1) {
            singletonLenSum += consLen;
            singletonCount++;
            // Singleton 不印出 Identity，也不列入 Global 計算
        } else if (blockSegCount > 1 && consLen > 0) {
            // 計算並印出非 Singleton 的 Identity
            uint64_t blockDenominator = (uint64_t)blockSegCount * consLen;
            double blockIdentity = 1.0 - ((double)blockVarLen / blockDenominator);
            
            if (verbose) {
                std::cout << "  ├─ Block ID: " << std::setw(6) << std::left << blk->getId() 
                          << " | Len: " << std::setw(7) << consLen 
                          << " | Segs: " << std::setw(3) << blockSegCount
                          << " | Identity: " << std::fixed << std::setprecision(4) << blockIdentity << "\n";
            }    
            // 累積至全局計算
            globalVarLen += blockVarLen;
            globalDenominator += blockDenominator;
        }
    }

    // ---------------------------------------------------------
    // 3. 計算 N50
    // ---------------------------------------------------------
    std::sort(allBlockLengths.rbegin(), allBlockLengths.rend()); // 降序排列
    uint64_t runningSum = 0;
    int n50 = 0;
    for (int len : allBlockLengths) {
        runningSum += len;
        if (runningSum >= totalConsensusLen / 2) {
            n50 = len;
            break;
        }
    }

    // ---------------------------------------------------------
    // 4. 計算衍生指標與輸出
    // ---------------------------------------------------------
    double lenIncreaseRatio = (maxSeqLen > 0) ? ((double)totalConsensusLen / maxSeqLen - 1.0) * 100.0 : 0.0;
    double singletonRatio = (totalConsensusLen > 0) ? ((double)singletonLenSum / totalConsensusLen) * 100.0 : 0.0;
    double globalIdentity = (globalDenominator > 0) ? (1.0 - ((double)globalVarLen / globalDenominator)) : 0.0;

    std::cout << "\n------------------------------------------------------------\n";
    std::cout << ">>> GRAPH METRICS SUMMARY <<<\n\n";

    std::cout << "[1. Sequence & Graph Size]\n";
    std::cout << "  - Total Sequences        : " << totalSequenceCount << "\n";
    std::cout << "  - Total Blocks           : " << blocks_.size() << "\n";
    std::cout << "  - Longest Input Sequence : " << maxSeqLen << " bp (" << longestSeqName << ")\n";
    std::cout << "  - Total Graph Length     : " << totalConsensusLen << " bp\n";
    std::cout << "  - Graph Size Inflation   : +" << std::fixed << std::setprecision(2) << lenIncreaseRatio << " %\n";
    
    std::cout << "\n[2. Fragmentation & Contiguity]\n";
    std::cout << "  - Block N50              : " << n50 << " bp\n";
    std::cout << "  - Singleton Blocks       : " << singletonCount << " blocks\n";
    std::cout << "  - Singleton Length Ratio : " << std::fixed << std::setprecision(2) << singletonRatio << " %\n";

    std::cout << "\n[3. Evolution & Conservation]\n";
    std::cout << "  - Strict Core (100%)     : " << coreBlocksCount << " blocks (" << coreLenSum << " bp)\n";
    std::cout << "  - Soft Core (>= 95%)     : " << softCore95BlocksCount << " blocks (" << softCore95LenSum << " bp)\n";
    std::cout << "  - Soft Core (>= 90%)     : " << softCore90BlocksCount << " blocks (" << softCore90LenSum << " bp)\n";
    // std::cout << "  - Accessory (< 100%)     : " << accessoryBlocksCount << " blocks (" << accessoryLenSum << " bp)\n";

    std::cout << "\n[4. Alignment Quality]\n";
    if (globalDenominator > 0) {
        std::cout << "  - Global Average Identity: " << std::fixed << std::setprecision(4) << globalIdentity << " (Excluded singletons)\n";
    } else {
        std::cout << "  - Global Average Identity: N/A (No valid multi-segment blocks found)\n";
    }

    std::cout << "============================================================\n\n";
}


void BlockSet::refine() {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet Refine: Optimizing Graph Topology ===\n"
                              << "============================================================\n";

    // ==========================================
    // 內部神器：基於真實基因體座標的全局指標重建！
    // ==========================================
    auto rebuildAllPointers = [&]() {
        for (auto blk : this->getAllBlocks()) {
            blk->clearLinkages();
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    segPair.second.setPrevBlock(std::shared_ptr<Block>(nullptr));
                    segPair.second.setNextBlock(std::shared_ptr<Block>(nullptr));
                }
            }
        }
        struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
        std::map<std::string, std::vector<SegRef>> seqTracks;
        for (auto blk : this->getAllBlocks()) {
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    seqTracks[seqPair.first].push_back({ &segPairInner.second, blk });
                }
            }
        }
        for (auto& trackPair : seqTracks) {
            auto& track = trackPair.second;
            std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
                return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
            });
            for (size_t i = 0; i < track.size(); ++i) {
                if (i > 0) {
                    auto& prevRef = track[i-1];
                    auto& currRef = track[i];
                    if (!prevRef.seg->isReverse()) prevRef.seg->setNextBlock(currRef.blk);
                    else prevRef.seg->setPrevBlock(currRef.blk); 
                    if (!currRef.seg->isReverse()) currRef.seg->setPrevBlock(prevRef.blk);
                    else currRef.seg->setNextBlock(prevRef.blk); 
                }
            }
        }
    };

    // ==========================================
    // 階段 1: 找出長 Gap (>100) 並利用 splitSingleBlock 獨立出來
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Refine Phase 1] Isolating long gaps (>100bp)...\n";
    bool splitOccurred = true;
    while (splitOccurred) {
        splitOccurred = false;
        auto currentBlocks = this->getAllBlocks(); 
        for (auto blk : currentBlocks) {
            int cutPos = -1;
            int consLen = blk->getConsensus().length();
            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    for (auto& v : segPair.second.getVariations()) {
                        if (v.getType() == Variation::GAP && (v.getEnd() - v.getStart() > 100)) {
                            if (v.getStart() > 0 && v.getStart() < consLen) { cutPos = v.getStart(); break; }
                            if (v.getEnd() > 0 && v.getEnd() < consLen) { cutPos = v.getEnd(); break; }
                        }
                    }
                    if (cutPos != -1) break;
                }
                if (cutPos != -1) break;
            }
            if (cutPos != -1) {
                if (DEBUG_MODE) std::cout << "  -> Found long gap in Block " << blk->getId() << ", cutting at " << cutPos << "\n";
                this->splitSingleBlock(blk->getId(), cutPos);
                splitOccurred = true; 
                break; 
            }
        }
    }

    // ==========================================
    // 階段 2: 移除 Pure Gap Segments
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Refine Phase 2] Removing pure gap segments...\n";
    int removedSegs = 0;
    for (auto blk : this->getAllBlocks()) {
        std::vector<std::string> seqsToRemove;
        for (auto& seqPair : blk->getSequences()) {
            std::string seqID = seqPair.first;
            std::vector<int> segsToRemove;
            for (auto& segPairInner : seqPair.second.getSegments()) {
                if (segPairInner.second.getStart() == segPairInner.second.getEnd()) {
                    segsToRemove.push_back(segPairInner.first);
                    removedSegs++;
                }
            }
            for (int sCoord : segsToRemove) seqPair.second.getSegments().erase(sCoord);
            if (seqPair.second.getSegments().empty()) seqsToRemove.push_back(seqID);
        }
        for (const auto& s : seqsToRemove) blk->getSequences().erase(s);
    }
    if (DEBUG_MODE) std::cout << "  -> Removed " << removedSegs << " pure gap segments.\n";

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
                        
                        // 【修復核心】：針對正反股分別擴張正確的邊界
                        if (isContiguousFwd) {
                            newSeg.setEnd(rSeg.getEnd()); 
                        } else {
                            newSeg.setStart(rSeg.getStart()); 
                        }

                        for (auto v : rSeg.getVariations()) { 
                            v.shift(leftLen); 
                            newSeg.getVariations().push_back(v); 
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
                    newSeg.getVariations().push_back(Variation::createGap(leftLen, leftLen + rightLen));
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            for (size_t i = 0; i < rightSegs.size(); ++i) {
                if (!rUsed[i]) {
                    Segment newSeg = rightSegs[i];
                    newSeg.getVariations().clear(); 
                    newSeg.getVariations().push_back(Variation::createGap(0, leftLen));
                    for (auto v : rightSegs[i].getVariations()) {
                        v.shift(leftLen);
                        newSeg.getVariations().push_back(v);
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }
            if (!newSeqInfo.getSegments().empty()) newBlock->addSequence(newSeqInfo);
        }
        
        for (auto& blockPair : this->blocks_) {
            auto currentBlock = blockPair.second;
            if (currentBlock == left || currentBlock == right) continue;
            for (auto& seqPair : currentBlock->getSequences()) {
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

    // ==========================================
    // 太極迴圈：Phase 3 (Unzip) 與 Phase 4 (Concat) 互相制衡，直到圖形結晶穩定
    // ==========================================
    bool topologyChanged = true;
    while (topologyChanged) {
        topologyChanged = false;

        // ------------------------------------------
        // Phase 3: 解壓縮 (Unzip) - 門檻 <100bp
        // ------------------------------------------
        if (DEBUG_MODE) std::cout << "\n[Refine Phase 3] Unzipping highly-branched small blocks (<100bp)...\n";
        bool unzippedOccurred = true;
        while (unzippedOccurred) {
            unzippedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            
            for (auto blk : currentBlocks) {
                if (!this->getBlock(blk->getId())) continue;
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

                struct Bucket {
                    Block* pBlk;
                    Block* nBlk;
                    std::map<std::string, Segment> segs;
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
                                bucket.segs[seqID] = seg;
                                placed = true;
                                break;
                            }
                        }
                        if (!placed) {
                            Bucket newB;
                            newB.pBlk = pBlk;
                            newB.nBlk = nBlk;
                            newB.segs[seqID] = seg;
                            buckets.push_back(newB);
                        }
                    }
                }

                if (buckets.size() > 1) {
                    if (DEBUG_MODE) std::cout << "  -> Unzipping entangled Block " << blk->getId() << " into " << buckets.size() << " independent paths.\n";
                    for (auto& bucket : buckets) {
                        auto newBlock = this->createBlock(blk->getConsensus());
                        for (auto& kv : bucket.segs) {
                            SequenceInfo newSeqInfo(kv.first);
                            Segment cleanSeg = kv.second;
                            cleanSeg.setPrevBlock(std::shared_ptr<Block>(nullptr));
                            cleanSeg.setNextBlock(std::shared_ptr<Block>(nullptr));
                            newSeqInfo.getSegments()[cleanSeg.getStart()] = cleanSeg;
                            newBlock->addSequence(newSeqInfo);
                        }
                    }
                    this->deleteBlock(blk->getId());
                    unzippedOccurred = true;
                    topologyChanged = true; 
                    break;
                }
            }
            if (unzippedOccurred) rebuildAllPointers();
        }

        // ------------------------------------------
        // Phase 4: 合併 (Concat) - <100bp
        // ------------------------------------------
        if (DEBUG_MODE) std::cout << "\n[Refine Phase 4] Concatenating small blocks (<100bp)...\n";
        bool mergedOccurred = true;
        while (mergedOccurred) {
            mergedOccurred = false;
            auto currentBlocks = this->getAllBlocks();
            
            for (auto blk : currentBlocks) {
                if (!this->getBlock(blk->getId())) continue; 
                if (blk->getConsensus().length() >= 100 || blk->getSequences().empty()) continue;

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

                if (allSamePrev && sharedPrev && sharedPrev != blk) {
                    if (DEBUG_MODE) std::cout << "  -> Merging small Block " << blk->getId() << " into Prev Block " << sharedPrev->getId() << "\n";
                    concatBlocks(sharedPrev, blk);
                    mergedOccurred = true;
                    topologyChanged = true;
                    break;
                }

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

                if (allSameNext && sharedNext && sharedNext != blk) {
                    if (DEBUG_MODE) std::cout << "  -> Merging small Block " << blk->getId() << " into Next Block " << sharedNext->getId() << "\n";
                    concatBlocks(blk, sharedNext);
                    mergedOccurred = true;
                    topologyChanged = true;
                    break;
                }
            }
        }
    } // End of Outer While Loop

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    rebuildAllPointers();

    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== Refine Complete! ===\n"
                              << "============================================================\n\n";
}
/*
void BlockSet::generateRepresentativeConsensus(std::vector<std::pair<std::string, std::string>>& consensus, std::vector<std::pair<std::string, std::string>>& remaining) {

    // Use Dynamic Programming to find the optimal path of the DAG
    std::unordered_map<Block::ID, double> memo;
    std::unordered_map<Block::ID, Block::ID> path_tracker;
    
    // Objective Function
    auto calculate_weight = [](const std::shared_ptr<Block>& block) -> double {
        if (!block) return 0.0;
        double length = block->getConsensus().length();
        double count = block->getSequences().size();
        return length * count; 
    };

    std::function<double(Block::ID)> fill_dp_matrix = 
        [&](Block::ID current_id) -> double {
        
        if (memo.count(current_id)) {
            return memo[current_id];
        }

        auto current_block = getBlock(current_id);
        if (!current_block) return 0.0;

        double my_weight = calculate_weight(current_block);
        double max_prev_score = 0.0;
        Block::ID best_prev_id = -1;

        auto prev_blocks = current_block->getPrevBlocks();
        for (auto& weak_prev : prev_blocks) {
            if (auto prev = weak_prev.lock()) {
                double prev_score = fill_dp_matrix(prev->getId());
                if (prev_score > max_prev_score) {
                    max_prev_score = prev_score;
                    best_prev_id = prev->getId();
                }
            }
        }

        if (best_prev_id != -1) {
            path_tracker[current_id] = best_prev_id;
        }
        
        double total_score = my_weight + max_prev_score;
        memo[current_id] = total_score;
        return total_score;
    };

    double global_max_score = -1.0;
    Block::ID end_block_id = -1;

    for (const auto& pair : blocks_) {
        Block::ID id = pair.first;
        double score = fill_dp_matrix(id);
        if (score > global_max_score) {
            global_max_score = score;
            end_block_id = id;
        }
    }

    std::unordered_set<Block::ID> path_block_ids;
    
    consensus_path_.clear();
    Block::ID curr = end_block_id;
    while (curr != -1) {
        consensus_path_.push_back(curr);
        path_block_ids.insert(curr);
        
        if (path_tracker.find(curr) != path_tracker.end()) {
            curr = path_tracker[curr];
        } else {
            curr = -1;
        }
    }
    
    std::reverse(consensus_path_.begin(), consensus_path_.end());
    
    // Construct main consensus of the BlockSet
    std::string main_consensus_seq = "";
    for (Block::ID id : consensus_path_) {
        auto block = getBlock(id);
        if (block) {
            main_consensus_seq += block->getConsensus();
        }
    }
    
    // First element: BlockSet Consensus
    std::string main_name = id_ + "_main";
    consensus.push_back({main_name, main_consensus_seq});

    // 2. Blocks that did not be selected
    std::vector<Block::ID> all_ids;
    for(const auto& pair : blocks_) {
        all_ids.push_back(pair.first);
    }
    std::sort(all_ids.begin(), all_ids.end()); // For fix outupt

    for (Block::ID id : all_ids) {
        if (path_block_ids.find(id) == path_block_ids.end()) {
            auto block = getBlock(id);
            if (block) {
                // name: <blockSetID>_<blockID>
                std::string name = id_ + "_" + std::to_string(id);
                remaining.push_back({name, block->getConsensus()});
            }
        }
    }
}
*/