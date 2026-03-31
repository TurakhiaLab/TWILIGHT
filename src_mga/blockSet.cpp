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
    return new_block;
}

std::shared_ptr<Block> BlockSet::getBlock(Block::ID id) const {
    auto it = blocks_.find(id);
    if (it != blocks_.end()) {
        return it->second;
    }
    return nullptr;
}

std::vector<std::shared_ptr<Block>> BlockSet::getAllBlocks() const {
    std::vector<std::shared_ptr<Block>> all_blocks;
    all_blocks.reserve(blocks_.size());
    for (const auto& pair : blocks_) {
        all_blocks.push_back(pair.second);
    }
    return all_blocks;
}

std::vector<Block::ID> BlockSet::getConsensusBlocks() const {
    return consensus_path_;
}

std::vector<std::shared_ptr<Block>> BlockSet::getRemainingBlocks() const {
    std::vector<std::shared_ptr<Block>> remaining_blocks;
    std::set<Block::ID> consensus_ids;
    for (const auto& id : consensus_path_) {
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
    return blocks_.erase(id) > 0;
}


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

    std::string temp_dir = option.tempDir;
    std::string seqName = blocks_.begin()->second->getSequences().begin()->first;
    std::string seq = blocks_.begin()->second->getConsensus();
    std::string seqFile = temp_dir + "/" +  seqName + ".fa";
    std::string pafFile = temp_dir + "/" +  seqName + "-self.fa";
    mga::stringPairVec seqPair (1, {seqName, seq});
    std::cout << "[SelfMapping] 0. Writing and Self-aligning [" << seqName << "]...\n";
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
    // 假設你已經有一個函數可以讀取 PAF 並回傳 std::vector<mga::Alignment>
    std::vector<mga::Alignment> rawAlignments = mga::parser::parseMinimap2PAF(pafFile);
    std::vector<mga::Alignment> uniqueAlignments;

    // A. 處理 "兩兩一對" 的重複 (Dual mappings) 並排序
    // 簡單的作法：如果兩條 Alignment 的座標交叉重疊，只取其中一條 (例如取 idx 為偶數的)
    for (size_t i = 0; i < rawAlignments.size(); ++i) {
        if (i % 2 == 0) { // 你提到的取單數/偶數去重法
            uniqueAlignments.push_back(rawAlignments[i]);
        }
    }

    std::sort(uniqueAlignments.begin(), uniqueAlignments.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        return a.alnScore > b.alnScore; // 照 Score 由高到低處理
    });

    std::cout << "[SelfMapping] 2. Iterative Splitting & Merging...\n";

    for (const auto& aln : uniqueAlignments) {
        // 為了簡單起見，我們直接蒐集這條 Alignment 的切點
        std::set<int> cuts;
        cuts.insert(aln.refIdx.first); cuts.insert(aln.refIdx.second);
        cuts.insert(aln.qryIdx.first); cuts.insert(aln.qryIdx.second);

        // B. 切割現有 Block (橫跨多個 Block 時自動切開)
        // 這裡會沿用你寫好的 splitBlocksByCuts
        auto blocksMap = this->splitBlocksByCuts(cuts, true); 

        // C. 將 Alignment 切割成多個 Fragment (如果它橫跨了多個 Block)
        // (這裡簡化處理：透過找出 Ref 和 Qry 覆蓋的 Blocks 來逐一比對)
        int currentRefPos = aln.refIdx.first;
        int currentQryPos = aln.qryIdx.first;

        while (currentRefPos < aln.refIdx.second && currentQryPos < aln.qryIdx.second) {
            
            // 找出這段碎片對應的 Block
            auto rIt = blocksMap.upper_bound(currentRefPos); rIt--;
            auto qIt = blocksMap.upper_bound(currentQryPos); qIt--;

            int rBlkStart = rIt->first;
            int qBlkStart = qIt->first;
            std::shared_ptr<Block> refBlk = this->getBlock(rIt->second);
            std::shared_ptr<Block> qryBlk = this->getBlock(qIt->second);

            if (!refBlk || !qryBlk) break;
            if (refBlk == qryBlk) {
                // 如果已經在同一個 Block 裡（之前 Merge 過了），就跳過
                currentRefPos += refBlk->getConsensus().length();
                currentQryPos += qryBlk->getConsensus().length();
                continue;
            }

            // D. 計算這段碎片的長度，如果小於 100 則丟棄
            int fragLen = std::min({
                (int)refBlk->getConsensus().length(), 
                (int)qryBlk->getConsensus().length(), 
                aln.refIdx.second - currentRefPos
            });

            if (fragLen < 100) {
                currentRefPos += fragLen;
                currentQryPos += fragLen;
                continue; 
            }

            // E. 尋找對應的 SequenceInfo (這裡利用 original coordinate 來比對)
            SequenceInfo* targetRefSeq = nullptr;
            for (auto& seq : refBlk->getSequences()) {
                if (seq.start_coordinate <= currentRefPos && seq.end_coordinate >= currentRefPos) {
                    targetRefSeq = &seq; break;
                }
            }
            SequenceInfo* targetQrySeq = nullptr;
            for (auto& seq : qryBlk->getSequences()) {
                if (seq.start_coordinate <= currentQryPos && seq.end_coordinate >= currentQryPos) {
                    targetQrySeq = &seq; break;
                }
            }

            if (targetRefSeq && targetQrySeq) {
                // F. 從原本的巨型 CIGAR 中截取這段 Fragment 的 CIGAR
                mga::Cigar fragOrigCigar = extractSubCigar(aln.CIGAR, currentRefPos - aln.refIdx.first, fragLen);

                // G. 核心魔法：根據之前留下的 Variations 動態調整 CIGAR！
                mga::Cigar adjustedCigar = adjustCigarWithVariations(fragOrigCigar, *targetRefSeq, *targetQrySeq);

                // H. 物理合併！
                auto mergedBlock = this->mergeTwoBlocks(refBlk, qryBlk, adjustedCigar, aln.inverse);
                
                // (因為都在同一個 Set 裡，這裡記得把舊的 Block 從你的內部資料結構中刪除，替換成 mergedBlock)
                this->deleteBlock(refBlk->getId());
                this->deleteBlock(qryBlk->getId());
                this->addBlock(mergedBlock);
                
                // 重新接線 (Rewiring) - 沿用你原本 splitSingleBlock 裡的接線邏輯
                // ...
            }

            currentRefPos += fragLen;
            currentQryPos += fragLen;
        }
    }
    std::cout << "[SelfMapping] Done. All duplications collapsed!\n";
}



/*
std::shared_ptr<Block> BlockSet::mergeTwoBlocks(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, const mga::Cigar& cigar, bool inverse) 
{
    // 1. 取得舊 Consensus 與 SequenceInfo
    std::string refSeq = refBlock->getConsensus();
    std::string qrySeq = qryBlock->getConsensus();
    
    std::vector<SequenceInfo> refSeqs = refBlock->getSequences();
    std::vector<SequenceInfo> qrySeqs = qryBlock->getSequences();

    // 2. 處理反股 (Reverse Complement)
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
        for (auto& seq : qrySeqs) {
            seq.reverseComplement(qrySeq.length());
        }
    }

    // 3. 走訪 CIGAR，建構新的 Consensus 與座標對應表
    std::string mergedConsensus = "";
    mergedConsensus.reserve(refSeq.length() + qrySeq.length()); 

    std::vector<int> refOldToNew(refSeq.length() + 1, 0);
    std::vector<int> qryOldToNew(qrySeq.length() + 1, 0);

    std::vector<Variation> newRefGaps;
    std::vector<Variation> newQryGaps;

    int rPos = 0, qPos = 0, mPos = 0; 

    // 【優化】：使用 Binary Search (O(log V)) 快速尋找 SNV
    auto getSnvBase = [](const std::vector<Variation>& vars, int pos, char defaultBase) -> char {
        auto it = std::lower_bound(vars.begin(), vars.end(), pos, 
            [](const Variation& v, int p) { return v.start < p; });
        if (it != vars.end() && it->start == pos && it->type == Variation::SNV) {
            return it->alt;
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
                } 
                else {
                    // Mismatch: 重新計算 Consensus，使用二元搜尋優化效能
                    std::map<char, int> baseFreq;
                    baseFreq[rBase] = 0; // 確保 rBase 在 Map 中
                    
                    for (const auto& seq : refSeqs) {
                        baseFreq[getSnvBase(seq.variations, rPos, rBase)]++;
                    }
                    for (const auto& seq : qrySeqs) {
                        baseFreq[getSnvBase(seq.variations, qPos, qBase)]++;
                    }
                
                    // 【優化】：找出數量最多的鹼基，平手時優先保留 rBase
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
    // 補上最後的尾部邊界映射
    refOldToNew[rPos] = mPos;
    qryOldToNew[qPos] = mPos;

    // 4. 建立新的 Merged Block
    // (這裡假設 BlockSet 類別有一個 createBlock 方法，或者你可以換成 std::make_shared<Block>(...))
    auto mergedBlock = this->createBlock(mergedConsensus);

    // 5. Lambda: 更新 SequenceInfo 並寫入新 Block
    auto updateAndAddSeqs = [&](std::vector<SequenceInfo>& seqs, const std::vector<int>& oldToNew, const std::vector<Variation>& inducedGaps) {
        for (auto& seq : seqs) {
            
            // A. 將原有的 Variations 轉換到新座標
            for (auto& var : seq.variations) {
                var.start = oldToNew[var.start];
                if (var.type == Variation::SNV) {
                    // 【修復 Bug 1】：SNV 的長度必須永遠維持 1，不能被 Insertion 拉長！
                    var.end = var.start + 1;
                } else {
                    // GAP 則正常使用 mapped 的 end，讓它自然擴張
                    var.end = oldToNew[var.end];
                }
            }
            
            // B. 加入因為比對產生的新 Gaps
            seq.variations.insert(seq.variations.end(), inducedGaps.begin(), inducedGaps.end());

            // C. 排序
            std::sort(seq.variations.begin(), seq.variations.end(), [](const Variation& a, const Variation& b) {
                if (a.start != b.start) return a.start < b.start;
                return a.type > b.type; // 若起點相同，GAP 排在 SNV 前面
            });
            
            // D. 合併相連/重疊的 Gaps
            std::vector<Variation> cleanedVars;
            for (const auto& var : seq.variations) {
                if (cleanedVars.empty()) {
                    cleanedVars.push_back(var);
                } else {
                    auto& last = cleanedVars.back();
                    // 【修復 Bug 2】：只有 GAP 可以被合併！相鄰的 SNV 絕對不能合併。
                    if (last.type == Variation::GAP && var.type == Variation::GAP && last.end >= var.start) {
                        last.end = std::max(last.end, var.end);
                    } else {
                        cleanedVars.push_back(var);
                    }
                }
            }
            seq.variations = std::move(cleanedVars);

            // E. 把更新完畢的 sequence 加進新 Block
            mergedBlock->addSequence(seq);
        }
    };

    // 6. 套用更新並回傳
    updateAndAddSeqs(refSeqs, refOldToNew, newRefGaps);
    updateAndAddSeqs(qrySeqs, qryOldToNew, newQryGaps);

    return mergedBlock;
}
*/