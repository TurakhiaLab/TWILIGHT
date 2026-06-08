#include "block.hpp"
#include "mga.hpp"
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <functional>
#include <unordered_set>
#include <chrono>

// Helper to truncate strings for printing
std::string truncate(const std::string& str, size_t width) {
    if (str.length() > width) {
        return str.substr(0, width - 3) + "...";
    }
    return str;
}



/**
 * 根據 minimap2 的 alignments，判斷每個 Query 基因組是否需要進行 Reverse Complement。
 * * @param alignments 來自 minimap2 解析後的 alignment 陣列
 * @return std::map<std::string, bool> 回傳每個 qryName 是否需要 RC (true 代表需要 RC)
 */
std::map<std::string, bool> decideReverseComplement(const std::vector<Alignment>& alignments) {
    
    
}

void BlockManager::orientCircularGenomes(Option& option, std::string refSequenceName) {
    bool debug = true;
    bool write_rotate = option.writeOriented;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto log_step_time = [&](const std::string& step_name) {
        if (!debug) return;
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = current_time - start_time;
        std::cout << "[DEBUG] " << std::left << std::setw(30) << step_name 
                  << " | Elapsed: " << elapsed.count() << " ms\n";
    };

    if (debug) std::cout << "\n[DEBUG] --- Start orientCircularGenomes ---\n";

    if (blocksets.empty() || sequence_lengths.empty()) {
        std::cerr << "Warning: No blocksets or sequences to orient.\n";
        return;
    }

    // 1. 決定 Reference Sequence Name
    std::string bestRefName = refSequenceName;
    if (bestRefName.empty()) { 
        int maxLen = -1;
        for (const auto& pair : sequence_lengths) {
            if (pair.second > maxLen) {
                maxLen = pair.second;
                bestRefName = pair.first;
            }
        }
    }

    if (blocksets.find(bestRefName) == blocksets.end()) {
        std::cerr << "Error: Chosen reference " << bestRefName << " not found in blocksets.\n";
        return;
    }

    if (debug) std::cout << "[DEBUG] Selected Reference: " << bestRefName 
                         << " (Length: " << sequence_lengths[bestRefName] << ")\n";
    log_step_time("Reference Selected");

    // 2. 擷取 Reference 和 Query 的序列
    BlockSet* refBlockSet = blocksets[bestRefName].get();
    auto reference_genome = refBlockSet->getRepresentativeConsensus();
    
    StringPairs query_genomes(this->blocksets.size() - 1);
    int qryCount = 0;
    for (const auto& pair : this->blocksets) {
        if (pair.first != bestRefName) {
            BlockSet* qryBlockSet = pair.second.get();
            query_genomes[qryCount] = qryBlockSet->getRepresentativeConsensus()[0];
            qryCount++;
        }
    }
    log_step_time("Prepared Queries for minimap2");

    // 3. 執行 Minimap2
    if (debug) std::cout << "[DEBUG] Running minimap2 on " << qryCount << " queries...\n";
    auto allAlignments = runMinimap2(reference_genome, query_genomes, "reference_genome", "query_genomes", option);
    log_step_time("Minimap2 Completed");

    // 4. 統計方向 & 尋找最佳定錨點 (Hybrid Approach)
    struct StrandStats {
        long long forwardLength = 0;
        long long reverseLength = 0;
    };
    std::map<std::string, StrandStats> qryStatsMap;
    std::map<std::string, int> maxAlnLenMap; // 記錄每個 query 的最長 alignment 長度
    
    for (const auto& aln : allAlignments) {
        if (!aln.valid) continue; 
        
        if (aln.inverse) {
            qryStatsMap[aln.qryName].reverseLength += aln.alnLength;
        } else {
            qryStatsMap[aln.qryName].forwardLength += aln.alnLength;
        }

        // 更新該 Query 的最大長度
        maxAlnLenMap[aln.qryName] = std::max(maxAlnLenMap[aln.qryName], aln.alnLength);
    }
    
    std::map<std::string, bool> rcDecisionMap;
    for (const auto& pair : qryStatsMap) {
        rcDecisionMap[pair.first] = (pair.second.reverseLength > pair.second.forwardLength);
    }

    // 尋找「夠長」且「最接近 Ref 起點」的 Alignment 當作 Shift 基準
    std::map<std::string, const Alignment*> bestAlnMap; 
    std::map<std::string, int> minRefStartMap; // 記錄最接近 0 的距離
    
    for (const auto& aln : allAlignments) {
        if (!aln.valid) continue;
        bool needsRC = rcDecisionMap[aln.qryName];
        if (aln.inverse != needsRC) continue; // 方向錯誤的雜訊不要理

        // 門檻：長度必須大於該 Query 最長 Alignment 的 50% (可調參數)
        int lengthThreshold = 10000;
        if (aln.alnLength < lengthThreshold) continue;

        // 在符合長度條件的候選者中，找 ref.start 最接近 0 的
        auto it = minRefStartMap.find(aln.qryName);
        if (it == minRefStartMap.end() || aln.refIdx.first < it->second) {
            minRefStartMap[aln.qryName] = aln.refIdx.first;
            bestAlnMap[aln.qryName] = &aln;
        }
    }
    log_step_time("Strand Voting & Hybrid Anchor Found");
    
    // 5. 執行 Reverse Complement 與 Shift (Rotation)
    int rcCount = 0;
    int shiftCount = 0;

    for (const auto& pair : this->blocksets) {
        std::string qryName = pair.first;
        if (qryName == bestRefName) continue; // 略過 Reference 本身

        bool needsRC = rcDecisionMap[qryName];
        std::string seq = this->sequences[qryName];
        int L = sequence_lengths[qryName];
        int raw_shift = 0;

        // 計算偏移量 (如果有找到對應的主 Alignment)
        if (bestAlnMap.count(qryName)) {
            const Alignment* bestAln = bestAlnMap[qryName];
            if (needsRC) {
                // RC 狀態下的位移數學轉換
                raw_shift = (L - bestAln->qryIdx.second) - bestAln->refIdx.first;
            } else {
                raw_shift = bestAln->qryIdx.first - bestAln->refIdx.first;
            }
        }

        // 把 shift 處理為標準的正數 (環狀模除)
        int true_shift = ((raw_shift % L) + L) % L;

        // --- 應用轉換 ---
        if (needsRC) {
            seq = getReverseComplement(seq);
            rcCount++;
        }
        
        if (true_shift != 0) {
            // 字串環狀位移 (將前面 true_shift 長度搬到最後面)
            seq = seq.substr(true_shift) + seq.substr(0, true_shift);
            shiftCount++;
        }

        if (debug) {
            std::cout << "[DEBUG] " << std::left << std::setw(15) << qryName 
                      << " | RC: " << (needsRC ? "YES" : "NO ") 
                      << " (Forward: " << qryStatsMap[qryName].forwardLength << " bp)"
                      << " vs (Reverse: " << qryStatsMap[qryName].reverseLength << " bp"
                      << " | Shifted: " << true_shift << " bp\n";
        }

        // 寫回原本的資料結構
        BlockSet* qryBlockSet = pair.second.get();
        auto first_block = qryBlockSet->getAllBlocks()[0];
        if (auto sp = first_block.lock()) {
            sp->setConsensus(seq, true); 
        }
        this->sequences[qryName] = seq;
    }

    if (debug) std::cout << "[DEBUG] Inverted: " << rcCount << " | Shifted: " << shiftCount << " genomes.\n";
    log_step_time("Transformations Applied");

    // ==========================================================
    // 🌟 6. 輸出轉換後的完整序列 (Reference + Oriented Queries)
    // ==========================================================

    StringPairs output_seqs;
    
    // 6-1. 將 Reference 放在第一筆，維持對齊基準的直覺性
    output_seqs.push_back({bestRefName, this->sequences[bestRefName]});
    
    // 6-2. 收集所有其他的 Query 序列
    for (const auto& pair : this->sequences) {
        if (pair.first != bestRefName) {
            output_seqs.push_back({pair.first, pair.second});
        }
    }

    // 6-3. 呼叫 IO 寫出檔案 (檔名可以依需求改為 option 裡的變數)
    std::string outFileName = option.tempDir + "/oriented_genomes.fasta"; 
    bool isCompressed = false; // 依需求調整
    bool appendMode = false;   // 覆寫模式
    
    if (write_rotate) mga::io::writeAlignment(outFileName, output_seqs, isCompressed, appendMode);
    
    if (debug && write_rotate) std::cout << "[DEBUG] Wrote " << output_seqs.size() << " oriented sequences to '" << outFileName << "'\n";
    if (debug && write_rotate) log_step_time("Output Sequences Written");
    // ==========================================================

    if (debug) std::cout << "[DEBUG] --- End orientCircularGenomes ---\n\n";
    
    /*
    // ==========================================================
    // 🌟 6. 輸出轉換後的完整序列 (Reference + Oriented Queries)
    // ==========================================================
    const size_t SAMPLE_SIZE = 50; // 🌟 設定最大輸出數量 (包含 Reference)
    StringPairs output_seqs;
    
    // 6-1. 將 Reference 放在第一筆，維持對齊基準的直覺性
    output_seqs.push_back({bestRefName, this->sequences[bestRefName]});
    
    // 6-2. 收集所有其他的 Query 序列 (套用 Quality Pass 與 Sample Size 限制)
    int passed_quality_count = 0;
    for (const auto& pair : this->sequences) {
        if (pair.first != bestRefName) {
            std::string qryName = pair.first;
            
            // 取得該序列的正反向比對長度
            long long f_len = qryStatsMap[qryName].forwardLength;
            long long r_len = qryStatsMap[qryName].reverseLength;
            
            // 🛡️ Quality Pass 條件：長度必須壓倒性地偏向某一邊 (> 3倍)
            bool pass_quality = (f_len > 3 * r_len) || (r_len > 3 * f_len);

            if (pass_quality) {
                passed_quality_count++;
                // 檢查是否還沒達到 Sample 數量上限
                if (output_seqs.size() < SAMPLE_SIZE) {
                    output_seqs.push_back({qryName, pair.second});
                }
            } else {
                // 沒有通過 Quality Pass (可能是 50-50 倒位，或是完全沒有 alignment)
                if (debug && write_rotate) {
                    std::cout << "[DEBUG] 🛑 [FILTERED] " << std::left << std::setw(15) << qryName 
                              << " excluded from output (F: " << f_len << " vs R: " << r_len << "). Strand ambiguous.\n";
                }
            }
        }
    }

    // 6-3. 呼叫 IO 寫出檔案
    std::string outFileName = option.tempDir + "/oriented_genomes.fasta"; 
    bool isCompressed = false; 
    bool appendMode = false;   
    
    if (write_rotate) mga::io::writeAlignment(outFileName, output_seqs, isCompressed, appendMode);
    
    if (debug && write_rotate) {
        std::cout << "[DEBUG] Total queries passed quality : " << passed_quality_count << "\n";
        std::cout << "[DEBUG] Wrote " << output_seqs.size() << " oriented sequences (capped at " << SAMPLE_SIZE << ") to '" << outFileName << "'\n";
        log_step_time("Output Sequences Written");
    }
    */
    // ==========================================================
}


// =======================
// BlockManager Implementation
// =======================

/*
void BlockManager::updateLongestSequences() {
    for (const auto& pair : blocksets) {
        if (pair.second) {
            pair.second->updateLongestSequence(this->sequence_lengths);
        }
    }
}

void BlockManager::print(std::ostream& os) const {
    os << "\n";
    os << "############################################################\n";
    os << "               BLOCK MANAGER SYSTEM STATUS                  \n";
    os << "############################################################\n\n";
    os << "> Total BlockSets Active: " << blocksets.size() << "\n";

    for (const auto& pair : blocksets) {
        if (pair.second) {
            pair.second->print(os);
        }
    }
    os << "############################################################\n";
}
*/

/*
void mergeAdjacentBlocks(std::list<std::shared_ptr<Block>>& list, const std::set<int>& cuts) {
    std::list<std::shared_ptr<Block>> stitchedList;
    std::vector<std::shared_ptr<Block>> buffer;
        
    int currentGlobalPos = 0;

    for (auto& block : list) {
        int len = block->getConsensus().length();
        
        // 先將當前 Block 加入緩衝區
        buffer.push_back(block);
        currentGlobalPos += len;
        // 檢查當前位置是否是一個「切點」
        // 如果是切點，或是列表最後一個，就執行縫合
        if (cuts.count(currentGlobalPos) || (!cuts.empty() && currentGlobalPos == *cuts.rbegin())) {
            if (buffer.size() == 1) {
                stitchedList.push_back(buffer[0]);
            } else {
                stitchedList.push_back(stitchBlocks(buffer));
            }
            buffer.clear();
        }
    }
    
    if (!buffer.empty()) {
        stitchedList.push_back(stitchBlocks(buffer));
    }
    list = std::move(stitchedList);
}

// --- [NEW] Helper: Stitch multiple blocks into one ---
std::shared_ptr<Block> stitchBlocks(const std::vector<std::shared_ptr<Block>>& blocks) {
    if (blocks.empty()) return nullptr;
    std::cerr << "Merge " << blocks.size() << " blocks\n";
    // 1. 合併 Consensus
    std::string newConsensus = "";
    std::vector<int> offsets; // 記錄每一塊的起始偏移量
    int totalLen = 0;
    for (const auto& b : blocks) {
        offsets.push_back(totalLen);
        newConsensus += b->getConsensus();
        totalLen += b->getConsensus().length();
    }

    auto newID = blocks[0]->getId();

    // 產生新 ID (使用 internal counter)
    auto newBlock = std::make_shared<Block>(newID, newConsensus);
        
    // 2. 收集所有涉及的 Sequence IDs (Union)
    std::map<std::string, std::vector<const SequenceInfo*>> seqMap;
    for (size_t i = 0; i < blocks.size(); ++i) {
        for (const auto& seq : blocks[i]->getSequences()) {
            // 記錄 seqID -> 在第 i 塊的 SequenceInfo 指標
            // 我們需要知道它在哪一塊出現，哪一塊沒出現
            // 為了處理 "沒出現"，我們可以用一個固定大小的 vector
            if (seqMap.find(seq.sequence_id) == seqMap.end()) {
                seqMap[seq.sequence_id].resize(blocks.size(), nullptr);
            }
            seqMap[seq.sequence_id][i] = &seq;
        }
    }

    // 3. 建構合併後的 Sequences
    for (const auto& pair : seqMap) {
        const std::string& seqID = pair.first;
        const auto& fragments = pair.second; // vector of SequenceInfo* (可以是 null)
        SequenceInfo newSeq;
        newSeq.sequence_id = seqID;
        
        // 尋找這一長串中的 Start 和 End Coordinate
        // Start: 第一個非 null fragment 的 start
        // End: 最後一個非 null fragment 的 end
        bool foundStart = false;
        int lastValidEnd = -1;
        for (size_t i = 0; i < blocks.size(); ++i) {
            int offset = offsets[i];
            int blockLen = blocks[i]->getConsensus().length();
            if (fragments[i] != nullptr) {
                // --- Case A: Sequence exists in this block ---
                const auto* sub = fragments[i];
                // Coordinate Logic
                if (!foundStart) {
                    newSeq.start_coordinate = sub->start_coordinate;
                    foundStart = true;
                }
                lastValidEnd = sub->end_coordinate;
                // Copy Variations (並加上 Offset)
                for (auto var : sub->variations) {
                    var.shift(offset); // start += offset, end += offset
                    newSeq.variations.push_back(var);
                }
            } else {
                // --- Case B: Sequence missing in this block (GAP) ---
                // 在 Consensus [offset, offset + blockLen) 這段區間是 Gap
                
                // Variation: Create a GAP covering this whole block
                newSeq.variations.push_back(Variation::createGap(offset, offset + blockLen));
                // Raw Sequence: 不增加長度 (Gap 不佔用 raw bases)
                
                // Coordinate: 不更新 start/end，保持原樣
            }
        }
        newSeq.end_coordinate = lastValidEnd;
        newBlock->addSequence(newSeq);
    }

    return newBlock;
}


/*
void BlockManager::applySplits(std::list<std::shared_ptr<Block>>& list, const std::set<int>& cuts) {
    size_t currentGlobalPos = 0; 
    int max_id = -1;
    
    for (auto it = list.begin(); it != list.end(); ++it) {
        if ((*it)->getId() > max_id) {
            max_id = (*it)->getId();
        }
    }
    
    auto it = list.begin();
    while (it != list.end()) {
        auto block = *it;
        size_t blkLen = block->getConsensus().length();
        size_t blkStart = currentGlobalPos;
        size_t blkEnd = currentGlobalPos + blkLen;
        
        // Use upper_bound to find the first cut that is larger than start
        auto cutIt = cuts.upper_bound(blkStart);
        
        std::vector<int> internalCuts;
        while (cutIt != cuts.end() && *cutIt < blkEnd) {
            internalCuts.push_back(*cutIt - blkStart); 
            cutIt++;
        }
        
        if (!internalCuts.empty()) {
            std::reverse(internalCuts.begin(), internalCuts.end());
            
            std::shared_ptr<Block> curr = block;
            std::vector<std::shared_ptr<Block>> fragments;
            
            for (int offset : internalCuts) {
                auto id1 = ++max_id; 
                auto id2 = ++max_id; 
                
                auto res = curr->split(offset, id1, id2);
                if (res.first && res.second) {
                    fragments.push_back(res.second);
                    curr = res.first;
                } else {
                    std::cerr << "[Error] Split failed for block ID: " 
                              << curr->getId() << " at offset: " << offset << '\n';
                }
            }
            fragments.push_back(curr); 
            
            it = list.erase(it); 
            
            for (auto r_it = fragments.rbegin(); r_it != fragments.rend(); ++r_it) {
                list.insert(it, *r_it); 
            }
            
            currentGlobalPos += blkLen; 
        } else {
            currentGlobalPos += blkLen;
            ++it;
        }
    }
}

std::list<std::shared_ptr<Block>> BlockManager::cloneConsensusPath(BlockSet* set) {
    std::list<std::shared_ptr<Block>> list;
    for (Block::ID id : set->getConsensusBlocks()) {
        auto original = set->getBlock(id);
        if (original) list.push_back(original->clone(id));
    }
    return list;
} 

std::list<std::shared_ptr<Block>> BlockManager::prepareBlocks(BlockSet* set, const std::vector<mga::Alignment>& alns, bool isRef) {
    auto list = cloneConsensusPath(set);
    std::set<int> cuts;
    
    for(const auto& aln : alns) {
        if(!aln.valid) continue;
        // Collect cuts based on type
        if (isRef) {
            if (aln.refIdx.first != -1) {
                cuts.insert(aln.refIdx.first);
                cuts.insert(aln.refIdx.second);
            }
        } else {
            if (aln.qryIdx.first != -1) {
                cuts.insert(std::min(aln.qryIdx.first, aln.qryIdx.second));
                cuts.insert(std::max(aln.qryIdx.first, aln.qryIdx.second));
            }
        }
    }
    
    applySplits(list, cuts);
    mergeAdjacentBlocks(list, cuts); // Stitching
    return list;
}


BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, const std::vector<mga::Alignment>& alignments) {
        
    std::cout << "[Merger] 1. Preparing Blocks (Clone -> Split -> Stitch)...\n";

    auto refWorkList = prepareBlocks(refSet, alignments, true); // true = collect ref cuts
    auto qryWorkList = prepareBlocks(qrySet, alignments, false); // false = collect qry cuts
        
    // 建立索引 (GlobalPos -> BlockPtr) 以便快速查找
    auto refIndex = indexBlocks(refWorkList);
    auto qryIndex = indexBlocks(qryWorkList);

    // 2. 初始化結果 BlockSet
    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* resultSet = createBlockSet(newID);

    // Map: Old Block ID (from WorkList) -> New Merged Block ID
    // 用來追蹤每個碎塊最後跑去哪裡了
    std::map<Block::ID, std::shared_ptr<Block>> refBlockMapping; // Old Ref ID -> New Block
    std::map<Block::ID, std::shared_ptr<Block>> qryBlockMapping; // Old Qry ID -> New Block
    std::set<Block::ID> consumedRefBlocks;
    std::set<Block::ID> consumedQryBlocks;

    std::cout << "[Merger] 2. Processing Case 1: Primary Merges...\n";

        // =========================================================
        // Phase A: Handle Primary Alignments (Case 1)
        // =========================================================
        for (const auto& aln : alignments) {
            if (!aln.valid || aln.type != mga::PRIMARY) continue;

            mga::Alignment currAln = aln;

            while (currAln.valid && !currAln.CIGAR.empty()) {
                
                auto rBlocks = getBlocksInRange(refIndex, {currAln.refIdx.first, currAln.refIdx.first + 1});
                if (rBlocks.empty()) break; 
                auto rBlk = rBlocks[0];

                // 2. 找出這個 Ref Block 的絕對座標邊界
                // (你需要從 index 中找回它的座標，假設有 helper 或是你能從 map 反查)
                int rStart = getBlockGlobalStart(refIndex, rBlk); 
                int rEnd = rStart + rBlk->getConsensus().length();

                // 3. 切割 Alignment！取得剛好對應這個 Ref Block 的 sub-alignment
                mga::SplitResult parts = mga::splitAlignment(currAln, rStart, rEnd, true );
                
                if (parts.head.valid) std::cerr << "ERROR: Head Valid.\n";
                

                mga::Alignment subAln = parts.middle;
                currAln = parts.tail; // 剩下的留到下一個 while 迴圈處理

                if (!subAln.valid) continue;
                
                std::cerr << subAln.refIdx.first << " " << subAln.refIdx.second << '\n';
                std::cerr << subAln.qryIdx.first << " " << subAln.qryIdx.second << '\n';


                // 4. 找出對應的 Qry Block
                int qMin = std::min(subAln.qryIdx.first, subAln.qryIdx.second);
                int qMax = std::max(subAln.qryIdx.first, subAln.qryIdx.second);
                
                auto qBlocks = getBlocksInRange(qryIndex, {qMin, qMax});
                if (qBlocks.empty()) continue;
                
                auto qBlk = qBlocks[0]; // 假設已做好預切割與 Stitching，這裡應是一對一

                // 5. 呼叫新的核心函數：執行真正的 CIGAR Projection 與 Sequence Merge
                auto newBlock = mergeBlockPair(resultSet, rBlk, qBlk, subAln);

                std::cerr << "Finish Merging\n";
                // 6. 記錄映射關係
                refBlockMapping[rBlk->getId()] = newBlock;
                consumedRefBlocks.insert(rBlk->getId());

                for (auto qb : qBlocks) { 
                    qryBlockMapping[qb->getId()] = newBlock;
                    consumedQryBlocks.insert(qb->getId());
                }
            }
        }
        // =========================================================
        // Phase B: Handle Unaligned Blocks (Remainders)
        // =========================================================
        auto addRemainder = [&](const std::list<std::shared_ptr<Block>>& list, 
                                std::set<Block::ID>& consumed,
                                std::map<Block::ID, std::shared_ptr<Block>>& mapping) {
            for (const auto& blk : list) {
                if (consumed.find(blk->getId()) == consumed.end()) {
                    auto newBlock = resultSet->createBlock(blk->getConsensus());
                    newBlock->setFromPrimary(false); 
                    for (const auto& seq : blk->getSequences()) newBlock->addSequence(seq);

                    mapping[blk->getId()] = newBlock; // 存入對應的 Map
                    consumed.insert(blk->getId());
                }
            }
        };

        addRemainder(refWorkList, consumedRefBlocks, refBlockMapping);
        addRemainder(qryWorkList, consumedQryBlocks, qryBlockMapping);

        std::cout << "[Merger] 3. Processing Case 2 & 3: Secondary Alignments...\n";

        // =========================================================
        // Phase C: Handle Secondaries (Case 2 & Case 3)
        // =========================================================
        for (const auto& aln : alignments) {
            if (!aln.valid || aln.type != mga::SECONDARY) continue;

            // 1. 找出這條 Secondary Alignment 連接了哪兩個 "New Block"
            // 注意：我們找的是 Ref 區間對應的 Block 和 Qry 區間對應的 Block
            // 因為所有的 Block (Primary or Unaligned) 都已經建立好了，這裡只是查表
            
            auto rBlocks = getBlocksInRange(refIndex, aln.refIdx);
            auto qBlocks = getBlocksInRange(qryIndex, aln.qryIdx);

            if (rBlocks.empty() || qBlocks.empty()) continue;

            // 針對每一對涉及的 Block 處理
            // 通常 Secondary Alignment 涵蓋的區間應該只對應到一個或少數幾個 Block
            for (auto rBlkOriginal : rBlocks) {
                for (auto qBlkOriginal : qBlocks) {
                    
                    auto targetBlock = refBlockMapping[rBlkOriginal->getId()]; // Ref side block
                    auto sourceBlock = qryBlockMapping[qBlkOriginal->getId()]; // Qry side block

                    if (!targetBlock || !sourceBlock || targetBlock == sourceBlock) continue;

                    // 關鍵判斷邏輯
                    bool targetIsPrimary = targetBlock->isFromPrimary();
                    bool sourceIsPrimary = sourceBlock->isFromPrimary();

                    // --- Case 3: Paralog (Both are Primary) ---
                    // "同 Case 2 但那對 duplication 在另一個方向又有 match 到"
                    // 意思是：Source Block 其實已經跟別的 Ref Block 合併了 (所以它是 Primary)
                    if (targetIsPrimary && sourceIsPrimary) {
                        // 標記 Paralog，然後丟掉 Alignment 片段 (不合併 sequence)
                        targetBlock->addParalog(sourceBlock->getId());
                        sourceBlock->addParalog(targetBlock->getId());
                        // Done. Discard duplication segment.
                    }
                    
                    // --- Case 2: Duplication (One is Primary, One is Unaligned) ---
                    // "只 merge 在 primary 上的那對，然後把 duplication 保留"
                    // 意思是：Source Block 沒有參與 Primary Merge (它是孤兒/Unaligned)，
                    // 但它透過 Secondary Alignment 連到了 Target。
                    else {
                        // 標記 Duplication
                        // 通常是 Primary 指向 Secondary (Copy)
                        if (targetIsPrimary && !sourceIsPrimary) {
                            targetBlock->addDuplication(sourceBlock->getId());
                        } else if (!targetIsPrimary && sourceIsPrimary) {
                            sourceBlock->addDuplication(targetBlock->getId());
                        } else {
                            targetBlock->addDuplication(sourceBlock->getId());
                        }
                    }
                }
            }
        }

        // 4. 重建 BlockSet 的拓樸 (Prev/Next)
        // rewireGraph(resultSet, blockMapping, refWorkList); // 根據 Ref 路徑連線
        // 注意：對於 Qry 特有的路徑，可能需要額外的連線邏輯，
        // 但因為我們主要是 Merge 到 Ref backbone，通常依賴 Ref 的連線就夠了，
        // 除非處理 insertion bubbles。

        std::cout << "[Merger] 4. Rewiring Graph Topology...\n";

        auto addEdgeUnique = [](std::shared_ptr<Block> from, std::shared_ptr<Block> to) {
        if (!from || !to || from == to) return; 
        
        bool alreadyConnected = false;
        for (const auto& weak_next : from->getNextBlocks()) {
            if (auto next_ptr = weak_next.lock()) {
                if (next_ptr->getId() == to->getId()) {
                    alreadyConnected = true;
                    break;
                }
            }
        }
        
        if (!alreadyConnected) {
            from->addNextBlock(to);
            to->addPrevBlock(from);
        }
    };

    // [修改重點 5] 傳入獨立的 Map 進行對應
    auto traceAndConnect = [&](const std::list<std::shared_ptr<Block>>& workList, 
                               const std::map<Block::ID, std::shared_ptr<Block>>& mapping) {
        if (workList.empty()) return;
        
        auto it = workList.begin();
        auto prevOldBlock = *it;
        ++it;
        
        for (; it != workList.end(); ++it) {
            auto currOldBlock = *it;
            
            // 安全地透過 at 或 find 來取值，避免 map 自動插入空指標
            auto prevNewIt = mapping.find(prevOldBlock->getId());
            auto currNewIt = mapping.find(currOldBlock->getId());
            
            if (prevNewIt != mapping.end() && currNewIt != mapping.end()) {
                addEdgeUnique(prevNewIt->second, currNewIt->second);
            }
            
            prevOldBlock = currOldBlock;
        }
    };

    // 分別用對應的 Map 去走訪 Ref 和 Qry 的原始路徑
    traceAndConnect(refWorkList, refBlockMapping);
    traceAndConnect(qryWorkList, qryBlockMapping);
        
    return resultSet;
}

*/
/*
std::shared_ptr<Block> BlockManager::mergeBlockPair(BlockSet* resultSet, std::shared_ptr<Block> rBlk, std::shared_ptr<Block> qBlk, const mga::Alignment& subAln) {
        

    // TODO: Update Consensus, leave it for now
    if (subAln.inverse) {
        std::cout << "Inverse Alignment.\n";
        qBlk->reverseComplement(); 
    }

    std::string rCons = rBlk->getConsensus();
    std::string qCons = qBlk->getConsensus();
    std::string newCons = "";

    // std::cout << rCons.size() << " " << qCons.size() << '\n';

    std::map<int,int> rMap;
    std::map<int,int> qMap;
        
    std::vector<Variation> refExtraGaps, qryExtraGaps;

    int rPos = 0; 
    int qPos = 0; 
    int nPos = 0; 

    for (const auto& op : subAln.CIGAR) {
        int len = op.first;
        char type = op.second;
        std::cout << std::to_string(len) << type;
    }
    // 1. 走訪 CIGAR，構建 New Consensus 與 Mapping
    for (const auto& op : subAln.CIGAR) {
        int len = op.first;
        char type = op.second;
        if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                rMap[rPos] = nPos;
                qMap[qPos] = nPos;
                // Mismatch: Take reference
                if (rPos >= rCons.length() || qPos >= qCons.length()) {
                    std::cerr << "[Fatal Error] CIGAR length exceeds consensus length!\n";
                    std::cerr << "rPos: " << rPos << " / " << rCons.length() << "\n";
                    std::cerr << "qPos: " << qPos << " / " << qCons.length() << "\n";
                    abort(); // 提早發現是哪條 Alignment 搞鬼
                }
                newCons += rCons[rPos]; 
                
                rPos++; qPos++; nPos++;
            }
        } 
        else if (type == 'I' || type == 'S') { 
            // Gap on Ref
            refExtraGaps.push_back(Variation::createGap(nPos, nPos + len));
            for (int i = 0; i < len; ++i) {
                qMap[qPos] = nPos;
                newCons += qCons[qPos];
                qPos++; nPos++;
            }
        } 
        else if (type == 'D' || type == 'N') { 
            // Gap on Qry
            qryExtraGaps.push_back(Variation::createGap(nPos, nPos + len));
            for (int i = 0; i < len; ++i) {
                rMap[rPos] = nPos;
                newCons += rCons[rPos]; 
                rPos++; nPos++;
            }
        }
        
    }
    
    rMap[rCons.length()] = nPos;
    qMap[qCons.length()] = nPos;

    // 2. 建立新 Block
    auto newBlock = resultSet->createBlock(newCons);
    newBlock->setFromPrimary(true);

    // 簡化後的 Mapping Lambda (不用管 inverse 了)
    auto mapVariation = [&](const Variation& v, std::map<int,int>& map) {
        Variation newV = v;
        newV.start = map[v.start];
        newV.end = map[v.end];
        return newV;
    };

    // 3. 遷移並更新 Ref Sequences
    for (const auto& seq : rBlk->getSequences()) {
        SequenceInfo newSeq(seq);
        newSeq.variations.clear();
        
        for (const auto& v : seq.variations) newSeq.variations.push_back(mapVariation(v, rMap));
        for (const auto& gap : refExtraGaps) newSeq.variations.push_back(gap);
        newBlock->addSequence(newSeq);
    }

    // 4. 遷移並更新 Qry Sequences
    for (const auto& seq : qBlk->getSequences()) {
        SequenceInfo newSeq (seq);
        newSeq.variations.clear();
        
        for (const auto& v : seq.variations) newSeq.variations.push_back(mapVariation(v, qMap));
        for (const auto& gap : qryExtraGaps) newSeq.variations.push_back(gap);
        newBlock->addSequence(newSeq);
    }

    return newBlock;
}


void BlockManager::integrateRemainingBlocks(
    BlockSet* refSet, 
    BlockSet* qrySet, 
    BlockSet* mergedSet, 
    const std::vector<mga::Alignment>& remainingAlns,
    int mergedConsensusLen) 
{
    std::cout << "\n[Integration] Starting Remaining Blocks Integration...\n";

    // 常數設定 (你可以視後續測試結果調整)
    const float ALPHA = 10.0f; // Hanging segment penalty
    const float BETA = 5.0f;   // Mismatch/Gap penalty
    const int HANG_TOLERANCE = 5; // 容許的末端未對齊長度 (bp)，超過才算 Hanging

    // =========================================================
    // 1. Scoring & Sorting: 找出每個 Remaining Block 分數最負 (最好) 的 Alignment
    // =========================================================
    std::map<std::string, mga::Alignment> bestAlignments;
    std::map<std::string, float> minScores;

    for (const auto& aln : remainingAlns) {
        if (!aln.valid) continue;

        // 從 qryName (e.g., "Set1_4") 反查原本的 Block 取得長度
        std::string setPrefix = aln.qryName.substr(0, aln.qryName.find('_'));
        Block::ID origBlockId = std::stoi(aln.qryName.substr(aln.qryName.find('_') + 1));
        
        BlockSet* origSet = (setPrefix == refSet->getId()) ? refSet : qrySet;
        auto origBlock = origSet->getBlock(origBlockId);
        if (!origBlock) continue;

        int qryLen = origBlock->getConsensus().length();
        int refLen = mergedConsensusLen;

        // 計算 l (Alignment Length)
        int l = aln.refIdx.second - aln.refIdx.first; 

        // 計算 N (Hanging Segments)
        int N = 0;
        if (aln.refIdx.first > HANG_TOLERANCE) N++;
        if (refLen - aln.refIdx.second > HANG_TOLERANCE) N++;
        if (aln.qryIdx.first > HANG_TOLERANCE) N++;
        if (qryLen - aln.qryIdx.second > HANG_TOLERANCE) N++;

        // 計算 M (Mismatches + Gaps)
        int matchCount = 0;
        int M = 0;
        for (const auto& op : aln.CIGAR) {
            if (op.second == 'M' || op.second == '=') matchCount += op.first;
            else if (op.second == 'X' || op.second == 'I' || op.second == 'D') M += op.first;
        }
        // 如果有 NM tag 也可以直接拿來用，這裡假設從 CIGAR 和 aln 推算
        
        // Energy Score Formula: E = -l + alpha*N + beta*M
        float score = -l + (ALPHA * N) + (BETA * M);

        // 取最負 (Minimum) 的 Score
        if (minScores.find(aln.qryName) == minScores.end() || score < minScores[aln.qryName]) {
            minScores[aln.qryName] = score;
            bestAlignments[aln.qryName] = aln;
        }
    }

    // =========================================================
    // 2. 切點收集與 Backbone Splitting
    // =========================================================
    // 這裡我們只收集有效 (< 0) 的 Alignments 的切點
    std::vector<mga::Alignment> validIntegrations;
    std::set<int> backboneCuts;

    for (const auto& pair : bestAlignments) {
        float score = minScores[pair.first];
        if (score < 0) { // 只有分數是負的才 apply merging
            validIntegrations.push_back(pair.second);
            backboneCuts.insert(pair.second.refIdx.first);
            backboneCuts.insert(pair.second.refIdx.second);
        }
    }

    std::cout << "[Integration] Found " << validIntegrations.size() << " valid remaining blocks to integrate.\n";

    // 取得 Backbone Blocks 並進行切割
    auto backboneBlocks = mergedSet->getAllBlocks();
    std::list<std::shared_ptr<Block>> backboneList(backboneBlocks.begin(), backboneBlocks.end());
    applySplits(backboneList, backboneCuts);
    
    // 建立 Backbone Index 以便快速查找 (GlobalPos -> BlockPtr)
    auto backboneIndex = indexBlocks(backboneList);

    // =========================================================
    // 3. Cycle Prevention & Merging
    // =========================================================
    for (const auto& aln : validIntegrations) {
        std::string setPrefix = aln.qryName.substr(0, aln.qryName.find('_'));
        Block::ID origBlockId = std::stoi(aln.qryName.substr(aln.qryName.find('_') + 1));
        
        BlockSet* origSet = (setPrefix == refSet->getId()) ? refSet : qrySet;
        auto remBlock = origSet->getBlock(origBlockId);
        
        // 找到對應的 Backbone Block(s)
        auto targetBlocks = getBlocksInRange(backboneIndex, {aln.refIdx.first, aln.refIdx.second});
        if (targetBlocks.empty()) continue;

        // 為了簡化，假設這裡已經是 1 對 1 的 mapping (因為我們剛剛切過了)
        auto targetBlock = targetBlocks[0]; 

        // -----------------------------------------------------
        // Cycle Detection (防環機制)
        // -----------------------------------------------------
        bool createsCycle = false;
        std::unordered_set<std::string> backboneSeqIDs;
        
        for (const auto& seq : targetBlock->getSequences()) {
            backboneSeqIDs.insert(seq.sequence_id);
        }

        // 檢查 Remaining Block 裡面的序列是否已經存在於 Backbone Block 中
        for (const auto& seq : remBlock->getSequences()) {
            if (backboneSeqIDs.count(seq.sequence_id)) {
                createsCycle = true;
                break;
            }
        }

        if (createsCycle) {
            // [Case A] 發生 Cycle -> 不 Merge Sequence，改掛 Duplication Link
            std::cout << "  -> Cycle detected for " << aln.qryName << ". Adding Duplication Link.\n";
            targetBlock->addDuplication(remBlock->getId());
            
            // 將 Remaining Block 加進 MergedSet 當作獨立節點 (如果是第一次加)
            if (!mergedSet->getBlock(remBlock->getId())) {
                auto standaloneBlock = mergedSet->createBlock(remBlock->getConsensus());
                standaloneBlock->setFromPrimary(false);
                for (const auto& seq : remBlock->getSequences()) standaloneBlock->addSequence(seq);
            }
        } else {
            // [Case B] 無 Cycle -> 執行 Merge Sequence
            std::cout << "  -> Merging sequence from " << aln.qryName << " into Backbone.\n";
            
            // 這裡呼叫你的 CIGAR Mapping 邏輯 (類似 mergeBlockPair 裡面的 seq map)
            // 將 remBlock 的 Sequences 根據 aln.CIGAR 映射到 targetBlock 上
            // mergeSequencesIntoBackbone(targetBlock, remBlock, aln); 
        }
    }
}
*/