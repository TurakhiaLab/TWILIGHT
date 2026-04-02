#include "mga.hpp"


#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

std::vector<mga::Alignment> mga::splitSingleAlignment(const mga::Alignment& aln,const std::set<int>& refCuts,const std::set<int>& qryCuts)  {
    std::vector<mga::Alignment> frags; // 用來裝切碎的子片段

    // 1. 篩選並排序切點
    std::vector<int> rCuts;
    for (int c : refCuts) {
        if (c > aln.refIdx.first && c < aln.refIdx.second) rCuts.push_back(c);
    }

    std::vector<int> qCuts;
    for (int c : qryCuts) {
        if (c > aln.qryIdx.first && c < aln.qryIdx.second) qCuts.push_back(c);
    }

    if (aln.inverse) {
        std::sort(qCuts.rbegin(), qCuts.rend());
    } else {
        std::sort(qCuts.begin(), qCuts.end());
    }

    if (rCuts.empty() && qCuts.empty()) {
        frags.push_back(aln);
        return frags;
    }

    // 2. 準備走訪 CIGAR 進行動態切割
    int rCutIdx = 0;
    int qCutIdx = 0;

    int rPos = aln.refIdx.first;
    int qPos = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first; 
    int qDir = aln.inverse ? -1 : 1;

    int currRStart = rPos;
    int currQStart = qPos;

    mga::Alignment currAln = aln;
    currAln.CIGAR.clear(); 

    // 3. 逐一消耗 CIGAR Operations
    for (const auto& op : aln.CIGAR) {
        int len = op.first;
        char type = op.second;

        while (len > 0) {
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            int step = len;

            if (consumesRef && rCutIdx < rCuts.size()) {
                int distR = rCuts[rCutIdx] - rPos;
                if (distR > 0 && distR < step) step = distR;
            }

            if (consumesQry && qCutIdx < qCuts.size()) {
                int distQ = std::abs(qCuts[qCutIdx] - qPos);
                if (distQ > 0 && distQ < step) step = distQ;
            }

            if (!currAln.CIGAR.empty() && currAln.CIGAR.back().second == type) {
                currAln.CIGAR.back().first += step; 
            } else {
                currAln.CIGAR.push_back({step, type});
            }

            if (consumesRef) rPos += step;
            if (consumesQry) qPos += step * qDir;
            len -= step;

            // 4. 檢查是否精準踩到切點
            bool hitRef = (consumesRef && rCutIdx < rCuts.size() && rPos == rCuts[rCutIdx]);
            bool hitQry = (consumesQry && qCutIdx < qCuts.size() && qPos == qCuts[qCutIdx]);

            if (hitRef || hitQry) {
                currAln.refIdx.first = currRStart;
                currAln.refIdx.second = rPos;

                if (aln.inverse) {
                    currAln.qryIdx.first = qPos;        
                    currAln.qryIdx.second = currQStart; 
                } else {
                    currAln.qryIdx.first = currQStart;  
                    currAln.qryIdx.second = qPos;       
                }

                if (!currAln.CIGAR.empty()) {
                    frags.push_back(currAln); // 改存進 frags
                }

                currRStart = rPos;
                currQStart = qPos;
                currAln = aln; 
                currAln.CIGAR.clear();

                if (hitRef) rCutIdx++;
                if (hitQry) qCutIdx++;
            }
        }
    }

    // 5. 收尾
    if (!currAln.CIGAR.empty()) {
        currAln.refIdx.first = currRStart;
        currAln.refIdx.second = rPos;

        if (aln.inverse) {
            currAln.qryIdx.first = qPos;
            currAln.qryIdx.second = currQStart;
        } else {
            currAln.qryIdx.first = currQStart;
            currAln.qryIdx.second = qPos;
        }
        frags.push_back(currAln); // 改存進 frags
    }

    return frags;
}

std::vector<mga::Alignment> mga::splitAlignmentsByCuts( const std::vector<mga::Alignment>& alignments, const std::set<int>& refCuts, const std::set<int>& qryCuts) 
{
    std::vector<mga::Alignment> newAlignments;

    for (const auto& aln : alignments) {
        // 保留你原本的過濾邏輯
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) {
            newAlignments.push_back(aln);
            continue;
        }

        // 直接呼叫單一處理函數！
        auto frags = splitSingleAlignment(aln, refCuts, qryCuts);
        
        // 將切碎的片段接上 newAlignments 的尾端
        newAlignments.insert(newAlignments.end(), frags.begin(), frags.end());
    }

    return newAlignments;
}

/*
std::pair<SequenceInfo, SequenceInfo> SequenceInfo::split(int cut_point) {
    SequenceInfo left(*this);
    SequenceInfo right(*this);
    
    left.variations.clear();
    right.variations.clear();
    
    int gap_bases_before_cut = 0;
    
    for (const auto& var : variations) {
        if (var.end <= cut_point) {
            left.variations.push_back(var);
            if (var.type == Variation::GAP) {
                gap_bases_before_cut += (var.end - var.start);
            }
        } else if (var.start >= cut_point) {
            Variation right_var = var;
            right_var.shift(-cut_point); 
            right.variations.push_back(right_var);
        } else {
            Variation left_gap = Variation::createGap(var.start, cut_point);
            left.variations.push_back(left_gap);
            gap_bases_before_cut += (cut_point - var.start);
            
            Variation right_gap = Variation::createGap(cut_point, var.end);
            right_gap.shift(-cut_point);
            right.variations.push_back(right_gap);
        }
    }
    
    int consumed_seq_len = cut_point - gap_bases_before_cut;
    
    // 【修復】：解開註解，因為如果該 Sequence 處於反股 (Reverse Complement)，
    // 它的 Consensus 從左往右，其實是對應基因體座標的「由右往左 (大到小)」遞減！
    if (!is_reverse) {
        left.end_coordinate = start_coordinate + consumed_seq_len;
        right.start_coordinate = left.end_coordinate;
    } else {
        left.start_coordinate = end_coordinate - consumed_seq_len;
        right.end_coordinate = left.start_coordinate;
    }
    
    return {left, right};
}
*/

// Return {LeftBlock, RightBlock}
/*
std::pair<int, int> BlockSet::splitSingleBlock(int parentID, int localCut) {
    auto parent = this->getBlock(parentID);
    
    // 防呆：如果找不到 Block 就不做任何事，避免 Crash
    if (!parent) return {-1, -1};

    // Split blocks
    auto left = this->createBlock(parent->getConsensus().substr(0, localCut));
    auto right = this->createBlock(parent->getConsensus().substr(localCut));

    for (auto& seqInfo : parent->getSequences()) {
        auto splitSeqInfo = seqInfo.split(localCut);
        
        // 【新增】：檢查 left 是否包含真實序列 (非全 GAP)
        if (splitSeqInfo.first.start_coordinate != splitSeqInfo.first.end_coordinate) {
            left->addSequence(splitSeqInfo.first);
        }
        
        // 【新增】：檢查 right 是否包含真實序列 (非全 GAP)
        if (splitSeqInfo.second.start_coordinate != splitSeqInfo.second.end_coordinate) {
            right->addSequence(splitSeqInfo.second);
        }
    }

    left->addNextBlock(right);
    right->addPrevBlock(left);
    
    // 【修復】：接上新邊緣前，必須先呼叫 remove 拔除指向舊 Parent 的邊
    for(auto wp : parent->getPrevBlocks()) {
        if(auto p = wp.lock()) { 
            p->removeNextBlock(parent); // 拔除舊邊
            p->addNextBlock(left); 
            left->addPrevBlock(p); 
        }
    }
    for(auto wn : parent->getNextBlocks()) {
        if(auto n = wn.lock()) { 
            n->removePrevBlock(parent); // 拔除舊邊
            right->addNextBlock(n); 
            n->addPrevBlock(right); 
        }
    }

    this->deleteBlock(parent->getId());
    return {left->getId(), right->getId()};
}

std::map<int, int> BlockSet::splitBlocksByCuts(std::set<int>& cuts, bool isRef) 
{
    std::map<int, int> blocksMap;
    auto consensusBlocks = this->getConsensusBlocks();
    int start = 0, end = 0;

    for (auto blkID : consensusBlocks) {
        auto blk = this->getBlock(blkID);
        start = end;
        end = start + blk->getConsensus().length();
        
        auto currentBlk = blk;
        int currentStart = start;

        for (int cut : cuts) {
            if (cut > currentStart && cut < end) {
                // 【修復】：必須使用 currentBlk->getId()，因為原始的 blkID 可能已經被切碎刪除了！
                auto parts = splitSingleBlock(currentBlk->getId(), cut - currentStart);
                blocksMap[currentStart] = parts.first;  
                
                currentBlk = this->getBlock(parts.second);
                currentStart = cut;
            }
        }
        blocksMap[currentStart] = currentBlk->getId();
    }
    return blocksMap;
}
*/
#include <numeric> // for std::iota

// 輔助結構：Disjoint Set Union (用於 Grouping)
struct UnionFind {
    std::map<Block::ID, Block::ID> parent;
    Block::ID find(Block::ID i) {
        if (parent.find(i) == parent.end()) parent[i] = i;
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    }
    void unite(Block::ID i, Block::ID j) {
        Block::ID rootI = find(i);
        Block::ID rootJ = find(j);
        if (rootI != rootJ) parent[rootI] = rootJ;
    }
};

/*
BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, std::vector<mga::Alignment>& alignments) {
    std::cout << "[Merger] 1. Concatenating Involved Blocks into Super-Blocks...\n";

    // ==========================================
    // Phase 0: 串接所有 Consensus Blocks 成為單一 Super-Block (保留你原本的邏輯)
    // ==========================================
    auto concatenateBlocks = [&](BlockSet* bSet, Block::ID superId) -> std::shared_ptr<Block> {
        auto consensusBlocks = bSet->getConsensusBlocks();
        if (consensusBlocks.empty()) return nullptr;

        std::string super_consensus = "";
        std::vector<SequenceInfo> super_seqs;

        for (auto blkID : consensusBlocks) {
            auto blk = bSet->getBlock(blkID);
            int current_super_len = super_consensus.length();
            int blk_len = blk->getConsensus().length();
            
            // 拼接 Consensus
            super_consensus += blk->getConsensus();
            
            std::vector<SequenceInfo> next_super_seqs;
            std::vector<bool> blk_seq_used(blk->getSequences().size(), false);
            
            // 掃描目前累積的 Sequences，嘗試與新 Block 的 Sequences 進行配對
            for (auto& old_seq : super_seqs) {
                bool extended = false;
                for (size_t j = 0; j < blk->getSequences().size(); ++j) {
                    if (blk_seq_used[j]) continue;
                    auto& new_seq = blk->getSequences()[j];
                    
                    // Case 1-1: ID 相同且 Coordinate 連續 (無縫接軌)
                    if (old_seq.sequence_id == new_seq.sequence_id && 
                        old_seq.end_coordinate == new_seq.start_coordinate &&
                        old_seq.is_reverse == new_seq.is_reverse) {
                        
                        SequenceInfo merged_seq = old_seq;
                        merged_seq.end_coordinate = new_seq.end_coordinate;
                        
                        // 將新 Block 內的 Variations 平移後加入
                        for (auto var : new_seq.variations) {
                            var.shift(current_super_len);
                            merged_seq.variations.push_back(var);
                        }
                        next_super_seqs.push_back(merged_seq);
                        blk_seq_used[j] = true;
                        extended = true;
                        break;
                    }
                }
                
                // Case 1-2: 如果沒有找到延續的 Sequence -> 代表在這個新 Block 區段是缺失的
                if (!extended) {
                    SequenceInfo gap_extended_seq = old_seq;
                    // 加入一段長度等同於 Block2 consensus 的 GAP
                    gap_extended_seq.variations.push_back(Variation::createGap(current_super_len, current_super_len + blk_len));
                    next_super_seqs.push_back(gap_extended_seq);
                }
            }
            
            // Case 1-3 & 1-4: 新 Block 裡有，但前面沒配對到的 Sequence (全新或是不連續的)
            for (size_t j = 0; j < blk->getSequences().size(); ++j) {
                if (!blk_seq_used[j]) {
                    auto& new_seq = blk->getSequences()[j];
                    SequenceInfo padded_seq = new_seq;
                    
                    // 先平移它本身的 Variations
                    for (auto& var : padded_seq.variations) {
                        var.shift(current_super_len);
                    }
                    
                    // 因為它前面沒有出現過，必須在前面補上一段長度等同於目前 Super-Block 的 GAP
                    if (current_super_len > 0) {
                        padded_seq.variations.insert(padded_seq.variations.begin(), Variation::createGap(0, current_super_len));
                    }
                    next_super_seqs.push_back(padded_seq);
                }
            }
            
            super_seqs = std::move(next_super_seqs);
        }

        auto super_block = std::make_shared<Block>(superId, super_consensus);
        for (const auto& seq : super_seqs) {
            super_block->addSequence(seq);
        }
        return super_block;
    };
    auto refSuperBlock = concatenateBlocks(refSet, 9999991); 
    auto qrySuperBlock = concatenateBlocks(qrySet, 9999992);


    // ==========================================
    // Phase 1: 切割 Backbone (Primary + Secondary 的切點)
    // ==========================================
    std::cout << "[Merger] 2. Extracting Cut Points from Alignments...\n";
    std::set<int> refCuts, qryCuts;
    for (const auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;
        if (aln.refIdx.first != -1) { refCuts.insert(aln.refIdx.first); refCuts.insert(aln.refIdx.second); }
        if (aln.qryIdx.first != -1) { qryCuts.insert(aln.qryIdx.first); qryCuts.insert(aln.qryIdx.second); }
    }
    
    auto refBlocksMap = refSet->splitBlocksByCuts(refCuts, true);
    auto qryBlocksMap = qrySet->splitBlocksByCuts(qryCuts, false);

    // 建立全局 Block 查找池，方便跨 Set 拿取 Block
    std::map<Block::ID, std::shared_ptr<Block>> globalBlockPool;
    for (const auto& kv : refBlocksMap) globalBlockPool[kv.second] = refSet->getBlock(kv.second);
    for (const auto& kv : qryBlocksMap) globalBlockPool[kv.second] = qrySet->getBlock(kv.second);


    // ==========================================
    // Phase 2: Grouping (Primary & Secondary)
    // ==========================================
    std::cout << "[Merger] 3. Grouping Homologous Blocks (Primary & Secondary)...\n";
    UnionFind uf;

    // 先把所有的細碎 Block 加進 DSU (自己就是一個 Group)
    for (const auto& kv : globalBlockPool) uf.find(kv.first);

    for (const auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;

        int rBlkId = refBlocksMap[aln.refIdx.first];
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qBlkId = qryBlocksMap[qMin];

        // 不管是 1-to-1 matching (Primary), 還是 Duplication/Paralog (Secondary)
        // 只要同源，全部歸類到同一個 Group！
        uf.unite(rBlkId, qBlkId); 
    }


    // ==========================================
    // Phase 3: 處理 Remaining Blocks (Split & Join Groups)
    // ==========================================
    std::cout << "[Merger] 4. Projecting and Grouping Remaining Blocks...\n";
    
    std::vector<std::shared_ptr<Block>> remainingBlocks = refSet->getRemainingBlocks();
    auto qryRemainingBlocks = qrySet->getRemainingBlocks();
    remainingBlocks.insert(remainingBlocks.end(), qryRemainingBlocks.begin(), qryRemainingBlocks.end());
    
    // 2. 找出每個 Remaining Block 分數最高的 Alignment
    std::map<Block::ID, mga::Alignment> bestRemainingAln;
    
    // 輔助變數：用來記錄目前存下來的 best alignment，它對齊的對象是誰（方便做 main vs remaining 的優先權判斷）
    std::map<Block::ID, std::string> bestTargetName;

    // 輔助函數：從名稱 (ex: "Set1_123") 萃取出 Block::ID (123)
    auto extractBlockIdFromName = [](const std::string& name) -> Block::ID {
        size_t pos = name.find_last_of('_');
        if (pos != std::string::npos && pos + 1 < name.length()) {
            std::string idStr = name.substr(pos + 1);
            if (idStr != "main") {
                return std::stoull(idStr); // 假設 Block::ID 對應的是 unsigned long long
            }
        }
        return -1; 
    };

    // 輔助函數：判斷這個名字是否為 main backbone
    auto isMain = [](const std::string& name) -> bool {
        return name.find("_main") != std::string::npos;
    };
    /*

    // 輔助函數：判斷 newAln 是否比 oldAln 更好
    auto isBetterAlignment = [&isMain](const mga::Alignment& newAln, const std::string& newTarget,
                                       const mga::Alignment& oldAln, const std::string& oldTarget) -> bool {
        bool newIsMain = isMain(newTarget);
        bool oldIsMain = isMain(oldTarget);

        // 規則 1: Main 優先權絕對大於 Remaining
        if (newIsMain && !oldIsMain) return true;  // 新的對到 Main，舊的對到 Remaining -> 新的贏
        if (!newIsMain && oldIsMain) return false; // 舊的對到 Main，新的對到 Remaining -> 舊的贏

        // 規則 2: 如果都對到 Main，或是都對到 Remaining，則比較 alnScore
        return newAln.alnScore > oldAln.alnScore;
    };

    // 遍歷所有 Alignment 進行過濾
    for (const auto& aln : alignments) {
        if (!aln.valid || aln.type != mga::REMAINING_ALN) continue;

        bool refIsMain = isMain(aln.refName);
        bool qryIsMain = isMain(aln.qryName);

        // 情境 A: 如果 Ref 端是 Remaining Block (它對齊到了 Qry)
        if (!refIsMain) {
            Block::ID refId = extractBlockIdFromName(aln.refName);
            if (refId != (Block::ID)-1) {
                // 如果是第一次遇到，或者新的 Alignment 條件更好
                if (bestRemainingAln.find(refId) == bestRemainingAln.end() ||
                    isBetterAlignment(aln, aln.qryName, bestRemainingAln[refId], bestTargetName[refId])) {
                    
                    bestRemainingAln[refId] = aln;
                    bestTargetName[refId] = aln.qryName; // 記錄它對齊的目標名稱
                }
            }
        }

        // 情境 B: 如果 Qry 端是 Remaining Block (它對齊到了 Ref)
        // (注意：有可能這條 Alignment 兩端都是 Remaining，所以兩個 if 都要獨立檢查！)
        if (!qryIsMain) {
            Block::ID qryId = extractBlockIdFromName(aln.qryName);
            if (qryId != (Block::ID)-1) {
                // 如果是第一次遇到，或者新的 Alignment 條件更好
                if (bestRemainingAln.find(qryId) == bestRemainingAln.end() ||
                    isBetterAlignment(aln, aln.refName, bestRemainingAln[qryId], bestTargetName[qryId])) {
                    
                    bestRemainingAln[qryId] = aln;
                    bestTargetName[qryId] = aln.refName; // 記錄它對齊的目標名稱
                }
            }
        }
    }

    const int FUZZY_TOLERANCE = 15;
    
    for (auto& remBlock : remainingBlocks) {
        Block::ID remId = remBlock->getId();
        globalBlockPool[remId] = remBlock;

        if (bestRemainingAln.count(remId)) {
            auto aln = bestRemainingAln[remId];
            bool targetIsBackbone = (aln.refName.find("_main") != std::string::npos);

            if (targetIsBackbone) {
                // 3-1: Project cuts & Split
                std::set<int> projectedCuts = projectCutsViaCigar(refCuts, aln); 
                
                // Fuzzy 邊界過濾
                std::set<int> validCuts;
                int remLen = remBlock->getConsensus().length();
                for(int c : projectedCuts) {
                    if (c >= FUZZY_TOLERANCE && (remLen - c) >= FUZZY_TOLERANCE) validCuts.insert(c);
                }

                // 切割 Remaining Block
                auto fragMap = qrySet->splitBlocksByCuts(validCuts, false); // 假設 remBlock 在 qrySet

                // 將碎片加入對應的 Group
                for (const auto& fragKv : fragMap) {
                    Block::ID fragId = fragKv.second;
                    globalBlockPool[fragId] = qrySet->getBlock(fragId);

                    int targetRefPos = mapQryPosToRefPos(fragKv.first, aln); 
                    if (refBlocksMap.count(targetRefPos)) {
                        Block::ID targetGroupBlkId = refBlocksMap[targetRefPos];
                        uf.unite(targetGroupBlkId, fragId); // 碎片加入 Backbone 群組
                    } else {
                        // 如果超出 Alignment 範圍 (Overhang)，自己成為一個 Group
                        uf.find(fragId);
                    }
                }
            } else {
                // 3-2: Remaining 互相 Align (Floating Islands)
                Block::ID targetRemId = extractBlockIdFromName(aln.refName);
                uf.unite(targetRemId, remId); 
            }
        } else {
            // 3-3: 孤兒自己成組
            uf.find(remId);
        }
    }


    // ==========================================
    // Phase 4: Iterative Physical Merging per Group
    // ==========================================
    std::cout << "[Merger] 5. Iterative Physical Merging within Groups...\n";
    
    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* resultSet = createBlockSet(newID);

    std::map<Block::ID, std::vector<Block::ID>> groupedBlocks;
    for (auto& kv : uf.parent) groupedBlocks[uf.find(kv.first)].push_back(kv.first);

    std::set<std::shared_ptr<Block>> workingPool;
    std::map<Block::ID, std::shared_ptr<Block>> oldToNewBlockMap;

    for (const auto& group : groupedBlocks) {
        const auto& members = group.second;
        if (members.empty()) continue;

        // A. 挑選 Base Block (優先挑選 ref_main 的片段)
        std::shared_ptr<Block> baseBlock = globalBlockPool[members[0]]; 
        for (Block::ID id : members) {
            // 這裡可以寫一個條件：如果 globalBlockPool[id] 來自 ref_main，就設為 baseBlock
        }

        // B. 建立「座標變換追蹤器」 (Tracking Coordinates)
        // 初始狀態：1-to-1 映射
        std::vector<int> currentCoordMap(baseBlock->getConsensus().length() + 1);
        std::iota(currentCoordMap.begin(), currentCoordMap.end(), 0);

        // C. 依序將其他 Member 併入 Base Block
        for (Block::ID memberId : members) {
            if (memberId == baseBlock->getId()) continue;
            auto memberBlock = globalBlockPool[memberId];

            // 1. 取得 Base 和 Member 之間的原始 Alignment CIGAR
            // (這需要你從 alignments 裡面撈出對應的 aln)
            mga::Cigar origCigar = getOriginalCigar(baseBlock->getId(), memberId, alignments); 
            bool inverse = false; // 取得原本的 strand 資訊

            // 2. 轉換 CIGAR！利用 currentCoordMap 把舊的 Ref 座標轉換成現在長胖後的座標
            mga::Cigar adjustedCigar = adjustCigarWithMap(origCigar, currentCoordMap);

            // 3. 執行物理合併！(你需要稍微修改 mergeTwoBlocks，讓它回傳這次合併的 refOldToNew 陣列)
            std::vector<int> stepCoordMap;
            baseBlock = resultSet->mergeTwoBlocksWithTracking(baseBlock, memberBlock, adjustedCigar, inverse, stepCoordMap);

            // 4. 更新全局的座標追蹤器
            // currentCoordMap[i] 現在指向一個舊座標 X，經過這次合併，X 變成了 stepCoordMap[X]
            for (size_t i = 0; i < currentCoordMap.size(); ++i) {
                currentCoordMap[i] = stepCoordMap[currentCoordMap[i]];
            }
        }

        // 完成這個群組的合併，登記到 workingPool
        workingPool.insert(baseBlock);
        for (Block::ID memberId : members) {
            oldToNewBlockMap[memberId] = baseBlock;
        }
    }


    // ==========================================
    // Phase 5: 拓撲重建 (Topology Reconstruction)
    // ==========================================
    std::cout << "[Merger] 6. Rewiring Graph Edges...\n";
    
    // (將 workingPool 放進 resultSet, 清除舊連結)
    for (auto& blk : workingPool) {
        resultSet->addBlock(blk);
        blk->clearLinkages();
    }

    // 利用 oldToNewBlockMap 重建 prev/next 連線 
    // (你原本的 rebuildEdges 邏輯完全相容！)
    // ...

    return resultSet;
    
    return nullptr;
}
*/

