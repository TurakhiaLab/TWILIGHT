#include "mga.hpp"


#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric> // for std::iota
#include <tuple>

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

// ==========================================
// Helper 1: 根據 Map 動態補償 Gap (絕對數學嚴謹版)
// ==========================================
mga::Cigar adjustCigarWithMap(const mga::Cigar& origCigar, const std::vector<int>& coordMap, int superBaseLen) {
    mga::Cigar adjCigar;
    int origRefPos = 0; 

    auto addOp = [&](int len, char op) {
        if (len <= 0) return;
        if (!adjCigar.empty() && adjCigar.back().second == op) adjCigar.back().first += len;
        else adjCigar.push_back({len, op});
    };

    // 1. 補齊 Prefix Gaps (在 Member Block 之前的 SuperBase 序列)
    if (coordMap.size() > 0 && coordMap[0] > 0) {
        addOp(coordMap[0], 'D');
    }

    // 2. 走訪並轉換原始 CIGAR
    for (auto op : origCigar) {
        int len = op.first; char type = op.second;

        if (type == 'I' || type == 'S' || type == 'H') {
            addOp(len, 'I'); // 純 Qry 消耗
        }
        else if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                if (origRefPos + 1 < coordMap.size()) {
                    int p1 = coordMap[origRefPos];
                    int p2 = coordMap[origRefPos + 1];

                    if (p2 == p1) {
                        // 座標重疊：代表這個鹼基在 SuperBase 中不存在，轉為 Insertion
                        addOp(1, 'I');
                    } else {
                        // 正常匹配
                        addOp(1, type);
                        // 補齊兩點之間的 Gap (SuperBase 中多出來的序列)
                        int gaps = p2 - p1 - 1;
                        if (gaps > 0) addOp(gaps, 'D');
                    }
                    origRefPos++;
                }
            }
        }
        else if (type == 'D') {
            for (int i = 0; i < len; ++i) {
                if (origRefPos + 1 < coordMap.size()) {
                    int p1 = coordMap[origRefPos];
                    int p2 = coordMap[origRefPos + 1];

                    if (p2 > p1) {
                        addOp(1, 'D'); // Qry 本來就沒有，SuperBase 有，所以是 D
                        int gaps = p2 - p1 - 1;
                        if (gaps > 0) addOp(gaps, 'D');
                    }
                    origRefPos++;
                }
            }
        }
    }

    // 3. 補齊 Suffix Gaps (在 Member Block 之後的 SuperBase 序列)
    int currentEnd = coordMap.empty() ? 0 : coordMap.back();
    if (superBaseLen > currentEnd) {
        addOp(superBaseLen - currentEnd, 'D');
    }

    return adjCigar;
}

// ==========================================
// Helper 2: 根據 Consensus 座標裁切 Sub-CIGAR 並強制對齊
// ==========================================
std::tuple<bool, mga::Cigar, bool> extractSubCigar(
    Block::ID baseId, Block::ID memId, 
    const std::vector<mga::Alignment>& alignments,
    const std::map<Block::ID, std::pair<int, int>>& blockConsensusCoords,
    const std::map<Block::ID, bool>& isRefBlockMap) 
{
    bool baseFromRef = isRefBlockMap.at(baseId);
    bool memFromRef  = isRefBlockMap.at(memId);

    int tgtRStart = baseFromRef ? blockConsensusCoords.at(baseId).first  : blockConsensusCoords.at(memId).first;
    int tgtREnd   = baseFromRef ? blockConsensusCoords.at(baseId).second : blockConsensusCoords.at(memId).second;
    int tgtQStart = !baseFromRef ? blockConsensusCoords.at(baseId).first  : blockConsensusCoords.at(memId).first;
    int tgtQEnd   = !baseFromRef ? blockConsensusCoords.at(baseId).second : blockConsensusCoords.at(memId).second;

    int targetBaseLen = blockConsensusCoords.at(baseId).second - blockConsensusCoords.at(baseId).first;
    int targetMemLen  = blockConsensusCoords.at(memId).second - blockConsensusCoords.at(memId).first;

    for (const auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;

        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        if (tgtRStart >= rMax || tgtREnd <= rMin) continue;
        if (tgtQStart >= qMax || tgtQEnd <= qMin) continue;

        mga::Cigar subCigar;
        int currR = rMin;

        for (auto op : aln.CIGAR) {
            int len = op.first; char type = op.second;
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            if (!consumesRef) {
                if (currR >= tgtRStart && currR <= tgtREnd) subCigar.push_back(op);
                continue;
            }
            if (currR >= tgtREnd) break; 

            int overlapStart = std::max(tgtRStart, currR);
            int overlapEnd = std::min(tgtREnd, currR + len);
            if (overlapStart < overlapEnd) {
                subCigar.push_back({overlapEnd - overlapStart, type});
            }
            currR += len;
        }

        // 翻轉方向
        if (!baseFromRef) {
            for (auto& op : subCigar) {
                if (op.second == 'I') op.second = 'D';
                else if (op.second == 'D') op.second = 'I';
            }
        }

        // 【核心修復】：1-bp 微步進對齊引擎，強制讓 CIGAR 貼合真實長度！
        mga::Cigar finalCigar;
        int curBase = 0, curMem = 0;

        auto addOp = [&](char t) {
            if (!finalCigar.empty() && finalCigar.back().second == t) finalCigar.back().first++;
            else finalCigar.push_back({1, t});
        };

        for (auto op : subCigar) {
            int l = op.first; char t = op.second;
            if (t == 'S' || t == 'H') t = 'I';

            bool bCons = (t == 'M' || t == '=' || t == 'X' || t == 'D');
            bool mCons = (t == 'M' || t == '=' || t == 'X' || t == 'I');

            for (int i = 0; i < l; ++i) {
                bool useB = false, useM = false;
                
                // 走訪並控制不超出上限
                if (bCons && curBase < targetBaseLen) { useB = true; curBase++; }
                if (mCons && curMem < targetMemLen) { useM = true; curMem++; }

                if (useB && useM) addOp((t == 'M' || t == '=' || t == 'X') ? t : 'M');
                else if (useB) addOp('D');
                else if (useM) addOp('I');
            }
        }

        // 如果不足，強制補尾刀 Gap
        if (curBase < targetBaseLen) {
            addOp('D'); finalCigar.back().first += (targetBaseLen - curBase - 1);
        }
        if (curMem < targetMemLen) {
            addOp('I'); finalCigar.back().first += (targetMemLen - curMem - 1);
        }

        return {true, finalCigar, aln.inverse};
    }
    
    return {false, {}, false};
}

// ==========================================
// 主函數：Graph Merge
// ==========================================
BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, std::vector<mga::Alignment>& alignments) {
    bool DEBUG_MODE = true;

    if (DEBUG_MODE) std::cout << "\n========================================================\n"
                              << "=== GRAPH MERGE START: " << refSet->getId() << " + " << qrySet->getId() << " ===\n"
                              << "========================================================\n";

    // ==========================================
    // Phase 0: 串接 Consensus Blocks (加入 Debug 驗證)
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 0] Concatenating Involved Blocks into Super-Blocks...\n";
    
    auto concatenateBlocks = [&](BlockSet* bSet, Block::ID superId) -> std::shared_ptr<Block> {
        auto consensusBlocks = bSet->getRepresentativeBlocks(); 
        if (consensusBlocks.empty()) return nullptr;
        
        std::string super_consensus = "";
        
        struct Track {
            std::string seqID;
            Segment seg;
        };
        std::vector<Track> super_tracks;
    
        for (auto blkID : consensusBlocks) {
            auto blk = bSet->getBlock(blkID);
            int current_super_len = super_consensus.length();
            int blk_len = blk->getConsensus().length();
            
            super_consensus += blk->getConsensus();
            std::vector<Track> next_super_tracks;
            std::map<std::string, std::set<int>> used_blk_segs;
            
            for (auto& old_track : super_tracks) {
                bool extended = false;
                std::string id = old_track.seqID;
                
                auto it = blk->getSequences().find(id);
                if (it != blk->getSequences().end()) {
                    SequenceInfo blk_info = it->second; 
                    for (auto& kv : blk_info.getSegments()) {
                        int seg_start = kv.first;
                        Segment& new_seg = kv.second;
                        
                        if (used_blk_segs[id].count(seg_start)) continue;
                        
                        if (old_track.seg.getEnd() == new_seg.getStart() && 
                            old_track.seg.isReverse() == new_seg.isReverse()) {
                            
                            Track merged_track = old_track;
                            merged_track.seg.setEnd(new_seg.getEnd());
                            
                            for (auto var : new_seg.getVariations()) {
                                var.shift(current_super_len);
                                merged_track.seg.getVariations().push_back(var);
                            }
                            
                            next_super_tracks.push_back(merged_track);
                            used_blk_segs[id].insert(seg_start);
                            extended = true;
                            break;
                        }
                    }
                }
                
                if (!extended) {
                    Track gap_extended_track = old_track;
                    gap_extended_track.seg.getVariations().push_back(Variation::createGap(current_super_len, current_super_len + blk_len));
                    next_super_tracks.push_back(gap_extended_track);
                }
            }
            
            for (auto& seqPair : blk->getSequences()) {
                std::string id = seqPair.first;
                SequenceInfo blk_info = seqPair.second;
                
                for (auto& kv : blk_info.getSegments()) {
                    int seg_start = kv.first;
                    Segment new_seg = kv.second;
                    
                    if (!used_blk_segs[id].count(seg_start)) {
                        Track padded_track;
                        padded_track.seqID = id;
                        padded_track.seg = new_seg;
                        
                        for (auto& var : padded_track.seg.getVariations()) {
                            var.shift(current_super_len);
                        }
                        
                        if (current_super_len > 0) {
                            padded_track.seg.getVariations().insert(
                                padded_track.seg.getVariations().begin(), 
                                Variation::createGap(0, current_super_len)
                            );
                        }
                        next_super_tracks.push_back(padded_track);
                    }
                }
            }
            super_tracks = std::move(next_super_tracks);
        }
    
        auto super_block = std::make_shared<Block>(superId, super_consensus);
        std::map<std::string, SequenceInfo> final_seq_map;
        
        for (auto& track : super_tracks) {
            if (final_seq_map.find(track.seqID) == final_seq_map.end()) {
                final_seq_map[track.seqID] = SequenceInfo(track.seqID);
            }
            final_seq_map[track.seqID].addSegments(track.seg);
        }
        
        for (auto& kv : final_seq_map) {
            super_block->addSequence(kv.second);
        }
        
        if (DEBUG_MODE) {
            int expected_len = super_block->getConsensus().length();
            std::cout << "[DEBUG-CHECK] Validating SuperBlock " << superId << " (Target Len: " << expected_len << ")...\n";
            bool all_passed = true;
            for (auto& seqPair : super_block->getSequences()) {
                for (auto& segPair : seqPair.second.getSegments()) {
                    Segment& seg = segPair.second;
                    int coord_diff = std::abs(seg.getEnd() - seg.getStart());
                    int gap_len = 0;
                    for (auto& var : seg.getVariations()) {
                        if (var.getType() == Variation::GAP) {
                            gap_len += (var.getEnd() - var.getStart());
                        }
                    }
                    int total_calculated_len = coord_diff + gap_len;
                    if (total_calculated_len != expected_len) {
                        std::cerr << "  ❌ [WARNING] Length mismatch in Seq: " << seqPair.first 
                                  << " | Coord Diff: " << coord_diff << " + Gaps: " << gap_len 
                                  << " = " << total_calculated_len << " (Expected: " << expected_len << ")\n";
                        all_passed = false;
                    }
                }
            }
            if (all_passed) std::cout << "  ✅ All segments dynamically sum to " << expected_len << " bp perfectly!\n";
        }
        
        return super_block;
    };

    DEBUG_MODE = false;
    auto refSuperBlock = concatenateBlocks(refSet, 9999991); 
    auto qrySuperBlock = concatenateBlocks(qrySet, 9999992);
    BlockSet refSuperSet ("ref_super");
    BlockSet qrySuperSet ("qry_super");
    refSuperBlock = refSuperSet.addBlock(refSuperBlock);
    qrySuperBlock = qrySuperSet.addBlock(qrySuperBlock);
    DEBUG_MODE = true;

    // ==========================================
    // Phase 1: 使用 Dictionary + splitSingleBlock 切割
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 1] Extracting Cuts and Splitting SuperBlocks...\n";
    
    std::set<int> refCuts, qryCuts;
    std::vector<mga::Alignment> validAlignments; // [新增] 用來取代原本的 alignments

    for (auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;
        
        // 【新增條件】：如果 ref 端或 qry 端小於 100 base，直接丟棄這個 Alignment
        int rLen = std::abs(aln.refIdx.second - aln.refIdx.first);
        int qLen = std::abs(aln.qryIdx.second - aln.qryIdx.first);
        if (rLen < 100 || qLen < 100) {
            aln.setValid2False();
            continue;
        }

        refCuts.insert(aln.refIdx.first); refCuts.insert(aln.refIdx.second);
        qryCuts.insert(aln.qryIdx.first); qryCuts.insert(aln.qryIdx.second);
        
        validAlignments.push_back(aln); // 保留合格的 Alignment
    }

    // if (DEBUG_MODE) {
    //     std::cout << "Reference Cut Points: \n";
    //     for (auto cut: refCuts) std::cout << cut << " "; std::cout << "\n";
    //     std::cout << "Query Cut Points: \n";
    //     for (auto cut: qryCuts) std::cout << cut << " "; std::cout << "\n";
    // }
    
    
    std::map<int, BlockSet::SegNode> refDict;
    refDict[0] = {0, (int)refSuperBlock->getConsensus().length(), refSuperBlock->getId()};

    std::map<int, BlockSet::SegNode> qryDict;
    qryDict[0] = {0, (int)qrySuperBlock->getConsensus().length(), qrySuperBlock->getId()};

    auto splitDictBlock = [&](std::map<int, BlockSet::SegNode>& dict, BlockSet* bSet, int cutPos) {
        if (DEBUG_MODE) std::cout << "  [DEBUG-SPLIT] Requested cut at " << cutPos << " -> ";
        
        auto it = dict.upper_bound(cutPos);
        if (it == dict.begin()) {
            if (DEBUG_MODE) std::cout << "Ignored (Out of lower bounds)\n";
            return; 
        }
        it--;
        
        if (it->first == cutPos) {
            if (DEBUG_MODE) std::cout << "Ignored (Already a boundary)\n";
            return; 
        }

        int start = it->second.start;
        int end = it->second.end;
        Block::ID targetBlkId = it->second.blkId;
        std::shared_ptr<Block> targetBlk = bSet->getBlock(targetBlkId);
        
        if (!targetBlk) {
            if (DEBUG_MODE) std::cout << "Failed (Target Block ID " << targetBlkId << " not found)\n";
            return;
        }

        int localCut = cutPos - start;

        if (localCut <= 0 || localCut >= targetBlk->getConsensus().length()) {
            if (DEBUG_MODE) std::cout << "Ignored (Invalid localCut: " << localCut << " for Block " << targetBlkId << " length " << targetBlk->getConsensus().length() << ")\n";
            return; 
        }

        if (DEBUG_MODE) std::cout << "Cutting Block " << targetBlkId << " [" << start << ", " << end << "] at local idx " << localCut << " ... ";

        auto parts = bSet->splitSingleBlock(targetBlkId, localCut);
        if (parts.first == (uint64_t)-1) {
            if (DEBUG_MODE) std::cout << "Failed (splitSingleBlock returned -1)\n";
            return; 
        }

        if (DEBUG_MODE) std::cout << "Success! New Blocks: " << parts.first << " & " << parts.second << "\n";

        dict.erase(it);
        dict[start] = {start, cutPos, parts.first};
        dict[cutPos] = {cutPos, end, parts.second};
    };

    DEBUG_MODE = false;
    for (int cut : refCuts) splitDictBlock(refDict, &refSuperSet, cut);
    for (int cut : qryCuts) splitDictBlock(qryDict, &qrySuperSet, cut);
    DEBUG_MODE = true;

    if (DEBUG_MODE) {
        std::cout << "  -> Ref SuperBlock split into " << refDict.size() << " atomic blocks.\n";
        std::cout << "  -> Qry SuperBlock split into " << qryDict.size() << " atomic blocks.\n";
    }

    std::map<int, Block::ID> refBlocksMap;
    for (auto& kv : refDict) {
        refBlocksMap[kv.first] = kv.second.blkId;
    }

    std::map<int, Block::ID> qryBlocksMap;
    // 遍歷 qryDict，把切好的 qry block 逐一加進 refSuperSet
    for (auto& kv : qryDict) {
        int qryStart = kv.first;              // Qry 的 Genomic Coordinate
        Block::ID oldQryId = kv.second.blkId; // 在 qrySuperSet 裡的舊 ID

        // 從 qrySuperSet 拿出切好的 Block
        std::shared_ptr<Block> qryBlock = qrySuperSet.getBlock(oldQryId);

        if (qryBlock) {
            // 搬家：加進 refSuperSet，取得全新 ID
            auto newUnifiedBlock = refSuperSet.addBlock(qryBlock);

            // 1. 記錄到供 Grouping 使用的快速查詢表
            qryBlocksMap[qryStart] = newUnifiedBlock->getId();

            // 2. 【關鍵新增】：直接更新 qryDict 本身裡面的記錄！
            // 因為 kv 是 auto& (參照)，所以這裡改了，map 裡面的值就會跟著改
            kv.second.blkId = newUnifiedBlock->getId();
        }
    }

    // 建立全局 Block 查找池
    std::map<Block::ID, std::shared_ptr<Block>> globalBlockPool;
    
    for (const auto& kv : refBlocksMap) {
        globalBlockPool[kv.second] = refSuperSet.getBlock(kv.second);
    }
    
    for (const auto& kv : qryBlocksMap) {
        // 【注意】：因為 Qry 已經搬家了，所以這裡也是從 refSuperSet 拿積木！
        globalBlockPool[kv.second] = refSuperSet.getBlock(kv.second); 
    }
    
    // ==========================================
    // Phase 2: Grouping Homologous Blocks
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 2] Grouping Homologous Blocks (Primary & Secondary)...\n";
    DEBUG_MODE = false;
    UnionFind uf; 
    for (const auto& kv : globalBlockPool) uf.find(kv.first); 

    int matchCount = 0;
    int alnCounter = 0;
    for (const auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;
        alnCounter++;

        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        if (DEBUG_MODE) {
            std::cout << "  [DEBUG-GROUP] Aln #" << alnCounter << " (Type " << aln.type << ") | Ref: [" << rMin << ", " << rMax << "] | Qry: [" << qMin << ", " << qMax << "]\n";
        }

        // 【修復核心】：抓取此 Alignment 範圍內涵蓋到的 所有 Ref 與 Qry 積木
        std::vector<Block::ID> rBlocks;
        auto rIt = refDict.upper_bound(rMin); 
        if (rIt != refDict.begin()) rIt--;
        while (rIt != refDict.end() && rIt->second.start < rMax) {
            rBlocks.push_back(rIt->second.blkId);
            rIt++;
        }

        std::vector<Block::ID> qBlocks;
        auto qIt = qryDict.upper_bound(qMin); 
        if (qIt != qryDict.begin()) qIt--;
        while (qIt != qryDict.end() && qIt->second.start < qMax) {
            qBlocks.push_back(qIt->second.blkId);
            qIt++;
        }

        // 把這個範圍內的所有 Ref 積木與 Qry 積木，全部拉進同一個 UnionFind 群組！
        for (Block::ID rId : rBlocks) {
            for (Block::ID qId : qBlocks) {
                uf.unite(rId, qId);
                matchCount++;
            }
        }
    }
    
    if (DEBUG_MODE) std::cout << "  -> Grouped " << matchCount << " homologous pairs.\n";



    // ==========================================
    // 建立供 Phase 4 提取 CIGAR 使用的座標與來源查詢表
    // ==========================================
    std::map<Block::ID, std::pair<int, int>> blockConsensusCoords;
    std::map<Block::ID, bool> isRefBlockMap;

    for (const auto& kv : refDict) {
        blockConsensusCoords[kv.second.blkId] = {kv.second.start, kv.second.end};
        isRefBlockMap[kv.second.blkId] = true; 
    }
    for (const auto& kv : qryDict) {
        blockConsensusCoords[kv.second.blkId] = {kv.second.start, kv.second.end};
        isRefBlockMap[kv.second.blkId] = false; 
    }

    DEBUG_MODE = true;
    // ==========================================
    // Phase 4: Iterative Spanning-Tree Merging per Group
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 4] Iterative Spanning-Tree Merging within Groups...\n";
    
    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* resultSet = createBlockSet(newID); 

    DEBUG_MODE = false;
    std::map<Block::ID, std::vector<Block::ID>> groupedBlocks;
    for (auto& kv : uf.parent) {
        groupedBlocks[uf.find(kv.first)].push_back(kv.first);
    }

    if (DEBUG_MODE) std::cout << "  -> Total unique groups to process: " << groupedBlocks.size() << "\n\n";

    std::set<std::shared_ptr<Block>> workingPool;
    std::map<Block::ID, std::shared_ptr<Block>> oldToNewBlockMap;
    int groupCounter = 1;

    for (const auto& group : groupedBlocks) {
        const auto& members = group.second;
        if (members.empty()) continue;

        if (members.size() == 1) {
            workingPool.insert(globalBlockPool[members[0]]);
            oldToNewBlockMap[members[0]] = globalBlockPool[members[0]];
            continue;
        }

        if (DEBUG_MODE) std::cout << "[Group " << groupCounter++ << "] Members: " << members.size() << "\n";

        // 1. 決定群組的 起始 Hub (優先選擇來自 Ref 或是 _main 的積木)
        Block::ID origBaseId = members[0];
        for (Block::ID id : members) {
            if (isRefBlockMap[id]) { origBaseId = id; break; }
        }

        std::shared_ptr<Block> baseBlock = globalBlockPool[origBaseId];

        if (DEBUG_MODE) std::cout << "  ├─ Initial Hub Block ID: " << origBaseId 
                                  << " (Initial Len: " << baseBlock->getConsensus().length() << " bp)\n";

        // 2. 初始化群組的 Spanning Tree 狀態
        std::map<Block::ID, std::vector<int>> coordMaps; // 記錄每一個 Original Member 的座標變化
        coordMaps[origBaseId] = std::vector<int>(baseBlock->getConsensus().length() + 1);
        std::iota(coordMaps[origBaseId].begin(), coordMaps[origBaseId].end(), 0);

        std::set<Block::ID> mergedMembers = {origBaseId};
        std::set<Block::ID> unmergedMembers;
        for (Block::ID id : members) {
            if (id != origBaseId) unmergedMembers.insert(id);
        }

        // 3. 核心迴圈：利用圖的邊緣 (Valid Alignments) 依序拉攏未合併的積木
        while (!unmergedMembers.empty()) {
            bool foundEdge = false;
            Block::ID targetMergedId = 0;
            Block::ID targetUnmergedId = 0;
            mga::Cigar bestOrigCigar;
            bool bestInverse = false;

            // 尋找任一個 "已合併積木" 與 "未合併積木" 之間的合法 Alignment
            for (Block::ID u : unmergedMembers) {
                for (Block::ID m : mergedMembers) {
                    auto alnData = extractSubCigar(m, u, alignments, blockConsensusCoords, isRefBlockMap);
                    if (std::get<0>(alnData) == true) { // 如果找到有效的 Alignment
                        targetMergedId = m;
                        targetUnmergedId = u;
                        bestOrigCigar = std::get<1>(alnData);
                        bestInverse = std::get<2>(alnData);
                        foundEdge = true;
                        break;
                    }
                }
                if (foundEdge) break;
            }

            if (!foundEdge) {
                if (DEBUG_MODE) std::cout << "  │  [WARNING] Group disconnected! Dropping " << unmergedMembers.size() << " orphaned members.\n";
                break; // 斷圖，結束此 Group
            }

            auto memberBlock = globalBlockPool[targetUnmergedId];
            if (DEBUG_MODE) std::cout << "  ├─ Merging Member " << targetUnmergedId << " (via edge from " << targetMergedId << ")\n";

            // A. 將 CIGAR 補償為指向「當前 SuperBase 的 Consensus」
            // 【核心修改】：取得當前 baseBlock (Hub) 的長度，傳給 adjustCigarWithMap 讓它自動在頭尾補齊 D
            int superBaseLen = baseBlock->getConsensus().length();
            mga::Cigar adjustedCigar = adjustCigarWithMap(bestOrigCigar, coordMaps[targetMergedId], superBaseLen);

            if (DEBUG_MODE) {
                std::cout << "  │    - Orig CIGAR: ";
                for (auto op : bestOrigCigar) std::cout << op.first << op.second;
                std::cout << "\n  │    - Adj. CIGAR: ";
                for (auto op : adjustedCigar) std::cout << op.first << op.second;
                std::cout << "\n";
            }
            
            // B. 建立即將加入的 Member 的專屬 coordMap
            int uLen = memberBlock->getConsensus().length();
            std::vector<int> uCoordMap(uLen + 1, 0);
            int tmpR = 0, tmpQ = 0;
            for (auto op : adjustedCigar) {
                int len = op.first; char type = op.second;
                if (type == 'M' || type == '=' || type == 'X') {
                    for(int i=0; i<len; ++i) { if (tmpQ < uLen) uCoordMap[tmpQ++] = tmpR++; else tmpR++; }
                } else if (type == 'D') {
                    tmpR += len;
                } else if (type == 'I') {
                    for(int i=0; i<len; ++i) { if (tmpQ < uLen) uCoordMap[tmpQ++] = tmpR; }
                } else if (type == 'S' || type == 'H') {
                    tmpQ += len;
                }
            }
            uCoordMap[uLen] = tmpR;
            coordMaps[targetUnmergedId] = uCoordMap; // 加入追蹤池

            // C. 準備 stepCoordMap 計算即將發生的共識長胖
            int adjRefLen = 0;
            for (auto op : adjustedCigar) {
                if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') adjRefLen += op.first;
            }
            std::vector<int> stepCoordMap(adjRefLen + 1, 0);
            tmpR = 0; int mPos = 0;
            for (auto op : adjustedCigar) {
                int len = op.first; char type = op.second;
                if (type == 'M' || type == '=' || type == 'X' || type == 'D') {
                    for (int i=0; i<len; ++i) { if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR++] = mPos++; }
                } else if (type == 'I') {
                    mPos += len;
                }
            }
            if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR] = mPos;

            // D. 執行真實的物理合併
            baseBlock = refSuperSet.mergeTwoBlocks(baseBlock, memberBlock, adjustedCigar, bestInverse); 

            // E. 聯動更新【所有】已在池內的 Member 的座標系！
            for (auto& kv : coordMaps) {
                for (size_t i = 0; i < kv.second.size(); ++i) {
                    if (kv.second[i] < stepCoordMap.size()) {
                        kv.second[i] = stepCoordMap[kv.second[i]];
                    }
                }
            }

            // F. 狀態更新
            mergedMembers.insert(targetUnmergedId);
            unmergedMembers.erase(targetUnmergedId);
        }

        if (DEBUG_MODE) std::cout << "  └─ Final Merged Block ID: " << baseBlock->getId() << "\n\n";

        workingPool.insert(baseBlock);
        for (Block::ID memberId : mergedMembers) {
            oldToNewBlockMap[memberId] = baseBlock;
        }
    }

    DEBUG_MODE = true;
    // ==========================================
    // Phase 5: 拓撲重建 (Topology Reconstruction)
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 5] Rewiring Graph Edges and Finalizing...\n";
    
    // ==========================================
    // Phase 5: 基於生物座標的自動拓撲重建 
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 5] Re-wiring Pangenome Graph Edges based on genomic coordinates...\n";
    
    for (auto& blk : workingPool) {
        resultSet->addBlock(blk);
        blk->clearLinkages(); 
    }

    // 【修改點 3】：利用真實序列座標排序來重建圖拓撲
    struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
    std::map<std::string, std::vector<SegRef>> seqTracks;
    
    auto allBlocks = resultSet->getAllBlocks();
    for (auto& blk : allBlocks) {
        for (auto& seqPair : blk->getSequences()) {
            for (auto& segPair : seqPair.second.getSegments()) {
                seqTracks[seqPair.first].push_back({ &segPair.second, blk });
            }
        }
    }

    for (auto& trackPair : seqTracks) {
        auto& track = trackPair.second;
        
        // 依照原始基因體座標排序
        std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
            return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
        });
        
        for (size_t i = 0; i < track.size(); ++i) {
            track[i].seg->setPrevBlock(nullptr);
            track[i].seg->setNextBlock(nullptr);
            
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

    for (auto& seq: refSet->getSequences()) resultSet->addSequence(seq);
    for (auto& seq: qrySet->getSequences()) resultSet->addSequence(seq);


    if (DEBUG_MODE) std::cout << "========================================================\n"
                              << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
                              << "========================================================\n\n";

    return resultSet;
}