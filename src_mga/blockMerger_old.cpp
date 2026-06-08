
#include "type.hpp"
#include "block.hpp"
#include "alignment.hpp"



#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <numeric> // for std::iota
#include <tuple>
#include <chrono>
#include <tbb/parallel_invoke.h>







// 輔助結構：Disjoint Set Union (用於 Grouping)
struct UnionFind {
    std::map<BlockID, BlockID> parent;
    BlockID find(BlockID i) {
        if (parent.find(i) == parent.end()) parent[i] = i;
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    }
    void unite(BlockID i, BlockID j) {
        BlockID rootI = find(i);
        BlockID rootJ = find(j);
        if (rootI != rootJ) parent[rootI] = rootJ;
    }
};

// ==========================================
// Helper 1: 根據 Map 動態補償 Gap (絕對數學嚴謹版)
// ==========================================
CigarString adjustCigarWithMap(const CigarString& origCigar, const std::vector<int>& coordMap, int superBaseLen, bool isRefTarget = true) {
    CigarString adjCigar;
    int origRefPos = 0; 

    // 定義「消耗 map 座標」與「不消耗」的運算子
    // RefTarget 模式：D 消耗 map (Deletion), I 不消耗
    // QryTarget 模式：I 消耗 map (Insertion), D 不消耗
    char consume_op = isRefTarget ? 'D' : 'I';
    char pass_op    = isRefTarget ? 'I' : 'D';

    auto addOp = [&](int len, char op) {
        if (len <= 0) return;
        if (!adjCigar.empty() && adjCigar.back().second == op) adjCigar.back().first += len;
        else adjCigar.push_back({len, op});
    };

    // 處理開頭的 Gap
    if (!coordMap.empty() && coordMap[0] > 0) {
        addOp(coordMap[0], consume_op);
    }

    for (const auto& op : origCigar) {
        int len = op.first; char type = op.second;

        if (type == 'S' || type == 'H') continue; 
        
        // 1. 不消耗 map 的運算子 (例如 RefTarget 下的 I)
        if (type == pass_op) {
            addOp(len, type); 
        } 
        // 2. 消耗 map 的運算子 (M, =, X 以及對應的 Deletion/Insertion)
        else if (type == 'M' || type == '=' || type == 'X' || type == consume_op) {
            for (int i = 0; i < len; ++i) {
                if (origRefPos + 1 < (int)coordMap.size()) {
                    int p1 = coordMap[origRefPos];
                    int p2 = coordMap[origRefPos + 1];

                    // 如果 p2 == p1，代表在目標塊上這是個 Insertion (QryTarget 下是 D 的變體)
                    if (p2 == p1) {
                        addOp(1, consume_op);
                    } else {
                        addOp(1, type);
                        int gaps = p2 - p1 - 1;
                        if (gaps > 0) addOp(gaps, consume_op);
                    }
                    origRefPos++;
                }
            }
        }
    }

    // 🚨 移除掉原本強制補全結尾 D 的邏輯，這就是導致 190D 的元兇
    return adjCigar;
}

// ==========================================
// Helper 2: 根據 Consensus 座標裁切 Sub-CIGAR 並強制對齊
// ==========================================
std::tuple<bool, CigarString, bool> extractSubCigar(
    BlockID baseId, BlockID memId, 
    const Alignments& alignments,
    const std::map<BlockID, std::pair<int, int>>& blockConsensusCoords,
    const std::map<BlockID, bool>& isRefBlockMap) 
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
        if (!aln.valid) continue;

        int rMin = std::min(aln.refIdx.first, aln.refIdx.second);
        int rMax = std::max(aln.refIdx.first, aln.refIdx.second);
        int qMin = std::min(aln.qryIdx.first, aln.qryIdx.second);
        int qMax = std::max(aln.qryIdx.first, aln.qryIdx.second);

        if (tgtRStart >= rMax || tgtREnd <= rMin) continue;
        if (tgtQStart >= qMax || tgtQEnd <= qMin) continue;

        CigarString subCigar;
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
        CigarString finalCigar;
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
BlockSet* BlockManager::merge(BlockSet* refSet, BlockSet* qrySet, std::vector<Alignment>& alignments) {
    bool DEBUG_MODE = false;
    auto time0 = std::chrono::high_resolution_clock::now();
    if (DEBUG_MODE) std::cout << "\n========================================================\n"
                              << "=== GRAPH MERGE START: " << refSet->getId() << " + " << qrySet->getId() << " ===\n"
                              << "========================================================\n";

    // ==========================================
    // Phase 0: 串接 Consensus Blocks (加入 Debug 驗證)
    // ==========================================
    if (DEBUG_MODE) std::cout << "[Phase 0] Concatenating Involved Blocks into Super-Blocks...\n";
    
    auto refSuperBlock = refSet->concatenateBlocks(9999991); 
    auto qrySuperBlock = qrySet->concatenateBlocks(9999992);
    BlockSet refSuperSet ("ref_super");
    BlockSet qrySuperSet ("qry_super");
    refSuperBlock = refSuperSet.addBlock(refSuperBlock);
    qrySuperBlock = qrySuperSet.addBlock(qrySuperBlock);

    auto time1 = std::chrono::high_resolution_clock::now();

    // ==========================================
    // Phase 1: 使用 Dictionary + splitSingleBlock 切割
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 1] Extracting Cuts and Splitting SuperBlocks...\n";
    
    std::set<int> refCuts, qryCuts;
    std::vector<mga::Alignment> validAlignments; 

    for (auto& aln : alignments) {
        if (!aln.valid || (aln.type != mga::PRIMARY && aln.type != mga::SECONDARY)) continue;
        
        // 【新增條件】：如果 ref 端或 qry 端小於 100 base，直接丟棄這個 Alignment
        int rLen = std::abs(aln.refIdx.second - aln.refIdx.first);
        int qLen = std::abs(aln.qryIdx.second - aln.qryIdx.first);
        // if (rLen < 100 || qLen < 100) {
        //     aln.setValid2False();
        //     continue;
        // }

        refCuts.insert(aln.refIdx.first); refCuts.insert(aln.refIdx.second);
        qryCuts.insert(aln.qryIdx.first); qryCuts.insert(aln.qryIdx.second);
        
        validAlignments.push_back(aln); 
    }

    if (DEBUG_MODE) {
        std::cout << "Reference Cut Points: \n";
        for (auto cut: refCuts) std::cout << cut << " "; std::cout << "\n";
        std::cout << "Query Cut Points: \n";
        for (auto cut: qryCuts) std::cout << cut << " "; std::cout << "\n";
    }

    auto batchSplitSuperBlock = [&](BlockSet* bSet, std::shared_ptr<Block> superBlock, const std::set<int>& cuts) {
        auto start_t = std::chrono::high_resolution_clock::now();
        std::map<int, BlockSet::SegNode> finalDict;
        int superLen = superBlock->getConsensus().length();
        
        // 1. 過濾合法切點並轉為 Vector (保留 std::set 的升序特性)
        std::vector<int> validCuts;
        for (int c : cuts) {
            if (c > 0 && c < superLen) validCuts.push_back(c);
        }
        
        if (validCuts.empty()) {
            finalDict[0] = {0, superLen, superBlock->getId()};
            return finalDict;
        }

        // 2. 呼叫終極特化函數，一刀將 SuperBlock 切成 N+1 塊！
        std::vector<BlockID> newBlockIDs = bSet->splitMultiBlocks(superBlock->getId(), validCuts);

        // 防呆：如果因為某些原因切割失敗或沒切，退回原狀
        if (newBlockIDs.size() != validCuts.size() + 1) {
            if (DEBUG_MODE) std::cerr << "[ERROR] splitMultiBlock mismatch!\n";
            finalDict[0] = {0, superLen, superBlock->getId()};
            return finalDict;
        }
        
        // 3. 直接將結果映射回 Dict (因為 splitMultiBlock 保證是由左至右生成)
        int prevCut = 0;
        for (size_t i = 0; i < validCuts.size(); ++i) {
            finalDict[prevCut] = {prevCut, validCuts[i], newBlockIDs[i]};
            prevCut = validCuts[i];
        }
        // 補上最後一段尾巴
        finalDict[prevCut] = {prevCut, superLen, newBlockIDs.back()};
        
        auto end_t = std::chrono::high_resolution_clock::now();
        
        if (DEBUG_MODE) {
            // 改回 milliseconds，或者把字串改成 us
            auto totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count();
            // 由於被 parallel_invoke 呼叫，簡單印一行就好，避免嚴重的交錯
            std::cout << "[DEBUG] Split into " << newBlockIDs.size() << " blocks. Time: " << totalTime << " ms.\n";
        }
        
        return finalDict;
    };

    // 啟動雙管齊下的平行處理！
    std::map<int, BlockSet::SegNode> refDict, qryDict;
    tbb::parallel_invoke(
        [&] { refDict = batchSplitSuperBlock(&refSuperSet, refSuperBlock, refCuts); },
        [&] { qryDict = batchSplitSuperBlock(&qrySuperSet, qrySuperBlock, qryCuts); }
    );

    if (DEBUG_MODE) {
        std::cout << "  -> Ref SuperBlock split into " << refDict.size() << " atomic blocks.\n";
        std::cout << "  -> Qry SuperBlock split into " << qryDict.size() << " atomic blocks.\n";
    }

    auto time1_5 = std::chrono::high_resolution_clock::now();
    

    refSuperSet.debugValidateSegments(false);
    qrySuperSet.debugValidateSegments(false);

    auto time2 = std::chrono::high_resolution_clock::now();


    std::map<int, BlockID> refBlocksMap;
    for (auto& kv : refDict) {
        refBlocksMap[kv.first] = kv.second.blkId;
    }

    std::map<int, BlockID> qryBlocksMap;
    // 遍歷 qryDict，把切好的 qry block 逐一加進 refSuperSet
    for (auto& kv : qryDict) {
        int qryStart = kv.first;              // Qry 的 Genomic Coordinate
        BlockID oldQryId = kv.second.blkId; // 在 qrySuperSet 裡的舊 ID

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
    std::map<BlockID, std::shared_ptr<Block>> globalBlockPool;
    
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
        std::vector<BlockID> rBlocks;
        auto rIt = refDict.upper_bound(rMin); 
        if (rIt != refDict.begin()) rIt--;
        while (rIt != refDict.end() && rIt->second.start < rMax) {
            rBlocks.push_back(rIt->second.blkId);
            rIt++;
        }

        std::vector<BlockID> qBlocks;
        auto qIt = qryDict.upper_bound(qMin); 
        if (qIt != qryDict.begin()) qIt--;
        while (qIt != qryDict.end() && qIt->second.start < qMax) {
            qBlocks.push_back(qIt->second.blkId);
            qIt++;
        }

        // 把這個範圍內的所有 Ref 積木與 Qry 積木，全部拉進同一個 UnionFind 群組！
        for (BlockID rId : rBlocks) {
            for (BlockID qId : qBlocks) {
                uf.unite(rId, qId);
                matchCount++;
            }
        }
    }
    
    if (DEBUG_MODE) std::cout << "  -> Grouped " << matchCount << " homologous pairs.\n";



    auto time3 = std::chrono::high_resolution_clock::now();

    // ==========================================
    // 建立供 Phase 4 提取 CIGAR 使用的座標與來源查詢表
    // ==========================================
    std::map<BlockID, std::pair<int, int>> blockConsensusCoords;
    std::map<BlockID, bool> isRefBlockMap;

    for (const auto& kv : refDict) {
        blockConsensusCoords[kv.second.blkId] = {kv.second.start, kv.second.end};
        isRefBlockMap[kv.second.blkId] = true; 
    }
    for (const auto& kv : qryDict) {
        blockConsensusCoords[kv.second.blkId] = {kv.second.start, kv.second.end};
        isRefBlockMap[kv.second.blkId] = false; 
    }

    // ==========================================
    // Phase 4: Iterative Spanning-Tree Merging per Group
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n[Phase 4] Iterative Spanning-Tree Merging within Groups...\n";
    
    std::string newID = "Merged_" + refSet->getId() + "_" + qrySet->getId();
    BlockSet* resultSet = createBlockSet(newID); 

    std::map<BlockID, std::vector<BlockID>> groupedBlocks;
    for (auto& kv : uf.parent) {
        groupedBlocks[uf.find(kv.first)].push_back(kv.first);
    }

    if (DEBUG_MODE) std::cout << "  -> Total unique groups to process: " << groupedBlocks.size() << "\n\n";

    std::set<std::shared_ptr<Block>> workingPool;
    std::map<BlockID, std::shared_ptr<Block>> oldToNewBlockMap;
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
        BlockID origBaseId = members[0];
        for (BlockID id : members) {
            if (isRefBlockMap[id]) { origBaseId = id; break; }
        }

        std::shared_ptr<Block> baseBlock = globalBlockPool[origBaseId];

        if (DEBUG_MODE) std::cout << "  ├─ Initial Hub Block ID: " << origBaseId 
                                  << " (Initial Len: " << baseBlock->getConsensus().length() << " bp)\n";

        // ==========================================
        // 2. 初始化群組的 Spanning Tree 狀態
        // ==========================================
        std::map<BlockID, std::vector<int>> coordMaps; // 記錄每一個 Original Member 的座標變化
        coordMaps[origBaseId] = std::vector<int>(baseBlock->getConsensus().length() + 1);
        std::iota(coordMaps[origBaseId].begin(), coordMaps[origBaseId].end(), 0);

        std::set<BlockID> mergedMembers = {origBaseId};
        std::set<BlockID> unmergedMembers;
        for (BlockID id : members) {
            if (id != origBaseId) unmergedMembers.insert(id);
        }

        // 【新增】：追蹤每個積木加入 Hub 時的「絕對反轉狀態」
        std::map<BlockID, bool> isReversedInHub;
        isReversedInHub[origBaseId] = false;

        // 3. 核心迴圈：利用圖的邊緣 (Valid Alignments) 依序拉攏未合併的積木
        while (!unmergedMembers.empty()) {
            bool foundEdge = false;
            BlockID targetMergedId = 0;
            BlockID targetUnmergedId = 0;
            mga::Cigar bestOrigCigar;
            bool bestInverse = false;

            // 尋找任一個 "已合併積木" 與 "未合併積木" 之間的合法 Alignment
            for (BlockID u : unmergedMembers) {
                for (BlockID m : mergedMembers) {
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

            // ==========================================
            // A. 處理 Inverse 邏輯與 CIGAR 方向性 (負正得負)
            // ==========================================
            // 1. 取得橋樑積木 (M) 當初加入 Hub 時的方向狀態
            bool mIsRev = isReversedInHub[targetMergedId];
            
            // 2. 計算新積木 (U) 應該套用的真實反轉狀態
            bool effectiveInverse = mIsRev ^ bestInverse;

            // 3. 【核心修復】：原本的 bestOrigCigar 是基於原始 M 到原始 U。
            // 如果 M 已經在 Hub 中反轉，我們必須將 CIGAR 左右對調，才能讓 coordMaps 正確映射！
            mga::Cigar inputCigar = bestOrigCigar;
            if (mIsRev) {
                std::reverse(inputCigar.begin(), inputCigar.end());
            }

            // 4. 將 CIGAR 補償為指向「當前 SuperBase 的 Consensus」
            int superBaseLen = baseBlock->getConsensus().length();
            mga::Cigar adjustedCigar = adjustCigarWithMap(inputCigar, coordMaps[targetMergedId], superBaseLen);

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

            // D. 執行真實的物理合併 (傳入算好的 effectiveInverse)
            baseBlock = refSuperSet.mergeTwoBlocks(baseBlock, memberBlock, adjustedCigar, effectiveInverse); 

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
            isReversedInHub[targetUnmergedId] = effectiveInverse; // 【新增】：記錄這塊積木的最終反轉狀態
        }

        if (DEBUG_MODE) std::cout << "  └─ Final Merged Block ID: " << baseBlock->getId() << "\n\n";

        workingPool.insert(baseBlock);
        for (BlockID memberId : mergedMembers) {
            oldToNewBlockMap[memberId] = baseBlock;
        }
    }

    auto time4 = std::chrono::high_resolution_clock::now();

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

    auto time5 = std::chrono::high_resolution_clock::now();

    if (DEBUG_MODE) {
        std::cout << "Time 1  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << " ms\n";
        std::cout << "Time 1.5: " << std::chrono::duration_cast<std::chrono::milliseconds>(time1_5 - time1).count() << " ms\n";
        std::cout << "Time 2  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1_5).count() << " ms\n";
        std::cout << "Time 3  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << " ms\n";
        std::cout << "Time 4  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time4 - time3).count() << " ms\n";
        std::cout << "Time 5  : " << std::chrono::duration_cast<std::chrono::milliseconds>(time5 - time4).count() << " ms\n";
    }


    if (DEBUG_MODE) std::cout << "========================================================\n"
                              << "=== GRAPH MERGE COMPLETED SUCCESSFULLY ===\n"
                              << "========================================================\n\n";


                

    return resultSet;
}



/*
Alignments splitAlignmentsByCuts( const Alignments& alignments, const std::set<int>& refCuts, const std::set<int>& qryCuts) 
{
    Alignments newAlignments;

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
*/