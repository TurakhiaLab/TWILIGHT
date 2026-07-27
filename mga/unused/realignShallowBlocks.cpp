#include "mga.hpp"

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <numeric>

void BlockSet::realignShallowBlocks(std::string tempDir) {
    bool DEBUG_MODE = true;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet: Realigning & Absorbing Shallow Blocks ===\n"
                              << "============================================================\n";

    const int DEEP_THRESHOLD = 3;
    
    std::string refFasta = tempDir + "/" + this->id_ + "_shallow_ref.fa";
    std::string qryFasta = tempDir + "/" + this->id_ + "_shallow_qry.fa";
    std::string pafFile  = tempDir + "/" + this->id_ + "_shallow_aln.paf";

    std::ofstream refOut(refFasta);
    std::ofstream qryOut(qryFasta);

    int refCount = 0, qryCount = 0;

    // 1. 分類 Blocks 並輸出 FASTA
    for (auto blk : this->getAllBlocks()) {
        if (!blk || blk->getSequences().empty()) continue;
        int seqCount = blk->getSequences().size();
        if (seqCount > DEEP_THRESHOLD) {
            refOut << ">" << blk->getId() << "\n" << blk->getConsensus() << "\n";
            refCount++;
        } else {
            qryOut << ">" << blk->getId() << "\n" << blk->getConsensus() << "\n";
            qryCount++;
        }
    }
    refOut.close(); qryOut.close();

    if (refCount == 0 || qryCount == 0) {
        std::cout << "No shallow blocks!\n";
        return;
    }

    std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
   

    // 2. 執行 minimap2 進行對齊 (-c 取得 cg:Z CIGAR)
    // std::string cmd = minimap2_path + " -cx asm5 " + refFasta + " " + qryFasta + " > " + pafFile + " 2> /dev/null";
    std::string cmd = minimap2_path + " -cx asm20 " + refFasta + " " + qryFasta + " > " + pafFile;
    
    if (std::system(cmd.c_str()) != 0) {
        std::cerr << "[ERROR] minimap2 execution failed!\n";
        return;
    }

    // 3. 解析 PAF (使用你寫好的 parser)，挑選最高分 Alignment
    mga::alnVec parsedAlignments = mga::parser::parseMinimap2PAF(pafFile);
    std::unordered_map<Block::ID, mga::Alignment> bestAlignments;
    
    for (const auto& aln : parsedAlignments) {
        Block::ID qId = std::stoull(aln.qryName);
        if (bestAlignments.find(qId) == bestAlignments.end() || aln.alnScore > bestAlignments[qId].alnScore) {
            bestAlignments[qId] = aln;
        }
    }

    std::map<Block::ID, std::vector<Block::ID>> refToQrys;
    for (const auto& kv : bestAlignments) {
        Block::ID rId = std::stoull(kv.second.refName);
        refToQrys[rId].push_back(kv.first);
    }

    // ==========================================
    // 4. 強制合併 (完美還原你的 Phase 4 核心邏輯)
    // ==========================================
    for (auto& group : refToQrys) {
        Block::ID refId = group.first;
        auto baseBlk = this->getBlock(refId);
        if (!baseBlk) continue;

        if (DEBUG_MODE) std::cout << "\n[Target Ref] Block " << refId << " absorbing " << group.second.size() << " blocks.\n";

        // 【核心機制】：初始化座標映射表 (Tracking Map)
        int origRefLen = baseBlk->getConsensus().length();
        std::vector<int> coordMap(origRefLen + 1);
        std::iota(coordMap.begin(), coordMap.end(), 0);

        for (Block::ID qryId : group.second) {
            auto qryBlk = this->getBlock(qryId);
            if (!qryBlk) continue;
            auto& aln = bestAlignments[qryId];

            // A. Padding: 把 Local Alignment 補成 Full-length CIGAR
            int padRFront = aln.refIdx.first;
            int padRTail  = origRefLen - aln.refIdx.second;
            int padQFront, padQTail;

            int qBlkLen = qryBlk->getConsensus().length();
            if (!aln.inverse) {
                padQFront = aln.qryIdx.first;
                padQTail  = qBlkLen - aln.qryIdx.second;
            } else {
                // PAF 反向鏈座標轉換陷阱：CIGAR 是配對 Reverse Complement Qry
                padQFront = qBlkLen - aln.qryIdx.second;
                padQTail  = aln.qryIdx.first;
            }

            mga::Cigar paddedCigar;
            if (padRFront > 0) paddedCigar.push_back({padRFront, 'D'});
            if (padQFront > 0) paddedCigar.push_back({padQFront, 'I'});
            paddedCigar.insert(paddedCigar.end(), aln.CIGAR.begin(), aln.CIGAR.end());
            if (padRTail > 0) paddedCigar.push_back({padRTail, 'D'});
            if (padQTail > 0) paddedCigar.push_back({padQTail, 'I'});
            paddedCigar = mga::compressCigar(paddedCigar);

            // B. 【核心機制】：依據 coordMap，將 CIGAR 補償為「當前變胖後的 Consensus」
            mga::Cigar adjustedCigar;
            int origR = 0, currR = 0;
            
            for (auto op : paddedCigar) {
                char t = op.second; int l = op.first;
                if (t == 'S' || t == 'H') t = 'I';

                if (t == 'M' || t == '=' || t == 'X' || t == 'D') {
                    for (int i = 0; i < l; ++i) {
                        int targetR = (origR < coordMap.size()) ? coordMap[origR] : currR;
                        // 如果別人在這裡先 Insert 了，我們就必須補上 D，跳過那些新鹼基！
                        while (currR < targetR) {
                            adjustedCigar.push_back({1, 'D'});
                            currR++;
                        }
                        adjustedCigar.push_back({1, t});
                        origR++; currR++;
                    }
                } else if (t == 'I') {
                    adjustedCigar.push_back({l, 'I'});
                }
            }
            int finalTargetR = coordMap.back();
            while (currR < finalTargetR) {
                adjustedCigar.push_back({1, 'D'});
                currR++;
            }
            adjustedCigar = mga::compressCigar(adjustedCigar);

            // C. 執行真實的物理合併 (由於 Ref 永遠是正向，直接傳入 aln.inverse 即可)
            baseBlk = this->mergeTwoBlocks(baseBlk, qryBlk, adjustedCigar, aln.inverse);

            // D. 【核心機制】：準備 stepCoordMap，計算這次合併造成的 Consensus 長胖
            int adjRefLen = 0;
            for (auto op : adjustedCigar) {
                if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') adjRefLen += op.first;
            }
            std::vector<int> stepCoordMap(adjRefLen + 1, 0);
            int tmpR = 0, mPos = 0;
            for (auto op : adjustedCigar) {
                char t = op.second; int l = op.first;
                if (t == 'M' || t == '=' || t == 'X' || t == 'D') {
                    for (int i = 0; i < l; ++i) { if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR++] = mPos++; }
                } else if (t == 'I') {
                    mPos += l; // Insertion 會撐大新的 Consensus
                }
            }
            if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR] = mPos;

            // E. 【核心機制】：套用聯動更新
            for (size_t i = 0; i < coordMap.size(); ++i) {
                if (coordMap[i] < stepCoordMap.size()) {
                    coordMap[i] = stepCoordMap[coordMap[i]];
                }
            }
            
            this->deleteBlock(qryId); // 清除已吸收的積木
        }
    }

    // ==========================================
    // 5. 掃尾：重建圖譜拓撲與清垃圾
    // ==========================================
    this->rebuildAllPointers();

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    if (DEBUG_MODE) std::cout << "=== Realign Shallow Blocks: Completed ===\n\n";
}

void BlockSet::realignAllToAll(std::string tempDir) {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet: All-to-All Realignment & Absorption ===\n"
                              << "============================================================\n";

    std::string allFasta = tempDir + "/" + this->id_ + "_all_blocks.fa";
    std::string pafFile  = tempDir + "/" + this->id_ + "_all_aln.paf";

    std::ofstream allOut(allFasta);
    int blockCount = 0;

    // 1. 輸出所有 Block 到 FASTA
    for (auto blk : this->getAllBlocks()) {
        if (!blk || blk->getSequences().empty()) continue;
        allOut << ">" << blk->getId() << "\n" << blk->getConsensus() << "\n";
        blockCount++;
    }
    allOut.close();

    if (blockCount < 2) {
        if (DEBUG_MODE) std::cout << "[INFO] Not enough blocks for All-to-All realignment. Skipping.\n";
        return;
    }

    if (DEBUG_MODE) std::cout << "  -> Identified " << blockCount << " Blocks for All-to-All Battle Royale.\n";

    // 2. 執行 minimap2 進行 All-to-All 對齊 (-X 避免 self-hits 和 dual-hits)
    if (DEBUG_MODE) std::cout << "  -> Running minimap2...\n";
    std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
   

    // 2. 執行 minimap2 進行對齊 (-c 取得 cg:Z CIGAR)
    // std::string cmd = minimap2_path + " -cx asm5 " + refFasta + " " + qryFasta + " > " + pafFile + " 2> /dev/null";
    std::string cmd = minimap2_path + " -cx asm5 -s 50 " + allFasta + " " + allFasta + " > " + pafFile + " 2> /dev/null"; 
    
    // std::string cmd = "minimap2 -X -x asm5 -c " + allFasta + " " + allFasta + " > " + pafFile + " 2> /dev/null";
    if (std::system(cmd.c_str()) != 0) {
        std::cerr << "[ERROR] minimap2 execution failed!\n";
        return;
    }

    // 3. 翻轉 Alignment 工具 (強制讓 大積木=Ref, 小積木=Qry)
    auto flipAlignment = [](mga::Alignment aln) -> mga::Alignment {
        std::swap(aln.refName, aln.qryName);
        std::swap(aln.refIdx, aln.qryIdx);
        for (auto& op : aln.CIGAR) {
            if (op.second == 'D') op.second = 'I';
            else if (op.second == 'I') op.second = 'D';
        }
        return aln;
    };

    // 4. 解析 PAF 並建立階級制度 (Hierarchy)
    mga::alnVec parsedAlignments = mga::parser::parseMinimap2PAF(pafFile);
    std::unordered_map<Block::ID, mga::Alignment> bestAlignments;
    
    for (auto aln : parsedAlignments) {
        Block::ID idA = std::stoull(aln.qryName);
        Block::ID idB = std::stoull(aln.refName);
        if (idA == idB) continue; // 防呆：略過自我比對

        auto blkA = this->getBlock(idA);
        auto blkB = this->getBlock(idB);
        if (!blkA || !blkB) continue;

        // 定義大小 (Score = 序列數量)
        int sizeA = blkA->getSequences().size();
        int sizeB = blkB->getSequences().size();

        Block::ID hubId, spokeId;
        
        // 強制：序列多的人當 Hub (Ref)，序列少的當 Spoke (Qry)
        // 如果數量一樣，ID 小的當 Hub (保證穩定性)
        if (sizeB > sizeA || (sizeB == sizeA && idB < idA)) {
            // PAF 剛好是 B=Ref, A=Qry，不需翻轉
            hubId = idB; spokeId = idA;
        } else {
            // PAF 是 A=Qry, B=Ref，但 A 比較大！所以要翻轉 Alignment
            hubId = idA; spokeId = idB;
            aln = flipAlignment(aln);
        }

        // 每個 Spoke 只挑選一個最高分、最穩的 Hub 投靠
        if (bestAlignments.find(spokeId) == bestAlignments.end() || aln.alnScore > bestAlignments[spokeId].alnScore) {
            bestAlignments[spokeId] = aln;
        }
    }

    // 5. 【防崩潰機制】：強制星狀拓撲 (Enforce Star Topology)
    std::set<Block::ID> activeSpokes;
    for (const auto& kv : bestAlignments) activeSpokes.insert(kv.first);

    std::map<Block::ID, std::vector<Block::ID>> hubToSpokes;
    int droppedChains = 0;

    for (const auto& kv : bestAlignments) {
        Block::ID spoke = kv.first;
        Block::ID hub = std::stoull(kv.second.refName);

        // 如果 Hub 自己也準備被別人吃掉，那就剝奪它吃別人的權力 (避免 Chain Reaction)
        if (activeSpokes.count(hub)) {
            droppedChains++;
            continue; // 這個 Spoke 這回合先放生，等 Hub 被吃完後，下一回合它自然會直接對齊到最終的老大
        }
        hubToSpokes[hub].push_back(spoke);
    }

    // if (DEBUG_MODE) {
    if (true) {
        std::cout << "  -> Found " << bestAlignments.size() << " valid absorption candidates.\n";
        if (droppedChains > 0) std::cout << "  -> Dropped " << droppedChains << " chained interactions to strictly enforce Star Topology.\n";
    }

    // ==========================================
    // 6. 強制合併 (含動態 CIGAR 座標校正機制)
    // ==========================================
    for (auto& group : hubToSpokes) {
        Block::ID refId = group.first;
        auto baseBlk = this->getBlock(refId);
        if (!baseBlk) continue;

        if (DEBUG_MODE) std::cout << "\n[Hub] Block " << refId << " absorbing " << group.second.size() << " spokes.\n";

        int origRefLen = baseBlk->getConsensus().length();
        std::vector<int> coordMap(origRefLen + 1);
        std::iota(coordMap.begin(), coordMap.end(), 0);

        for (Block::ID qryId : group.second) {
            auto qryBlk = this->getBlock(qryId);
            if (!qryBlk) continue;
            auto& aln = bestAlignments[qryId];

            int padRFront = aln.refIdx.first;
            int padRTail  = origRefLen - aln.refIdx.second;
            int padQFront, padQTail;

            int qBlkLen = qryBlk->getConsensus().length();
            if (!aln.inverse) {
                padQFront = aln.qryIdx.first;
                padQTail  = qBlkLen - aln.qryIdx.second;
            } else {
                padQFront = qBlkLen - aln.qryIdx.second;
                padQTail  = aln.qryIdx.first;
            }

            mga::Cigar paddedCigar;
            if (padRFront > 0) paddedCigar.push_back({padRFront, 'D'});
            if (padQFront > 0) paddedCigar.push_back({padQFront, 'I'});
            paddedCigar.insert(paddedCigar.end(), aln.CIGAR.begin(), aln.CIGAR.end());
            if (padRTail > 0) paddedCigar.push_back({padRTail, 'D'});
            if (padQTail > 0) paddedCigar.push_back({padQTail, 'I'});
            paddedCigar = mga::compressCigar(paddedCigar);

            mga::Cigar adjustedCigar;
            int origR = 0, currR = 0;
            
            for (auto op : paddedCigar) {
                char t = op.second; int l = op.first;
                if (t == 'S' || t == 'H') t = 'I';

                if (t == 'M' || t == '=' || t == 'X' || t == 'D') {
                    for (int i = 0; i < l; ++i) {
                        int targetR = (origR < coordMap.size()) ? coordMap[origR] : currR;
                        while (currR < targetR) {
                            adjustedCigar.push_back({1, 'D'});
                            currR++;
                        }
                        adjustedCigar.push_back({1, t});
                        origR++; currR++;
                    }
                } else if (t == 'I') {
                    adjustedCigar.push_back({l, 'I'});
                }
            }
            int finalTargetR = coordMap.back();
            while (currR < finalTargetR) {
                adjustedCigar.push_back({1, 'D'});
                currR++;
            }
            adjustedCigar = mga::compressCigar(adjustedCigar);

            if (DEBUG_MODE) std::cout << "  ├─ Absorbing Spoke " << qryId << " (Score: " << aln.alnScore << ")... ";
            
            Block::ID oldBaseId = baseBlk->getId(); 
            baseBlk = this->mergeTwoBlocks(baseBlk, qryBlk, adjustedCigar, aln.inverse);
            this->deleteBlock(oldBaseId);
            this->deleteBlock(qryId);
            if (DEBUG_MODE) std::cout << "New Block ID: " << baseBlk->getId() << "... Done.\n";

            int adjRefLen = 0;
            for (auto op : adjustedCigar) {
                if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') adjRefLen += op.first;
            }
            std::vector<int> stepCoordMap(adjRefLen + 1, 0);
            int tmpR = 0, mPos = 0;
            for (auto op : adjustedCigar) {
                char t = op.second; int l = op.first;
                if (t == 'M' || t == '=' || t == 'X' || t == 'D') {
                    for (int i = 0; i < l; ++i) { if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR++] = mPos++; }
                } else if (t == 'I') {
                    mPos += l; 
                }
            }
            if (tmpR < stepCoordMap.size()) stepCoordMap[tmpR] = mPos;

            for (size_t i = 0; i < coordMap.size(); ++i) {
                if (coordMap[i] < stepCoordMap.size()) {
                    coordMap[i] = stepCoordMap[coordMap[i]];
                }
            }
            
            this->deleteBlock(qryId);
        }
    }

    // ==========================================
    // 7. 掃尾：重建圖譜拓撲與清垃圾
    // ==========================================
    if (DEBUG_MODE) std::cout << "\n  -> Rebuilding graph pointers...\n";
    this->rebuildAllPointers();

    std::vector<Block::ID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk->getSequences().empty()) emptyBlocks.push_back(blk->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    std::remove(allFasta.c_str());
    std::remove(pafFile.c_str());

    if (DEBUG_MODE) std::cout << "=== All-to-All Realignment: Completed ===\n\n";
}