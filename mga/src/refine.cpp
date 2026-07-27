#include "block_manager.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

/*
void BlockSet::refineGraph() {
    std::cout << "\n[Pipeline] Starting Graph Refinement...\n";
    // 1. Zip identical branch
    auto time0 = std::chrono::high_resolution_clock::now();

    this->splitBlocksByLongGaps();

    auto time1 = std::chrono::high_resolution_clock::now();

    this->extractMicroSegments();

    auto time2 = std::chrono::high_resolution_clock::now();

   // this->absorbMicroBlocksIndividual();

    auto time3 = std::chrono::high_resolution_clock::now();

    std::cout << "refine 1: " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << "ms\n";
    std::cout << "refine 2: " << std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() << "ms\n";
    std::cout << "refine 3: " << std::chrono::duration_cast<std::chrono::milliseconds>(time3 - time2).count() << "ms\n";

    this->cleanGappedColumns();
    this->rebuildAllPointers();
    

    std::cout << "[Pipeline] Graph Refinement Completed!\n";
}


void BlockSet::cleanGappedColumns() {
    bool DEBUG_MODE = false;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << "=== BlockSet: Cleaning All-Gap Columns (Consensus Shrink) ===\n"
                              << "============================================================\n";

    auto allBlocks = this->getAllBlocks();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t idx = r.begin(); idx != r.end(); ++idx) {
                auto blk = allBlocks[idx].lock();
                if (!blk) continue;

                int consLen = blk->getConsensus().length();
                if (consLen == 0) continue;

                // ==========================================
                // 1. 精準計數法：計算每個 Column 出現過幾次 GAP
                // ==========================================
                std::vector<int> gapCount(consLen, 0);
                int numSegments = 0;

                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPairInner : seqPair.second.getSegments()) {
                        numSegments++;
                        Segment& seg = segPairInner.second;
                        
                        for (auto& var : seg.getVariants()) {
                            if (var.getType() == VariantType::GAP) {
                                for (int i = var.getStart(); i < var.getEnd(); ++i) {
                                    if (i < consLen) gapCount[i]++;
                                }
                            }
                        }
                    }
                }

                if (numSegments == 0) continue; // 空積木不處理

                // ==========================================
                // 2. 判斷是否為純空氣欄位 (所有 Segment 都在這有 Gap)
                // ==========================================
                std::vector<bool> hasBase(consLen, true);
                bool needsCleaning = false;
                
                for (int i = 0; i < consLen; ++i) {
                    if (gapCount[i] == numSegments) {
                        hasBase[i] = false; // 確定是 All-Gap Column
                        needsCleaning = true;
                    }
                }

                if (!needsCleaning) continue; 

                // ==========================================
                // 3. 建立映射表並縮短 Consensus
                // ==========================================
                std::vector<int> oldToNew(consLen + 1, 0);
                std::string oldCons = blk->getConsensus();
                std::string newCons = "";
                newCons.reserve(consLen);

                int newPos = 0;
                for (int i = 0; i < consLen; ++i) {
                    oldToNew[i] = newPos;
                    if (hasBase[i]) {
                        newCons += oldCons[i];
                        newPos++;
                    }
                }
                oldToNew[consLen] = newPos;

                blk->setConsensus(newCons);

                // ==========================================
                // 4. 更新 Segment 的 Variations (相對座標平移)
                // ==========================================
                for (auto& seqPair : blk->getSequences()) {
                    for (auto& segPairInner : seqPair.second.getSegments()) {
                        Segment& seg = segPairInner.second;
                        
                        std::vector<Variant> newVars;
                        for (auto& var : seg.getVariants()) {
                            if (var.getType() == VariantType::GAP) {
                                int ns = oldToNew[var.getStart()];
                                int ne = oldToNew[var.getEnd()];
                                if (ne > ns) {
                                    newVars.push_back(Variant::createGap(ns, ne));
                                }
                            } else if (var.getType() == VariantType::SNV) {
                                int pos = var.getStart();
                                if (hasBase[pos]) { 
                                    newVars.push_back(Variant(oldToNew[pos], var.getAlt()));
                                }
                            }
                        }
                        seg.getVariants() = std::move(newVars);
                    }
                }
            }
        }
    );

    if (DEBUG_MODE) std::cout << "=== Clean Gapped Columns Completed ===\n\n";
}



struct GapInterval {
    int start;
    int end;
    bool operator<(const GapInterval& other) const {
        if (start != other.start) return start < other.start;
        return end < other.end;
    }
};


void BlockSet::splitBlocksByLongGaps() {
    const int min_gap_length = 50;
    const int merge_distance = 10; 
    bool debug = false;
    
    // --- [Debug 工具保持不變] ---
    auto printBlockMSA = [&](std::shared_ptr<Block> b, const std::string& title) {
        if (!b) return;
        std::string cons = b->getConsensus();
        std::cout << "\n>>> " << title << " (Block ID: " << b->getId() << ", Length: " << cons.length() << ") <<<\n";
        std::cout << "Consensus :\t" << cons << "\n";
        
        for (auto& seq_pair : b->getSequences()) {
            const std::string& seq_name = seq_pair.first;
            for (auto& seg_pair : seq_pair.second.getSegments()) {
                Segment& seg = seg_pair.second;
                std::string aligned_seq = cons; 
                for (Variant& var : seg.getVariants()) {
                    if (var.getType() == VariantType::GAP) {
                        for (int i = var.getStart(); i < var.getEnd() && i < aligned_seq.length(); ++i) {
                            aligned_seq[i] = '-';
                        }
                    } else if (var.getType() == VariantType::SNV) {
                        if (var.getStart() < aligned_seq.length()) {
                            aligned_seq[var.getStart()] = var.getAlt();
                        }
                    }
                }
                std::cout << seq_name << " :\t" << aligned_seq << "\n";
            }
        }
        std::cout << "--------------------------------------------------\n";
    };
    // ---------------------------------------------------------

    // 取得所有目標 Block ID
    std::vector<BlockID> target_blocks;
    target_blocks.reserve(blocks.size());
    for (const auto& pair : blocks) {
        target_blocks.push_back(pair.first);
    }

    // 定義任務結構：記錄「哪個 Block」要在「哪些位置」切斷
    struct SplitTask {
        BlockID blkId;
        std::vector<int> cuts;
    };
    
    // 無鎖安全陣列，用來接收多執行緒探勘到的切斷任務
    tbb::concurrent_vector<SplitTask> pendingSplits;

    // ==========================================
    // 階段 1: [平行化] 尋找與計算最佳切點
    // ==========================================
    tbb::parallel_for(tbb::blocked_range<size_t>(0, target_blocks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                BlockID blkId = target_blocks[i];
                auto block = getBlock(blkId);
                if (!block) continue;

                int consensus_len = block->getConsensus().length();
                
                // 💡 Early Rejection: 長度不夠直接跳過，省去巨量迴圈！
                if (consensus_len < min_gap_length) continue;

                std::vector<int> raw_cut_points;

                // 收集所有符合條件的 Gap 邊界
                for (auto& seq_pair : block->getSequences()) {
                    for (auto& seg_pair : seq_pair.second.getSegments()) {
                        for (Variant& var : seg_pair.second.getVariants()) {
                            if (var.getType() == VariantType::GAP) {
                                int gap_length = var.getEnd() - var.getStart();
                                if (gap_length >= min_gap_length) {
                                    raw_cut_points.push_back(var.getStart());
                                    raw_cut_points.push_back(var.getEnd());
                                }
                            }
                        }
                    }
                }

                if (raw_cut_points.empty()) continue;

                // 排序與一維分群
                std::sort(raw_cut_points.begin(), raw_cut_points.end());
                std::set<int> final_cut_points;
                std::vector<int> current_cluster;
                current_cluster.push_back(raw_cut_points[0]);

                for (size_t j = 1; j < raw_cut_points.size(); ++j) {
                    if (raw_cut_points[j] - current_cluster.back() <= merge_distance) {
                        current_cluster.push_back(raw_cut_points[j]);
                    } else {
                        int rep_cut = current_cluster[current_cluster.size() / 2];
                        if (rep_cut > 0 && rep_cut < consensus_len) {
                            final_cut_points.insert(rep_cut);
                        }
                        current_cluster.clear();
                        current_cluster.push_back(raw_cut_points[j]);
                    }
                }

                if (!current_cluster.empty()) {
                    int rep_cut = current_cluster[current_cluster.size() / 2];
                    if (rep_cut > 0 && rep_cut < consensus_len) {
                        final_cut_points.insert(rep_cut);
                    }
                }

                // 如果有算到切點，打包成任務丟進安全陣列
                if (!final_cut_points.empty()) {
                    std::vector<int> cuts(final_cut_points.begin(), final_cut_points.end());
                    pendingSplits.push_back({blkId, cuts});
                }
            }
        }
    );

    // ==========================================
    // 階段 2: [循序] 執行安全的結構切斷
    // ==========================================
    auto time0 = std::chrono::high_resolution_clock::now();
    std::cout << "Tasks: " << pendingSplits.size() << "\n";
    for (const auto& task : pendingSplits) {
        if (debug) {
            auto block = getBlock(task.blkId);
            if (block && block->getConsensus().size() < 1000) {
                std::cout << "\n==================================================\n";
                std::cout << "[DEBUG] Preparing to split Block " << task.blkId << "\n";
                std::cout << "Target Cut Points: ";
                for (int c : task.cuts) std::cout << c << " ";
                std::cout << "\n";
                printBlockMSA(block, "BEFORE SPLIT");
            }
        }

        // 執行實際切斷動作
        std::vector<BlockID> new_block_ids = splitMultiBlocks(task.blkId, task.cuts);

        if (debug) {
            std::cout << "[DEBUG] Split successful! Generated " << new_block_ids.size() << " new blocks.\n";
            for (BlockID new_id : new_block_ids) {
                auto new_block = getBlock(new_id);
                if (new_block && new_block->getConsensus().size() < 1000) {
                    printBlockMSA(new_block, "AFTER SPLIT (New Block)");
                }
            }
            std::cout << "==================================================\n";
        }
    }
    auto time1 = std::chrono::high_resolution_clock::now();
    std::cout << "Split Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() << "ms\n";
    // 只有在真的有切斷發生時，才需要花時間重建指標
    if (!pendingSplits.empty()) {
        this->rebuildAllPointers();
    }
}


void BlockSet::extractMicroSegments() {
    bool DEBUG_MODE = true;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << ">>> Executing extractMicroSegments (Parallel <= 15bp isolation)...\n"
                              << "============================================================\n";

    const int MICRO_THRESHOLD = 15;
    std::atomic<int> extracted_count{0};
    std::atomic<int> deleted_zeros{0};

    auto getRawSeq = [](const std::string& cons, Segment& seg) {
        std::string raw = "";
        int pos = 0;
        for (auto& var : seg.getVariants()) {
            if (var.getStart() > pos) raw += cons.substr(pos, var.getStart() - pos);
            if (var.getType() == VariantType::SNV) {
                raw += var.getAlt();
                pos = var.getStart() + 1;
            } else if (var.getType() == VariantType::GAP) {
                pos = var.getEnd();
            }
        }
        if (pos < cons.length()) raw += cons.substr(pos);
        return raw;
    };

    // 建立一個結構來暫存需要建立的新 Block 資訊
    struct MicroTask {
        std::string rawSeq;
        std::string seqID;
        Segment cleanSeg;
    };

    // 無鎖安全陣列，所有 Thread 都可以安全地把拔下來的資料丟進來
    tbb::concurrent_vector<MicroTask> pendingCreations;

    auto allBlocks = this->getAllBlocks();

    // ==========================================
    // 階段 1: [平行化] 尋找、拔除並暫存微小 Segment
    // ==========================================
    tbb::parallel_for(tbb::blocked_range<size_t>(0, allBlocks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto blk = allBlocks[i].lock();
                if (!blk) continue;
                if (blk->getConsensus().size() <= MICRO_THRESHOLD) continue;

                std::vector<std::pair<std::string, int>> segsToRemove;
                int local_extracted = 0;
                int local_deleted = 0;

                // 尋找要拔除的 Segment
                for (auto& seqPair : blk->getSequences()) {
                    std::string seqID = seqPair.first;
                    for (auto& segPairInner : seqPair.second.getSegments()) {
                        Segment& seg = segPairInner.second;
                        int actualLen = std::abs(seg.getEnd() - seg.getStart());

                        if (actualLen <= MICRO_THRESHOLD) {
                            segsToRemove.push_back({seqID, seg.getStart()});

                            if (actualLen > 0) {
                                std::string rawSeq = getRawSeq(blk->getConsensus(), seg);
                                Segment cleanSeg = seg;
                                cleanSeg.getVariants().clear(); 
                                
                                // 把建立新 Block 需要的素材丟進無鎖陣列
                                pendingCreations.push_back({rawSeq, seqID, cleanSeg});
                                local_extracted++;
                            } else {
                                local_deleted++;
                            }
                        }
                    }
                }

                // 執行安全移除 (每個 Thread 只改自己的 blk，絕對安全)
                for (auto& rmPair : segsToRemove) {
                    auto it = blk->getSequences().find(rmPair.first);
                    if (it != blk->getSequences().end()) {
                        it->second.getSegments().erase(rmPair.second);
                        if (it->second.getSegments().empty()) {
                            blk->getSequences().erase(it);
                        }
                    }
                }

                // 更新全域計數器
                if (local_extracted > 0) extracted_count += local_extracted;
                if (local_deleted > 0) deleted_zeros += local_deleted;
            }
        }
    );

    // ==========================================
    // 階段 1.5: [循序] 安全地建立新 Block
    // (這段非常快，因為不涉及字串比對，只是單純 new Object)
    // ==========================================
    for (auto& task : pendingCreations) {
        auto isolatedBlk = this->createBlock(task.rawSeq);
        Sequence newSeqInfo(task.seqID);
        newSeqInfo.getSegments()[task.cleanSeg.getStart()] = task.cleanSeg;
        isolatedBlk->addSequence(newSeqInfo);
    }

    // ==========================================
    // 階段 2: 刪除被徹底掏空的 Block
    // ==========================================
    std::vector<BlockID> emptyBlocks;
    for (auto blk : this->getAllBlocks()) {
        if (blk.lock()->getSequences().empty()) emptyBlocks.push_back(blk.lock()->getId());
    }
    for (auto id : emptyBlocks) this->deleteBlock(id);

    if (DEBUG_MODE) {
        std::cout << "  -> Extracted " << extracted_count << " micro segments (1-" << MICRO_THRESHOLD << "bp) into isolated blocks.\n";
        std::cout << "  -> Purged " << deleted_zeros << " pure gap (0-char) segments.\n";
    }

    // ==========================================
    // 階段 3: 自動擠壓與圖譜重建
    // ==========================================
    if (extracted_count > 0 || deleted_zeros > 0) {
        if (DEBUG_MODE) std::cout << "  -> Triggering cleanGappedColumns() to shrink host blocks...\n";
        this->cleanGappedColumns();
        this->rebuildAllPointers();
    }

    if (DEBUG_MODE) std::cout << "============================================================\n\n";
}

std::vector<BlockID> BlockSet::absorbMicroBlocksIndividual() {
    bool DEBUG_MODE = false; // 強制開啟 DEBUG 以印出你要求的資訊
    std::unordered_set<BlockID> modified_hosts;
    bool absorptionOccurred = true;

    int total_segments_absorbed = 0;
    int total_blocks_deleted = 0;
    int total_length_absorbed = 0;

    auto reverseComplement = [](const std::string& seq) {
        std::string rc = seq;
        std::reverse(rc.begin(), rc.end());
        for (char& c : rc) {
            if (c == 'A') c = 'T'; else if (c == 'T') c = 'A';
            else if (c == 'C') c = 'G'; else if (c == 'G') c = 'C';
        }
        return rc;
    };

    auto getGenomePrev = [](Segment& s) { return s.isReverse() ? s.getNextBlock() : s.getPrevBlock(); };
    auto getGenomeNext = [](Segment& s) { return s.isReverse() ? s.getPrevBlock() : s.getNextBlock(); };
    auto setGenomePrev = [](Segment& s, std::weak_ptr<Block> b) { if (s.isReverse()) s.setNextBlock(b); else s.setPrevBlock(b); };
    auto setGenomeNext = [](Segment& s, std::weak_ptr<Block> b) { if (s.isReverse()) s.setPrevBlock(b); else s.setNextBlock(b); };

    auto needFrontGap = [](bool isRightNeighbor, bool isReverse) {
        return (isRightNeighbor && !isReverse) || (!isRightNeighbor && isReverse);
    };

    while (absorptionOccurred) {
        absorptionOccurred = false;
        auto currentBlocks = this->getAllBlocks();

        for (auto mB_weak : currentBlocks) {
            auto mB = mB_weak.lock();
            if (!this->getBlock(mB->getId())) continue;
            
            int consLen = mB->getConsensus().length();
            if (consLen == 0 || consLen >= 30) continue; 

            std::vector<std::string> seqIDs;
            for (auto& kv : mB->getSequences()) seqIDs.push_back(kv.first);

            for (const auto& seqID : seqIDs) {
                if (mB->getSequences().count(seqID) == 0) continue;
                auto& mSegsMap = mB->getSequences().at(seqID).getSegments();
                
                std::vector<int> segKeys;
                for (auto& kv : mSegsMap) segKeys.push_back(kv.first);

                for (int mSegStartKey : segKeys) {
                    if (mB->getSequences().at(seqID).getSegments().count(mSegStartKey) == 0) continue;
                    Segment mSeg = mB->getSequences().at(seqID).getSegments().at(mSegStartKey);

                    std::string mSeq = "";
                    std::string mCons = mB->getConsensus();
                    for (int i = 0; i < mCons.length(); ++i) {
                        bool isGap = false;
                        char c = mCons[i];
                        for (auto& v : mSeg.getVariants()) {
                            if (v.getType() == VariantType::GAP && i >= v.getStart() && i < v.getEnd()) {
                                isGap = true; break;
                            } else if (v.getType() == VariantType::SNV && i == v.getStart()) {
                                c = v.getAlt();
                            }
                        }
                        if (!isGap) mSeq += c;
                    }
                    
                    int mLen = mSeq.length();
                    
                    if (mLen == 0) {
                        if (mSeg.getStart() == mSeg.getEnd()) {
                            auto bLeft = getGenomePrev(mSeg).lock();
                            auto bRight = getGenomeNext(mSeg).lock();
                            
                            if (bLeft && this->getBlock(bLeft->getId()) && bLeft->getSequences().count(seqID)) {
                                for (auto& pair : bLeft->getSequences().at(seqID).getSegments()) {
                                    if (getGenomeNext(pair.second).lock() == mB) setGenomeNext(pair.second, bRight);
                                }
                            }
                            if (bRight && this->getBlock(bRight->getId()) && bRight->getSequences().count(seqID)) {
                                for (auto& pair : bRight->getSequences().at(seqID).getSegments()) {
                                    if (getGenomePrev(pair.second).lock() == mB) setGenomePrev(pair.second, bLeft);
                                }
                            }
                            mB->getSequences().at(seqID).getSegments().erase(mSegStartKey);
                        }
                        continue; 
                    }

                    // 🌟 1. 探勘左邊與右邊的鄰居，並確認可用 Gap 大小
                    std::shared_ptr<Block> bLeft = nullptr; Segment sLeft; bool leftValid = false; int leftAvailGap = 0; bool leftIsFront = false;
                    std::shared_ptr<Block> bRight = nullptr; Segment sRight; bool rightValid = false; int rightAvailGap = 0; bool rightIsFront = false;

                    if (auto pB = getGenomePrev(mSeg).lock()) {
                        if (this->getBlock(pB->getId()) && pB != mB && pB->getSequences().count(seqID)) {
                            for (auto& hPair : pB->getSequences().at(seqID).getSegments()) {
                                if (getGenomeNext(hPair.second).lock() == mB && hPair.second.getEnd() == mSeg.getStart()) {
                                    bLeft = pB; sLeft = hPair.second; leftValid = true;
                                    leftIsFront = needFrontGap(false, sLeft.isReverse());
                                    int H = bLeft->getConsensus().length();
                                    if (leftIsFront && !sLeft.getVariants().empty() && sLeft.getVariants().front().getType() == VariantType::GAP && sLeft.getVariants().front().getStart() == 0)
                                        leftAvailGap = sLeft.getVariants().front().getEnd();
                                    else if (!leftIsFront && !sLeft.getVariants().empty() && sLeft.getVariants().back().getType() == VariantType::GAP && sLeft.getVariants().back().getEnd() == H)
                                        leftAvailGap = H - sLeft.getVariants().back().getStart();
                                    break;
                                }
                            }
                        }
                    }

                    if (auto nB = getGenomeNext(mSeg).lock()) {
                        if (this->getBlock(nB->getId()) && nB != mB && nB->getSequences().count(seqID)) {
                            for (auto& hPair : nB->getSequences().at(seqID).getSegments()) {
                                if (getGenomePrev(hPair.second).lock() == mB && hPair.second.getStart() == mSeg.getEnd()) {
                                    bRight = nB; sRight = hPair.second; rightValid = true;
                                    rightIsFront = needFrontGap(true, sRight.isReverse());
                                    int H = bRight->getConsensus().length();
                                    if (rightIsFront && !sRight.getVariants().empty() && sRight.getVariants().front().getType() == VariantType::GAP && sRight.getVariants().front().getStart() == 0)
                                        rightAvailGap = sRight.getVariants().front().getEnd();
                                    else if (!rightIsFront && !sRight.getVariants().empty() && sRight.getVariants().back().getType() == VariantType::GAP && sRight.getVariants().back().getEnd() == H)
                                        rightAvailGap = H - sRight.getVariants().back().getStart();
                                    break;
                                }
                            }
                        }
                    }

                    // 🌟 2. 智慧挑選最佳宿主 (Priority: 有足夠 Gap > 長度比較長)
                    std::shared_ptr<Block> hB = nullptr; Segment hSeg; bool isRightNeighbor = false; bool hostIsGenomeLeft = false; bool hostIsGenomeRight = false;
                    int available_gap = 0; bool isFrontGapNeeded = false;

                    bool leftHasGap = (leftValid && leftAvailGap >= mLen);
                    bool rightHasGap = (rightValid && rightAvailGap >= mLen);

                    if (leftHasGap && !rightHasGap) {
                        hB = bLeft; hSeg = sLeft; isRightNeighbor = false; hostIsGenomeLeft = true; available_gap = leftAvailGap; isFrontGapNeeded = leftIsFront;
                    } else if (!leftHasGap && rightHasGap) {
                        hB = bRight; hSeg = sRight; isRightNeighbor = true; hostIsGenomeRight = true; available_gap = rightAvailGap; isFrontGapNeeded = rightIsFront;
                    } else if (leftValid && rightValid) {
                        if (bLeft->getConsensus().length() >= bRight->getConsensus().length()) {
                            hB = bLeft; hSeg = sLeft; isRightNeighbor = false; hostIsGenomeLeft = true; available_gap = leftAvailGap; isFrontGapNeeded = leftIsFront;
                        } else {
                            hB = bRight; hSeg = sRight; isRightNeighbor = true; hostIsGenomeRight = true; available_gap = rightAvailGap; isFrontGapNeeded = rightIsFront;
                        }
                    } else if (leftValid) {
                        hB = bLeft; hSeg = sLeft; isRightNeighbor = false; hostIsGenomeLeft = true; available_gap = leftAvailGap; isFrontGapNeeded = leftIsFront;
                    } else if (rightValid) {
                        hB = bRight; hSeg = sRight; isRightNeighbor = true; hostIsGenomeRight = true; available_gap = rightAvailGap; isFrontGapNeeded = rightIsFront;
                    } else continue; // 沒有任何合法相鄰積木則跳過

                    if (hB) {
                        // 🌟 Debug 輸出：合併前的狀態
                        if (DEBUG_MODE) {
                            std::string pm = getGenomePrev(mSeg).lock() ? std::to_string(getGenomePrev(mSeg).lock()->getId()) : "NULL";
                            std::string nm = getGenomeNext(mSeg).lock() ? std::to_string(getGenomeNext(mSeg).lock()->getId()) : "NULL";
                            std::string ph = getGenomePrev(hSeg).lock() ? std::to_string(getGenomePrev(hSeg).lock()->getId()) : "NULL";
                            std::string nh = getGenomeNext(hSeg).lock() ? std::to_string(getGenomeNext(hSeg).lock()->getId()) : "NULL";
                            std::cout << "--------------------------------------------------\n";
                            std::cout << "[DEBUG] SHORT SEG: Prev=" << pm << "\tNext=" << nm << "\tLen=" << (mSeg.getEnd() - mSeg.getStart()) << "\t(Block " << mB->getId() << ")\n";
                            std::cout << "[DEBUG] MERGE TO : Prev=" << ph << "\tNext=" << nh << "\tLen=" << (hSeg.getEnd() - hSeg.getStart()) << "\t(Block " << hB->getId() << ")\n";
                        }

                        std::string mSeq_h = (mSeg.isReverse() == hSeg.isReverse()) ? mSeq : reverseComplement(mSeq);
                        int H = hB->getConsensus().length();
                        int hSegStartKey = hSeg.getStart();
                        auto& realHSeg = hB->getSequences().at(seqID).getSegments().at(hSegStartKey);
                        std::string hostCons = hB->getConsensus();

                        int overflow = std::max(0, mLen - available_gap);
                        int consume_gap = std::min(mLen, available_gap);

                        // 🌟 3. 變異與共識序列合併處理
                        if (overflow > 0) {
                            if (isFrontGapNeeded) {
                                hostCons = mSeq_h.substr(0, overflow) + hostCons;
                                hB->setConsensus(hostCons);
                                for (auto& seqPair : hB->getSequences()) {
                                    for (auto& segPair : seqPair.second.getSegments()) {
                                        Segment& s = segPair.second;
                                        auto& s_vars = s.getVariants();
                                        for (auto& v : s_vars) v.shift(overflow); // 全部平移
                                        
                                        if (&s != &realHSeg) {
                                            // 其他分支沒有這段擴張，補上 GAP
                                            if (!s_vars.empty() && s_vars.front().getType() == VariantType::GAP && s_vars.front().getStart() == overflow) s_vars.front().setStart(0);
                                            else s_vars.insert(s_vars.begin(), Variant::createGap(0, overflow));
                                        } else {
                                            if (available_gap > 0 && !s_vars.empty() && s_vars.front().getType() == VariantType::GAP && s_vars.front().getStart() == overflow) {
                                                s_vars.erase(s_vars.begin());
                                            }
                                            for (int i = 0; i < consume_gap; ++i) {
                                                int pos = overflow + i; char c = mSeq_h[pos];
                                                if (c != hostCons[pos]) s_vars.push_back(Variant(pos, c));
                                            }
                                        }
                                    }
                                }
                            } else {
                                hostCons = hostCons + mSeq_h.substr(consume_gap, overflow);
                                hB->setConsensus(hostCons);
                                for (auto& seqPair : hB->getSequences()) {
                                    for (auto& segPair : seqPair.second.getSegments()) {
                                        Segment& s = segPair.second;
                                        auto& s_vars = s.getVariants();
                                        if (&s != &realHSeg) {
                                            if (!s_vars.empty() && s_vars.back().getType() == VariantType::GAP && s_vars.back().getEnd() == H) s_vars.back().setEnd(H + overflow);
                                            else s_vars.push_back(Variant::createGap(H, H + overflow));
                                        } else {
                                            if (available_gap > 0 && !s_vars.empty() && s_vars.back().getType() == VariantType::GAP && s_vars.back().getEnd() == H) {
                                                s_vars.pop_back();
                                            }
                                            for (int i = 0; i < consume_gap; ++i) {
                                                int pos = H - available_gap + i; char c = mSeq_h[i];
                                                if (c != hostCons[pos]) s_vars.push_back(Variant(pos, c));
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            // Gap 夠大，純粹消耗 Gap，不擴張 Consensus
                            auto& s_vars = realHSeg.getVariants();
                            if (isFrontGapNeeded) {
                                if (mLen == available_gap) s_vars.erase(s_vars.begin());
                                else s_vars.front().setStart(mLen); // 消耗掉前面的 mLen
                                for (int i = 0; i < mLen; ++i) if (mSeq_h[i] != hostCons[i]) s_vars.push_back(Variant(i, mSeq_h[i]));
                            } else {
                                if (mLen == available_gap) s_vars.pop_back();
                                else s_vars.back().setEnd(H - mLen); // 消耗掉後面的 mLen
                                for (int i = 0; i < mLen; ++i) {
                                    int pos = H - mLen + i;
                                    if (mSeq_h[i] != hostCons[pos]) s_vars.push_back(Variant(pos, mSeq_h[i]));
                                }
                            }
                        }

                        // 變異重排序確保安全
                        for (auto& seqPair : hB->getSequences()) {
                            for (auto& segPair : seqPair.second.getSegments()) {
                                std::sort(segPair.second.getVariants().begin(), segPair.second.getVariants().end(), [](Variant& a, Variant& b){
                                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                                    return a.getType() > b.getType();
                                });
                            }
                        }

                        // 🌟 4. 指標重新縫合與座標更新
                        Segment finalSegCopy; // 複製一份拿來 Debug 印出
                        if (hostIsGenomeLeft) { 
                            realHSeg.setEnd(mSeg.getEnd());
                            setGenomeNext(realHSeg, getGenomeNext(mSeg)); 

                            auto bRight = getGenomeNext(mSeg).lock();
                            if (bRight && bRight->getSequences().count(seqID)) {
                                for (auto& pair : bRight->getSequences().at(seqID).getSegments()) {
                                    if (getGenomePrev(pair.second).lock() == mB) setGenomePrev(pair.second, hB);
                                }
                            }
                            finalSegCopy = realHSeg;
                        } else if (hostIsGenomeRight) { 
                            int oldStart = realHSeg.getStart();
                            realHSeg.setStart(mSeg.getStart());
                            setGenomePrev(realHSeg, getGenomePrev(mSeg)); 

                            auto bLeft = getGenomePrev(mSeg).lock();
                            if (bLeft && bLeft->getSequences().count(seqID)) {
                                for (auto& pair : bLeft->getSequences().at(seqID).getSegments()) {
                                    if (getGenomeNext(pair.second).lock() == mB) setGenomeNext(pair.second, hB);
                                }
                            }

                            if (realHSeg.getStart() != oldStart) {
                                Segment movedSeg = std::move(realHSeg);
                                hB->getSequences().at(seqID).getSegments().erase(oldStart);
                                hB->getSequences().at(seqID).getSegments()[movedSeg.getStart()] = std::move(movedSeg);
                                finalSegCopy = hB->getSequences().at(seqID).getSegments()[movedSeg.getStart()];
                            } else {
                                finalSegCopy = realHSeg;
                            }
                        }

                        // 🌟 Debug 輸出：合併後的結果 (驗證指標縫合是否成功)
                        if (DEBUG_MODE) {
                            std::string pf = getGenomePrev(finalSegCopy).lock() ? std::to_string(getGenomePrev(finalSegCopy).lock()->getId()) : "NULL";
                            std::string nf = getGenomeNext(finalSegCopy).lock() ? std::to_string(getGenomeNext(finalSegCopy).lock()->getId()) : "NULL";
                            std::cout << "[DEBUG] AFTER MRG: Prev=" << pf << "\tNext=" << nf << "\tLen=" << (finalSegCopy.getEnd() - finalSegCopy.getStart()) << "\n";
                        }

                        mB->getSequences().at(seqID).getSegments().erase(mSegStartKey);
                        if (mB->getSequences().at(seqID).getSegments().empty()) mB->getSequences().erase(seqID);

                        modified_hosts.insert(hB->getId());
                        absorptionOccurred = true;

                        total_segments_absorbed++;
                        total_length_absorbed += mLen;
                    }
                } 
            } 

            if (mB->getSequences().empty()) {
                total_blocks_deleted++;
                this->deleteBlock(mB->getId());
            }
        }
    }

    this->rebuildAllPointers();

    if (total_segments_absorbed > 0 || total_blocks_deleted > 0) {
        std::cout << "\n==================================================\n";
        std::cout << "          [ Micro-Block Absorption Summary ]        \n";
        std::cout << "==================================================\n";
        std::cout << " - Total Segments Absorbed : " << total_segments_absorbed << "\n";
        std::cout << " - Total Micro-Blocks Deleted : " << total_blocks_deleted << "\n";
        std::cout << " - Total Sequence Length Absorbed : " << total_length_absorbed << " bp\n";
        std::cout << " - Modified Host Blocks : " << modified_hosts.size() << "\n";
        std::cout << "==================================================\n\n";
    }

    return std::vector<BlockID>(modified_hosts.begin(), modified_hosts.end());
}


BlockIDs BlockSet::splitMultiBlocks(BlockID parentID, const std::vector<int>& cuts) {
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
        std::vector<Sequence> chunkSeqs;
    };
    std::vector<SplitResult> splitResults(seqIDs.size());

    tbb::parallel_for(tbb::blocked_range<size_t>(0, seqIDs.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                const std::string& seqID = seqIDs[i];
                auto& parentSeqInfo = parent->getSequences().at(seqID);

                // 預先分配記憶體 (避免 chunkSeqs 內部動態擴容)
                std::vector<Sequence> chunkSeqs(newBlocks.size(), Sequence(seqID));

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

    std::vector<BlockID> resultIDs;
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

void BlockSet::reconnectBlocks() {
    bool debug = true; // 如果你有全域的 DEBUG_MODE，可以直接替換掉這行
    
    if (debug) std::cout << "\n[BlockSet] 🔍 Starting Reconnect Blocks pass...\n";

    bool merged_any = false;
    int reconnectCount = 0; // 🌟 新增：追蹤總共合併了幾次
    auto linear_path = getLinearizeBlocks(); 
    
    if (linear_path.size() < 2) {
        if (debug) std::cout << "  -> Only " << linear_path.size() << " block(s) present. Skipping reconnect.\n";
        return;
    }

    // 輔助 Lambda：檢查兩個 Block 內的所有 Segments 是否完美頭尾相接
    auto isPerfectlyContiguous = [&](std::shared_ptr<Block> left, std::shared_ptr<Block> right) -> bool {
        if (!left || !right) return false;

        if (!left->isCoreBlock() || !right->isCoreBlock()) {
            return false;
        }

        auto& leftSeqs = left->getSequences();
        auto& rightSeqs = right->getSequences();
        
        if (leftSeqs.size() != rightSeqs.size()) return false;
        
        for (auto& [seqName, leftSeqInfo] : leftSeqs) {
            auto it = rightSeqs.find(seqName);
            if (it == rightSeqs.end()) return false;
            
            auto& lSegs = leftSeqInfo.getSegments();
            auto& rSegs = it->second.getSegments();
            
            if (lSegs.size() != rSegs.size()) return false;
            
            std::vector<bool> rUsed(rSegs.size(), false);
            
            for (auto& [lStart, lSeg] : lSegs) {
                bool matched = false;
                int rIdx = 0;
                
                for (auto& [rStart, rSeg] : rSegs) {
                    if (!rUsed[rIdx] && (lSeg.isReverse() == rSeg.isReverse())) {
                        bool isContiguousFwd = !lSeg.isReverse() && (lSeg.getEnd() == rSeg.getStart());
                        bool isContiguousRev = lSeg.isReverse()  && (lSeg.getStart() == rSeg.getEnd());
                        
                        if (isContiguousFwd || isContiguousRev) {
                            matched = true;
                            rUsed[rIdx] = true;
                            break;
                        }
                    }
                    rIdx++;
                }
                if (!matched) return false;
            }
        }
        return true;
    };

    auto connectBlocks = [&](std::shared_ptr<Block> left, std::shared_ptr<Block> right) {
        if (!left || !right) return std::shared_ptr<Block>(nullptr);

        int leftLen = left->getConsensus().length();
        int rightLen = right->getConsensus().length();
        
        // 1. 建立新的合併 Block
        auto newBlock = this->createBlock(left->getConsensus() + right->getConsensus());
        
        // 2. 收集左右兩邊涵蓋的所有 Sequence ID
        std::set<std::string> allSeqs;
        for (auto const& [seqName, info] : left->getSequences()) allSeqs.insert(seqName);
        for (auto const& [seqName, info] : right->getSequences()) allSeqs.insert(seqName);
        
        // 3. 遍歷每一條 Sequence，處理 Segment 的完美縫合或補 Gap
        for (auto& seqID : allSeqs) {
            Sequence newSeqInfo(seqID);
            std::vector<Segment> leftSegs, rightSegs;

            if (left->getSequences().count(seqID)) {
                for (auto const& [start, seg] : left->getSequences().at(seqID).getSegments()) {
                    leftSegs.push_back(seg);
                }
            }
            if (right->getSequences().count(seqID)) {
                for (auto const& [start, seg] : right->getSequences().at(seqID).getSegments()) {
                    rightSegs.push_back(seg);
                }
            }

            std::vector<bool> rUsed(rightSegs.size(), false);

            // 處理左側的 Segments，嘗試與右側對接
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
                            v.shift(leftLen); // 右側的 Variant 必須向右推移左側的長度

                            // 檢查 Gap 連續性，如果可以縫合就合併成一個大 Gap
                            if (!vars.empty() && vars.back().getType() == VariantType::GAP && 
                                v.getType() == VariantType::GAP && vars.back().getEnd() == v.getStart()) {
                                int oldStart = vars.back().getStart();
                                vars.pop_back();
                                vars.push_back(Variant::createGap(oldStart, v.getEnd()));
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

                // 如果左邊這個 Segment 找不到右邊的接續，代表右邊斷尾了，補上尾部 Gap
                if (!matched) {
                    Segment newSeg = lSeg;
                    auto& vars = newSeg.getVariants();
                    if (!vars.empty() && vars.back().getType() == VariantType::GAP && vars.back().getEnd() == leftLen) {
                        int oldStart = vars.back().getStart();
                        vars.pop_back();
                        vars.push_back(Variant::createGap(oldStart, leftLen + rightLen));
                    } else {
                        vars.push_back(Variant::createGap(leftLen, leftLen + rightLen));
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }

            // 處理右側剩餘的未配對 Segments，代表左邊斷頭了，補上頭部 Gap
            for (size_t i = 0; i < rightSegs.size(); ++i) {
                if (!rUsed[i]) {
                    Segment newSeg = rightSegs[i];
                    newSeg.getVariants().clear(); 
                    newSeg.getVariants().push_back(Variant::createGap(0, leftLen));

                    auto& vars = newSeg.getVariants();
                    for (auto v : rightSegs[i].getVariants()) {
                        v.shift(leftLen);
                        if (!vars.empty() && vars.back().getType() == VariantType::GAP && 
                            v.getType() == VariantType::GAP && vars.back().getEnd() == v.getStart()) {
                            int oldStart = vars.back().getStart();
                            vars.pop_back();
                            vars.push_back(Variant::createGap(oldStart, v.getEnd()));
                        } else {
                            vars.push_back(v);
                        }
                    }
                    newSeqInfo.getSegments()[newSeg.getStart()] = newSeg;
                }
            }

            if (!newSeqInfo.getSegments().empty()) {
                newBlock->addSequence(newSeqInfo);
            }
        }

        // 4. 尋找所有相鄰的 Neighbors (包含 Left 與 Right 的外圍鄰居)
        std::set<std::shared_ptr<Block>> neighbors;
        auto addNeighbors = [&](std::shared_ptr<Block> b) {
            for (auto& [seqName, seqInfo] : b->getSequences()) {
                for (auto& [segStart, seg] : seqInfo.getSegments()) {
                    if (auto sp = seg.getPrevBlock().lock()) neighbors.insert(sp);
                    if (auto sp = seg.getNextBlock().lock()) neighbors.insert(sp);
                }
            }
        };

        addNeighbors(left); 
        addNeighbors(right);
        neighbors.erase(left); 
        neighbors.erase(right);

        // 5. 🌟 關鍵修復：確保取得的是 Reference (auto&)，否則無法寫入記憶體
        for (auto& neighbor : neighbors) {
            // 注意：Block::getSequences() 必須回傳 std::unordered_map<...>& (Reference)
            for (auto& [seqName, seqInfo] : neighbor->getSequences()) {
                for (auto& [segStart, seg] : seqInfo.getSegments()) {

                    // 將所有原本指向 Left 或 Right 的連線，轉接到 newBlock
                    if (seg.getPrevBlock().lock() == left || seg.getPrevBlock().lock() == right) {
                        seg.setPrevBlock(newBlock);
                    }
                    if (seg.getNextBlock().lock() == left || seg.getNextBlock().lock() == right) {
                        seg.setNextBlock(newBlock);
                    }
                }
            }
        }

        // 6. 刪除舊的碎塊並回傳新合併的 Block
        this->deleteBlock(left->getId());
        this->deleteBlock(right->getId());

        return newBlock;
    };

    // 🌟 修改：使用 auto 兼容 VBlockID (即 std::pair<BlockID, int>)
    auto it = linear_path.begin();
    while (it != linear_path.end() && std::next(it) != linear_path.end()) {
        
        // 解析 VBlockID
        auto leftVId = *it;
        auto rightVId = *std::next(it);
        
        BlockID leftId = leftVId.first;
        BlockID rightId = rightVId.first;
        
        auto leftBlock = getBlock(leftId);
        auto rightBlock = getBlock(rightId);
        
        if (leftBlock && rightBlock && isPerfectlyContiguous(leftBlock, rightBlock)) {
            
            auto newBlock = connectBlocks(leftBlock, rightBlock); 
            
            if (debug) {
                std::cout << "  -> 🔗 [RECONNECT] Perfectly contiguous! Merged Core Block " 
                          << leftId << " + Core Block " << rightId 
                          << " => New Core Block " << newBlock->getId() << "\n";
            }
            
            // 🌟 將新的實體 BlockID 包裝回 VBlockID (因為是 Core Block，Copy 固定延續)
            *it = {newBlock->getId(), leftVId.second}; 
            linear_path.erase(std::next(it));
            
            merged_any = true;
            reconnectCount++; 
            
        } else {
            ++it; 
        }
    }

    if (merged_any) {
        invalidateRepCache();
    }

    this->rebuildAllPointers();

    if (debug) {
        if (reconnectCount > 0) {
            std::cout << "  ✅ [RECONNECT SUMMARY] Successfully reconnected " << reconnectCount << " core block pairs.\n";
        } else {
            std::cout << "  ✅ [RECONNECT SUMMARY] No contiguous core blocks needed reconnection.\n";
        }
    }
}


/*
void BlockSet::popFuzzyBubbles() {
    bool DEBUG_MODE = true;
    if (DEBUG_MODE) std::cout << "\n============================================================\n"
                              << ">>> Executing Fuzzy Bubble Popping (Alignment-based)...\n"
                              << "============================================================\n";

    const int MAX_BUBBLE_LEN = 500;
    const int MAX_LEN_DIFF = 20;
    const double MIN_IDENTITY = 0.85;

    int total_merged = 0;
    bool mergedOccurred = true;

    // 太極迴圈：只要有發生合併，圖譜拓撲就改變了，必須重新找一次 Bubble
    while (mergedOccurred) {
        mergedOccurred = false;
        
        // 1. 尋找氣泡：使用 (Genomic Left ID, Genomic Right ID) 作為氣泡空間特徵
        std::map<std::pair<BlockID, BlockID>, std::vector<BlockID>> bubble_groups;

        for (auto blk : this->getAllBlocks()) {
            if (!blk) continue;

            BlockID left_id = 0, right_id = 0;
            bool is_consistent = true;
            bool first = true;

            for (auto& seqPair : blk->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    Segment& seg = segPairInner.second;
                    // 根據方向性取得 Genomic 絕對左邊和右邊的 Hub
                    auto l_ptr = !seg.isReverse() ? seg.getPrevBlock().lock() : seg.getNextBlock().lock();
                    auto r_ptr = !seg.isReverse() ? seg.getNextBlock().lock() : seg.getPrevBlock().lock();

                    BlockID cur_l = l_ptr ? l_ptr->getId() : 0;
                    BlockID cur_r = r_ptr ? r_ptr->getId() : 0;

                    if (first) {
                        left_id = cur_l; right_id = cur_r; first = false;
                    } else {
                        if (cur_l != left_id || cur_r != right_id) {
                            is_consistent = false; break; // 分支跑去別的地方了，不是單純的氣泡
                        }
                    }
                }
                if (!is_consistent) break;
            }

            if (is_consistent && left_id != 0 && right_id != 0 && left_id != right_id) {
                // 為了無向圖的穩定性，永遠把較小的 ID 放前面作為 Key
                BlockID min_id = std::min(left_id, right_id);
                BlockID max_id = std::max(left_id, right_id);
                bubble_groups[{min_id, max_id}].push_back(blk->getId());
            }
        }

        // 2. 針對同一個 Bubble 空間內的分支，嘗試進行 Global Alignment
        for (auto& kv : bubble_groups) {
            auto& branches = kv.second;
            if (branches.size() < 2) continue; // 沒有平行分支

            for (size_t i = 0; i < branches.size(); ++i) {
                for (size_t j = i + 1; j < branches.size(); ++j) {
                    auto b1 = this->getBlock(branches[i]);
                    auto b2 = this->getBlock(branches[j]);
                    if (!b1 || !b2) continue;

                    std::string seq1 = b1->getConsensus();
                    std::string seq2 = b2->getConsensus();
                    int len1 = seq1.length();
                    int len2 = seq2.length();

                    // 【條件 1】：長度都小於 500，且長度差異小於 20
                    if (len1 < MAX_BUBBLE_LEN && len2 < MAX_BUBBLE_LEN && std::abs(len1 - len2) <= MAX_LEN_DIFF) {
                        
                        // 【條件 2】：執行 Global Alignment
                        auto res = runGlobalAlignment(seq1, seq2); 

                        // 【條件 3】：Identity 達標
                        if (res.success && res.identity >= MIN_IDENTITY) {
                            if (DEBUG_MODE) {
                                std::cout << "  [Fuzzy Merge] Anchor (" << kv.first.first << ", " << kv.first.second << ")\n"
                                          << "    -> Merging Block " << branches[i] << " (Len: " << len1 << ") and Block " 
                                          << branches[j] << " (Len: " << len2 << ")\n"
                                          << "    -> Identity: " << res.identity << " | Score: " << res.score << "\n";
                            }

                            // 呼叫你已經寫好的 Merge 函數，把 b2 揉進 b1 (注意：你需要確定你的 mergeTwoBlocks 怎麼傳參數)
                            // 假設你的 mergeTwoBlocks 是 this->mergeTwoBlocks(hostBlock, donorBlock, alignmentResult);
                            this->mergeTwoBlocks(b1, b2, res.cigar, false); 
                            vim
                            total_merged++;
                            mergedOccurred = true;
                            break; // 只要發生合併，指標就亂了，立刻 break 讓迴圈重來
                        }
                    }
                }
                if (mergedOccurred) break; 
            }
            if (mergedOccurred) break;
        }
    }

    if (total_merged > 0) {
        if (DEBUG_MODE) std::cout << "  -> Success! Popped " << total_merged << " fuzzy bubbles.\n";
        // 如果你的 mergeTwoBlocks 裡面沒有呼叫重建，記得在這裡呼叫
        this->rebuildAllPointers(); 
    } else {
        if (DEBUG_MODE) std::cout << "  -> No eligible fuzzy bubbles found.\n";
    }
    if (DEBUG_MODE) std::cout << "============================================================\n\n";
}
*/