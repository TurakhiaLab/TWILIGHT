
#include "block.hpp"

#include <cmath>
#include <iomanip>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>


void BlockSet::debugValidateSegments(bool verbose) {
    
    if (verbose) {
        std::cout << "\n============================================================\n"
                  << "=== BlockSet Debug Validation: " << ID << " ===\n"
                  << "============================================================\n";
    }

    int totalBlocks = 0;
    int totalSegments = 0;
    int errorCount = 0;

    for (const auto& blockPair : blocks) {
        std::shared_ptr<Block> blk = blockPair.second;
        int consLen = blk->getConsensus().length();
        totalBlocks++;

        // How many segments are in the block
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
                for (auto& var : seg.getVariants()) {
                    if (var.getType() == VariantType::GAP) {
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
    
    if (verbose || errorCount > 0) {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "Validation Complete. Checked " << totalBlocks << " blocks and " << totalSegments << " segments.\n";
    }
    if (errorCount == 0) {
        if (verbose) std::cout << "🎉 PERFECT! All segment lengths match their block consensus length perfectly.\n";
    } else {
        std::cout << "🚨 CRITICAL: FOUND " << errorCount << " LENGTH MISMATCH ERROR(S)!\n";
        exit(1);
    }
    if (verbose || errorCount > 0) {
        std::cout << "============================================================\n\n";
    }
}

void BlockSet::debugValidateLinkages(bool verbose) {
    if (verbose) {
        std::cout << "\n============================================================\n"
                  << "=== Topology & Linkage Debug Validation: " << ID << " ===\n"
                  << "============================================================\n";
    }

    struct SegRef { Segment* seg; std::shared_ptr<Block> blk; };
    std::map<std::string, std::vector<SegRef>> seqTracks;
    
    // 1. Collect all segment
    for (const auto& blockPair : blocks) {
        for (auto& seqPair : blockPair.second->getSequences()) {
            for (auto& segPairInner : seqPair.second.getSegments()) {
                seqTracks[seqPair.first].push_back({ &segPairInner.second, blockPair.second });
            }
        }
    }

    int pointerErrorCount = 0;
    int coordinateGapCount = 0;

    // 2. Check all sequences
    for (auto& trackPair : seqTracks) {
        auto& track = trackPair.second;
        
        // Sort by genome coordinate
        std::sort(track.begin(), track.end(), [](const SegRef& a, const SegRef& b) {
            return std::min(a.seg->getStart(), a.seg->getEnd()) < std::min(b.seg->getStart(), b.seg->getEnd());
        });

        if (verbose) std::cout << "Checking Sequence Track: " << trackPair.first << " (" << track.size() << " segments)\n";

        for (size_t i = 0; i < track.size(); ++i) {
            auto& currRef = track[i];

            if (i > 0) {
                auto& prevRef = track[i-1];

                int prevEnd = std::max(prevRef.seg->getStart(), prevRef.seg->getEnd());
                int currStart = std::min(currRef.seg->getStart(), currRef.seg->getEnd());
                
                if (prevEnd != currStart) {
                    std::cerr << "  ⚠️ [WARNING] Coordinate Gap: Block " << prevRef.blk->getId() 
                              << " end(" << prevEnd << ") != Block " << currRef.blk->getId() 
                              << " start(" << currStart << ")\n";
                    coordinateGapCount++;
                }

                std::shared_ptr<Block> expectedNextForPrev = currRef.blk;
                std::shared_ptr<Block> expectedPrevForCurr = prevRef.blk;
                std::shared_ptr<Block> actualNextForPrev;
                std::shared_ptr<Block> actualPrevForCurr;

                if (this->sequence_names.size() > 1) { // Pangenome graph
                    actualNextForPrev = (!prevRef.seg->isReverse()) ? prevRef.seg->getNextBlock().lock() : prevRef.seg->getPrevBlock().lock();
                    actualPrevForCurr = (!currRef.seg->isReverse()) ? currRef.seg->getPrevBlock().lock() : currRef.seg->getNextBlock().lock();
                } 
                else { // Self-mapping
                    actualNextForPrev = prevRef.seg->getNextBlock().lock();
                    actualPrevForCurr = currRef.seg->getPrevBlock().lock();
                }
                
                if (actualNextForPrev != expectedNextForPrev) {
                    std::cerr << "  ❌ [ERROR] Broken Pointer (Forward): Block " << prevRef.blk->getId() 
                              << " does NOT point to Block " << currRef.blk->getId() << " (instead: Block " 
                              << actualNextForPrev->getId() << ")\n";
                    pointerErrorCount++;
                }

                if (actualPrevForCurr != expectedPrevForCurr) {
                    std::cerr << "  ❌ [ERROR] Broken Pointer (Backward): Block " << currRef.blk->getId() 
                              << " does NOT point back to Block " << prevRef.blk->getId() << " (instead: Block " 
                              << actualPrevForCurr->getId() << ")\n";
                    pointerErrorCount++;
                }
            }
        }
    }

    
    if (pointerErrorCount == 0 && coordinateGapCount == 0) {
        if (verbose) std::cout << "🎉 PERFECT! All Graph linkages and coordinates are contiguous and sound.\n";
    } else {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "🚨 SUMMARY: Found " << pointerErrorCount << " Pointer Error(s) and " 
                  << coordinateGapCount << " Coordinate Gap(s)!\n";
        std::cout << "============================================================\n\n";
    }
    
}

void BlockSet::debugValidateQuality(bool verbose) {
    std::cout << "\n============================================================\n"
              << "=== Pangenome Graph Quality Report (NEW): " << ID << " ===\n"
              << "============================================================\n";

    if (blocks.empty()) {
        std::cout << "  [Warning] Graph is empty. No metrics to calculate.\n";
        return;
    }

    std::map<std::string, int> seqRealLengths;
    int totalSequenceCount = 0;
    
    for (auto& blkPair : blocks) {
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

    int threshold95 = std::ceil(totalSequenceCount * 0.95);
    int threshold90 = std::ceil(totalSequenceCount * 0.90);

    uint64_t totalConsensusLen = 0;
    uint64_t singletonLenSum = 0;
    int singletonCount = 0;
    
    uint64_t globalVarLen = 0;
    uint64_t globalDenominator = 0;

    uint64_t softCore90VarLen = 0;
    uint64_t softCore90Denominator = 0;
    uint64_t narrowCoreVarLen = 0;
    uint64_t narrowCoreDenominator = 0;

    std::vector<int> allBlockLengths;
    
    int coreBlocksCount = 0;
    int softCore95BlocksCount = 0;
    int softCore90BlocksCount = 0;
    int accessoryBlocksCount = 0;
    int narrowCoreBlocksCount = 0; 
    
    uint64_t coreLenSum = 0;
    uint64_t softCore95LenSum = 0;
    uint64_t softCore90LenSum = 0;
    uint64_t accessoryLenSum = 0;
    uint64_t narrowCoreLenSum = 0;

    if (verbose) std::cout << "[Block-level Identity Info]\n";

    for (const auto& blkPair : blocks) {
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
                for (auto& var : segPair.second.getVariants()) {
                    if (var.getType() == VariantType::SNV) {
                        blockVarLen += 1;
                    } else if (var.getType() == VariantType::GAP) {
                        blockVarLen += (var.getEnd() - var.getStart());
                    }
                }
            }
        }

        // 統計 Core, Soft-core (95%, 90%), Accessory
        int uniqueSeqCount = uniqueSeqsInBlock.size();
        
        if (uniqueSeqCount == totalSequenceCount && blockSegCount == totalSequenceCount && totalSequenceCount > 1) {
            narrowCoreBlocksCount++;
            narrowCoreLenSum += consLen;
        }
        
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

            if (uniqueSeqCount >= threshold90 && totalSequenceCount > 1) {
                softCore90VarLen += blockVarLen;
                softCore90Denominator += blockDenominator;
            }
            if (uniqueSeqCount == totalSequenceCount && blockSegCount == totalSequenceCount && totalSequenceCount > 1) {
                narrowCoreVarLen += blockVarLen;
                narrowCoreDenominator += blockDenominator;
            }
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
    double softCore90Identity = (softCore90Denominator > 0) ? (1.0 - ((double)softCore90VarLen / softCore90Denominator)) : 0.0;
    double narrowCoreIdentity = (narrowCoreDenominator > 0) ? (1.0 - ((double)narrowCoreVarLen / narrowCoreDenominator)) : 0.0;

    std::cout << "\n------------------------------------------------------------\n";
    std::cout << ">>> GRAPH METRICS SUMMARY <<<\n\n";

    std::cout << "[1. Sequence & Graph Size]\n";
    std::cout << "  - Total Sequences        : " << totalSequenceCount << "\n";
    std::cout << "  - Total Blocks           : " << blocks.size() << "\n";
    std::cout << "  - Longest Input Sequence : " << maxSeqLen << " bp (" << longestSeqName << ")\n";
    std::cout << "  - Total Graph Length     : " << totalConsensusLen << " bp\n";
    std::cout << "  - Graph Size Inflation   : +" << std::fixed << std::setprecision(2) << lenIncreaseRatio << " %\n";
    
    std::cout << "\n[2. Fragmentation & Contiguity]\n";
    std::cout << "  - Block N50              : " << n50 << " bp\n";
    std::cout << "  - Singleton Blocks       : " << singletonCount << " blocks\n";
    std::cout << "  - Singleton Length       : " << singletonLenSum << " bp\n";
    std::cout << "  - Singleton Length Ratio : " << std::fixed << std::setprecision(2) << singletonRatio << " %\n";

    std::cout << "\n[3. Evolution & Conservation]\n";
    std::cout << "  - Narrow Core (1-to-1)   : " << narrowCoreBlocksCount << " blocks (" << narrowCoreLenSum << " bp)\n";
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
    if (softCore90Denominator > 0) {
        std::cout << "  - Broad Core (>=90%) Identity: " << std::fixed << std::setprecision(4) << softCore90Identity << "\n";
    } else {
        std::cout << "  - Broad Core (>=90%) Identity: N/A\n";
    }
    if (narrowCoreDenominator > 0) {
        std::cout << "  - Narrow Core (1-to-1) Identity: " << std::fixed << std::setprecision(4) << narrowCoreIdentity << "\n";
    } else {
        std::cout << "  - Narrow Core (1-to-1) Identity: N/A\n";
    }

    std::cout << "============================================================\n\n";
}

void BlockSet::debugValidateSequences(BlockManager* manager, bool verbose) {
    if (verbose) {
        std::cout << "\n============================================================\n"
                  << "=== Sequence Debug Validation: " << ID << " ===\n"
                  << "============================================================\n";
    }
    for (auto& seqName: this->getSequences()) {
        auto seq_after = this->reconstructSequence(seqName);
        auto seq_before = manager->getSequence(seqName);
        bool pass = (seq_after == seq_before);
        if (pass) {
            if (verbose) {
                std::cout << "Validate Sequence " << seqName << '\n';
                std::cout << "  ✅ [PERFECT] Validation Passed! Reconstructed sequence perfectly matches the raw sequence. (Len: " << seq_after.length() << " bp)\n";    
            }
        }
        else {
            std::cout << "Validate Sequence " << seqName << '\n';
            std::cerr << "  ❌ [CRITICAL ERROR] Validation Failed! Sequence mismatch. (Len: " << seq_before.length() << " bp -> " <<  seq_after.length() << " bp)\n"; 
        }
    }
}

void BlockSet::debugValidateLinearizedBlocks(bool verbose) {
    auto linearBlocks = this->getLinearizeBlocks();
    if (linearBlocks.empty()) return;

    std::unordered_map<BlockID, int> id_to_index;
    for (int i = 0; i < (int)linearBlocks.size(); ++i) {
        id_to_index[linearBlocks[i]] = i;
    }

    int N = (int)linearBlocks.size();
    bool all_valid = true;

    // 2. 收集所有序列在圖中的完整路徑 (Path)
    std::unordered_map<std::string, std::vector<int>> seq_path_indices;
    for (int i = 0; i < N; ++i) {
        auto blk = this->getBlock(linearBlocks[i]);
        for (auto& [seqName, seqData] : blk->getSequences()) {
            
            // 🌟 核心修正：必須深入遍歷每一個 Segment！
            // 這樣如果同一個 Block 有兩個片段，`i` 就會被 push_back 兩次。
            for (auto& segPair : seqData.getSegments()) {
                seq_path_indices[seqName].push_back(i);
            }
            
        }
    }

    // 3. 檢查每一條 Sequence 的順序是否連貫
    for (auto& [seqName, path] : seq_path_indices) {
        if (path.size() <= 1) continue;

        for (size_t k = 0; k < path.size() - 1; ++k) {
            int curr_idx = path[k];
            int next_idx = path[k+1];

            // 正常情況：嚴格遞增
            if (next_idx > curr_idx) {
                continue; 
            } 
            // 環狀情況：從最後面跳回最前面
            else if (curr_idx == N - 1 && next_idx == 0) {
                if (verbose) std::cout << "  [Info] Seq " << seqName << " circular wrap-around detected (" << curr_idx << " -> " << next_idx << ").\n";
                continue;
            }
            // 🌟 異常情況：順序倒退，或者是「自我迴圈 (Self-loop)」！
            else {
                std::cerr << "  🚨 [ERROR] Order violation in Seq " << seqName 
                          << ": Block " << linearBlocks[curr_idx] << " (index " << curr_idx 
                          << ") -> Block " << linearBlocks[next_idx] << " (index " << next_idx << ")\n";
                          
                // 順便幫你加上專屬的 Self-loop 提示，方便你 Debug
                if (curr_idx == next_idx) {
                    std::cerr << "      -> 🔥 Reason: Self-loop detected! Multiple segments exist in the same block.\n";
                }
                
                all_valid = false;
            }
        }
    }

    if (all_valid) {
        if (verbose) std::cout << "  ✅ [Validation] All sequences follow a valid linearized order.\n";
    } else {
        std::cerr << "  ❌ [Validation] Found order violations in BlockSet " << ID << "!\n";
    }
}

void BlockSet::debugValidateBlocks(bool verbose) {
    auto linearBlocks = this->getLinearizeBlocks();
    if (linearBlocks.empty()) {
        std::cout << "\n  [Warning] Graph is empty. No blocks to validate.\n";
        return;
    }

    std::cout << "\n=========================================================================================================\n"
              << "=== Linearized Block Distribution & Quality Report ===\n"
              << "=========================================================================================================\n";

    // 定義統計級距 (Bins)
    struct SizeBin {
        std::string label;
        int min_val;
        int max_val;
        int block_count = 0;
        double sum_identity = 0.0;
        int sum_sequences = 0;
        uint64_t sum_snvs = 0;       // 🌟 新增：紀錄該級距的總 SNV
        uint64_t sum_gaps = 0;       // 🌟 新增：紀錄該級距的總 Gap 長度
    };

    std::vector<SizeBin> bins = {
        {"1-10", 1, 10},
        {"11-50", 11, 50},
        {"51-100", 51, 100},
        {"101-200", 101, 200},
        {"201-500", 201, 500},
        {"501-1000", 501, 1000},
        {"1001-10000", 1001, 10000},
        {"10000+", 10001, INT32_MAX}
    };

    if (verbose) {
        std::cout << "[Block-level Details (Linear Order)]\n";
    }

    for (BlockID id : linearBlocks) {
        auto blk = this->getBlock(id);
        if (!blk) continue;

        int consLen = blk->getConsensus().length();
        int famId = blk->getFamilyId();
        int seqCount = blk->getSequences().size();

        // 🌟 分開統計 SNV 數量與 Total Gap 長度
        uint64_t snv_count = 0;
        uint64_t total_gap_len = 0;
        int segCount = 0;

        for (auto& [seqName, seqData] : blk->getSequences()) {
            for (auto& [segId, segNode] : seqData.getSegments()) {
                segCount++;
                for (auto& var : segNode.getVariants()) {
                    if (var.getType() == VariantType::SNV) {
                        snv_count += 1;
                    } else if (var.getType() == VariantType::GAP) {
                        total_gap_len += (var.getEnd() - var.getStart());
                    }
                }
            }
        }

        // 計算總變異長度供 Identity 使用
        uint64_t varLen = snv_count + total_gap_len;

        double identity = 1.0; // 預設 1.0 (包含 Singleton)
        if (segCount > 1 && consLen > 0) {
            uint64_t denom = (uint64_t)segCount * consLen;
            identity = 1.0 - ((double)varLen / denom);
        }

        // 🌟 單行詳細輸出
        bool low_identity = (identity < 0.995);
        if (verbose && low_identity) {
            std::cout << "  ├─ Block ID: " << std::setw(6) << std::left << id 
                      << " | Fam: " << std::setw(6) << famId
                      << " | Len: " << std::setw(7) << consLen
                      << " | Seqs: " << std::setw(3) << seqCount
                      << " | Identity: " << std::fixed << std::setprecision(4) << identity
                      << " | SNVs: " << std::setw(5) << snv_count
                      << " | TotalGapLen: " << std::setw(6) << total_gap_len;
                      
            if (segCount == 1) std::cout << " (Singleton)";
            std::cout << "\n";
        }

        // 分類到對應的級距中
        for (auto& bin : bins) {
            if (consLen >= bin.min_val && consLen <= bin.max_val) {
                bin.block_count++;
                bin.sum_sequences += seqCount;
                bin.sum_identity += identity;
                bin.sum_snvs += snv_count;           // 🌟 累積 SNV 到對應 Bin
                bin.sum_gaps += total_gap_len;       // 🌟 累積 Gap 到對應 Bin
                break;
            }
        }
    }

    // ==========================================
    // 印出統計結果表格
    // ==========================================
    std::cout << "\n[Block Size Distribution & Statistics]\n";
    std::cout << std::string(105, '-') << "\n"; // 🌟 拉寬表格分隔線
    std::cout << std::setw(15) << std::left << "Size Range"
              << "| " << std::setw(10) << "Blocks"
              << "| " << std::setw(12) << "Avg Seqs"
              << "| " << std::setw(15) << "Avg Identity" 
              << "| " << std::setw(15) << "Total SNVs"       // 🌟 表格新增欄位
              << "| " << std::setw(15) << "Total GapLen"     // 🌟 表格新增欄位
              << "\n";
    std::cout << std::string(105, '-') << "\n";

    int total_blocks = 0;
    uint64_t grand_total_snvs = 0; // 🌟 總和變數
    uint64_t grand_total_gaps = 0; // 🌟 總和變數

    for (const auto& bin : bins) {
        if (bin.block_count == 0) continue; // 跳過沒有積木的級距
        
        double avg_seqs = (double)bin.sum_sequences / bin.block_count;
        double avg_identity = bin.sum_identity / bin.block_count; // 直接無權重平均

        std::cout << std::setw(15) << std::left << bin.label
                  << "| " << std::setw(10) << bin.block_count
                  << "| " << std::setw(12) << std::fixed << std::setprecision(2) << avg_seqs
                  << "| " << std::setw(15) << std::fixed << std::setprecision(4) << avg_identity
                  << "| " << std::setw(15) << bin.sum_snvs       // 🌟 印出該級距的 SNVs
                  << "| " << std::setw(15) << bin.sum_gaps       // 🌟 印出該級距的 Gaps
                  << "\n";
                  
        total_blocks += bin.block_count;
        grand_total_snvs += bin.sum_snvs; // 累加至全域總和
        grand_total_gaps += bin.sum_gaps; // 累加至全域總和
    }
    
    // 🌟 全域最終數據 Summary
    std::cout << std::string(105, '-') << "\n";
    std::cout << "Total Linear Blocks Evaluated : " << total_blocks << "\n";
    std::cout << "Total SNVs in Graph           : " << grand_total_snvs << "\n";
    std::cout << "Total Gap Length in Graph     : " << grand_total_gaps << "\n";
    std::cout << "=========================================================================================================\n\n";
}

/*
void BlockSet::debugValidateBubble(bool verbose) {
    std::cout << "\n============================================================\n";
    std::cout << "=== BlockSet: Bubble Topology Validator (Raw Consensus) ===\n";
    std::cout << "============================================================\n";

    // 尋找氣泡：使用 Genomic 左右兩端的 Hub Block ID 作為 Key
    std::map<std::pair<Block::ID, Block::ID>, std::vector<Block::ID>> bubble_groups;

    for (auto blk : this->getAllBlocks()) {
        if (!blk) continue;

        Block::ID left_id = 0, right_id = 0;
        bool is_consistent = true;
        bool first = true;

        for (auto& seqPair : blk->getSequences()) {
            for (auto& segPairInner : seqPair.second.getSegments()) {
                Segment& seg = segPairInner.second;
                
                // 根據方向性，精準取得 Genomic 上的左邊與右邊積木
                auto l_ptr = !seg.isReverse() ? seg.getPrevBlock().lock() : seg.getNextBlock().lock();
                auto r_ptr = !seg.isReverse() ? seg.getNextBlock().lock() : seg.getPrevBlock().lock();

                Block::ID cur_l = l_ptr ? l_ptr->getId() : 0;
                Block::ID cur_r = r_ptr ? r_ptr->getId() : 0;

                if (first) {
                    left_id = cur_l; right_id = cur_r; first = false;
                } else {
                    if (cur_l != left_id || cur_r != right_id) {
                        is_consistent = false; break;
                    }
                }
            }
            if (!is_consistent) break;
        }

        // 如果這個積木完美夾在兩個 Hub 之間，就把它加入該氣泡群組
        if (is_consistent && left_id != 0 && right_id != 0 && left_id != right_id) {
            Block::ID min_id = std::min(left_id, right_id);
            Block::ID max_id = std::max(left_id, right_id);
            bubble_groups[{min_id, max_id}].push_back(blk->getId());
        }
    }

    int singletonBubbleLen = 0;

    int bubble_count = 0;
    for (auto& kv : bubble_groups) {
        if (kv.second.size() < 2) continue; // 只有一條路徑，不是氣泡
        
        bubble_count++;
        Block::ID left_id = kv.first.first;
        Block::ID right_id = kv.first.second;
        auto leftBlk = this->getBlock(left_id);
        auto rightBlk = this->getBlock(right_id);

        int left_seg = 0, right_seg = 0;
        for (auto& seq: leftBlk->getSequences()) 
            for (auto& seg: seq.second.getSegments()) left_seg++;

        for (auto& seq: rightBlk->getSequences()) 
            for (auto& seg: seq.second.getSegments()) right_seg++;

        std::cout << "\n[Bubble #" << bubble_count << "]\n";
        std::cout << "  ├─ Flanking Hub A (ID: " << left_id << ") | Len: " 
                  << (leftBlk ? leftBlk->getConsensus().length() : 0) << " bp | Seqs: " 
                  << (leftBlk ? left_seg : 0) << "\n";
        std::cout << "  ├─ Flanking Hub B (ID: " << right_id << ") | Len: " 
                  << (rightBlk ? rightBlk->getConsensus().length() : 0) << " bp | Seqs: " 
                  << (rightBlk ? right_seg : 0) << "\n";
        std::cout << "  └─ Branches (" << kv.second.size() << " parallel paths):\n\n";

        for (Block::ID bId : kv.second) {
            auto bBlk = this->getBlock(bId);
            if (!bBlk) continue;

            int bSeg = 0;
            for (auto& seq: bBlk->getSequences()) 
                for (auto& seg: seq.second.getSegments()) bSeg++;

            // 判斷方向：如果這條分支的序列在 Genomic 上是反向，我們印出來前先把它轉正
            bool is_rev = false;
            if (!bBlk->getSequences().empty() && !bBlk->getSequences().begin()->second.getSegments().empty()) {
                is_rev = bBlk->getSequences().begin()->second.getSegments().begin()->second.isReverse();
            }

            std::string seqStr = bBlk->getConsensus();
            if (is_rev) {
                std::reverse(seqStr.begin(), seqStr.end());
                for (char& c : seqStr) {
                    switch (c) {
                        case 'A': c='T'; break; case 'T': c='A'; break;
                        case 'C': c='G'; break; case 'G': c='C'; break;
                        case 'a': c='t'; break; case 't': c='a'; break;
                        case 'c': c='g'; break; case 'g': c='c'; break;
                    }
                }
            }

            // 直接印出這塊分支的 ID、長度、包含序列數，以及最純粹的 Consensus
            std::cout << "       [Block " << std::left << std::setw(6) << bId << "] "
                      << "(Len: " << std::setw(4) << seqStr.length() << " bp, Seqs: " << std::setw(2) << bSeg << ")\n"
                      << "       Seq: " << seqStr << "\n\n";

            if (bSeg == 1) singletonBubbleLen+= seqStr.length();
        }
        std::cout << "------------------------------------------------------------\n";
    }
    
    if (bubble_count == 0) {
        std::cout << "  -> No simple bubbles found!\n";
    } else {
        std::cout << "  -> Total " << bubble_count << " bubbles identified and validated.\n";
        std::cout << "  -> Total singleton bubble length: " << singletonBubbleLen << " bp\n";
    }
    std::cout << "============================================================\n\n";
}
*/