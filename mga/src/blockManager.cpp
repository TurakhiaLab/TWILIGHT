#include "block_manager.hpp"
#include "mga.hpp"
#include "type.hpp"
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <chrono>

// Helper to truncate strings for printing
std::string truncate(const std::string& str, size_t width) {
    if (str.length() > width) {
        return str.substr(0, width - 3) + "...";
    }
    return str;
}





void BlockManager::orientCircularGenomes(Option& option, std::string refSequenceName) {
    bool debug = true;
    bool filter = false;
    int target_num [2] = {20, 30};
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

    // 2. 擷取 Reference 和 Query 的序列 (使用 SequenceRefs 避免複製序列)
    BlockSet* refBlockSet = blocksets[bestRefName].get();
    SequenceRefs reference_genome;
    reference_genome.push_back({refBlockSet->getId(), refBlockSet->getAncestralSequence()});
    
    SequenceRefs query_genomes;
    query_genomes.reserve(this->blocksets.size() - 1);
    for (const auto& pair : this->blocksets) {
        if (pair.first != bestRefName) {
            BlockSet* qryBlockSet = pair.second.get();
            query_genomes.push_back({qryBlockSet->getId(), qryBlockSet->getAncestralSequence()});
        }
    }
    log_step_time("Prepared Queries for minimap2");

    // 3. Run Minimap2 (orientation only requires chaining length, CIGAR is disabled and max_gap is reduced for speed)
    if (debug) std::cout << "[DEBUG] Running chaining-only minimap2 on " << query_genomes.size() << " queries...\n";
    
    mm_config_orientation fastConfig(option);
    
    auto allAlignments = runMinimap2(reference_genome, query_genomes, "reference_genome", "query_genomes", option, fastConfig, true);
    log_step_time("Minimap2 Completed");

    // 4. 統計方向 & 尋找最佳定錨點 (Hybrid Approach)
    struct StrandStats {
        long long forwardLength = 0;
        long long reverseLength = 0;
    };
    std::map<std::string, StrandStats> qryStatsMap;
    for (const auto& aln : allAlignments) {
        if (!aln.valid) continue; 
        if (aln.inverse) {
            qryStatsMap[aln.qryName].reverseLength += aln.alnLength;
        } else {
            qryStatsMap[aln.qryName].forwardLength += aln.alnLength;
        }
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
    
    // 備份原始 sequence (轉向與旋轉前)
    std::unordered_map<std::string, std::string> original_sequences = this->sequences;

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
                      << " vs (Reverse: " << qryStatsMap[qryName].reverseLength << " bp)"
                      << " | Shifted: " << true_shift << " bp\n";
        }

        // 寫回原本的資料結構
        BlockSet* qryBlockSet = pair.second.get();
        auto first_block = qryBlockSet->getAllBlocks()[0];
        if (auto sp = first_block.lock()) {
            sp->setConsensus(Consensus(seq)); 
        }
        qryBlockSet->rebuildAllPointers();
        this->sequences[qryName] = seq;
    }

    if (debug) std::cout << "[DEBUG] Inverted: " << rcCount << " | Shifted: " << shiftCount << " genomes.\n";
    log_step_time("Transformations Applied");

    // ==========================================================
    // 🌟 6. 過濾序列 (filter == true 時僅保留與 ref 有明確 strand 且覆蓋率 > 90% 之序列) 並根據 target_num 輸出 FASTA
    // ==========================================================

    long long refLen = sequence_lengths[bestRefName];

    // 收集符合條件的 sequence 名稱 (Reference 本身為第一筆)
    std::vector<std::string> qualified_names;
    qualified_names.push_back(bestRefName);

    for (const auto& pair : this->sequences) {
        std::string qryName = pair.first;
        if (qryName == bestRefName) continue;

        if (filter) {
            long long f_len = qryStatsMap[qryName].forwardLength;
            long long r_len = qryStatsMap[qryName].reverseLength;

            bool pass = (refLen > 0) && (((double)f_len / refLen > 0.9) || ((double)r_len / refLen > 0.9));
            if (pass) {
                qualified_names.push_back(qryName);
            } else if (debug) {
                std::cout << "[DEBUG] [FILTERED OUT] " << qryName 
                          << " (Forward: " << f_len << " bp, Reverse: " << r_len << " bp, RefLen: " << refLen << " bp)\n";
            }
        } else {
            qualified_names.push_back(qryName);
        }
    }

    if (debug) {
        std::cout << "[DEBUG] Filter: " << (filter ? "ON" : "OFF") 
                  << " | Total sequences: " << this->sequences.size() 
                  << " | Qualified sequences: " << qualified_names.size() << "\n";
    }

    // 依據 target_num 分別截取前 N 筆序列寫出 original 與 orientated FASTA (僅在 filter 開啟時才寫入檔案)
    if (filter) {
        std::string outDir = option.tempDir.empty() ? "." : option.tempDir;
        for (int N : target_num) {
            size_t count = std::min(static_cast<size_t>(N), qualified_names.size());
            StringPairs orig_seqs;
            StringPairs orient_seqs;
            orig_seqs.reserve(count);
            orient_seqs.reserve(count);

            for (size_t i = 0; i < count; ++i) {
                const std::string& name = qualified_names[i];
                orig_seqs.push_back({name, original_sequences[name]});
                orient_seqs.push_back({name, this->sequences[name]});
            }

            std::string origFileName = outDir + "/seq_" + std::to_string(N) + ".original.fa";
            std::string orientFileName = outDir + "/seq_" + std::to_string(N) + ".orientated.fa";

            bool isCompressed = false;
            bool appendMode = false;

            mga::io::writeAlignment(origFileName, orig_seqs, isCompressed, appendMode);
            mga::io::writeAlignment(orientFileName, orient_seqs, isCompressed, appendMode);

            if (debug) {
                std::cout << "[DEBUG] Target N=" << N 
                          << ": Wrote " << orig_seqs.size() << " sequences to '" << origFileName << "' and '" << orientFileName << "'\n";
            }
        }
        log_step_time("Output FASTA Files Written");
    }

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
