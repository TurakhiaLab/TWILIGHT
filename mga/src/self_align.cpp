#include "block_set.hpp"
#include "alignment.hpp"
#include "minimap2_util.hpp"
#include "minimap2_config.hpp"
#include "cigar_util.hpp"
#include "option.hpp"
#include "coordinate_manager.hpp"
#include "global_alignment.hpp"

#include <boost/mpl/minus.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>

// 輔助函式：根據累積的 Reference Insertions 調整後續 Query 的 CIGAR
static CigarString adjustCigarForRefInsertions(const CigarString &origCigar, 
                                               const std::map<int, int> &ref_insertions, 
                                               int startRefPos = 0) 
{
    if (ref_insertions.empty()) return origCigar;

    CigarString newCigar;
    int currRefPos = startRefPos;
    std::unordered_set<int> appliedInsertions;

    for (const auto &op : origCigar) {
        int len = op.first;
        char type = op.second;

        if (consumesRef(type)) {
            int remainingLen = len;

            while (remainingLen > 0) {
                if (appliedInsertions.find(currRefPos) == appliedInsertions.end()) {
                    auto it = ref_insertions.find(currRefPos);
                    if (it != ref_insertions.end()) {
                        newCigar.push_back({ it->second, 'D' });
                        appliedInsertions.insert(currRefPos);
                    }
                }

                auto itNext = ref_insertions.upper_bound(currRefPos);
                int nextInsPos = (itNext != ref_insertions.end()) ? itNext->first : -1;

                int step = remainingLen;
                if (nextInsPos != -1 && nextInsPos < currRefPos + remainingLen) {
                    step = nextInsPos - currRefPos;
                }

                newCigar.push_back({ step, type });
                currRefPos += step;
                remainingLen -= step;
            }

            if (appliedInsertions.find(currRefPos) == appliedInsertions.end()) {
                auto it = ref_insertions.find(currRefPos);
                if (it != ref_insertions.end()) {
                    newCigar.push_back({ it->second, 'D' });
                    appliedInsertions.insert(currRefPos);
                }
            }
        } else {
            if (appliedInsertions.find(currRefPos) == appliedInsertions.end()) {
                auto it = ref_insertions.find(currRefPos);
                if (it != ref_insertions.end()) {
                    newCigar.push_back({ it->second, 'D' });
                    appliedInsertions.insert(currRefPos);
                }
            }
            newCigar.push_back({ len, type });
        }
    }

    return compressCigar(newCigar);
}

Alignments BlockSet::selfAlign(Option &option, CoordinateManager *coordMgr, double min_ratio) {
    bool DEBUG_MODE = false;

    // 1. 分類 Blocks：
    // - distant_blocks (Queries)：完全 distant 的 block (所有 copy 均為 distant)
    // - ref_blocks (References)：至少有一個 copy 不是 distant 的 block
    uint64_t total_distant_len = 0;
    std::vector<BlockPtr> distant_blocks;
    std::vector<BlockPtr> ref_blocks;

    std::unordered_map<std::string, BlockPtr> name_to_block;

    for (const auto &pair : blocks) {
        BlockPtr blk = pair.second;
        if (!blk) continue;
        std::string blkName = "Block_" + std::to_string(blk->getId());
        name_to_block[blkName] = blk;

        if (blk->isAllDistant()) {
            uint64_t len = blk->getConsensus().length();
            total_distant_len += len;
            distant_blocks.push_back(blk);
        } else {
            ref_blocks.push_back(blk);
        }
    }

    // 2. 取得 reference (非 distant) blocks 長度總和用於門檻確認 (不需依賴 Ancestral Cache)
    size_t anc_len = 0;
    for (const auto &blk : ref_blocks) {
        if (blk) anc_len += blk->getConsensus().length();
    }

    if (DEBUG_MODE) {
        std::cout << "\n=================== [SelfAlign DEBUG MODE] ===================\n";
        std::cout << "[DEBUG SelfAlign] Before Alignment:\n";
        std::cout << "  - Ancestral sequence length : " << anc_len << " bp\n";
        std::cout << "  - Distant blocks count      : " << distant_blocks.size() << " blocks\n";
        std::cout << "  - Distant blocks total len  : " << total_distant_len << " bp\n";
        std::cout << "  - Reference blocks count    : " << ref_blocks.size() << " blocks\n";
    }

    if (anc_len == 0 || distant_blocks.empty() || ref_blocks.empty()) {
        if (option.verbose || DEBUG_MODE) {
            std::cout << "[SelfAlign] Ancestral sequence is empty, or no distant/reference blocks found. Skipping.\n";
            std::cout << "=============================================================\n\n";
        }
        return {};
    }

    min_ratio = 0;
    // 3. 確認 distant block 總長度是否有 >= ancestral sequence 的 min_ratio (預設 10%)
    if (total_distant_len < static_cast<uint64_t>(anc_len * min_ratio)) {
        if (option.verbose || DEBUG_MODE) {
            std::cout << "[SelfAlign] Distant blocks total length (" << total_distant_len
                      << " bp) < " << (min_ratio * 100.0) << "% of ancestral sequence length ("
                      << anc_len << " bp). Skipping self alignment.\n";
            std::cout << "=============================================================\n\n";
        }
        return {};
    }

    if (option.verbose || DEBUG_MODE) {
        std::cout << "[SelfAlign] Distant blocks total length (" << total_distant_len
                  << " bp) >= " << (min_ratio * 100.0) << "% of ancestral sequence length ("
                  << anc_len << " bp). Proceeding with self alignment ("
                  << distant_blocks.size() << " queries vs " << ref_blocks.size() << " references)...\n";
    }

    // 4. 準備 reference (非全 distant 的 blocks) 與 query (全 distant 的 blocks) 序列
    std::vector<std::string> ref_seqs;
    ref_seqs.reserve(ref_blocks.size());
    SequenceRefs ref;
    ref.reserve(ref_blocks.size());

    for (const auto &blk : ref_blocks) {
        ref_seqs.push_back(blk->getConsensus().getConsensusString());
    }
    for (size_t i = 0; i < ref_blocks.size(); ++i) {
        std::string r_name = "Block_" + std::to_string(ref_blocks[i]->getId());
        ref.push_back({ r_name, ref_seqs[i] });
    }

    std::vector<std::string> query_seqs;
    query_seqs.reserve(distant_blocks.size());
    SequenceRefs qry;
    qry.reserve(distant_blocks.size());

    for (const auto &blk : distant_blocks) {
        query_seqs.push_back(blk->getConsensus().getConsensusString());
    }
    for (size_t i = 0; i < distant_blocks.size(); ++i) {
        std::string q_name = "Block_" + std::to_string(distant_blocks[i]->getId());
        qry.push_back({ q_name, query_seqs[i] });
    }

    // 5. 設定 minimap2 參數
    Minimap2Config config("asm5", true);
    config.setBestN(50);      // 允許搜尋多個候選 reference block
    config.setPriRatio(0.0f); // -p 0

    // 6. 執行 minimap2 alignment
    Alignments raw_alignments = runMinimap2(ref, qry, "NonDistantBlocks_" + ID, "DistantBlocks", option, config);

    // 7. 過濾結果：每個 query 選擇 score 最高的 alignment 留下來
    std::unordered_map<std::string, Alignment> best_alignment_per_query;

    for (const auto &aln : raw_alignments) {
        auto it = best_alignment_per_query.find(aln.qryName);
        if (it == best_alignment_per_query.end()) {
            best_alignment_per_query[aln.qryName] = aln;
        } else {
            if (aln.alnScore > it->second.alnScore ||
               (aln.alnScore == it->second.alnScore && aln.alnLength > it->second.alnLength)) {
                it->second = aln;
            }
        }
    }

    Alignments final_alignments;
    final_alignments.reserve(best_alignment_per_query.size());
    for (const auto &pair : best_alignment_per_query) {
        final_alignments.push_back(pair.second);
    }

    if (DEBUG_MODE) {
        std::cout << "\n[DEBUG SelfAlign] Query to Reference Alignments (" << final_alignments.size() << " query blocks mapped):\n";
        for (const auto &aln : final_alignments) {
            auto qIt = name_to_block.find(aln.qryName);
            auto rIt = name_to_block.find(aln.refName);
            size_t qLen = (qIt != name_to_block.end() && qIt->second) ? qIt->second->getConsensus().length() : 0;
            size_t rLen = (rIt != name_to_block.end() && rIt->second) ? rIt->second->getConsensus().length() : 0;

            std::cout << "  - " << aln.qryName << " (len: " << qLen << ") -> " << aln.refName << " (len: " << rLen << ")"
                      << " | Score: " << aln.alnScore
                      << " | AlnLen: " << aln.alnLength
                      << " | Inverse: " << (aln.inverse ? "true" : "false")
                      << " | CIGAR: " << cigarToStr(aln.CIGAR) << "\n";
        }
    }

    // 8. 依據 Reference Block 群組 Alignments
    std::unordered_map<std::string, std::vector<Alignment>> alignments_by_ref;
    for (const auto &aln : final_alignments) {
        alignments_by_ref[aln.refName].push_back(aln);
    }

    // 9. 針對每個 Reference Block，只選擇 Align Score 最高的單一 Query Block 進行 merge (Rank 1)
    int total_merged_count = 0;
    std::unordered_set<std::string> merged_queries;

    for (auto &pair : alignments_by_ref) {
        const std::string &refName = pair.first;
        auto &group = pair.second;

        // 依據 alnScore 降序排序 (Score 越高越優先 merge，長度作為同分比對)
        std::sort(group.begin(), group.end(), [](const Alignment &a, const Alignment &b) {
            if (a.alnScore == b.alnScore) return a.alnLength > b.alnLength;
            return a.alnScore > b.alnScore;
        });

        auto refIt = name_to_block.find(refName);
        if (refIt == name_to_block.end() || !refIt->second) continue;

        BlockPtr currRefBlock = refIt->second;

        // 每個 Reference Block 在本輪中僅 Merge 分數最高的單一 Query Block (Rank 1)
        for (size_t i = 0; i < std::min<size_t>(1, group.size()); ++i) {
            const Alignment &aln = group[i];
            if (merged_queries.count(aln.qryName)) continue;

            auto qryIt = name_to_block.find(aln.qryName);
            if (qryIt == name_to_block.end() || !qryIt->second) continue;

            BlockPtr qryBlock = qryIt->second;

            int origRefLen = currRefBlock ? (int)currRefBlock->getConsensus().length() : 0;
            int origQryLen = qryBlock ? (int)qryBlock->getConsensus().length() : 0;

            int rStart = aln.refIdx.first;
            int rEnd = aln.refIdx.second;

            int qStart = aln.qryIdx.first;
            int qEnd = aln.qryIdx.second;

            int padRHead = rStart;
            int padRTail = std::max(0, origRefLen - rEnd);

            int padQHead = aln.inverse ? std::max(0, origQryLen - qEnd) : qStart;
            int padQTail = aln.inverse ? qStart : std::max(0, origQryLen - qEnd);

            bool cutRHead = (padRHead > 100);
            bool cutRTail = (padRTail > 100);

            bool cutQHead = (padQHead > 100);
            bool cutQTail = (padQTail > 100);

            // ==========================================
            // A. Reference Block 切割 (> 100 bp)
            // ==========================================
            // 先切 Tail（避免切 Head 導致 Tail 的相對座標位移）
            if (cutRTail && rEnd > 0 && currRefBlock && rEnd < (int)currRefBlock->getConsensus().length()) {
                BlockID oldRefId = currRefBlock->getId();
                auto splitPair = this->splitSingleBlock(oldRefId, rEnd);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldRefId, splitPair.first, splitPair.second, rEnd);
                    }
                    currRefBlock = this->getBlock(splitPair.first);
                    padRTail = 0;
                }
            }

            if (cutRHead && rStart > 0 && currRefBlock && rStart < (int)currRefBlock->getConsensus().length()) {
                BlockID oldRefId = currRefBlock->getId();
                auto splitPair = this->splitSingleBlock(oldRefId, rStart);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldRefId, splitPair.first, splitPair.second, rStart);
                    }
                    currRefBlock = this->getBlock(splitPair.second);
                    padRHead = 0;
                }
            }

            // ==========================================
            // B. Query Block 切割 (> 100 bp) 且繼承 distant 屬性
            // ==========================================
            bool qryIsAllDistant = qryBlock ? qryBlock->isAllDistant() : false;
            auto preserveDistant = [&](BlockPtr blk) {
                if (blk && qryIsAllDistant) {
                    for (int c = 0; c <= blk->getMaxCopy(); ++c) {
                        blk->setDistant(true, c);
                    }
                }
            };

            bool cutQ3Prime = aln.inverse ? cutQHead : cutQTail;
            if (cutQ3Prime && qEnd > 0 && qryBlock && qEnd < (int)qryBlock->getConsensus().length()) {
                BlockID oldQryId = qryBlock->getId();
                auto splitPair = this->splitSingleBlock(oldQryId, qEnd);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldQryId, splitPair.first, splitPair.second, qEnd);
                    }
                    auto leftQ = this->getBlock(splitPair.first);
                    auto rightQ = this->getBlock(splitPair.second);
                    preserveDistant(leftQ);
                    preserveDistant(rightQ);
                    qryBlock = leftQ;
                    if (aln.inverse) padQHead = 0;
                    else padQTail = 0;
                }
            }

            bool cutQ5Prime = aln.inverse ? cutQTail : cutQHead;
            if (cutQ5Prime && qStart > 0 && qryBlock && qStart < (int)qryBlock->getConsensus().length()) {
                BlockID oldQryId = qryBlock->getId();
                auto splitPair = this->splitSingleBlock(oldQryId, qStart);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldQryId, splitPair.first, splitPair.second, qStart);
                    }
                    auto leftQ = this->getBlock(splitPair.first);
                    auto rightQ = this->getBlock(splitPair.second);
                    preserveDistant(leftQ);
                    preserveDistant(rightQ);
                    qryBlock = rightQ;
                    if (aln.inverse) padQTail = 0;
                    else padQHead = 0;
                }
            }

            // ==========================================
            // C. 頭尾剩餘 <= 100 bp 的 Overhang：Linear Gap Tiling 動態對齊
            // ==========================================
            CigarString headCigar;
            if (padRHead > 0 || padQHead > 0) {
                std::string rHeadSeq = (padRHead > 0 && currRefBlock) ? 
                    currRefBlock->getConsensus().getConsensusString().substr(0, padRHead) : "";
                std::string qHeadSeq = "";
                if (padQHead > 0 && qryBlock) {
                    std::string fullQ = qryBlock->getConsensus().getConsensusString();
                    if (!aln.inverse) {
                        qHeadSeq = fullQ.substr(0, padQHead);
                    } else {
                        std::string sub = fullQ.substr(fullQ.length() - padQHead);
                        qHeadSeq = getReverseComplement(sub);
                    }
                }

                if (!rHeadSeq.empty() && !qHeadSeq.empty()) {
                    headCigar = runTilingAlignmentLinearGap(rHeadSeq, qHeadSeq);
                } else if (!rHeadSeq.empty()) {
                    headCigar = {{ (int)rHeadSeq.length(), 'D' }};
                } else if (!qHeadSeq.empty()) {
                    headCigar = {{ (int)qHeadSeq.length(), 'I' }};
                }
            }

            CigarString tailCigar;
            if (padRTail > 0 || padQTail > 0) {
                std::string rTailSeq = "";
                if (padRTail > 0 && currRefBlock) {
                    std::string fullR = currRefBlock->getConsensus().getConsensusString();
                    rTailSeq = fullR.substr(fullR.length() - padRTail);
                }
                std::string qTailSeq = "";
                if (padQTail > 0 && qryBlock) {
                    std::string fullQ = qryBlock->getConsensus().getConsensusString();
                    if (!aln.inverse) {
                        qTailSeq = fullQ.substr(fullQ.length() - padQTail);
                    } else {
                        std::string sub = fullQ.substr(0, padQTail);
                        qTailSeq = getReverseComplement(sub);
                    }
                }

                if (!rTailSeq.empty() && !qTailSeq.empty()) {
                    tailCigar = runTilingAlignmentLinearGap(rTailSeq, qTailSeq);
                } else if (!rTailSeq.empty()) {
                    tailCigar = {{ (int)rTailSeq.length(), 'D' }};
                } else if (!qTailSeq.empty()) {
                    tailCigar = {{ (int)qTailSeq.length(), 'I' }};
                }
            }

            CigarString paddedCigar;
            for (const auto &op : headCigar) paddedCigar.push_back(op);
            for (const auto &op : aln.CIGAR) paddedCigar.push_back(op);
            for (const auto &op : tailCigar) paddedCigar.push_back(op);
            paddedCigar = compressCigar(paddedCigar);

            if (DEBUG_MODE) {
                size_t qLen = qryBlock ? qryBlock->getConsensus().length() : 0;
                size_t rLen = currRefBlock ? currRefBlock->getConsensus().length() : 0;
                std::cout << "  [DEBUG Merge] Merging " << aln.qryName << " (len: " << qLen << ") into " << refName << " (len: " << rLen << ")"
                          << " [CIGAR Padded: " << cigarToStr(paddedCigar) << "]\n";
            }

            // 執行 mergeTwoBlocks (Mode 3: Scenario B, Ref 已存在而 Query 為新 Block)
            BlockPtr nextRefBlock = mergeTwoBlocks(currRefBlock, qryBlock, paddedCigar, aln.inverse, 3);
            if (nextRefBlock) {
                currRefBlock = nextRefBlock;
                std::string newRefName = "Block_" + std::to_string(nextRefBlock->getId());
                name_to_block[refName] = nextRefBlock;
                name_to_block[newRefName] = nextRefBlock;
                merged_queries.insert(aln.qryName);
                total_merged_count++;
            }
        }
    }

    if (DEBUG_MODE) {
        int after_distant_count = 0;
        uint64_t after_distant_len = 0;
        for (const auto &pair : blocks) {
            BlockPtr blk = pair.second;
            if (blk && blk->isAllDistant()) {
                after_distant_count++;
                after_distant_len += blk->getConsensus().length();
            }
        }

        std::cout << "\n[DEBUG SelfAlign] After Alignment & Merge:\n";
        std::cout << "  - Merged query blocks count : " << total_merged_count << " blocks\n";
        std::cout << "  - Remaining distant blocks  : " << after_distant_count << " blocks (was " << distant_blocks.size() << ")\n";
        std::cout << "  - Remaining distant length  : " << after_distant_len << " bp (was " << total_distant_len << " bp)\n";
        std::cout << "=============================================================\n\n";
    }

    return final_alignments;
}

Alignments BlockSet::selfAlignDistant(Option &option, CoordinateManager *coordMgr, double min_ratio) {
    bool DEBUG_MODE = false;

    // 1. 分類出所有 distant blocks (Queries & References 均為 distant blocks)
    uint64_t total_distant_len = 0;
    std::vector<BlockPtr> distant_blocks;
    std::unordered_map<std::string, BlockPtr> name_to_block;

    for (const auto &pair : blocks) {
        BlockPtr blk = pair.second;
        if (!blk) continue;
        std::string blkName = "Block_" + std::to_string(blk->getId());
        name_to_block[blkName] = blk;

        if (blk->isAllDistant()) {
            uint64_t len = blk->getConsensus().length();
            total_distant_len += len;
            distant_blocks.push_back(blk);
        }
    }

    if (DEBUG_MODE) {
        std::cout << "\n=================== [SelfAlignDistant DEBUG MODE] ===================\n";
        std::cout << "[DEBUG SelfAlignDistant] Before Alignment:\n";
        std::cout << "  - Distant blocks count     : " << distant_blocks.size() << " blocks\n";
        std::cout << "  - Distant blocks total len : " << total_distant_len << " bp\n";
    }

    if (distant_blocks.size() < 2) {
        if (option.verbose || DEBUG_MODE) {
            std::cout << "[SelfAlignDistant] Less than 2 distant blocks found. Skipping.\n";
            std::cout << "=============================================================\n\n";
        }
        return {};
    }

    // 2. 準備 distant vs distant 的 minimap2 序列
    std::vector<std::string> distant_seqs;
    distant_seqs.reserve(distant_blocks.size());
    SequenceRefs refs;
    refs.reserve(distant_blocks.size());

    for (const auto &blk : distant_blocks) {
        distant_seqs.push_back(blk->getConsensus().getConsensusString());
    }
    for (size_t i = 0; i < distant_blocks.size(); ++i) {
        std::string r_name = "Block_" + std::to_string(distant_blocks[i]->getId());
        refs.push_back({ r_name, distant_seqs[i] });
    }

    // 3. 設定 minimap2 參數
    Minimap2Config config("asm5", true);
    config.setBestN(50);
    config.setPriRatio(0.0f);

    // 4. 執行 minimap2 alignment (distant blocks vs distant blocks)
    Alignments raw_alignments = runMinimap2(refs, refs, "DistantBlocks_Ref", "DistantBlocks_Qry", option, config);

    // 5. 使用 DSU (Disjoint Set Union) 根據 Coverage >= 90% 將 Distant Blocks 分組 (Group / Cluster)
    struct DSU {
        std::unordered_map<std::string, std::string> parent;
        std::string find(const std::string &i) {
            if (parent.find(i) == parent.end()) parent[i] = i;
            if (parent[i] == i) return i;
            return parent[i] = find(parent[i]);
        }
        void unite(const std::string &i, const std::string &j) {
            std::string rootI = find(i);
            std::string rootJ = find(j);
            if (rootI != rootJ) parent[rootI] = rootJ;
        }
    } dsu;

    int candidate_edge_count = 0;
    for (const auto &aln : raw_alignments) {
        if (aln.refName == aln.qryName) continue; // 排除自比自
        if (aln.alnLength < 100) continue; // 排除 alignment 長度 < 100 bp

        auto refIt = name_to_block.find(aln.refName);
        auto qryIt = name_to_block.find(aln.qryName);
        if (refIt == name_to_block.end() || !refIt->second) continue;
        if (qryIt == name_to_block.end() || !qryIt->second) continue;

        double rLen = refIt->second->getConsensus().length();
        double qLen = qryIt->second->getConsensus().length();
        if (rLen <= 0 || qLen <= 0) continue;

        double covR = (double)(aln.refIdx.second - aln.refIdx.first) / rLen;
        double covQ = (double)(aln.qryIdx.second - aln.qryIdx.first) / qLen;

        if (DEBUG_MODE) {
            std::cout << "  [DEBUG SelfAlignDistant Candidate] " << aln.qryName << " (len:" << (int)qLen << ") -> "
                      << aln.refName << " (len:" << (int)rLen << ") | alnLen:" << aln.alnLength
                      << " | covR:" << (covR * 100.0) << "% | covQ:" << (covQ * 100.0) << "%\n";
        }

        // 放寬 Coverage 條件：只要任一邊對齊覆蓋率 >= 50% 且 alnLength >= 100 bp 即可劃入同一 Group
        if (covR >= 0.50 || covQ >= 0.50) {
            dsu.unite(aln.refName, aln.qryName);
            candidate_edge_count++;
        }
    }

    // 6. 收集各大群組 (Clusters)
    std::unordered_map<std::string, std::vector<BlockPtr>> cluster_map;
    for (const auto &blk : distant_blocks) {
        std::string blkName = "Block_" + std::to_string(blk->getId());
        std::string root = dsu.find(blkName);
        cluster_map[root].push_back(blk);
    }

    if (DEBUG_MODE) {
        int valid_clusters = 0;
        for (const auto &pair : cluster_map) {
            if (pair.second.size() >= 2) valid_clusters++;
        }
        std::cout << "[DEBUG SelfAlignDistant] Accepted " << candidate_edge_count << " alignment edges. Formed " 
                  << valid_clusters << " clusters (with size >= 2) out of "
                  << cluster_map.size() << " total groups.\n";
    }

    // 7. 針對每個 Cluster 的 Blocks 進行一條一條對齊與合併 (Sequential Merge)
    int total_merged_count = 0;
    Alignments final_selected_alignments;

    for (auto &pair : cluster_map) {
        auto &cluster = pair.second;
        if (cluster.size() < 2) continue;

        // 依據 consensus 長度降序排序，選擇最長者作為初始 anchor reference block
        std::sort(cluster.begin(), cluster.end(), [](const BlockPtr &a, const BlockPtr &b) {
            return a->getConsensus().length() > b->getConsensus().length();
        });

        BlockPtr currRefBlock = cluster[0];

        for (size_t i = 1; i < cluster.size(); ++i) {
            BlockPtr qryBlock = cluster[i];
            if (!currRefBlock || !qryBlock) continue;

            std::string r_name = "Block_" + std::to_string(currRefBlock->getId());
            std::string q_name = "Block_" + std::to_string(qryBlock->getId());
            std::string r_str = currRefBlock->getConsensus().getConsensusString();
            std::string q_str = qryBlock->getConsensus().getConsensusString();

            SequenceRefs single_ref = {{ r_name, r_str }};
            SequenceRefs single_qry = {{ q_name, q_str }};

            Alignments pair_alns = runMinimap2(single_ref, single_qry, "DistantGroup_Ref", "DistantGroup_Qry", option, config);
            if (pair_alns.empty()) continue;

            const Alignment &aln = pair_alns[0];
            final_selected_alignments.push_back(aln);

            int origRefLen = (int)currRefBlock->getConsensus().length();
            int origQryLen = (int)qryBlock->getConsensus().length();

            int rStart = aln.refIdx.first;
            int rEnd = aln.refIdx.second;

            int qStart = aln.qryIdx.first;
            int qEnd = aln.qryIdx.second;

            int padRHead = rStart;
            int padRTail = std::max(0, origRefLen - rEnd);

            int padQHead = aln.inverse ? std::max(0, origQryLen - qEnd) : qStart;
            int padQTail = aln.inverse ? qStart : std::max(0, origQryLen - qEnd);

            bool cutRHead = (padRHead > 100);
            bool cutRTail = (padRTail > 100);

            bool cutQHead = (padQHead > 100);
            bool cutQTail = (padQTail > 100);

            // A. Reference Block 切割 (> 100 bp)
            bool refIsAllDistant = currRefBlock ? currRefBlock->isAllDistant() : false;
            auto preserveDistantRef = [&](BlockPtr blk) {
                if (blk && refIsAllDistant) {
                    for (int c = 0; c <= blk->getMaxCopy(); ++c) {
                        blk->setDistant(true, c);
                    }
                }
            };

            if (cutRTail && rEnd > 0 && currRefBlock && rEnd < (int)currRefBlock->getConsensus().length()) {
                BlockID oldRefId = currRefBlock->getId();
                auto splitPair = this->splitSingleBlock(oldRefId, rEnd);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldRefId, splitPair.first, splitPair.second, rEnd);
                    }
                    auto leftR = this->getBlock(splitPair.first);
                    auto rightR = this->getBlock(splitPair.second);
                    preserveDistantRef(leftR);
                    preserveDistantRef(rightR);
                    currRefBlock = leftR;
                    padRTail = 0;
                }
            }

            if (cutRHead && rStart > 0 && currRefBlock && rStart < (int)currRefBlock->getConsensus().length()) {
                BlockID oldRefId = currRefBlock->getId();
                auto splitPair = this->splitSingleBlock(oldRefId, rStart);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldRefId, splitPair.first, splitPair.second, rStart);
                    }
                    auto leftR = this->getBlock(splitPair.first);
                    auto rightR = this->getBlock(splitPair.second);
                    preserveDistantRef(leftR);
                    preserveDistantRef(rightR);
                    currRefBlock = rightR;
                    padRHead = 0;
                }
            }

            // B. Query Block 切割 (> 100 bp)
            bool qryIsAllDistant = qryBlock ? qryBlock->isAllDistant() : false;
            auto preserveDistantQry = [&](BlockPtr blk) {
                if (blk && qryIsAllDistant) {
                    for (int c = 0; c <= blk->getMaxCopy(); ++c) {
                        blk->setDistant(true, c);
                    }
                }
            };

            bool cutQ3Prime = aln.inverse ? cutQHead : cutQTail;
            if (cutQ3Prime && qEnd > 0 && qryBlock && qEnd < (int)qryBlock->getConsensus().length()) {
                BlockID oldQryId = qryBlock->getId();
                auto splitPair = this->splitSingleBlock(oldQryId, qEnd);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldQryId, splitPair.first, splitPair.second, qEnd);
                    }
                    auto leftQ = this->getBlock(splitPair.first);
                    auto rightQ = this->getBlock(splitPair.second);
                    preserveDistantQry(leftQ);
                    preserveDistantQry(rightQ);
                    qryBlock = leftQ;
                    if (aln.inverse) padQHead = 0;
                    else padQTail = 0;
                }
            }

            bool cutQ5Prime = aln.inverse ? cutQTail : cutQHead;
            if (cutQ5Prime && qStart > 0 && qryBlock && qStart < (int)qryBlock->getConsensus().length()) {
                BlockID oldQryId = qryBlock->getId();
                auto splitPair = this->splitSingleBlock(oldQryId, qStart);
                if (splitPair.first != (uint64_t)-1 && splitPair.second != (uint64_t)-1) {
                    if (coordMgr) {
                        coordMgr->updateAfterSplit(oldQryId, splitPair.first, splitPair.second, qStart);
                    }
                    auto leftQ = this->getBlock(splitPair.first);
                    auto rightQ = this->getBlock(splitPair.second);
                    preserveDistantQry(leftQ);
                    preserveDistantQry(rightQ);
                    qryBlock = rightQ;
                    if (aln.inverse) padQTail = 0;
                    else padQHead = 0;
                }
            }

            // C. 頭尾剩餘 <= 100 bp 的 Overhang：Linear Gap Tiling 動態對齊
            CigarString headCigar;
            if (padRHead > 0 || padQHead > 0) {
                std::string rHeadSeq = (padRHead > 0 && currRefBlock) ? 
                    currRefBlock->getConsensus().getConsensusString().substr(0, padRHead) : "";
                std::string qHeadSeq = "";
                if (padQHead > 0 && qryBlock) {
                    std::string fullQ = qryBlock->getConsensus().getConsensusString();
                    if (!aln.inverse) {
                        qHeadSeq = fullQ.substr(0, padQHead);
                    } else {
                        std::string sub = fullQ.substr(fullQ.length() - padQHead);
                        qHeadSeq = getReverseComplement(sub);
                    }
                }

                if (!rHeadSeq.empty() && !qHeadSeq.empty()) {
                    headCigar = runTilingAlignmentLinearGap(rHeadSeq, qHeadSeq);
                } else if (!rHeadSeq.empty()) {
                    headCigar = {{ (int)rHeadSeq.length(), 'D' }};
                } else if (!qHeadSeq.empty()) {
                    headCigar = {{ (int)qHeadSeq.length(), 'I' }};
                }
            }

            CigarString tailCigar;
            if (padRTail > 0 || padQTail > 0) {
                std::string rTailSeq = "";
                if (padRTail > 0 && currRefBlock) {
                    std::string fullR = currRefBlock->getConsensus().getConsensusString();
                    rTailSeq = fullR.substr(fullR.length() - padRTail);
                }
                std::string qTailSeq = "";
                if (padQTail > 0 && qryBlock) {
                    std::string fullQ = qryBlock->getConsensus().getConsensusString();
                    if (!aln.inverse) {
                        qTailSeq = fullQ.substr(fullQ.length() - padQTail);
                    } else {
                        std::string sub = fullQ.substr(0, padQTail);
                        qTailSeq = getReverseComplement(sub);
                    }
                }

                if (!rTailSeq.empty() && !qTailSeq.empty()) {
                    tailCigar = runTilingAlignmentLinearGap(rTailSeq, qTailSeq);
                } else if (!rTailSeq.empty()) {
                    tailCigar = {{ (int)rTailSeq.length(), 'D' }};
                } else if (!qTailSeq.empty()) {
                    tailCigar = {{ (int)qTailSeq.length(), 'I' }};
                }
            }

            CigarString paddedCigar;
            for (const auto &op : headCigar) paddedCigar.push_back(op);
            for (const auto &op : aln.CIGAR) paddedCigar.push_back(op);
            for (const auto &op : tailCigar) paddedCigar.push_back(op);
            paddedCigar = compressCigar(paddedCigar);

            if (DEBUG_MODE) {
                size_t qLen = qryBlock ? qryBlock->getConsensus().length() : 0;
                size_t rLen = currRefBlock ? currRefBlock->getConsensus().length() : 0;
                std::cout << "  [DEBUG MergeDistantGroup] Merging " << aln.qryName << " (len: " << qLen << ") into " << aln.refName << " (len: " << rLen << ")"
                          << " [CIGAR Padded: " << cigarToStr(paddedCigar) << "]\n";
            }

            BlockPtr nextRefBlock = mergeTwoBlocks(currRefBlock, qryBlock, paddedCigar, aln.inverse, 3);
            if (nextRefBlock) {
                currRefBlock = nextRefBlock;
                total_merged_count++;
            }
        }
    }

    if (DEBUG_MODE) {
        int after_distant_count = 0;
        uint64_t after_distant_len = 0;
        for (const auto &pair : blocks) {
            BlockPtr blk = pair.second;
            if (blk && blk->isAllDistant()) {
                after_distant_count++;
                after_distant_len += blk->getConsensus().length();
            }
        }

        std::cout << "\n[DEBUG SelfAlignDistant] After Grouping & Sequential Merge:\n";
        std::cout << "  - Total merged distant blocks : " << total_merged_count << " blocks\n";
        std::cout << "  - Remaining distant blocks    : " << after_distant_count << " blocks (was " << distant_blocks.size() << ")\n";
        std::cout << "  - Remaining distant length    : " << after_distant_len << " bp (was " << total_distant_len << " bp)\n";
        std::cout << "=============================================================\n\n";
    }

    return final_selected_alignments;
}
