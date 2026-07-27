#include "mga.hpp"

#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace fs = boost::filesystem;

void BlockSet::selfMapping(Option& option) {

    const int MIN_CHUNK_LEN = 15;

    bool debug = false;

    std::string temp_dir = option.tempDir;
    std::string seqName = blocks_.begin()->second->getSequences().begin()->first;
    std::string seq = blocks_.begin()->second->getConsensus();
    std::string seqFile = temp_dir + "/" +  seqName + ".fa";
    std::string pafFile = temp_dir + "/" +  seqName + "-self.paf";
    mga::stringPairVec seqPair (1, {seqName, seq});
    
    int total_original_length = seq.length();
    
    if (debug) std::cout << "[SelfMapping] 0. Writing and Self-aligning [" << seqName << "] (Len: " << total_original_length << " bp)...\n";
    mga::io::writeAlignment(seqFile, seqPair, false, false);

    std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
    std::string command;
    int system_ret; 

    auto minimap2_start = std::chrono::high_resolution_clock::now();
    command = minimap2_path + " -cx asm5 -g 500 -r 500,500 -N 20 -X -s 100 " + seqFile + " " + seqFile + " > " + pafFile + " 2> /dev/null";
    system_ret = system(command.c_str());
    if (system_ret != 0) {
        std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
    }
    auto minimap2_end = std::chrono::high_resolution_clock::now();
    option.minimap2_time += std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count();


    if (debug) std::cout << "[SelfMapping] 1. Parsing PAF and Deduplicating...\n";
    std::vector<mga::Alignment> rawAlignments = mga::parser::parseMinimap2PAF(pafFile);
    std::vector<mga::Alignment> uniqueAlignments;

    for (size_t i = 0; i < rawAlignments.size(); ++i) {
        if (i % 2 == 0) { 
            if (rawAlignments[i].alnScore > 0) uniqueAlignments.push_back(rawAlignments[i]);
        }
    }

    std::sort(uniqueAlignments.begin(), uniqueAlignments.end(), [](const mga::Alignment& a, const mga::Alignment& b) {
        return a.alnScore > b.alnScore; 
    });

    std::remove(seqFile.c_str());
    std::remove(pafFile.c_str());

    if (debug) std::cout << "[SelfMapping] 2. Initializing Coordinate Dictionary...\n";

    //Key: start_coordinate, Value: Segment info.
    std::map<int, SegNode> dict;

    auto initialBlock = blocks_.begin()->second;
    int seqLen = total_original_length; 
    dict[0] = {0, seqLen, initialBlock->getId()};

    auto splitDictSegment = [&](int cutPos) {
        auto it = dict.upper_bound(cutPos); it--;
        if (it->first == cutPos) return;

        int start = it->second.start;
        Block::ID targetBlkId = it->second.blkId;
        std::shared_ptr<Block> targetBlk = this->getBlock(targetBlkId);
        if (!targetBlk) return;

        Segment* targetSeg = nullptr;
        for (auto& seqPair : targetBlk->getSequences()) {
            if (seqPair.second.getSegments().count(start)) {
                targetSeg = &(seqPair.second.getSegments()[start]);
                break;
            }
        }
        if (!targetSeg) return;

        int origLocalCut = targetSeg->isReverse() ? (targetSeg->getEnd() - cutPos) : (cutPos - targetSeg->getStart());
        int consensusCut = origLocalCut;
        for (auto& v : targetSeg->getVariants()) {
            if (v.getType() == Variation::GAP && v.getStart() <= consensusCut) {
                consensusCut += (v.getEnd() - v.getStart());
            }
        }

        if (consensusCut <= 0 || consensusCut >= targetBlk->getConsensus().length()) return; 

        auto parts = this->splitSingleBlock(targetBlkId, consensusCut);
        if (parts.first == -1) return; 

        this->rebuildDictionary(dict, seqName); 
    };

    if (debug) std::cout << "[SelfMapping] 3. Greedy Interval Merging...\n";

    int total_collapsed_bases = 0;
    int merge_operations_count = 0;

    int alnCount = 0;

    for (const auto& aln : uniqueAlignments) {

        alnCount++;

        int subAlnCount = 0;

        int rCursor = aln.refIdx.first;
        int qCursor = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first;

        int cigarIdx = 0;
        if (aln.CIGAR.empty()) continue;
        int opRemain = aln.CIGAR[0].first;
        char opType = aln.CIGAR[0].second;

        while (rCursor < aln.refIdx.second && 
              (aln.inverse ? qCursor > aln.qryIdx.first : qCursor < aln.qryIdx.second) && 
              cigarIdx < aln.CIGAR.size()) {

            subAlnCount++;
            // if (subAlnCount > 20) exit(1);
            
            // 1. 查字典，找出目前的 rCursor 和 qCursor 屬於哪兩個 Segment 區間
            auto rIt = dict.upper_bound(rCursor); rIt--;
            SegNode rSeg = rIt->second;
            int rRemain = rSeg.end - rCursor; // Ref 永遠向右走，看離右邊界多遠

            SegNode qSeg;
            int qRemain = 0;
            if (!aln.inverse) {
                auto qIt = dict.upper_bound(qCursor); qIt--;
                qSeg = qIt->second;
                qRemain = qSeg.end - qCursor; // 向右走，看離右邊界多遠
            } else {
                // 向左走！要尋找包含 (qCursor - 1) 的區間
                auto qIt = dict.upper_bound(qCursor - 1); qIt--;
                qSeg = qIt->second;
                qRemain = qCursor - qSeg.start; // 向左走，看離左邊界多遠
            }

            // 2. 解析 CIGAR 碎片
            mga::Cigar fragCigar;
            int chunkRefLen = 0;
            int chunkQryLen = 0;

            while (cigarIdx < aln.CIGAR.size()) {
                bool consumesRef = (opType == 'M' || opType == '=' || opType == 'X' || opType == 'D');
                bool consumesQry = (opType == 'M' || opType == '=' || opType == 'X' || opType == 'I');

                int maxStep = opRemain;
                if (consumesRef && maxStep > rRemain - chunkRefLen) maxStep = rRemain - chunkRefLen;
                if (consumesQry && maxStep > qRemain - chunkQryLen) maxStep = qRemain - chunkQryLen;

                if (maxStep == 0) break; // 撞到邊界了

                if (!fragCigar.empty() && fragCigar.back().second == opType) {
                    fragCigar.back().first += maxStep;
                } else {
                    fragCigar.push_back({maxStep, opType});
                }
                
                if (consumesRef) chunkRefLen += maxStep;
                if (consumesQry) chunkQryLen += maxStep;

                opRemain -= maxStep;
                if (opRemain == 0) {
                    cigarIdx++;
                    if (cigarIdx < aln.CIGAR.size()) {
                        opRemain = aln.CIGAR[cigarIdx].first;
                        opType = aln.CIGAR[cigarIdx].second;
                    }
                }
            }

            // ==========================================
            // 【核心修復】：精準還原這段 Chunk 的真實左右邊界
            // ==========================================
            int rChunkStart = rCursor;
            int rChunkEnd = rCursor + chunkRefLen;
            int qChunkStart = aln.inverse ? (qCursor - chunkQryLen) : qCursor;
            int qChunkEnd = aln.inverse ? qCursor : (qCursor + chunkQryLen);

            if (debug) {
                std::cout << "\n------------------------------------------------------------\n";
                std::cout << "[DEBUG] Alignment [" << alnCount << "-" << subAlnCount << "/" << uniqueAlignments.size() << "]\n";
                std::cout << "[DEBUG] Direction: " << (aln.inverse ? "Reverse Complement (-)" : "Forward (+)") << "\n";
                std::cout << "[DEBUG] Ref Raw Chunk: [" << rChunkStart << " -> " << rChunkEnd << ") | Length: " << chunkRefLen << " | Curr Blk ID: " << rSeg.blkId << "\n";
                std::cout << "[DEBUG] Qry Raw Chunk: [" << qChunkStart << " -> " << qChunkEnd << ") | Length: " << chunkQryLen << " | Curr Blk ID: " << qSeg.blkId << "\n";
                std::cout << "[DEBUG] Fragment CIGAR: ";
                for (auto op : fragCigar) std::cout << op.first << op.second;
                std::cout << "\n";
            }

            // 如果已經在同一個 Block (之前 Merge 過了)
            if (rSeg.blkId == qSeg.blkId) {
                bool skip = false;
                if (rChunkStart <= qChunkEnd && qChunkStart <= rChunkEnd) {
                    skip = true;
                } else {
                    int dist = (rChunkStart < qChunkStart) ? (qChunkStart - rChunkEnd) : (rChunkStart - qChunkEnd);
                    if (dist < 100) skip = true;
                }
                if (!skip && rSeg.start != qSeg.start) skip = true;
                
                if (skip) {
                    rCursor = rChunkEnd;
                    qCursor = aln.inverse ? qChunkStart : qChunkEnd; // 更新游標
                    continue;
                }
            }

            // 修改這裡：放寬 Chunk 過濾門檻，拯救被截斷的合法尾巴
            if (chunkRefLen < MIN_CHUNK_LEN || chunkQryLen < MIN_CHUNK_LEN) {
                rCursor = rChunkEnd;
                qCursor = aln.inverse ? qChunkStart : qChunkEnd; 
                continue;
            }

            // 3. 進行【絕對嚴格】的 Overhang 切割與字典更新
            // 不管差幾 bp，只要沒有切齊，就強制切斷！保證後續合併是 100% 滿版合併。
            if (rChunkStart > rSeg.start) splitDictSegment(rChunkStart);
            if (rChunkEnd < rSeg.end) splitDictSegment(rChunkEnd);
            
            // 重新取得切齊後的 Ref SegNode
            auto rIt2 = dict.upper_bound(rChunkStart); rIt2--; 
            rSeg = rIt2->second;

            if (qChunkStart > qSeg.start) splitDictSegment(qChunkStart);
            if (qChunkEnd < qSeg.end) splitDictSegment(qChunkEnd);
            
            // 重新取得切齊後的 Qry SegNode
            auto qIt2 = dict.upper_bound(qChunkStart); qIt2--; 
            qSeg = qIt2->second;

            // 4. 取得原始 Block 與 Segment (此時它們已經是完美切齊的了！)
            std::shared_ptr<Block> refBlk = this->getBlock(rSeg.blkId);
            std::shared_ptr<Block> qryBlk = this->getBlock(qSeg.blkId);

            Segment* targetRefSeg = nullptr;
            for (auto& seqPair : refBlk->getSequences()) {
                if (seqPair.second.getSegments().count(rSeg.start)) {
                    targetRefSeg = &(seqPair.second.getSegments()[rSeg.start]);
                    break;
                }
            }
            Segment* targetQrySeg = nullptr;
            for (auto& seqPair : qryBlk->getSequences()) {
                if (seqPair.second.getSegments().count(qSeg.start)) {
                    targetQrySeg = &(seqPair.second.getSegments()[qSeg.start]);
                    break;
                }
            }

            if (!targetRefSeg || !targetQrySeg) {
                rCursor = rChunkEnd;
                qCursor = aln.inverse ? qChunkStart : qChunkEnd;
                continue;
            }

            // ==========================================
            // 5. 【終極修復】：廢除 Overhang Padding！
            // 既然已經完美切齊，CIGAR 只需要反映對齊區段，不需要補 D 和 I
            // ==========================================
            mga::Cigar paddedCigar = fragCigar; // 直接使用原始的對齊 CIGAR！

            // 轉換為 Consensus Space 的方向處理
            if (targetRefSeg->isReverse()) {
                std::reverse(paddedCigar.begin(), paddedCigar.end());
            }

            // ==========================================
            // 6. 計算真實反轉關係並合併
            // ==========================================
            // effectiveInverse 告訴 mergeTwoBlocks：Qry Consensus 是否需要相對於 Ref Consensus 進行翻轉
            bool effectiveInverse = targetRefSeg->isReverse() ^ aln.inverse ^ targetQrySeg->isReverse();

            mga::Cigar adjustedCigar = adjustCigarWithVariations(
                paddedCigar, 
                *targetRefSeg, 
                *targetQrySeg, 
                effectiveInverse, 
                qryBlk->getConsensus().length()
            );
            
            if (debug) {
                std::cout << "[DEBUG] Adjusted CIGAR (Consensus space): ";
                for (auto op : adjustedCigar) std::cout << op.first << op.second;
                std::cout << "\n";
            }

            auto mergedBlock = this->mergeTwoBlocks(refBlk, qryBlk, adjustedCigar, effectiveInverse);

            // ==========================================
            // === 新增的 Debug 輸出：顯示結果 ===
            // ==========================================
            if (debug) std::cout << "[DEBUG] -> Merged into New Block ID: " << mergedBlock->getId() << " (Consensus Len: " << mergedBlock->getConsensus().length() << ")\n";

            this->updateSegmentLinks(refBlk, mergedBlock);
            this->updateSegmentLinks(qryBlk, mergedBlock);
            this->deleteBlock(refBlk->getId());
            this->deleteBlock(qryBlk->getId());

            this->rebuildDictionary(dict, seqName);

            total_collapsed_bases += chunkRefLen; // 用實際配對長度統計
            merge_operations_count++;

            // 推進下一輪座標
            rCursor = rChunkEnd;
            qCursor = aln.inverse ? qChunkStart : qChunkEnd;
        }
        
    }

    this->absorbMicroBlocks();
    // ==========================================
    // Phase 4: 列印統計數據與最終 Block 內容
    // ==========================================
    double collapsed_ratio = (total_original_length > 0) ? 
        ((double)total_collapsed_bases / total_original_length) * 100.0 : 0.0;

    if (debug) {
        std::cout << "\n[SelfMapping] === Self-Mapping Summary ===\n";
        std::cout << "  - Total Original Length: " << total_original_length << " bp\n";
        std::cout << "  - Duplications Merged:   " << total_collapsed_bases << " bp (" 
                  << std::fixed << std::setprecision(2) << collapsed_ratio << "% collapsed)\n";
        std::cout << "  - Total Merge Operations: " << merge_operations_count << "\n";
        std::cout << "==========================================\n\n";
    }

    if (debug) {
        std::cout << "[SelfMapping] === Final Block Topology & Segments ===\n";
        for (const auto& blockPair : this->blocks_) {
            auto blk = blockPair.second;

            std::cout << "Block ID: " << blk->getId() << "\t| Consensus Len: " << blk->getConsensus().length() << "\n";

            for (auto& seqPair : blk->getSequences()) {
                std::cout << "  ├─ Sequence: " << seqPair.first << "\n";
                for (auto& segPair : seqPair.second.getSegments()) {
                    Segment& seg = segPair.second;
                    std::cout << "  │    └─ Segment Range: [" << seg.getStart() << ", " << seg.getEnd() << "]"
                              << "\t Strand: " << (seg.isReverse() ? "(-)" : "(+)");
                    if (seg.getNextBlock().lock()) std::cout << "    Next: " << seg.getNextBlock().lock()->getId();              
                    std::cout << "\n";
                }
            }
            std::cout << "  └--------------------------------------------------\n";
        }
        std::cout << "\n";
    }
}
