#include "block.hpp"
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>



// =========================================================
// 原始的 Merge (建立新 Block，搬移 Sequences，並刪除舊 Block)
// =========================================================
BlockPtr BlockSet::mergeTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse) 
{
    // 1. 扁平化 Segment 以供 TBB 使用
    std::vector<Segment*> refSegsFlat, qrySegsFlat;
    for (auto& seqPair : refBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) refSegsFlat.push_back(&segPairInner.second);
    }
    for (auto& seqPair : qryBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) qrySegsFlat.push_back(&segPairInner.second);
    }

    // 2. 核心運算：取得 Mappings 與計算出新的 Variations
    MergeMappingData mapData = calculateMergedConsensusAndMappings(refBlock, qryBlock, cigar, inverse, refSegsFlat, qrySegsFlat);
    
    // 【修正】：這個函數現在只會更新相對座標 (Variants)，絕對不會動到 Segment 的 Start/End
    updateSegmentVariations(refSegsFlat, qrySegsFlat, mapData, inverse, 
                            refBlock->getConsensus().length(), qryBlock->getConsensus().length());

    // 3. 建立全新的 Block
    auto mergedBlock = this->createBlock(mapData.mergedConsensus);

    // 4. 將 Sequence 搬移掛載到新的 Block 身上
    auto& mergedSeqs = mergedBlock->getSequences();
    for (auto& seqPair : refBlock->getSequences()) {
        if (mergedSeqs.find(seqPair.first) == mergedSeqs.end()) mergedSeqs[seqPair.first] = Sequence(seqPair.first);
        for (auto& segPairInner : seqPair.second.getSegments()) {
            mergedSeqs[seqPair.first].getSegments()[segPairInner.second.getStart()] = std::move(segPairInner.second);
        }
    }
    for (auto& seqPair : qryBlock->getSequences()) {
        if (mergedSeqs.find(seqPair.first) == mergedSeqs.end()) mergedSeqs[seqPair.first] = Sequence(seqPair.first);
        for (auto& segPairInner : seqPair.second.getSegments()) {
            mergedSeqs[seqPair.first].getSegments()[segPairInner.second.getStart()] = std::move(segPairInner.second);
        }
    }

    // ========================================================
    // 🛠️ 核心修復：解決「產生 3 個 Blocks」的問題
    // 將 Sequences 抽乾後，必須把原本的 refBlock 和 qryBlock 從 BlockSet 刪除
    // ========================================================
    this->deleteBlock(refBlock->getId());
    this->deleteBlock(qryBlock->getId());

    return mergedBlock;
}

// =========================================================
// 全新的 Link (不建立新 Block，就地同步 Consensus 與 Variations)
// =========================================================
void BlockSet::linkTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse) 
{
    std::string refSeq = refBlock->getConsensus();
    std::string qrySeq = qryBlock->getConsensus();

    // ========================================================
    // 🚨 終極防護網：CIGAR 長度嚴格校驗 (防止記憶體越界與 Core Dump)
    // ========================================================
    int cRefLen = 0, cQryLen = 0;
    for (const auto& op : cigar) {
        if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'D') cRefLen += op.first;
        if (op.second == 'M' || op.second == '=' || op.second == 'X' || op.second == 'I') cQryLen += op.first;
    }

    if (cRefLen != refSeq.length() || cQryLen != qrySeq.length()) {
        std::cerr << "\n[CRITICAL ERROR] CIGAR length mismatch in linkTwoBlocks!\n"
                  << "  Ref Block ID: " << refBlock->getId() << " (Len: " << refSeq.length() << ") vs CIGAR Ref: " << cRefLen << "\n"
                  << "  Qry Block ID: " << qryBlock->getId() << " (Len: " << qrySeq.length() << ") vs CIGAR Qry: " << cQryLen << "\n"
                  << "  -> [ACTION] Aborting link to prevent HEAP CORRUPTION (Core Dump)!\n";
        printCIGAR(cigar);
        exit(1);
    }

    // 1. 扁平化 Segment 以供 TBB 使用
    std::vector<Segment*> refSegsFlat, qrySegsFlat;
    for (auto& seqPair : refBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) refSegsFlat.push_back(&segPairInner.second);
    }
    for (auto& seqPair : qryBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) qrySegsFlat.push_back(&segPairInner.second);
    }

    // 2. 核心運算：取得 Mappings 與計算出新的 Variations
    MergeMappingData mapData = calculateMergedConsensusAndMappings(refBlock, qryBlock, cigar, inverse, refSegsFlat, qrySegsFlat);
    
    // 就地修改 Variants 的相對座標
    updateSegmentVariations(refSegsFlat, qrySegsFlat, mapData, inverse, 
                            refBlock->getConsensus().length(), qryBlock->getConsensus().length());

    // 3. 同步雙方的 Consensus
    std::cout << mapData.mergedConsensus.size() << '\n';
    refBlock->setConsensus(mapData.mergedConsensus, true);
    std::cout << mapData.mergedConsensus.size() << '\n';
    qryBlock->setConsensus(mapData.mergedConsensus);
    
    FamilyID f_id = static_cast<FamilyID>(refBlock->getId());
    refBlock->setFamilyId(f_id);
    qryBlock->setFamilyId(f_id);


    // 4. Update Map (確保兩個 Block 都歸屬於目前的 BlockSet)
    // if (this->getBlock(refBlock->getId()) == nullptr) {
    //     this->addBlock(refBlock); 
    // }
    // if (this->getBlock(qryBlock->getId()) == nullptr) {
    //     this->addBlock(qryBlock); 
    // }
}


// =========================================================
// Helper 1: 計算 Merged Consensus 與所有座標映射
// =========================================================
BlockSet::MergeMappingData BlockSet::calculateMergedConsensusAndMappings(
    std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, 
    const CigarString& cigar, bool inverse, 
    const std::vector<Segment*>& refSegsFlat, const std::vector<Segment*>& qrySegsFlat) 
{
    MergeMappingData data;
    std::string refSeq = refBlock->getConsensus();
    std::string qrySeq = qryBlock->getConsensus();
    int refConsLen = refSeq.length();
    int qryConsLen = qrySeq.length();

    if (inverse) {
        qrySeq = getReverseComplement(qrySeq);
    }

    data.mergedConsensus.reserve(refSeq.length() + qrySeq.length());
    data.refOldToNew.assign(refSeq.length() + 1, 0);
    data.qryOldToNew.assign(qrySeq.length() + 1, 0);

    size_t numRefSegs = refSegsFlat.size();
    size_t totalSegs = numRefSegs + qrySegsFlat.size();

    auto getBaseFromSeg = [](Segment& seg, int pos, char defaultBase, bool needRc, int consLen) -> char {
        int lookupPos = needRc ? (consLen - 1 - pos) : pos; 
        auto& vars = seg.getVariants();
        auto it = std::lower_bound(vars.begin(), vars.end(), lookupPos, 
            [](Variant& v, int p) { return v.getStart() < p; });
        if (it != vars.end() && it->getStart() == lookupPos && it->getType() == VariantType::SNV) {
            char alt = it->getAlt();
            if (needRc) return complement(alt);
            return alt;
        }
        return defaultBase;
    };

    using FreqArray = std::array<int, 256>;
    int rPos = 0, qPos = 0, mPos = 0; 

    for (const auto& op : cigar) {
        int len = op.first;
        char type = op.second;

        if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                char rBase = refSeq[rPos];
                char qBase = qrySeq[qPos];
            
                if (rBase == qBase) {
                    data.mergedConsensus += rBase;
                } else {
                    FreqArray finalCounts = tbb::parallel_reduce(
                        tbb::blocked_range<size_t>(0, totalSegs),
                        FreqArray{}, 
                        [&](const tbb::blocked_range<size_t>& r, FreqArray localCounts) -> FreqArray {
                            for (size_t idx = r.begin(); idx != r.end(); ++idx) {
                                if (idx < numRefSegs) {
                                    char b = getBaseFromSeg(*(refSegsFlat[idx]), rPos, rBase, false, refConsLen);
                                    localCounts[(unsigned char)b]++;
                                } else {
                                    char b = getBaseFromSeg(*(qrySegsFlat[idx - numRefSegs]), qPos, qBase, inverse, qryConsLen);
                                    localCounts[(unsigned char)b]++;
                                }
                            }
                            return localCounts;
                        },
                        [](FreqArray a, const FreqArray& b) -> FreqArray {
                            for(int k=0; k<256; ++k) a[k] += b[k];
                            return a;
                        }
                    );
                
                    char bestBase = rBase; 
                    int maxFreq = -1; 
                    for (int k = 0; k < 256; ++k) {
                        if (finalCounts[k] > maxFreq) {
                            maxFreq = finalCounts[k];
                            bestBase = (char)k;
                        }
                    }
                    data.mergedConsensus += bestBase;
                }

                if (rBase != data.mergedConsensus.back()) data.refConsensusChanges.push_back({rPos, rBase});
                if (qBase != data.mergedConsensus.back()) data.qryConsensusChanges.push_back({qPos, qBase});
                
                data.refOldToNew[rPos] = mPos;
                data.qryOldToNew[qPos] = mPos;
                rPos++; qPos++; mPos++;
            }
        }
        else if (type == 'I') { 
            data.mergedConsensus += qrySeq.substr(qPos, len);
            for(int i = 0; i < len; ++i) data.qryOldToNew[qPos + i] = mPos + i; 
            data.newRefGaps.push_back(Variant::createGap(mPos, mPos + len));
            qPos += len; mPos += len;
        } 
        else if (type == 'D') { 
            data.mergedConsensus += refSeq.substr(rPos, len);
            for(int i = 0; i < len; ++i) data.refOldToNew[rPos + i] = mPos + i; 
            data.newQryGaps.push_back(Variant::createGap(mPos, mPos + len));
            rPos += len; mPos += len;
        }
    }
    data.refOldToNew[rPos] = mPos;
    data.qryOldToNew[qPos] = mPos;

    return data;
}

// =========================================================
// Helper 2: 平行更新 Segment Variations (純更新相對座標)
// =========================================================
void BlockSet::updateSegmentVariations(
    std::vector<Segment*>& refSegsFlat, std::vector<Segment*>& qrySegsFlat,
    const MergeMappingData& mapData, bool inverse, int refConsLen, int qryConsLen) 
{
    struct SegUpdateTask {
        Segment* seg;
        bool isQrySide;
        int originalConsLen;
        const std::vector<int>* oldToNew;
        const std::vector<Variant>* inducedGaps;
        const std::vector<std::pair<int, char>>* consensusChanges;
    };

    std::vector<SegUpdateTask> updateTasks;
    updateTasks.reserve(refSegsFlat.size() + qrySegsFlat.size());
    for (auto* seg : refSegsFlat) updateTasks.push_back({seg, false, refConsLen, &mapData.refOldToNew, &mapData.newRefGaps, &mapData.refConsensusChanges});
    for (auto* seg : qrySegsFlat) updateTasks.push_back({seg, true, qryConsLen, &mapData.qryOldToNew, &mapData.newQryGaps, &mapData.qryConsensusChanges});

    tbb::parallel_for(tbb::blocked_range<size_t>(0, updateTasks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto& task = updateTasks[i];
                Segment& seg = *(task.seg); 
                
                if (task.isQrySide && inverse) seg.reverseVariants(task.originalConsLen);
                
                std::vector<Variant> segmentGaps;
                std::vector<Variant> candidateSnvs;
                
                // 1. 更新原有的 Variations 到新的共識座標系
                for (auto& var : seg.getVariants()) {
                    if (var.getType() == VariantType::GAP) {
                        if (var.getStart() >= task.oldToNew->size() || var.getEnd() >= task.oldToNew->size()) continue;
                        int newStart = (*task.oldToNew)[var.getStart()];
                        int newEnd = (*task.oldToNew)[var.getEnd()];
                        segmentGaps.push_back(Variant::createGap(newStart, newEnd));
                    } else {
                        if (var.getStart() >= task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[var.getStart()];
                        if (var.getAlt() != mapData.mergedConsensus[newPos]) {
                            candidateSnvs.push_back(Variant(newPos, var.getAlt()));
                        }
                    }
                }
                
                // 2. 加入因為比對而新產生的 Induced Gaps
                segmentGaps.insert(segmentGaps.end(), task.inducedGaps->begin(), task.inducedGaps->end());
                std::sort(segmentGaps.begin(), segmentGaps.end(), [](Variant& a, Variant& b) { return a.getStart() < b.getStart(); });

                // 3. 合併重疊的 Gaps
                std::vector<Variant> mergedGaps;
                for (auto& gap : segmentGaps) {
                    if (mergedGaps.empty()) { mergedGaps.push_back(gap); } 
                    else {
                        auto& lastGap = mergedGaps.back();
                        if (lastGap.getEnd() >= gap.getStart()) { 
                            int mStart = lastGap.getStart();
                            int mEnd = std::max(lastGap.getEnd(), gap.getEnd());
                            mergedGaps.pop_back();
                            mergedGaps.push_back(Variant::createGap(mStart, mEnd));
                        } else {
                            mergedGaps.push_back(gap);
                        }
                    }
                }

                // 4. 檢查舊共識被投票蓋掉所產生的新 SNV
                for (auto& change : *(task.consensusChanges)) {
                    int oldPos = change.first; char oldConsBase = change.second;
                    bool hasOldSnv = false;
                    for (auto& v : seg.getVariants()) {
                        if (v.getType() == VariantType::SNV && v.getStart() == oldPos) { hasOldSnv = true; break; }
                    }
                    if (!hasOldSnv) {
                        if (oldPos >= task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[oldPos];
                        if (oldConsBase != mapData.mergedConsensus[newPos]) { 
                            candidateSnvs.push_back(Variant(newPos, oldConsBase));
                        }
                    }
                }
                
                // 5. 過濾掉掉進 Gap 裡面的 SNV
                std::vector<Variant> finalSnvs;
                for (auto& snv : candidateSnvs) {
                    bool inGap = false;
                    for (auto& gap : mergedGaps) {
                        if (snv.getStart() >= gap.getStart() && snv.getStart() < gap.getEnd()) { inGap = true; break; }
                    }
                    if (!inGap) finalSnvs.push_back(snv);
                }
                
                // 6. 整合並排序所有的 Variants
                std::vector<Variant> finalVars;
                finalVars.reserve(mergedGaps.size() + finalSnvs.size());
                finalVars.insert(finalVars.end(), mergedGaps.begin(), mergedGaps.end());
                finalVars.insert(finalVars.end(), finalSnvs.begin(), finalSnvs.end());

                std::sort(finalVars.begin(), finalVars.end(), [](Variant& a, Variant& b) {
                    if (a.getStart() != b.getStart()) return a.getStart() < b.getStart();
                    return a.getType() > b.getType(); 
                });

                std::vector<Variant> cleanedVars;
                for (auto& var : finalVars) {
                    if (cleanedVars.empty()) cleanedVars.push_back(var);
                    else {
                        auto& last = cleanedVars.back();
                        if (last.getStart() == var.getStart() && last.getType() == VariantType::SNV && var.getType() == VariantType::SNV) continue;
                        cleanedVars.push_back(var);
                    }
                }

                // 【修正】：只就地覆寫 Variants，絕對不動 Segment 的 Start/End 絕對座標！
                seg.getVariants() = std::move(cleanedVars);
            }
        }
    );
}