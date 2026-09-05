#include "block_set.hpp"
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// =========================================================
// 靜態輔助函數
// =========================================================
static CigarString invertCigar(const CigarString& cigar) {
    CigarString inverted = cigar;
    for (auto& op : inverted) {
        if (op.second == 'I') op.second = 'D';
        else if (op.second == 'D') op.second = 'I';
    }
    return inverted;
}

#include <unordered_set>

// =========================================================
// 原始的 Merge (建立新 Block，搬移 Sequences，並刪除舊 Block)
// 加入 mode 控制 Copy ID 的分發邏輯
// =========================================================
BlockPtr BlockSet::mergeTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse, int mode, int r_probe_start, int q_probe_start, int r_probe_copy, int q_probe_copy) 
{
    if (!refBlock || !qryBlock || refBlock == qryBlock) return nullptr;

    // 1. 扁平化 Segment 以供 TBB 使用
    std::vector<Segment*> refSegsFlat, qrySegsFlat;
    for (auto& seqPair : refBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) refSegsFlat.push_back(&segPairInner.second);
    }
    for (auto& seqPair : qryBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) qrySegsFlat.push_back(&segPairInner.second);
    }

    // ========================================================
    // 🌟 新增：根據 Mode 設定或平移 Copy ID
    // ========================================================
    int maxRefCopy = -1;
    for (auto seg : refSegsFlat) maxRefCopy = std::max(maxRefCopy, seg->getCopyCount());
    
    int maxQryCopy = -1;
    for (auto seg : qrySegsFlat) maxQryCopy = std::max(maxQryCopy, seg->getCopyCount());

    switch (mode) {
        case 1:
        {
            int target_r_copy = (r_probe_copy >= 0) ? r_probe_copy : 0;
            int target_q_copy = (q_probe_copy >= 0) ? q_probe_copy : 0;

            // 1. Ref 端：所有 segment 保持原本的 copyCount（同源的維持 target_r_copy，額外的維持其原 ID）
            // 2. Qry 端：
            //    - 屬於 target_q_copy 的同源 segment 合併至 target_r_copy (前提是 Ref 端該 sequence 尚未佔用 target_r_copy)
            //    - 其餘額外的 copy 或已存在衝突的 sequence 則平移至 (maxRefCopy + 1) 起跳，避免與 Ref 端衝突
            int qry_shift = std::max(0, maxRefCopy + 1);

            std::unordered_set<std::string> ref_target_seqs;
            for (auto& seqPair : refBlock->getSequences()) {
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    if (segPairInner.second.getCopyCount() == target_r_copy) {
                        ref_target_seqs.insert(seqPair.first);
                        break;
                    }
                }
            }

            for (auto& seqPair : qryBlock->getSequences()) {
                bool ref_has_target = (ref_target_seqs.count(seqPair.first) > 0);
                for (auto& segPairInner : seqPair.second.getSegments()) {
                    Segment& seg = segPairInner.second;
                    if (seg.getCopyCount() == target_q_copy && !ref_has_target) {
                        seg.setCopyCount(target_r_copy);
                    } else {
                        seg.setCopyCount(seg.getCopyCount() + qry_shift);
                    }
                }
            }
            break;
        }
            
        case 2:
        {
            // Ref 端保持原有的 copyCount 分配不變
            // Qry 端每個 segment 則平移 maxRefCopy + 1，確保兩邊 copy ID 互不衝突且各自內部 copy 區隔完整保留
            int shift = std::max(0, maxRefCopy + 1);
            for (auto seg : qrySegsFlat) {
                seg->setCopyCount(seg->getCopyCount() + shift);
            }
            break;
        }
            
        case 3:
        {
            int shift = std::max(0, maxRefCopy + 1);
            for (auto seg : qrySegsFlat) seg->setCopyCount(seg->getCopyCount() + shift);
            break;
        }
            
        case 4:
        {
            int shift = std::max(0, maxQryCopy + 1);
            for (auto seg : refSegsFlat) seg->setCopyCount(seg->getCopyCount() + shift);
            break;
        }
            
        case 5:
            if (refSegsFlat.size() >= qrySegsFlat.size()) {
                int shift = std::max(0, maxRefCopy + 1);
                for (auto seg : qrySegsFlat) seg->setCopyCount(seg->getCopyCount() + shift);
            } else {
                int shift = std::max(0, maxQryCopy + 1);
                for (auto seg : refSegsFlat) seg->setCopyCount(seg->getCopyCount() + shift);
            }
            break;
    }

    // 2. 核心運算：取得 Mappings 與計算出新的 Variations
    MergeMappingData mapData = calculateMergedConsensusAndMappings(refBlock, qryBlock, cigar, inverse, refSegsFlat, qrySegsFlat);
    
    updateSegmentVariations(refSegsFlat, qrySegsFlat, mapData, inverse, refBlock->getConsensus().getLength(), qryBlock->getConsensus().getLength());

    // 3. 建立全新的 Block
    // std::cout << refBlock->getConsensus().getLength() << '\n';
    // std::cout << qryBlock->getConsensus().getLength() << '\n';
    // std::cout << mapData.mergedCons.getConsensusString().size() << '\n';
    auto mergedBlock = this->createBlock(mapData.mergedCons);

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

    // std::cout << "Create " << mergedBlock->getId() << ", deleted " << refBlock->getId() << " and " << qryBlock->getId() << "\n";

    this->deleteBlock(refBlock->getId());
    this->deleteBlock(qryBlock->getId());

    return mergedBlock;
}

// =========================================================
// Helper 1: 計算 Merged Consensus 與所有座標映射
// =========================================================
BlockSet::MergeMappingData BlockSet::calculateMergedConsensusAndMappings(
    std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, 
    const CigarString& cigar, bool inverse, 
    const std::vector<Segment*>& refSegsFlat, const std::vector<Segment*>& qrySegsFlat) 
{
    bool DEBUG_MODE = false;
    MergeMappingData data;
    int refLen = refBlock->getConsensus().getLength();
    int qryLen = qryBlock->getConsensus().getLength();

    bool refIsA = inverse ? true : (refLen >= qryLen);
    
    if (DEBUG_MODE) {
        std::cout << "--- [Before Merge] ---\n";
        std::cout << "Reference Block:\n";
        refBlock->getConsensus().print(std::cout);
        std::cout << "Query Block:\n";
        qryBlock->getConsensus().print(std::cout);
        std::cout << "----------------------\n";
    }

    if (refIsA) {
        data.mergedCons = refBlock->getConsensus();
        data.mergedCons.merge(qryBlock->getConsensus(), cigar, inverse);
    } else {
        data.mergedCons = qryBlock->getConsensus();
        CigarString invCigar = invertCigar(cigar);
        data.mergedCons.merge(refBlock->getConsensus(), invCigar, inverse);
    }
    const Consensus& consA = refIsA ? refBlock->getConsensus() : qryBlock->getConsensus();
    const Consensus& consB = refIsA ? qryBlock->getConsensus() : refBlock->getConsensus();

    std::string seqA = consA.getConsensusString();
    std::string seqB = consB.getConsensusString();

    int lenA = seqA.length();
    int lenB = seqB.length();

    data.refOldToNew.assign(refLen + 1, 0);
    data.qryOldToNew.assign(qryLen + 1, 0);

    CigarString activeCigar = refIsA ? cigar : invertCigar(cigar);

    int aPos = 0, bPos = 0, mPos = 0;

    for (const auto& op : activeCigar) {
        int len = op.first;
        char type = op.second;

        if (type == 'M' || type == '=' || type == 'X') {
            for (int i = 0; i < len; ++i) {
                if (refIsA) {
                    if (aPos < data.refOldToNew.size()) data.refOldToNew[aPos] = mPos;
                    if (bPos < data.qryOldToNew.size()) data.qryOldToNew[bPos] = mPos;

                    if (aPos < lenA && bPos < lenB) {
                        char refBase = seqA[aPos];
                        char qryBase = inverse ? complement(seqB[lenB - 1 - bPos]) : seqB[bPos];
                        if (refBase != qryBase) {
                            data.qryConsensusChanges.push_back({bPos, qryBase});
                        }
                    }
                } else {
                    if (bPos < data.refOldToNew.size()) data.refOldToNew[bPos] = mPos;
                    if (aPos < data.qryOldToNew.size()) data.qryOldToNew[aPos] = mPos;

                    if (aPos < lenA && bPos < lenB) {
                        char qryBase = seqA[aPos];
                        char refBase = inverse ? complement(seqB[lenB - 1 - bPos]) : seqB[bPos];
                        if (refBase != qryBase) {
                            data.refConsensusChanges.push_back({bPos, refBase});
                        }
                    }
                }

                aPos++; bPos++; mPos++;
            }
        }
        else if (type == 'D') { // A 側有字元，B 側是 Gap
            for (int i = 0; i < len; ++i) {
                if (refIsA) {
                    if (aPos + i < data.refOldToNew.size()) data.refOldToNew[aPos + i] = mPos + i;
                } else {
                    if (aPos + i < data.qryOldToNew.size()) data.qryOldToNew[aPos + i] = mPos + i;
                }
            }
            if (refIsA) {
                data.newQryGaps.push_back(Variant::createGap(mPos, mPos + len));
            } else {
                data.newRefGaps.push_back(Variant::createGap(mPos, mPos + len));
            }
            aPos += len; mPos += len;
        }
        else if (type == 'I') { // A 側是 Gap，B 側有字元
            for (int i = 0; i < len; ++i) {
                if (refIsA) {
                    if (bPos + i < data.qryOldToNew.size()) data.qryOldToNew[bPos + i] = mPos + i;
                } else {
                    if (bPos + i < data.refOldToNew.size()) data.refOldToNew[bPos + i] = mPos + i;
                }
            }
            if (refIsA) {
                data.newRefGaps.push_back(Variant::createGap(mPos, mPos + len));
            } else {
                data.newQryGaps.push_back(Variant::createGap(mPos, mPos + len));
            }
            bPos += len; mPos += len;
        }
    }

    if (refIsA) {
        if (aPos < data.refOldToNew.size()) data.refOldToNew[aPos] = mPos;
        if (bPos < data.qryOldToNew.size()) data.qryOldToNew[bPos] = mPos;
    } else {
        if (bPos < data.refOldToNew.size()) data.refOldToNew[bPos] = mPos;
        if (aPos < data.qryOldToNew.size()) data.qryOldToNew[aPos] = mPos;
    }

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

    const std::string& mergedStr = mapData.mergedCons.getConsensusString();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, updateTasks.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto& task = updateTasks[i];
                Segment& seg = *(task.seg); 
                
                if (task.isQrySide && inverse) {
                    seg.reverseVariants(task.originalConsLen);
                }
                
                std::vector<Variant> segmentGaps;
                std::vector<Variant> candidateSnvs;
                
                for (auto& var : seg.getVariants()) {
                    if (var.getType() == VariantType::GAP) {
                        int s = var.getStart();
                        int e = var.getEnd();
                        if (s < 0 || s >= (int)task.oldToNew->size()) continue;
                        int rawStart = (*task.oldToNew)[s];
                        int rawEnd = (e >= (int)task.oldToNew->size()) 
                                     ? (task.oldToNew->empty() ? rawStart : task.oldToNew->back() + 1)
                                     : (e < 0 ? rawStart : (*task.oldToNew)[e]);
                        int newStart = std::min(rawStart, rawEnd);
                        int newEnd = std::max(rawStart, rawEnd);
                        segmentGaps.push_back(Variant::createGap(newStart, newEnd));
                    } else {
                        int s = var.getStart();
                        if (s < 0 || s >= (int)task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[s];
                        if (newPos >= 0 && newPos < (int)mergedStr.length() && var.getAlt() != mergedStr[newPos]) {
                            candidateSnvs.push_back(Variant(newPos, var.getAlt()));
                        }
                    }
                }
                
                segmentGaps.insert(segmentGaps.end(), task.inducedGaps->begin(), task.inducedGaps->end());
                std::sort(segmentGaps.begin(), segmentGaps.end(), [](Variant& a, Variant& b) { return a.getStart() < b.getStart(); });

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

                for (auto& change : *(task.consensusChanges)) {
                    int oldPos = change.first; 
                    char changeBase = change.second;
                    bool hasOldSnv = false;
                    for (auto& v : seg.getVariants()) {
                        if (v.getType() == VariantType::SNV && v.getStart() == oldPos) { hasOldSnv = true; break; }
                    }
                    if (!hasOldSnv) {
                        if (oldPos >= (int)task.oldToNew->size()) continue;
                        int newPos = (*task.oldToNew)[oldPos];
                        if (newPos < (int)mergedStr.length() && changeBase != mergedStr[newPos]) { 
                            candidateSnvs.push_back(Variant(newPos, changeBase));
                        }
                    }
                }
                
                std::vector<Variant> finalSnvs;
                for (auto& snv : candidateSnvs) {
                    bool inGap = false;
                    for (auto& gap : mergedGaps) {
                        if (snv.getStart() >= gap.getStart() && snv.getStart() < gap.getEnd()) { inGap = true; break; }
                    }
                    if (!inGap) finalSnvs.push_back(snv);
                }
                
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

                seg.getVariants() = std::move(cleanedVars);
            }
        }
    );
}