#include "block_set.hpp"
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// =========================================================
// Static Helper Functions
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
// Merge Two Blocks (Create new Block, migrate Sequences, and delete old Blocks)
// Include mode parameter to control Copy ID assignment logic
// =========================================================
BlockPtr BlockSet::mergeTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse, int mode, int r_probe_start, int q_probe_start, int r_probe_copy, int q_probe_copy) 
{
    if (!refBlock || !qryBlock || refBlock == qryBlock) return nullptr;

    // 1. Flatten Segments for TBB parallel processing
    std::vector<Segment*> refSegsFlat, qrySegsFlat;
    for (auto& seqPair : refBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) refSegsFlat.push_back(&segPairInner.second);
    }
    for (auto& seqPair : qryBlock->getSequences()) {
        for (auto& segPairInner : seqPair.second.getSegments()) qrySegsFlat.push_back(&segPairInner.second);
    }

    // ========================================================
    // Set or shift Copy ID based on mode
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

            // 1. Ref side: Keep original copyCount for all segments (homologous stays target_r_copy, extras keep their original ID)
            // 2. Qry side:
            //    - Homologous segments belonging to target_q_copy are merged into target_r_copy (provided that the sequence in Ref does not already occupy target_r_copy)
            //    - Other extra copies or sequences with conflicting copies are shifted starting from (maxRefCopy + 1) to avoid collision with Ref
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
            // Keep original copyCount allocation on Ref side
            // Shift each segment on Qry side by (maxRefCopy + 1) to ensure no copy ID collision while preserving internal copy distinctions
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

    // 2. Core computation: Obtain mappings and calculate updated variations
    MergeMappingData mapData = calculateMergedConsensusAndMappings(refBlock, qryBlock, cigar, inverse, refSegsFlat, qrySegsFlat);
    
    updateSegmentVariations(refSegsFlat, qrySegsFlat, mapData, inverse, refBlock->getConsensus().getLength(), qryBlock->getConsensus().getLength());

    // 3. Create a new Block
    // std::cout << refBlock->getConsensus().getLength() << '\n';
    // std::cout << qryBlock->getConsensus().getLength() << '\n';
    // std::cout << mapData.mergedCons.getConsensusString().size() << '\n';
    auto mergedBlock = this->createBlock(mapData.mergedCons);

    // 4. Move and attach sequences to the new Block
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
// Helper 1: Calculate Merged Consensus and all coordinate mappings
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
        else if (type == 'D') { // Deletion: char on A side, gap on B side
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
        else if (type == 'I') { // Insertion: gap on A side, char on B side
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
// Helper 2: Parallel update of Segment Variations (relative coordinates only)
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