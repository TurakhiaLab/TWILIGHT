#include "msa.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <atomic>
#include <mutex>

#include <tbb/parallel_for.h>

bool msa::accurate::ResidueKey::operator==(const ResidueKey& other) const
{
    return seqId == other.seqId && pos == other.pos;
}

std::size_t msa::accurate::ResidueKeyHash::operator()(const ResidueKey& key) const
{
    std::size_t seed = static_cast<std::size_t>(key.seqId);
    seed = seed * 1315423911u + static_cast<std::size_t>(key.pos);
    return seed;
}

msa::accurate::ResiduePairKey::ResiduePairKey(int firstSeq, int firstPos, int secondSeq, int secondPos)
{
    if (firstSeq < secondSeq || (firstSeq == secondSeq && firstPos <= secondPos)) {
        this->seqA = firstSeq;
        this->posA = firstPos;
        this->seqB = secondSeq;
        this->posB = secondPos;
    }
    else {
        this->seqA = secondSeq;
        this->posA = secondPos;
        this->seqB = firstSeq;
        this->posB = firstPos;
    }
}

bool msa::accurate::ResiduePairKey::operator==(const ResiduePairKey& other) const
{
    return seqA == other.seqA && posA == other.posA && seqB == other.seqB && posB == other.posB;
}

std::size_t msa::accurate::ResiduePairKeyHash::operator()(const ResiduePairKey& key) const
{
    std::size_t seed = static_cast<std::size_t>(key.seqA);
    seed = seed * 1315423911u + static_cast<std::size_t>(key.posA);
    seed = seed * 1315423911u + static_cast<std::size_t>(key.seqB);
    seed = seed * 1315423911u + static_cast<std::size_t>(key.posB);
    return seed;
}


std::string getCurrentSequence(const msa::SequenceDB::SequenceInfo* sequence)
{
    return std::string(sequence->alnStorage[sequence->storage], sequence->len);
}

float computePairNew(msa::accurate::DirectPairLibrary& directLib, int seqA, int posA, int seqB, int posB)
{
    float directSupport = 0.0f;
    if (directLib.is_aligned(seqA, seqB, posA, posB)) {
        directSupport = directLib.getWeight(seqA, seqB);
    }

    float numerator = directSupport;
    float supportA = directSupport;
    float supportB = directSupport;

    bool hasNeighborA = false;
    bool hasNeighborB = false;

    int totalSeq = directLib.getSequence();

    for (int seqC = 0; seqC < totalSeq; ++seqC) {
        if (seqC == seqA || seqC == seqB) continue;

        int alignedPosC_A = directLib.getAlignedPos(seqA, seqC, posA);

        int alignedPosC_B = directLib.getAlignedPos(seqB, seqC, posB);

        if (alignedPosC_A != -1) {
            hasNeighborA = true;
            supportA += directLib.getWeight(seqA, seqC);
        }
        if (alignedPosC_B != -1) {
            hasNeighborB = true;
            supportB += directLib.getWeight(seqB, seqC);
        }

        if (alignedPosC_A != -1 && alignedPosC_B != -1 && alignedPosC_A == alignedPosC_B) {
            float wA = directLib.getWeight(seqA, seqC);
            float wB = directLib.getWeight(seqB, seqC);
            if (wA > 0.0f && wB > 0.0f) {
                numerator += std::min(wA, wB);
            }
        }
    }

    if (!hasNeighborA && !hasNeighborB) return directSupport;

    const float denominator = supportA + supportB;
    if (denominator <= 0.0f) return 0.0f;

    return std::clamp((2.0f * numerator) / denominator, 0.0f, 1.0f);
}

msa::accurate::DirectPairLibrary::DirectPairLibrary(int sequenceCount, int length) {
    this->length = length;
    this->totalSequence = sequenceCount;
    this->totalPairs = (sequenceCount + 1) * sequenceCount / 2 - sequenceCount;
    this->weights.resize(totalPairs, -1.0f);
    this->forward.resize(totalPairs * length, -1);
    this->backward.resize(totalPairs * length, -1);
    this->pairWeight.resize(totalPairs * length, 0.0f);
}

// --- Getter ---
int msa::accurate::DirectPairLibrary::getSequence() {return totalSequence;}
int msa::accurate::DirectPairLibrary::getPairs() {return totalPairs;}
inline uint64_t msa::accurate::DirectPairLibrary::getOffset(int i, int j) { return static_cast<uint64_t>(length) * static_cast<uint64_t>(idx(i,j));}
// --- Modifier ---
void msa::accurate::DirectPairLibrary::addLocalAlignmentResult(int refID, int qryID, LocalAlignmentResult& result) {
    uint64_t offset = getOffset(refID, qryID);
    int id = idx(refID, qryID);
    this->weights[id] = result.identity;
    if (refID < qryID) { // Ref-based
        for (auto& alignedPair: result.alignedPairs) {
            int base = alignedPair.refIndex;
            int align = alignedPair.qryIndex;
            this->forward[offset+base] = align;
            this->backward[offset+align] = base;
        }
    }
    else {
        for (auto& alignedPair: result.alignedPairs) {
            int base = alignedPair.qryIndex;
            int align = alignedPair.refIndex;
            this->forward[offset+base] = align;
            this->backward[offset+align] = base;
        }
    }
}
// --- lookup ---
float msa::accurate::DirectPairLibrary::getWeight(int refIdx, int qryIdx) {
    return this->weights[idx(refIdx, qryIdx)];
};
bool msa::accurate::DirectPairLibrary::is_aligned(int refIdx, int qryIdx, int refPos, int qryPos) {
    int offset = getOffset(refIdx, qryIdx);
    if (refIdx < qryIdx) return (this->forward[offset+refPos] == qryPos);
    else                 return (this->backward[offset+refPos] == qryPos);
}
int msa::accurate::DirectPairLibrary::getAlignedPos(int baseIdx, int alnIdx, int basePos) {
    int offset = getOffset(baseIdx, alnIdx);
    if (baseIdx < alnIdx) return this->forward[offset+basePos];
    else                  return this->backward[offset+basePos];
}
float msa::accurate::DirectPairLibrary::getPairWeight(int refIdx, int qryIdx, int refPos) {
    int offset = getOffset(refIdx, qryIdx);
    if (refIdx < qryIdx) return this->pairWeight[offset+refPos];
    else {
        int qryPos = this->backward[offset+refPos];
        if (qryPos == -1) return 0.0f;
        return this->pairWeight[offset+qryPos];
    }
}

// --- Mapping function ---
inline int msa::accurate::DirectPairLibrary::idx(int i, int j) {
    if (i == j) {
        std::cerr << "ERROR: Lookup pairwise alignment of same sequence\n.";
        return -1;
    }
    if (i > j) std::swap(i, j);
    return i * totalSequence - i * (i + 1) / 2 + (j - i - 1);
}

std::shared_ptr<msa::accurate::SubtreeAccurateState> msa::accurate::buildSubtreeAccurateState(SequenceDB* database, Option* option, int subtreeIdx, Params& params)
{
    auto time0 = std::chrono::high_resolution_clock::now();
    auto accurateState = std::make_shared<SubtreeAccurateState>();
    accurateState->subtreeIdx = subtreeIdx; 
    SmithWatermanAligner aligner;

    const auto& sequences = database->sequences;
    // -------
    std::vector<std::size_t> activeSeqIdx;
    activeSeqIdx.reserve(sequences.size());
    for (std::size_t seqIdx = 0; seqIdx < sequences.size(); ++seqIdx) {
        if (!sequences[seqIdx]->lowQuality || option->noFilter) activeSeqIdx.push_back(seqIdx);
    }

    // activeSeqIdx: vectorIdx -> seqIdx
    // currentSequences: seqIdx -> raw sequence

    int maxLength = -1;
    int maxSeqID = -1;

    std::vector<std::string> currentSequences(sequences.size());
    for (std::size_t idx = 0; idx < activeSeqIdx.size(); ++idx) {
        const std::size_t seqIdx = activeSeqIdx[idx];
        auto seq = getCurrentSequence(sequences[seqIdx]);
        currentSequences[seqIdx] = seq;
        if (static_cast<int>(seq.size()) > maxLength) maxLength = static_cast<int>(seq.size());
        if (static_cast<int>(seqIdx) > maxSeqID) maxSeqID = static_cast<int>(seqIdx);
    }
    
    // Construct Lookup Table
    DirectPairLibrary ConsistencyLibrary(maxSeqID+1, maxLength);

    std::vector<std::pair<std::size_t, std::size_t>> pairJobs;
    pairJobs.reserve((activeSeqIdx.size() + 1) * activeSeqIdx.size() / 2 - activeSeqIdx.size());
    
    for (std::size_t refPos = 0; refPos < activeSeqIdx.size(); ++refPos) {
        for (std::size_t qryPos = refPos + 1; qryPos < activeSeqIdx.size(); ++qryPos) {
            pairJobs.push_back({activeSeqIdx[refPos], activeSeqIdx[qryPos]});
        }
    }

    std::atomic<size_t> progress{0};
    std::mutex cout_mutex;
    size_t total = pairJobs.size();

    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, total), [&](const tbb::blocked_range<std::size_t>& range) {
        for (std::size_t pairIdx = range.begin(); pairIdx < range.end(); ++pairIdx) {
            const auto [refIdx, qryIdx] = pairJobs[pairIdx];
            const auto* refSeq = sequences[refIdx];
            const auto* qrySeq = sequences[qryIdx];
            auto alignment = aligner.align_affine(
                currentSequences[refIdx],
                currentSequences[qryIdx],
                option->type,
                params
            );
            ConsistencyLibrary.addLocalAlignmentResult(
                refSeq->id, qrySeq->id, alignment
            );

            // print progress
            size_t current = ++progress;
            if ((current % 10 == 0 || current == total) && pairIdx == range.begin()) {
                double percent = 100.0 * current / total;
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "1. All-to-all Pairwise Alignment: ["
                          << current << "/" << total
                          << " pairs aligned] ("
                          << std::fixed << std::setprecision(1)
                          << percent << "%)\r"
                          << std::flush;
            }
        }
    });

    std::cout << std::endl;

    accurateState.get()->directLib = ConsistencyLibrary;
    
    auto time1 = std::chrono::high_resolution_clock::now();
    accurateState.get()->directLib.computeResidueSupport();
    auto time2 = std::chrono::high_resolution_clock::now();
    accurateState.get()->directLib.computePairWeights();
    auto time3 = std::chrono::high_resolution_clock::now();

    if (option->printDetail) {
        auto ms = [](auto a, auto b) {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1000000;
        };
        std::cerr << "Built accurate-mode direct library for subtree " << subtreeIdx
                  << " with " << ConsistencyLibrary.getPairs() 
                  << " pairs in " << ms(time0, time3) << " ms.\n"
                  << "Time breakdown: \n"
                  << "  1. Build pairwise library: " << ms(time0, time1) << " ms\n"
                  << "  2. Compute residue support: " << ms(time1, time2) << " ms\n"
                  << "  3. Compute pair weights: " << ms(time2, time3) << " ms\n";
    }

    return accurateState;
}

// -------

float msa::accurate::getOrComputePairNew(SubtreeAccurateState& accurateState, int seqA, int posA, int seqB, int posB)
{
    return computePairNew(accurateState.directLib, seqA, posA, seqB, posB);
}

// -------




// -------
void msa::accurate::DirectPairLibrary::computePairWeights() {
    pairWeight.assign(totalPairs * length, 0.0f);

    std::atomic<int> progress{0};
    std::mutex cout_mutex;
    int total = totalSequence;

    tbb::parallel_for(
        tbb::blocked_range<int>(0, totalSequence),
        [&](const tbb::blocked_range<int>& range) {

        for (int seqA = range.begin(); seqA < range.end(); ++seqA) {

            for (int seqB = seqA + 1; seqB < totalSequence; ++seqB) {
                float w_AB = getWeight(seqA, seqB);
                uint64_t offset = getOffset(seqA, seqB);

                for (int posA = 0; posA < length; ++posA) {
                    int posB = forward[offset + posA];
                    if (posB == -1) continue;

                    float numerator = w_AB > 0.0f ? w_AB : 0.0f;

                    for (int seqC = 0; seqC < totalSequence; ++seqC) {
                        if (seqC == seqA || seqC == seqB) continue;

                        int posC_A = getAlignedPos(seqA, seqC, posA);
                        if (posC_A == -1) continue;

                        int posC_B = getAlignedPos(seqB, seqC, posB);
                        if (posC_B == -1) continue;

                        if (posC_A == posC_B) {
                            float w_AC = getWeight(seqA, seqC);
                            float w_BC = getWeight(seqB, seqC);
                            if (w_AC > 0.0f && w_BC > 0.0f) {
                                numerator += std::min(w_AC, w_BC);
                            }
                        }
                    }

                    pairWeight[offset + posA] = numerator;
                }
            }

            int current = ++progress;

            if ((current % 5 == 0 || current == total) && seqA == range.begin()) {
                double percent = 100.0 * current / total;

                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "3. Compute Pair Weights: "
                          << current << "/" << total
                          << "] ("
                          << std::fixed << std::setprecision(1)
                          << percent << "%)\r"
                          << std::flush;
            }
        }
    });

    std::cout << std::endl;
}

// -------
void msa::accurate::DirectPairLibrary::computeResidueSupport() {
    std::cerr << "2. Compute Residue Support: ";
    residueSupport.assign(totalSequence * length, 0.0f);
    for (int i = 0; i < totalSequence; ++i) {
        for (int j = i + 1; j < totalSequence; ++j) {
            float w = getWeight(i, j);
            if (w <= 0.0f) continue;
            uint64_t offset = getOffset(i, j);
            for (int posI = 0; posI < length; ++posI) {
                int posJ = forward[offset + posI];
                if (posJ != -1) {
                    residueSupport[i * length + posI] += w;
                    residueSupport[j * length + posJ] += w;
                }
            }
        }
    }
    std::cerr << "Done.\n";
}


void msa::accurate::removeColumns(ColumnProvenance& provenance, const IntPairVec& removedColumns)
{
    if (removedColumns.empty()) return;
    ColumnProvenance filtered;
    filtered.reserve(provenance.size());
    int removeIdx = 0;
    int nextRemoveStart = removedColumns[removeIdx].first;
    int nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
    for (int col = 0; col < static_cast<int>(provenance.size()); ++col) {
        while (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveEnd) {
            ++removeIdx;
            if (removeIdx < static_cast<int>(removedColumns.size())) {
                nextRemoveStart = removedColumns[removeIdx].first;
                nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
            }
        }
        const bool removed = (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveStart && col < nextRemoveEnd);
        if (!removed) filtered.push_back(std::move(provenance[col]));
    }
    provenance = std::move(filtered);
}

std::vector<std::vector<float>> msa::accurate::buildConsistencyTable(
    const ColumnProvenance& refProvenance,
    const ColumnProvenance& qryProvenance,
    SubtreeAccurateState& accurateState)
{
    int refLen = refProvenance.size();
    int qryLen = qryProvenance.size();
    std::vector<std::vector<float>> consistencyTable(refLen, std::vector<float>(qryLen, 0.0f));

    int totalSeq = accurateState.directLib.getSequence();
    int maxLen = accurateState.directLib.length;

    struct QryResInfo {
        int qryCol;
        int seqB;
        int posB;
        float weight;
        float S_B;
    };

    std::vector<QryResInfo> qryResList;
    std::vector<int> qryMap(totalSeq * maxLen, -1);
    std::vector<float> qryColWeightSum(qryLen, 0.0f);
    std::vector<float> refColWeightSum(refLen, 0.0f);

    for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
        for (const auto& qryRes : qryProvenance[qryCol]) {
            int seqB = qryRes.seqId;
            int posB = qryRes.residueIndex;
            int id = qryResList.size();
            float S_B = accurateState.directLib.residueSupport[seqB * maxLen + posB];
            qryResList.push_back({qryCol, seqB, posB, qryRes.weight, S_B});
            qryMap[seqB * maxLen + posB] = id;
            qryColWeightSum[qryCol] += qryRes.weight;
        }
    }

    for (int refCol = 0; refCol < refLen; ++refCol) {
        for (const auto& refRes : refProvenance[refCol]) {
            refColWeightSum[refCol] += refRes.weight;
        }
    }

    int qryResCount = qryResList.size();

    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, refLen), [&](const tbb::blocked_range<std::size_t>& range) {
        std::vector<float> num_B(qryResCount, 0.0f);
        std::vector<int> active_B;
        active_B.reserve(totalSeq);

        for (std::size_t refCol = range.begin(); refCol < range.end(); ++refCol) {
            if (refProvenance[refCol].empty()) continue;

            std::vector<float> col_scores(qryLen, 0.0f);
            for (const auto& refRes : refProvenance[refCol]) {
                int seqA = refRes.seqId;
                int posA = refRes.residueIndex;
                float wA = refRes.weight;
                float S_A = accurateState.directLib.residueSupport[seqA * maxLen + posA];
                if (S_A <= 0.0f) continue;

                for (int seqB = 0; seqB < totalSeq; ++seqB) {
                    if (seqB == seqA) continue;
                    int posB = accurateState.directLib.getAlignedPos(seqA, seqB, posA);
                    if (posB != -1) {
                        int qID = qryMap[seqB * maxLen + posB];
                        if (qID != -1) {
                            float extendedWeight = accurateState.directLib.getPairWeight(seqA, seqB, posA);
                            if (extendedWeight > 0.0f) {
                                if (num_B[qID] == 0.0f) active_B.push_back(qID);
                                num_B[qID] += extendedWeight;
                            }
                        }
                    }
                }
                
                for (int qID : active_B) {
                    float num = num_B[qID];
                    num_B[qID] = 0.0f;
                    const auto& qInfo = qryResList[qID];
                    float S_B = qInfo.S_B;
                    float tcs = (S_A + S_B > 0.0f) ? std::clamp((2.0f * num) / (S_A + S_B), 0.0f, 1.0f) : 0.0f;
                    col_scores[qInfo.qryCol] += wA * qInfo.weight * tcs;
                }
                active_B.clear();
            }
            
            float denom_ref = refColWeightSum[refCol];
            if (denom_ref <= 0.0f) continue;
            for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
                float denom = denom_ref * qryColWeightSum[qryCol];
                if (denom > 0.0f) {
                    consistencyTable[refCol][qryCol] = col_scores[qryCol] / denom;
                }
            }
        }
    });

    return consistencyTable;
}
