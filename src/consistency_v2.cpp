#include "consistency_v2.hpp"

#include "local-align.hpp"
#include "msa.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include <tbb/enumerable_thread_specific.h>
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

namespace {

std::string getCurrentSequence(const msa::SequenceDB::SequenceInfo* sequence)
{
    return std::string(sequence->alnStorage[sequence->storage], sequence->len);
}

/*
float computePairTcs(const msa::accurate::DirectPairLibrary& directLib, int seqA, int posA, int seqB, int posB)
{
    const float directSupport = directLib.getSupport(seqA, posA, seqB, posB);
    const auto& neighborsA = directLib.neighbors(seqA, posA);
    const auto& neighborsB = directLib.neighbors(seqB, posB);

    if (neighborsA.empty() && neighborsB.empty()) return directSupport;

    float numerator = directSupport;
    float supportA = directSupport;
    float supportB = directSupport;

    for (const auto& neighborA : neighborsA) supportA += neighborA.second;
    for (const auto& neighborB : neighborsB) supportB += neighborB.second;

    for (const auto& neighborA : neighborsA) {
        const float supportThroughB = directLib.getSupport(
            neighborA.first.seqId,
            neighborA.first.pos,
            seqB,
            posB
        );
        if (supportThroughB > 0.0f) {
            numerator += std::min(neighborA.second, supportThroughB);
        }
    }

    const float denominator = supportA + supportB;
    if (denominator <= 0.0f) return 0.0f;

    return std::clamp((2.0f * numerator) / denominator, 0.0f, 1.0f);
}
*/

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

} // namespace

std::shared_ptr<msa::accurate::SubtreeAccurateState> msa::accurate::buildSubtreeAccurateState(SequenceDB* database, Option* option, int subtreeIdx)
{
    auto time0 = std::chrono::high_resolution_clock::now();
    auto accurateState = std::make_shared<SubtreeAccurateState>();
    accurateState->subtreeIdx = subtreeIdx;
    auto aligner = makeDefaultLocalAligner();

    const auto& sequences = database->sequences;
    // -------
    std::vector<std::size_t> activeSeqIdx;
    activeSeqIdx.reserve(sequences.size());
    for (std::size_t seqIdx = 0; seqIdx < sequences.size(); ++seqIdx) {
        if (!sequences[seqIdx]->lowQuality) activeSeqIdx.push_back(seqIdx);
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

    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, pairJobs.size()), [&](const tbb::blocked_range<std::size_t>& range) {
        for (std::size_t pairIdx = range.begin(); pairIdx < range.end(); ++pairIdx) {
            const auto [refIdx, qryIdx] = pairJobs[pairIdx];
            const auto* refSeq = sequences[refIdx];
            const auto* qrySeq = sequences[qryIdx];
            auto alignment = aligner->align(currentSequences[refIdx], currentSequences[qryIdx], option->type);
            ConsistencyLibrary.addLocalAlignmentResult(refSeq->id, qrySeq->id, alignment);
        }
    });

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

// -------
void msa::accurate::DirectPairLibrary::computePairWeights() {
    pairWeight.assign(totalPairs * length, 0.0f);
    
    tbb::parallel_for( tbb::blocked_range<int>(0, totalSequence), [&](const tbb::blocked_range<int>& range) {
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
        }
    });
}

// -------
void msa::accurate::DirectPairLibrary::computeResidueSupport() {
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
}
