#include "consistency.hpp"

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

void msa::accurate::DirectPairLibrary::addSupport(int seqA, int posA, int seqB, int posB, float weight)
{
    ResiduePairKey key(seqA, posA, seqB, posB);
    auto it = weights.find(key);
    if (it == weights.end()) {
        weights.emplace(key, weight);
        adjacency[{seqA, posA}].push_back({{seqB, posB}, weight});
        adjacency[{seqB, posB}].push_back({{seqA, posA}, weight});
        ++alignedResiduePairs;
    }
    else {
        it->second = std::max(it->second, weight);
        for (auto& neighbor : adjacency[{seqA, posA}]) {
            if (neighbor.first.seqId == seqB && neighbor.first.pos == posB) {
                neighbor.second = std::max(neighbor.second, weight);
                break;
            }
        }
        for (auto& neighbor : adjacency[{seqB, posB}]) {
            if (neighbor.first.seqId == seqA && neighbor.first.pos == posA) {
                neighbor.second = std::max(neighbor.second, weight);
                break;
            }
        }
    }
}

float msa::accurate::DirectPairLibrary::getSupport(int seqA, int posA, int seqB, int posB) const
{
    auto it = weights.find(ResiduePairKey(seqA, posA, seqB, posB));
    return (it == weights.end()) ? 0.0f : it->second;
}

const std::vector<std::pair<msa::accurate::ResidueKey, float>>& msa::accurate::DirectPairLibrary::neighbors(int seqId, int pos) const
{
    static const std::vector<std::pair<ResidueKey, float>> empty;
    auto it = adjacency.find({seqId, pos});
    return (it == adjacency.end()) ? empty : it->second;
}

// -------
void msa::accurate::DirectPairLibrary::mergeFrom(const DirectPairLibrary& other)
{
    for (const auto& [key, weight] : other.weights) {
        addSupport(key.seqA, key.posA, key.seqB, key.posB, weight);
    }
}
// -------

namespace {

std::string getCurrentSequence(const msa::SequenceDB::SequenceInfo* sequence)
{
    return std::string(sequence->alnStorage[sequence->storage], sequence->len);
}

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

} // namespace

std::shared_ptr<msa::accurate::SubtreeAccurateState> msa::accurate::buildSubtreeAccurateState(SequenceDB* database, Option* option, int subtreeIdx)
{
    auto buildStart = std::chrono::high_resolution_clock::now();
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

    std::vector<std::string> currentSequences(sequences.size());
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, activeSeqIdx.size()),
        [&](const tbb::blocked_range<std::size_t>& range) {
            for (std::size_t idx = range.begin(); idx < range.end(); ++idx) {
                const std::size_t seqIdx = activeSeqIdx[idx];
                currentSequences[seqIdx] = getCurrentSequence(sequences[seqIdx]);
            }
        }
    );

    std::vector<std::pair<std::size_t, std::size_t>> pairJobs;
    pairJobs.reserve((activeSeqIdx.size() * (activeSeqIdx.size() - 1)) / 2);
    for (std::size_t refPos = 0; refPos < activeSeqIdx.size(); ++refPos) {
        for (std::size_t qryPos = refPos + 1; qryPos < activeSeqIdx.size(); ++qryPos) {
            pairJobs.push_back({activeSeqIdx[refPos], activeSeqIdx[qryPos]});
        }
    }

    tbb::enumerable_thread_specific<DirectPairLibrary> partialLibraries;
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, pairJobs.size()),
        [&](const tbb::blocked_range<std::size_t>& range) {
            auto& localLibrary = partialLibraries.local();
            for (std::size_t pairIdx = range.begin(); pairIdx < range.end(); ++pairIdx) {
                const auto [refIdx, qryIdx] = pairJobs[pairIdx];
                const auto* refSeq = sequences[refIdx];
                const auto* qrySeq = sequences[qryIdx];
                const auto alignment = aligner->align(currentSequences[refIdx], currentSequences[qryIdx], option->type);
                for (const auto& residuePair : alignment.alignedPairs) {
                    localLibrary.addSupport(
                        refSeq->id,
                        residuePair.refIndex,
                        qrySeq->id,
                        residuePair.qryIndex,
                        alignment.identity
                    );
                }
            }
        }
    );

    for (const auto& partialLibrary : partialLibraries) {
        accurateState->directLib.mergeFrom(partialLibrary);
    }
    // -------

    if (option->printDetail) {
        auto buildEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds buildTime = buildEnd - buildStart;
        std::cerr << "Built accurate-mode direct library for subtree " << subtreeIdx
                  << " with " << accurateState->directLib.alignedResiduePairs
                  << " residue pairs in " << buildTime.count() / 1000000 << " ms.\n";
    }

    return accurateState;
}

// -------
float msa::accurate::getOrComputePairTcs(SubtreeAccurateState& accurateState, int seqA, int posA, int seqB, int posB)
{
    ResiduePairKey key(seqA, posA, seqB, posB);
    {
        tbb::spin_rw_mutex::scoped_lock readLock(accurateState.pairTcsCache.mutex, false);
        auto it = accurateState.pairTcsCache.values.find(key);
        if (it != accurateState.pairTcsCache.values.end()) {
            ++accurateState.pairTcsCache.cacheHits;
            return it->second;
        }
    }

    const float pairTcs = computePairTcs(accurateState.directLib, seqA, posA, seqB, posB);
    {
        tbb::spin_rw_mutex::scoped_lock writeLock(accurateState.pairTcsCache.mutex);
        auto [it, inserted] = accurateState.pairTcsCache.values.emplace(key, pairTcs);
        if (inserted) {
            ++accurateState.pairTcsCache.cacheMisses;
            return pairTcs;
        }
        ++accurateState.pairTcsCache.cacheHits;
        return it->second;
    }
}
// -------
