#ifndef CONSISTENCY_HPP
#define CONSISTENCY_HPP

#include "local-align.hpp"

#include <atomic>
#include <cstddef>
#include <memory>
#include <vector>
#include <unordered_map>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_rw_mutex.h>
#include <iostream>

namespace msa {

struct Option;
struct SequenceDB;

namespace accurate {

struct ResidueKey
{
    int seqId;
    int pos;

    bool operator==(const ResidueKey& other) const;
};

struct ResidueKeyHash
{
    std::size_t operator()(const ResidueKey& key) const;
};

struct ResiduePairKey
{
    int seqA;
    int posA;
    int seqB;
    int posB;

    ResiduePairKey(int firstSeq, int firstPos, int secondSeq, int secondPos);
    bool operator==(const ResiduePairKey& other) const;
};

struct ResiduePairKeyHash
{
    std::size_t operator()(const ResiduePairKey& key) const;
};

struct DirectPairLibrary
{
    // --- Constructor ---
    DirectPairLibrary() {};
    DirectPairLibrary(int sequenceCount, int length) {
        this->length = length;
        this->totalSequence = sequenceCount;
        this->totalPairs = (sequenceCount + 1) * sequenceCount / 2 - sequenceCount;
        this->weights.resize(totalPairs, -1.0f);
        this->forward.resize(totalPairs * length, -1);
        this->backward.resize(totalPairs * length, -1);
        this->pairWeight.resize(totalPairs * length, 0.0f);
    }

    // --- Getter ---
    int getSequence() {return totalSequence;}
    int getPairs() {return totalPairs;}
    inline uint64_t getOffset(int i, int j) {
        return static_cast<uint64_t>(length) * static_cast<uint64_t>(idx(i,j));
    }

    // --- Modifier ---
    void addLocalAlignmentResult(int refID, int qryID, LocalAlignmentResult& result);

    // --- lookup ---
    float getWeight(int refIdx, int qryIdx) {
        return this->weights[idx(refIdx, qryIdx)];
    };
    bool is_aligned(int refIdx, int qryIdx, int refPos, int qryPos) {
        int offset = getOffset(refIdx, qryIdx);
        if (refIdx < qryIdx) return (this->forward[offset+refPos] == qryPos);
        else                 return (this->backward[offset+refPos] == qryPos);
    }
    int getAlignedPos(int baseIdx, int alnIdx, int basePos) {
        int offset = getOffset(baseIdx, alnIdx);
        if (baseIdx < alnIdx) return this->forward[offset+basePos];
        else                  return this->backward[offset+basePos];
    }
    float getPairWeight(int refIdx, int qryIdx, int refPos) {
        int offset = getOffset(refIdx, qryIdx);
        if (refIdx < qryIdx) return this->pairWeight[offset+refPos];
        else {
            int qryPos = this->backward[offset+refPos];
            if (qryPos == -1) return 0.0f;
            return this->pairWeight[offset+qryPos];
        }
    }
    
    // --- Mapping function ---
    inline int idx(int i, int j) {
        if (i == j) {
            std::cerr << "ERROR: Lookup pairwise alignment of same sequence\n.";
            return -1;
        }
        if (i > j) std::swap(i, j);
        return i * totalSequence - i * (i + 1) / 2 + (j - i - 1);
    }
    std::vector<float> weights;
    std::vector<int32_t> forward;
    std::vector<int32_t> backward;
    std::vector<float> pairWeight;
    int totalSequence;
    int totalPairs;
    int length;
    
    // --- Pre-compute supportive values
    std::vector<float> residueSupport;
    void computeResidueSupport();
    void computePairWeights();

};

struct PairTcsCache
{
    tbb::spin_rw_mutex mutex;
    std::unordered_map<ResiduePairKey, float, ResiduePairKeyHash> values;
    std::atomic<std::size_t> cacheHits{0};
    std::atomic<std::size_t> cacheMisses{0};
};

struct SubtreeAccurateState
{
    int subtreeIdx = -1;
    DirectPairLibrary directLib;
    // -------
    PairTcsCache pairTcsCache;
    // -------
};

std::shared_ptr<SubtreeAccurateState> buildSubtreeAccurateState(SequenceDB* database, Option* option, int subtreeIdx, Params& param);
// -------
// float getOrComputePairTcs(SubtreeAccurateState& accurateState, int seqA, int posA, int seqB, int posB);
float getOrComputePairNew(SubtreeAccurateState& accurateState, int seqA, int posA, int seqB, int posB);
// -------

} // namespace accurate
} // namespace msa

#endif
