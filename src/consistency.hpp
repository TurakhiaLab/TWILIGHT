#ifndef CONSISTENCY_HPP
#define CONSISTENCY_HPP

#include <atomic>
#include <cstddef>
#include <memory>
#include <vector>
#include <unordered_map>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/spin_rw_mutex.h>

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
    std::unordered_map<ResiduePairKey, float, ResiduePairKeyHash> weights;
    std::unordered_map<ResidueKey, std::vector<std::pair<ResidueKey, float>>, ResidueKeyHash> adjacency;
    std::size_t alignedResiduePairs = 0;

    void addSupport(int seqA, int posA, int seqB, int posB, float weight = 1.0f);
    float getSupport(int seqA, int posA, int seqB, int posB) const;
    const std::vector<std::pair<ResidueKey, float>>& neighbors(int seqId, int pos) const;
    // -------
    void mergeFrom(const DirectPairLibrary& other);
    // -------
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

std::shared_ptr<SubtreeAccurateState> buildSubtreeAccurateState(SequenceDB* database, Option* option, int subtreeIdx);
// -------
float getOrComputePairTcs(SubtreeAccurateState& accurateState, int seqA, int posA, int seqB, int posB);
// -------

} // namespace accurate
} // namespace msa

#endif
