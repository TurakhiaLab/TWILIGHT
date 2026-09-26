#ifndef LOCAL_ALIGN_HPP
#define LOCAL_ALIGN_HPP

#include "msa.hpp"

#include <memory>
#include <string>
#include <vector>

namespace msa {

    struct Params;

namespace accurate {

struct AlignedResiduePair
{
    int refIndex;
    int qryIndex;
};

struct AlignmentResult
{
    int score = 0;
    float identity = 0.0f;
    std::vector<AlignedResiduePair> alignedPairs;
};

struct Aligner
{
    // Reusable DP scratchpad buffers to eliminate per-pair heap allocations
    std::vector<int> score_0;
    std::vector<int> score_1;
    std::vector<int> score_2;
    std::vector<int> E_buf;
    std::vector<uint8_t> tb_buf;
    std::vector<int> qry_idx_buf;
    std::vector<AlignedResiduePair> aligned_pairs_buf;
    std::vector<std::pair<int, int>> pair_idx_buf;

    AlignmentResult align(const std::string& reference, const std::string& query, char type, Params& params);
    AlignmentResult align_affine(const std::string& reference, const std::string& query, char type, Params& params);
    AlignmentResult align_affine_local(const std::string& reference, const std::string& query, char type, Params& params);
    AlignmentResult align_linear_local (const std::string& reference, const std::string& query, char type, msa::Params& params);
    LocalHomTable::PairResult align_affine_local_segments(const std::string& reference, const std::string& query, char type, msa::Params& params);
    LocalHomTable::PairResult align_affine_local_segments_banded(const std::string& reference, const std::string& query, char type, msa::Params& params, int bandWidth);
};

} // namespace accurate
} // namespace msa




#endif
