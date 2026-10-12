#ifndef ALIGN_METHOD_HPP
#define ALIGN_METHOD_HPP

#include "msa.hpp"

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <utility>

namespace msa {

namespace accurate {

    // =========================================================================
    // 1. Pairwise Alignment Declarations
    // =========================================================================
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
        std::vector<int> score_mat;
        std::vector<uint8_t> used_ref_buf;
        std::vector<uint8_t> used_qry_buf;
        std::vector<std::pair<int, int>> sub_pair_buf;

        AlignmentResult align(const std::string& reference, const std::string& query, char type, Params& params);
        AlignmentResult align_affine(const std::string& reference, const std::string& query, char type, Params& params);
        AlignmentResult align_affine_local(const std::string& reference, const std::string& query, char type, Params& params);
        AlignmentResult align_linear_local(const std::string& reference, const std::string& query, char type, Params& params);
        LocalHomTable::PairResult align_affine_local_segments(const std::string& reference, const std::string& query, char type, Params& params);
        LocalHomTable::PairResult align_affine_local_segments_banded(const std::string& reference, const std::string& query, char type, Params& params, int bandWidth);
    };

} // namespace accurate

    // =========================================================================
    // 2. Profile-Profile Alignment Declarations
    // =========================================================================
    std::vector<int8_t> alignProfile(
        const AlnProfile& profile1,
        const AlnProfile& profile2,
        const Params& param,
        const std::vector<std::vector<float>>* consistencyTable = nullptr,
        float consistencyWeight = 0.0f,
        int bandWidth = -1
    );

    // Compatibility alias
    inline std::vector<int8_t> alignProfile_FAMSA(
        const AlnProfile& profile1,
        const AlnProfile& profile2,
        const Params& param,
        const std::vector<std::vector<float>>* consistencyTable = nullptr,
        float consistencyWeight = 0.0f,
        int bandWidth = -1)
    {
        return alignProfile(profile1, profile2, param, consistencyTable, consistencyWeight, bandWidth);
    }

} // namespace msa

#endif // ALIGN_METHOD_HPP
