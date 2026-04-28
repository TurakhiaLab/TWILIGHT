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
    AlignmentResult align(const std::string& reference, const std::string& query, char type, Params& params);
    AlignmentResult align_affine(const std::string& reference, const std::string& query, char type, Params& params);
    AlignmentResult align_affine_local(const std::string& reference, const std::string& query, char type, Params& params);
};

} // namespace accurate
} // namespace msa




#endif
