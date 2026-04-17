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

struct LocalAlignmentResult
{
    int score = 0;
    float identity = 0.0f;
    std::vector<AlignedResiduePair> alignedPairs;
};

class LocalAligner
{
public:
    virtual ~LocalAligner() = default;
    virtual LocalAlignmentResult align(const std::string& reference, const std::string& query, char type, Params& params) const = 0;
    virtual LocalAlignmentResult align_affine(const std::string& reference, const std::string& query, char type, Params& params) const = 0;
};

std::shared_ptr<LocalAligner> makeDefaultLocalAligner();

} // namespace accurate
} // namespace msa

std::vector<int8_t> alignProfile(
    const std::vector<std::vector<float>>& refProfile,
    const std::vector<std::vector<float>>& qryProfile,
    const std::vector<std::vector<float>>& gapOp,
    const std::vector<std::vector<float>>& gapEx,
    const std::pair<float, float>& num,
    msa::Params& param,
    const std::vector<std::vector<float>>* consistencyTable, // 允許為空指標
    float consistencyWeight
);


#endif
