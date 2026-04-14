#ifndef LOCAL_ALIGN_HPP
#define LOCAL_ALIGN_HPP

#include <memory>
#include <string>
#include <vector>

namespace msa {
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
    virtual LocalAlignmentResult align(const std::string& reference, const std::string& query, char type) const = 0;
};

std::shared_ptr<LocalAligner> makeDefaultLocalAligner();

} // namespace accurate
} // namespace msa

#endif
