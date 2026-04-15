#include "local-align.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>

namespace msa {
namespace accurate {
namespace {

class SmithWatermanAligner final : public LocalAligner
{
public:
    LocalAlignmentResult align(const std::string& reference, const std::string& query, char type, msa::Params& params) const override
    {
        // const int matchScore = (type == 'p') ? 3 : 2;
        // const int mismatchScore = (type == 'p') ? -2 : -1;
        // const int gapPenalty = -2;

        const int gapPenalty = static_cast<int>(params.gapExtend) * 5;
        

        const int refLen = static_cast<int>(reference.size());
        const int qryLen = static_cast<int>(query.size());
        std::vector<std::vector<int>> score(refLen + 1, std::vector<int>(qryLen + 1, 0));
        std::vector<std::vector<uint8_t>> traceback(refLen + 1, std::vector<uint8_t>(qryLen + 1, 0));

        int bestScore = 0;
        int bestRef = 0;
        int bestQry = 0;

        for (int i = 1; i <= refLen; ++i) {
            for (int j = 1; j <= qryLen; ++j) {
                const int diag = score[i - 1][j - 1] + residueScore(reference[i - 1], query[j - 1], type, params);
                const int up = score[i - 1][j] + gapPenalty;
                const int left = score[i][j - 1] + gapPenalty;

                int best = 0;
                uint8_t direction = 0;
                if (diag > best) {
                    best = diag;
                    direction = 1;
                }
                if (up > best) {
                    best = up;
                    direction = 2;
                }
                if (left > best) {
                    best = left;
                    direction = 3;
                }

                score[i][j] = best;
                traceback[i][j] = direction;

                if (best > bestScore) {
                    bestScore = best;
                    bestRef = i;
                    bestQry = j;
                }
            }
        }

        LocalAlignmentResult result;
        result.score = bestScore;
        for (int i = bestRef, j = bestQry; i > 0 && j > 0 && score[i][j] > 0; ) {
            const uint8_t direction = traceback[i][j];
            if (direction == 1) {
                result.alignedPairs.push_back({i - 1, j - 1});
                --i;
                --j;
            }
            else if (direction == 2) {
                --i;
            }
            else if (direction == 3) {
                --j;
            }
            else {
                break;
            }
        }

        std::reverse(result.alignedPairs.begin(), result.alignedPairs.end());
        int identicalPairs = 0;
        for (const auto& alignedPair : result.alignedPairs) {
            const char ref = static_cast<char>(std::toupper(static_cast<unsigned char>(reference[alignedPair.refIndex])));
            const char qry = static_cast<char>(std::toupper(static_cast<unsigned char>(query[alignedPair.qryIndex])));
            if (ref == qry) ++identicalPairs;
        }
        if (!result.alignedPairs.empty()) {
            result.identity = static_cast<float>(identicalPairs) / static_cast<float>(result.alignedPairs.size());
        }
        return result;
    }

private:
    static int residueScore(char referenceBase, char queryBase, char type, Params& params)
    {
        const char ref = static_cast<char>(std::toupper(static_cast<unsigned char>(referenceBase)));
        const char qry = static_cast<char>(std::toupper(static_cast<unsigned char>(queryBase)));
        int refIndex = letterIdx(type, toupper(referenceBase));
        int qryIndex = letterIdx(type, toupper(queryBase));
        return params.scoringMatrix[refIndex][qryIndex];
    }
};

} // namespace

std::shared_ptr<LocalAligner> makeDefaultLocalAligner()
{
    return std::make_shared<SmithWatermanAligner>();
}

} // namespace accurate
} // namespace msa
