#include "alignment.hpp"
#include "cigar_util.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <set>

CigarString extractSubCigar(const CigarString& origCigar, int refOffset, int refLen) {
    CigarString subCigar;
    int currentRef = 0;
    int extractedRef = 0;

    auto appendOp = [&](int len, char type) {
        if (len <= 0) return;
        if (!subCigar.empty() && subCigar.back().second == type) {
            subCigar.back().first += len;
        } else {
            subCigar.push_back({len, type});
        }
    };

    for (const auto& op : origCigar) {
        if (extractedRef >= refLen) break; 

        int opLen = op.first;
        char type = op.second;

        bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
        int refOpLen = consumesRef ? opLen : 0;
        
        if (consumesRef) {
            if (currentRef + refOpLen <= refOffset) {
                currentRef += refOpLen;
                continue;
            }
        } else {
            if (currentRef < refOffset) {
                continue;
            }
        }

        int useLen = opLen;

        if (consumesRef && currentRef < refOffset) {
            int trim = refOffset - currentRef;
            useLen -= trim;
            currentRef += trim;
        }

        if (consumesRef && (extractedRef + useLen > refLen)) {
            useLen = refLen - extractedRef;
        }

        appendOp(useLen, type);

        if (consumesRef) {
            currentRef += useLen;
            extractedRef += useLen;
        }
    }
    
    return subCigar;
}

Alignments splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts, const std::set<int>& qryCuts) {
    Alignments frags; // Stores split sub-fragments

    // 1. Filter and sort cut points (WGA Edge-Cut Fix)
    std::vector<int> rCuts;
    for (int c : refCuts) {
        // Reference is always forward-scanned, so cutting at the endpoint (second) is allowed
        if (c > aln.refIdx.first && c <= aln.refIdx.second) rCuts.push_back(c);
    }

    std::vector<int> qCuts;
    for (int c : qryCuts) {
        if (aln.inverse) {
            // Reverse scan (from second down to first), allow cutting at endpoint (first)
            if (c >= aln.qryIdx.first && c < aln.qryIdx.second) qCuts.push_back(c);
        } else {
            // Forward scan, allow cutting at endpoint (second)
            if (c > aln.qryIdx.first && c <= aln.qryIdx.second) qCuts.push_back(c);
        }
    }

    if (aln.inverse) {
        std::sort(qCuts.rbegin(), qCuts.rend());
    } else {
        std::sort(qCuts.begin(), qCuts.end());
    }

    if (rCuts.empty() && qCuts.empty()) {
        frags.push_back(aln);
        return frags;
    }

    // 2. Prepare CIGAR traversal for dynamic splitting
    int rCutIdx = 0;
    int qCutIdx = 0;

    int rPos = aln.refIdx.first;
    int qPos = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first; 
    int qDir = aln.inverse ? -1 : 1;

    int currRStart = rPos;
    int currQStart = qPos;

    Alignment currAln = aln;
    currAln.CIGAR.clear(); 

    // 3. Consume CIGAR operations step-by-step
    for (const auto& op : aln.CIGAR) {
        int len = op.first;
        char type = op.second;

        while (len > 0) {
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            int step = len;

            if (consumesRef && rCutIdx < rCuts.size()) {
                int distR = rCuts[rCutIdx] - rPos;
                if (distR > 0 && distR < step) step = distR;
            }

            if (consumesQry && qCutIdx < qCuts.size()) {
                int distQ = std::abs(qCuts[qCutIdx] - qPos);
                if (distQ > 0 && distQ < step) step = distQ;
            }

            if (!currAln.CIGAR.empty() && currAln.CIGAR.back().second == type) {
                currAln.CIGAR.back().first += step; 
            } else {
                currAln.CIGAR.push_back({step, type});
            }

            if (consumesRef) rPos += step;
            if (consumesQry) qPos += step * qDir;
            len -= step;

            // 4. Check if the current position hits a cut point exactly
            bool hitRef = (consumesRef && rCutIdx < rCuts.size() && rPos == rCuts[rCutIdx]);
            bool hitQry = (consumesQry && qCutIdx < qCuts.size() && qPos == qCuts[qCutIdx]);

            if (hitRef || hitQry) {
                currAln.refIdx.first = currRStart;
                currAln.refIdx.second = rPos;

                if (aln.inverse) {
                    currAln.qryIdx.first = qPos;        
                    currAln.qryIdx.second = currQStart; 
                } else {
                    currAln.qryIdx.first = currQStart;  
                    currAln.qryIdx.second = qPos;       
                }

                if (!currAln.CIGAR.empty()) {
                    frags.push_back(currAln); // Store into frags
                }

                currRStart = rPos;
                currQStart = qPos;
                currAln = aln; 
                currAln.CIGAR.clear();

                if (hitRef) rCutIdx++;
                if (hitQry) qCutIdx++;
            }
        }
    }

    // 5. Finalize the remaining fragment
    if (!currAln.CIGAR.empty()) {
        currAln.refIdx.first = currRStart;
        currAln.refIdx.second = rPos;

        if (aln.inverse) {
            currAln.qryIdx.first = qPos;
            currAln.qryIdx.second = currQStart;
        } else {
            currAln.qryIdx.first = currQStart;
            currAln.qryIdx.second = qPos;
        }
        frags.push_back(currAln); // Store into frags
    }

    // 6. Update alignment length
    for (auto& frag : frags) {
        frag.updateAlnLength();
    }

    return frags;
}

void snapAlignment(Alignment& aln, int r_pad_left, int r_pad_right, int q_pad_left, int q_pad_right) {
    std::vector<std::pair<int, char>> prepend_ops;
    std::vector<std::pair<int, char>> append_ops;

    // 1. Reference padding (reference is always forward)
    if (r_pad_left > 0)  prepend_ops.push_back({r_pad_left, 'D'});
    if (r_pad_right > 0) append_ops.push_back({r_pad_right, 'D'});

    // 2. Query padding (prepend or append depending on strand orientation)
    if (aln.inverse) {
        if (q_pad_right > 0) prepend_ops.push_back({q_pad_right, 'I'});
        if (q_pad_left > 0)  append_ops.push_back({q_pad_left, 'I'});
    } else {
        if (q_pad_left > 0)  prepend_ops.push_back({q_pad_left, 'I'});
        if (q_pad_right > 0) append_ops.push_back({q_pad_right, 'I'});
    }

    // 3. Apply D/I padding operations to CIGAR
    if (!prepend_ops.empty()) aln.CIGAR.insert(aln.CIGAR.begin(), prepend_ops.begin(), prepend_ops.end());
    if (!append_ops.empty())  aln.CIGAR.insert(aln.CIGAR.end(), append_ops.begin(), append_ops.end());

    // 4. Update coordinates (first is guaranteed to be <= second, so adjust directly)
    aln.refIdx.first  -= r_pad_left;
    aln.refIdx.second += r_pad_right;

    aln.qryIdx.first  -= q_pad_left;
    aln.qryIdx.second += q_pad_right;

    // 5. Update total alignment length
    aln.updateAlnLength(); 
}