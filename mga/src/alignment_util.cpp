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

Alignments splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts,const std::set<int>& qryCuts)  {
    Alignments frags; // 用來裝切碎的子片段

    // 1. 篩選並排序切點 (WGA Edge-Cut Fix)
    std::vector<int> rCuts;
    for (int c : refCuts) {
        // Ref 永遠是正向掃描，所以允許在終點 (second) 切割
        if (c > aln.refIdx.first && c <= aln.refIdx.second) rCuts.push_back(c);
    }

    std::vector<int> qCuts;
    for (int c : qryCuts) {
        if (aln.inverse) {
            // 反向掃描 (從 second 往 down 走到 first)，允許在終點 (first) 切割
            if (c >= aln.qryIdx.first && c < aln.qryIdx.second) qCuts.push_back(c);
        } else {
            // 正向掃描，允許在終點 (second) 切割
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

    // 2. 準備走訪 CIGAR 進行動態切割
    int rCutIdx = 0;
    int qCutIdx = 0;

    int rPos = aln.refIdx.first;
    int qPos = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first; 
    int qDir = aln.inverse ? -1 : 1;

    int currRStart = rPos;
    int currQStart = qPos;

    Alignment currAln = aln;
    currAln.CIGAR.clear(); 

    // 3. 逐一消耗 CIGAR Operations
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

            // 4. 檢查是否精準踩到切點
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
                    frags.push_back(currAln); // 改存進 frags
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

    // 5. 收尾
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
        frags.push_back(currAln); // 改存進 frags
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

    // 1. Ref 的補丁 (Ref 永遠是正向)
    if (r_pad_left > 0)  prepend_ops.push_back({r_pad_left, 'D'});
    if (r_pad_right > 0) append_ops.push_back({r_pad_right, 'D'});

    // 2. Qry 的補丁 (根據 Inverse 決定加在頭還是尾)
    if (aln.inverse) {
        if (q_pad_right > 0) prepend_ops.push_back({q_pad_right, 'I'});
        if (q_pad_left > 0)  append_ops.push_back({q_pad_left, 'I'});
    } else {
        if (q_pad_left > 0)  prepend_ops.push_back({q_pad_left, 'I'});
        if (q_pad_right > 0) append_ops.push_back({q_pad_right, 'I'});
    }

    // 3. 將 D/I 補丁貼上 CIGAR
    if (!prepend_ops.empty()) aln.CIGAR.insert(aln.CIGAR.begin(), prepend_ops.begin(), prepend_ops.end());
    if (!append_ops.empty())  aln.CIGAR.insert(aln.CIGAR.end(), append_ops.begin(), append_ops.end());

    // 4. 更新座標 (因為你保證了 first 永遠小於 second，所以直接加減即可)
    aln.refIdx.first  -= r_pad_left;
    aln.refIdx.second += r_pad_right;

    aln.qryIdx.first  -= q_pad_left;
    aln.qryIdx.second += q_pad_right;

    // 5. 更新對齊總長度
    aln.updateAlnLength(); 
}