#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include <tbb/parallel_for.h>

void msa::alignment_helper::calculateProfile(float *profile, NodePair &nodes, SequenceDB *database, Option *option, int32_t memLen)
{
    int profileSize = (option->type == 'n') ? 6 : 22;
    int refNum = nodes.first->getAlnNum(database->currentTask), qryNum = nodes.second->getAlnNum(database->currentTask);
    int refLen = nodes.first->getAlnLen(database->currentTask), qryLen = nodes.second->getAlnLen(database->currentTask);
    float refWeight = nodes.first->alnWeight, qryWeight = nodes.second->alnWeight;
    bool storeFreq = ((refNum >= _CAL_PROFILE_TH || qryNum >= _CAL_PROFILE_TH) || (!nodes.first->msaFreq.empty() || !nodes.second->msaFreq.empty()));
    { // Ref
        if (!nodes.first->msaFreq.empty()) {
            for (int t = 0; t < refLen; ++t) {
                for (int v = 0; v < profileSize; ++v)
                    profile[profileSize * t + v] = nodes.first->msaFreq[t][v] / refWeight * refNum;
            }
        }
        else {
            for (auto sIdx : nodes.first->seqsIncluded) {
                float w = database->sequences[sIdx]->weight / refWeight * refNum;;
                int storage = database->sequences[sIdx]->storage;
                tbb::this_task_arena::isolate([&] { 
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int t = r.begin(); t < r.end(); ++t) {
                    int letterIndex = letterIdx(option->type, toupper(database->sequences[sIdx]->alnStorage[storage][t]));
                    profile[profileSize*t+letterIndex] += 1.0 * w;  
                } 
                }); 
                });
            }
            if (storeFreq) {
                nodes.first->msaFreq.resize(refLen, std::vector<float>(profileSize, 0.0));
                for (int t = 0; t < refLen; ++t) 
                    for (int v = 0; v < profileSize; ++v)
                        nodes.first->msaFreq[t][v] = profile[profileSize * t + v] / refNum * refWeight;
            }
        }
    }
    { // Query
        if (!nodes.second->msaFreq.empty()) {
            for (int t = 0; t < qryLen; ++t) {
                for (int v = 0; v < profileSize; ++v)
                    profile[profileSize * (memLen+t) + v] = nodes.second->msaFreq[t][v] / qryWeight * qryNum;
            }
        }
        else {
            for (auto sIdx : nodes.second->seqsIncluded) {
                float w = database->sequences[sIdx]->weight / qryWeight * qryNum;;
                int storage = database->sequences[sIdx]->storage;
                tbb::this_task_arena::isolate([&] { 
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int t = r.begin(); t < r.end(); ++t) {
                    int letterIndex = letterIdx(option->type, toupper(database->sequences[sIdx]->alnStorage[storage][t]));
                    profile[profileSize*(memLen+t)+letterIndex] += 1.0 * w;
                } 
                }); 
                });
            }
            if (storeFreq) {
                nodes.second->msaFreq.resize(qryLen, std::vector<float>(profileSize, 0.0));
                for (int t = 0; t < qryLen; ++t)
                    for (int v = 0; v < profileSize; ++v)
                        nodes.second->msaFreq[t][v] = profile[profileSize * (memLen + t) + v] / qryNum * qryWeight;
            }
        }
    }
    return;
}

void msa::alignment_helper::removeGappyColumns(float *hostFreq, NodePair &nodes, Option *option, std::pair<IntPairVec, IntPairVec> &gappyColumns, int32_t memLen, IntPair &lens, int currentTask)
{
    float gap_threshold = option->gappyVertical;
    if (gap_threshold == 1.0) return;
    int32_t profileSize = (option->type == 'n') ? 6 : 22;
    int refLen = nodes.first->getAlnLen(currentTask), qryLen = nodes.second->getAlnLen(currentTask);
    int refNum = nodes.first->getAlnNum(currentTask), qryNum = nodes.second->getAlnNum(currentTask);
    int start = -1, length = 0;
    // Reference
    for (int i = 0; i < lens.first; ++i) {
        if ((hostFreq[profileSize * i + profileSize - 1]) / refNum > gap_threshold) {
            if (start == -1) {
                start = i; // new region starts
                length = 1;
            }
            else {
                ++length; // continue region
            }
        }
        else if (start != -1) {
            gappyColumns.first.push_back({start, length});
            start = -1;
            length = 0;
        }
    }
    if (start != -1) { // handle case where sequence ends inside a region
        gappyColumns.first.push_back({start, length});
    }
    // Query
    start = -1, length = 0;
    for (int i = 0; i < lens.second; ++i) {
        if ((hostFreq[profileSize * (memLen + i) + profileSize - 1]) / qryNum > gap_threshold)
        {
            if (start == -1) {
                start = i; // new region starts
                length = 1;
            }
            else {
                ++length; // continue region
            }
        }
        else if (start != -1) {
            gappyColumns.second.push_back({start, length});
            start = -1;
            length = 0;
        }
    }
    if (start != -1) {
        gappyColumns.second.push_back({start, length});
    }
    // Remove gappy columns
    int orgIdx = 0, newIdx = 0, gapIdx = 0, orgLen = lens.first;
    if (!gappyColumns.first.empty()) {
        while (orgIdx < orgLen) {
            bool gapCols = (gapIdx >= gappyColumns.first.size()) ? false : (orgIdx == gappyColumns.first[gapIdx].first);
            if (gapCols) {
                orgIdx += gappyColumns.first[gapIdx].second;
                gapIdx++;
            }
            else {
                for (int t = 0; t < profileSize; ++t) hostFreq[profileSize * newIdx + t] = hostFreq[profileSize * orgIdx + t];
                ++newIdx;
                ++orgIdx;
            }
        }
        lens.first = newIdx;
        while (newIdx < orgLen) {
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize * newIdx + t] = 0;
            newIdx++;
        }
    }
    orgIdx = 0, newIdx = 0, gapIdx = 0, orgLen = lens.second;
    if (!gappyColumns.second.empty()) {
        while (orgIdx < orgLen) {
            bool gapCols = (gapIdx >= gappyColumns.second.size()) ? false : (orgIdx == gappyColumns.second[gapIdx].first);
            if (gapCols) {
                orgIdx += gappyColumns.second[gapIdx].second;
                gapIdx++;
            }
            else {
                for (int t = 0; t < profileSize; ++t) hostFreq[profileSize * (memLen + newIdx) + t] = hostFreq[profileSize * (memLen + orgIdx) + t];
                ++newIdx;
                ++orgIdx;
            }
        }
        lens.second = newIdx;
        while (newIdx < orgLen) {
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize * (memLen + newIdx) + t] = 0;
            newIdx++;
        }
    }
    return;
}

void msa::alignment_helper::calculatePSGP(float *hostFreq, float *hostGapOp, float *hostGapEx, NodePair &nodes, SequenceDB* database, Option *option, int memLen, IntPair offset, IntPair lens, Params &param)
{
    int32_t refLen = lens.first;
    int32_t qryLen = lens.second;
    int32_t offsetf = offset.first;
    int32_t offsetg = offset.second;
    int32_t refNum = nodes.first->getAlnNum(database->currentTask);
    int32_t qryNum = nodes.second->getAlnNum(database->currentTask);
    int32_t profileSize = (option->type == 'n') ? 6 : 22;

    // Clustalw's method
    float scale = (option->type == 'n') ? 0.5 : 1.0;
    float min_gapExtend = param.gapExtend * 0.2;
    float min_gapOpen = param.gapOpen * 0.1;
    tbb::this_task_arena::isolate([&]  { 
    tbb::parallel_for(tbb::blocked_range<int>(0, memLen), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        if (s < refLen) {
            float gapRatio = hostFreq[offsetf+profileSize*s+profileSize-1];
            if (gapRatio > 0) {
                hostGapOp[offsetg+s] = std::min(min_gapOpen,   static_cast<float>(param.gapOpen * scale * ((refNum-gapRatio)*1.0 / refNum)));
                hostGapEx[offsetg+s] = std::min(min_gapExtend, static_cast<float>(param.gapExtend * ((refNum-gapRatio)*1.0 / refNum)));
            }
            else {
                hostGapOp[offsetg+s] = param.gapOpen;
                hostGapEx[offsetg+s] = param.gapExtend;
            }
        }
        else {
            hostGapOp[offsetg+s] = 0.0;
            hostGapEx[offsetg+s] = 0.0;
        }
        if (s < qryLen) {
            float gapRatio = hostFreq[offsetf+profileSize*(memLen+s)+profileSize-1];
            if (gapRatio > 0) {
                hostGapOp[offsetg+memLen+s] = std::min(min_gapOpen,  static_cast<float>(param.gapOpen * scale * ((qryNum-gapRatio) * 1.0 / qryNum)));
                hostGapEx[offsetg+memLen+s] = std::min(min_gapExtend, static_cast<float>(param.gapExtend * ((qryNum-gapRatio) * 1.0 / qryNum)));
            }
            else {
                hostGapOp[offsetg+memLen+s] = param.gapOpen;
                hostGapEx[offsetg+memLen+s] = param.gapExtend;
            }
        }
        else {
            hostGapOp[offsetg+memLen+s] = 0.0;
            hostGapEx[offsetg+memLen+s] = 0.0;
        }     
    } 
    }); 
    });
    return;
}

void msa::alignment_helper::getConsensus(Option *option, float *profile, std::string &consensus, int len)
{
    const char bases[5] = {'A', 'C', 'G', 'T', 'N'};
    const char acids[21] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};
    const char *lut = (option->type == 'n') ? bases : acids;
    const int profileSize = (option->type == 'n') ? 6 : 22;
    consensus.reserve(len);
    for (int i = 0; i < len; ++i) {
        int max_idx = profileSize - 2;
        float max_count = 0;
        for (int j = 0; j < profileSize - 2; ++j) {
            if (profile[profileSize * i + j] > max_count)
            {
                max_count = profile[profileSize * i + j];
                max_idx = j;
            }
        }
        consensus.append(std::string(1, lut[max_idx]));
    }
    return;
}

void msa::alignment_helper::pairwiseGlobal(const std::string &seq1, const std::string &seq2, alnPath &alnPath, Params &param)
{
    int m = seq1.size();
    int n = seq2.size();
    char type = (param.matrixSize == 5) ? 'n' : 'p';
    float gap_open = param.gapOpen, gap_extend = param.gapExtend;
    // DP matrices: M = match/mismatch, X = gap in seq1 (deletion), Y = gap in seq2 (insertion)
    std::vector<std::vector<float>> M(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<float>> X(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<float>> Y(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<int8_t>> tb(m + 1, std::vector<int8_t>(n + 1, 0)); // traceback: 0/1/2
    // Initialization
    M[0][0] = 0;
    for (int i = 1; i <= m; ++i)
    {
        // M[i][0] = gap_open + (i - 1) * gap_extend;
        M[i][0] = 0;
        X[i][0] = M[i][0];
        Y[i][0] = -1e9; // impossible
        tb[i][0] = 2;   // gap in query
    }
    for (int j = 1; j <= n; ++j)
    {
        // M[0][j] = gap_open + (j - 1) * gap_extend;
        M[0][j] = 0;
        Y[0][j] = M[0][j];
        X[0][j] = -1e9;
        tb[0][j] = 1; // gap in reference
    }
    // Fill matrices
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            // Compute scores
            int idx_1 = letterIdx(type, toupper(seq1[i - 1])), idx_2 = letterIdx(type, toupper(seq2[j - 1]));
            float base_score = param.scoringMatrix[idx_1][idx_2];
            // Match/mismatch
            M[i][j] = base_score + std::max({M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]});
            // Gap in seq1 (deletion)
            X[i][j] = std::max(M[i - 1][j] + gap_open,
                               X[i - 1][j] + gap_extend);
            // Gap in seq2 (insertion)
            Y[i][j] = std::max(M[i][j - 1] + gap_open,
                               Y[i][j - 1] + gap_extend);
            // Choose best path for traceback
            float best = std::max({M[i][j], X[i][j], Y[i][j]});
            if (best == M[i][j])
                tb[i][j] = 0;
            else if (best == Y[i][j])
                tb[i][j] = 1; // gap in ref (insert in query)
            else
                tb[i][j] = 2; // gap in qry (delete in ref)
        }
    }
    // Traceback
    alnPath.clear();
    int i = m, j = n;
    float max_score = std::max({M[m][n], X[m][n], Y[m][n]});
    while (i > 0 || j > 0)
    {
        int8_t dir = tb[i][j];
        alnPath.push_back(dir);
        if (dir == 0)
        {
            i--;
            j--;
        }
        else if (dir == 1)
        {
            j--;
        }
        else if (dir == 2)
        {
            i--;
        }
    }
    std::reverse(alnPath.begin(), alnPath.end());
    return;
}

void msa::alignment_helper::addGappyColumnsBack(alnPath &aln_before, alnPath &aln_after, std::pair<IntPairVec, IntPairVec> &gappyColumns, Params &param, IntPair rgcLens, stringPair orgSeqs)
{
    int rIdx = 0, qIdx = 0, alnIdx = 0;
    int gc_rIdx = 0, gc_qIdx = 0;
    for (alnIdx = 0; alnIdx < aln_before.size()+1; ++alnIdx) {
        bool gapR = (gc_rIdx >= gappyColumns.first.size())  ? false : (rIdx == gappyColumns.first[gc_rIdx].first);
        bool gapQ = (gc_qIdx >= gappyColumns.second.size()) ? false : (qIdx == gappyColumns.second[gc_qIdx].first);
        if (gapR && gapQ) {
            int gc_len_r = gappyColumns.first[gc_rIdx].second;
            int gc_len_q = gappyColumns.second[gc_qIdx].second;
            std::string consensus_r = orgSeqs.first.substr(rIdx, gc_len_r);
            std::string consensus_q = orgSeqs.second.substr(qIdx, gc_len_q);
            alnPath aln_temp(0);
            pairwiseGlobal(consensus_r, consensus_q, aln_temp, param);
            for (auto a : aln_temp) aln_after.push_back(a);
            ++gc_rIdx; ++gc_qIdx;
            rIdx += gc_len_r; qIdx += gc_len_q;
        }
        else {
            if (gapR) {
                int gc_len = gappyColumns.first[gc_rIdx].second;
                aln_after.insert(aln_after.end(), gc_len, 2);
                rIdx += gc_len;
                gc_rIdx++;
            }
            if (gapQ) {
                int gc_len = gappyColumns.second[gc_qIdx].second;
                aln_after.insert(aln_after.end(), gc_len, 1);
                qIdx += gc_len;
                gc_qIdx++;
            }
        }
        if (alnIdx < aln_before.size()) {
            aln_after.push_back(aln_before[alnIdx]);
            switch (aln_before[alnIdx]) {
            case 0:
                ++rIdx; ++qIdx;
                break;
            case 1: 
                ++qIdx;
                break;
            case 2:
                ++rIdx;
                break;
            default:
                std::cerr << "ERROR: Undefined TB Path: " << aln_before[alnIdx] << '\n';
                break;
            }
        }
    }
    return;
}

void msa::alignment_helper::updateAlignment(NodePair &nodes, SequenceDB *database, Option* option, alnPath &aln)
{
    int totalLen = aln.size(), startLen = 0, endLen = 0;
    tbb::this_task_arena::isolate([&] { 
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.first->seqsIncluded.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = nodes.first->seqsIncluded[idx];
        if (database->currentTask != 2 && sIdx >= 0) {
            database->sequences[sIdx]->memCheck(totalLen);
            int storeFrom = database->sequences[sIdx]->storage;
            int storeTo = 1 - storeFrom;
            int orgIdx = 0;
            for (int k = 0; k < aln.size(); ++k) {
                int kto = k + startLen;
                if (aln[k] == 0 || aln[k] == 2) {
                    database->sequences[sIdx]->alnStorage[storeTo][kto] = database->sequences[sIdx]->alnStorage[storeFrom][orgIdx];
                    orgIdx++;
                }
                else {
                    database->sequences[sIdx]->alnStorage[storeTo][kto] = '-';
                }
            }
            database->sequences[sIdx]->len = totalLen;
            database->sequences[sIdx]->changeStorage();
        }
        else {
            int orgIdx = 0;
            alnPath orgAln = database->subtreeAln[sIdx];
            database->subtreeAln[sIdx].resize(aln.size());
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 2) {
                    database->subtreeAln[sIdx][k] = orgAln[orgIdx];
                    orgIdx++;
                }
                else {
                    database->subtreeAln[sIdx][k] = 1; // 1 for Gaps
                }
            }
        }
    } 
    }); 
    });
    tbb::this_task_arena::isolate([&] { 
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.second->seqsIncluded.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = nodes.second->seqsIncluded[idx];
        if (database->currentTask != 2 && sIdx >= 0) {
            database->sequences[sIdx]->memCheck(totalLen);
            int orgIdx = 0;
            int storeFrom = database->sequences[sIdx]->storage;
            int storeTo = 1 - storeFrom;
            for (int k = 0; k < aln.size(); ++k) {
                int kto = k + startLen;
                if (aln[k] == 0 || aln[k] == 1) {
                    database->sequences[sIdx]->alnStorage[storeTo][kto] = database->sequences[sIdx]->alnStorage[storeFrom][orgIdx];
                    orgIdx++;
                }
                else {
                    database->sequences[sIdx]->alnStorage[storeTo][kto] = '-';
                }
            }
            database->sequences[sIdx]->len = totalLen;
            database->sequences[sIdx]->changeStorage();
        }
        else {
            int orgIdx = 0;
            alnPath orgAln = database->subtreeAln[sIdx];
            database->subtreeAln[sIdx].resize(aln.size());
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 1) {
                    database->subtreeAln[sIdx][k] = orgAln[orgIdx];
                    orgIdx++;
                }
                else {
                    database->subtreeAln[sIdx][k] = 1; // 1 for Gaps
                }
            }
        }
    } 
    }); 
    });
    nodes.first->alnNum += nodes.second->alnNum;
    nodes.first->alnLen = totalLen;
    nodes.first->alnWeight += nodes.second->alnWeight;
    for (auto idx: nodes.second->seqsIncluded) nodes.first->seqsIncluded.push_back(idx);
    nodes.second->seqsIncluded.clear();
    if (nodes.first->seqsIncluded.size() > _UPDATE_SEQ_TH && !nodes.first->msaFreq.empty() && database->currentTask != 2) {
        int seqCount = 0, firstSeqID = 0;
        for (auto idx: nodes.first->seqsIncluded) {
            if (idx > 1) {
               if (firstSeqID == 0) firstSeqID = -idx;
               seqCount++;
            }
        }
        if (seqCount >= _UPDATE_SEQ_TH) {
            database->subtreeAln[firstSeqID] = alnPath(totalLen, 0);
            std::vector<int> new_seqsIncluded;
            new_seqsIncluded.push_back(firstSeqID);
            for (auto idx: nodes.first->seqsIncluded) {
                if (idx >= 0) database->sequences[idx]->subtreeIdx = firstSeqID;
                else new_seqsIncluded.push_back(idx);
            }
            nodes.first->seqsIncluded = new_seqsIncluded;
        }
    }
    
    return;
}

void msa::alignment_helper::updateFrequency(NodePair &nodes, SequenceDB *database, alnPath &aln, FloatPair weights)
{
    if (nodes.first->msaFreq.empty() || nodes.second->msaFreq.empty()) return;
    int profileSize = nodes.first->msaFreq[0].size();
    float refWeight = weights.first, qryWeight = weights.second;
    Profile mergeFreq(aln.size(), std::vector<float>(profileSize, 0.0));
    int rIdx = 0, qIdx = 0;
    for (int j = 0; j < aln.size(); ++j) {
        if (aln[j] == 0) {
            for (int k = 0; k < profileSize; ++k)
                mergeFreq[j][k] = (nodes.first->msaFreq[rIdx][k] + nodes.second->msaFreq[qIdx][k]);
            ++rIdx;
            ++qIdx;
        }
        else if (aln[j] == 1) {
            for (int k = 0; k < profileSize - 1; ++k)
                mergeFreq[j][k] = nodes.second->msaFreq[qIdx][k];
            mergeFreq[j][profileSize - 1] = (nodes.second->msaFreq[qIdx][profileSize - 1] + 1.0 * refWeight);
            ++qIdx;
        }
        else if (aln[j] == 2) {
            for (int k = 0; k < profileSize - 1; ++k)
                mergeFreq[j][k] = nodes.first->msaFreq[rIdx][k];
            mergeFreq[j][profileSize - 1] = (nodes.first->msaFreq[rIdx][profileSize - 1] + 1.0 * qryWeight);
            ++rIdx;
        }
    }
    // Free memory
    nodes.first->msaFreq.clear();
    nodes.second->msaFreq.clear();
    nodes.first->msaFreq = mergeFreq;
    nodes.first->alnLen = mergeFreq.size();
    return;
}

void msa::alignment_helper::fallback2cpu(std::vector<int>& fallbackPairs,  NodePairVec& nodes, SequenceDB* database, Option* option) {
    int totalSeqs = 0;
    bool filtering = !option->noFilter;
    std::sort(fallbackPairs.begin(), fallbackPairs.end());
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int32_t refNum = nodes[nIdx].first->alnNum;
        int32_t qryNum = nodes[nIdx].second->alnNum;
        bool lowQ_r = (refNum > 1) ? false : (database->sequences[nodes[nIdx].first->seqsIncluded[0]]->lowQuality);
        bool lowQ_q = (qryNum > 1) ? false : (database->sequences[nodes[nIdx].second->seqsIncluded[0]]->lowQuality);
        if ((refNum < qryNum) || lowQ_r) {
            bool fallback = (!filtering) || (!lowQ_r);
            if (fallback) {
                database->fallback_nodes.push_back(nodes[nIdx].second);
                if (lowQ_r) database->sequences[nodes[nIdx].first->seqsIncluded[0]]->lowQuality = false;
            }
            // swap reference and query
            int32_t refLen = nodes[nIdx].first->alnLen;
            int32_t qryLen = nodes[nIdx].second->alnLen;
            float refWeight = nodes[nIdx].first->alnWeight;
            float qryWeight = nodes[nIdx].second->alnWeight;
            nodes[nIdx].second->alnLen = refLen;
            nodes[nIdx].first->alnLen = qryLen;
            nodes[nIdx].second->alnNum = refNum;
            nodes[nIdx].first->alnNum = qryNum;
            nodes[nIdx].second->alnWeight = refWeight;
            nodes[nIdx].first->alnWeight = qryWeight;
            auto temp = nodes[nIdx].second->seqsIncluded;
            nodes[nIdx].second->seqsIncluded = nodes[nIdx].first->seqsIncluded;
            nodes[nIdx].first->seqsIncluded = temp;
            auto freq_temp = nodes[nIdx].second->msaFreq;
            nodes[nIdx].second->msaFreq = nodes[nIdx].first->msaFreq;
            nodes[nIdx].first->msaFreq = freq_temp;
            totalSeqs += refNum;
        }
        else {
            bool fallback = (!filtering) || (!lowQ_q);
            if (fallback) {
                database->fallback_nodes.push_back(nodes[nIdx].second);
                if (lowQ_q) database->sequences[nodes[nIdx].second->seqsIncluded[0]]->lowQuality = false;
            }
            totalSeqs += qryNum;
        }
    }
    if (option->printDetail) {
        if (fallbackPairs.size() == 1 && totalSeqs == 1) std::cout << "Deferring/excluding 1 pair (1 sequence).\n";
        else if (fallbackPairs.size() == 1 && totalSeqs > 1) printf("Deferring/excluding 1 pair (%d sequences).\n", totalSeqs); 
        else printf("Deferring/excluding %lu pair (%d sequences).\n", fallbackPairs.size(), totalSeqs); 
    }
    return;
}

void msa::alignment_helper::mergeInsertions(SequenceDB *database, Node* root) {
    
    std::vector<std::unordered_map<int,int>> insertions (database->subtreeAln.size());
    int refLen = database->subtreeAln[-1].size();
    tbb::parallel_for(tbb::blocked_range<int>(0, database->subtreeAln.size()), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        if (database->subtreeAln.find(s) != database->subtreeAln.end() && !database->sequences[s]->lowQuality) {
            int refIdx = 0,  start = -1, len = 0;
            for (int i = 0; i < database->subtreeAln[s].size(); ++i) {
                if (database->subtreeAln[s][i] == 1) {
                    if (start == -1) start = refIdx;
                    len++;
                }
                else {
                    if (start != -1) {
                        insertions[s][start] = len;
                        start = -1;
                        len = 0;
                    }
                    refIdx++;
                }
            }
            if (start != -1) insertions[s][start] = len;
            if (refIdx != refLen) {
                std::cerr << "ERROR: Length not match.\t" <<database->subtreeAln[s].size() << ","<< refLen << '-' << refIdx << '\n';
            }
        }
    }
    });

    std::vector<int> longest_insertions (refLen+1, 0);
    tbb::parallel_for(tbb::blocked_range<int>(0, refLen+1), [&](tbb::blocked_range<int> r) {
    for (int i = r.begin(); i < r.end(); ++i) {
        int maxLen = 0;
        for (int s = 0; s < database->subtreeAln.size(); ++s) {
            if (insertions[s].find(i) != insertions[s].end()) {
                maxLen = std::max(maxLen, insertions[s][i]);
            }
        }
        longest_insertions[i] = maxLen;
    }
    });

    int totalLen = refLen;
    for (auto l: longest_insertions) totalLen += l;
    alnPath refAln;
    refAln.reserve(totalLen);
    for (int i = 0; i < longest_insertions.size(); ++i) {
        int l = longest_insertions[i];
        for (int i = 0; i < l; ++i) refAln.push_back(3);
        if (i < longest_insertions.size() - 1) refAln.push_back(0);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, database->sequences.size()), [&](tbb::blocked_range<int> range) {
    for (int sIdx = range.begin(); sIdx < range.end(); ++sIdx) {
        if (database->sequences[sIdx]->lowQuality) continue;
        database->sequences[sIdx]->memCheck(totalLen);
        int orgIdx = 0, alnIdx = 0;
        int storeFrom = database->sequences[sIdx]->storage;
        int storeTo = 1 - storeFrom;
        for (int k = 0; k < totalLen; ++k) {
            if (refAln[k] == 0) {
                if (database->subtreeAln[sIdx][alnIdx] == 0) {
                    database->sequences[sIdx]->alnStorage[storeTo][k] = database->sequences[sIdx]->alnStorage[storeFrom][orgIdx];
                    ++alnIdx; ++orgIdx;
                }
                else if (database->subtreeAln[sIdx][alnIdx] == 2) {
                    database->sequences[sIdx]->alnStorage[storeTo][k] = '-';
                    ++alnIdx;
                }
                else {
                    std::cerr << "ERROR: " << k << ": " << database->subtreeAln[sIdx][alnIdx] << '\n';
                }
            }
            else { // refAln[k] == 3
                if (alnIdx >= database->subtreeAln[sIdx].size()) {
                    database->sequences[sIdx]->alnStorage[storeTo][k] = '.';
                }
                else {
                    if (database->subtreeAln[sIdx][alnIdx] == 1) {
                        database->sequences[sIdx]->alnStorage[storeTo][k] = database->sequences[sIdx]->alnStorage[storeFrom][orgIdx];
                        alnIdx++; ++orgIdx;
                    }
                    else {
                        database->sequences[sIdx]->alnStorage[storeTo][k] = '.';
                    }
                }
                
            }
            
        }
        database->sequences[sIdx]->changeStorage();
    }
    });
    database->subtreeAln[-1] = refAln;
    root->alnLen = totalLen;
    return;
}




/*
void msa::alignment_helper::getStartingAndEndingPoints(std::string& ref, std::string& qry, IntPair& starting, IntPair& ending) {
    const float ratio = 0.1;
    int segLenRef = ref.size() * ratio, segLenQry = qry.size() * ratio;
    std::string rStart = ref.substr(0, segLenRef), rEnd = ref.substr(ref.size()-segLenRef, segLenRef);
    std::string qStart = qry.substr(0, segLenQry), qEnd = qry.substr(qry.size()-segLenQry, segLenQry);
    std::reverse(rStart.begin(),rStart.end());
    std::reverse(qStart.begin(), qStart.end());
    starting = smith_waterman(rStart, qStart);
    starting.first = segLenRef - starting.first;
    starting.second = segLenQry - starting.second;
    ending = smith_waterman(rEnd,qEnd);
    ending.first = ref.size()-segLenRef + ending.first;
    ending.second = qry.size()-segLenQry + ending.second;
    return;
}

void msa::alignment_helper::alignEnds (NodePair &nodes, SequenceDB *database, char type, IntPair& starting, IntPair& ending, std::map<int, std::string>& startAln, std::map<int, std::string>& endAln) 
{
    // Extract Sequence and remove gaps
    int refNum = nodes.first->getAlnNum(database->currentTask);
    int qryNum = nodes.second->getAlnNum(database->currentTask);
    int refLen = nodes.first->getAlnLen(database->currentTask);
    int qryLen = nodes.second->getAlnLen(database->currentTask);
    std::vector<std::pair<int, std::string>> startingSeqs (refNum+qryNum), endingSeqs(refNum+qryNum);
    // Reference
    tbb::this_task_arena::isolate([&] { 
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.first->seqsIncluded.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = nodes.first->seqsIncluded[idx];
        int storage = database->sequences[sIdx]->storage;
        startingSeqs[idx] = {sIdx, std::string(database->sequences[sIdx]->alnStorage[storage], starting.first)};
        // std::cout << starting.first << ',' << std::string(database->sequences[sIdx]->alnStorage[storage], starting.first) << ',' << std::string(database->sequences[sIdx]->alnStorage[storage], 10) << '\n';
        endingSeqs[idx]   = {sIdx, std::string(database->sequences[sIdx]->alnStorage[storage]+ending.first, refLen-ending.first)};
        size_t k = 0;
        for (size_t i = 0; i < startingSeqs[idx].second.size(); ++i) {
            if (startingSeqs[idx].second[i] != '-') startingSeqs[idx].second[k++] = startingSeqs[idx].second[i];
        }
        startingSeqs[idx].second.resize(k);
        k = 0;
        for (size_t i = 0; i < endingSeqs[idx].second.size(); ++i) {
            if (endingSeqs[idx].second[i] != '-') endingSeqs[idx].second[k++] = endingSeqs[idx].second[i];
        }
        endingSeqs[idx].second.resize(k);
    }
    });
    });
    tbb::this_task_arena::isolate([&] { 
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.second->seqsIncluded.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = nodes.second->seqsIncluded[idx];
        int storage = database->sequences[sIdx]->storage;
        startingSeqs[idx+refNum] = {sIdx, std::string(database->sequences[sIdx]->alnStorage[storage], starting.second)};
        endingSeqs[idx+refNum]   = {sIdx, std::string(database->sequences[sIdx]->alnStorage[storage]+ending.second, qryLen-ending.second)};
        size_t k = 0;
        for (size_t i = 0; i < startingSeqs[idx+refNum].second.size(); ++i) {
            if (startingSeqs[idx+refNum].second[i] != '-') startingSeqs[idx+refNum].second[k++] = startingSeqs[idx+refNum].second[i];
        }
        startingSeqs[idx+refNum].second.resize(k);
        k = 0;
        for (size_t i = 0; i < endingSeqs[idx+refNum].second.size(); ++i) {
            if (endingSeqs[idx+refNum].second[i] != '-') endingSeqs[idx+refNum].second[k++] = endingSeqs[idx+refNum].second[i];
        }
        endingSeqs[idx+refNum].second.resize(k);
    }
    });
    });
    std::sort(startingSeqs.begin(), startingSeqs.end(), [&](const std::pair<int,std::string> &a, const std::pair<int,std::string> &b) {
        if (a.second.size() == b.second.size()) return a.first < b.first;
        return a.second.size() > b.second.size(); // descending
    });
    std::sort(endingSeqs.begin(), endingSeqs.end(), [&](const std::pair<int,std::string> &a, const std::pair<int,std::string> &b) {
        if (a.second.size() == b.second.size()) return a.first < b.first;
        return a.second.size() > b.second.size(); // descending
    });
    // Alignment of starting
    if (starting.first > 0 || starting.second > 0) {
        int consensusLen = startingSeqs[0].second.size();
        for (auto s: startingSeqs) startAln[s.first] = s.second;
        for (int i = 1; i < startingSeqs.size(); ++i) {
            std::string seqConsensus = "";
            for (int k = 0; k < consensusLen; ++k) {
                int charNum [6] = {0};
                for (int j = 0; j < i; ++j) {
                    int rsIdx = startingSeqs[j].first;
                    charNum[letterIdx('n', toupper(startAln[rsIdx][k]))] += 1;
                }
                int idx = std::max_element(charNum, charNum + 5) - charNum;
                if (idx == 0) seqConsensus += "A";
                else if (idx == 1) seqConsensus += "C";
                else if (idx == 2) seqConsensus += "G";
                else if (idx == 3) seqConsensus += "T";
                else               seqConsensus += "A";
            }
            auto seq = startingSeqs[i].second;
            // std::reverse(seq.begin(), seq.end());
            auto sIdx = startingSeqs[i].first;
            std::cout << "Consen: " << seqConsensus << '\n';
            std::cout << "Added : " << seq << '\n';
            alnPath tempAln = smith_waterman_tb(seqConsensus, seq);
            // std::reverse(tempAln.begin(), tempAln.end());
            for (auto a: tempAln) std::cout << (a & 0xFFFF);
            std::cout << '\n';
            // Update aln path and consensus
            std::string seqTemp = "";
            seqTemp.reserve(tempAln.size());
            // Update the consensus
            // tbb::this_task_arena::isolate([&] { 
            // tbb::parallel_for(tbb::blocked_range<int>(0, i), [&](tbb::blocked_range<int> range) {
            // for (int j = range.begin(); j < range.end(); ++j) {
            for (int j = 0; j < i; ++j) {
                int rsIdx = startingSeqs[j].first, rsSeqIdx = 0;
                std::string temp = "";
                temp.reserve(tempAln.size());
                for (auto a: tempAln) {
                    if (a == 0 || a == 2) {
                        temp += startAln[rsIdx][rsSeqIdx];
                        rsSeqIdx++;
                    }
                    else if (a == 1) {
                        temp += '-';
                    }
                    // std::cout << (a & 0xFFFF) << ',' << temp << '\n';
                }
                startAln[rsIdx] = temp;
            }
            // });
            // });
            // Update the newly aligned seq
            {
                int sqIdx = 0;
                std::string temp = "";
                temp.reserve(tempAln.size());
                for (auto a: tempAln) {
                    if (a == 0 || a == 1) {
                        temp += startAln[sIdx][sqIdx];
                        sqIdx++; 
                    }
                    else if (a == 2) {
                        temp += '-';
                    }
                    // std::cout << (a & 0xFFFF) << ',' << temp << '\n';
                }
                startAln[sIdx] = temp;
                consensusLen = tempAln.size();
            }
        }
        for (auto a: startAln) {
            std::cout << a.first << '\t' << a.second << '\n';
        }
    }
    
    if (ending.first < refLen || ending.second < qryLen) {
        // Alignment of ending
        int consensusLen = endingSeqs[0].second.size();
        for (auto s: endingSeqs) endAln[s.first] = s.second;
        for (int i = 1; i < endingSeqs.size(); ++i) {
            std::string seqConsensus = "";
            for (int k = 0; k < consensusLen; ++k) {
                int charNum [6] = {0};
                for (int j = 0; j < i; ++j) {
                    int rsIdx = startingSeqs[j].first;
                    charNum[letterIdx('n', toupper(startAln[rsIdx][k]))] += 1;
                }
                int idx = std::max_element(charNum, charNum + 5) - charNum;
                if (idx == 0) seqConsensus += "A";
                else if (idx == 1) seqConsensus += "C";
                else if (idx == 2) seqConsensus += "G";
                else if (idx == 3) seqConsensus += "T";
                else               seqConsensus += "A";
            }
            auto seq = endingSeqs[i].second;
            auto sIdx = endingSeqs[i].first;
            std::cout << "Consen: " << seqConsensus << '\n';
            std::cout << "Added : " << seq << '\n';
            alnPath tempAln = smith_waterman_tb(seqConsensus, seq);
            for (auto a: tempAln) std::cout << (a & 0xFFFF);
            std::cout << '\n';
            // Update aln path and consensus
            std::string seqTemp = "";
            seqTemp.reserve(tempAln.size());
            // Update the consensus
            // tbb::this_task_arena::isolate([&] { 
            // tbb::parallel_for(tbb::blocked_range<int>(0, i), [&](tbb::blocked_range<int> range) {
            // for (int j = range.begin(); j < range.end(); ++j) {
            for (int j = 0; j < i; ++j) {
                int rsIdx = endingSeqs[j].first, rsAlnIdx = 0;
                std::string temp = "";
                temp.reserve(tempAln.size());
                for (auto a: tempAln) {
                    if (a == 0 || a == 2) {
                        temp += endAln[rsIdx][rsAlnIdx];
                        rsAlnIdx++;
                    }
                    else if (a == 1) {
                        temp += '-';
                    }
                    // std::cout << (a & 0xFFFF) << ',' << temp << '\n';
                }
                endAln[rsIdx] = temp;
            }
            // });
            // });
            // Update the newly aligned seq
            {
                int sqIdx = 0;
                std::string temp;
                temp.reserve(tempAln.size());
                for (auto a: tempAln) {
                    if (a == 0 || a == 1) {
                        temp += endAln[sIdx][sqIdx];
                        sqIdx++; 
                    }
                    else if (a == 2) {
                        temp += '-';
                    }
                }
                endAln[sIdx] = temp;
                consensusLen = tempAln.size();
            }
        }
        for (auto a: endAln) {
            std::cout << a.first << '\t' << a.second << '\n';
        }
    }
    
    // Clean msaFreq
    if (!nodes.first->msaFreq.empty() || !nodes.second->msaFreq.empty()) {
        nodes.first->msaFreq.clear();
        nodes.second->msaFreq.clear();
    }

}

msa::IntPair msa::alignment_helper::smith_waterman(const std::string &seq1, const std::string &seq2) 
{
    int m = seq1.size();
    int n = seq2.size();
    // Scoring parameters
    const int MATCH = 4;
    const int MISMATCH = -1;
    const int GAP = -2;
    std::vector<std::vector<int>> H(m+1, std::vector<int>(n+1, 0));
    int max_score = 0;
    int max_i = 0, max_j = 0;
    // Fill scoring matrix
    for(int i=1;i<=m;++i){
        for(int j=1;j<=n;++j){
            int match = (seq1[i-1]==seq2[j-1]) ? MATCH : MISMATCH;
            int score_diag = H[i-1][j-1] + match;
            int score_up = H[i-1][j] + GAP;
            int score_left = H[i][j-1] + GAP;
            H[i][j] = std::max({0, score_diag, score_up, score_left});
            if(H[i][j] > max_score){
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    return {max_i, max_j};
}

msa::alnPath msa::alignment_helper::smith_waterman_tb(const std::string &seq1, const std::string &seq2)
{
    int m = seq1.size();
    int n = seq2.size();
    const int MATCH = 4;
    const int MISMATCH = -2;
    const int GAP_OPEN = -5;
    const int GAP_EXTEND = -1;
    const int NSCORE = 1;
    // DP matrices: M = match/mismatch, X = gap in seq1 (deletion), Y = gap in seq2 (insertion)
    std::vector<std::vector<float>> M(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<float>> X(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<float>> Y(m + 1, std::vector<float>(n + 1, 0));
    std::vector<std::vector<int8_t>> tb(m + 1, std::vector<int8_t>(n + 1, 0)); // traceback: 0/1/2
    alnPath aln;
    // Initialization
    M[0][0] = 0;
    for (int i = 1; i <= m; ++i)
    {
        // M[i][0] = gap_open + (i - 1) * gap_extend;
        M[i][0] = 0;
        X[i][0] = M[i][0];
        Y[i][0] = -1e9; // impossible
        tb[i][0] = 2;   // gap in query
    }
    for (int j = 1; j <= n; ++j)
    {
        // M[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND;
        M[0][j] = 0;
        Y[0][j] = M[0][j];
        X[0][j] = -1e9;
        tb[0][j] = 1; // gap in reference
    }
    // Fill matrices
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= n; ++j)
        {
            // Compute scores
            int base_score = ((toupper(seq1[i - 1]) == 'N') || (toupper(seq2[j - 1]) == 'N')) ? NSCORE : 
                             (toupper(seq1[i - 1]) == toupper(seq2[j - 1])) ? MATCH : MISMATCH;
            // Match/mismatch
            M[i][j] = base_score + std::max({M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]});
            // Gap in seq1 (deletion)
            X[i][j] = std::max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXTEND);
            // if (i < m) X[i][j] = std::max(M[i - 1][j] + GAP_OPEN, X[i - 1][j] + GAP_EXTEND);
            // else       X[i][j] = std::max(M[i - 1][j] + 0, X[i - 1][j] + GAP_EXTEND);
            // Gap in seq2 (insertion)
            Y[i][j] = std::max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXTEND);
            // if (j < n) Y[i][j] = std::max(M[i][j - 1] + GAP_OPEN, Y[i][j - 1] + GAP_EXTEND);
            // else       Y[i][j] = std::max(M[i][j - 1] + 0, Y[i][j - 1] + GAP_EXTEND);
            // Choose best path for traceback
            float best = std::max({M[i][j], X[i][j], Y[i][j]});
            if (best == M[i][j])
                tb[i][j] = 0;
            else if (best == X[i][j])
                tb[i][j] = 2; // gap in ref (insert in query)
            else
                tb[i][j] = 1; // gap in qry (delete in ref)
        }
    }
    // Traceback
    aln.clear();
    int i = m, j = n;
    float max_score = std::max({M[m][n], X[m][n], Y[m][n]});
    while (i > 0 || j > 0)
    {
        int8_t dir = tb[i][j];
        aln.push_back(dir);
        if (dir == 0)
        {
            i--;
            j--;
        }
        else if (dir == 1)
        {
            j--;
        }
        else if (dir == 2)
        {
            i--;
        }
    }
    std::reverse(aln.begin(), aln.end());
    return aln;
}
*/