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
                float w = database->id_map[sIdx]->weight / refWeight * refNum;;
                int storage = database->id_map[sIdx]->storage;
                tbb::this_task_arena::isolate([&] { 
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int t = r.begin(); t < r.end(); ++t) {
                    int letterIndex = letterIdx(option->type, toupper(database->id_map[sIdx]->alnStorage[storage][t]));
                    profile[profileSize*t+letterIndex] += 1.0 * w;  
                } 
                }); 
                });
            }
            if (storeFreq) {
                nodes.first->msaFreq.resize(refLen, std::vector<float>(6, 0.0));
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
                float w = database->id_map[sIdx]->weight / qryWeight * qryNum;;
                int storage = database->id_map[sIdx]->storage;
                tbb::this_task_arena::isolate([&] { 
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int t = r.begin(); t < r.end(); ++t) {
                    int letterIndex = letterIdx(option->type, toupper(database->id_map[sIdx]->alnStorage[storage][t]));
                    profile[profileSize*(memLen+t)+letterIndex] += 1.0 * w;
                } 
                }); 
                });
            }
            if (storeFreq) {
                nodes.second->msaFreq.resize(qryLen, std::vector<float>(6, 0.0));
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
        M[i][0] = gap_open + (i - 1) * gap_extend;
        X[i][0] = M[i][0];
        Y[i][0] = -1e9; // impossible
        tb[i][0] = 2;   // gap in query
    }
    for (int j = 1; j <= n; ++j)
    {
        M[0][j] = gap_open + (j - 1) * gap_extend;
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
            int idx_1 = letterIdx(type, seq1[i - 1]), idx_2 = letterIdx(type, seq2[j - 1]);
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
    for (alnIdx = 0; alnIdx < aln_before.size(); ++alnIdx) {
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
    return;
}

void msa::alignment_helper::updateAlignment(NodePair &nodes, SequenceDB *database, alnPath &aln)
{
    tbb::this_task_arena::isolate([&] { 
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.first->seqsIncluded.size()), [&](tbb::blocked_range<int> range) {
    for (int idx = range.begin(); idx < range.end(); ++idx) {
        int sIdx = nodes.first->seqsIncluded[idx];
        bool updateSeq = (database->currentTask != 2) || 
                         (database->currentTask == 2 && database->id_map.find(sIdx) != database->id_map.end());
        if (updateSeq) {
            database->id_map[sIdx]->memCheck(aln.size());
            int orgIdx = 0;
            int storeFrom = database->id_map[sIdx]->storage;
            int storeTo = 1 - storeFrom;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 2) {
                    database->id_map[sIdx]->alnStorage[storeTo][k] = database->id_map[sIdx]->alnStorage[storeFrom][orgIdx];
                    orgIdx++;
                }
                else {
                    database->id_map[sIdx]->alnStorage[storeTo][k] = '-';
                }
            }
            database->id_map[sIdx]->len = aln.size();
            database->id_map[sIdx]->changeStorage();
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
        bool updateSeq = (database->currentTask != 2) || 
                         (database->currentTask == 2 && database->id_map.find(sIdx) != database->id_map.end());
        if (updateSeq) {
            database->id_map[sIdx]->memCheck(aln.size());
            int orgIdx = 0;
            int storeFrom = database->id_map[sIdx]->storage;
            int storeTo = 1 - storeFrom;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 1) {
                    database->id_map[sIdx]->alnStorage[storeTo][k] = database->id_map[sIdx]->alnStorage[storeFrom][orgIdx];
                    orgIdx++;
                }
                else {
                    database->id_map[sIdx]->alnStorage[storeTo][k] = '-';
                }
            }
            database->id_map[sIdx]->len = aln.size();
            database->id_map[sIdx]->changeStorage();
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
    for (auto idx: nodes.second->seqsIncluded) nodes.first->seqsIncluded.push_back(idx);
    nodes.first->alnNum += nodes.second->alnNum;
    nodes.second->seqsIncluded.clear();
    nodes.first->alnLen = aln.size();
    nodes.first->alnWeight += nodes.second->alnWeight;
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
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int32_t refNum = nodes[nIdx].first->seqsIncluded.size();
        int32_t qryNum = nodes[nIdx].second->seqsIncluded.size();
        bool lowQ_r = (refNum == 1 && database->id_map[nodes[nIdx].first->seqsIncluded[0]]->lowQuality);
        bool lowQ_q = (qryNum == 1 && database->id_map[nodes[nIdx].second->seqsIncluded[0]]->lowQuality);
        if ((refNum < qryNum) || lowQ_r) {
            bool fallback = (!filtering) || (!lowQ_r);
            if (fallback) {
                database->fallback_nodes.push_back(nodes[nIdx].second);
                if (lowQ_r) database->id_map[nodes[nIdx].first->seqsIncluded[0]]->lowQuality = false;
            }
            // swap reference and query
            int32_t refLen = nodes[nIdx].first->alnLen;
            int32_t qryLen = nodes[nIdx].second->alnLen;
            nodes[nIdx].second->alnLen = refLen;
            nodes[nIdx].first->alnLen = qryLen;
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
                if (lowQ_q) database->id_map[nodes[nIdx].second->seqsIncluded[0]]->lowQuality = false;
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