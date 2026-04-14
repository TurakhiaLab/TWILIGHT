#ifndef MSA_HPP
#include "msa.hpp"
#endif

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#include "consistency_v2.hpp"

#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>
#include <chrono>
#include <algorithm>

namespace {

void removeColumns(msa::ColumnProvenance& provenance, const msa::IntPairVec& removedColumns)
{
    if (removedColumns.empty()) return;
    msa::ColumnProvenance filtered;
    filtered.reserve(provenance.size());
    int removeIdx = 0;
    int nextRemoveStart = removedColumns[removeIdx].first;
    int nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
    for (int col = 0; col < static_cast<int>(provenance.size()); ++col) {
        while (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveEnd) {
            ++removeIdx;
            if (removeIdx < static_cast<int>(removedColumns.size())) {
                nextRemoveStart = removedColumns[removeIdx].first;
                nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
            }
        }
        const bool removed = (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveStart && col < nextRemoveEnd);
        if (!removed) filtered.push_back(std::move(provenance[col]));
    }
    provenance = std::move(filtered);
}

std::vector<std::vector<float>> buildConsistencyTable(
    const msa::ColumnProvenance& refProvenance,
    const msa::ColumnProvenance& qryProvenance,
    msa::accurate::SubtreeAccurateState& accurateState)
{
    int refLen = refProvenance.size();
    int qryLen = qryProvenance.size();
    std::vector<std::vector<float>> consistencyTable(refLen, std::vector<float>(qryLen, 0.0f));

    int totalSeq = accurateState.directLib.getSequence();
    int maxLen = accurateState.directLib.length;

    struct QryResInfo {
        int qryCol;
        int seqB;
        int posB;
        float weight;
        float S_B;
    };

    std::vector<QryResInfo> qryResList;
    std::vector<int> qryMap(totalSeq * maxLen, -1);
    std::vector<float> qryColWeightSum(qryLen, 0.0f);
    std::vector<float> refColWeightSum(refLen, 0.0f);

    for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
        for (const auto& qryRes : qryProvenance[qryCol]) {
            int seqB = qryRes.seqId;
            int posB = qryRes.residueIndex;
            int id = qryResList.size();
            float S_B = accurateState.directLib.residueSupport[seqB * maxLen + posB];
            qryResList.push_back({qryCol, seqB, posB, qryRes.weight, S_B});
            qryMap[seqB * maxLen + posB] = id;
            qryColWeightSum[qryCol] += qryRes.weight;
        }
    }

    for (int refCol = 0; refCol < refLen; ++refCol) {
        for (const auto& refRes : refProvenance[refCol]) {
            refColWeightSum[refCol] += refRes.weight;
        }
    }

    int qryResCount = qryResList.size();

    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, refLen), [&](const tbb::blocked_range<std::size_t>& range) {
        std::vector<float> num_B(qryResCount, 0.0f);
        std::vector<int> active_B;
        active_B.reserve(totalSeq);

        for (std::size_t refCol = range.begin(); refCol < range.end(); ++refCol) {
            if (refProvenance[refCol].empty()) continue;

            std::vector<float> col_scores(qryLen, 0.0f);
            for (const auto& refRes : refProvenance[refCol]) {
                int seqA = refRes.seqId;
                int posA = refRes.residueIndex;
                float wA = refRes.weight;
                float S_A = accurateState.directLib.residueSupport[seqA * maxLen + posA];
                if (S_A <= 0.0f) continue;

                for (int seqB = 0; seqB < totalSeq; ++seqB) {
                    if (seqB == seqA) continue;
                    int posB = accurateState.directLib.getAlignedPos(seqA, seqB, posA);
                    if (posB != -1) {
                        int qID = qryMap[seqB * maxLen + posB];
                        if (qID != -1) {
                            if (num_B[qID] == 0.0f) active_B.push_back(qID);
                            num_B[qID] += accurateState.directLib.getWeight(seqA, seqB);
                        }
                    }
                }

                for (int seqC = 0; seqC < totalSeq; ++seqC) {
                    if (seqC == seqA) continue;
                    int posC = accurateState.directLib.getAlignedPos(seqA, seqC, posA);
                    if (posC == -1) continue;
                    float w_AC = accurateState.directLib.getWeight(seqA, seqC);
                    for (int seqB = 0; seqB < totalSeq; ++seqB) {
                        if (seqB == seqA || seqB == seqC) continue;
                        int posB = accurateState.directLib.getAlignedPos(seqC, seqB, posC);
                        if (posB != -1) {
                            int qID = qryMap[seqB * maxLen + posB];
                            if (qID != -1) {
                                float w_BC = accurateState.directLib.getWeight(seqC, seqB);
                                if (num_B[qID] == 0.0f) active_B.push_back(qID);
                                num_B[qID] += std::min(w_AC, w_BC);
                            }
                        }
                    }
                }
                
                for (int qID : active_B) {
                    float num = num_B[qID];
                    num_B[qID] = 0.0f;
                    const auto& qInfo = qryResList[qID];
                    float S_B = qInfo.S_B;
                    float tcs = (S_A + S_B > 0.0f) ? std::clamp((2.0f * num) / (S_A + S_B), 0.0f, 1.0f) : 0.0f;
                    col_scores[qInfo.qryCol] += wA * qInfo.weight * tcs;
                }
                active_B.clear();
            }
            
            float denom_ref = refColWeightSum[refCol];
            if (denom_ref <= 0.0f) continue;
            for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
                float denom = denom_ref * qryColWeightSum[qryCol];
                if (denom > 0.0f) {
                    consistencyTable[refCol][qryCol] = col_scores[qryCol] / denom;
                }
            }
        }
    });

    return consistencyTable;
}

} // namespace

void msa::progressive::cpu::allocateMemory_and_Initialize(float*& freq, float*& gapOp, float*& gapEx, int memLen, int profileSize)
{
    freq =  new float [profileSize * 2 * memLen];
    gapOp = new float [              2 * memLen];
    gapEx = new float [              2 * memLen];
    for (int n = 0; n < profileSize * 2 * memLen; ++n) freq[n] = 0;
    for (int n = 0; n <               2 * memLen; ++n) gapOp[n] = 0;
    for (int n = 0; n <               2 * memLen; ++n) gapEx[n] = 0;
    return;
}

void msa::progressive::cpu::freeMemory(float*& freq, float*& gapOp, float*& gapEx)
{
    delete [] freq;
    delete [] gapOp;
    delete [] gapEx;
    return;
}

void msa::progressive::cpu::alignmentKernel_CPU(Tree* T, NodePairVec& alnPairs, SequenceDB* database, Option* option, Params &param) {
    parallelAlignmentCPU(T, alnPairs, database, option, param);
}

void msa::progressive::cpu::parallelAlignmentCPU(Tree *tree, NodePairVec &nodes, SequenceDB *database, Option *option, Params &param)
{
    int profileSize = param.matrixSize + 1;
    tbb::spin_rw_mutex fallbackMutex;
    std::vector<int> fallbackPairs;

    std::atomic<uint64_t> preprocessTime, talcoTime;
    preprocessTime.store(0);
    talcoTime.store(0);

    std::cout << "Total Pairs: " << nodes.size() << '\n';
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
    // for (int nIdx = 0; nIdx <  nodes.size(); ++nIdx) {
        // allocate memory
        auto preStart = std::chrono::high_resolution_clock::now();
        int32_t refLen = nodes[nIdx].first->getAlnLen(database->currentTask);
        int32_t qryLen = nodes[nIdx].second->getAlnLen(database->currentTask);
        int32_t refNum = nodes[nIdx].first->getAlnNum(database->currentTask);
        int32_t qryNum = nodes[nIdx].second->getAlnNum(database->currentTask);
        int32_t memLen = std::max(refLen,qryLen);
        float* hostFreq, * hostGapOp, * hostGapEx; 
        allocateMemory_and_Initialize(hostFreq, hostGapOp, hostGapEx, memLen, profileSize);
        std::pair<IntPairVec, IntPairVec> gappyColumns;
        stringPair consensus({"", ""});
        IntPair lens = {refLen, qryLen};
        // -------
        ColumnProvenance refProvenance, qryProvenance;
        std::vector<std::vector<float>> consistencyTable;
        // -------
        

        alignment_helper::calculateProfile(hostFreq, nodes[nIdx], database, option, memLen);
        // -------
        if (option->accurate && database->currentTask == 0) {
            alignment_helper::extractColumnProvenance(nodes[nIdx].first, database, refProvenance);
            alignment_helper::extractColumnProvenance(nodes[nIdx].second, database, qryProvenance);
            if (option->printDetail) {
                size_t refResidues = 0, qryResidues = 0;
                for (const auto& column : refProvenance) refResidues += column.size();
                for (const auto& column : qryProvenance) qryResidues += column.size();
                // std::cerr << "Accurate mode provenance for pair " << nIdx
                //           << ": ref columns=" << refProvenance.size() << ", residues=" << refResidues
                //           << "; qry columns=" << qryProvenance.size() << ", residues=" << qryResidues << '\n';
            }
        }

        // -------
        alignment_helper::getConsensus(option, hostFreq,                    consensus.first,  refLen);
        alignment_helper::getConsensus(option, hostFreq+profileSize*memLen, consensus.second, qryLen);
        alignment_helper::removeGappyColumns(hostFreq, nodes[nIdx], option, gappyColumns, memLen, lens, database->currentTask);
        // -------
        if (option->accurate && database->currentTask == 0 && database->accurateState) {
            
            removeColumns(refProvenance, gappyColumns.first);
            removeColumns(qryProvenance, gappyColumns.second);
            auto consistStart = std::chrono::high_resolution_clock::now();
            consistencyTable = buildConsistencyTable(refProvenance, qryProvenance, *database->accurateState);
            if (option->printDetail) {
                std::size_t nonZeroEntries = 0;
                for (const auto& row : consistencyTable) {
                    for (float score : row) {
                        if (score > 0.0f) ++nonZeroEntries;
                    }
                }
                auto consistEnd = std::chrono::high_resolution_clock::now();
                std::chrono::nanoseconds consistTime = consistEnd - consistStart;

                std::cerr << "Accurate mode consistency table for pair " << nIdx
                          << ": " << consistencyTable.size() << "x"
                          << (consistencyTable.empty() ? 0 : consistencyTable.front().size())
                          << ", non-zero entries=" << nonZeroEntries << " Time: " << consistTime.count() / 1000 << " us\n";
            }
            
        }
        // -------
        alignment_helper::calculatePSGP(hostFreq, hostGapOp, hostGapEx, nodes[nIdx], database, option, memLen, {0,0}, lens, param);
        auto preEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds pTime = preEnd - preStart;
        int pt = preprocessTime.fetch_add(pTime.count());
        // Start alignment
        std::vector<int8_t> aln_wo_gc;
        Profile freqRef (lens.first, std::vector<float>(profileSize, 0.0));
        Profile freqQry (lens.second, std::vector<float>(profileSize, 0.0));
        Profile gapOp (2), gapEx (2);    
        for (int s = 0; s < lens.first; s++) for (int t = 0; t < profileSize; ++t)  freqRef[s][t] = hostFreq[profileSize*s+t];
        for (int s = 0; s < lens.second; s++) for (int t = 0; t < profileSize; ++t) freqQry[s][t] = hostFreq[profileSize*(memLen+s)+t];     
        for (int r = 0; r < lens.first; ++r) {
            gapOp[0].push_back(hostGapOp[r]);
            gapEx[0].push_back(hostGapEx[r]);
        }
        for (int q = 0; q < lens.second; ++q) {
            gapOp[1].push_back(hostGapOp[memLen+q]);
            gapEx[1].push_back(hostGapEx[memLen+q]);
        }           
        freeMemory(hostFreq, hostGapOp, hostGapEx);
        std::pair<float, float> num = std::make_pair(static_cast<float>(refNum), static_cast<float>(qryNum));
        Talco_xdrop::Params* talco_params = new Talco_xdrop::Params(param);
        if (database->currentTask == 1 || database->currentTask == 2 || refNum > 10000 || qryNum > 10000) talco_params->gapCharScore = 0;
        if (refLen == 0) for (int j = 0; j < qryLen; ++j) aln_wo_gc.push_back(1);
        if (qryLen == 0) for (int j = 0; j < refLen; ++j) aln_wo_gc.push_back(2);
        bool lowQ_r = (option->alnMode == MERGE_MSA) ? false : ((refNum > 1) ? false : database->sequences[nodes[nIdx].first->seqsIncluded[0]]->lowQuality);
        bool lowQ_q = (option->alnMode == MERGE_MSA) ? false : ((qryNum > 1) ? false : database->sequences[nodes[nIdx].second->seqsIncluded[0]]->lowQuality);
       if (!lowQ_r && !lowQ_q) {
            auto talcoStart = std::chrono::high_resolution_clock::now();
            while (aln_wo_gc.empty()) {
                int16_t errorType = 0;
                aln_wo_gc.clear();
                Talco_xdrop::Align_freq (
                    talco_params,
                    freqRef,
                    freqQry,
                    gapOp,
                    gapEx,
                    num,
                    // -------
                    (option->accurate && database->currentTask == 0 && database->accurateState && !consistencyTable.empty()) ? &consistencyTable : nullptr,
                    (option->accurate && database->currentTask == 0 && database->accurateState) ? option->consistencyWeight : 0.0f,
                    // -------
                    aln_wo_gc,
                    errorType
                );
                if (database->currentTask == 0 && errorType != 0) {
                    aln_wo_gc.clear();
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                        fallbackPairs.push_back(nIdx);
                    }
                    break;
                }
                if (errorType == 2) {
                    if (option->printDetail) std::cout << "Updated anti-diagonal limit on No. " << nIdx << '\n';
                    talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->fLen * 1.2) << 1, std::min(lens.first, lens.second)));
                }
                else if (errorType == 3) {
                    std::cout << "There might be some bugs in the code!\n";
                    exit(1);
                }
                else if (errorType == 1) {
                    talco_params->updateXDrop(static_cast<int32_t>(talco_params->xdrop * 2));
                    talco_params->updateFLen(std::min(static_cast<int32_t>(talco_params->xdrop * 4) << 1, std::min(lens.first, lens.second)));
                    if (option->printDetail) std::cout << "Updated x-drop value on No. " << nIdx << "\tNew Xdrop: " << talco_params->xdrop << '\n';

                }
            }
            auto talcoEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds tTime = talcoEnd - talcoStart;
            int ct = talcoTime.fetch_add(tTime.count());
        }
        delete talco_params;
        if (database->currentTask == 0 && (refNum == 1 || qryNum == 1)) {
            if (lowQ_r || lowQ_q) {
                aln_wo_gc.clear();
                {
                    tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                    fallbackPairs.push_back(nIdx);
                }
            }
        }
        
        if (!aln_wo_gc.empty()) {
            alnPath aln_w_gc;
            int alnRef = 0, alnQry = 0;
            // Alignment path without gappy columns
            for (auto a: aln_wo_gc) {
                if (a == 0) {alnRef += 1; alnQry += 1;}
                if (a == 1) {alnQry += 1;}
                if (a == 2) {alnRef += 1;}
            }
            // Add gappy columns back
            alignment_helper::addGappyColumnsBack(aln_wo_gc, aln_w_gc, gappyColumns, param, {alnRef,alnQry}, consensus);
            alnRef = 0, alnQry = 0;
            for (auto a: aln_w_gc) {
                if (a == 0) {alnRef += 1; alnQry += 1;}
                if (a == 1) {alnQry += 1;}
                if (a == 2) {alnRef += 1;}
            }
            float refWeight = nodes[nIdx].first->alnWeight, qryWeight = nodes[nIdx].second->alnWeight;
            if (alnRef != refLen) std::cout << "R: Post " << nodes[nIdx].first->identifier << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
            if (alnQry != qryLen) std::cout << "Q: Post " << nodes[nIdx].second->identifier << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";
            std::map<int,std::string> startAln, endAln;
            
            if (option->alnMode != PLACE_WO_TREE) {
                alignment_helper::updateFrequency(nodes[nIdx], database, aln_w_gc, {refWeight, qryWeight});
                alignment_helper::updateAlignment(nodes[nIdx], database, option, aln_w_gc);
            }
            else {
                database->subtreeAln[nodes[nIdx].second->seqsIncluded[0]] = aln_w_gc;
            }
        }
        // std::cerr << "Alignment Kernel Time: " << talcoTime / 1000 << " us\n"
        //           << "Preprocessing Time:    " << preprocessTime / 1000 << " us\n";
    }    
    });
    // std::cerr << "Alignment Kernel Time: " << talcoTime / 1000 << " us\n"
    //           << "Preprocessing Time:    " << preprocessTime / 1000 << " us\n";
    if (fallbackPairs.empty()) return;
    alignment_helper::fallback2cpu(fallbackPairs, nodes, database, option);
    return;
}
