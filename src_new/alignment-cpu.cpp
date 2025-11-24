#ifndef MSA_HPP
#include "msa.hpp"
#endif

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>
#include <chrono>

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
    
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){ 
    for (int nIdx = range.begin(); nIdx < range.end(); ++nIdx) {
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
        IntPair lens = {refLen, qryLen}, starting = {0,0}, ending = {refLen,qryLen};

        const int refineTh = 1000;
        // bool refine = (refNum < refineTh && qryNum < refineTh && (refNum+qryNum > refineTh));
        bool refine = false;
        alignment_helper::calculateProfile(hostFreq, nodes[nIdx], database, option, memLen);
        alignment_helper::getConsensus(option, hostFreq,                    consensus.first,  refLen);
        alignment_helper::getConsensus(option, hostFreq+profileSize*memLen, consensus.second, qryLen);
        if (refine) {
            alignment_helper::getStartingAndEndingPoints(consensus.first, consensus.second, starting, ending);
            lens.first = (ending.first-starting.first);
            lens.second = (ending.second-starting.second);
            std::cout << refNum << ',' << qryNum << '\n';
            std::cout << starting.first << ',' << starting.second << ',' << ending.first << ',' << ending.second << '\n';
            std::cout << refLen << ',' << qryLen << ',' << lens.first << ',' << lens.second << '\n';
            for (int i = 0; i < profileSize * lens.first; ++i) hostFreq[i] = hostFreq[profileSize*starting.first+i];
            for (int i = profileSize * lens.first; i < profileSize * memLen; ++i) hostFreq[i] = 0.0;
            for (int i = 0; i < profileSize * lens.second; ++i) hostFreq[profileSize*memLen+i] = hostFreq[profileSize*(memLen+starting.second)+i];
            for (int i = profileSize * lens.second; i < profileSize * memLen; ++i) hostFreq[profileSize*memLen+i] = 0.0;
        }
        alignment_helper::removeGappyColumns(hostFreq, nodes[nIdx], option, gappyColumns, memLen, lens, database->currentTask);
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
        bool lowQ_r = (refNum == 1) && database->id_map[nodes[nIdx].first->seqsIncluded[0]]->lowQuality;
        bool lowQ_q = (qryNum == 1) && database->id_map[nodes[nIdx].second->seqsIncluded[0]]->lowQuality;
        if (option->noFilter || (!lowQ_r && !lowQ_q)) {
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
                    talco_params->updateXDrop(static_cast<int32_t>(talco_params->xdrop * 1.2));
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
            if (alnRef != ending.first-starting.first)   std::cout << "R: Post " << nodes[nIdx].first->identifier << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
            if (alnQry != ending.second-starting.second) std::cout << "Q: Post " << nodes[nIdx].second->identifier << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";
            std::map<int,std::string> startAln, endAln;
            if (refine) database->debug();
            if (refine) alignment_helper::alignEnds(nodes[nIdx], database, option->type, starting, ending, startAln, endAln);
            alignment_helper::updateAlignment(nodes[nIdx], database, aln_w_gc, starting, ending, startAln, endAln);
            alignment_helper::updateFrequency(nodes[nIdx], database, aln_w_gc, {refWeight, qryWeight});
            if (refine) database->debug();
        }
    }    
    });
    // std::cerr << "Alignment Kernel Time: " << talcoTime / 1000 << " us\n"
    //           << "Preprocessing Time:    " << preprocessTime / 1000 << " us\n";
    if (fallbackPairs.empty()) return;
    alignment_helper::fallback2cpu(fallbackPairs, nodes, database, option);
    return;
}
