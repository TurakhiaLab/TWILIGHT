#ifndef MSA_HPP

#include <hip/hip_runtime.h>
#include "../msa.hpp"
#endif

#ifndef DEVICE_HPP
#include "device-function.cuh"
#endif

#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>

void msa::progressive::gpu::alignmentKernel_GPU(Tree* T, NodePairVec& alnPairs, SequenceDB* database, Option* option, Params &param) {
    int cpuThres = option->cpuNum * 3;
    if (option->cpuOnly || alnPairs.size() < cpuThres || database->currentTask == 2) cpu::parallelAlignmentCPU(T, alnPairs, database, option, param);
    else                                                                                  parallelAlignmentGPU(T, alnPairs, database, option, param);
}

void msa::progressive::gpu::GPU_pointers::allocateMemory(Option *option, int len, int numBlocks) {
    hipSetDevice(this->gpu_index);
    int profileSize = option->type == 'n' ? 6 : 22;
    int paramSize   = option->type == 'n' ? 5*5+4 : 21*21+4;
    // For Profile pointers
    this->hostPointers = (float**)malloc(this->memBlock * sizeof(float*));
    hipMalloc((void **)&(this->devicePointers), this->memBlock * sizeof(float*));
    // Allocate memory
    int pairPerMemBlock = (numBlocks % this->memBlock == 0) ? numBlocks/this->memBlock : numBlocks/this->memBlock+1;
    for (int i = 0; i < this->memBlock; ++i) {
        this->hostFreq[i]  = (float*)malloc(        profileSize * 2 * len * pairPerMemBlock * sizeof(float));
        hipMalloc((void **)&(this->deviceFreq[i]), profileSize * 2 * len * pairPerMemBlock * sizeof(float));
    }
    this->hostGapOp =     (float *)malloc(2 * len * numBlocks * sizeof(float));
    this->hostGapEx =     (float *)malloc(2 * len * numBlocks * sizeof(float));
    this->hostLen   =   (int32_t *)malloc(2 * numBlocks * sizeof(int32_t));
    this->hostNum   =   (int32_t *)malloc(2 * numBlocks * sizeof(int32_t));
    this->hostAln   =    (int8_t *)malloc(2 * len * numBlocks * sizeof(int8_t));
    this->hostAlnLen =  (int32_t *)malloc(numBlocks * sizeof(int32_t));
    this->hostSeqInfo = (int32_t *)malloc(4 * sizeof(int32_t));
    this->hostParam =     (float *)malloc(paramSize * sizeof(float));
    hipMalloc((void **)&(this->deviceGapOp), 2 * len * numBlocks * sizeof(float));
    hipMalloc((void **)&(this->deviceGapEx), 2 * len * numBlocks * sizeof(float));
    hipMalloc((void **)&(this->deviceLen), 2 * numBlocks * sizeof(int32_t));
    hipMalloc((void **)&(this->deviceNum), 2 * numBlocks * sizeof(int32_t));
    hipMalloc((void **)&(this->deviceAln), 2 * len * numBlocks * sizeof(int8_t));
    hipMalloc((void **)&(this->deviceAlnLen), numBlocks * sizeof(int32_t));
    hipMalloc((void **)&(this->deviceSeqInfo), 4 * sizeof(int32_t));
    hipMalloc((void **)&(this->deviceParam), paramSize * sizeof(float));
    std::string error = hipGetErrorString(hipGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to allocate memory [%s] on GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::GPU_pointers::initializeHostMemory(Option *option, int len, int numBlocks, Params &param)
{
    int profileSize = (option->type == 'n') ? 5 + 1 : 21 + 1;
    int paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    for (int i = 0; i < 2 * len * numBlocks; ++i) this->hostAln[i] = 0;
    for (int i = 0; i < 2 * numBlocks; ++i)       this->hostLen[i] = 0;
    for (int i = 0; i < 2 * numBlocks; ++i)       this->hostNum[i] = 0;
    for (int i = 0; i < 2 * len * numBlocks; ++i) this->hostGapOp[i] = 0;
    for (int i = 0; i < 2 * len * numBlocks; ++i) this->hostGapEx[i] = 0;
    for (int i = 0; i < param.matrixSize; ++i)
        for (int j = 0; j < param.matrixSize; ++j)
            this->hostParam[i * param.matrixSize + j] = param.scoringMatrix[i][j];
    this->hostParam[paramSize - 4] = param.gapOpen;
    this->hostParam[paramSize - 3] = param.gapExtend;
    this->hostParam[paramSize - 2] = param.gapBoundary;
    this->hostParam[paramSize - 1] = param.xdrop;
    // For pointers
    int pairPerMemBlock = (numBlocks % this->memBlock == 0) ? numBlocks/this->memBlock : numBlocks/this->memBlock+1;
    for (int i = 0; i < this->memBlock; ++i) {
        this->hostPointers[i] = this->deviceFreq[i];
        for (int j = 0; j < profileSize * 2 * len * pairPerMemBlock; ++j) this->hostFreq[i][j] = 0.0;
    }
    int idx = 0;
    for (int i = 0; i < numBlocks; ++i) {
        this->hostAlnLen[i] = idx;
        if ((i+1) % pairPerMemBlock == 0) idx += 1;
    }
    return;
}

void msa::progressive::gpu::GPU_pointers::freeMemory()
{
    // free memory
    hipSetDevice(this->gpu_index);
    for (int i = 0; i < this->memBlock; ++i) hipFree(this->deviceFreq[i]);
    hipFree(this->deviceGapOp);
    hipFree(this->deviceGapEx);
    hipFree(this->deviceLen);
    hipFree(this->deviceNum);
    hipFree(this->deviceAln);
    hipFree(this->deviceAlnLen);
    hipFree(this->deviceSeqInfo);
    hipFree(this->deviceParam);
    hipFree(this->devicePointers);
    std::string freeErr = hipGetErrorString(hipGetLastError());
    if (freeErr != "no error") fprintf(stderr, "CUDA Error: Failed to free memory [%s] on GPU %d\n", freeErr.c_str(), this->gpu_index);
    for (int i = 0; i < this->memBlock; ++i) free(this->hostFreq[i]);
    free(this->hostGapOp);
    free(this->hostGapEx);
    free(this->hostLen);
    free(this->hostNum);
    free(this->hostAln);
    free(this->hostAlnLen);
    free(this->hostSeqInfo);
    free(this->hostParam);
    free(this->hostPointers);
    return;
}

void msa::progressive::gpu::GPU_pointers::memcpyHost2Device(Option *option, int len, int alnPairs, int numBlocks)
{
    hipSetDevice(this->gpu_index);
    int paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    int profileSize = (option->type == 'n') ? 5 + 1 : 21 + 1;
    int pairPerMemBlock = (numBlocks % this->memBlock == 0) ? numBlocks/this->memBlock : numBlocks/this->memBlock+1;
    int usedMemBlock = (alnPairs % pairPerMemBlock == 0) ? alnPairs/pairPerMemBlock : alnPairs/pairPerMemBlock+1;
    for (int i = 0; i < usedMemBlock; ++i) 
        hipMemcpy(this->deviceFreq[i], this->hostFreq[i], profileSize * 2 * len * pairPerMemBlock * sizeof(float),   hipMemcpyHostToDevice);
    hipMemcpy(this->devicePointers, this->hostPointers,              this->memBlock * sizeof(float*),   hipMemcpyHostToDevice);
    hipMemcpy(this->deviceGapOp,    this->hostGapOp,              2 * len * alnPairs * sizeof(float),   hipMemcpyHostToDevice);
    hipMemcpy(this->deviceGapEx,    this->hostGapEx,              2 * len * alnPairs * sizeof(float),   hipMemcpyHostToDevice);
    hipMemcpy(this->deviceAln,      this->hostAln,                2 * len * alnPairs * sizeof(int8_t),  hipMemcpyHostToDevice);
    hipMemcpy(this->deviceLen,      this->hostLen,                2       * alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
    hipMemcpy(this->deviceNum,      this->hostNum,                2       * alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
    hipMemcpy(this->deviceAlnLen,   this->hostAlnLen,                       alnPairs * sizeof(int32_t), hipMemcpyHostToDevice);
    hipMemcpy(this->deviceSeqInfo,  this->hostSeqInfo,            4 *                  sizeof(int32_t), hipMemcpyHostToDevice);
    hipMemcpy(this->deviceParam,    this->hostParam,                       paramSize * sizeof(float),   hipMemcpyHostToDevice);
    std::string error = hipGetErrorString(hipGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to copy memory [%s] to GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::GPU_pointers::memcpyDevice2Host(Option *option, int len, int alnPairs)
{
    hipSetDevice(this->gpu_index);
    hipMemcpy(this->hostAln,    this->deviceAln,    2 * len * alnPairs * sizeof(int8_t),  hipMemcpyDeviceToHost);
    hipMemcpy(this->hostAlnLen, this->deviceAlnLen,           alnPairs * sizeof(int32_t), hipMemcpyDeviceToHost);
    std::string error = hipGetErrorString(hipGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to copy memory [%s] from GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::parallelAlignmentGPU(Tree *tree, NodePairVec &nodes, SequenceDB *database, Option *option, Params &param)
{
    auto startT = std::chrono::high_resolution_clock::now();
    int numBlocks = device_function::_BLOCKSIZE;
    int blockSize = device_function::_THREAD_NUM;
    int gpuNum = option->gpuNum;
    int profileSize = param.matrixSize + 1;
    const float gpuMemUsage = 0.7;

    // get maximum sequence/profile length
    int32_t seqLen = 0;
    for (auto n : nodes) seqLen = std::max(seqLen, std::max(n.first->getAlnLen(database->currentTask), n.second->getAlnLen(database->currentTask)));
    auto minMem = std::min_element(option->gpuMem.begin(), option->gpuMem.end());
    if ((static_cast<float>(numBlocks * profileSize * 2 * seqLen) / (*minMem)*gpuMemUsage) > 1.0) {
        float larger = (static_cast<float>(numBlocks * profileSize * 2 * seqLen) / (*minMem)*gpuMemUsage);
        if (larger < 2.0)       numBlocks = 1024;
        else if (larger >= 2.0) numBlocks = 512;
    }
    int memBlock = seqLen / (_MAX_LENGTH_N * ( device_function::_BLOCKSIZE/numBlocks)) + 1;
    int pairPerMemBlock = (numBlocks % memBlock == 0) ? numBlocks/memBlock : numBlocks/memBlock+1;

    // Divide into batches
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size() % numBlocks == 0) roundGPU -= 1;
    if (roundGPU < gpuNum) gpuNum = roundGPU;
    // Adjust Length
    int gpuMemLen = seqLen;

    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::spin_rw_mutex fallbackMutex;
    std::atomic<uint64_t> kernelTime, copyTime;
    kernelTime.store(0);
    copyTime.store(0);

    std::vector<int> fallbackPairs;

    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range) {
    for (int gn = range.begin(); gn < range.end(); ++gn) {
        GPU_pointers* gp = new GPU_pointers(option->gpuIdx[gn], memBlock);
        gp->allocateMemory(option, gpuMemLen, numBlocks);
        while (nowRound < roundGPU) {
            int rn = nowRound.fetch_add(1);
            int alnPairs = (nodes.size() - rn * numBlocks > numBlocks) ? numBlocks : nodes.size() - rn * numBlocks;
            // Initialization
            alnPathVec aln(alnPairs);
            stringPairVec consensus(alnPairs, {"", ""});
            IntPairVec starting(alnPairs, {0,0});
            IntPairVec ending(alnPairs, {0,0});
            std::vector<std::pair<IntPairVec, IntPairVec>> gappyColumns(alnPairs);
            gp->initializeHostMemory(option, gpuMemLen, numBlocks, param);
            // Compute Profile and Gap Penalties
            tbb::this_task_arena::isolate([&] { 
            tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> r) {
            for (int n = r.begin(); n < r.end(); ++n) {
                int32_t nIdx = n + rn*numBlocks;
                int32_t refLen = nodes[nIdx].first->getAlnLen(database->currentTask);
                int32_t qryLen = nodes[nIdx].second->getAlnLen(database->currentTask);
                int32_t refNum = nodes[nIdx].first->getAlnNum(database->currentTask);
                int32_t qryNum = nodes[nIdx].second->getAlnNum(database->currentTask);
                int memBlockIdx = gp->hostAlnLen[n];  
                int32_t offset_f = profileSize * 2 * gpuMemLen * (n%pairPerMemBlock);
                int32_t offset_g =               2 * gpuMemLen * (n);
                IntPair lens = {refLen, qryLen};
                ending[n] = {refLen, qryLen};
                float* rawProfile = new float [profileSize * 2 * gpuMemLen];
                for (int i = 0; i < profileSize * 2 * gpuMemLen; ++i) rawProfile[i] = 0.0;
                alignment_helper::calculateProfile(rawProfile, nodes[nIdx], database, option, gpuMemLen); 
                // Get consensus sequences
                alignment_helper::getConsensus(option, rawProfile,                       consensus[n].first,  lens.first);
                alignment_helper::getConsensus(option, rawProfile+profileSize*gpuMemLen, consensus[n].second, lens.second);
                // alignment_helper::getStartingAndEndingPoints(consensus[n].first, consensus[n].second, starting[n], ending[n]);
                lens.first = (ending[n].first-starting[n].first);
                lens.second = (ending[n].second-starting[n].second);
                // for (int i = 0; i < profileSize * lens.first; ++i) rawProfile[i] = rawProfile[i+starting[n].first*profileSize];
                // for (int i = profileSize * lens.first; i < profileSize * gpuMemLen; ++i) rawProfile[i] = 0.0;
                // for (int i = 0; i < profileSize * lens.second; ++i) rawProfile[profileSize*gpuMemLen+i] = rawProfile[profileSize*gpuMemLen+i+starting[n].second*profileSize];
                // for (int i = profileSize * lens.second; i < profileSize * gpuMemLen; ++i) rawProfile[profileSize*gpuMemLen+i] = 0.0;
                alignment_helper::removeGappyColumns(rawProfile, nodes[nIdx], option, gappyColumns[n], gpuMemLen, lens, database->currentTask);
                // Copy profile to host memory
                for (int i = 0; i < profileSize * gpuMemLen; ++i) gp->hostFreq[memBlockIdx][offset_f+i] =                       (i < profileSize * lens.first)  ? rawProfile[i] : 0.0;
                for (int i = 0; i < profileSize * gpuMemLen; ++i) gp->hostFreq[memBlockIdx][offset_f+profileSize*gpuMemLen+i] = (i < profileSize * lens.second) ? rawProfile[profileSize*gpuMemLen+i] : 0.0;
                std::pair<int32_t, int32_t> offset = std::make_pair(offset_f, offset_g);
                alignment_helper::calculatePSGP(gp->hostFreq[memBlockIdx], gp->hostGapOp, gp->hostGapEx, nodes[nIdx], database, option, gpuMemLen, offset, lens, param);
                gp->hostLen[2*n] = lens.first; gp->hostLen[2*n+1] = lens.second;
                gp->hostNum[2*n] = refNum;     gp->hostNum[2*n+1] = qryNum;
                delete [] rawProfile;
            } 
            }); 
            });
            gp->hostSeqInfo[0] = alnPairs;
            gp->hostSeqInfo[1] = gpuMemLen;
            gp->hostSeqInfo[2] = profileSize;
            gp->hostSeqInfo[3] = pairPerMemBlock;
            auto copyStart = std::chrono::high_resolution_clock::now();
            gp->memcpyHost2Device(option, gpuMemLen, alnPairs, numBlocks);
            auto copyEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds cTime = copyEnd - copyStart;
            int ct = copyTime.fetch_add(cTime.count());
            std::string berr = hipGetErrorString(hipGetLastError());
            if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
            auto kernelStart = std::chrono::high_resolution_clock::now();
            device_function::parallelProfileAlignment<<<numBlocks, blockSize>>>(
                gp->devicePointers,
                gp->deviceAln,
                gp->deviceLen,
                gp->deviceNum,
                gp->deviceAlnLen,
                gp->deviceSeqInfo,
                gp->deviceGapOp,
                gp->deviceGapEx,
                gp->deviceParam
            );
            hipDeviceSynchronize();
            std::string aerr = hipGetErrorString(hipGetLastError());
            if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
            auto kernelEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds kTime = kernelEnd - kernelStart;
            int kt = kernelTime.fetch_add(kTime.count());
            gp->memcpyDevice2Host(option, gpuMemLen, alnPairs);
            tbb::this_task_arena::isolate([&] { 
            tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
            for (int n = range.begin(); n < range.end(); ++n) {
                int32_t nIdx = n + rn*numBlocks;
                int32_t refLen = nodes[nIdx].first->getAlnLen(database->currentTask);
                int32_t qryLen = nodes[nIdx].second->getAlnLen(database->currentTask);
                int32_t refNum = nodes[nIdx].first->getAlnNum(database->currentTask);
                int32_t qryNum = nodes[nIdx].second->getAlnNum(database->currentTask);
                if ((database->currentTask == 0) && (gp->hostAlnLen[n] == -1) && (refLen > 0 && qryLen > 0) && (option->alnMode != 2)) { 
                    // Unable to complete the alignment on GPU
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex);
                        fallbackPairs.push_back(nIdx);
                    }
                    aln[n].clear();
                }
                else {
                    if ((database->currentTask == 0 && option->alnMode != 2) && 
                        ((refNum == 1 && database->id_map[nodes[nIdx].first->seqsIncluded[0]]->lowQuality) || 
                         (qryNum == 1 && database->id_map[nodes[nIdx].second->seqsIncluded[0]]->lowQuality))) {
                        aln[n].clear();
                        {
                            tbb::spin_rw_mutex::scoped_lock lock(fallbackMutex, true);
                            fallbackPairs.push_back(nIdx);
                        }
                    }    
                    else {
                        alnPath aln_wo_gc;
                        int alnRef = 0, alnQry = 0;
                        if (gp->hostAlnLen[n] > 0) for (int j = 0; j < gp->hostAlnLen[n]; ++j) aln_wo_gc.push_back(gp->hostAln[n*2*gpuMemLen+j]);
                        else if (refLen == 0 || qryLen == 0) {
                            if (refLen == 0) for (int j = 0; j < qryLen; ++j) aln_wo_gc.push_back(1);
                            if (qryLen == 0) for (int j = 0; j < refLen; ++j) aln_wo_gc.push_back(2);
                        }
                        // Alignment path without gappy columns
                        for (auto a: aln_wo_gc) {
                            if (a == 0) {alnRef += 1; alnQry += 1;}
                            if (a == 1) {alnQry += 1;}
                            if (a == 2) {alnRef += 1;}
                        }
                        if (alnRef != gp->hostLen[2*n])   std::cout << "R: Pre " << nodes[nIdx].first->identifier << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
                        if (alnQry != gp->hostLen[2*n+1]) std::cout << "Q: Pre " << nodes[nIdx].second->identifier << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";

                        // Add gappy columns back
                        IntPair debugIdx = std::make_pair(alnRef,alnQry);
                        msa::alignment_helper::addGappyColumnsBack(aln_wo_gc, aln[n], gappyColumns[n], param, {alnRef,alnQry}, consensus[n]);
                        alnRef = 0, alnQry = 0;
                        for (auto a: aln[n]) {
                            if (a == 0) {alnRef += 1; alnQry += 1;}
                            if (a == 1) {alnQry += 1;}
                            if (a == 2) {alnRef += 1;}
                        }
                        if (alnRef != refLen) {
                            std::cout << "R: Post " << nodes[nIdx].first->identifier <<  '\t' << refNum << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
                            for (auto gp: gappyColumns[n].first) std::cout << '(' << gp.first << ',' << gp.second << ')' << ',';
                            std::cout << '\n';
                        }
                        if (alnQry != qryLen) {
                            std::cout << "Q: Post " << nodes[nIdx].second->identifier << '\t' << qryNum << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";
                            for (auto gp: gappyColumns[n].second) std::cout << '(' << gp.first << ',' << gp.second << ')' << ',';
                            std::cout << '\n';
                        }

                    }
                }
            } 
            });
            });
            tbb::this_task_arena::isolate([&] { 
            tbb::parallel_for(tbb::blocked_range<int>(0, alnPairs), [&](tbb::blocked_range<int> range) {
            for (int n = range.begin(); n < range.end(); ++n) {
                if (!aln[n].empty()) {
                    int32_t nIdx = n + rn*numBlocks;
                    std::map<int,std::string> startAln, endAln;
                    float refWeight = nodes[nIdx].first->alnWeight, qryWeight = nodes[nIdx].second->alnWeight;
                    // alignment_helper::alignEnds(nodes[nIdx], database, option->type, starting[n], ending[n], startAln, endAln);
                    alignment_helper::updateAlignment(nodes[nIdx], database, aln[n], starting[n], ending[n], startAln, endAln);
                    alignment_helper::updateFrequency(nodes[nIdx], database, aln[n], {refWeight, qryWeight});
                }
            } 
            }); 
            });
        }
        gp->freeMemory();
        delete gp;
    }
    });

    // free memory
    if (fallbackPairs.empty()) return;
    // fallback to CPU
    alignment_helper::fallback2cpu(fallbackPairs, nodes, database, option);
    return;
}
