#ifndef MSA_HPP
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
    cudaSetDevice(this->gpu_index);
    int profileSize = option->type == 'n' ? 6 : 22;
    int paramSize   = option->type == 'n' ? 5*5+4 : 21*21+4;
    std::string error;
    // Pointers (profile, gapOp, gapEx, Aln)
    this->hostFreqPointers  = (float**)malloc(this->memBlock * sizeof(float*));
    this->hostGapOpPointers = (float**)malloc(this->memBlock * sizeof(float*));
    this->hostGapExPointers = (float**)malloc(this->memBlock * sizeof(float*));
    this->hostAlnPointers   = (int8_t**)malloc(this->memBlock * sizeof(int8_t*));

    cudaMalloc((void **)&(this->deviceFreqPointers),  this->memBlock * sizeof(float*));
    cudaMalloc((void **)&(this->deviceGapOpPointers), this->memBlock * sizeof(float*));
    cudaMalloc((void **)&(this->deviceGapExPointers), this->memBlock * sizeof(float*));
    cudaMalloc((void **)&(this->deviceAlnPointers),   this->memBlock * sizeof(int8_t*));

    // Allocate memory
    int pairPerMemBlock = numBlocks/this->memBlock;
    for (int i = 0; i < this->memBlock; ++i) {
        this->hostFreq[i]  = (float*)malloc(        profileSize * 2 * len * pairPerMemBlock * sizeof(float));
        cudaMalloc((void **)&(this->deviceFreq[i]), profileSize * 2 * len * pairPerMemBlock * sizeof(float));
        error = cudaGetErrorString(cudaGetLastError());
        if (error != "no error") fprintf(stderr, "CUDA Error: Failed to allocate memory Freq [%s] on GPU %d\n", error.c_str(), this->gpu_index);
        this->hostGapOp[i] = (float*)malloc(                      2 * len * pairPerMemBlock * sizeof(float));
        cudaMalloc((void **)&(this->deviceGapOp[i]),              2 * len * pairPerMemBlock * sizeof(float));
        error = cudaGetErrorString(cudaGetLastError());
        if (error != "no error") fprintf(stderr, "CUDA Error: Failed to allocate memory GapOp [%s] on GPU %d\n", error.c_str(), this->gpu_index);
    
        this->hostGapEx[i] = (float*)malloc(                      2 * len * pairPerMemBlock * sizeof(float));
        cudaMalloc((void **)&(this->deviceGapEx[i]),              2 * len * pairPerMemBlock * sizeof(float));
        error = cudaGetErrorString(cudaGetLastError());
        if (error != "no error") fprintf(stderr, "CUDA Error: Failed to allocate memory GapEx [%s] on GPU %d\n", error.c_str(), this->gpu_index);
    
        this->hostAln[i]   = (int8_t*)malloc(                     2 * len * pairPerMemBlock * sizeof(int8_t));
        cudaMalloc((void **)&(this->deviceAln[i]),                2 * len * pairPerMemBlock * sizeof(int8_t));
    }
    this->hostLen   =   (int32_t *)malloc(2 * numBlocks * sizeof(int32_t));
    this->hostNum   =   (int32_t *)malloc(2 * numBlocks * sizeof(int32_t));
    this->hostAlnLen =  (int32_t *)malloc(    numBlocks * sizeof(int32_t));
    this->hostSeqInfo = (int32_t *)malloc(4 * sizeof(int32_t));
    this->hostParam =     (float *)malloc(paramSize * sizeof(float));
    cudaMalloc((void **)&(this->deviceLen), 2 * numBlocks * sizeof(int32_t));
    cudaMalloc((void **)&(this->deviceNum), 2 * numBlocks * sizeof(int32_t));
    cudaMalloc((void **)&(this->deviceAlnLen),  numBlocks * sizeof(int32_t));
    cudaMalloc((void **)&(this->deviceSeqInfo), 4 * sizeof(int32_t));
    cudaMalloc((void **)&(this->deviceParam), paramSize * sizeof(float));
    error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to allocate memory [%s] on GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::GPU_pointers::initializeHostMemory(Option *option, int len, int numBlocks, Params &param)
{
    int profileSize = param.matrixSize+1;
    int paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    for (int i = 0; i < 2 * numBlocks; ++i)       this->hostLen[i] = 0;
    for (int i = 0; i < 2 * numBlocks; ++i)       this->hostNum[i] = 0;
    for (int i = 0; i < param.matrixSize; ++i)
        for (int j = 0; j < param.matrixSize; ++j)
            this->hostParam[i * param.matrixSize + j] = param.scoringMatrix[i][j];
    this->hostParam[paramSize - 4] = param.gapOpen;
    this->hostParam[paramSize - 3] = param.gapExtend;
    this->hostParam[paramSize - 2] = param.gapBoundary;
    this->hostParam[paramSize - 1] = param.xdrop;
    // For pointers
    int pairPerMemBlock = numBlocks/(this->memBlock);
    for (int i = 0; i < this->memBlock; ++i) {
        this->hostFreqPointers[i]  = this->deviceFreq[i];
        this->hostGapOpPointers[i] = this->deviceGapOp[i];
        this->hostGapExPointers[i] = this->deviceGapEx[i];
        this->hostAlnPointers[i]   = this->deviceAln[i];
        for (int j = 0; j < profileSize * 2 * len * pairPerMemBlock; ++j) this->hostFreq[i][j] = 0.0;
        for (int j = 0; j < 2 * len * pairPerMemBlock; ++j) {
            this->hostGapOp[i][j] = 0;
            this->hostGapEx[i][j] = 0;
            this->hostAln[i][j] = 0;
        }
    }
    for (int i = 0; i < numBlocks; ++i) this->hostAlnLen[i] = (i/pairPerMemBlock);
    return;
}

void msa::progressive::gpu::GPU_pointers::freeMemory()
{
    // free memory
    cudaSetDevice(this->gpu_index);
    for (int i = 0; i < this->memBlock; ++i) {
        cudaFree(this->deviceFreq[i]);
        cudaFree(this->deviceGapOp[i]);
        cudaFree(this->deviceGapEx[i]);
        cudaFree(this->deviceAln[i]);

    }
    cudaFree(this->deviceLen);
    cudaFree(this->deviceNum);
    cudaFree(this->deviceAlnLen);
    cudaFree(this->deviceSeqInfo);
    cudaFree(this->deviceParam);
    cudaFree(this->deviceFreqPointers);
    cudaFree(this->deviceGapOpPointers);
    cudaFree(this->deviceGapExPointers);
    cudaFree(this->deviceAlnPointers);
    std::string freeErr = cudaGetErrorString(cudaGetLastError());
    if (freeErr != "no error") fprintf(stderr, "CUDA Error: Failed to free memory [%s] on GPU %d\n", freeErr.c_str(), this->gpu_index);
    for (int i = 0; i < this->memBlock; ++i) {
        free(this->hostFreq[i]);
        free(this->hostGapOp[i]);
        free(this->hostGapEx[i]);
        free(this->hostAln[i]);
    }
    free(this->hostLen);
    free(this->hostNum);
    free(this->hostAlnLen);
    free(this->hostSeqInfo);
    free(this->hostParam);
    free(this->hostFreqPointers);
    free(this->hostGapOpPointers);
    free(this->hostGapExPointers);
    free(this->hostAlnPointers);
    return;
}

void msa::progressive::gpu::GPU_pointers::memcpyHost2Device(Option *option, int len, int alnPairs, int numBlocks)
{
    cudaSetDevice(this->gpu_index);
    int paramSize = (option->type == 'n') ? 5 * 5 + 4 : 21 * 21 + 4;
    int profileSize = (option->type == 'n') ? 5 + 1 : 21 + 1;
    int pairPerMemBlock = numBlocks/(this->memBlock);
    int usedMemBlock = (alnPairs % pairPerMemBlock == 0) ? alnPairs/pairPerMemBlock : alnPairs/pairPerMemBlock+1;
    for (int i = 0; i < usedMemBlock; ++i) {
        cudaMemcpy(this->deviceFreq[i],  this->hostFreq[i], profileSize * 2 * len * pairPerMemBlock * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(this->deviceGapOp[i], this->hostGapOp[i],              2 * len * pairPerMemBlock * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(this->deviceGapEx[i], this->hostGapEx[i],              2 * len * pairPerMemBlock * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(this->deviceAln[i],   this->hostAln[i],                2 * len * pairPerMemBlock * sizeof(int8_t), cudaMemcpyHostToDevice);
    }
    cudaMemcpy(this->deviceFreqPointers,  this->hostFreqPointers,  this->memBlock * sizeof(float*),  cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceGapOpPointers, this->hostGapOpPointers, this->memBlock * sizeof(float*),  cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceGapExPointers, this->hostGapExPointers, this->memBlock * sizeof(float*),  cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceAlnPointers,   this->hostAlnPointers,   this->memBlock * sizeof(int8_t*), cudaMemcpyHostToDevice);

    cudaMemcpy(this->deviceLen,      this->hostLen,                2       * alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceNum,      this->hostNum,                2       * alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceAlnLen,   this->hostAlnLen,                       alnPairs * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceSeqInfo,  this->hostSeqInfo,            4 *                  sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(this->deviceParam,    this->hostParam,                       paramSize * sizeof(float),   cudaMemcpyHostToDevice);
    std::string error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to copy memory [%s] to GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::GPU_pointers::memcpyDevice2Host(Option *option, int len, int alnPairs, int numBlocks)
{
    cudaSetDevice(this->gpu_index);
    int pairPerMemBlock = numBlocks/this->memBlock;
    int usedMemBlock = (alnPairs % pairPerMemBlock == 0) ? alnPairs/pairPerMemBlock : alnPairs/pairPerMemBlock+1;
    for (int i = 0; i < usedMemBlock; ++i) {
        cudaMemcpy(this->hostAln[i], this->deviceAln[i], 2 * len * pairPerMemBlock * sizeof(int8_t),  cudaMemcpyDeviceToHost);
    }
    cudaMemcpy(this->hostAlnLen, this->deviceAlnLen,           alnPairs * sizeof(int32_t), cudaMemcpyDeviceToHost);
    std::string error = cudaGetErrorString(cudaGetLastError());
    if (error != "no error") fprintf(stderr, "CUDA Error: Failed to copy memory [%s] from GPU %d\n", error.c_str(), this->gpu_index);
    return;
}

void msa::progressive::gpu::parallelAlignmentGPU(Tree *tree, NodePairVec &nodes, SequenceDB *database, Option *option, Params &param)
{
    auto startT = std::chrono::high_resolution_clock::now();
    int blockSize = device_function::_THREAD_NUM;
    int gpuNum = option->gpuNum;
    int profileSize = param.matrixSize + 1;
    const float gpuMemUsage = 0.7;

    // get maximum sequence/profile length
    int32_t seqLen = 0;
    const int32_t KB = 1024;
    int32_t neededBytes = (profileSize * 2 + 2) * sizeof(float) + 2 * sizeof(int8_t);
    for (auto n : nodes) seqLen = std::max(seqLen, std::max(n.first->getAlnLen(database->currentTask), n.second->getAlnLen(database->currentTask)));
    auto minMem = std::min_element(option->gpuMem.begin(), option->gpuMem.end());
    int32_t requiredMem = (((neededBytes * seqLen)/KB));
    int32_t availableMem = static_cast<int>((*minMem)*gpuMemUsage)*KB;
    int numBlocks = std::min(device_function::_BLOCKSIZE, availableMem/requiredMem);
    if (numBlocks < device_function::_BLOCKSIZE) {
        int new_numBlocks = 1;
        while (new_numBlocks <= numBlocks) new_numBlocks <<= 1;
        new_numBlocks >>= 1;
        numBlocks = new_numBlocks;
    }
    // std::cout << numBlocks << ',' << availableMem << ',' << (availableMem/requiredMem) << '\n';
    int memBlock = (option->type == 'n') ? 
                   ((seqLen * numBlocks) / (_MAX_LENGTH_N * device_function::_BLOCKSIZE) + 1):
                   ((seqLen * numBlocks) / (_MAX_LENGTH_P * device_function::_BLOCKSIZE) + 1);
    switch (memBlock) {
        case 1:  memBlock=1; break;
        case 2:  memBlock=2; break;
        case 3:  memBlock=4; break;
        case 4:  memBlock=4; break;
        case 5:  memBlock=8; break;
        case 6:  memBlock=8; break;
        case 7:  memBlock=8; break;
        case 8:  memBlock=8; break;
        default: memBlock=1; break;
    }
    int pairPerMemBlock = numBlocks / memBlock;
    if (pairPerMemBlock > nodes.size()) memBlock = 1;

    // std::cout << pairPerMemBlock << ';' << memBlock << '\t' << seqLen << '\n'; 

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
                int32_t offset_g =               2 * gpuMemLen * (n%pairPerMemBlock);
                IntPair lens = {refLen, qryLen};
                float* rawProfile = new float [profileSize * 2 * gpuMemLen];
                for (int i = 0; i < profileSize * 2 * gpuMemLen; ++i) rawProfile[i] = 0.0;
                alignment_helper::calculateProfile(rawProfile, nodes[nIdx], database, option, gpuMemLen); 
                // Get consensus sequences
                alignment_helper::getConsensus(option, rawProfile,                       consensus[n].first,  lens.first);
                alignment_helper::getConsensus(option, rawProfile+profileSize*gpuMemLen, consensus[n].second, lens.second);
                // Copy profile to host memory
                for (int i = 0; i < profileSize * gpuMemLen; ++i) gp->hostFreq[memBlockIdx][offset_f+i] =                       (i < profileSize * lens.first)  ? rawProfile[i] : 0.0;
                for (int i = 0; i < profileSize * gpuMemLen; ++i) gp->hostFreq[memBlockIdx][offset_f+profileSize*gpuMemLen+i] = (i < profileSize * lens.second) ? rawProfile[profileSize*gpuMemLen+i] : 0.0;
                std::pair<int32_t, int32_t> offset = std::make_pair(offset_f, offset_g);
                alignment_helper::calculatePSGP(gp->hostFreq[memBlockIdx], gp->hostGapOp[memBlockIdx], gp->hostGapEx[memBlockIdx], nodes[nIdx], database, option, gpuMemLen, offset, lens, param);
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
            std::string berr = cudaGetErrorString(cudaGetLastError());
            if (berr != "no error") printf("ERROR: Before kernel %s!\n", berr.c_str());
            auto kernelStart = std::chrono::high_resolution_clock::now();
            device_function::parallelProfileAlignment<<<numBlocks, blockSize>>>(
                gp->deviceFreqPointers,
                gp->deviceAlnPointers,
                gp->deviceLen,
                gp->deviceNum,
                gp->deviceAlnLen,
                gp->deviceSeqInfo,
                gp->deviceGapOpPointers,
                gp->deviceGapExPointers,
                gp->deviceParam
            );
            cudaDeviceSynchronize();
            std::string aerr = cudaGetErrorString(cudaGetLastError());
            if (aerr != "no error") printf("ERROR: After kernel %s!\n", aerr.c_str());
            auto kernelEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds kTime = kernelEnd - kernelStart;
            int kt = kernelTime.fetch_add(kTime.count());
            gp->memcpyDevice2Host(option, gpuMemLen, alnPairs, numBlocks);
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
                        int offset = 2 * gpuMemLen * (n%pairPerMemBlock);
                        int memBlockIdx = (n/pairPerMemBlock);
                        if (gp->hostAlnLen[n] > 0) for (int j = 0; j < gp->hostAlnLen[n]; ++j) aln_wo_gc.push_back(gp->hostAln[memBlockIdx][offset+j]);
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
                        // if (alnRef != gp->hostLen[2*n])   std::cout << "R: Pre " << nodes[nIdx].first->identifier << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
                        // if (alnQry != gp->hostLen[2*n+1]) std::cout << "Q: Pre " << nodes[nIdx].second->identifier << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";

                        // Add gappy columns back
                        IntPair debugIdx = std::make_pair(alnRef,alnQry);
                        msa::alignment_helper::addGappyColumnsBack(aln_wo_gc, aln[n], gappyColumns[n], param, {alnRef,alnQry}, consensus[n]);
                        alnRef = 0, alnQry = 0;
                        for (auto a: aln[n]) {
                            if (a == 0) {alnRef += 1; alnQry += 1;}
                            if (a == 1) {alnQry += 1;}
                            if (a == 2) {alnRef += 1;}
                        }
                        // if (alnRef != refLen) {
                        //     std::cout << "R: Post " << nodes[nIdx].first->identifier <<  '\t' << refNum << "(" << alnRef << "/" << nodes[nIdx].first->getAlnLen(database->currentTask) << ")\n";
                        //     for (auto gp: gappyColumns[n].first) std::cout << '(' << gp.first << ',' << gp.second << ')' << ',';
                        //     std::cout << '\n';
                        // }
                        // if (alnQry != qryLen) {
                        //     std::cout << "Q: Post " << nodes[nIdx].second->identifier << '\t' << qryNum << "(" << alnQry << "/" << nodes[nIdx].second->getAlnLen(database->currentTask) << ")\n";
                        //     for (auto gp: gappyColumns[n].second) std::cout << '(' << gp.first << ',' << gp.second << ')' << ',';
                        //     std::cout << '\n';
                        // }

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
                    float refWeight = nodes[nIdx].first->alnWeight, qryWeight = nodes[nIdx].second->alnWeight;
                    alignment_helper::updateAlignment(nodes[nIdx], database, option, aln[n]);
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

    std::cerr << "Kernel Time: " << kernelTime / 1000000 << "ms\n";

    // free memory
    if (fallbackPairs.empty()) return;
    // fallback to CPU
    alignment_helper::fallback2cpu(fallbackPairs, nodes, database, option);
    return;
}
