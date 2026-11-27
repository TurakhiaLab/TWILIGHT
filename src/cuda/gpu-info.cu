#ifndef MSA_HPP
#include "../msa.hpp"
#endif


void msa::Option::getGpuInfo(const CliOptions &opt) {

    int maxGpuNum = 0;
    cudaError_t err = cudaGetDeviceCount(&maxGpuNum);
    if (err != cudaSuccess) {
        maxGpuNum = 0; // or handle accordingly
    }

    if (this->accurate) {
        msg::info("Maximum available GPUs: %d. Using 0 GPUs (Current accurate mode only has CPU version; no GPUs will be used).\n", maxGpuNum);
        this->cpuOnly = true;
        this->gpuNum = 0;
        return;
    }
    int gpuNum = (opt.hasGpu) ? opt.gpuNum : maxGpuNum;
    if (gpuNum < 0 || gpuNum > maxGpuNum)
    {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid number of GPUs. Please request between 0 and " + std::to_string(maxGpuNum) + ".");
    }
    
    if (this->cpuOnly) gpuNum = 0;
    else if (gpuNum > this->cpuNum) {
        if (this->cpuNum == 1) msg::warn("Requesting more GPUs than requested CPU threads, and will only execute on %d GPU.\n", this->cpuNum);
        else msg::warn("Requesting more GPUs than requested CPU threads, and will only execute on %d GPUs.\n", this->cpuNum);
        gpuNum = this->cpuNum;
    }
    std::vector<int> gpuIdx;
    if (!opt.gpuIndex.empty()) {
        std::string gpuIdxString = opt.gpuIndex;
        std::string id = "";
        for (int i = 0; i < gpuIdxString.size(); ++i) {
            if (gpuIdxString[i] == ',') {
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else if (i == gpuIdxString.size() - 1) {
                id += gpuIdxString[i];
                gpuIdx.push_back(std::atoi(id.c_str()));
                id = "";
            }
            else id += gpuIdxString[i];
        }
    }
    else {
        for (int i = 0; i < gpuNum; ++i)
            gpuIdx.push_back(i);
    }

    if (gpuIdx.size() != gpuNum) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "The number of requested GPUs does not match the number of specified gpu indexes.");
    }
    for (auto id : gpuIdx) {
        if (id >= maxGpuNum)
        {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Specified gpu index >= the number of GPUs.");
        }
    }
    msg::info("Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);
    if (gpuNum == 0) {
        msg::info("Requested 0 GPU or no GPU resources detected; switching to CPU version.\n");
        this->cpuOnly = true;
    }

    std::vector<size_t> gpuMem;
    std::vector<size_t> gpuTotalMem;
    const size_t MB = 1024 *1024;

    for (auto idx: gpuIdx) {
        cudaSetDevice(idx);
        size_t free_bytes, total_bytes;
        cudaError_t status = cudaMemGetInfo(&free_bytes, &total_bytes);
        if (status != cudaSuccess) {
            throw msa::TwilightError(msa::ErrorCode::INTERNAL,
                std::string("CUDA error: ") + cudaGetErrorString(status));
        }
        gpuTotalMem.push_back(total_bytes/MB);
        gpuMem.push_back(free_bytes/MB);
    }
    
    this->gpuNum = gpuNum;
    this->gpuIdx = gpuIdx;
    this->gpuMem = gpuMem;

    if (this->gpuNum == 0) return;

    
    msg::info("========= GPU Info =========\n");
    for (int i = 0; i < this->gpuIdx.size(); ++i) {
        msg::info("GPU %d, Available Memory: %zu / %zu MB.\n", this->gpuIdx[i], this->gpuMem[i], gpuTotalMem[i]);
    }
    return;
}