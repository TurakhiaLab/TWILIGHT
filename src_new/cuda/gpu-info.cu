#ifndef MSA_HPP
#include "../msa.hpp"
#endif


void msa::Option::getGpuInfo(po::variables_map &vm) {

    int maxGpuNum = 0;
    cudaError_t err = cudaGetDeviceCount(&maxGpuNum);
    if (err != cudaSuccess) {
        std::cerr << "No GPU detected." << std::endl;
        maxGpuNum = 0; // or handle accordingly
    }
    int gpuNum = (vm.count("gpu")) ? vm["gpu"].as<int>() : maxGpuNum;
    if (gpuNum < 0 || gpuNum > maxGpuNum)
    {
        std::cerr << "ERROR: Invalid number of GPUs. Please request between 0 and " << maxGpuNum << ".\n";;
        exit(1);
    }
    
    if (this->cpuOnly) gpuNum = 0;
    else if (gpuNum > this->cpuNum) {
        if (this->cpuNum == 1) std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << this->cpuNum << " GPU.\n";
        else std::cerr << "WARNING: Requesting more GPUs than requested CPU threads, and will only execute on " << this->cpuNum << " GPUs.\n";
        gpuNum = this->cpuNum;
    }
    std::vector<int> gpuIdx;
    if (vm.count("gpu-index")) {
        std::string gpuIdxString = vm["gpu-index"].as<std::string>();
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
        std::cerr << "ERROR: the number of requested GPUs does not match the number of specified gpu indexes.\n";
        exit(1);
    }
    for (auto id : gpuIdx) {
        if (id >= maxGpuNum)
        {
            std::cerr << "ERROR: specified gpu index >= the number of GPUs\n";
            exit(1);
        }
    }
    fprintf(stderr, "Maximum available GPUs: %d. Using %d GPUs.\n", maxGpuNum, gpuNum);
    if (gpuNum == 0) {
        std::cerr << "Requested 0 GPU or no GPU resources detected; switching to CPU version.\n";
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
            std::cerr << "Error: " << cudaGetErrorString(status) << std::endl;
            exit(1);
        }
        gpuTotalMem.push_back(total_bytes/MB);
        gpuMem.push_back(free_bytes/MB);
    }
    
    this->gpuNum = gpuNum;
    this->gpuIdx = gpuIdx;
    this->gpuMem = gpuMem;

    if (this->gpuNum == 0) return;

    
    std::cerr << "========= GPU Info =========\n";
    for (int i = 0; i < this->gpuIdx.size(); ++i) {
        std::cerr << "GPU " << this->gpuIdx[i] << ", Available Memory: "<< (this->gpuMem[i]) << " / " << (gpuTotalMem[i]) << " MB.\n";
    }
    return;
}