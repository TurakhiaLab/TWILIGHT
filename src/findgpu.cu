int main (){
    // A simple code to detect cuda-capable GPUs
    int deviceCount;
    cudaError_t e = cudaGetDeviceCount(&deviceCount);
    return e == cudaSuccess ? 0 : 1;
}