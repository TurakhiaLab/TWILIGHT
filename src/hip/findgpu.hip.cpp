
#include <hip/hip_runtime.h>
#include <iostream>

int main (){
    // A simple code to detect cuda-capable GPUs
    int device;
    hipDeviceProp_t prop;
    if (hipGetDevice(&device) != hipSuccess) return 1;
    if (hipGetDeviceProperties(&prop, device) != hipSuccess) return 1;
    std::cout << prop.major << prop.minor;
    return 0;
}