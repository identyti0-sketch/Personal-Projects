#include <omp.h>
#include <iostream>
#include <omp_llvm.h>

int main() {

    std::cout << "Num devices: " << omp_get_num_devices() << "\n";
    #pragma omp target
    {
        std::cout << "Hello from device " << omp_get_device_num() << "\n";
    } 
    return 0;
}