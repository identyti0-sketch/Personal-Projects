#include <omp.h>
#include <iostream>

int main() {

    std::cout << "Num devices: " << omp_get_num_devices() << "\n";
    return 0;
}