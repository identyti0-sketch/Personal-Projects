#include "valuation.hpp"

int main() {
    string filename = "example1.txt";
    ValuationContinuousSMC smc = ValuationContinuousSMC(filename);
    for (int i = 1; i < 200; i*= 2) {
        smc.updateSteps(i);
        smc.repeatExperiment(50);
    }
    return 0;
}