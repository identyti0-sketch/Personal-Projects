#include "valuation.hpp"

int main() {
    string filename = "example1.txt";
    ValuationSemiContinuousSMC smc = ValuationSemiContinuousSMC(filename);
    for (int i = 1; i < 200; i*= 2) {
        smc.updateSteps(i);
        smc.repeatExperiment(50);
    }
    /*
    ValuationMC mc = ValuationMC(filename);
    
    
    for (int i = 1; i < 200; i*= 2) {
        mc.updateSteps(i);
        mc.repeatExperiment(50);
    }
    */
    return 0;
}