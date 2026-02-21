#include "valuation.hpp"

int main() {
    string filename = "example1.txt";
    ValuationSMC smc = ValuationSMC(filename);
    smc.updateSteps(100000);
    smc.repeatExperiment(1);
    /*
    ValuationMC mc = ValuationMC(filename);
    for (int i = 1; i < 200; i*= 2) {
        smc.updateSteps(i);
        smc.repeatExperiment(50);
    }
    
    for (int i = 1; i < 200; i*= 2) {
        mc.updateSteps(i);
        mc.repeatExperiment(50);
    }
    */
    return 0;
}