#include "valuation.hpp"

int main() {
    string filename = "example1.txt";
    ValuationContinuousSMC smc = ValuationContinuousSMC(filename);
    ofstream out("output.csv");
    if (!out) {
        cerr << "Error creating file!\n";
        return 1;
    }
    out << "N,Average Price,Average Survival Rate,Standard Error,Time (s)" << endl;
    out.close();
    for (int i = 1; i < 200; i*= 2) {
        smc.updateSteps(i);
        smc.repeatExperiment(50);
    }
    return 0;
}