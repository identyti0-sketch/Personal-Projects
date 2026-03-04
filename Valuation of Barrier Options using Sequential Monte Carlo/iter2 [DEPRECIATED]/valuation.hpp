#ifndef VALUATION
#define VALUATION

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>
#include <numeric>

using namespace std;


class Valuation {
public:
    Valuation(string& filename);

    void repeatExperiment(int n);

    /// @brief Update the number of steps taken for the option to mature 
    /// and recalculate the step volatility
    /// @param steps 
    void updateSteps(int steps);

protected: 
    double S0;      // initial asset price
    double K;       // strike price
    double U;       // Upper barrier
    double L;       // Lower barrier
    double r;       // risk-neutral drift
    double sigmaYearly; // yearly volatility
    double T;       // time to maturity
    int N;          // number of time steps
    int M;          // number of particles
    int n;          // times to repeat experiment
    mt19937 gen;
    normal_distribution<double> dist;

    struct result {
        vector<double> price;
        vector<double> survivalRate;
    };

    virtual void runSimulation(struct result& res);

    double stdErr(struct result& res, double mean);

    // returns true if the particle hits either barrier, false otherwise
    bool ifHitBarrier(double price);

    double getDrift();
};


class ValuationSMC : public Valuation {
public:

    using Valuation::Valuation;

private: 


    vector<double> particles;

    // runs M simulations of the stock price in parallel
    void runSimulation(struct result& res);

    // moves the particles one step forward, updates the survival rate and resamples the particles
    void step(double& survivalRate);

    void increment(double& survivalRate, vector<int>& alive, vector<int>& dead);

    
    // randomly (uniform distribution) copy living particles onto dead particles.
    void resample(vector<int>& alive, vector<int>& dead);

    // calculates the payoff of the option based on the average of the final particle states
    // does not account for the discount factor or the survival rate
    double getPayoff();
};

class ValuationMC : public Valuation {
public:
    using Valuation::Valuation;
private:
    void runSimulation(struct result& res);
    
    double simulatePrice(bool& survived);

};

#endif