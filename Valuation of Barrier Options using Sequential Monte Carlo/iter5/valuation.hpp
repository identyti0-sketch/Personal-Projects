#ifndef VALUATION
#define VALUATION

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>
#include <numeric>
#include <omp.h>
#include "openrand/philox.h"
#include <memory>

using namespace std;


class Rng {
public: 
    Rng(uint64_t seed, uint32_t counter, double sigma);
    double getRand();
    normal_distribution<double> dist;
private:
    openrand::Philox rng;
};

class Valuation {
public:
    Valuation(const string& filename);

    void repeatExperiment(int n);

    /// @brief Update the number of steps taken for the option to mature 
    /// and recalculate the step volatility
    /// @param steps 
    void updateSteps(int steps);

protected: 
    double S0;      // log of initial asset price
    double K;       // strike price
    double U;       // log of Upper barrier
    double L;       // log of Lower barrier
    double r;       // risk-neutral drift
    double sigmaYearly; // yearly volatility
    double T;       // time to maturity
    double sigmaStep; // volatility of the time frame of a step;
    int N;          // number of time steps
    int M;          // number of particles
    int n;          // times to repeat experiment
    unsigned int seed;
    openrand::Philox gen;
    normal_distribution<double> dist;

    struct result {
        vector<double> price;
        vector<double> survivalRate;
    };
    static thread_local unique_ptr<Rng> rngEngine;
    static thread_local bool rngInitialized;
    virtual void runSimulation(struct result& res);

    double stdErr(struct result& res, double mean);
    double updateweight(double price, double prevPrice);
    // returns true if the particle hits either barrier, false otherwise
    bool ifHitBarrier(double price);

    double getDrift();

    inline void initRng();
};


class ValuationSMC : public Valuation {
public:

    using Valuation::Valuation;

protected: 


    vector<double> particles;

    // runs M simulations of the stock price in parallel
    void runSimulation(struct result& res);

    // moves the particles one step forward, updates the survival rate and resamples the particles
    bool step(double& survivalRate);

    virtual void increment(vector<int>& alive, vector<int>& dead);

    
    // randomly (uniform distribution) copy living particles onto dead particles.
    void resample(vector<int>& alive, vector<int>& dead);

    // calculates the payoff of the option based on the average of the final particle states
    // does not account for the discount factor or the survival rate
    double getPayoff();
};



class ValuationMC : public Valuation {
public:
    using Valuation::Valuation;
protected:
    void runSimulation(struct result& res);
    
    virtual double simulatePrice(double& survivedg);

};

class ValuationContinuousMC : public ValuationMC {
public:
    using ValuationMC::ValuationMC;
private:
    double simulatePrice(double& survived) override;
};




class ValuationContinuousSMC : public Valuation {
public:

    using Valuation::Valuation;

private: 
    class Particle
    {
    public:
        double p; // price
        double w; // weight
    };
    double totalWeight;
    vector<Particle> particles;
    

    // runs M simulations of the stock price in parallel
    void runSimulation(struct result& res);

    // moves the particles one step forward, updates the survival rate and resamples the particles
    bool step();

    int increment(vector<char>& deadMask);

    void unmask(vector<Particle *>& alive, vector<Particle *>& dead, vector<char>& deadMask, vector<double>& prefixSum);    
    // randomly (proportionally weighted distribution) copy living particles onto dead particles.
    // returns totalweight/newtotalweight
    void resample(vector<Particle *>& alive, vector<Particle *>& dead, vector<double>& prefixSum);

    void partitionResample(vector<Particle *>& alive, vector<Particle *>& dead, vector<double>& randVect, vector<double>& prefixSum, int startEnd[4]);
    
    // calculates the payoff of the option based on the average of the final particle states
    // does not account for the discount factor or the survival rate
    double getPayoff();
};



vector<double> getRandomVector(int n, double max);

void getRandomHelper(vector<double>& randoms, int n, double max, double min,
    openrand::Philox& rng, uniform_real_distribution<double>& dist);

int binSearchSmallestAbove(double target, vector<double>& arr, int low, int high);
int binSearchGreatestBelow(double target, vector<double>& arr, int low, int high);
#endif