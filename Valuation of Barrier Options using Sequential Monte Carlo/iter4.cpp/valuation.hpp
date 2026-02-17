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
    double S0;      // log of initial asset price
    double K;       // strike price
    double U;       // log of Upper barrier
    double L;       // log of Lower barrier
    double r;       // risk-neutral drift
    double sigmaYearly; // yearly volatility
    double T;       // time to maturity
    int N;          // number of time steps
    int M;          // number of particles
    int n;          // times to repeat experiment
    mt19937_64 gen;
    normal_distribution<double> dist;

    struct result {
        vector<double> price;
        vector<double> survivalRate;
    };

    virtual void runSimulation(struct result& res);

    double stdErr(struct result& res, double mean);
    double updateweight(double price, double prevPrice);
    // returns true if the particle hits either barrier, false otherwise
    bool ifHitBarrier(double price);

    double getDrift();
};


class ValuationSMC : public Valuation {
public:

    using Valuation::Valuation;

protected: 


    vector<double> particles;

    // runs M simulations of the stock price in parallel
    void runSimulation(struct result& res);

    // moves the particles one step forward, updates the survival rate and resamples the particles
    void step(double& survivalRate);

    virtual void increment(vector<int>& alive, vector<int>& dead);

    
    // randomly (uniform distribution) copy living particles onto dead particles.
    void resample(vector<int>& alive, vector<int>& dead);

    // calculates the payoff of the option based on the average of the final particle states
    // does not account for the discount factor or the survival rate
    double getPayoff();
};


class ValuationSemiContinuousSMC : public ValuationSMC {
    public:
    using ValuationSMC::ValuationSMC;
    private:
    void increment(vector<int>& alive, vector<int>& dead) override;

    bool ifHitBarrier(double price, double prevPrice);

};

class ValuationMC : public Valuation {
public:
    using Valuation::Valuation;
private:
    void runSimulation(struct result& res);
    
    double simulatePrice(bool& survived);

};

#endif


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
    void step();

    void increment(vector<int>& alive, vector<int>& dead);

    
    // randomly (proportionally weighted distribution) copy living particles onto dead particles.
    // returns totalweight/newtotalweight
    void resample(vector<int>& alive, vector<int>& dead);


    
    // calculates the payoff of the option based on the average of the final particle states
    // does not account for the discount factor or the survival rate
    double getPayoff();
};



vector<double> getRandomVector(int n, double max);

void getRandomHelper(vector<double>& randoms, int n, double max, double min,
    mt19937_64& rng, uniform_real_distribution<double>& dist);