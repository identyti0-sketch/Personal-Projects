#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>

using namespace std;



struct Setting {
    double S0;      // initial asset price
    double K;       // strike price
    double U;       // Upper barrier
    double L;       // Lower barrier
    double r;       // risk-neutral drift
    double sigmaYearly; // yearly volatility
    double sigma;   // volatility of each iteration
    double T;       // time to maturity
    int N;          // number of time steps
    int M;          // number of particles
    int n;          // times to repeat experiment
    mt19937 gen;
    normal_distribution<double> dist;
};


struct result {
    vector<double> price;
    vector<double> survivalRate;
};

// struct Setting setSettings();
void runSimulation(struct Setting& setting, struct result& res);
void step(vector<double>& particles, struct Setting& setting, double& survivalRate);
void increment(vector<double>& particles, struct Setting& setting, double& survivalRate, vector<int>& alive, vector<int>& dead);
bool ifHitBarrier(double price, struct Setting& setting);
void resample(vector<double>& particles, vector<int>& alive, vector<int>& dead);
double getPayoff(struct Setting& setting, vector<double>& particles);
struct Setting getSettings();
void repeatExperiment(struct Setting& setting);

int main() {
    struct Setting setting = getSettings();
    // Initialize particles
    for (int i = 1; i < 200; i *= 2) {
        setting.N = i;
        setting.sigma = setting.sigmaYearly * sqrt(setting.T / setting.N);
        setting.dist = normal_distribution<double>(0, setting.sigma);
        repeatExperiment(setting);
    }
    return 0;

}
/*
struct Setting setSettings() {
    struct Setting setting;
    cout << "Enter initial asset price S0: ";
    cin >> setting.S0;
    cout << "Enter strike price K: ";
    cin >> setting.K;
    cout << "Enter upper barrier U: ";
    cin >> setting.U;
    cout << "Enter lower barrier L: ";
    cin >> setting.L;
    double i;
    cout << "Enter interest rate r: ";
    cin >> i;
    double q;
    cout << "Enter dividend yield q: ";
    cin >> q;
    setting.r = i - q;
    cout << "Enter volatility sigma: ";
    cin >> setting.sigma;
    cout << "Enter time (years) to maturity T: ";
    cin >> setting.T;
    cout << "Enter number of time steps N: ";
    cin >> setting.N;
    cout << "Enter number of particles M: ";
    cin >> setting.M;
    setting.sigma = setting.sigma * sqrt(setting.T / setting.N);
    setting.gen = mt19937(random_device{}());
    setting.dist = normal_distribution<double>(0, setting.sigma);
    return setting;
}
*/
struct Setting getSettings() {
    struct Setting setting;
    setting.S0 = 100.0;
    setting.K = 100.0;
    setting.U = 110.0;
    setting.L = 90.0;
    setting.r = 0.1;
    setting.sigmaYearly = 0.3;
    setting.T = 0.5;
    setting.N = 1;
    setting.sigma = setting.sigmaYearly * sqrt(setting.T / setting.N);
    setting.M = 100000;
    setting.n = 50;
    setting.gen = mt19937(random_device{}());
    setting.dist = normal_distribution<double>(0, setting.sigma);
    return setting;
}

void repeatExperiment(struct Setting& setting) {
    struct result res;
    for (int i = 0; i < setting.n; i++) {
        runSimulation(setting, res);
    }
    double averagePrice = accumulate(res.price.begin(), res.price.end(), 0.0) / res.price.size();
    double averageSurvivalRate = accumulate(res.survivalRate.begin(), res.survivalRate.end(), 0.0) / res.survivalRate.size();
    cout <<"N: " << setting.N << " Average Price: " << averagePrice << " Average Survival Rate: " << averageSurvivalRate << endl; 
}

void runSimulation(struct Setting& setting, struct result& res) {
    vector<double> particles(setting.M, setting.S0);
    double survivalRate = 1.0;
    for (int i = 0; i < setting.N; i++) {
        step(particles, setting, survivalRate);
    }
    double discountFactor = exp(-setting.r * setting.T);
    double payoff = getPayoff(setting, particles);
    res.price.push_back(payoff * discountFactor * survivalRate);
    res.survivalRate.push_back(survivalRate);
}

void step(vector<double>& particles, struct Setting& setting, double& survivalRate) {
    vector<int> alive;
    vector<int> dead;

    increment(particles, setting, survivalRate, alive, dead);    
    survivalRate *= (double)alive.size() / particles.size();
    resample(particles, alive, dead);
}


void increment(vector<double>& particles, struct Setting& setting, double& survivalRate, vector<int>& alive, vector<int>& dead) {

    for (int i = 0; i < particles.size(); i++) {
        particles[i] *= 1 + setting.dist(setting.gen);   
        if (ifHitBarrier(particles[i], setting)) {
            dead.push_back(i);
        } else {
            alive.push_back(i);
        }
    }
}

bool ifHitBarrier(double price, struct Setting& setting) {
    if (price >= setting.U || price <= setting.L) {
        return true;
    }
    return false;
}


void resample(vector<double>& particles, vector<int>& alive, vector<int>& dead) {
    for (int i = 0; i < dead.size(); i++) {
        int index = rand() % alive.size();
        particles[dead[i]] = particles[alive[index]];
    }
}

double getPayoff(struct Setting& setting, vector<double>& particles) {
    double sum = 0.0;
    for (double price : particles) {
        sum += max(0.0, price - setting.K);
    }
    return sum / particles.size();
}