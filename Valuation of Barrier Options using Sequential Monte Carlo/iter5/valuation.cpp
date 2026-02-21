#include <ctime>
#include "valuation.hpp"

Valuation::Valuation(string& filename) {
    ifstream file(filename);
    
    file >> S0 >> K >> U >> L >> r >> sigmaYearly >> T >> N >> M;
    S0 = log(S0);
    U = log(U);
    L = log(L);
    double sigma = sigmaYearly * sqrt(T / N);
    gen = mt19937_64(random_device{}());  
    dist = normal_distribution<double>(0, sigma);
}


void Valuation::repeatExperiment(int n) {
    struct result res;
    clock_t startTime = clock();
    for (int i = 0; i < n; i++) {
        runSimulation(res);
    }
    clock_t endTime = clock();
    double elapsedTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
    double averagePrice = accumulate(res.price.begin(), res.price.end(), 0.0) / res.price.size();
    double averageSurvivalRate = accumulate(res.survivalRate.begin(), res.survivalRate.end(), 0.0) / res.survivalRate.size();
    double standardErr = stdErr(res, averagePrice);
    ofstream out("output.csv", ios::app);
    if (!out) {
        cerr << "Error opening file!\n";
        return;
    }
    out << N << "," << averagePrice << "," << averageSurvivalRate << "," << standardErr << "," << elapsedTime << endl;
    out.close();
}

void Valuation::updateSteps(int steps) {
    N = steps;
    double sigma = sigmaYearly * sqrt(T / N);
    dist = normal_distribution<double>(0, sigma);
}

double Valuation::stdErr(struct result& res, double mean) {
    size_t n = res.price.size();
    if (n < 2 || mean == 0.0) {
        return 0.0;
    }
    double var = 0.0;
    for (double price : res.price) {
        double diff = price - mean;
        var += diff * diff;
    }
    var /= (n - 1);
    return sqrt(var / n) / mean * 100.0;
}

bool Valuation::ifHitBarrier(double price) {
    if (price >= U || price <= L) {
        return true;
    }
    return false;
}

void Valuation::runSimulation(struct result& res) {}

double Valuation::getDrift()  {
    return (r - 0.5 * sigmaYearly * sigmaYearly) * (T / N) + dist(gen);
}

double Valuation::updateweight(double price, double prevPrice) {
    double doubleVar = 2* (sigmaYearly*sigmaYearly*T/N);
    double priceDiff = 2*(price - prevPrice);
    double barrierDist = 2*(U - L);
    double distToUpper = 2*(U - prevPrice);
    double distToLower = 2*(prevPrice - L);
    
    auto brownianBridge = [](double z, double x, double doubleVar) {
        return exp(- z*(z - x)/doubleVar);
    };

    double sum = 1.0;
    for (int i = 1; i < 4; i++) {
        double bdi = barrierDist*i;
        sum -= brownianBridge(bdi - distToLower, priceDiff, doubleVar)
                + brownianBridge(-bdi + distToUpper, priceDiff, doubleVar);
        sum += brownianBridge(bdi, priceDiff, doubleVar)
                + brownianBridge(-bdi, priceDiff, doubleVar);
    }
    return sum;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ValuationSMC::runSimulation(struct result& res) {
    particles = vector<double>(M, S0);
    double survivalRate = 1.0;
    for (int i = 0; i < N; i++) {
        step(survivalRate);
    }
    double discountFactor = exp(-r * T);
    double payoff = getPayoff();
    res.price.push_back(payoff * discountFactor * survivalRate);
    res.survivalRate.push_back(survivalRate);
}

void ValuationSMC::step(double& survivalRate) {
    vector<int> alive;
    vector<int> dead;

    increment(alive, dead);    
    survivalRate *= (double)alive.size() / particles.size();
    resample(alive, dead);
}


void ValuationSMC::increment(vector<int>& alive, vector<int>& dead) {
    for (int i = 0; i < particles.size(); i++) {
        particles[i] += getDrift();   
        if (ifHitBarrier(particles[i])) {
            dead.push_back(i);
        } else {
            alive.push_back(i);
        }
    }
}


void ValuationSMC::resample(vector<int>& alive, vector<int>& dead) {
    for (int i : dead) {
        int index = rand() % alive.size();
        particles[i] = particles[alive[index]];
    }
}

double ValuationSMC::getPayoff() {
    double sum = 0.0;
    for (double price : particles) {
        sum += max(0.0, exp(price) - K);
    }
    return sum / particles.size();
}


//////////////////////////////////////////////////////////////////////////////////////////////

void ValuationMC::runSimulation(struct result& res) {
    double averagePrice = 0.0;
    double survivalCount = 0.0;
    for (int i = 0; i < M; i++) {
        double survived = 1.0;
        averagePrice += simulatePrice(survived);
        survivalCount += survived;
    }
    double discountFactor = exp(-r * T);
    res.price.push_back(averagePrice / M * discountFactor);
    res.survivalRate.push_back(survivalCount / M);
}

double ValuationMC::simulatePrice(double& survived) {
    double value = S0;
    for (int i = 0; i < N; i++) {
        value += getDrift();     
        if (ifHitBarrier(value)) {
            survived = 0.0;
            return 0.0;
        }
    }
    return max(exp(value) - K, 0.0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////


void ValuationSemiContinuousSMC::increment(vector<int>& alive, vector<int>& dead) {

    for (int i = 0; i < particles.size(); i++) {
        double prev = particles[i];
        particles[i] += getDrift();   
        if (ifHitBarrier(particles[i], prev)) {
            dead.push_back(i);
        } else {
            alive.push_back(i);
        }
    }
}

bool ValuationSemiContinuousSMC::ifHitBarrier(double price, double prevPrice) {
    if (price >= U || price <= L) {
        return true;
    }
    return updateweight(price, prevPrice) < (double)rand() / RAND_MAX;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ValuationContinuousSMC::runSimulation(struct result& res) {
    particles = vector<Particle>(M, Particle{S0, 1.0});
    totalWeight = M;
    for (int i = 0; i < N; i++) {
        step();
    }
    double discountFactor = exp(-r * T);
    double payoff = getPayoff();
    res.price.push_back(payoff * discountFactor);
    res.survivalRate.push_back(totalWeight / M);
}

void ValuationContinuousSMC::step() {
    vector<Particle *> alive;
    vector<Particle *> dead;

    increment(alive, dead);    
    resample(alive, dead);
}


void ValuationContinuousSMC::increment(vector<Particle *>& alive, vector<Particle *>& dead) {
    totalWeight = 0.0;
    for (Particle& particle : particles) {
        double prev = particle.p;
        particle.p += getDrift();   
        if (ifHitBarrier(particle.p)) {
            particle.w = 0.0;
            dead.push_back(&particle);
        } else {
            particle.w *= updateweight(particle.p, prev);
            alive.push_back(&particle);
            totalWeight += particle.w;
        }
    }
}


void ValuationContinuousSMC::resample(vector<Particle *>& alive, vector<Particle *>& dead) {
    vector<double> randVec = getRandomVector(dead.size(), totalWeight);
    int i = 0;
    double low = 0.0;
    for (Particle *p : alive) {
        low += p->w;
        int copies = 0;
        while (i < randVec.size() && randVec[i] <= low) {
            *dead[i] = *p;
            copies++;
            i++;
        }
        if (copies > 0) {
            double newWeight = p->w / (copies + 1);
            p->w = newWeight;
            for (int j = i - copies; j < i; j++) {
                dead[j]->w = newWeight;
            }
        }
    }
}

vector<double> getRandomVector(int n, double max) {
    vector<double> vec;
    mt19937_64 rng(std::random_device{}());
    uniform_real_distribution<double> dist(0.0, 1.0);
    getRandomHelper(vec, n, max, 0.0, rng, dist);
    return vec;
}


void getRandomHelper(vector<double>& vec, int n, double max, double min,
      mt19937_64& rng, uniform_real_distribution<double>& dist) {
    if (n <= 0) {
        return;
    }

    double val = min + (max - min)*(1.0 - pow(dist(rng), 1.0 / n));
    vec.push_back(val);
    getRandomHelper(vec, n - 1, max, val, rng, dist);
}

double ValuationContinuousSMC::getPayoff() {
    double sum = 0.0;
    for (Particle& particle : particles) {
        sum += max(0.0, exp(particle.p) - K)*particle.w;
    }
    return sum / particles.size();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ValuationContinuousMC::simulatePrice(double& survived) {
    double value = S0;
    vector<double> pricePath = vector<double>(N + 1);
    for (int i = 0; i < N; i++) {
        pricePath[i] = value;
        value += getDrift();     
        if (ifHitBarrier(value)) {
            survived = 0.0;
            return 0.0;
        }
    }
    pricePath[N] = value;
    for (int i = 0; i < N; i++) {
        survived *= updateweight(pricePath[i + 1], pricePath[i]);
    }
    return max(exp(value) - K, 0.0) * survived;
}
