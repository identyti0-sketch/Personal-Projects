#include <ctime>
#include <chrono>
#include "valuation.hpp"

thread_local unique_ptr<Rng> Valuation::rngEngine = nullptr;
thread_local bool Valuation::rngInitialized = false;

Valuation::Valuation(const string& filename) :gen(random_device{}(), 0) {
    ifstream file(filename);
    
    if (!(file >> S0 >> K >> U >> L >> r >> sigmaYearly >> T >> N >> M)) {
        cerr << "Input format error\n";
        return;
    }
    S0 = log(S0);
    U = log(U);
    L = log(L);
    sigmaStep = sigmaYearly * sqrt(T / N);
    seed = random_device{}();
    dist = normal_distribution<double>(0, sigmaStep);
}


void Valuation::repeatExperiment(int n) {
    struct result res;
    
    omp_set_dynamic(0);
    
    auto startTime = chrono::steady_clock::now();
    for (int i = 0; i < n; i++) {
        seed = random_device{}();
        runSimulation(res);
    }
    auto endTime = chrono::steady_clock::now();
    double elapsedTime = chrono::duration<double>(endTime - startTime).count();
    double averagePrice = accumulate(res.price.begin(), res.price.end(), 0.0) / res.price.size();
    double averageSurvivalRate = accumulate(res.survivalRate.begin(), res.survivalRate.end(), 0.0) / res.survivalRate.size();
    double standardErr = stdErr(res, averagePrice);
    double efficiency = standardErr * standardErr * elapsedTime;
    ofstream out("output.csv", ios::app);
    if (!out) {
        cerr << "Error opening file!\n";
        return;
    }    
    out << N << "," << averagePrice << "," << averageSurvivalRate << "," << standardErr << "," << efficiency << endl;
    out.close();
}

void Valuation::updateSteps(int steps) {
    N = steps;
    sigmaStep = sigmaYearly * sqrt(T / N);
    dist = normal_distribution<double>(0, sigmaStep);
    #pragma omp parallel
    {
        if (rngInitialized) {
            rngEngine->dist = normal_distribution<double>(0, sigmaStep);
        }
    }
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
    return (r - 0.5 * sigmaYearly * sigmaYearly) * (T / N);
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


inline void Valuation::initRng() {
    if (!rngInitialized) {
        uint32_t threadNum = static_cast<uint32_t>(omp_get_thread_num());
        rngEngine = make_unique<Rng>(seed, threadNum, sigmaStep);
        rngInitialized = true;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Rng::Rng(uint64_t seed, uint32_t counter, double sigma) : rng(seed, counter) {
    dist = normal_distribution<double>(0, sigma);
}

double Rng::getRand() {
    return dist(rng);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ValuationSMC::runSimulation(struct result& res) {
    particles = vector<double>(M, S0);
    double survivalRate = 1.0;
    for (int i = 0; i < N; i++) {
        if (step(survivalRate)) {
            break;
        }
    }
    double discountFactor = exp(-r * T);
    double payoff = getPayoff();
    res.price.push_back(payoff * discountFactor * survivalRate);
    res.survivalRate.push_back(survivalRate);
}

bool ValuationSMC::step(double& survivalRate) {
    vector<int> alive;
    vector<int> dead;

    increment(alive, dead);    
    survivalRate *= (double)alive.size() / particles.size();
    if (alive.size() == 0) {
        return true;
    } 
    resample(alive, dead);
    return false;
}


void ValuationSMC::increment(vector<int>& alive, vector<int>& dead) {
    for (int i = 0; i < particles.size(); i++) {
        particles[i] += getDrift() + dist(gen);   
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
    #pragma omp parallel reduction(+:averagePrice, survivalCount) 
    {   
        initRng();
        
        #pragma omp for
        for (int i = 0; i < M; i++) {
            double survived = 1.0;
            double price = simulatePrice(survived);
            averagePrice += price;
            survivalCount += survived;
        }
    }

    double discountFactor = exp(-r * T);
    res.price.push_back(averagePrice / M * discountFactor);
    res.survivalRate.push_back(survivalCount / M);
}

double ValuationMC::simulatePrice(double& survived) {
    double value = S0;
    vector<double> pricePath = vector<double>(N);
    for (int i = 0; i < N; i++) {
        value += getDrift() + rngEngine->getRand();
        pricePath[i] = value;
        if (ifHitBarrier(value)) {
            survived = 0.0;
            return 0.0;
        }
    }
    double weight = 1.0;
    if (survived > 0.0) {
        weight *= updateweight(pricePath[0], S0);
        for (int i = 1; i < N; i++) {
            weight *= updateweight(pricePath[i], pricePath[i - 1]);
        }
    }
    return max(exp(value) - K, 0.0) * weight;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////


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

bool ValuationContinuousSMC::step() {
    vector<char> deadMask(particles.size(), false);
    int aliveCount = increment(deadMask);
    if (aliveCount == 0) {
        return true;
    }
    vector<Particle *> alive(aliveCount);
    vector<Particle *> dead(M - aliveCount);
    vector<double> prefixSum(aliveCount);
    unmask(alive, dead, deadMask, prefixSum);
    totalWeight = prefixSum.back();
    resample(alive, dead, prefixSum);
    return false;
}


int ValuationContinuousSMC::increment(vector<char>& deadMask) {
    int aliveCount = 0;
    #pragma omp parallel reduction(+:aliveCount) 
    {
        initRng();
        #pragma omp for
        for (int i = 0; i < particles.size(); i++) {
            Particle& particle = particles[i];
            double prev = particle.p;
            particle.p += getDrift() + rngEngine->getRand();   
            if (ifHitBarrier(particle.p)) {
                particle.w = 0.0;
                deadMask[i] = true;
            } else {
                particle.w *= updateweight(particle.p, prev);
                aliveCount++;
            }
        }
    }
    return aliveCount;
}

void ValuationContinuousSMC::unmask(vector<Particle *>& alive, vector<Particle *>& dead, vector<char>& deadMask, vector<double>& prefixSum) {
    int aliveIndex = 0;
    int deadIndex = 0;
    double cumulativeWeight = 0.0;
    for (int i = 0; i < particles.size(); i++) {
        if (deadMask[i]) {
            dead[deadIndex++] = &particles[i];
        } else {
            alive[aliveIndex] = &particles[i];
            cumulativeWeight += particles[i].w;
            prefixSum[aliveIndex] = cumulativeWeight;
            aliveIndex++;
        }
    }
}

void ValuationContinuousSMC::resample(vector<Particle *>& alive, vector<Particle *>& dead, vector<double>& prefixSum) {
    if (dead.size() == 0 || alive.size() == 0) {
        return;
    }
    vector<double> randVec = getRandomVector(dead.size(), prefixSum.back());
    vector<int> partitionCeilDeadSetup(omp_get_max_threads());
    vector<int> partitionCeilDead(omp_get_max_threads());
    vector<int> partitionCeilAlive(omp_get_max_threads());
    int baseSize = dead.size() / partitionCeilDead.size();
    int remainder = dead.size() % partitionCeilDead.size();
    #pragma omp parallel
    {   
        initRng();
        int i = omp_get_thread_num();
        int ceil = (i + 1) * baseSize + min(i + 1, remainder) - 1;
        if ( i ==  partitionCeilDead.size() - 1) {
            ceil = dead.size() - 1;
            partitionCeilAlive[i] = prefixSum.size() - 1;
        } else {
            partitionCeilAlive[i] = binSearchSmallestAbove(randVec[ceil], prefixSum, 0, prefixSum.size() - 1);
        }
        partitionCeilDeadSetup[i] = ceil;
        #pragma omp barrier
        int upper = ceil;
        int j = i;
        while ( j < partitionCeilDeadSetup.size() - 1 && upper <= partitionCeilAlive[j]) {
            j++;
            upper = partitionCeilDeadSetup[j];
        }
        partitionCeilDead[i] = binSearchGreatestBelow(prefixSum[partitionCeilAlive[i]], randVec, ceil, upper);
        #pragma omp barrier
        int startEnd[4] = {(i == 0) ? -1 : partitionCeilAlive[i - 1], partitionCeilAlive[i], i == 0 ? -1 : partitionCeilDead[i - 1], partitionCeilDead[i]};
        partitionResample(alive, dead, randVec, prefixSum, startEnd);
    }
}

void ValuationContinuousSMC::partitionResample(vector<Particle *>& alive, vector<Particle *>& dead, vector<double>& randVect, vector<double>& prefixSum, int startEnd[4]) {
    int deadIndex = startEnd[2] + 1;
    for (int i = startEnd[0] + 1; i <= startEnd[1]; i++) {
        int copies = 0;
        while (deadIndex <= startEnd[3] && randVect[deadIndex] <= prefixSum[i]) {
            *dead[deadIndex] = *alive[i];
            copies++;
            deadIndex++;
        }

        if (copies > 0) {
            double newWeight = alive[i]->w / (copies + 1);
            alive[i]->w = newWeight;
            for (int j = deadIndex - copies; j < deadIndex; j++) {
                dead[j]->w = newWeight;
            }
        }
    }
}


vector<double> getRandomVector(int n, double max) {
    vector<double> vec(n);
    openrand::Philox rng(random_device{}(), 0);
    uniform_real_distribution<double> dist(0.0, 1.0);
    getRandomHelper(vec, n, max, 0.0, rng, dist);
    return vec;
}


void getRandomHelper(vector<double>& vec, int n, double max, double min,
      openrand::Philox& rng, uniform_real_distribution<double>& dist) {
    if (n <= 0) {
        return;
    }

    int i = 0;
    while (n > 0) {
        double val = min + (max - min)*(1.0 - pow(dist(rng), 1.0 / n));
        vec[i] = val;
        i++;
        n--;
        min = val;
    }
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
        value += getDrift() + rngEngine->getRand();     
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int binSearchSmallestAbove(double target, vector<double>& arr, int low, int high) {
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (arr[mid] < target) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

int binSearchGreatestBelow(double target, vector<double>& arr, int low, int high) {
    while (low < high) {
        int mid = low + (high - low) / 2;
        if (arr[mid] > target)
            high = mid - 1;
        else
            low = mid + 1;
    }
    return high;
}

/*
void ValuationContinuousSMC::runSimulation(struct result& res) {
    particles = vector<Particle>(M, Particle{S0, 1.0});
    totalWeight = M;
    int aliveCount = 0;
    vector<char> deadMask(particles.size(), false);
    #pragma omp parallel
    for (int i = 0; i < N; i++) {
        initRng();
        
        #pragma omp for reduction(+:aliveCount)
        for (int i = 0; i < particles.size(); i++) {
            Particle& particle = particles[i];
            double prev = particle.p;
            particle.p += getDrift() + rngEngine->getRand();   
            if (ifHitBarrier(particle.p)) {
                particle.w = 0.0;
                deadMask[i] = true;
            } else {
                particle.w *= updateweight(particle.p, prev);
                aliveCount++;
            }
        }

        #pragma omp single
        {
            vector<Particle *> alive(aliveCount);
            vector<Particle *> dead(M - aliveCount);
            vector<double> prefixSum(aliveCount);
            unmask(alive, dead, deadMask, prefixSum);
            totalWeight = prefixSum.back();
            resample(alive, dead, prefixSum);
            aliveCount = 0;
        }

        #pragma omp for
        for (int i = 0; i < deadMask.size(); i++) {
            deadMask[i] = false;
        }
        #pragma omp barrier
    }
    double discountFactor = exp(-r * T);
    double payoff = getPayoff();
    res.price.push_back(payoff * discountFactor);

    res.survivalRate.push_back(totalWeight / M);
}

*/