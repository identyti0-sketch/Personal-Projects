# Iteration 1: Basic Proof of concept 
Implemented a basic Sequential Monte Carlo (SMC) simulation with discrete monitoring, where the the parameters are set to the scenario described in chapter 6.

The process is as given:

1. The sampled prices are all propogated forward in time. These prices follow geometric Browian Motion:
  - d(log Sₜ) = (r − ½σ²) dt + σ dWₜ
  - Where Sₜ is the price of the sample asset, r is the risk-neutral drift (interest - dividend), σ is the volatility, and dWₜ is a random number generated from a normal distribution.
2. If a sampled price crosses the given barrier, the instance of asset prices is resampled in a uniform distribution from instances which have not crossed the barrier.
3. The final value of the option price is given by discount-factor*accumlative-survival-rate*average-value, where the discount factor is given by exp(-r*T), where T is the length of the option's maturity.

In the implemented SMC the return values roughly match the average value, surival rate and stderr given in Table 2, page 26 of the paper, though the it appears that my implementation of SMC is less efficient relative to the standard monte carlo implementation.
