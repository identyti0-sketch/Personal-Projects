# Valuation of Barrier Options using Sequential Monte Carlo

This repository contains a C++ implementation of barrier option pricing using **Sequential Monte Carlo (SMC)** methods. The project replicates and optimises techniques described in:

Shevchenko, P. V., & Del Moral, P. (2014). *Valuation of Barrier Options using Sequential Monte Carlo*.

The implementation simulates asset price paths under risk-neutral geometric Brownian motion and improves accuracy over standard Monte Carlo using Brownian bridge correction and asset instance resampling.

---

## Overview

Barrier options are path-dependent derivatives whose payoff depends on whether the underlying asset crosses specified barrier levels. Accurate pricing under continuous monitoring is challenging with standard Monte Carlo methods due to missed barrier crossings between discrete time steps.

This project:

- Implements Sequential Monte Carlo (SMC)
- Incorporates Brownian bridge correction for continuous monitoring
- Uses particle resampling to reduce estimator variance
- Studies convergence behaviour under different discretisations
- Explores performance optimisations

### Requirements
- C++17 compatible compiler
