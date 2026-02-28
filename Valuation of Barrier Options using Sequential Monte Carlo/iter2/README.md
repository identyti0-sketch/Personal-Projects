# Iteration 2

Changes:
 - Refactored code using OOP to improve readability.
 - Reads input parameters off example1.txt. Reads parameters in follwoing order:
    1. Initial Price
    2. Strike Price
    3. Upper Knock-out Barrier
    4. Lower Knock-out Barrier
    5. Risk-Neutral Drift
    6. Volatility
    7. Time to Maturity
    8. Number of steps
    9. Number of samples of asset prices
  - Added a standard Monte Carlo simulation where samples are propogated individually instead of together such that when a sample price crosses a barrier its entire simulation path is discarded rather than being resampled.
