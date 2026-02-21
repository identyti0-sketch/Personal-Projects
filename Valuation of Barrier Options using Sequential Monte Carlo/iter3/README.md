# Iteration 3

Changes:
  - To better facilitate Browian calculations prices for sample instances and barriers are stored in log form.
  - Added a "Semi Continuous" and a "Continuous" monitoring system:
    - If a sampled price instance doesn't explicitly cross the barrier, the probability of barrier crossing is approximated using the Brownian Bridge correction for single barrier under the assumption that the crossing the upper and lower barrier is independent with the probability of crossing each barrier being given by:
      - P<sub>B</sub> = exp(-2*(log Sₜ - log B)*(log S<sub>t-1</sub> - log B)), where B is the Barrier price and S<sub>t-1</sub> is the previous sampled price, and the probability of crossing both barriers is given by:
      - P<sub>U</sub> + P<sub>L</sub> - P<sub>U</sub>P<sub>L</sub>.
    - For the "Semi Continuous" implementation, the instances which do not explicitly pass the barrier may be resampled based on the probability given by the Brownian Bridge correction.
    - For the "Continuous" implementation, instances are not resampled if the barrier isn't explicitly crossed, instead they have their weights reduced proportional to the probabliity. Should the barrier be explicitly crossed, the instances are resampled based on weight of other instances rather than uniformly.
      The weight of the resampling instance is distributed evenly among the resampled particles.
    - After empirical testing, it was determined that the Continuous implementation was more precise (lower stderr) and more time efficient (shorter execution times).
    - Due to the probability calculation being an approximation under the assumption that the probability of crossing either barrier is independent of each other, it takes multiple steps (N > 10) for the simulation to reach around a 95% accuracy for the given scenario.

 
