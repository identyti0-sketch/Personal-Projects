# Iteration 5

Changes:
- Added parallel processing for the standard Monte Carlo stimulation using OpenMP.
- Restructured Continuous Sequential Monte Carlo:
  - Incrementing phase for the sampled instances has been made parallel, to reduce the load of the heavy brownian bridge calculations. Samples are mapped onto a mask based on whether they directly cross the barrier instead of being pushed into a vector to avoid race conditions.
  -  Sequentially, the mask is converted into two vectors which represent samples which crossed (dead) or didn't cross the barrier (alive) in the previous step. In preparation to the resamping step, a prefixsum vector, representing the cumulative distribution of the weights of samples which did not cross the barrier and a ordered uniformly distributed vector of cumulative weight values with a size equal to the number of "dead" samples are created. During the resampling process the values of these vectors will be matched to proportionally distribute the alive sample values to dead samples to resample them.
  -   The prefix sum  and ordered uniformly distributed vectors are partitioned such that there are roughly equal numbers of elements in each partition and that each thread will resample one and only one partition.
- Due to the sequential nature of SMC, where certain steps had to be done in sequence, and at certain points in the code the threads had to be synced up, parallelisation was more complex to implement and less efficient for SMC than it was for standard MC.
