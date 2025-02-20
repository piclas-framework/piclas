# DSMC - Macroscopic restart
* Test case based on Flow_N2_70degCone
* Restart the simulation with a different weighting factor while using the sampled macroscopic result to insert particles
* Comparing only the total number of particles in the simulation domain as a comparison of the DSMCState is not possible due to strong statistical fluctuations
* Cell-local weighting: do not change the global weighting factor but determine a cell-local weighting based on DSMCState result