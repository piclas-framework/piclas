# DSMC - Macroscopic restart
* Test case based on Flow_N2_70degCone
* Restart the simulation with a different weighting factor while using the sampled macroscopic result to insert particles
* Cell-local weighting: do not change the global weighting factor but determine a cell-local weighting based on DSMCState result
* Enable variable/local time step to check combination with macroscopic restart
* Comparing the total number of particles in the simulation domain and the number density in the DSMCState