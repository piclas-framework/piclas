# 2D: vmpf initial particle insertion
* Testing the initial particle inserting routine for 2D, and axisymmetric (radial weighting=T/F)
* Inserts particles in different spaceICs: planar - cuboid and cylinder; axisymmetric - cylinder and sphere
* Half of the tests is in the "Symmetry_Initial_Particle_Emission" check due to too many crosscombinations
* Comparison of PartAnalyze: number density and temperature
* Desired result of DSMC output: constant number density of at n=1E16 in the choosen spaceIC (integral density is lower due to empty regions in the computational domain) and all temperatures at T=1000K
* Comparison of DSMCState not possible due to strong fluctuations
