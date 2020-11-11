# VHS/VSS Self-Diffusion of N2
* Testing the collision models, Variable Hard Sphere (VHS) and Variable Soft Sphere (VSS), comparison with the analytical solution of Fick's law
* Channel with a square cross-section (dz=10nm) and a length of L = 4000nm in x
  * Left side at x=0 is an open inlet boundary condition at 80%/20% N2 (Species1/Species2)
  * Right side at x=4000nm is a perfectly reflecting wall
* Test case is based on J. A. H. Dreyer et al., “Simulation of gas diffusion in highly porous nanostructures by direct simulation Monte Carlo,” Chem. Eng. Sci., vol. 105, pp. 69–76, 2014.
* VHS-VSS-SelfDiffusion.ods shows spatial distribution at t=40ns (produced with MPF=1) and the temporal evolution integrated number density in the channel (corresponds to the number density shown in PartAnalyze.csv)
  * Simulation result of the spatial distribution is normalized with the number density at x=0 to avoid influence due to inlet
* Produced DSMCState result can be used to compare the spatial distribution at t=40ns (however, not suitable for the regression test due to statistical noise)