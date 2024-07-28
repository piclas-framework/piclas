# PIC-axisymmetric HDG - O'Connor Benchmark Expanding electron beam
* A beam of electrons is injected into a cylindrical drift tube and evolved over time
* Taken from S.O'Connor et al., A Set of Benchmark Tests for Validation of 3D Particle In Cell Methods (https://doi.org/10.48550/arXiv.2101.09299)
* Beam is expanding over time due to the mutual repulsion (analytical values for beam radius over time are given)
* Axisymmetric depositions of HDG are tested using
* Three different depositions are tested with PETSc and a direct solver (PrecondType = 10)
  * cell_mean
  * cell_volweight
  * cell_volweight_mean
* Comparing the kinetic (PartAnalyze.csv) and potential energy (FieldAnalyze.csv) of the electrons over time
