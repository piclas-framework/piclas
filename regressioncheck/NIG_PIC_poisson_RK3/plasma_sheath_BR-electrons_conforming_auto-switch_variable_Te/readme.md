# Plasma sheath with Boltzmann relation electrons
- automatic switch between BR and fully kinetic
- multiple-switch (kin->BR->kin->BR->kin....) continued for ever
- temperature for Te is variable 
  - from initially 8.617332E-02 eV to 8.186461742254e-2, i.e., from 1000 K to 950 K just for testing the functionality
  - this model the electron temperature relaxation towards the BGGas properties over time
- automatic load balance after the first time step for MPI>1
