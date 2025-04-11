# Case 5: MCC with ioniation at high background pressure
This can only considers the ionization process and deactivates all other cross-sections.
Note that the cross-section data for ionization has been altered by a factor of 1/10 in order
to shift the ionization process to lover energies. These cross-sections cannot be used for actual
simulations, where physical results are sought.

Changes in XSec_Database_H_Parodi_Benchmark_only-ionization.h5 as compared to case 4
- H-HIon1
  - Removed BACKSCATTER data completely
  - Changed ELASTIC data to [[1e-4, 0.0],[1e4, 0.0]], effectively setting the elastic cross-section to zero
- H-electron
  - Removed ELECTRONIC cross-section data for the levels 10.20001 and 10.20002 completely
  - Changed ELASTIC data to [[...],[0.0, ..., 0.0]], effectively setting the elastic cross-section to zero

In piclas, the energy after the ionization process is distributed among all product species, which includes ions.
In Parodi, only electrons are assigned the remaining energy after the chemical reaction.
This results in different densities between piclas and VKI/LPP simulation results.

The energy level in Electronic-State-Database-Parodi-Benchmark-ionization-energy.h5 has been changed
so that the resulting enthalpy of reaction matches the value in Parodi, which is 1.43 eV.
Therefore, the dataset has been changed for H to 16594.4703702144 K.

