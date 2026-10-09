# 1D Plasma wave
- Testing different shape function deposition methods
- 1D plasma wave test case (1D periodic) with N=5 polynomial degree
- Desired result is a sinus function of the electric potential
- The electric energy is integrated (FieldAnalyze.csv) over time and compared with a reference value (origin unkown)
- Particles have been inserted into the restart file via SpaceIC = sin_deviation
- For this quick test only 50 ions 
- HDG solver was originally set to epcCG=1e-6 but is now set to 1.0 as it does not converge else wise (the system is fully periodic)
