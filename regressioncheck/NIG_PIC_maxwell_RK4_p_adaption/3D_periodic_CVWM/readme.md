# Fully 3D periodic with CVWM deposition
- 2 species with positive charge and 2000 particles in total
- cell-volweight mean (CVWM) deposition across periodic boundaries (fully 3D periodic)
- the mesh is built on-the-fly using hopr and ./pre-hopr/hopr.ini
- create restart + reference file in ./pre-piclas/parameter.ini with MPI=1 (single-core run)
  - restart from /pre-piclas/plasma_wave_State_000.00000000000000000.h5 using 1,2,...,30 cores
  - compare output with /pre-piclas/plasma_wave_State_000.00000000010000000.h5: this does not work for random N between 1 and 5
- Test PIC temporal averaging for sampling the field data on Nloc

      CalcTimeAverage = T
      VarNameAvg      = ElectricFieldX
      VarNameAvg      = ElectricFieldY
      VarNameAvg      = ElectricFieldZ
      VarNameAvg      = ElectricFieldMagnitude
      VarNameAvg      = ChargeDensity-Spec01
      VarNameFluc     = ElectricFieldX
      VarNameFluc     = ElectricFieldY
      VarNameFluc     = ElectricFieldZ
      VarNameFluc     = ElectricFieldMagnitude
      VarNameFluc     = ChargeDensity-Spec01
