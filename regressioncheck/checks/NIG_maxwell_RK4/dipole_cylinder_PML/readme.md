# 3D Electromagnetic dipole with PMLs
- Tests Maxwell only simulation in combination with PMLs
- cylindrical domain with curved elements
- compares
  - L2 error (indirectly checks the correct mesh volume)
  - volume data for PML layers (Dipole_PMLZetaGlobal_000.00000000000000000.h5) that is written to h5 at the beginning of the simulation
