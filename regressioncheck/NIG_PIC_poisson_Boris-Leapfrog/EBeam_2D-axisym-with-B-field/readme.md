# PIC - Electron beam acceleration
* Inserting electrons with a zero velocity and at a low temperature
* Acceleration with a potential difference of 120 kV
* 2D axisymmetric with === axial direction in x ===
* 2D magnetic field is read from reggie-rot-symmetry-B-field.h5 (Br, Bz, r and z) and is used for the axisymmetric setting with the
  axial direction in x and the radial direction in y (z is the symmetry dimension)
* Comparing the current at the inlet (should correspond to 2E-3) and the kinetic energy of the exiting particles
  * The total kinetic energy of the leaving particles must be 2.4E-9 Joule (= V I dt AnalyzeStep = (120000eV) (2E-3A) (1E-12s) (100))
