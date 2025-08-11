# PIC - Electron beam acceleration
* Inserting electrons with a zero velocity and at a low temperature
* Acceleration with a potential difference of 120 kV
* ATTENTION: two different meshes are tested
  * 3D cylinder ... with === axial direction in z ===
  * 2D axisymmetric with === axial direction in x ===
* Comparing the current at the inlet (should correspond to 2E-3) and the kinetic energy of the exiting particles
  * The total kinetic energy of the leaving particles must be 2.4E-9 Joule (= V I dt AnalyzeStep = (120000eV) (2E-3A) (1E-12s) (100))
