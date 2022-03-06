# plasma ball with save CVWM deposition
- A ball of particles with, each with charge=1.60217653E-5 is initialized and deposited on the grid
    in order to check the deposition method "cell_volweight_mean"
- two different meshes are simulated, one Cartesian and one heavily deformed mesh (with a total of 2 elements per mesh)
- the fallback in the CVWM algorithm is triggered by particles that are tracked via TriaTracking, but the RefMapping algorithm that
is used for deposition fails to find the particle in the same element as TriaTracking due to the conflicting convex/concave side
because one element interface is strongly deformed
- deposition error calculation using NAnalyze=12 (polynomial degree is N=1, but polynomial for error is calculated with N=12)
  - this feature is useful if the Jacobians are determined good enough with N and the integration error is then reduced using a much
  higher NAnalyze (in this case the Jacobians are bad from the beginning, which is why the calculation of the error does not improve
  much)
- also considers superB magnetic field to test functionality of CVWM fallback in combination with B-field
