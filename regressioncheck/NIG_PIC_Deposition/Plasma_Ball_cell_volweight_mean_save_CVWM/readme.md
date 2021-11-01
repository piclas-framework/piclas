# plasma ball example
- A ball of particles with, each with charge=1.60217653E-5 is initialized and deposited on the grid
    in order to check the deposition method "cell_volweight_mean"
- two different meshes are compared with a total of 2 elements per mesh
- the fallback in the CVWM algorithm is triggered by particles that are tracked via TriaTracking, but the RefMapping algorithm that
is used for deposition fails to find the particle in the same element as TriaTracking due to the conflicting convex/concave side
because one element interface is strongly deformed
