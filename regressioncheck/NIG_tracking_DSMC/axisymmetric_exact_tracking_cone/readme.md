# Exact axisymmetric tracking on a conical wall
* Conical domain (45 deg outer wall, radius grows from 0.5 mm at x- to 1.5 mm at x+) with a directed particle beam
* Three runs (VeloVecIC sweep, velocityDistribution=constant) exercise the degenerate geometries in ParticleThroughSideCheck2DRotSym:
  * (/1,0,0/): purely axial beam, regular (non-degenerate) reference case
  * (/0,1,0/): purely radial beam (no axial velocity, sx -> 0); checking a constant-x side hits case 3, where S = (xNode1-x_pos_start)/sx is a division by zero
  * (/1,1,0/): beam exactly parallel to the cone wall (dr/dx = vy/vx = 1); the intersection quadratic collapses (leading coefficient a -> 0, trajectory parallel to the cone surface)
* Without the guards against sx -> 0 (case 3) and a -> 0 (cases 1/2), these intersections divide by zero, producing NaN/Inf and aborting the run in debug mode
* The test passes if all runs complete without aborting, i.e. the degenerate cases are handled correctly
