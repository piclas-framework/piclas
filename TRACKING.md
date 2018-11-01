## PICLas Tracking Routines

PICLas uses two different tracking routines. The switch between both is

		DoRefMapping=T/F

Only particle located inside of the elments are emitted, hence particle located at boundary faces are not taken 
into account for the primary emission, neglecting the SurfaceFlux.

### DoRefMapping

This method is the slowest implemented method for linear grids and large particle movements. A particle is mapped into 
a element to compute the particle position
in the reference space. This test determines in which element a particle is located. Each element has a slightly larger
reference space due to tolerance. Starting from reference values >=1. the best element is found and used for the 
hosting element. In order to take boundary interactions into account, all BC faces in the halo vicinity of the element
are checked for boundary interactions and a boundary condition is performed accordingly. This algorithm has a 
inherent self check. If a boundary condition is not detected, the particle position is located outside of all elements.
A fall-back algorithm is used to recompute the position and boundary interaction. Periodic domains are only possible
for Cartesian meshes. The particle position is used for periodic displacements.

| Option                 | Values     |  Notes                                                  |
|:----------------------:|:----------:|:-------------------------------------------------------:|
| CartesianPeriodic      | T/F        | If a fully periodic box (all 6 sides) is computed, the  |
|                        |            | intersections do not have to be computed. Instead, each |
|                        |            | particle can be simply shifted by the periodic vector.  |
| FastPeriodic           | T/F        | Moves particle the whole periodic distance once, which  |
|                        |            | can be several times the mesh size in this direction.   |


### Tracing

This method traces the particles throughout the domain. The initial element is determined by computing the intersection
between the particle-element-origin vector and each element face. If non of the six element faces are hit, the particle is 
located inside of this element. Next, the particle trajectory is traced throughout the domain. Hence, each face is checked
for an intersection and a particle mapped accordingly into the neighbor element or perform a boundary condition. This 
algorithm has no inherent self-consistency check. For critical intersections (beginning,end of particle path or close to 
edges of faces) an additional safety check is performed by recomputing the element check and if it fails a re-localization of 
the particle. Particles traveling parallel to faces are in a undefined state and a currently removed. This prints a warning
message. Note, the tracing on periodic meshes works only for non-mpi computations. Periodic displacement requires 
additional coding.


| Option                 | Values     |  Notes                                                  |
|:----------------------:|:----------:|:-------------------------------------------------------:|
| CountNbOfLostParts     | T/F        | Count number of lost particles due to tolerance issues. |
|                        |            | This number is a global number, summed over the full t. |



### Options Available For Both Methods

Following parameters can be used for both schemes.

| Option                 | Values     |  Notes                                                  | 
|:----------------------:|:----------:|:-------------------------------------------------------:|
| MeasureTrackTime       | T/F        | Measure the time required for tracking and init local.  |
| RefMappingGuess        | 1-4        | Prediction of particle position in reference space:     |
|                        | 1          | Assumption of a linear element coord system.            |
|                        | 2          | Gauss point which is closest to the particle.           |
|                        | 3          | CL point which is closest to the particle.              |
|                        | 4          | Trival guess: element origin                            |
| RefMappingEps          | 1e-4       | Tolerance of the Newton algorithm for mapping in ref.   |
|                        |            | space. It is the L2 norm of the delta Xi in ref space.  |
| BezierElevation        | 0-50       | Increase polinomial degree of BezierControlPoints to    |
|                        |            | construct a thighter bounding box for each side.        |
| BezierSampleN          | NGeo       | Polynomial degree to sample sides for SurfaceFlux and   |
|                        |            | Sampling of DSMC surface data.                          |
| BezierNewtonAngle      | <PI/2      | Angle to switch between Clipping and a Newton algorithm.|
| BezierClipTolerance    | 1e-8       | Tolerance of Bezier-Clipping and Bezier-Newton          |
| BezierClipHit          | 1e-6       | Tolerance to increase sides and path during Bezier-Algo.|
| BezierSplitLimit       | 0.6        | Minimum degrees of side during clipping. A larger       |
|                        |            | surface is spit in two to increase convergence rate and |
|                        |            | predict several intersections.                          |
| BezierClipMaxIntersec  | 2*NGeo     | Maximum number of roots for curvilinear faces.          |
| epsilontol             | 100*epsM   | Tolerance for linear and bilinear algorithm.            |

### Possible outdated

| Option                 | Values     |  Notes                                                  | 
|:----------------------:|:----------:|:-------------------------------------------------------:|
| BezierEpsilonBilinear  | T/F        | Tolerance for linear-bilinear side. Obsolet.            |
