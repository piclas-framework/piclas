# Particle Tracking

Three different particle tracking methods are implemented in PICLas and are selected via

    TrackingMethod = triatracking ! Define Method that is used for tracking of
                                  ! particles:
                                  ! refmapping (1): reference mapping of particle
                                  ! position with (bi-)linear and bezier (curved)
                                  ! description of sides.
                                  ! tracing (2): tracing of particle path with
                                  ! (bi-)linear and bezier (curved) description of
                                  ! sides.
                                  ! triatracking (3): tracing of particle path with
                                  ! triangle-aproximation of (bi-)linear sides.

For conventional computations on (bi-, tri-) linear meshes, the following tracking algorithm is recommended:

    TrackingMethod = triatracking

Following options are available to get more information about the tracking, e.g. number of lost particles:

| Option               | Values | Notes                                                                                                                    |
| :------------------- | :----: | :----------------------------------------------------------------------------------------------------------------------- |
| DisplayLostParticles |  F/T   | Display position, velocity, species and host element of                                                                  |
|                      |        | particles lost during particle tracking (TrackingMethod = triatracking, tracing) in the std.out                          |
| CountNbrOfLostParts  |  T/F   | Count number of lost particles due to tolerance issues.                                                                  |
|                      |        | This number is a global number, summed over the full simulation duration and includes particles lost during the restart. |
|                      |        | The lost particles are output in a separate `*_PartStateLost*.h5` file.                                                  |

The two alternative tracking routines and their options are described in the following.

## DoRefMapping

    TrackingMethod = refmapping

This method is the slowest implemented method for linear grids and large particle displacements.
A particle is mapped into a element to compute the particle position in the reference space.
This test determines in which element a particle is located. Each element has a slightly larger
reference space due to tolerance. Starting from reference values >=1. the best element is found and used for the
hosting element. In order to take boundary interactions into account, all BC faces in the halo vicinity of the element
are checked for boundary interactions and a boundary condition is performed accordingly. This algorithm has an
inherent self-check. If a boundary condition is not detected, the particle position is located outside of all elements.
A fall-back algorithm is then used to recompute the position and boundary interaction. Periodic domains are only possible
for Cartesian meshes. The particle position is used for periodic displacements.

|      Option       | Values |                          Notes                          |
| :---------------: | :----: | :-----------------------------------------------------: |
| CartesianPeriodic |  T/F   | If a fully periodic box (all 6 sides) is computed, the  |
|                   |        | intersections do not have to be computed. Instead, each |
|                   |        | particle can be simply shifted by the periodic vector.  |
|   FastPeriodic    |  T/F   | Moves particle the whole periodic distance once, which  |
|                   |        |  can be several times the mesh size in this direction.  |


## Tracing

    TrackingMethod = tracing

This method traces the particle trajectory throughout the domain. The initial element is determined by computing the intersection
between the particle-element-origin vector and each element face. If none of the six element faces are hit, the particle is
located inside of this element. Next, the particle trajectory is traced throughout the domain. Hence, each face is checked
for an intersection and a particle assigned accordingly to neighbor elements or the interaction with boundary conditions occur. This
algorithm has no inherent self-consistency check. For critical intersections (beginning or end of a particle path or if a particle is located close to
the edges of element faces) an additional safety check is performed by recomputing the element check and if it fails a re-localization of
the particle is required. Particles traveling parallel to element faces are in an undefined state and are currently removed from the computation.
This leads to a warning message.

## Parameters for DoRefMapping and Tracing  (NEEDS UPDATING)

Following parameters can be used for both schemes.

|        Option         |  Values  |                          Notes                           |
| :-------------------: | :------: | :------------------------------------------------------: |
|   MeasureTrackTime    |   T/F    |  Measure the time required for tracking and init local.  |
|    RefMappingGuess    |   1-4    |   Prediction of particle position in reference space:    |
|                       |    1     |       Assumption of a linear element coord system.       |
|                       |    2     |      Gauss point which is closest to the particle.       |
|                       |    3     |        CL point which is closest to the particle.        |
|                       |    4     |              Trivial guess: element origin               |
|     RefMappingEps     |   1e-4   |  Tolerance of the Newton algorithm for mapping in ref.   |
|                       |          |  space. It is the L2 norm of the delta Xi in ref space.  |
|    BezierElevation    |   0-50   |   Increase polynomial degree of BezierControlPoints to   |
|                       |          |     construct a tighter bounding box for each side.      |
|     BezierSampleN     |   NGeo   |  Polynomial degree to sample sides for SurfaceFlux and   |
|                       |          |              Sampling of DSMC surface data.              |
|   BezierNewtonAngle   | $<PI/2$  | Angle to switch between Clipping and a Newton algorithm. |
|  BezierClipTolerance  |   1e-8   |      Tolerance of Bezier-Clipping and Bezier-Newton      |
|     BezierClipHit     |   1e-6   | Tolerance to increase sides and path during Bezier-Algo. |
|   BezierSplitLimit    |   0.6    |    Minimum degrees of side during clipping. A larger     |
|                       |          | surface is spit in two to increase convergence rate and  |
|                       |          |              predict several intersections.              |
| BezierClipMaxIntersec |  2*NGeo  |      Maximum number of roots for curvilinear faces.      |
|      epsilontol       | 100*epsM |       Tolerance for linear and bilinear algorithm.       |
