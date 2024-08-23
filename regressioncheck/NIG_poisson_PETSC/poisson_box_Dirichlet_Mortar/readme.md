# Pure Poisson Solver
* Test the restart functionality of the HDG solver with different number of MPI processes
    * the number of iterations should be zero but can vary depending on
      * the machine precision
      * the number of MPI processes
      * the compiler version
* Dirichlet boundary conditions left/right
* Periodic boundary conditions in the other two directions
* 13 elements and one Mortar interface 1:4
* Required compiled executable with
  * PICLAS_CODE_ANALYZE = ON
  * PICLAS_EQNSYSNAME = poisson
* This test works as follows
  * Single execution from t=0 to create a state file at t=5ns
  * Restart from this state file with a different number of procs (within the post-piclas-restart directory) as defined in externals.ini
  * Each run stores its output in FieldAnalyze.csv (with the same time=0.5000000000000000E-008) that is compared in the end with HDGIterations.csv
