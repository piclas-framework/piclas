# Pure Poisson Solver
* Test the restart functionality of the HDG solver with different number of MPI processes
    * the number of iterations should be zero but can vary depending on the machine precision and number of MPI processes
* Dirichlet boundary conditions left/right
* Periodic boundary conditions in the other two directions
* 13 elements and one Mortar interface 1:4
* required compiled executable with PICLAS_CODE_ANALYZE = ON and PICLAS_EQNSYSNAME = poisson to force FieldAnalyze.ini output during restart
