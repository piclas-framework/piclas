Setup:
  - simple setup that solves Maxwell's equations for code coverage (no physical results are created)

IMPORTANT NOTICE: 
  - This setup requires PICLAS_READIN_CONSTANTS = ON (see ../builds.ini) because
    otherwise time step will be extremely small and the simulation will run forever!
  - The setup considers (see parameter.ini)
    c0  = 1.0
    eps = 1.0
    mu  = 1.0
