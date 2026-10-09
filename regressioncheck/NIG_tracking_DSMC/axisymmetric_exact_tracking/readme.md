# Exact axisymmetric tracking
* Isothermal cylinder test case with initial particle insertion and surface flux at zero velocity
* Temperature should be constant and velocity remain zero
* At large timesteps, the old axisymmetric tracking would produce a temperature hotspot / velocity above zero
* Comparing the sampled velocity in y-direction, result without exact tracking shows a distinct negative velocity in y-direction in the top cell layer
