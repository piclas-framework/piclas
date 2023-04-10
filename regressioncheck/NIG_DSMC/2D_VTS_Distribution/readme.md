# 2D axisymmetric: Variable Time Step - Distribution
* Testing the calculation of a variable time step distribution, based on a read-in DSMCState file
* Utilized criteria is the number of particles and the time step is adapted to achieve approximately 20 particles per cell
* Parameters:
  * Macroscopic restart should be performed to immediately insert the appropriate number of particles:
    * Particles-MacroscopicRestart = T
    * Particles-MacroscopicRestart-Filename = 2D_CellLocal_Insert_DSMCHOState_000.00000000010000000.h5
  * Variable Time Step
    * Part-VariableTimeStep-Distribution = T
    * Part-VariableTimeStep-Distribution-Adapt = T
    * Part-VariableTimeStep-Distribution-MaxFactor = 1.0
    * Part-VariableTimeStep-Distribution-MinFactor = 0.1
    * Part-VariableTimeStep-Distribution-MinPartNum = 20
