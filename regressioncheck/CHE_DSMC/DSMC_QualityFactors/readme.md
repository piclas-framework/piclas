# DSMC: Quality Factors
* Testing the calculation of the mean/maximum collision probability as well as the mean collision separation distrance over mean free path
* Comparing with the theoretical values from gas kinetic theory
  * Mean collision probability is expected to be 0.06105 (chosen time step divided by the mean collision time, which is mean free path over average thermal velocity)
  * Mean free path is expected to be 2.98E-04 meters (Variable hard sphere model)
  * Mean collision separation distance over mean free path (only in DSMCState) is 0.0001199, compared to the approximated value of 0.0001558, which assumes that the particles are distributed equidistantly in the cell. The difference is due to a 20% lower mean collision separation distance compared to the equidistant distribution (most likely due to nearest neighbour routine)

