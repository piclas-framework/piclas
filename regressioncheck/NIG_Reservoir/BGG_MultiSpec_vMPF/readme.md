# Multi-species background gas with VHS and basic vMPF
* Testing the multi-species background gas with a CO2-N2-He reservoir using the Variable Hard Sphere model
* Background gas: 6E23 1/m3 (1/3 each), Electrons: n_e = 1E21 1/m3
* Testing basic variable weighting factor functionality, load balance, and the background gas trace species split feature
  * Part-vMPF = T with equal species MPF, constant number of particles
  * Part-vMPF = T with a different particle MPF compared to background gas, half the number of particles (weighting factor of background gas is irrelevant)
  * Part-vMPF = T with a trace background species, where particles are split to test more often with the trace background species, number of particles is increasing as electrons with lower weighting factor are introduced after collisions with the trace species
* Comparison of collision rate with Part-vMPF = F reference