# FPFlow Reservoir
* Test the relaxation from thermal non-equilibrium with ESFP (CollModel=2), with its exact (-ESFP-Model=1) and approximative method (-ESFP-Model=2)
* Not testing the cubic model, since its only suitable for atomic species in the current state
* Tests also the octree cell refinement and continuous treatment of vibrational energy
* T_trans = 10000K, T_rot = 7500K, T_vib = 5000K
* Relevant parameter
  * Particles-FP-CollModel = 2
  * Particles-ESFP-Model=1,2
  * Particles-FP-DoVibRelaxation=true
  * Particles-FP-UseQuantVibEn=true,false
  * Particles-FP-DoCellAdaptation=F,T
  * Particles-FP-MinPartsPerCell=20
