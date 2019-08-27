# BGK Reservoir - Relaxation of N2
* Test the relaxation from thermal non-equilibrium with ESBGK (CollModel=1) and SBGK (CollModel=2)
* Tests also the octree cell refinement and continuous treatment of vibrational energy
* T_trans = 10000K, T_rot = 7500K, T_vib = 5000K
* Relevant parameter
  * Particles-BGK-CollModel = 1,2
  * Particles-SBGK-EnergyConsMethod=1 (only relevant for CollModel=2)
  * Particles-BGK-DoVibRelaxation=true
  * Particles-BGK-UseQuantVibEn=true,false
  * Particles-BGK-DoCellAdaptation=F,T
  * Particles-BGK-MinPartsPerCell=20
