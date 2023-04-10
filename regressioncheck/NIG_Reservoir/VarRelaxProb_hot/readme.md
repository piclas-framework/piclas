# Relaxation of a hot reservoir of N2 and O2 with variable relaxation probabilities 
* Test of thermal relaxation 
* Settings
  * PartDensity           = 2*1E20
  * MacroParticleFactor   = 5E9
  * T_trans               = 20000
  * T_Rot_N2              = 35000
  * T_Vib_N2              = 5000
  * T_Rot_O2              = 15000
  * T_Vib_O2              = 25000
  * initial autorestart at the middle of the simulation time.
  * First, translational temperature rises, due to the higher rotational temperature and its faster relaxation. After total relaxation of the rotation, the translational temperature decreases, due to slower relaxation of the vibrational degrees of freedom. The position and and temperature of the plateau of translational temperature is very sensitive to the relaxation probabilities. The initial autorestart shows the capability to save and load the average vibrational relaxation probability to and from the state file, respectively. Dissociation is neglected, due to short simulation time.
