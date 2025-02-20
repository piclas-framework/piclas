# BGK Reservoir - Relaxation of a Ar-He mixture
* Test the relaxation of a mixture from thermal non-equilibrium with ESBGK (CollModel=1)
* Equal parts of Argon and Helium at different initial temperatures, T_Ar = 10000K and T_He=1000K
* Comparison with DSMC is shown in Figure_Comparison_DSMC.png (dashed line: ESBGK)
  * Time is normalized with the relaxation time (time at which the temperature relaxed to 1/e, T=6471.5 K for Argon, interpolated between time steps)
    * tau_c DSMC: 3.88548E-08 s
    * tau_c ESBGK: 1.23369E-08 s