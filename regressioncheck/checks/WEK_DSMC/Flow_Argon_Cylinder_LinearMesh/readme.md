# DSMC - Hypersonic flow around a cylinder (2D)
* Simulation of a hypersonic Argon flow around a cylinder at a Knudsen number of 0.25
* Test case based on Lofthouse, A. J., Boyd, I. D., & Wright, M. J. (2007). Effects of continuum breakdown on hypersonic aerothermodynamics. Physics of Fluids, 19(2). https://doi.org/10.1063/1.2710289
* 2D simulation of the half of simulation domain: Linear mesh with surface flux particle emission, only front part of the cylinder modelled to avoid the wake (otherwise the comparison of the DSMC surface state is prevented by high statistical fluctuations)
* Comparison of the heat flux and force per area with a reference surface state file, which was compared to the publication mentioned above
  * Maximum heat flux (Deviation < 4%):
    * PICLas: q = 5760 W/m²
    * Publication: q = 5984 W/m²
  * Drag (Deviation = 2%, determined with VisIt: Sum of ForcePerArea[0]*area(mesh)/0.002 [0.002 = Length in z-direction])
    * PICLas: 2.04 N/m (half model, drag multiplied by 2, no wake)
    * Publication: 2.08 N/m (full model)