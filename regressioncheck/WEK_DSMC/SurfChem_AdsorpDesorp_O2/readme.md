# DSMC - Adsorption/Desorption of O2 on a catalytic wall
* Free-molecular gap (25 nm) between an inflow plane (BC_Xminus, surface flux) and a catalytic wall (BC_Xplus), symmetric in y/z
* The only case with `Part-nSurfSample = 2` (4 sides x 4 sub-surfaces): the surface coverage is indexed per sub-surface, so this guards the sub-surface mapping
* Literature reaction parameters for O2 on palladium; wall temperature (476.0 K) and the inflow density are chosen, so that desorption is resolved within the simulated 10 ns
* Adsorption and desorption are second order (`DissOrder = 2`): one O2 needs two adjacent free sites and creates two adsorbate units, two units recombine on desorption
* The surface starts clean and fills towards the equilibrium at 0.125, with the sticking coefficient falling from 0.42 to 0.08
* Comparison of the coverage after 10 ns (DSMCSurfState), number density and temperature (PartAnalyze) and the number of impacts (SurfaceAnalyze)
* `verify_analytic.py` checks coverage and impacts against a 0D solution of the surface model; pass a `*_ref.csv` to re-check the stored reference files
