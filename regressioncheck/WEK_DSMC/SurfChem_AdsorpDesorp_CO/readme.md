# DSMC - Adsorption/Desorption of CO on a catalytic wall
* Free-molecular gap (25 nm) between an inflow plane (BC_Xminus, surface flux) and a catalytic wall (BC_Xplus, 16 sub-surfaces), symmetric in y/z
* Literature reaction parameters for CO on palladium; wall temperature (1116.8 K) and the inflow density are chosen, so that desorption is resolved within the simulated 10 ns
* The surface starts at a coverage of 0.09 and desorbs towards the equilibrium at 0.030
* Comparison of the coverage after 10 ns (DSMCSurfState), number density and temperature (PartAnalyze) and the number of impacts (SurfaceAnalyze)
* `verify_analytic.py` checks coverage and impacts against a 0D solution of the surface model; pass a `*_ref.csv` to re-check the stored reference files
