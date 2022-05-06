# Background gas using regions
* Reservoir with 3 species and 3 background gas regions
  * 1 region with N2 (1E23/m3), 1 region with He (1E24/m3) and 1 region with N2/He (1E25/m3 in total)
* Corresponding number density and temparte is verified with DSMCState output in the respective species container
  * Compared as part of the regression test
* Interaction with particle species was verified by the collision probability
  * Comparison with DSMCState not possible due to strong fluctuations, third species is not inserted