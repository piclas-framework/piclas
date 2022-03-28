# Multi-species background gas with MCC
* Testing the multi-species background gas with a CO2-N2-He reservoir using effective cross-section data (Phelps database, www.lxcat.net, retrieved on February 18, 2020)
* Background gas: 6E23 1/m3 (1/3 each), Electrons: n_e = 1E21 1/m3
* Adapting the constant velocity of the electrons from 0.01 eV to 100 eV
* Two additional values at 1500 eV and 2500 eV to test for extrapolation and potential cross-sections below 0, respectively (only for CO2 and He, N2 has values up to 10000 eV)
* Collisional rates compared with the theoretical results (see Figure)
* "get_values.py" script can be utilized to collect the collision rates from the reference databases