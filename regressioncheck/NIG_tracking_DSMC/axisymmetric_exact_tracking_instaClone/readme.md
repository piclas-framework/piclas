# DSMC - Massflow driven axisymmetric pipe flow
* Simulation of a subsonic channel flow with O2
* Testing the adaptive surface flux boundary conditions (Type=4: Inflow, constant mass flow and temperature, Type=2: Outflow with constant pressure)
* Inlet: 5E-12 kg/s, 300K, Outlet: 0.02 Pa, 300K
* Testing the instantaneuos cloning feature with more than one clone: delayed cloning with only one allowed clone will abort
* Additional check whether particles leave the simulation domain due to a large timestep: Particles-Symmetry2DAxisymmetricExact = F will fail, since it neglects the curvature of the outer boundary