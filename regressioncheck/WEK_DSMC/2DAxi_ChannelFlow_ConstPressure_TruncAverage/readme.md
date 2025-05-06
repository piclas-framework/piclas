# DSMC - Pressure driven axisymmetric channel flow
* Simulation of a subsonic channel flow with O2
* Testing the adaptive surface flux boundary conditions (Type=1: Inflow, constant pressure and temperature and Type=2: Outflow with constant pressure)
* Inlet: 5 Pa, 300K, Outlet: 2.5 Pa, 300K
* Temporal evolution of the channel is sensitive to the chosen sampling method for the adaptive BC
* Comparing only the pressure with reference (absolute comparison to exclude very low mass flow values)
* Mass flow is not yet stationary and thus very different at the BCs