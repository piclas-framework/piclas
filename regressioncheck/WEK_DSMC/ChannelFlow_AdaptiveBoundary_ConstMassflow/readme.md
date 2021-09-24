# DSMC - Massflow driven channel flow
* Simulation of a subsonic channel flow with O2
* Testing the adaptive surface flux boundary conditions (Type=4: Inflow, constant massflow and temperature and Type=2: Outflow with constant pressure)
* Massflow is taken from ChannelFlow_AdaptiveBoundary_ConstPressure regression test, compared to the same reference
* Inlet: 3.5E-14 kg/s, 300K, Outlet: 2.5 Pa, 300K
* Temporal evolution of the channel is sensitive to the chosen sampling method for the adaptive BC