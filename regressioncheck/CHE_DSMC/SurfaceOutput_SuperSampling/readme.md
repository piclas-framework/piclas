# Flow inside cube: Testing of surface output with super-sampling
* Constant inflow with a fixed number density and circular inflow defined by a circle
* Circle should be visible in the sampled heat flux using nSurfSample = 6 despite only having 4 volume cells
* Comparison of heat flux in HDF5 and in converted VTK file