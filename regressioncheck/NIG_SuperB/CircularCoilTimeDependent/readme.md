# SuperB: Circular coil (time-dependent)
- Magnetic field of a single circular coil with a current of 1 A that changes in time with a sin function
- Output of DG_Solution field to SuperB_CircularCoil_BGField.h5 with additional dimension that accounts for temporal change
- When running this reggie with piclas and MPI=4, the field is created and the read-in when the automatic load balance is performed
after the first time step
