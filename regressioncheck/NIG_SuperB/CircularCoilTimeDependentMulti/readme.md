# SuperB: Circular coils (time-dependent)
- Magnetic field of a three circular coils with a current of 1e6 A that changes in time with a sin function
- The three different coils must have the same frequency but have different phases in this example (shifted by 120 and 240 degrees)
- Output of DG_Solution field to SuperB_CircularCoil_BGField.h5 with additional dimension that accounts for temporal change
- When running this reggie with piclas and MPI=4, the field is created and the read-in when the automatic load balance is performed
after the first time step
