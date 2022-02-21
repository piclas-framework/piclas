# Integer Kind=8
- Test the functionality of using more than INTEGER(KIND=4), which is 2 147 483 647, particles 
  by setting all required arrays and variables to INTEGER(KIND=8)
- Cannot use 2.2 billion or more particles in this regression test, because this requires more RAM than available on the reggie server
  - if this is possible in the future, increase the number of particles above the INTEGER(KIND=4) limit
