# Change Reference Frame
* The frame of references is localy changed along the z-axis. The reference frame is changed from a stationary frame to a rotating frame of reference and vise versa several times.
* A single particle is inserted at (x,y,z) = (0.25,0.25,0.5) with a constant speed of 1 m/s in z-direction. 
* As a consequence the particle path should be a straight line in the stationary frame of reference.
* In contrast the particle path should be circular in the rotating frame of reference.
* MPI test failed randomly during regressioncheck with following error output: cannot perform analyze Analyze_compare_column, because the shape of the data in file=[ParticlePosition.csv] is (998,) and that of the reference=[ParticlePosition_ref.csv] is (111,). They cannot be different!
