# Rotational frame of reference: Subcycling
* Frame of reference is rotating with 5 revolutions per second
* A single particle is initially placed at coordinates (x,y) = (-0.25,-0.0) with a velocity vector of (10,10,0). It is anticipated to traverse a circular-like trajectory within the rotating frame of reference due to fictitious forces.
* A wall, rotating synchronously with the frame of reference, induces a specular reflection upon collision.
* A relatively large time step is employed (Delta_t=2E-3). Reference positions are calculated using a smaller time step (1E-5). Utilization of a subcycling with 200 substeps ensures identical positions in the regression test as in the reference case. Disabling subcycling results in failure of the regression test.
