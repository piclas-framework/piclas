# Implicit Maxwell and Particles
- The reggie is currently deactivated via Part-Species1-nInits = 0 as particle is lost during simulation
- The implicit solver cannot run without PIC-DoInterpolation=T as the array FieldAtParticle(:,:) is not allocated (this has changed,
as originally it simply was all zeros)
