#ifdef MPI
#ifdef PARTICLES
module mod_readIMD_vars

implicit none
public
save

! ==============================================================================
logical                       :: useIMDresults
character(len=300)            :: filenameIMDresults

!Debugging vars
logical                       :: killPIClas
end module mod_readIMD_vars
#endif /*MPI*/
#endif /*PARTICLES*/
