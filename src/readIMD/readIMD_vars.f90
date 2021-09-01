module mod_readIMD_vars
#if USE_MPI && defined(PARTICLES)

implicit none
public
save

! ==============================================================================
logical                       :: useIMDresults
character(len=300)            :: filenameIMDresults

!Debugging vars
logical                       :: killPIClas
#endif /*USE_MPI && defined(PARTICLES)*/
end module mod_readIMD_vars
