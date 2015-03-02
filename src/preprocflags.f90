MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"mpif90 -fpp -assume bscc -r8 -i4 -traceback -warn all -O2 -xHost -vec-report0 -DPP_TimeDiscMethod=5 -DPP_NodeType=1 -DP"//&
"P_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ -DPP_N=N -DPARTICLES=T -DMPI -DINTEL -I. -I/opt/hdf5/hdf5-1.8.14/intel14/incl"//&
"ude/ "
END MODULE MOD_PreProcFlags
