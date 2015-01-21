MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"ifort -fpp -assume bscc -r8 -i4 -traceback -warn all -O2 -xHost -vec-report0 -DPP_TimeDiscMethod=2 -DPP_NodeType=1 -DPP"//&
"_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ -DPP_N=N -DPARTICLES=T -DINTEL -I. -I../share/INTEL-SINGLE/hdf5-1.8.14/hdf5/in"//&
"clude/ "
END MODULE MOD_PreProcFlags
