MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"ftn -fpp -assume bscc -r8 -i4 -traceback -warn all -O2 -xCORE-AVX2 -vec-report0 -DPP_TimeDiscMethod=2 -DPP_NodeType=1 -"//&
"DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ -DPP_N=N -DLUSTRE -DMPI -DINTEL -I. -I/opt/cray/hdf5-parallel/1.8.13/INTEL/"//&
"140/include/ -mkl=sequential"
END MODULE MOD_PreProcFlags
