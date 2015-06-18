MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"ftn -fpp -assume bscc -r8 -i4 -traceback -warn all -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,noout"//&
"put_conversion,pointer,bounds,uninit -DPP_TimeDiscMethod=2 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ"//&
" -DPP_N=N -DPARTICLES=T -DLUSTRE -DMPI -DINTEL -I. -I/opt/cray/hdf5-parallel/1.8.13/INTEL/140/include/ "
END MODULE MOD_PreProcFlags
