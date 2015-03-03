MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"mpif90 -fpp -assume bscc -r8 -i4 -traceback -warn all -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,no"//&
"output_conversion,pointer,bounds,uninit -DPP_TimeDiscMethod=5 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NEL"//&
"EMZ -DPP_N=N -DPARTICLES=T -DMPI -DINTEL -I. -I/opt/hdf5/hdf5-1.8.14/intel14/include/ "
END MODULE MOD_PreProcFlags
