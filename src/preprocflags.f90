MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"ifort -fpp -assume bscc -r8 -i4 -traceback -warn all -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,noo"//&
"utput_conversion,pointer,bounds,uninit -DPP_TimeDiscMethod=2 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELE"//&
"MZ -DPP_N=N -DPARTICLES=T -DINTEL -I. -I../share/INTEL-SINGLE/hdf5-1.8.14/hdf5/include/ "
END MODULE MOD_PreProcFlags
