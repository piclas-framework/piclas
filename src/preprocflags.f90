MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"ifort -fpp -assume bscc -r8 -i4 -traceback -warn all -g -O0 -fpe0 -traceback -c"//&
"heck all,noarg_temp_created,noformat,nooutput_conversion,pointer,bounds,uninit -"//&
"DPP_TimeDiscMethod=5 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELE"//&
"MZ -DPP_N=N -DPARTICLES=T -DINTEL -I. -I../share/INTEL-SINGLE/hdf5-1.8.13/hdf5/i"//&
"nclude "
END MODULE MOD_PreProcFlags
