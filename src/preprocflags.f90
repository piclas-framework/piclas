MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"mpif90 -xf95-cpp-input -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU -Wall -Wextra -g -O0 -ggdb3 -fbounds-"//&
"check -finit-real=nan -fbacktrace -DPP_TimeDiscMethod=2 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ -D"//&
"PP_N=N -DPARTICLES=T -DMPI -DGNU -I. -I/opt/hdf5/hdf5-1.8.14/gnu48/include/ "
END MODULE MOD_PreProcFlags
