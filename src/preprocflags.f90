MODULE MOD_PreProcFlags
IMPLICIT NONE
CHARACTER(LEN=1000) :: PREPROC_FLAGS = &
"gfortran -xf95-cpp-input -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU -Wall -O3 -march=native -DPP_TimeDi"//&
"scMethod=2 -DPP_NodeType=1 -DPP_Riemann=1 -DPP_nVar=8 -DPP_nElems=NELEMZ -DPP_N=N -DPARTICLES=T -DGNU -I. -I../share/GNU"//&
"-SINGLE/hdf5-1.8.14/hdf5/include/ "
END MODULE MOD_PreProcFlags
