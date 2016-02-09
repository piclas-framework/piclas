!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbrevbiations for specific expressions are defined
!===================================================================================================================================

! Abbrevations
#ifdef SUN
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SUN COMPILER'
#  define IEEE_ISNAN
#elif SX
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SX COMPILER'
#elif PGI
#  define NO_ISNAN
#endif
#define __STAMP__ __FILE__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#ifdef MPI
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#else
#  define SWRITE WRITE
#  define IPWRITE(a,b) WRITE(a,b)0,
#endif
#define ERRWRITE(a,b) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)


#ifdef OPTIMIZED
#define PP_IJK     i,0,0
#define PP_ij      i,0
#else
#define PP_IJK     i,j,k
#define PP_ij      i,j
#endif
! Predefined "PARAMETER-like" variables
#define XI_MINUS   5
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

! Entry position in SideToElem
#define S2E_ELEM_ID        1
#define S2E_NB_ELEM_ID     2
#define S2E_LOC_SIDE_ID    3
#define S2E_NB_LOC_SIDE_ID 4
#define S2E_FLIP           5

! Entry position in SideToElem2
#define S2E2_ELEM_ID        1
#define S2E2_SIDE_ID        2
#define S2E2_LOC_SIDE_ID    3
#define S2E2_FLIP           4

! Entry position in ElemToSide
#define E2S_SIDE_ID 1
#define E2S_FLIP    2

! Entry position in ElemToElem
#define E2E_NB_ELEM_ID      1
#define E2E_NB_LOC_SIDE_ID  2

! Entry position in BC
#define BC_TYPE  1
#define BC_STATE 2
#define BC_ALPHA 3

!#define DEBUGMESH

#define PLANAR   0
#define BILINEAR 1
#define CURVED   2

! entries for PartHaloToProc
#define NATIVE_ELEM_ID  1
#define NATIVE_PROC_ID  2
#define LOCAL_PROC_ID   3
!#define NATIVE_SIDE_ID  1
#define LOCAL_SEND_ID   4
