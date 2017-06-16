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
#ifndef __FILENAME__ 
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  k).OR.x<-HUGE(1_  k))       CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ k).OR.x<-HUGE(1._ k))       CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#endif

! Test for equality: read description in src/globals/globals.f90 for further infos
! for variable relative tolerance
#define ALMOSTEQUALRELATIVE(x,y,tol)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(tol))
! for fixed relative tolerance (for double precision use twice the machine precision 2E-52 ~ 2.22e-16 -> 2*2.22E-16=4.44E-16)
#define ALMOSTEQUAL(x,y)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(4.44E-16))
#define ALMOSTZERO(x) (ABS(x).LE.(2.22e-16))

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

! Entry position in BC
#define MI_SIDEID 1
#define MI_FLIP   2

!#define DEBUGMESH

#define PLANAR_RECT   0
#define PLANAR_NONRECT   1
#define BILINEAR 2
#define PLANAR_CURVED   3
#define CURVED   4

! entries for PartHaloToProc
#define NATIVE_ELEM_ID  1
#define NATIVE_PROC_ID  2
#define LOCAL_PROC_ID   3
!#define NATIVE_SIDE_ID  1
#define LOCAL_SEND_ID   4

! format
#define OUTPUTFORMAT '(E25.14E3)'
