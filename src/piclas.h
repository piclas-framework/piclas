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

#define SIZEOF_F(x) (STORAGE_SIZE(x)/8)

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  k).OR.x<-HUGE(1_  k))       CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ k).OR.x<-HUGE(1._ k))       CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#elif CRAY
#define CHECKSAFEINT(x,k)
#define CHECKSAFEREAL(x,k)
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#endif

! Test for equality: read description in src/globals/globals.f90 for further infos
! for variable relative tolerance
#define ALMOSTEQUALRELATIVE(x,y,tol)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(tol))
! for fixed relative tolerance (for double precision use twice the machine precision 2E-52 ~ 2.22e-16 -> 2*2.22E-16=4.44E-16)
#define ALMOSTEQUAL(x,y)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(4.441E-16))
#define ALMOSTZERO(x) (ABS(x).LE.(2.22e-16))

! Check for charged particles: x = iPart
#define CHARGEDPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)

! Check for particles to be interpolated or deposited: x = iPart
#if (PP_TimeDiscMethod==300) /*FP-Flow*/ 
#define PUSHPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define DEPOSITPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define INTERPOLATEPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#elif (PP_TimeDiscMethod==400) /*BGK*/
#define PUSHPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define DEPOSITPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define INTERPOLATEPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#else /*all other methods, mainly PIC*/
#define PUSHPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define DEPOSITPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#define INTERPOLATEPARTICLE(x) (ABS(Species(PartSpecies(x))%ChargeIC).GT.0.0)
#endif

#if USE_MPI
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#else
#  define SWRITE WRITE
#  define IPWRITE(a,b) WRITE(a,b)0,
#endif
#define ERRWRITE(a,b) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b)  IF(Logging) WRITE(UNIT_logOut,b)
#define LOGWRITE_BARRIER  IF(Logging) CALL ReOpenLogFile()
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

#ifdef OPTIMIZED
#define PP_IJK     i,0,0
#define PP_ij      i,0
#else
#define PP_IJK     i,j,k
#define PP_ij      i,j
#endif

#ifdef INTKIND8
#define MPI_INTEGER_INT_KIND MPI_INTEGER8
#else
#define MPI_INTEGER_INT_KIND MPI_INTEGER
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

! define side indeces used in sidetype array
#define PLANAR_RECT    0
#define PLANAR_NONRECT 1
#define BILINEAR       2
#define PLANAR_CURVED  3
#define CURVED         4

! entries for PartHaloToProc
#define NATIVE_ELEM_ID  1
#define NATIVE_PROC_ID  2
#define LOCAL_PROC_ID   3
!#define NATIVE_SIDE_ID  1
#define LOCAL_SEND_ID   4

! Entry position for interface type for selecting the corresponding Riemann solver
#define RIEMANN_VACUUM            0
#define RIEMANN_PML               1
#define RIEMANN_DIELECTRIC        2
#define RIEMANN_DIELECTRIC2VAC    3
#define RIEMANN_VAC2DIELECTRIC    4
#define RIEMANN_DIELECTRIC2VAC_NC 5
#define RIEMANN_VAC2DIELECTRIC_NC 6

! formats
! print to std out like  "    1.41421356237310E+000   -1.41421356237310E+000   -1.41421356237310E+000"
! (looks good and prevents the first digit of being a zero)
#define WRITEFORMAT '(ES25.14E3)'
! print to csv file like "0.1414213562373095E+001,-.1414213562373095E+001,-.1414213562373095E+001"
! (does not look that good but it saves disk space)
#define CSVFORMAT '(A1,E23.16E3)'

! Load Balance (LB) position in array for measuring the time that is spent on specific operations
#define LB_DG            1
#define LB_DGANALYZE     2
#define LB_DGCOMM        3
#define LB_PML           4
#define LB_EMISSION      5
#define LB_TRACK         6
#define LB_INTERPOLATION 7
#define LB_DEPOSITION    8
#define LB_CARTMESHDEPO  9
#define LB_PUSH          10
#define LB_PARTANALYZE   11
#define LB_PARTCOMM      12
#define LB_DSMC          13
#define LB_DSMCANALYZE   14
#define LB_SPLITMERGE    15
#define LB_UNFP          16
#define LB_SURF          17
#define LB_SURFFLUX      18
#define LB_SURFCOMM      19
#define LB_ADAPTIVE      20

#define LB_NTIMES        20

! DSMC_analyze indeces used in arrays
#define DSMC_VELOX       1
#define DSMC_VELOY       2
#define DSMC_VELOZ       3
#define DSMC_TEMPX       4
#define DSMC_TEMPY       5
#define DSMC_TEMPZ       6
#define DSMC_NUMDENS     7
#define DSMC_TVIB        8
#define DSMC_TROT        9
#define DSMC_TELEC       10
#define DSMC_SIMPARTNUM  11
#define DSMC_TEMPMEAN    12

#define DSMC_NVARS       12

! Sampwall_analyze indeces used in arrays
#define SAMPWALL_ETRANSOLD        1
#define SAMPWALL_ETRANSWALL       2
#define SAMPWALL_ETRANSNEW        3
#define SAMPWALL_EROTOLD          4
#define SAMPWALL_EROTWALL         5
#define SAMPWALL_EROTNEW          6
#define SAMPWALL_EVIBOLD          7
#define SAMPWALL_EVIBWALL         8
#define SAMPWALL_EVIBNEW          9
#define SAMPWALL_DELTA_MOMENTUMX  10
#define SAMPWALL_DELTA_MOMENTUMY  11
#define SAMPWALL_DELTA_MOMENTUMZ  12

#define SAMPWALL_NVARS            12
