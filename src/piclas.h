!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!===================================================================================================================================

! From include/petsc/finclude/petscsys.h: #define PetscCallA(func) call func; CHKERRA(ierr)
#if USE_PETSC_FIX317
#define PetscCallA(a) CALL a; PetscCall(ierr)
#endif

! Abbrevations
#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

! Calculate GCC version
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

! Size of data types
#define SIZEOF_F(x) (STORAGE_SIZE(x)/8)
#define SIZE_LOG  KIND(.TRUE.)
#define SIZE_INT  KIND(INT(1))
#define SIZE_INT4 4
#define SIZE_INT8 8
#define SIZE_REAL KIND(REAL(1))
#define SIZE_CHAR KIND('a')

#ifdef DEBUG_MEMORY
#define Allocate_Shared(a,b,c)   Allocate_Shared_DEBUG(a,b,c,'b',TRIM(__FILE__),__LINE__)
#endif

#ifdef MEASURE_MPI_WAIT
! Field solver
#if USE_HDG
#define MPIW8SIZEFIELD 4
#else
#define MPIW8SIZEFIELD 2
#endif
! Particle solver
#ifdef PARTICLES
#define MPIW8SIZEPART 6
#else
#define MPIW8SIZEPART 0
#endif
! Combination
#define MPIW8SIZE (2+MPIW8SIZEFIELD+MPIW8SIZEPART)
#endif

! Deactivate PURE subroutines/functions when using DEBUG
#if USE_DEBUG
#define PPURE
#else
#define PPURE PURE
#endif

! Replace function with dummy function and additional argument
#define UNLOCK_AND_FREE(a)   UNLOCK_AND_FREE_DUMMY(a,'a')

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(INT( 1,KIND=k)).OR.x<-HUGE(INT( 1,KIND=k))) CALL Abort(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(REAL(1,KIND=k)).OR.x<-HUGE(REAL(1,KIND=k))) CALL Abort(__STAMP__,'Real conversion failed: out of range!')
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
#define ALMOSTALMOSTEQUAL(x,y)  (ABS((x)-(y)).LE.MAX(ABS(x),ABS(y))*(1E-10))
#define ALMOSTZERO(x) (ABS(x).LE.(2.22e-16))

! Check if the exponent is within the range of machine precision, RANGE() gives the maximum exponent of the given variable type
#define CHECKEXP(x) (ABS(x).LT.REAL(RANGE(1.)))

#if USE_MPI
#  define SWRITE IF(MPIRoot) WRITE
#if USE_LOADBALANCE
#  define LBWRITE IF(MPIRoot.AND.(.NOT.PerformLoadBalance)) WRITE
#else /*USE_LOADBALANCE*/
#  define LBWRITE IF(MPIRoot) WRITE
#endif /*USE_LOADBALANCE*/
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#  define LWRITE IF(myComputeNodeRank.EQ.0) WRITE
#  define GETTIME(a) a=MPI_WTIME()
#else
#  define SWRITE WRITE
#  define LBWRITE WRITE
#  define IPWRITE(a,b) WRITE(a,b)0,
#  define LWRITE WRITE
#  define GETTIME(a) CALL CPU_TIME(a)
#endif
#define ERRWRITE(a,b) WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b)  IF(Logging) WRITE(UNIT_logOut,b)
#define LOGWRITE_BARRIER  IF(Logging) CALL ReOpenLogFile()
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)
#define SNULLIFY(A)    IF(ASSOCIATED(A)) NULLIFY(A)

#if USE_MPI
#define ALLOCPOINT POINTER
#define ADEALLOCATE(A) IF(ASSOCIATED(A)) NULLIFY(A)
#else
#define ALLOCPOINT ALLOCATABLE
#define ADEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)
#endif

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

! number of entry in each line of ElemInfo
#define ELEMINFOSIZE_H5   6
#if USE_MPI
#define ELEMINFOSIZE      8
#else
#define ELEMINFOSIZE      6
#endif /* USE_MPI*/
! ElemInfo in H5 file
#define ELEM_TYPE         1
#define ELEM_ZONE         2
#define ELEM_FIRSTSIDEIND 3
#define ELEM_LASTSIDEIND  4
#define ELEM_FIRSTNODEIND 5
#define ELEM_LASTNODEIND  6
! ElemInfo for shared array
#define ELEM_RANK         7
#define ELEM_HALOFLAG     8
! number of entries in each line of SideInfo
#define SIDEINFOSIZE_H5   5
#define SIDEINFOSIZE      8
#define SIDE_TYPE         1
#define SIDE_ID           2
#define SIDE_NBELEMID     3
#define SIDE_FLIP         4
#define SIDE_BCID         5
#define SIDE_ELEMID       6
#define SIDE_LOCALID      7
#define SIDE_NBSIDEID     8
#define SIDE_NBELEMTYPE   9
! surface sampling entries
#define SURF_SIDEID       1
#define SURF_RANK         2
#define SURF_LEADER       3

! Predefined "PARAMETER-like" variables
#define XI_MINUS   5
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

! Entry position in ElemBCSides
#define ELEM_NBR_BCSIDES   1
#define ELEM_FIRST_BCSIDE  2

! Entry position in FIBGMToProc
#define FIBGM_FIRSTPROCIND 1
#define FIBGM_NPROCS       2
#define FIBGM_NLOCALPROCS  3

! Entry position in SideBCMetrics
#define BCSIDE_SIDEID      1
#define BCSIDE_ELEMID      2
#define BCSIDE_DISTANCE    3
#define BCSIDE_RADIUS      4
#define BCSIDE_ORIGINX     5
#define BCSIDE_ORIGINY     6
#define BCSIDE_ORIGINZ     7

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

! entries for PartExchange
#define EXCHANGE_PROC_SIZE 2
#define EXCHANGE_PROC_TYPE 1
#define EXCHANGE_PROC_RANK 2

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
#define SAMPWALL_ETRANSNEW        2
#define SAMPWALL_EROTOLD          3
#define SAMPWALL_EROTNEW          4
#define SAMPWALL_EVIBOLD          5
#define SAMPWALL_EVIBNEW          6
#define SAMPWALL_EELECOLD         7
#define SAMPWALL_EELECNEW         8
#define SAMPWALL_DELTA_MOMENTUMX  9
#define SAMPWALL_DELTA_MOMENTUMY  10
#define SAMPWALL_DELTA_MOMENTUMZ  11

#define SAMPWALL_NVARS            11

#define MACROSURF_NVARS           6

! Tracking method
#define REFMAPPING    1
#define TRACING       2
#define TRIATRACKING  3

! Time Step Minimum: dt_Min
#define DT_MIN        1
#define DT_ANALYZE    2
#define DT_END        3
#define DT_BR_SWITCH  4

! Secondary electron emission
#define SEE_MODELS_ID 4,5,6,7,8,9,10,11

#if USE_HDG
! HDG Dirichlet BC Side IDs
#define HDGDIRICHLETBCSIDEIDS 2,4,5,6,7,8,50,51,52,60
#endif