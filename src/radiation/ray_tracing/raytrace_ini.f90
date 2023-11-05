!==================================================================================================================================
! Copyright (c) 2023 - 2023 Marcel Pfeiffer, Stephen Copplestone
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_RayTracing_Init
!===================================================================================================================================
! Initialization of Radiation Transport
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC :: InitRayTracing, DefineParametersRayTracing, FinalizeRayTracing
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for FP-Flow
!==================================================================================================================================
SUBROUTINE DefineParametersRayTracing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Ray Tracing")
CALL prms%CreateIntOption(       'RayTracing-NumRays'         , 'Number of emitted rays from particle boundary with index [RayTracing-PartBound]')
CALL prms%CreateRealArrayOption( 'RayTracing-RayDirection'    , 'Direction vector for ray emission. Will be normalized after read-in.' , no=3)
CALL prms%CreateIntOption(       'RayTracing-PartBound'       , 'Particle boundary ID where rays are emitted from' , '0')
CALL prms%CreateRealOption(      'RayTracing-PulseDuration'   , 'Pulse duration tau for a Gaussian-type pulse with I~exp(-(t/tau)^2) [s]'                  )
CALL prms%CreateIntOption(       'RayTracing-NbrOfPulses'     , 'Number of pulses [-]'                                                                     , '1')
CALL prms%CreateRealOption(      'RayTracing-WaistRadius'     , 'Beam waist radius (in focal spot) w_b for Gaussian-type pulse with I~exp(-(r/w_b)^2) [m]' , '0.0')
CALL prms%CreateRealOption(      'RayTracing-WaveLength'      , 'Beam wavelength [m]'                                                                      )
CALL prms%CreateRealOption(      'RayTracing-RepetitionRate'  , 'Pulse repetition rate (pulses per second) [Hz]'                                           )
CALL prms%CreateRealOption(      'RayTracing-PowerDensity'    , 'Average pulse power density (power per area) [W/m2]')
CALL prms%CreateLogicalOption(   'RayTracing-ForceAbsorption' , 'Surface photon sampling is performed independent of the actual absorption/reflection outcome (default=T)', '.TRUE.')
CALL prms%CreateIntOption(       'RayTracing-NMax'            , 'Maximum polynomial degree within refined volume elements for photon tracking (p-adaption)')
CALL prms%CreateIntOption(       'RayTracing-VolRefineMode'   , 'High-order ray tracing volume sampling refinement method:\n'//&
                                                               ' 0: do nothing (default)\n'//&
                                                               ' 1: refine below user-defined z-coordinate with NMax\n'//&
                                                               ' 2: scale N according to the mesh element volume between NMin>=1 and NMax>=2\n'//&
                                                               ' 3: refine below user-defined z-coordinate and scale N according to the mesh element volume between NMin>=1 and NMax>=2\n'//&
                                                               '    (consider only elements below the user-defined z-coordinate for the scaling)'&
                                                              , '0')
CALL prms%CreateRealOption(      'RayTracing-VolRefineModeZ'  , 'Z-coordinate for switching between NMin (pos>z) and NMax (pos<z) depending on element position for high-order ray tracing')
CALL prms%CreateIntOption(       'RayTracing-nSurfSample'     , 'Define polynomial degree of ray tracing BC sampling. Default: nSurfSample (which itself defaults to NGeo)')
CALL prms%CreateStringOption(    'RayTracing-NodeType'        , 'Node type for volume and surface ray tracing super sampling:  VISU, VISU_INNER or GAUSS', 'VISU_INNER')


END SUBROUTINE DefineParametersRayTracing


SUBROUTINE InitRayTracing()
!===================================================================================================================================
! Initialization of the radiation transport solver
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_RayTracing_Vars
USE MOD_ReadInTools            ,ONLY: GETREAL,GETREALARRAY,GETINT,GETLOGICAL,PrintOption,GETSTR
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_RadiationTrans_Vars    ,ONLY: RadiationAbsorptionModel,RadObservationPointMethod
USE MOD_Interpolation_Vars     ,ONLY: NodeType,NodeTypeVISU
USE MOD_Particle_Boundary_Vars ,ONLY: nSurfSample
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_Photon_Tracking        ,ONLY: InitPhotonSurfSample
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: factor,SurfaceNormal(3),alpha
CHARACTER(LEN=3)  :: hilf ! auxiliary variable for INTEGER -> CHARACTER conversion
!===================================================================================================================================
IF(.NOT.UseRayTracing) RETURN
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT RAY TRACING MODEL ...'

! Do not absorb rays within the volume!
RadiationAbsorptionModel = 0
RadObservationPointMethod = 0

! Get index of particle boundary from which rays are emitted
RayPartBound = GETINT('RayTracing-PartBound')
IF(RayPartBound.LE.0) CALL CollectiveStop(__STAMP__,'RayTracing-PartBound must be > 0 to activate ray tracing on this boundary!')

! Get ray parameters
Ray%PulseDuration  = GETREAL('RayTracing-PulseDuration')
Ray%NbrOfPulses    = GETINT('RayTracing-NbrOfPulses')
Ray%tShift         = SQRT(8.0) * Ray%PulseDuration
Ray%WaistRadius    = GETREAL('RayTracing-WaistRadius')
Ray%WaveLength     = GETREAL('RayTracing-WaveLength')
Ray%RepetitionRate = GETREAL('RayTracing-RepetitionRate')
Ray%Period         = 1./Ray%RepetitionRate
Ray%PowerDensity   = GETREAL('RayTracing-PowerDensity')
Ray%Direction      = GETREALARRAY('RayTracing-RayDirection',3)
Ray%Direction      = UNITVECTOR(Ray%Direction)
WRITE(UNIT=hilf,FMT='(I0)') nSurfSample
Ray%nSurfSample    = GETINT('RayTracing-nSurfSample',hilf)
Ray%NodeType       = TRIM(GETSTR('RayTracing-NodeType'))
SELECT CASE(TRIM(Ray%NodeType))
CASE('VISU','VISU_INNER','GAUSS')
  ! Do nothing
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown node type for ray tracing: '//TRIM(Ray%NodeType)//'. Select VISU, VISU_INNER or GAUSS')
END SELECT

! Parameters only require for actual ray tracing computation
IF(PerformRayTracing)THEN
  NumRays            = GETINT('RayTracing-NumRays')
  RayForceAbsorption = GETLOGICAL('RayTracing-ForceAbsorption')
  Ray%VolRefineMode  = GETINT('RayTracing-VolRefineMode')
#if ! (CORE_SPLIT==0)
  ! Sanity check: ElemVolume_Shared is only built for nComputeNodeElems and not nComputeNodeTotalElems. Maybe more containers are
  ! similarly not fully built when running multi-node
  CALL CollectiveStop(__STAMP__,'Ray tracing implemented for node-level splitting; all global elements must be on the compute-node')
#endif /*! (CORE_SPLIT==0)*/
END IF ! PerformRayTracing

! Output of high-order p-adaptive info
Ray%NMin = 1 ! GETINT('RayTracing-NMin')
WRITE(UNIT=hilf,FMT='(I3)') PP_N
Ray%Nmax = GETINT('RayTracing-Nmax',hilf)
IF(Ray%Nmax.LT.Ray%Nmax) CALL abort(__STAMP__,'RayTracing-Nmax cannot be smaller than Nmin=',IntInfoOpt=Ray%NMin)

! Build surface and volume containers
CALL InitPhotonSurfSample()
CALL InitHighOrderRaySampling()

ASSOCIATE( &
      E0      => Ray%Energy             ,&
      wb      => Ray%WaistRadius        ,&
      tau     => Ray%PulseDuration      ,&
      I0      => Ray%IntensityAmplitude ,&
      tShift  => Ray%tShift             ,&
      Period  => Ray%Period             ,&
      tActive => Ray%tActive            ,&
      A       => Ray%Area               )

  ! ATTENTION: Rectangle only and uses GEO min/max in x- and y-direction!
  ! TODO: Ray emission area from chosen boundary surface?
  A = (GEO%xmaxglob-GEO%xminglob) * (GEO%ymaxglob-GEO%yminglob)
  ! Normal vector of the ray emission area
  SurfaceNormal = (/ 0., 0., 1. /)
  ! Angle between emitted rays and emission area
  alpha = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(Ray%Direction,SurfaceNormal))))

  ! Derived quantities
  Ray%Power = Ray%PowerDensity * A ! adjust power from [W/m2] to [W]
  E0 = Ray%Power / Ray%RepetitionRate

  ! Generate two base vectors perpendicular to the ray direction
  CALL OrthoNormVec(Ray%Direction,Ray%BaseVector1IC,Ray%BaseVector2IC)

  ! Calculate the peak intensity (uncorrected)
  I0 = E0 / (SQRT(PI)*tau*A)

  ! Correction factor due to temporal cut-off of the Gaussian pulse
  ! no need for correction in space because the function is not cut-off in space
  ! just consider the temporal cut-off for the rectangle
  factor = ERF(tShift/tau)
  factor = SQRT(PI)*tau*A
  I0 = E0 / factor

  ! Sanity check: overlapping of pulses is not implemented (use multiple emissions for this)
  IF(2.0*tShift.GT.Period) CALL abort(__STAMP__,'Pulse length (2*tShift) is greater than the pulse period. This is not implemented!')

  ! Active pulse time
  tActive = REAL(Ray%NbrOfPulses - 1)*Period + 2.0*tShift
END ASSOCIATE

CALL PrintOption('Rectangular ray emission area: A [m2]'                             , 'CALCUL.' , RealOpt=Ray%Area)
CALL PrintOption('Angle between emission area normal and ray direction: alpha [deg]' , 'CALCUL.' , RealOpt=alpha)
! Increased energy in the volume due to the increased optical path (only 2D approximation)
IF(ABS(alpha)+1e-4.LT.90.0)THEN
  CALL PrintOption('Enhancement factor for energy deposited in the volume [-]'       , 'CALCUL.' , RealOpt=1.0/COS(alpha*PI/180.0))
ELSE
  CALL PrintOption('Enhancement factor for energy deposited in the volume [-]'       , 'CALCUL.' , RealOpt=1.0)
END IF ! ABS(alpha).GT.0.0
CALL PrintOption('Single pulse energy [J]'                                           , 'CALCUL.' , RealOpt=Ray%Energy)
CALL PrintOption('Intensity amplitude: I0 [W/m^2]'                                   , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Corrected Intensity amplitude: I0_corr [W/m^2]'                    , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Pulse period (Time between maximum of two pulses) [s]'             , 'CALCUL.' , RealOpt=Ray%Period)
CALL PrintOption('Temporal pulse width (pulse time 2x tShift) [s]'                   , 'CALCUL.' , RealOpt=2.0*Ray%tShift)
CALL PrintOption('Pulse will end at tActive (pulse final time) [s]'                  , 'CALCUL.' , RealOpt=Ray%tActive)

LBWRITE(UNIT_stdOut,'(A)')' INIT RAY TRACING MODEL DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRayTracing


!===================================================================================================================================
!> Build all high-order mappings required for ray trace sampling in the volume on a p-adaptive polynomial basis
!===================================================================================================================================
SUBROUTINE InitHighOrderRaySampling()
! MODULES
USE MOD_PreProc
USE MOD_Globals            ,ONLY: abort,IERROR,UNIT_StdOut
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems,nElems,offSetElem
!USE MOD_RayTracing_Vars    ,ONLY: N_VolMesh_Ray
USE MOD_RayTracing_Vars    ,ONLY: N_DG_Ray,Ray,U_N_Ray,U_N_Ray_loc,nVarRay,PerformRayTracing
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID,GetGlobalElemID
USE MOD_ReadInTools        ,ONLY: GETREAL
USE MOD_Particle_Mesh_Vars ,ONLY: ElemVolume_Shared,ElemBaryNGeo
#if USE_MPI
USE MOD_Globals            ,ONLY: MPI_COMM_PICLAS
USE MPI
!USE MOD_Particle_Mesh_Vars ,ONLY: NodeCoords_Shared
USE MOD_RayTracing_Vars    ,ONLY: N_DG_Ray_Shared,N_DG_Ray_Shared_Win
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars   ,ONLY: MPI_COMM_SHARED,myComputeNodeRank,nComputeNodeProcessors,nComputeNodeTotalElems
#else
!USE MOD_Mesh_Vars          ,ONLY: NodeCoords
USE MOD_Mesh_Vars          ,ONLY: nElems
#endif /*USE_MPI*/
#if defined(CODE_ANALYZE)
USE MOD_Globals            ,ONLY: nProcessors
#endif /*defined(CODE_ANALYZE)*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: Nloc,iCNElem,firstElem,lastElem,iGlobalElem,iElem
REAL              :: VolMin,VolMax,m,NReal
#if defined(CODE_ANALYZE)
INTEGER           :: CNElemID
LOGICAL,PARAMETER :: debugRay=.FALSE.
#endif /*defined(CODE_ANALYZE)*/
LOGICAL           :: FoundElem
CHARACTER(LEN=40) :: hilf
!===================================================================================================================================
! When ray tracing is performed, these arrays are created, otherwise (restart or no actual ray tracing) they are read from h5
IF(PerformRayTracing)THEN

#if USE_MPI
  CALL Allocate_Shared((/nGlobalElems/),N_DG_Ray_Shared_Win,N_DG_Ray_Shared)
  CALL MPI_WIN_LOCK_ALL(0,N_DG_Ray_Shared_Win,IERROR)
  N_DG_Ray => N_DG_Ray_Shared
  ! only CN root initializes
  IF (myComputeNodeRank.EQ.0) N_DG_Ray = Ray%NMin ! default
  ! This sync/barrier is required as it cannot be guaranteed that the zeros have been written to memory by the time the MPI_REDUCE
  ! is executed (see MPI specification). Until the Sync is complete, the status is undefined, i.e., old or new value or utter nonsense.
  CALL BARRIER_AND_SYNC(N_DG_Ray_Shared_Win,MPI_COMM_SHARED)
#else
  ALLOCATE(N_DG_Ray(nGlobalElems))
  N_DG_Ray = Ray%NMin ! default
#endif /*USE_MPI*/

  ! Select volumetric resolution
  IF(Ray%NMin.NE.Ray%NMax)THEN
    ! Only set variable N if NMin and NMax are not the same
    SELECT CASE(Ray%VolRefineMode)
    CASE(-1)
      ! -1: odd elements 1, even elements NMax
      DO iElem = 1, PP_nElems
        IF(MOD(iElem,2).EQ.0)THEN
          ! even
          N_DG_Ray(iElem) = Ray%Nmax
        ELSE
          ! odd
          N_DG_Ray(iElem) = 1
        END IF ! MOD(iElem,2).EQ.0
      END DO ! iElem = 1, PP_nElems
    CASE(0)
      ! 0: do nothing (default)
    CASE(1,2,3)
      ! Set first and last elem loop indices
#if USE_MPI
      firstElem = INT(REAL( myComputeNodeRank   )*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))+1
      lastElem  = INT(REAL((myComputeNodeRank+1))*REAL(nComputeNodeTotalElems)/REAL(nComputeNodeProcessors))
#else
      firstElem = 1
      lastElem  = nElems
#endif
      ! 1: refine below user-defined z-coordinate with NMax
      IF(Ray%VolRefineMode.NE.2)THEN
        WRITE(UNIT=hilf,FMT=WRITEFORMAT) 1.0E200!HUGE(1.0) -> HUGE produces IEEE overflow
        Ray%VolRefineModeZ = GETREAL('RayTracing-VolRefineModeZ',TRIM(hilf))
        DO iCNElem = firstElem, lastElem
          iGlobalElem = GetGlobalElemID(iCNElem)
          !IPWRITE(UNIT_StdOut,*) "iCNElem,iGlobalElem =", iCNElem,iGlobalElem
          IF(ElemBaryNGeo(3,iCNElem).LT.Ray%VolRefineModeZ)THEN
            N_DG_Ray(iGlobalElem) = Ray%Nmax
          ELSE
            N_DG_Ray(iGlobalElem) = Ray%Nmin
          END IF ! ElemBaryNGeo(3,iCNElem).LT.Ray%VolRefineModeZ
        END DO ! iCNElem = firstElem, lastElem
      ELSE
        Ray%VolRefineModeZ = 1.0E200 ! dummy
      END IF ! Ray%VolRefineMode.NE.2

      ! 2: scale N according to the mesh element volume between NMin>=1 and NMax>=2
      ! 3: refine below user-defined z-coordinate and scale N according to the mesh element volume between NMin>=1 and NMax>=2
      !    (consider only elements below the user-defined z-coordinate for the scaling)
      IF((Ray%VolRefineMode.EQ.2).OR.(Ray%VolRefineMode.EQ.3))THEN
        ! Get global min and max volume: Only consider elements below the z-coordinate
        VolMin = HUGE(1.)
        VolMax = -HUGE(1.)
        FoundElem = .FALSE.
        DO iCNElem = firstElem, lastElem
          iGlobalElem = GetGlobalElemID(iCNElem)
          IF((ElemBaryNGeo(3,iCNElem).LT.Ray%VolRefineModeZ).OR.(Ray%VolRefineMode.EQ.2))THEN
            VolMin = MIN(VolMin, ElemVolume_Shared(iCNElem))
            VolMax = MAX(VolMax, ElemVolume_Shared(iCNElem))
            FoundElem = .TRUE.
          END IF
        END DO ! iCNElem = firstElem, lastElem
#if USE_MPI
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, VolMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_PICLAS, IERROR)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, VolMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_PICLAS, IERROR)
#endif /*USE_MPI*/

        ! Loop over all elements again and scale the polynomial degree using NINT.
        ! Check if the volumes of the elements are almost equal
        IF((VolMax.GT.VolMin).AND.(.NOT.ALMOSTEQUALRELATIVE(VolMax,VolMin,1e-2)))THEN
          ! Get the slope of the linear interpolation function between the maximum and minimum of the element volumes
          m = REAL(Ray%Nmax-Ray%Nmin)/(VolMax-VolMin)
          IF(FoundElem)THEN
            ! Loop over process elements
            DO iCNElem = firstElem, lastElem
              iGlobalElem = GetGlobalElemID(iCNElem)
              IF(ElemBaryNGeo(3,iCNElem).LT.Ray%VolRefineModeZ)THEN
                NReal = m * (ElemVolume_Shared(iCNElem)-VolMin) + REAL(Ray%Nmin)
                N_DG_Ray(iGlobalElem) = NINT(NReal)
              END IF
            END DO ! iCNElem = firstElem, lastElem
          END IF ! FoundElem
        END IF ! (VolMax.GT.VolMin).AND.(.NOT.ALMOST)

      END IF ! Ray%VolRefineMode.EQ.3
    CASE DEFAULT
      ! Debugging:
#if defined(CODE_ANALYZE)
      IF(nProcessors.GT.1) CALL abort(__STAMP__,'This only works for single-core runs')
      ! 3D box test case with diagonal rays
      IF(debugRay)THEN
        N_DG_Ray = Ray%Nmax
        DO iElem = 1, PP_nElems
          CNElemID = GetCNElemID(iElem)
          ASSOCIATE( &
                x => ElemBaryNGeo(1,CNElemID),&
                y => ElemBaryNGeo(2,CNElemID),&
                z => ElemBaryNGeo(3,CNElemID))
            IF(y+z.GE.1.40)THEN
              N_DG_Ray(iElem) = 1
              CYCLE
            END IF ! y+z.GT.1.5

            IF(y+z.LE.0.6)THEN
              N_DG_Ray(iElem) = 1
              CYCLE
            END IF ! y+z.LT.0.5

            IF(y+z.LT.0.9)THEN
              N_DG_Ray(iElem) = Ray%Nmax-1
              CYCLE
            END IF ! y+z.GT.1.00001

            IF(y+z.GT.1.1)THEN
              N_DG_Ray(iElem) = Ray%Nmax-1
              CYCLE
            END IF ! y+z.GT.1.00001
          END ASSOCIATE
        END DO ! iElem = 1, PP_nElems
      END IF ! debugRay
#else
      CALL abort(__STAMP__,'RayTracing-VolRefineMode unknown: ',IntInfoOpt=Ray%VolRefineMode)
#endif /*defined(CODE_ANALYZE)*/
    END SELECT

#if USE_MPI
    CALL BARRIER_AND_SYNC(N_DG_Ray_Shared_Win,MPI_COMM_SHARED)
#endif /*USE_MPI*/
  END IF ! Ray%NMin.NE.Ray%NMax

  ! Sanity check
  IF(ANY(N_DG_Ray.LE.0)) CALL abort(__STAMP__,'N_DG_Ray cannot contain zeros!')

  !ALLOCATE(N_VolMesh_Ray(1:nGlobalElems))
#if USE_MPI
  !CALL BuildElem_xGP_RayTrace(NodeCoords_Shared)
#else
  !CALL BuildElem_xGP_RayTrace(NodeCoords)
#endif /*USE_MPI*/

  ! The global DG solution in physical space
  ALLOCATE(U_N_Ray(1:nGlobalElems))
  DO iGlobalElem = 1, nGlobalElems
    Nloc = N_DG_Ray(iGlobalElem)
    ALLOCATE(U_N_Ray(iGlobalElem)%U(nVarRay,0:Nloc,0:Nloc,0:Nloc))
    U_N_Ray(iGlobalElem)%U = 0.
  END DO ! iGlobalElem = 1, nGlobalElems

  ! The local DG solution in physical space
  ALLOCATE(U_N_Ray_loc(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Ray(iElem+offSetElem)
    ALLOCATE(U_N_Ray_loc(iElem)%U(nVarRay,0:Nloc,0:Nloc,0:Nloc))
    U_N_Ray_loc(iElem)%U = 0.
  END DO ! iElem = 1, nElems

ELSE

END IF ! PerformRayTracing


CALL BuildNInterAndVandermonde()


END SUBROUTINE InitHighOrderRaySampling


!!==================================================================================================================================
!!> This routine takes the equidistant node coordinates of the mesh (on NGeo+1 points) and uses them to build the coordinates
!!> of solution/interpolation points of type NodeType on polynomial degree Nloc (Nloc+1 points per direction).
!!> Output: Elem_xGP(:,:,:,:) for each element with variably N
!!==================================================================================================================================
!SUBROUTINE BuildElem_xGP_RayTrace(NodeCoords)
!! MODULES
!USE MOD_Globals
!USE MOD_PreProc
!USE MOD_Mesh_Vars          ,ONLY: NGeo,nGlobalElems
!USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType
!USE MOD_RayTracing_Vars    ,ONLY: Ray,N_VolMesh_Ray,N_DG_Ray
!USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
!USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ, ChangeBasis3D
!USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
!!----------------------------------------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!REAL,INTENT(IN)               :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nGlobalElems)         !< Equidistant mesh coordinates
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                       :: iGlobalElem,Nloc
!
!TYPE VdmType
!  REAL, ALLOCATABLE           :: Vdm_EQNGeo_CLNloc(:,:)
!  REAL, ALLOCATABLE           :: Vdm_CLNloc_Nloc  (:,:)
!END TYPE VdmType
!
!TYPE(VdmType), DIMENSION(:), ALLOCATABLE :: Vdm
!
!!==================================================================================================================================
!
!! Build Vdm for every degree
!ALLOCATE(Vdm(Ray%Nmin:Ray%Nmax))
!DO Nloc = Ray%Nmin, Ray%Nmax
!  ALLOCATE(Vdm(Nloc)%Vdm_EQNGeo_CLNloc(0:Nloc,0:NGeo))
!  ALLOCATE(Vdm(Nloc)%Vdm_CLNloc_Nloc(0:Nloc,0:Nloc))
!  CALL GetVandermonde(NGeo, NodeTypeVISU, NLoc, NodeTypeCL, Vdm(Nloc)%Vdm_EQNGeo_CLNloc,  modal=.FALSE.)
!  CALL GetVandermonde(Nloc, NodeTypeCL  , Nloc, NodeType  , Vdm(Nloc)%Vdm_CLNloc_Nloc,     modal=.FALSE.)
!
!  ! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!  !       Important for curved meshes if NGeo<N, no effect for N>=NGeo
!
!  !1.a) Transform from EQUI_NGeo to solution points on Nloc
!  Vdm(Nloc)%Vdm_EQNGeo_CLNloc=MATMUL(Vdm(Nloc)%Vdm_CLNloc_Nloc, Vdm(Nloc)%Vdm_EQNGeo_CLNloc)
!END DO ! Nloc = Ray%Nmin, Ray%Nmax
!
!! Set Elem_xGP for each element
!DO iGlobalElem=1,nGlobalElems
!  Nloc = N_DG_Ray(iGlobalElem)
!
!  ! TODO: Currently each process has all global xGP (maybe put unrolled into a shared array)
!  ALLOCATE(N_VolMesh_Ray(iGlobalElem)%Elem_xGP(3,0:Nloc,0:Nloc,0:Nloc))
!  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm(Nloc)%Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iGlobalElem),&
!                     N_VolMesh_Ray(iGlobalElem)%Elem_xGP(:,:,:,:))
!
!END DO
!
!END SUBROUTINE BuildElem_xGP_RayTrace


!===================================================================================================================================
!> Builds the interpolation basis N_Inter_Ray and the Vandermonde matrices PREF_VDM_Ray used for high-order volume sampling for the
!> reay tracing model
!===================================================================================================================================
SUBROUTINE BuildNInterAndVandermonde()
! MODULES
USE MOD_RayTracing_Vars    ,ONLY: Ray,N_Inter_Ray,PREF_VDM_Ray,PerformRayTracing
USE MOD_Interpolation      ,ONLY: InitInterpolationBasis,GetVandermonde
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,Nin,Nout,Nloc
REAL, DIMENSION(:), ALLOCATABLE :: MappedGauss(:)
!===================================================================================================================================
! Allocate interpolation variables
ALLOCATE(N_Inter_Ray(Ray%Nmin:Ray%Nmax))

DO Nloc=Ray%Nmin,Ray%Nmax
  ! Build basis for polynomial of degree Nloc
  CALL InitInterpolationBasis(Nloc , N_Inter_Ray(Nloc)%xGP     , N_Inter_Ray(Nloc)%wGP     , N_Inter_Ray(Nloc)%wBary , &
                                     N_Inter_Ray(Nloc)%L_Minus , N_Inter_Ray(Nloc)%L_Plus  , N_Inter_Ray(Nloc)%L_PlusMinus , &
                                     N_Inter_Ray(Nloc)%swGP    , N_Inter_Ray(Nloc)%wGPSurf , &
                                     N_Inter_Ray(Nloc)%Vdm_Leg , N_Inter_Ray(Nloc)%sVdm_Leg, NodeType_in = Ray%NodeType)

  ! Build variables for nearest Gauss-point (NGP) method
  ALLOCATE(N_Inter_Ray(Nloc)%GaussBorder(1:Nloc))
  ALLOCATE(MappedGauss(1:Nloc+1))

  DO i = 0, Nloc
    MappedGauss(i+1) = N_Inter_Ray(Nloc)%xGP(i)
  END DO ! i = 0, Nloc

  DO i = 1, Nloc
    N_Inter_Ray(Nloc)%GaussBorder(i) = (MappedGauss(i+1) + MappedGauss(i))/2
  END DO ! i = 1, Nloc

  DEALLOCATE(MappedGauss)

END DO

! Only allocate the following arrays when actual ray tracing is performed
IF(PerformRayTracing)THEN
  ! Allocate Vandermonde matrices for p-refinement
  ALLOCATE(PREF_VDM_Ray(Ray%Nmin:Ray%Nmax,Ray%Nmin:Ray%Nmax))

  ! Fill Vandermonde matrices for p-refinement
  DO Nin=Ray%Nmin,Ray%Nmax
    DO Nout=Ray%Nmin,Ray%Nmax
      ALLOCATE(PREF_VDM_Ray(Nin,Nout)%Vdm(0:Nin,0:Nout))
      IF(Nin.EQ.Nout) THEN
        DO i=0,Nin; DO j=0,Nin
          IF(i.EQ.j) THEN
            PREF_VDM_Ray(Nin,Nout)%Vdm(i,j) = 1.
          ELSE
            PREF_VDM_Ray(Nin,Nout)%Vdm(i,j) = 0.
          END IF
        END DO
      END DO
    ELSE IF(Nin.GT.Nout) THEN ! p-coarsening: Project from higher degree to lower degree
      CALL GetVandermonde(Nin, Ray%NodeType, Nout, Ray%NodeType, PREF_VDM_Ray(Nin,Nout)%Vdm, modal=.TRUE. )
    ELSE                   ! p-refinement: Interpolate lower degree to higher degree
      CALL GetVandermonde(Nin, Ray%NodeType, Nout, Ray%NodeType, PREF_VDM_Ray(Nin,Nout)%Vdm, modal=.FALSE.)
    END IF
  END DO;END DO
END IF ! PerformRayTracing

END SUBROUTINE BuildNInterAndVandermonde


!===================================================================================================================================
!> Deallocate arrays
!===================================================================================================================================
SUBROUTINE FinalizeRayTracing()
! MODULES
USE MOD_Globals
USE MOD_RayTracing_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars     ,ONLY: MPI_COMM_SHARED
#endif /*USE_MPI*/
USE MOD_Photon_TrackingVars
USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
USE MOD_Photon_Tracking     ,ONLY: FinalizePhotonSurfSample
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iGlobalElem
!===================================================================================================================================

! Check if ray tracing is used
IF(.NOT.UseRayTracing) RETURN

! Check if actual ray tracing through the domain is performed
IF(PerformRayTracing)THEN

  ! 1: after ray tracing is performed
  SDEALLOCATE(RayElemPassedEnergy)
  DO iGlobalElem = 1, nGlobalElems
    DEALLOCATE(U_N_Ray(iGlobalElem)%U)
  END DO ! iGlobalElem = 1, nGlobalElems
  DEALLOCATE(U_N_Ray)      ! ray tracing
  SDEALLOCATE(PREF_VDM_Ray) ! ray tracing

#if USE_MPI
  SDEALLOCATE(PhotonSampWallProc) ! ray tracing
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  CALL UNLOCK_AND_FREE(RayElemPassedEnergy_Shared_Win)

  IF(nProcessors.GT.1)THEN
    SDEALLOCATE(RayElemOffset)
    CALL UNLOCK_AND_FREE(RayElemPassedEnergyHO_Shared_Win)
  END IF ! nProcessors.GT.1

  CALL UNLOCK_AND_FREE(N_DG_Ray_Shared_Win)

  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

  ADEALLOCATE(RayElemPassedEnergy_Shared)
  ADEALLOCATE(RayElemPassedEnergyHO_Shared)
  ADEALLOCATE(N_DG_Ray_Shared)
#endif /*USE_MPI*/
ELSE
  ! 2: at the end of the simulation or during load balance
  SDEALLOCATE(U_N_Ray_loc)  ! ray tracing + plasma simulation
  SDEALLOCATE(N_DG_Ray_loc) ! ray tracing + plasma simulation
  SDEALLOCATE(N_Inter_Ray)  ! ray tracing + plasma simulation

  ! TODO: see above: deallocate these arrays after simulation end because otherwise these fields will be corrupt in the state file
  ! and that can cause confusion
  SDEALLOCATE(RayElemPassedEnergyLoc1st)
  SDEALLOCATE(RayElemPassedEnergyLoc2nd)
  SDEALLOCATE(RaySecondaryVectorX)
  SDEALLOCATE(RaySecondaryVectorY)
  SDEALLOCATE(RaySecondaryVectorZ)
  SDEALLOCATE(ElemVolume)

  SDEALLOCATE(RayElemEmission)

  SDEALLOCATE(PhotonSampWall_loc)  ! ray tracing + plasma simulation

#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  CALL UNLOCK_AND_FREE(PhotonSampWallHDF5_Shared_Win)
  CALL MPI_BARRIER(MPI_COMM_SHARED,iError)
  ADEALLOCATE(PhotonSampWallHDF5_Shared)
  ADEALLOCATE(PhotonSampWallHDF5)
#endif /*USE_MPI*/
  ! Deallocate surface variables
  CALL FinalizePhotonSurfSample()
END IF ! PerformRayTracing

END SUBROUTINE FinalizeRayTracing

END MODULE MOD_RayTracing_Init
