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

CALL prms%CreateLogicalOption(   'RayTracing-AdaptiveRays'   , 'TODO' , '.FALSE.')
CALL prms%CreateIntOption(       'RayTracing-NumRays'        , 'TODO' , '1')
CALL prms%CreateIntOption(       'RayTracing-RayPosModel'    , 'TODO' , '1')
CALL prms%CreateRealArrayOption( 'RayTracing-RayDirection'   , 'Direction vector for ray emission. Will be normalized after read-in.' , no=3)
CALL prms%CreateIntOption(       'RayTracing-PartBound'      , 'Particle boundary ID where rays are emitted from' , '0')
CALL prms%CreateRealOption(      'RayTracing-PulseDuration'  , 'Pulse duration tau for a Gaussian-type pulse with I~exp(-(t/tau)^2) [s]'                  )
CALL prms%CreateIntOption(       'RayTracing-NbrOfPulses'    , 'Number of pulses [-]'                                                                     , '1')
CALL prms%CreateRealOption(      'RayTracing-WaistRadius'    , 'Beam waist radius (in focal spot) w_b for Gaussian-type pulse with I~exp(-(r/w_b)^2) [m]' , '0.0')
CALL prms%CreateRealOption(      'RayTracing-WaveLength'     , 'Beam wavelength [m]'                                                                      )
CALL prms%CreateRealOption(      'RayTracing-RepetitionRate' , 'Pulse repetition rate (pulses per second) [Hz]'                                           )
CALL prms%CreateRealOption(      'RayTracing-Power'          , 'Average pulse power (energy of a single pulse times repetition rate) [W]'                 )
CALL prms%CreateLogicalOption(   'RayTracing-ForceAbsorption', 'Surface photon sampling is performed independent of the actual absorption/reflection outcome (default=T)', '.TRUE.')

CALL prms%CreateIntOption(      'RayTracing-NMax'            , 'Maximum polynomial degree within refined volume elements for photon tracking (p-adaption)')

END SUBROUTINE DefineParametersRayTracing


SUBROUTINE InitRayTracing()
!===================================================================================================================================
! Initialization of the radiation transport solver 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools
USE MOD_RayTracing_Vars
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_RadiationTrans_Vars    ,ONLY: RadiationAbsorptionModel,RadObservationPointMethod
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
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RAY TRACING SOLVER ...'

! Do not absorb rays within the volume!
RadiationAbsorptionModel = 0
RadObservationPointMethod = 0

RayPartBound = GETINT('RayTracing-PartBound')
IF(RayPartBound.EQ.0) RETURN
IF(RayPartBound.LT.0) CALL CollectiveStop(__STAMP__,'RayTracing-PartBound must be > 0 to activate ray tracing on this boundary!')

! Get ray parameters
Ray%PulseDuration  = GETREAL('RayTracing-PulseDuration')
Ray%NbrOfPulses    = GETINT('RayTracing-NbrOfPulses')
Ray%tShift         = SQRT(8.0) * Ray%PulseDuration
Ray%WaistRadius    = GETREAL('RayTracing-WaistRadius')
Ray%WaveLength     = GETREAL('RayTracing-WaveLength')
Ray%RepetitionRate = GETREAL('RayTracing-RepetitionRate')
Ray%Period         = 1./Ray%RepetitionRate
Ray%Power          = GETREAL('RayTracing-Power')
Ray%Direction      = GETREALARRAY('RayTracing-RayDirection',3)
Ray%Direction      = UNITVECTOR(Ray%Direction)

AdaptiveRays       = GETLOGICAL('RayTracing-AdaptiveRays')
NumRays            = GETINT('RayTracing-NumRays')
RayPosModel        = GETINT('RayTracing-RayPosModel')
RayForceAbsorption = GETLOGICAL('RayTracing-ForceAbsorption')

! Output of high-order p-adaptive info
Ray%NMin = 1 ! GETINT('RayTracing-NMin')
WRITE(UNIT=hilf,FMT='(I3)') PP_N
Ray%Nmax = GETINT('RayTracing-Nmax',hilf)

! Build all mappings
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
  ! Derived quantities
  E0 = Ray%Power / Ray%RepetitionRate

  ! Rectangle
  ! TODO: Ray emission area from chosen boundary surface?
  A = (GEO%xmaxglob-GEO%xminglob) * (GEO%ymaxglob-GEO%yminglob)
  ! Normal vector of the ray emission area
  SurfaceNormal = (/ 0., 0., 1. /)
  ! Angle between emitted rays and emission area
  alpha = (90.-ABS(90.-(180./PI)*ACOS(DOT_PRODUCT(Ray%Direction,SurfaceNormal))))

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
CALL PrintOption('Single pulse energy [J]'                                           , 'CALCUL.' , RealOpt=Ray%Energy)
CALL PrintOption('Intensity amplitude: I0 [W/m^2]'                                   , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Corrected Intensity amplitude: I0_corr [W/m^2]'                    , 'CALCUL.' , RealOpt=Ray%IntensityAmplitude)
CALL PrintOption('Pulse period (Time between maximum of two pulses) [s]'             , 'CALCUL.' , RealOpt=Ray%Period)
CALL PrintOption('Temporal pulse width (pulse time 2x tShift) [s]'                   , 'CALCUL.' , RealOpt=2.0*Ray%tShift)
CALL PrintOption('Pulse will end at tActive (pulse final time) [s]'                  , 'CALCUL.' , RealOpt=Ray%tActive)

SWRITE(UNIT_stdOut,'(A)')' INIT RAY TRACING SOLVER DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitRayTracing


!===================================================================================================================================
!> Build all high-order mappings required for ray trace sampling in the volume on a p-adaptive polynomial basis
!===================================================================================================================================
SUBROUTINE InitHighOrderRaySampling()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars       ,ONLY: NodeCoords,nElems,ElemBaryNGeo
USE MOD_RayTracing_Vars ,ONLY: N_VolMesh_Ray,N_DG_Ray,Ray,N_Inter_Ray,PREF_VDM_Ray,U_N_Ray,nVarRay
USE MOD_Mesh_Tools      ,ONLY: GetCNElemID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: Nloc,iElem,CNElemID
LOGICAL,PARAMETER :: debugRay=.FALSE.
!===================================================================================================================================
ALLOCATE(N_DG_Ray(nElems))
N_DG_Ray = PP_N
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

! Allocate interpolation variables
ALLOCATE(N_Inter_Ray(Ray%Nmin:Ray%Nmax))
! Allocate Vandermonde matrices for p-refinement
ALLOCATE(PREF_VDM_Ray(Ray%Nmin:Ray%Nmax,Ray%Nmin:Ray%Nmax))
CALL BuildNInterAndVandermonde()

ALLOCATE(N_VolMesh_Ray(1:nElems))
CALL BuildElem_xGP_RayTrace(NodeCoords)

! the local DG solution in physical and reference space
ALLOCATE(U_N_Ray(1:PP_nElems))
DO iElem = 1, PP_nElems
  Nloc = N_DG_Ray(iElem)
  ALLOCATE(U_N_Ray(iElem)%U(nVarRay,0:Nloc,0:Nloc,0:Nloc))
  U_N_Ray(iElem)%U = 0.
END DO ! iElem = 1, PP_nElems

END SUBROUTINE InitHighOrderRaySampling


!==================================================================================================================================
!> This routine takes the equidistant node coordinates of the mesh (on NGeo+1 points) and uses them to build the coordinates
!> of solution/interpolation points of type NodeType on polynomial degree Nloc (Nloc+1 points per direction).
!> Output: Elem_xGP(:,:,:,:) for each element with variably N
!==================================================================================================================================
SUBROUTINE BuildElem_xGP_RayTrace(NodeCoords)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,nElems
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU,NodeType
USE MOD_RayTracing_Vars    ,ONLY: Ray,N_VolMesh_Ray,N_DG_Ray,N_Inter_Ray
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ, ChangeBasis3D
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)               :: NodeCoords(3,0:NGeo,0:NGeo,0:NGeo,nElems)         !< Equidistant mesh coordinates
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,Nloc,i

TYPE VdmType
  REAL, ALLOCATABLE           :: Vdm_EQNGeo_CLNloc(:,:)
  REAL, ALLOCATABLE           :: Vdm_CLNloc_Nloc  (:,:)
END TYPE VdmType

TYPE(VdmType), DIMENSION(:), ALLOCATABLE :: Vdm

REAL, DIMENSION(:), ALLOCATABLE :: MappedGauss(:)
!==================================================================================================================================

! Build Vdm for every degree
ALLOCATE(Vdm(Ray%Nmin:Ray%Nmax))
DO Nloc = Ray%Nmin, Ray%Nmax
  ALLOCATE(Vdm(Nloc)%Vdm_EQNGeo_CLNloc(0:Nloc,0:NGeo))
  ALLOCATE(Vdm(Nloc)%Vdm_CLNloc_Nloc(0:Nloc,0:Nloc))
  CALL GetVandermonde(NGeo, NodeTypeVISU, NLoc, NodeTypeCL, Vdm(Nloc)%Vdm_EQNGeo_CLNloc,  modal=.FALSE.)
  CALL GetVandermonde(Nloc, NodeTypeCL  , Nloc, NodeType  , Vdm(Nloc)%Vdm_CLNloc_Nloc,     modal=.FALSE.)

  ! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
  !       Important for curved meshes if NGeo<N, no effect for N>=NGeo

  !1.a) Transform from EQUI_NGeo to solution points on Nloc
  Vdm(Nloc)%Vdm_EQNGeo_CLNloc=MATMUL(Vdm(Nloc)%Vdm_CLNloc_Nloc, Vdm(Nloc)%Vdm_EQNGeo_CLNloc)
END DO ! Nloc = Ray%Nmin, Ray%Nmax

! Set Elem_xGP for each element
DO iElem=1,nElems
  Nloc = N_DG_Ray(iElem)

  ALLOCATE(N_VolMesh_Ray(iElem)%Elem_xGP(3,0:Nloc,0:Nloc,0:Nloc))
  CALL ChangeBasis3D(3,NGeo,Nloc,Vdm(Nloc)%Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),N_VolMesh_Ray(iElem)%Elem_xGP(:,:,:,:))

  ! Build variables for nearest Gauss-point (NGP) method
  ALLOCATE(N_VolMesh_Ray(iElem)%GaussBorder(1:Nloc))
  ALLOCATE(MappedGauss(1:Nloc+1))

  DO i = 0, Nloc
    MappedGauss(i+1) = N_Inter_Ray(Nloc)%xGP(i)
  END DO ! i = 0, Nloc

  DO i = 1, Nloc
    N_VolMesh_Ray(iElem)%GaussBorder(i) = (MappedGauss(i+1) + MappedGauss(i))/2
  END DO ! i = 1, Nloc

  DEALLOCATE(MappedGauss)

END DO

END SUBROUTINE BuildElem_xGP_RayTrace


!===================================================================================================================================
!> Builds the interpolation basis N_Inter_Ray and the Vandermonde matrices PREF_VDM_Ray used for high-order volume sampling for the
!> reay tracing model
!===================================================================================================================================
SUBROUTINE BuildNInterAndVandermonde()
! MODULES
USE MOD_RayTracing_Vars    ,ONLY: Ray,N_Inter_Ray,PREF_VDM_Ray
USE MOD_Interpolation      ,ONLY: InitInterpolationBasis,GetVandermonde
USE MOD_Interpolation_Vars ,ONLY: NodeType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j,Nin,Nout,Nloc
!===================================================================================================================================
DO Nloc=Ray%Nmin,Ray%Nmax
  CALL InitInterpolationBasis(Nloc , N_Inter_Ray(Nloc)%xGP     , N_Inter_Ray(Nloc)%wGP     , N_Inter_Ray(Nloc)%wBary , &
                                     N_Inter_Ray(Nloc)%L_Minus , N_Inter_Ray(Nloc)%L_Plus  , N_Inter_Ray(Nloc)%L_PlusMinus , &
                                     N_Inter_Ray(Nloc)%swGP    , N_Inter_Ray(Nloc)%wGPSurf , &
                                     N_Inter_Ray(Nloc)%Vdm_Leg , N_Inter_Ray(Nloc)%sVdm_Leg)
END DO

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
    CALL GetVandermonde(Nin, NodeType, Nout, NodeType, PREF_VDM_Ray(Nin,Nout)%Vdm, modal=.TRUE. )
  ELSE                   ! p-refinement: Interpolate lower degree to higher degree
    CALL GetVandermonde(Nin, NodeType, Nout, NodeType, PREF_VDM_Ray(Nin,Nout)%Vdm, modal=.FALSE.)
  END IF
END DO;END DO

END SUBROUTINE BuildNInterAndVandermonde


!===================================================================================================================================
!> Deallocate arrays
!===================================================================================================================================
SUBROUTINE FinalizeRayTracing()
! MODULES
USE MOD_RayTracing_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(N_DG_Ray)
SDEALLOCATE(N_VolMesh_Ray)
SDEALLOCATE(N_Inter_Ray)
SDEALLOCATE(PREF_VDM_Ray)
SDEALLOCATE(U_N_Ray)
END SUBROUTINE FinalizeRayTracing

END MODULE MOD_RayTracing_Init
