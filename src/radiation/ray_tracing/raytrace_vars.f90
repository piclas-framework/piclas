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

MODULE MOD_RayTracing_Vars
!===================================================================================================================================
! Contains the tadiation transport variables
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL RAY TRACING VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL :: UseRayTracing        ! Activate ray tracing based emission (also required for plasma simulation)
LOGICAL :: PerformRayTracing    ! Activate actual ray tracing algorithms that track rays through the complete mesh (full mesh mode)

TYPE tRayTrace
  REAL    :: PulseDuration      !<
  REAL    :: tShift             !<
  REAL    :: tActive            !<
  REAL    :: Period             !<
  INTEGER :: NbrOfPulses        !<
  REAL    :: WaistRadius        !<
  REAL    :: WaveLength         !<
  REAL    :: RepetitionRate     !<
  REAL    :: PowerDensity       !<
  REAL    :: Power              !<
  REAL    :: Area               !<
  REAL    :: Energy             !<
  REAL    :: IntensityAmplitude !<
  REAL    :: Direction(3)       !<
  REAL    :: BaseVector1IC(3)   !<
  REAL    :: BaseVector2IC(3)   !<

  ! Output of high-order p-adaptive info
  INTEGER :: NMin               !< Minimum polynomial degree for the high-order volume sampling (p-adaption)
  INTEGER :: NMax               !< Maximum polynomial degree for the high-order volume sampling (p-adaption)

  INTEGER :: VolRefineMode      !< High-order ray tracing volume sampling refinement method:
                                !<  0: do nothing (default)
                                !<  1: refine below user-defined z-coordinate with NMax
                                !<  2: scale N according to the mesh element volume between NMin>=1 and NMax>=2
                                !<  3: refine below user-defined z-coordinate and scale N according to the mesh element volume between NMin>=1 and NMax>=2
                                !<     (consider only elements below the user-defined z-coordinate for the scaling)

  REAL    :: VolRefineModeZ     !< z-coordinate for switching between NMin and NMax

  CHARACTER(LEN=255) :: NodeType  !< equidistant or Gauss nodes [-1,1]

  INTEGER :: nSurfSample        !< polynomial degree of ray tracing or radiation BC sampling

END TYPE

TYPE (tRayTrace)     :: Ray                            !<

TYPE tRadTrans
  INTEGER            :: NumPhotonsPerCell              !<
  REAL               :: GlobalRadiationPower           !<
  REAL               :: ScaledGlobalRadiationPower     !<
  INTEGER            :: GlobalPhotonNum                !<
END TYPE

TYPE (tRadTrans)     :: RadTrans                       !<

LOGICAL              :: RayForceAbsorption             !< Surface photon sampling is performed independent of the actual absorption/reflection outcome (default=T)
LOGICAL, ALLOCATABLE :: RayElemEmission(:,:)           !< Flag elements that are relevant for volume or surface photoionization
INTEGER              :: NumRays                        !<
INTEGER              :: RayPartBound                   !< Particle boundary ID where rays are emitted from

! Output of low-order info
INTEGER,PARAMETER    :: RayElemSize=6
REAL, ALLOCATABLE    :: RayElemPassedEnergy(:,:)       !<
#if USE_MPI
INTEGER              :: RayElemPassedEnergy_Shared_Win    !<
REAL,POINTER         :: RayElemPassedEnergy_Shared(:,:)   !<
INTEGER              :: RayElemPassedEnergyHO_Shared_Win  !< high-order sampling
REAL,POINTER         :: RayElemPassedEnergyHO_Shared(:,:) !< high-order sampling
INTEGER,ALLOCATABLE  :: RayElemOffset(:)                  !< Entry offset for high-order sampling
#endif
REAL, ALLOCATABLE    :: RayElemPassedEnergyLoc1st(:)
REAL, ALLOCATABLE    :: RayElemPassedEnergyLoc2nd(:)
REAL, ALLOCATABLE    :: RaySecondaryVectorX(:)
REAL, ALLOCATABLE    :: RaySecondaryVectorY(:)
REAL, ALLOCATABLE    :: RaySecondaryVectorZ(:)
REAL,ALLOCATABLE     :: ElemVolume(:)

! Output of high-order p-adaptive info
INTEGER,PARAMETER    :: nVarRay=5                      !< Number of variables for higher-order sampling for volume ray tracing

INTEGER,ALLOCATABLE  :: N_DG_Ray_loc(:)                !< for output to ElemData and usage in emission routines
INTEGER,ALLOCPOINT   :: N_DG_Ray(:)                    !< polynomial degree inside DG element for higher-order sampling for volume ray tracing, size(nElems)
#if USE_MPI
INTEGER              :: N_DG_Ray_Shared_Win
INTEGER,ALLOCPOINT   :: N_DG_Ray_Shared(:)
#endif

! DG solution volume
TYPE N_U_Vol
  REAL,ALLOCATABLE  :: U(:,:,:,:)                      !< Polynomial solution of sampled data in each volume element (p-adaptive construct)
END TYPE N_U_Vol

! DG solution (JU or U) vectors
TYPE(N_U_Vol),ALLOCATABLE :: U_N_Ray(:)                !< Solution variable for each equation, node and element,
TYPE(N_U_Vol),ALLOCATABLE :: U_N_Ray_loc(:)            !< Solution variable for each equation, node and element,

!-----------------------------------------------------------------------------------------------------------------------------------
! Volume mesh variables
!-----------------------------------------------------------------------------------------------------------------------------------
!TYPE, PUBLIC :: VolMesh
  !REAL,ALLOCATABLE :: Elem_xGP(:,:,:,:) !< XYZ positions (first index 1:3) of the volume Gauss Point
!END TYPE VolMesh

!TYPE(VolMesh),ALLOCATABLE  :: N_VolMesh_Ray(:) !< Array to store Mesh metrics object "VolMesh"

!-----------------------------------------------------------------------------------------------------------------------------------
! Interpolation variables
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE, PUBLIC :: Interpolation
  ! reserved for Gauss Points with polynomial degree N, all allocated (0:N)
  REAL,ALLOCATABLE :: L_Plus(:)                   !< L for boundary flux computation at plus side  (1)
  REAL,ALLOCATABLE :: L_Minus(:)                  !< L for boundary flux computation at minus side (-1)
  REAL,ALLOCATABLE :: L_PlusMinus(:,:)            !< L for boundary flux computation at both sides (-1,1)
  REAL,ALLOCATABLE :: xGP(:)                      !< Gauss point coordinates
  REAL,ALLOCATABLE :: wGP(:)                      !< GP integration weights
  REAL,ALLOCATABLE :: swGP(:)                     !< 1.0/ GP integration weights
  REAL,ALLOCATABLE :: wBary(:)                    !< barycentric weights
  REAL,ALLOCATABLE :: wGPSurf(:,:)                !< wGPSurf(i,j)=wGP(i)*wGP(j)
  REAL,ALLOCATABLE :: NChooseK(:,:)               !< array n over n
  REAL,ALLOCATABLE :: Vdm_Leg(:,:), sVdm_Leg(:,:) !< Legendre Vandermonde matrix
  REAL,ALLOCATABLE :: GaussBorder(:)              !< Variable required for Nearest Gauss Point (NGP) assignment
END TYPE Interpolation

TYPE(Interpolation),ALLOCATABLE    :: N_Inter_Ray(:)      !< Array of prebuild interpolation matrices

TYPE, PUBLIC :: pVDM
  REAL,ALLOCATABLE   :: Vdm(:,:)                          !< Vandermonde matrix (PP_in,PP_out)
END TYPE pVDM

TYPE(pVDM),ALLOCATABLE             :: PREF_VDM_Ray(:,:)   !< Vandermonde matrices used for p-refinement and coarsening
!===================================================================================================================================
END MODULE MOD_RayTracing_Vars
