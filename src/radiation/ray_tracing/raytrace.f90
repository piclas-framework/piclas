!==================================================================================================================================
! Copyright (c) 2018 - 2019 Marcel Pfeiffer
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

MODULE MOD_RayTracing
!===================================================================================================================================
! Module for the main radiation transport routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE RayTracing
  MODULE PROCEDURE RayTracing
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: RayTracing
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Main routine for the Radiation Transport
!===================================================================================================================================
SUBROUTINE RayTracing()
! MODULES
USE MOD_Globals
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_RayTracing_Vars         ,ONLY: RayPartBound,NumRays,Ray,RayElemPassedEnergy
USE MOD_Photon_TrackingVars     ,ONLY: PhotonProps
USE MOD_Photon_Tracking         ,ONLY: PhotonTriaTracking
USE MOD_Mesh_Tools              ,ONLY: GetGlobalElemID
USE MOD_Output                  ,ONLY: PrintStatusLineRadiation
USE MOD_Mesh_Tools              ,ONLY: GetCNElemID
USE MOD_part_emission_tools     ,ONLY: InsideQuadrilateral
USE MOD_Particle_Boundary_Vars  ,ONLY: nComputeNodeSurfTotalSides,PartBound,SurfSide2GlobalSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Tools ,ONLY: StoreBoundaryParticleProperties
USE MOD_Particle_Boundary_Vars  ,ONLY: PartStateBoundary,nVarPartStateBoundary
USE MOD_Photon_TrackingOutput   ,ONLY: WritePhotonSurfSampleToHDF5,WritePhotonVolSampleToHDF5
#if USE_MPI
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared
USE MOD_Photon_TrackingVars     ,ONLY: PhotonSampWall_Shared, PhotonSampWall_Shared_Win
USE MOD_RayTracing_Vars         ,ONLY: RayElemPassedEnergy_Shared,RayElemPassedEnergy_Shared_Win
#endif /*USE_MPI*/
USE MOD_Photon_TrackingVars     ,ONLY: PhotonSampWall
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: NonUniqueGlobalSideID,iSurfSide,iBC
INTEGER :: CNElemID, iRay, photonCount, RayVisCount, LocRayNum, RayDisp
LOGICAL :: FoundComputeNodeSurfSide
INTEGER :: ALLOCSTAT
REAL    :: RectPower
REAL    :: StartT,EndT                                                                               ! Timer
!===================================================================================================================================

IF(RayPartBound.EQ.0) RETURN
GETTIME(StartT)
SWRITE(UNIT_stdOut,'(A)') ' Start Ray Tracing Calculation ...'

! Allocate arrays
ALLOCATE(RayElemPassedEnergy(1:nGlobalElems))
RayElemPassedEnergy=0.0
ALLOCATE(PhotonSampWall(2,1:nComputeNodeSurfTotalSides))
PhotonSampWall=0.0

#if USE_MPI
!> Shared arrays for volume sampling
CALL Allocate_Shared((/nGlobalElems/),RayElemPassedEnergy_Shared_Win,RayElemPassedEnergy_Shared)
CALL MPI_WIN_LOCK_ALL(0,RayElemPassedEnergy_Shared_Win,IERROR)
!> Shared arrays for boundary sampling
CALL Allocate_Shared((/2,nComputeNodeSurfTotalSides/),PhotonSampWall_Shared_Win,PhotonSampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,PhotonSampWall_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) RayElemPassedEnergy_Shared = 0.
IF (myComputeNodeRank.EQ.0) PhotonSampWall_Shared = 0.
CALL BARRIER_AND_SYNC(RayElemPassedEnergy_Shared_Win,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(PhotonSampWall_Shared_Win,MPI_COMM_SHARED)
#endif

photonCount = 0
RayVisCount = 0
IF(nProcessors.GT.NumRays) CALL abort(__STAMP__,'Use more rays! Number of processes > Number of rays')
LocRayNum = NumRays/nProcessors
IF(myrank.LT.MOD(NumRays,nProcessors)) LocRayNum = LocRayNum + 1
! Output to screen every 20 rays to show that the tool is still running
RayDisp = MAX(1,INT(LocRayNum/20)) ! This value cannot be zero

RectPower = Ray%IntensityAmplitude * Ray%Area / REAL(NumRays)

! This array is not de-allocated during load balance as it is only written to .h5 during WriteStateToHDF5()
IF(.NOT.ALLOCATED(PartStateBoundary))THEN
  ALLOCATE(PartStateBoundary(1:nVarPartStateBoundary,1:LocRayNum), STAT=ALLOCSTAT)
  IF (ALLOCSTAT.NE.0) CALL abort(__STAMP__,'ERROR in particle_init.f90: Cannot allocate PartStateBoundary array!')
  PartStateBoundary=0.
END IF ! .NOT.ALLOCATED(PartStateBoundary)

DO iRay = 1, LocRayNum
  IF(MPIroot.AND.(MOD(RayVisCount,RayDisp).EQ.0)) CALL PrintStatusLineRadiation(REAL(RayVisCount),0.0,REAL(LocRayNum),.TRUE.)
  RayVisCount = RayVisCount + 1
  PhotonProps%PhotonPos(1:3) = SetRayPos()
  PhotonProps%PhotonLastPos(1:3) = PhotonProps%PhotonPos(1:3)
  PhotonProps%PhotonDirection(1:3) = Ray%Direction ! (/0.0,0.0,-1.0/)! SetPhotonStartDirection(iCNElem, iPhotLoc, RandRot)

  ! Loop over all sides of a specific iPartBoundary and find the side where the ray enters the domain
  ! count number of nSides connected to iPartBoundary BCSideID

  PhotonProps%ElemID = -1 ! Initialize

  FoundComputeNodeSurfSide = .FALSE.
  SurfLoop: DO iSurfSide = 1,nComputeNodeSurfTotalSides
    NonUniqueGlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSurfSide)
    ! Check if the surface side has a neighbor (and is therefore an inner BCs)
    IF(SideInfo_Shared(SIDE_NBSIDEID,NonUniqueGlobalSideID).LE.0) THEN ! BC side
      ! Get field BC index of side and compare with BC index of the corresponding particle boundary index of the emission side
      iBC = SideInfo_Shared(SIDE_BCID,NonUniqueGlobalSideID) ! Get field BC index from non-unique global side index
      IF(iBC.NE.PartBound%MapToFieldBC(RayPartBound)) CYCLE SurfLoop ! Correct BC not found, check next side
      FoundComputeNodeSurfSide = .TRUE.
      ! Check if ray starting position is within the quadrilateral that is spanned by the four corner nodes of the side
      IF(InsideQuadrilateral(PhotonProps%PhotonPos(1:2),NonUniqueGlobalSideID))THEN
        ! Found CN element index
        CNElemID = GetCNElemID(SideInfo_Shared(SIDE_ELEMID,NonUniqueGlobalSideID))
        ! Set global element index for current ray
        PhotonProps%ElemID = GetGlobalElemID(CNElemID)
        PhotonProps%PhotonEnergy = RectPower
        EXIT SurfLoop
      END IF ! InsideQuadrilateral(X,NonUniqueGlobalSideID)
    END IF ! SideInfo_Shared(SIDE_NBSIDEID,NonUniqueGlobalSideID).LE.0
  END DO SurfLoop! iSurfSide = 1,nComputeNodeSurfTotalSides

  ! Sanity check: nComputeNodeSurfTotalSides > 0 and the correct PartBCIndex for those sides
  IF(.NOT.FoundComputeNodeSurfSide)THEN
    IPWRITE(UNIT_StdOut,*) ": nComputeNodeSurfTotalSides =", nComputeNodeSurfTotalSides
    IPWRITE(UNIT_StdOut,*) ": RayPartBound               =", RayPartBound
    !IPWRITE(UNIT_StdOut,'(I0,A,I0,A)') ": Set Part-Boundary",RayPartBound,"-BoundaryParticleOutput = T"
    IPWRITE(UNIT_StdOut,*) ": Set Particles-DSMC-CalcSurfaceVal = T"
    CALL abort(__STAMP__,'No boundary found in list of nComputeNodeSurfTotalSides for defined RayPartBound!')
  END IF ! FoundComputeNodeSurfSide

  ! Sanity check
  IF(PhotonProps%ElemID.LE.0)THEN
    IPWRITE(UNIT_StdOut,*) "PhotonProps%PhotonPos(1:3) =", PhotonProps%PhotonPos(1:3)
    CALL abort(__STAMP__,'Ray starting element not found!')
  ELSE
    ! Output ray starting position and direction vector to .h5 for debugging
    CALL StoreBoundaryParticleProperties(iRay,&
         999,&
         PhotonProps%PhotonPos(1:3),&
         UNITVECTOR(PhotonProps%PhotonDirection(1:3)),(/0.0,0.0,1.0/),&
         iPartBound=RayPartBound,&
         mode=2,&
         MPF_optIN=0.0,&
         Velo_optIN=PhotonProps%PhotonDirection(1:3))
  END IF !PhotonProps%ElemID.LE.0

  CALL PhotonTriaTracking()
END DO

! Print 100%
IF(MPIroot)THEN
  CALL PrintStatusLineRadiation(REAL(RayVisCount),0.0,REAL(LocRayNum),.TRUE.)
  WRITE(UNIT_StdOut,'(A)') " "
END IF

CALL WritePhotonSurfSampleToHDF5()

CALL WritePhotonVolSampleToHDF5()

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, ' DONE!', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

END SUBROUTINE RayTracing


FUNCTION SetPhotonEnergy(iCNElem, Point, iWave)
!===================================================================================================================================
!> Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: Pi
USE MOD_RadiationTrans_Vars ,ONLY: RadEmiAdaptPhotonNum, Radiation_Emission_Spec_Total, RadTrans, RadTransPhotPerCell
USE MOD_RadiationTrans_Vars ,ONLY: RadObservationPoint, RadObservationPointMethod,RadTransObsVolumeFrac,RadObservationPOI
USE MOD_Particle_Mesh_Vars  ,ONLY: ElemVolume_Shared
USE MOD_Radiation_Vars      ,ONLY: RadiationParameter,Radiation_Emission_spec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)       :: iCNElem
REAL, INTENT(IN)          :: Point(3)
INTEGER, INTENT(IN), OPTIONAL :: iWave
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                      :: SetPhotonEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL              :: ProjectedDist(3), Dist(3), ClosestPoint(3), FarthestPoint(3), Vec1(3), Vec2(3), fullangle
REAL               :: cosTheta, Dist(3), DistNorm(3), spaceangle, absdistnorm
!===================================================================================================================================
IF (RadEmiAdaptPhotonNum) THEN
  SetPhotonEnergy = Radiation_Emission_Spec_Total(iCNElem)*ElemVolume_Shared(iCNElem)*RadTransObsVolumeFrac(iCNElem) / RadTransPhotPerCell(iCNElem)
ELSE
  SetPhotonEnergy = Radiation_Emission_Spec_Total(iCNElem)*ElemVolume_Shared(iCNElem)*RadTransObsVolumeFrac(iCNElem) / (RadTrans%NumPhotonsPerCell)
END IF

IF (RadObservationPointMethod.EQ.1) THEN
   Dist(1:3) = Point(1:3) - RadObservationPoint%MidPoint(1:3)
   absdistnorm = VECNORM(Dist(1:3))
   DistNorm(1:3) = Dist(1:3)/absdistnorm
   cosTheta = DOT_PRODUCT(RadObservationPoint%ViewDirection(1:3),DistNorm(1:3))/(VECNORM(RadObservationPoint%ViewDirection(1:3))*VECNORM(DistNorm(1:3)))
   spaceangle = cosTheta * RadObservationPoint%Area/(absdistnorm*absdistnorm)
!  ProjectedDist(1:3) = Dist(1:3) - DOT_PRODUCT(RadObservationPoint%ViewDirection(1:3),Dist(1:3))*RadObservationPoint%ViewDirection(1:3)
!  ClosestPoint(1:3) = RadObservationPoint%MidPoint(1:3) + RadObservationPoint%Diameter/2.*ProjectedDist(1:3)/VECNORM(ProjectedDist(1:3))
!  FarthestPoint(1:3) = RadObservationPoint%MidPoint(1:3) - RadObservationPoint%Diameter/2.*ProjectedDist(1:3)/VECNORM(ProjectedDist(1:3))
!  Vec1(1:3) = ClosestPoint(1:3) - Point(1:3)
!  Vec2(1:3) = FarthestPoint(1:3) - Point(1:3)
!  fullangle = ACOS(DOT_PRODUCT(Vec1,Vec2)/(VECNORM(Vec1)*VECNORM(Vec2)))
  SetPhotonEnergy = SetPhotonEnergy * spaceangle/(4.*Pi)
ELSEIF (RadObservationPointMethod.EQ.2) THEN
  IF (RadObservationPoint%CalcFullSpectra) THEN
    SetPhotonEnergy = Radiation_Emission_Spec(iWave, iCNElem) * RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor &
        *ElemVolume_Shared(iCNElem)*RadTransObsVolumeFrac(iCNElem)
  ELSE
    SetPhotonEnergy = SetPhotonEnergy /(4.*Pi)
  END IF
  SetPhotonEnergy = SetPhotonEnergy / (ElemVolume_Shared(iCNElem)*RadTransObsVolumeFrac(iCNElem))*RadObservationPOI(7, iCNElem)
END IF

END FUNCTION SetPhotonEnergy


FUNCTION SetRayPos()
!===================================================================================================================================
!> Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
USE MOD_Particle_Mesh_Vars   ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: SetRayPos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: RandVal(2)
!===================================================================================================================================
CALL RANDOM_NUMBER(RandVal)
SetRayPos = (/RandVal(1)*(GEO%xmaxglob-GEO%xminglob)+GEO%xminglob,&
              RandVal(2)*(GEO%ymaxglob-GEO%yminglob)+GEO%yminglob,&
              GEO%zmaxglob/)
END FUNCTION SetRayPos

FUNCTION SetPhotonStartDirection(iCNElem, iPhot, RandRot)
!===================================================================================================================================
! modified particle emmission for LD case
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars        ,ONLY: Pi
USE MOD_RadiationTrans_Vars ,ONLY: RadiationDirectionModel, RadTransPhotPerCell, RadObservationPointMethod,RadObservationPoint
USE MOD_Photon_TrackingVars ,ONLY: PhotonProps
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iCNElem, iPhot
REAL, INTENT(IN)                :: RandRot(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                             :: SetPhotonStartDirection(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: iRan,RandomDirection(2), X_new, Y_new, start, incr, SpiralPos, SpiralStep
INTEGER                          :: RadMod
!===================================================================================================================================
SELECT CASE(RadiationDirectionModel)
CASE(1)
  RadMod = RadiationDirectionModel
CASE(2)
  IF (RadTransPhotPerCell(iCNElem).EQ.1) THEN
    RadMod = 1
  ELSE
    RadMod = RadiationDirectionModel
  END IF
CASE DEFAULT
  CALL abort(&
  __STAMP__&
  ,' ERROR: Radiation-DirectionModel not implemented!. (unknown case)')
END SELECT !PartBound%MapToPartBC(BC(SideID)
IF (RadObservationPointMethod.EQ.1) THEN
  CALL RANDOM_NUMBER(iRan)
  RandomDirection(1)  = RadObservationPoint%Diameter/2. * SQRT(iRan)
  CALL RANDOM_NUMBER(iRan)
  RandomDirection(2) = iRan * 2. * Pi
  SetPhotonStartDirection(1) = 0.0
  SetPhotonStartDirection(2) = RandomDirection(1) * COS(RandomDirection(2))
  SetPhotonStartDirection(3) = RandomDirection(1) * SIN(RandomDirection(2))
  SetPhotonStartDirection(1:3) = MATMUL(RadObservationPoint%OrthoNormBasis, SetPhotonStartDirection(1:3))
  SetPhotonStartDirection(1:3) = SetPhotonStartDirection(1:3) + RadObservationPoint%MidPoint(1:3)
  SetPhotonStartDirection(1:3) = SetPhotonStartDirection(1:3) - PhotonProps%PhotonPos(1:3)
  SetPhotonStartDirection(1:3) = SetPhotonStartDirection(1:3) / VECNORM(SetPhotonStartDirection(1:3))
ELSEIF (RadObservationPointMethod.EQ.2) THEN
!    SetPhotonStartDirection(1:3) = RadObservationPoint%MidPoint(1:3)
!    SetPhotonStartDirection(1:3) = SetPhotonStartDirection(1:3) - RadObservationPoint%ViewDirection(1:3)
  SetPhotonStartDirection(1:3) = -RadObservationPoint%ViewDirection(1:3)
  SetPhotonStartDirection(1:3) = SetPhotonStartDirection(1:3) / VECNORM(SetPhotonStartDirection(1:3))
ELSE
  SELECT CASE(RadMod)
  CASE(1)
    CALL RANDOM_NUMBER(iRan)
    RandomDirection(1) = 2.*iRan - 1.
    CALL RANDOM_NUMBER(iRan)
    RandomDirection(2) = 2.*Pi*iRan - Pi
    SetPhotonStartDirection(1)  = SIN(RandomDirection(2))*SQRT(1.-RandomDirection(1)**2.)
    SetPhotonStartDirection(2)  = COS(RandomDirection(2))*SQRT(1.-RandomDirection(1)**2.)
    SetPhotonStartDirection(3)  = RandomDirection(1)
  CASE(2)
    SpiralStep = 0.1+1.2*REAL(RadTransPhotPerCell(iCNElem))
    start = (-1. + 1./(REAL(RadTransPhotPerCell(iCNElem))-1.))
    incr = (2.-2./(REAL(RadTransPhotPerCell(iCNElem))-1.))/(REAL(RadTransPhotPerCell(iCNElem))-1.)
    SpiralPos = start + (REAL(iPhot)-1.) *incr
    X_new = SpiralPos * SpiralStep
    Y_new = Pi/2.*SIGN(1.,SpiralPos)*(1.-SQRT(1.-ABS(SpiralPos)))
    SetPhotonStartDirection(1)  = COS(X_new)*COS(Y_new)
    SetPhotonStartDirection(2)  = SIN(X_new)*COS(Y_new)
    SetPhotonStartDirection(3)  = SIN(Y_new)
    SetPhotonStartDirection(1:3)  = MATMUL(RandRot, SetPhotonStartDirection(1:3))
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,' ERROR: Radiation-DirectionModel not implemented!. (unknown case)')
  END SELECT !PartBound%MapToPartBC(BC(SideID)
END IF

END FUNCTION SetPhotonStartDirection

FUNCTION RandomRotMatrix()
!===================================================================================================================================
! Calculation of the vibrational temperature (zero-point search) for the TSHO (Truncated Simple Harmonic Oscillator)
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,         ONLY : Pi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                    :: RandomRotMatrix(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: alpha(3) , A(3,3)
!===================================================================================================================================
CALL RANDOM_NUMBER(alpha)
alpha(1:3) = 2.*alpha(1:3)*Pi
RandomRotMatrix = RESHAPE((/1.,0.,0.,0.,COS(alpha(1)),SIN(alpha(1)),0.,-SIN(alpha(1)), COS(alpha(1))/),(/3,3/))
A = RESHAPE((/COS(alpha(2)),0.,-SIN(alpha(2)),0.,1.,0.,SIN(alpha(2)),0.0, COS(alpha(2))/),(/3,3/))
RandomRotMatrix = MATMUL(A,RandomRotMatrix)
A = RESHAPE((/COS(alpha(3)),SIN(alpha(3)),0.,-SIN(alpha(3)),COS(alpha(3)),0.,0.,0.0, 1./),(/3,3/))
RandomRotMatrix = MATMUL(A, RandomRotMatrix)

END FUNCTION RandomRotMatrix


FUNCTION SetParticleWavelengthAR(iCNElem)
!===================================================================================================================================
! modified particle emmission for LD case
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : Pi
USE MOD_RadiationTrans_Vars,   ONLY : Radiation_Emission_Spec_Max
USE MOD_Radiation_Vars,        ONLY: Radiation_Emission_spec, RadiationParameter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iCNElem
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                         :: SetParticleWavelengthAR
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iWaveLength
REAL                             :: iRan, iRadPower
!===================================================================================================================================

CALL RANDOM_NUMBER(iRan)
iWaveLength = INT(RadiationParameter%WaveLenDiscrCoarse*iRan) + 1
IF ((RadiationParameter%WaveLenReductionFactor.GT.1).AND.(iWaveLength.EQ.RadiationParameter%WaveLenDiscrCoarse)) THEN
  IF (MOD(RadiationParameter%WaveLenDiscr,RadiationParameter%WaveLenDiscrCoarse).NE.0) THEN
    iRadPower = 4.*Pi*Radiation_Emission_Spec(RadiationParameter%WaveLenDiscrCoarse, iCNElem) * RadiationParameter%WaveLenIncr &
       * (1.+RadiationParameter%WaveLenReductionFactor)
  ELSE
    iRadPower = 4.*Pi*Radiation_Emission_Spec(iWaveLength,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
  END IF
ELSE
  iRadPower = 4.*Pi*Radiation_Emission_Spec(iWaveLength,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
END IF
CALL RANDOM_NUMBER(iRan)
DO WHILE (iRan.GT.(iRadPower/Radiation_Emission_Spec_Max(iCNElem)))
  CALL RANDOM_NUMBER(iRan)
  iWaveLength = INT(RadiationParameter%WaveLenDiscrCoarse*iRan) + 1
  IF ((RadiationParameter%WaveLenReductionFactor.GT.1).AND.(iWaveLength.EQ.RadiationParameter%WaveLenDiscrCoarse)) THEN
    IF (MOD(RadiationParameter%WaveLenDiscr,RadiationParameter%WaveLenDiscrCoarse).NE.0) THEN
      iRadPower = 4.*Pi*Radiation_Emission_Spec(RadiationParameter%WaveLenDiscrCoarse, iCNElem) * RadiationParameter%WaveLenIncr &
         * (1.+RadiationParameter%WaveLenReductionFactor)
    ELSE
      iRadPower = 4.*Pi*Radiation_Emission_Spec(iWaveLength,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
    END IF
  ELSE
    iRadPower = 4.*Pi*Radiation_Emission_Spec(iWaveLength,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
  END IF
  CALL RANDOM_NUMBER(iRan)
END DO
SetParticleWavelengthAR = iWaveLength

END FUNCTION SetParticleWavelengthAR


FUNCTION SetParticleWavelengthBiSec(iCNElem)
!===================================================================================================================================
! modified particle emmission for LD case
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,          ONLY : Pi
USE MOD_RadiationTrans_Vars,   ONLY : Radiation_Emission_Spec_Total
USE MOD_Radiation_Vars,        ONLY: Radiation_Emission_spec, RadiationParameter
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iCNElem
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                         :: SetParticleWavelengthBiSec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iWaveLength, iWave, iWaveOld, iWaveMin, iWaveMax
REAL                             :: iRan, iRadPower, iRadPower2
!===================================================================================================================================

CALL RANDOM_NUMBER(iRan)
iWaveOld = 1
iWaveLength = INT(RadiationParameter%WaveLenDiscrCoarse/2)
iWaveMin = 1
iWaveMax = RadiationParameter%WaveLenDiscrCoarse
IF (iWaveLength.EQ.RadiationParameter%WaveLenDiscrCoarse) THEN
  iRadPower = Radiation_Emission_Spec_Total(iCNElem)
ELSE
  iRadPower = 0.0
  DO iWave = 1,  iWaveLength
    iRadPower = iRadPower + 4.*Pi*Radiation_Emission_Spec(iWave,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
  END DO
END IF

DO
  IF (iRan.GT.(iRadPower/Radiation_Emission_Spec_Total(iCNElem)))THEN
    iWaveMin = iWaveLength
  ELSE
    iWaveMax = iWaveLength
  END IF
  iWaveOld = iWaveLength
  iWaveLength = INT((iWaveMax+iWaveMin)/2)
  IF (iWaveLength.EQ.RadiationParameter%WaveLenDiscrCoarse) THEN
    iRadPower = Radiation_Emission_Spec_Total(iCNElem)
  ELSE
    iRadPower = 0.0
    DO iWave = 1,  iWaveLength
      iRadPower = iRadPower + 4.*Pi*Radiation_Emission_Spec(iWave,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
    END DO
  END IF
  IF (ABS(iWaveOld-iWaveLength).LE.1) EXIT

END DO

iWaveOld = iWaveLength
IF (iRan.LT.(iRadPower/Radiation_Emission_Spec_Total(iCNElem))) THEN
  IF (iWaveLength.EQ.1) THEN
    iWaveLength = iWaveLength
  ELSE
    iWaveLength = iWaveLength - 1
    iRadPower2 = 0.0
    DO iWave = 1,  iWaveLength
      iRadPower2 = iRadPower2 + 4.*Pi*Radiation_Emission_Spec(iWave,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
    END DO
    IF (ABS(iRan-(iRadPower/Radiation_Emission_Spec_Total(iCNElem))).LT.ABS(iRan-(iRadPower2/Radiation_Emission_Spec_Total(iCNElem)))) THEN
      iWaveLength = iWaveOld
    END IF
  END IF
ELSE
  iWaveLength = iWaveLength + 1
  iRadPower2 = 0.0
  DO iWave = 1,  iWaveLength
    iRadPower2 = iRadPower2 + 4.*Pi*Radiation_Emission_Spec(iWave,iCNElem)*RadiationParameter%WaveLenIncr*RadiationParameter%WaveLenReductionFactor
  END DO
  IF (ABS(iRan-(iRadPower/Radiation_Emission_Spec_Total(iCNElem))).LT.ABS(iRan-(iRadPower2/Radiation_Emission_Spec_Total(iCNElem)))) THEN
    iWaveLength = iWaveOld
  END IF
END IF
SetParticleWavelengthBiSec = iWaveLength

END FUNCTION SetParticleWavelengthBiSec

END MODULE MOD_RayTracing
