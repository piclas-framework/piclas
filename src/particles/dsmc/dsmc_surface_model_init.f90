#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModelInit
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE InitDSMCSurfModel
  MODULE PROCEDURE InitDSMCSurfModel
END INTERFACE

INTERFACE FinalizeDSMCSurfModel
  MODULE PROCEDURE FinalizeDSMCSurfModel
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC                       :: InitDSMCSurfModel
PUBLIC                       :: FinalizeDSMCSurfModel
!===================================================================================================================================

CONTAINS


SUBROUTINE InitDSMCSurfModel()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals,             ONLY : abort
USE MOD_Mesh_Vars,           ONLY : nBCSides
USE MOD_DSMC_Vars,           ONLY : Adsorption, SurfMesh
USE MOD_PARTICLE_Vars,       ONLY : nSpecies, PDM
USE MOD_PARTICLE_Vars,       ONLY : KeepWallParticles, PEM
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                      !
  CHARACTER(32)                    :: hilf , hilf2  
  INTEGER                          :: iSpec, iSide, IDcounter
  REAL                             :: maxPart, SurfArea
#endif
!===================================================================================================================================
KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
IF (KeepWallParticles) THEN
  ALLOCATE(PDM%ParticleAtWall(1:PDM%maxParticleNumber)  , &
          PDM%PartAdsorbSideIndx(1:2,1:PDM%maxParticleNumber))
  PDM%ParticleAtWall(1:PDM%maxParticleNumber) = .FALSE.
  ALLOCATE(PEM%wNumber(1:nElems))
END IF
ALLOCATE(Adsorption%AdsorpInfo(1:nSpecies))
DO iSpec = 1,nSpecies
  ALLOCATE(Adsorption%AdsorpInfo(iSpec)%ProbAds(1:SurfMesh%nSurfaceBCSides),&
            Adsorption%AdsorpInfo(iSpec)%ProbDes(1:SurfMesh%nSurfaceBCSides))
#if (PP_TimeDiscMethod==42)
  ALLOCATE(Adsorption%AdsorpInfo(iSpec)%NumOfAds(1:SurfMesh%nSurfaceBCSides),&
            Adsorption%AdsorpInfo(iSpec)%NumOfDes(1:SurfMesh%nSurfaceBCSides))
#endif
END DO
DO iSpec = 1,nSpecies
  DO iSide = 1,SurfMesh%nSurfaceBCSides
    Adsorption%AdsorpInfo(iSpec)%ProbAds(iSide)=0
    Adsorption%AdsorpInfo(iSpec)%ProbDes(iSide)=0
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(iSpec)%NumOfAds(iSide)=0
    Adsorption%AdsorpInfo(iSpec)%NumOfDes(iSide)=0
    Adsorption%TPD = GETLOGICAL('Particles-DSMC-Adsorption-doTPD','.FALSE.')
    Adsorption%TPD_beta = GETREAL('Particles-DSMC-Adsorption-TPD-Beta','0.')
    Adsorption%TPD_Temp = 0.
#endif
  END DO
END DO

ALLOCATE(Adsorption%Coverage(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%SumDesorbPart(1:2,1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%SumAdsorbPart(1:2,1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%MaxCoverage(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%InitStick(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%PrefactorStick(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%Adsorbexp(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%Nu_a(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%Nu_b(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%DesorbEnergy(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%Intensification(1:SurfMesh%nSurfaceBCSides,1:nSpecies),&
          Adsorption%SurfSideToGlobSideMap(1:SurfMesh%nSurfaceBCSides),&
          Adsorption%DensSurfAtoms(1:SurfMesh%nSurfaceBCSides))
IDcounter = 0         
DO iSide=1, nBCSides 
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN
    IDcounter = IDcounter + 1
    Adsorption%SurfSideToGlobSideMap(IDcounter) = iSide
  END IF
END DO
DO iSurf = 1, SurfMesh%nSurfaceBCSides
  WRITE(UNIT=hilf,FMT='(I2)') iSurf
  Adsorption%DensSurfAtoms(iSurf) = GETREAL('Particles-Surface'//TRIM(hilf)//'-AtomsDensity','1.5E+19')
END DO

DO iSpec = 1, nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%Coverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialCoverage','0.')
  Adsorption%MaxCoverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
  Adsorption%InitStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
  Adsorption%PrefactorStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
  Adsorption%Adsorbexp(:,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
  Adsorption%Nu_a(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
  Adsorption%Nu_b(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
  Adsorption%DesorbEnergy(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
  Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
END DO
Adsorption%SumDesorbPart(:,:,:) = 0
Adsorption%SumAdsorbPart(:,:,:) = 0

END SUBROUTINE InitDSMCSurfModel

END MODULE MOD_DSMC_SurfModelInit