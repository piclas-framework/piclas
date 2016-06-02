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
USE MOD_Globals,                ONLY : abort
USE MOD_Mesh_Vars,              ONLY : nElems, nBCSides, BC
USE MOD_DSMC_Vars,              ONLY : Adsorption
USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PDM
USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, PEM
USE MOD_ReadInTools
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES                                                                      !
  CHARACTER(32)                    :: hilf , hilf2  
  INTEGER                          :: iSpec, iSide, iSurf, p, q, IDcounter
  REAL                             :: maxPart, SurfArea
!===================================================================================================================================
KeepWallParticles = GETLOGICAL('Particles-KeepWallParticles','.FALSE.')
IF (KeepWallParticles) THEN
  ALLOCATE(PDM%ParticleAtWall(1:PDM%maxParticleNumber)  , &
          PDM%PartAdsorbSideIndx(1:3,1:PDM%maxParticleNumber))
  PDM%ParticleAtWall(1:PDM%maxParticleNumber) = .FALSE.
  ALLOCATE(PEM%wNumber(1:nElems))
END IF
! allocate info and constants
#if (PP_TimeDiscMethod==42)
ALLOCATE( Adsorption%AdsorpInfo(1:nSpecies))
#endif
ALLOCATE( Adsorption%InitStick(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%PrefactorStick(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Adsorbexp(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Nu_a(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Nu_b(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%DesorbEnergy(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Intensification(1:SurfMesh%nSides,1:nSpecies))
ALLOCATE( Adsorption%HeatOfAdsZero(1:nSpecies),& 
!           Adsorption%DissResultsSpecNum(1:nSpecies),&
!           Adsorption%DissResultsSpec(1:2,1:(2*(nSpecies-1)-1),1:nSpecies),&
!           Adsorption%EDissBond(1:nSpecies,1:nSpecies,1:nSpecies),& 
          Adsorption%Coordination(1:nSpecies))
! initialize info and constants
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
  Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
  Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
#endif
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%InitStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialStick','0.')
  Adsorption%PrefactorStick(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-PrefactorStick','0.')
  Adsorption%Adsorbexp(:,iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Adsorbexp','1')
  Adsorption%Nu_a(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-a','0.')
  Adsorption%Nu_b(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Nu-b','0.')
  Adsorption%DesorbEnergy(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Desorption-Energy-K','1.')
  Adsorption%Intensification(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-Intensification-K','0.')
  Adsorption%Coordination(iSpec) = GETINT('Part-Species'//TRIM(hilf)//'-Coordination','0')
  Adsorption%HeatOfAdsZero(iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-HeatOfAdsorption-K','0.')
END DO
#if (PP_TimeDiscMethod==42)
  Adsorption%TPD = GETLOGICAL('Particles-DSMC-Adsorption-doTPD','.FALSE.')
  Adsorption%TPD_beta = GETREAL('Particles-DSMC-Adsorption-TPD-Beta','0.')
  Adsorption%TPD_Temp = 0.
#endif
! allocate and initialize adsorption variables
ALLOCATE( Adsorption%MaxCoverage(1:SurfMesh%nSides,1:nSpecies),& 
          Adsorption%Coverage(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbAds(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%ProbDes(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%SumDesorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SumAdsorbPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies),&
          Adsorption%SurfSideToGlobSideMap(1:SurfMesh%nSides),&
          Adsorption%Sigma(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%ProbSigma(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies,1:36*nSpecies),&
          Adsorption%DensSurfAtoms(1:SurfMesh%nSides))
IDcounter = 0         
DO iSide = 1,nBCSides 
  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN
    IDcounter = IDcounter + 1
    Adsorption%SurfSideToGlobSideMap(IDcounter) = iSide
  END IF
END DO
DO iSurf = 1,SurfMesh%nSides
  WRITE(UNIT=hilf,FMT='(I2)') iSurf
  Adsorption%DensSurfAtoms(iSurf) = GETREAL('Particles-Surface'//TRIM(hilf)//'-AtomsDensity','1.5E+19')
END DO
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I2)') iSpec
  Adsorption%MaxCoverage(:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-MaximumCoverage','0.')
  Adsorption%Coverage(:,:,:,iSpec) = GETREAL('Part-Species'//TRIM(hilf)//'-InitialCoverage','0.')
END DO
Adsorption%ProbAds(:,:,:,:) = 0.
Adsorption%ProbDes(:,:,:,:,:) = 0.
Adsorption%SumDesorbPart(:,:,:,:) = 0
Adsorption%SumAdsorbPart(:,:,:,:) = 0
Adsorption%Sigma(:,:,:,:,:) = 1.
Adsorption%ProbSigma(:,:,:,:,:) = 0.

CALL Init_SurfDist()

END SUBROUTINE InitDSMCSurfModel

SUBROUTINE Init_SurfDist()
!===================================================================================================================================
! Initializing surface distibution Model for calculating of coverage effects on heat of adsorption
!===================================================================================================================================
  USE MOD_Globals,                ONLY : abort
  USE MOD_ReadInTools
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species, KeepWallParticles
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
  INTEGER                          :: SurfSideID, subsurfxi, subsurfeta, iSpec
  INTEGER                          :: surfsquare, xi, eta, i, j, dist, Adsorbates
  INTEGER                          :: left, right, up, down, counter, first,second
  INTEGER,ALLOCATABLE              :: xSurfIndx1(:), ySurfIndx1(:), xSurfIndx2(:), ySurfIndx2(:), xSurfIndx3(:), ySurfIndx3(:)
  REAL                             :: RanNum!, RanNum2
  INTEGER                          :: leftxi,lefteta,rightxi,righteta,upxi,upeta,downxi,downeta
  INTEGER                          :: firstval,secondval,thirdval,fourthval,pos
  INTEGER                          :: Surfpos, Surfnum1, Surfnum2, Surfnum3, Indx, Indy
!===================================================================================================================================
ALLOCATE(SurfDistInfo(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides))
DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    ALLOCATE(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1:3))
  END DO
  END DO
END DO
IF (.NOT.KeepWallParticles) THEN
  surfsquare = GETINT('Particles-DSMC-AdsorptionSites','100')
END IF

DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    IF (KeepWallParticles) THEN ! does not work with vMPF
      surfsquare = INT(Adsorption%DensSurfAtoms(SurfSideID) &
                    * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID) &
                    / Species(1)%MacroParticleFactor)
    END IF
    surfsquare = INT(SQRT(REAL(surfsquare))) - 1
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1) = INT(surfsquare**2)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2) = INT(2*(surfsquare**2))
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3) = INT((surfsquare+1)**2)
    
    ALLOCATE( SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(1:nSpecies,1:surfsquare,1:surfsquare),&
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(1:nSpecies,1:(surfsquare*2),1:(surfsquare)),&
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(1:nSpecies,1:surfsquare+1,1:surfsquare+1),&
!               SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%TopSite(1:nSpecies,1:surfsquare,1:surfsquare),&
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(1:3))
              
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(1:nSpecies,1:surfsquare,1:surfsquare) = .FALSE.
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(1:nSpecies,1:(surfsquare*2),1:surfsquare) = .FALSE.
!     SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%TopSite(1:nSpecies,1:surfsquare+1,1:surfsquare+1) = .FALSE.
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(:,:,:) = 0

    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(:) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(:)
  END DO
  END DO
END DO


DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    DO iSpec = 1,nSpecies
      ! surfsqaure chosen for nSite(1) only for right xSurfInd definition (see IF cases below)
      surfsquare = INT(SQRT(REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1))))
      Adsorbates = INT(Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,iSpec) &
                  * SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec)))
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Adsorption%Coordination(iSpec)).LT.Adsorbates) THEN
        CALL abort(&
        __STAMP__&
        ,'Error in Init_SurfDist: Too many Adsorbates! - Choose lower Coverages')
      END IF
      
      ! allocate and define surface indexes for adsorbate distribution
      IF (iSpec.EQ.1) THEN
        SDEALLOCATE(xSurfIndx1)
        SDEALLOCATE(ySurfIndx1)
        SDEALLOCATE(xSurfIndx2)
        SDEALLOCATE(ySurfIndx2)
        SDEALLOCATE(xSurfIndx3)
        SDEALLOCATE(ySurfIndx3)
        ALLOCATE( xSurfIndx1(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1)),&
                  ySurfIndx1(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1)),&
                  xSurfIndx2(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2)),&
                  ySurfIndx2(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2)),&
                  xSurfIndx3(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3)),&
                  ySurfIndx3(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3)))
        Surfnum1 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(1)
        Surfnum2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(2)
        Surfnum3 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3)
        
        Indx = 1
        Indy = 1
        DO Surfpos = 1,Surfnum1
          IF (Indx.GT.surfsquare) THEN
            Indx = 1
            Indy = Indy + 1
          END IF
          xSurfIndx1(Surfpos) = Indx
          ySurfIndx1(Surfpos) = Indy
          Indx = Indx + 1
        END DO
        Indx = 1
        Indy = 1
        DO Surfpos = 1,Surfnum2
          IF (Indx.GT.(2*surfsquare)) THEN
            Indx = 1
            Indy = Indy + 1
          END IF
          xSurfIndx2(Surfpos) = Indx
          ySurfIndx2(Surfpos) = Indy
          Indx = Indx + 1
        END DO
        Indx = 1
        Indy = 1
        DO Surfpos = 1,Surfnum3
          IF (Indx.GT.surfsquare+1) THEN
            Indx = 1
            Indy = Indy + 1
          END IF
          xSurfIndx3(Surfpos) = Indx
          ySurfIndx3(Surfpos) = Indy
          Indx = Indx + 1
        END DO
      END IF
      
      ! distribute adsorbates randomly on the surface on the correct site
      dist = 1
      DO While (dist.LE.Adsorbates) 
        CALL RANDOM_NUMBER(RanNum)
        SELECT CASE(Adsorption%Coordination(iSpec))
        CASE(1)
          Surfpos = 1 + INT(Surfnum1 * RanNum)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(iSpec,xSurfIndx1(Surfpos),ySurfIndx1(Surfpos)) = .TRUE.
          xSurfIndx1(Surfpos) = xSurfIndx1(Surfnum1)
          ySurfIndx1(Surfpos) = ySurfIndx1(Surfnum1)
          Surfnum1 = Surfnum1 - 1
          dist = dist + 1
        CASE(2)
          Surfpos = 1 + INT(Surfnum2 * RanNum)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(iSpec,xSurfIndx2(Surfpos),ySurfIndx2(Surfpos)) = .TRUE.
          xSurfIndx2(Surfpos) = xSurfIndx2(Surfnum2)
          ySurfIndx2(Surfpos) = ySurfIndx2(Surfnum2)
          Surfnum2 = Surfnum2 - 1
          dist = dist + 1
        CASE(3)
          Surfpos = 1 + INT(Surfnum3 * RanNum)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%TopSite(iSpec,xSurfIndx1(Surfpos),ySurfIndx1(Surfpos)) = .TRUE.
          xSurfIndx3(Surfpos) = xSurfIndx3(Surfnum3)
          ySurfIndx3(Surfpos) = ySurfIndx3(Surfnum3)
          Surfnum3 = Surfnum3 - 1
          dist = dist + 1
        CASE DEFAULT
        END SELECT
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Adsorption%Coordination(iSpec)) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Adsorption%Coordination(iSpec)) - 1
      END DO
      
    END DO
  END DO
  END DO
END DO

SDEALLOCATE(xSurfIndx1)
SDEALLOCATE(ySurfIndx1)
SDEALLOCATE(xSurfIndx2)
SDEALLOCATE(ySurfIndx2)
SDEALLOCATE(xSurfIndx3)
SDEALLOCATE(ySurfIndx3)
    
! assign bond order of surface atoms in the surfacelattice
DO SurfSideID=1,SurfMesh%nSides
  DO subsurfeta = 1,nSurfSample
  DO subsurfxi = 1,nSurfSample
    DO iSpec = 1,nSpecies
! surfsqaure for nSite(3) exactly number of surface atoms
surfsquare = INT(SQRT(REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(3))))    
SELECT CASE(Adsorption%Coordination(iSpec))
CASE(1)
  DO xi = 1,surfsquare
    DO eta = 1,surfsquare
      IF (xi.EQ.1) THEN
        left = surfsquare-1
        right = xi
!           ELSE IF (xi.EQ.(surfsquare)) THEN
!             left = xi - 1
!             right = 1
      ELSE
        left = xi - 1
        right = xi
      END IF
      IF (eta.EQ.1) THEN
        up = surfsquare-1
        down = eta
!           ELSE IF (eta.EQ.(surfsquare)) THEN
!             up = eta - 1
!             down = 1 
      ELSE
        up = eta - 1
        down = eta
      END IF
      
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(iSpec,left,up) .AND. xi.NE.1 .AND. eta.NE.1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(iSpec,right,up) .AND. eta.NE.1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(iSpec,left,down) .AND. xi.NE.1) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%HollowSite(iSpec,right,down)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
    END DO
  END DO
CASE(2)
  DO xi = 1,surfsquare
    DO eta = 1,surfsquare
      IF (xi.EQ.1) THEN
        leftxi = surfsquare
      ELSE
        leftxi = xi - 1
      END IF
      lefteta = eta
      rightxi = xi
      righteta = eta
      IF (eta.EQ.1) THEN
        upeta = surfsquare
      ELSE
        upeta = eta - 1
      END IF
      upxi = xi + surfsquare
      downxi = xi + surfsquare
      downeta = eta
      IF (xi.EQ.surfsquare) THEN
        rightxi = rightxi - 1
        upxi = upxi - 1
        downxi = downxi - 1
      END IF
      IF (eta.EQ.surfsquare) THEN
        downeta = downeta - 1
        lefteta = lefteta - 1
        righteta = righteta - 1
      END IF
      
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(iSpec,leftxi,lefteta) &
          .AND. xi.NE.1 .AND. (eta.NE.surfsquare)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(iSpec,rightxi,righteta) &
          .AND. (xi.NE.surfsquare) .AND. (eta.NE.surfsquare)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(iSpec,upxi,upeta) &
          .AND. eta.NE.1 .AND. (xi.NE.surfsquare)) THEN 
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%BridgeSite(iSpec,downxi,downeta) &
          .AND. (eta.NE.surfsquare) .AND. (xi.NE.surfsquare)) THEN
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xi,eta) + 1
      END IF
      
    END DO
  END DO
CASE DEFAULT
END SELECT

    END DO
  END DO
  END DO
END DO
    
END SUBROUTINE Init_SurfDist

SUBROUTINE FinalizeDSMCSurfModel()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
USE MOD_PARTICLE_Vars,          ONLY : PDM, PEM
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(PDM%ParticleAtWall)
SDEALLOCATE(PDM%PartAdsorbSideIndx)
SDEALLOCATE(PEM%wNumber)
#if (PP_TimeDiscMethod==42)
SDEALLOCATE( Adsorption%AdsorpInfo)
#endif
SDEALLOCATE(Adsorption%InitStick)
SDEALLOCATE(Adsorption%PrefactorStick)
SDEALLOCATE(Adsorption%Adsorbexp)
SDEALLOCATE(Adsorption%Nu_a)
SDEALLOCATE(Adsorption%Nu_b)
SDEALLOCATE(Adsorption%DesorbEnergy)
SDEALLOCATE(Adsorption%Intensification)
SDEALLOCATE(Adsorption%HeatOfAdsZero)
SDEALLOCATE(Adsorption%DissResultsSpecNum)
SDEALLOCATE(Adsorption%DissResultsSpec)
SDEALLOCATE(Adsorption%EDissBond)
SDEALLOCATE(Adsorption%Coordination)

SDEALLOCATE(Adsorption%MaxCoverage)
SDEALLOCATE(Adsorption%Coverage)
SDEALLOCATE(Adsorption%ProbAds)
SDEALLOCATE(Adsorption%ProbDes)
SDEALLOCATE(Adsorption%SumDesorbPart)
SDEALLOCATE(Adsorption%SumAdsorbPart)
SDEALLOCATE(Adsorption%SurfSideToGlobSideMap)
SDEALLOCATE(Adsorption%Sigma)
SDEALLOCATE(Adsorption%ProbSigma)
SDEALLOCATE(Adsorption%DensSurfAtoms)

SDEALLOCATE(SurfDistInfo%HollowSite)
SDEALLOCATE(SurfDistInfo%BridgeSite)
SDEALLOCATE(SurfDistInfo%TopSite)
SDEALLOCATE(SurfDistInfo%SurfAtomBondOrder)
SDEALLOCATE(SurfDistInfo%SitesRemain)
SDEALLOCATE(SurfDistInfo%nSites)
SDEALLOCATE(SurfDistInfo)

END SUBROUTINE FinalizeDSMCSurfModel

END MODULE MOD_DSMC_SurfModelInit