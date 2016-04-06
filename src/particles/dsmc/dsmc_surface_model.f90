#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModelInit
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE WallAdsorption
  MODULE PROCEDURE WallAdsorption
END INTERFACE

INTERFACE WallDesorption
  MODULE PROCEDURE WallDesorption
END INTERFACE

INTERFACE CalcAdsorbProb
  MODULE PROCEDURE CalcAdsorbProb
END INTERFACE

INTERFACE CalcDesorbProb
  MODULE PROCEDURE CalcDesorbProb
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC                       :: DSMC_Update_Wall_Vars
! PUBLIC                       :: WallAdsorption
PUBLIC                       :: Particle_Wall_Adsorb
PUBLIC                       :: Particle_Wall_Desorb
PUBLIC                       :: Calc_PartNum_Wall_Desorb
PUBLIC                       :: CalcAdsorbProb
PUBLIC                       :: CalcDesorbProb
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_Update_Wall_Vars()   
!===================================================================================================================================
! Update and sample DSMC-values for adsorption, desorption and reactions on surfaces
!===================================================================================================================================
  USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PDM
  USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, PEM
  USE MOD_Mesh_Vars,              ONLY : SideData
  USE MOD_DSMC_Vars,              ONLY : DSMC, Adsorption, SurfMesh
  USE MOD_DSMC_Analyze,           ONLY : CalcWallSample
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
!===================================================================================================================================
! argument list declaration                                                                        !
! Local variable declaration                                                                       !  
  INTEGER                          :: iSpec, iSurfSide
  REAL                             :: maxPart, SurfArea
!===================================================================================================================================

  IF (DSMC%WallModel.GT.0) THEN
#if (PP_TimeDiscMethod==42)
    DO iSpec = 1,nSpecies
      Adsorption%AdsorpInfo(iSpec)%NumOfAds(:) = 0
      IF (KeepWallParticles) THEN
        Adsorption%AdsorpInfo(iSpec)%NumOfDes(:) = 0
      END IF
    END DO
#endif
    IF (.NOT.KeepWallParticles) THEN
      DO iSpec = 1,nSpecies
        DO iSurfSide = 1,SurfMesh%nSurfaceBCSides
          SurfArea = SideData(1,Adsorption%SurfSideToGlobSideMap(iSurfSide))%area &
                + SideData(2,Adsorption%SurfSideToGlobSideMap(iSurfSide))%area
          maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfArea
          Adsorption%Coverage(iSurfSide,iSpec) = Adsorption%Coverage(iSurfSide,iSpec) &
              - (Adsorption%SumDesorbPart(1,iSurfSide,iSpec)+Adsorption%SumDesorbPart(2,iSurfSide,iSpec)) &
              * Species(iSpec)%MacroParticleFactor / maxPart
        END DO
      END DO
      Adsorption%SumDesorbPart(:,:,:) = 0
      Adsorption%SumAdsorbPart(:,:,:) = 0
    END IF
    CALL CalcAdsorbProb()
    IF (KeepWallParticles) CALL CalcDesorbprob()
  END IF

END SUBROUTINE DSMC_Update_Wall_Vars()  
    
SUBROUTINE Particle_Wall_Adsorb(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,BCSideID,adsindex) 
!===================================================================================================================================
! Particle Adsorption after wall collision
!===================================================================================================================================
  USE MOD_DSMC_Analyze,           ONLY : CalcWallSample
  USE MOD_Particle_Vars
  USE MOD_Mesh_Vars,              ONLY : ElemToSide, SideData
  USE MOD_DSMC_Vars,              ONLY : SurfMesh, useDSMC
  USE MOD_DSMC_Vars,              ONLY : Adsorption
  USE MOD_DSMC_Vars,              ONLY : PartStateIntEn, SpecDSMC, DSMC
  USE MOD_DSMC_Vars,              ONLY : CollisMode
  USE MOD_TimeDisc_Vars,          ONLY : TEnd
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER,INTENT(INOUT)            :: adsindex                                                    !
  REAL,INTENT(INOUT)               :: PartTrajectory(1:3), lengthPartTrajectory, alpha
  REAL,INTENT(IN)                  :: xi, eta
  INTEGER,INTENT(IN)               :: PartID, GlobSideID
  INTEGER,INTENT(IN),OPTIONAL      :: BCSideID
! Local variable declaration            
  REAL                             :: RanNum, maxPart                                             !
  INTEGER                          :: SurfSide, locBCID                                           !
  REAL                             :: TransACC                                                    !
  REAL                             :: VibACC                                                      !
  REAL                             :: RotACC                                                      !
  INTEGER                          :: VibQuant,VibQuantWall,VibQuantNew                           !
  REAL                             :: VibQuantNewR                                                !
  REAL                             :: EtraOld, EtraNew                                            !
  REAL                             :: EtraWall, ErotWall, ErotNew, EvibNew                        !
  REAL                             :: VelXold, VelYold, VelZold, VeloReal                         !
  REAL, PARAMETER                  :: PI=3.14159265358979323846                                   !
  REAL                             :: IntersectionPos(1:3)                                        !
  REAL                             :: TransArray(1:6),IntArray(1:6)                               !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  ! additional states
  locBCID=PartBound%MapToPartBC(BC(SideID))
  ! get BC values
  WallVelo     = PartBound%WallVelo(1:3,locBCID)
  WallTemp     = PartBound%WallTemp(locBCID)
  TransACC     = PartBound%TransACC(locBCID)
  VibACC       = PartBound%VibACC(locBCID)
  RotACC       = PartBound%RotACC(locBCID)

  TransArray(:) = 0.0
  IntArray(:) = 0.0
  SurfSide = SurfMesh%SideIDToSurfID(GlobSideID)
  
  ! compute p and q
  ! correction of xi and eta, can only be applied if xi & eta are not used later!
  Xitild =MIN(MAX(-1.,XI ),0.99)
  Etatild=MIN(MAX(-1.,Eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1 
  
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption%AdsorpInfo(PartSpecies(PartID))%ProbAds(p,q,SurfSide).GE.RanNum) .AND. &
       (Adsorption%Coverage(p,q,SurfSide,PartSpecies(PartID)).LT.Adsorption%MaxCoverage(SurfSide,PartSpecies(PartID))) ) THEN  
    
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(p,q,SurfSide)
    Adsorption%Coverage(p,q,SurfSide,PartSpecies(PartID)) = Adsorption%Coverage(p,q,SurfSide,PartSpecies(PartID)) & 
                                                 + Species(PartSpecies(PartID))%MacroParticleFactor/maxPart
    adsindex = 1
!     sum over collisionenergy of adsorbed particles
    ! Sample Adsorbing atoms
!     Adsorption%SumAdsorbPart(TriNum,SurfSide,PartSpecies(PartID)) = Adsorption%SumAdsorbPart(TriNum,SurfSide,PartSpecies(PartID)) + 1
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(PartSpecies(PartID))%NumOfAds(SurfSide) = Adsorption%AdsorpInfo(PartSpecies(PartID))%NumOfAds(SurfSide) + 1
#endif
    ! allocate particle belonging adsorbing side index and Side-Triangle
    IF (KeepWallParticles) THEN
      PDM%PartAdsorbSideIndx(1,PartID) = GlobSideID
      PDM%PartAdsorbSideIndx(2,PartID) = p
      PDM%PartAdsorbSideIndx(3,PartID) = q
    END IF
  
    LastPartPos(PartID,1) = PartState(PartID,1)
    LastPartPos(PartID,2) = PartState(PartID,2)
    LastPartPos(PartID,3) = PartState(PartID,3)
    VelXold = PartState(PartID,4)
    VelYold = PartState(PartID,5)
    VelZold = PartState(PartID,6)
    PartState(PartID,1)   = IntersectionPos(1)
    PartState(PartID,2)   = IntersectionPos(2)
    PartState(PartID,3)   = IntersectionPos(3)
    PartState(PartID,4)   = WallVelo(1)
    PartState(PartID,5)   = WallVelo(2)
    PartState(PartID,6)   = WallVelo(3)
  
    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(PartSpecies(PartID))%MassIC * VeloReal**2
    EtraWall = 0.0
    EtraNew = 0.0
    
    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    TransArray(4) = PartState(PartID,4)-VelXold
    TransArray(5) = PartState(PartID,5)-VelYold
    TransArray(6) = PartState(PartID,6)-VelZold
  
    !---- Internal energy accommodation
    IF (CollisMode.GT.1) THEN
    IF (SpecDSMC(PartSpecies(PartID))%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
      ErotNew  = PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib) &
                  - DSMC%GammaQuant)
      CALL RANDOM_NUMBER(RanNum)
      VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
      DO WHILE (VibQuantWall.GE.SpecDSMC(PartSpecies(PartID))%MaxVibQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(PartSpecies(PartID))%CharaTVib)
      END DO
      VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
      VibQuantNew = INT(VibQuantNewR)
      CALL RANDOM_NUMBER(RanNum)
      IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
        EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
      ELSE
        EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(PartSpecies(PartID))%CharaTVib
      END IF
      IntArray(4) = VibQuant + DSMC%GammaQuant * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib
      IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(PartSpecies(PartID))%CharaTVib
      IntArray(6) = EvibNew
      PartStateIntEn(PartID,1) = EvibNew
    END IF
    END IF
    !End internal energy accomodation
    
!----  Sampling at walls
    IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
      CALL CalcWallSample(PartID,SurfSide,TriNum,Transarray,IntArray)
    END IF
    
  ELSE
    adsindex = 0
  END IF
  
END SUBROUTINE Particle_Wall_Adsorb

SUBROUTINE Calc_PartNum_Wall_Desorb()
!===================================================================================================================================
! Particle Desorption when particles deleted at adsorption and inserted at desorption
! (use particle tracking at wall instead, if low particle densities, surface coverages and high wall temperatures)
!===================================================================================================================================
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,      ONLY : SideData
USE MOD_DSMC_Vars,      ONLY : Adsorption, SurfMesh
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER                          :: i, PartAds2                                                 !
! Local variable declaration                                                                       !
   INTEGER                          :: iSurfSide, iSpec                                            !
   INTEGER                          :: TriNum, NPois                                               !
   REAL                             :: PartAds, PartDes, RanNum, maxPart, SurfArea, TPois          !
! argument list declaration                                                                        !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
#if (PP_TimeDiscMethod==42) 
  DO iSpec = 1,nSpecies
    Adsorption%AdsorpInfo(iSpec)%NumOfDes(:) = 0
  END DO
#endif
  DO iSpec = 1,nSpecies
    DO iSurfSide = 1,SurfMesh%nSurfaceBCSides
      DO TriNum = 1,2
        IF (Adsorption%Coverage(iSurfSide,iSpec).GT.0) THEN
          PartAds = Adsorption%Coverage(iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                    * SideData(TriNum,Adsorption%SurfSideToGlobSideMap(iSurfSide))%area &
                    / Species(iSpec)%MacroParticleFactor
          PartDes = PartAds * Adsorption%AdsorpInfo(iSpec)%ProbDes(iSurfSide)
!           PartAds2 = INT(Adsorption%Coverage(iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
!                     * SideData(TriNum,Adsorption%SurfSideToGlobSideMap(iSurfSide))%area &
!                     / Species(iSpec)%MacroParticleFactor) 
!           DO i = 1,PartAds2
          CALL RANDOM_NUMBER(RanNum)
!           IF (Adsorption%AdsorpInfo(iSpec)%ProbDes(iSurfSide).GT.RanNum) THEN
!           Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) = Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) + 1
!           END IF
!           END DO
!           Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) = INT(PartAds &
!                                                            * Adsorption%AdsorpInfo(iSpec)%ProbDes(iSurfSide) + RanNum)

          
        IF (EXP(-PartDes).LE.TINY(PartDes) &
#if (PP_TimeDiscMethod==42)
        .OR. Adsorption%TPD &
#endif
        ) THEN
          Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) = INT(PartDes + RanNum)
        ELSE !poisson-sampling instead of random rounding (reduces numerical non-equlibrium effects [Tysanner and Garcia 2004]
          Npois=0
          Tpois=1.0
          DO
            Tpois=RanNum*Tpois
            IF (Tpois.LT.TINY(Tpois)) THEN
              Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) = INT(PartDes + RanNum)
              EXIT
            END IF
            IF (Tpois.GT.EXP(-PartDes)) THEN
              Npois=Npois+1
              CALL RANDOM_NUMBER(RanNum)
            ELSE
              Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec) = Npois
              EXIT
            END IF
          END DO
        END IF
#if (PP_TimeDiscMethod==42)
          Adsorption%AdsorpInfo(iSpec)%NumOfDes(iSurfSide) = Adsorption%AdsorpInfo(iSpec)%NumOfDes(iSurfSide) &
                                                           + Adsorption%SumDesorbPart(TriNum,iSurfSide,iSpec)
#endif
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE Calc_PartNum_Wall_Desorb

SUBROUTINE Particle_Wall_Desorb(i)
!===================================================================================================================================
! Particle Desorption when particles are remembered (use only if low particle densities, surface coverages 
! and high wall temperatures)
!===================================================================================================================================
USE MOD_DSMC_Analyze,   ONLY : CalcWallSample
USE MOD_Particle_Vars
USE MOD_Mesh_Vars,      ONLY : BC, SideToElem, ElemToSide
USE MOD_DSMC_Vars,      ONLY : Adsorption, SurfMesh, useDSMC
USE MOD_DSMC_Vars,      ONLY : PartStateIntEn, SpecDSMC, DSMC
USE MOD_DSMC_Vars,      ONLY : CollisMode, SampWall
USE MOD_TimeDisc_Vars,  ONLY : dt, TEnd
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   INTEGER, INTENT(IN)              :: i                                                           !
! Local variable declaration                                                                       !
   INTEGER                          :: iLocSide, SurfSide, globSide                                !
   INTEGER                          :: Element                                                     !
   INTEGER                          :: TriNum 
   REAL                             :: maxPart                                                     !
   REAL                             :: TransACC                                                    !
   REAL                             :: VibACC                                                      !
   REAL                             :: RotACC                                                      !
   REAL                             :: WallTemp                                                    !
   REAL                             :: WallVelo(1:3)                                               !
   INTEGER                          :: Node1, Node2                                                !
   INTEGER                          :: VibQuant,VibQuantWall,VibQuantNew                           !
   REAL                             :: VibQuantNewR                                                !
   REAL                             :: nx, ny, nz, nVal                                            !
   REAL                             :: xNod, yNod, zNod, VeloReal, VeloCrad, VeloCx, VeloCy ,VeloCz!
   REAL                             :: EtraOld, EtraNew, RanNum, Cmr, Phi, Fak_D                   !
   REAL                             :: EtraWall, ErotWall, EvibNew, ErotNew                        !
   REAL                             :: VelX, VelY, VelZ, VecX, VecY, VecZ                          !
   REAL                             :: VelXold, VelYold, VelZold
   REAL                             :: Vector1(1:3), Vector2(1:3)                                  !
   REAL, PARAMETER                  :: PI=3.14159265358979323846                                   !
   REAL                             :: IntersectionPos(1:3)                                        !
   REAL                             :: TransArray(1:6),IntArray(1:6)                               !
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
  TransArray(:) = 0.0
  IntArray(:) = 0.0
  globSide = PDM%PartAdsorbSideIndx(1,i)
  SurfSide = SurfMesh%GlobSideToSurfSideMap(globSide)
    
  TransACC = 1
  VibACC = 1
  RotACC = 1
  
  CALL RANDOM_NUMBER(RanNum)
  IF (Adsorption%AdsorpInfo(PartSpecies(i))%ProbDes(SurfSide).GE.RanNum) THEN 
 
    PDM%ParticleAtWall(i)=.FALSE.
    ! Sample desorping particles
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(PartSpecies(i))%NumOfDes(SurfSide) = Adsorption%AdsorpInfo(PartSpecies(i))%NumOfDes(SurfSide) + 1
#endif
    
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(SurfSide) !* Adsorption%MaxCoverage(:,iSpec)
    Adsorption%Coverage(SurfSide,PartSpecies(i)) = Adsorption%Coverage(SurfSide,PartSpecies(i)) & 
                                                  - Species(PartSpecies(i))%MacroParticleFactor/maxPart
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
    WallVelo(1:3) = PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(globSide)))
    
    Element = SideToElem(1,globSide)
    IF (Element.LT.1) THEN
      Element = SideToElem(2,globSide)
      iLocSide = SideToElem(4,globSide)
    ELSE
      iLocSide = SideToElem(3,globSide)
    END IF 
      
    xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
    yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
    zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

    !---- Calculate normal vector in global coordinates : (points into Element)
    TriNum = PDM%PartAdsorbSideIndx(2,i)
    
    Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
    Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle

    Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
    Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
    Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod

    Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
    Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
    Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod

    nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
    ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
    nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

    nVal = SQRT(nx*nx + ny*ny + nz*nz)

    nx = nx/nVal
    ny = ny/nVal
    nz = nz/nVal
    
    !---- Calculate new velocity vector (Extended Maxwellian Model)
  !    VeloReal = SQRT(PartState(i,4) * PartState(i,4) + &
  !              PartState(i,5) * PartState(i,5) + &
  !              PartState(i,6) * PartState(i,6))
  !    EtraOld = 0.5 * Species(PartSpecies(i))%MassIC * VeloReal**2
    VeloReal = 0.0
    EtraOld = 0.0
    CALL RANDOM_NUMBER(RanNum)
    VeloCrad    = SQRT(-LOG(RanNum))
    CALL RANDOM_NUMBER(RanNum)
    VeloCz      = SQRT(-LOG(RanNum))
    Fak_D       = VeloCrad**2 + VeloCz**2
    EtraWall    = BoltzmannConst * WallTemp * Fak_D
  !    EtraNew = EtraOld + TransACC * (EtraWall - EtraOld)
    EtraNew = EtraWall
    Cmr     = SQRT(2.0 * EtraNew / (Species(PartSpecies(i))%MassIC * Fak_D))
    CALL RANDOM_NUMBER(RanNum)
    Phi     = 2 * PI * RanNum
    VeloCx  = Cmr * VeloCrad * COS(Phi)
    VeloCy  = Cmr * VeloCrad * SIN(Phi)
    VeloCz  = Cmr * VeloCz
    
    !---- Transformation local distribution -> global coordinates
    VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
    VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
    VecZ = Vector1(3) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
    
    VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
    VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
    VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz

    lastPartPos(i,1) = PartState(i,1)
    lastPartPos(i,2) = PartState(i,2)
    lastPartPos(i,3) = PartState(i,3)
    VelXold = PartState(i,4)
    VelYold = PartState(i,5)
    VelZold = PartState(i,6)
    
    CALL RANDOM_NUMBER(RanNum)
    PartState(i,1)   = Partstate(i,1) + RanNum * dt * VelX
    PartState(i,2)   = Partstate(i,2) + RanNum * dt * VelY
    PartState(i,3)   = Partstate(i,3) + RanNum * dt * VelZ
    PartState(i,4) = PartState(i,4) + VelX
    PartState(i,5) = PartState(i,5) + VelY
    PartState(i,6) = PartState(i,6) + VelZ
   
    TransArray(1) = EtraOld
    TransArray(2) = EtraWall
    TransArray(3) = EtraNew
    TransArray(4) = PartState(i,4)-VelXold
    TransArray(5) = PartState(i,5)-VelYold
    TransArray(6) = PartState(i,6)-VelZold
    
    !----  Sampling at walls
    IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
      CALL CalcWallSample(i,SurfSide,TriNum,Transarray,IntArray)
    END IF

 ELSE 
   LastPartPos(i,1) = PartState(i,1)
   LastPartPos(i,2) = PartState(i,2)
   LastPartPos(i,3) = PartState(i,3)
 END IF !End Desorption

END SUBROUTINE Particle_Wall_Desorb

SUBROUTINE CalcAdsorbProb()
! Adsorption Models
  USE MOD_Particle_Vars
  USE MOD_DSMC_Vars,              ONLY : SurfMesh, Adsorption
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! Local variable declaration                                                                       !
   INTEGER                          :: BC_Side, iSpec
   REAL                             :: Theta_req, Kfactor, S_0
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------! 

!----------------------------------------------------------------------------------------------------------------------------------- 
! Kisluik Sticking Model from Kolasinski's Surface Science (book)
  DO BC_Side=1,SurfMesh%nSurfaceBCSides
    DO iSpec=1,nSpecies
      ! enhance later to co-adsorption
      Theta_req = (1.0 - Adsorption%Coverage(BC_Side,iSpec)/Adsorption%MaxCoverage(BC_Side,iSpec)) &
                **Adsorption%Adsorbexp(BC_Side,iSpec)  
      !----- kann später auf von Wandtemperatur abhängige Werte erweitert werden          
      Kfactor = Adsorption%PrefactorStick(BC_Side,iSpec)
      S_0 = Adsorption%InitStick(BC_Side,iSpec)
      !-----
      IF (Theta_req.EQ.0) THEN
        Adsorption%AdsorpInfo(iSpec)%ProbAds(BC_Side) = 0.
      ELSE
        Adsorption%AdsorpInfo(iSpec)%ProbAds(BC_Side) = S_0 / (1.0 + Kfactor * ( 1.0/Theta_req - 1.0))
      END IF
    END DO
  END DO
!----------------------------------------------------------------------------------------------------------------------------------- 
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP) model
!   DO BC_Side=1,SurfMesh%nSurfaceBCSides
!     DO iSpec=1,nSpecies
!       CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct)
!       CALL PartitionFuncGas(iSpec, WallTemp, VarPartitionFuncGas)
!     END DO
!   END DO

END SUBROUTINE CalcAdsorbProb

SUBROUTINE CalcDesorbProb()
! Desorption Models
  USE MOD_Particle_Vars
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : SurfMesh, Adsorption
  USE MOD_TimeDisc_Vars,          ONLY : dt
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
#endif
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! Local variable declaration                                                                       !
   INTEGER                          :: BC_Side, iSpec, globSide
   REAL                             :: Theta, E_des, nu_des, rate, WallTemp                        !
!
   REAL                             :: Q_0A, D_A, sigma, Q, m
   INTEGER                          :: n
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------! 
  DO BC_Side=1,SurfMesh%nSurfaceBCSides
    globSide = Adsorption%SurfSideToGlobSideMap(BC_Side)
! special TPD (temperature programmed desorption) temperature adjustment routine    
#if (PP_TimeDiscMethod==42)
    IF (Adsorption%TPD) THEN
      WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide))) + (Adsorption%TPD_beta * iter * dt)
      Adsorption%TPD_Temp = Walltemp
    ELSE
#else
      WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
#endif
#if (PP_TimeDiscMethod==42)
    END IF
#endif
!----------------------------------------------------------------------------------------------------------------------------------- 
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP)
DO iSpec=1,nSpecies   
  Theta = Adsorption%Coverage(BC_Side,iSpec)  
  IF (Theta.GT.0) THEN
    Q_0A = 26312.
    D_A  = 59870.
    n = 2
    m = Theta * n
    IF (m.LT.1) THEN
!       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A)
      Q = (Q_0A**2)/(Q_0A/n+D_A)
    ELSE 
      sigma = 1/m*(2-1/m)
!       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A) *sigma
      Q = (Q_0A**2)/(Q_0A/n+D_A) * sigma
    END IF
    
    CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct)
    CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
    
    E_des = Q
    nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
    rate = nu_des * exp(-E_des/WallTemp)
    Adsorption%AdsorpInfo(iSpec)%ProbDes(BC_Side) = rate * dt
  ELSE
    Adsorption%AdsorpInfo(iSpec)%ProbDes(BC_Side) = 0.0
  END IF
      
!     write(*,*)Q,E_des,nu_des,rate,Theta
END DO
!----------------------------------------------------------------------------------------------------------------------------------- 
! absolute rate theory / lattice gas model (QCA-quasi-chemical-approximation)

!-----------------------------------------------------------------------------------------------------------------------------------  
! Polanyi-Wigner-eq. from Kolasinski's Surface Science (book)
!     DO iSpec=1,nSpecies   
!       Theta = Adsorption%Coverage(BC_Side,iSpec)! / Adsorption%MaxCoverage(BC_Side,iSpec)
!       !----- kann später auf von Wandtemperatur/Translationsenergie abhängige Werte erweitert werden          
!       E_des = Adsorption%DesorbEnergy(BC_Side,iSpec) + Adsorption%Intensification(BC_Side,iSpec) * Theta
!       nu_des = 10**(Adsorption%Nu_a(BC_Side,iSpec) + Adsorption%Nu_b(BC_Side,iSpec) * Theta)/10000
!       !-----
!       rate = nu_des *(Adsorption%DensSurfAtoms(BC_Side)**(Adsorption%Adsorbexp(BC_Side,iSpec)-1)) &
!                     * (Theta**Adsorption%Adsorbexp(BC_Side,iSpec)) * exp(-E_des/WallTemp)
!       IF (Theta.GT.0) THEN
!         Adsorption%AdsorpInfo(iSpec)%ProbDes(BC_Side) = rate * dt / Theta
!       ELSE
!         Adsorption%AdsorpInfo(iSpec)%ProbDes(BC_Side) = 0.0
!       END IF
!     END DO

  END DO

END SUBROUTINE CalcDesorbProb

! SUBROUTINE PartitionFuncGas(iSpec, Temp, VarPartitionFuncGas)
! !===================================================================================================================================
! ! Init of DSMC Vars
! !===================================================================================================================================
! ! MODULES
! USE MOD_Globals
! USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
! USE MOD_Particle_Vars,      ONLY: BoltzmannConst, PlanckConst, Species
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! INTEGER, INTENT(IN)           :: iSpec, Temp
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! !-----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! REAL, PARAMETER               :: Pi=3.14159265358979323846_8
! INTEGER                       :: iPolyatMole, iDOF
! REAL                          :: Qtra, Qrot, Qvib, Qelec, Temp
! !===================================================================================================================================
!   Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
!   IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
!     IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
!       iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!       IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
!         Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
!       ELSE
!         Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
!       END IF
!       Qvib = 1.
!       DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!         Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
!                 / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
!       END DO
!     ELSE
!       Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
!       Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
!     END IF
!   ELSE
!     Qrot = 1.
!     Qvib = 1.
!   END IF
!   Qelec = 0.
!   DO iDOF=1, SpecDSMC(iSpec)%NumElecLevels
!     Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
!   END DO
!   SpecDSMC(iSpec)%PartitionFunction(iInter) = Qtra * Qrot * Qvib * Qelec
! 
! END SUBROUTINE PartitionFuncGas

SUBROUTINE PartitionFuncAct(iSpec, Temp, VarPartitionFuncAct)
!===================================================================================================================================
! Partitionfunction of activated complex
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY: BoltzmannConst, PlanckConst, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
  Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))*1.5E-19
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
      IF(PolyatomMolDSMC(iPolyatMole)%LinearMolec) THEN
!         Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1))
        Qrot = Temp / (2 * 2.1)
      ELSE
!         Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(1)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(2)    &
!                                                                         * PolyatomMolDSMC(iPolyatMole)%CharaTRotDOF(3)))
        Qrot = SQRT(Pi) / SpecDSMC(iSpec)%SymmetryFactor * SQRT(Temp**3/( 2.1**3))
      END IF
      Qvib = 1.
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
      END DO
    ELSE
!       Qrot = Temp / (SpecDSMC(iSpec)%SymmetryFactor * SpecDSMC(iSpec)%CharaTRot)
      Qrot = Temp / (2 * 2.1)
      Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
    END IF
  ELSE
    Qrot = 1.
    Qvib = 1.
  END IF
  VarPartitionFuncAct = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncAct

SUBROUTINE PartitionFuncSurf(iSpec, Temp, VarPartitionFuncSurf)
!===================================================================================================================================
! partition function of adsorbates
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY: BoltzmannConst, PlanckConst, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncSurf
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
  Qtra = 1
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
        Qvib = Qvib * EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / (2. * Temp)) &
                / (1. - EXP(-PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) / Temp))
      END DO
    ELSE
      Qvib = EXP(-SpecDSMC(iSpec)%CharaTVib / (2. * Temp)) / (1. - EXP(-SpecDSMC(iSpec)%CharaTVib / Temp))
    END IF
  ELSE
    Qvib = 1.
  END IF
  Qrot = 1.
  VarPartitionFuncSurf = Qtra * Qrot * Qvib

END SUBROUTINE PartitionFuncSurf