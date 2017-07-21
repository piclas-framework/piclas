#include "boltzplatz.h"

MODULE MOD_DSMC_SurfModel_Tools
!===================================================================================================================================
! Initialization of DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC                       :: DSMC_Update_Wall_Vars
PUBLIC                       :: Particle_Wall_Adsorb
! PUBLIC                       :: Particle_Wall_Desorb
PUBLIC                       :: Calc_PartNum_Wall_Desorb
PUBLIC                       :: CalcAdsorbProb
PUBLIC                       :: CalcDesorbProb
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_Update_Wall_Vars()   
!===================================================================================================================================
! Update and sample DSMC-values for adsorption, desorption and reactions on surfaces
!===================================================================================================================================
  USE MOD_PARTICLE_Vars,          ONLY : nSpecies
  USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, Species
  USE MOD_DSMC_Vars,              ONLY : DSMC, Adsorption
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration 
  INTEGER                          :: iSpec, iSurfSide, p, q
  REAL                             :: maxPart
!===================================================================================================================================

  IF (DSMC%WallModel.GT.0) THEN
#if (PP_TimeDiscMethod==42)
    DO iSpec = 1,nSpecies
      Adsorption%AdsorpInfo(iSpec)%WallCollCount = 0
      Adsorption%AdsorpInfo(iSpec)%NumOfAds = 0
      Adsorption%AdsorpInfo(iSpec)%MeanProbAds = 0.
      IF (KeepWallParticles) THEN
        Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
      END IF
    END DO
#endif
    IF (.NOT.KeepWallParticles) THEN
      DO iSpec = 1,nSpecies
        DO iSurfSide = 1,SurfMesh%nSides
          DO q = 1,nSurfSample
            DO p = 1,nSurfSample
              maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
              Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  - Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)* Species(iSpec)%MacroParticleFactor / maxPart
            END DO
          END DO
        END DO
      END DO
      Adsorption%SumDesorbPart(:,:,:,:) = 0
      Adsorption%SumAdsorbPart(:,:,:,:) = 0
    END IF
    CALL CalcAdsorbProb()
    IF (KeepWallParticles) CALL CalcDesorbprob()
  END IF

END SUBROUTINE DSMC_Update_Wall_Vars  
    
SUBROUTINE Particle_Wall_Adsorb(PartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,adsindex) 
!===================================================================================================================================
! Particle Adsorption after wall collision
!===================================================================================================================================
  USE MOD_DSMC_Analyze,           ONLY : CalcWallSample
  USE MOD_Particle_Vars,          ONLY : WriteMacroValues, KeepWallParticles, PDM
  USE MOD_Particle_Vars,          ONLY : PartState,LastPartPos,Species,BoltzmannConst,PartSpecies
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : CollisMode, Adsorption
  USE MOD_DSMC_Vars,              ONLY : PartStateIntEn, SpecDSMC, DSMC, PolyatomMolDSMC, VibQuantsPar
  USE MOD_Particle_Boundary_Vars, ONLY : SurfMesh, dXiEQ_SurfSample, Partbound
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! argument list declaration
  INTEGER,INTENT(INOUT)            :: adsindex
  REAL,INTENT(INOUT)               :: PartTrajectory(1:3), alpha
  REAL,INTENT(IN)                  :: xi, eta
  INTEGER,INTENT(IN)               :: PartID, GlobSideID
  LOGICAL,INTENT(IN)               :: IsSpeciesSwap
! LOCAL VARIABLES
  INTEGER                              :: locBCID, VibQuant, VibQuantNew, VibQuantWall
  REAL                                 :: VibQuantNewR
  REAL                                 :: VeloReal, RanNum, EtraOld
  REAL                                 :: EtraWall, EtraNew
  REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
  REAL                                 :: ErotNew, ErotWall, EVibNew
  REAL                                 :: Xitild,EtaTild
  INTEGER                              :: p,q
! Polyatomic Molecules
  REAL, ALLOCATABLE                    :: RanNumPoly(:), VibQuantNewRPoly(:)
  INTEGER                              :: iPolyatMole, iDOF
  INTEGER, ALLOCATABLE                 :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
! Local variable declaration            
  REAL                             :: maxPart
  INTEGER                          :: SurfSide, iSpec
  REAL                             :: VelXold, VelYold, VelZold
  REAL, PARAMETER                  :: PI=3.14159265358979323846
  REAL                             :: IntersectionPos(1:3)
  REAL                             :: TransArray(1:6),IntArray(1:6)
!===================================================================================================================================
  ! additional states
  locBCID=PartBound%MapToPartBC(BC(GlobSideID))
  ! get BC values
  WallVelo     = PartBound%WallVelo(1:3,locBCID)
  WallTemp     = PartBound%WallTemp(locBCID)
  TransACC     = PartBound%TransACC(locBCID)
  VibACC       = PartBound%VibACC(locBCID)
  RotACC       = PartBound%RotACC(locBCID)

  TransArray(:) = 0.0
  IntArray(:) = 0.0
  SurfSide = SurfMesh%SideIDToSurfID(GlobSideID)
  iSpec = PartSpecies(PartID)
#if (PP_TimeDiscMethod==42)  
  ! Update wallcollision counter
  Adsorption%AdsorpInfo(iSpec)%WallCollCount = Adsorption%AdsorpInfo(iSpec)%WallCollCount + 1
#endif  
  ! compute p and q
  ! correction of xi and eta, can only be applied if xi & eta are not used later!
  Xitild =MIN(MAX(-1.,xi ),0.99)
  Etatild=MIN(MAX(-1.,eta),0.99)
  p=INT((Xitild +1.0)/dXiEQ_SurfSample)+1
  q=INT((Etatild+1.0)/dXiEQ_SurfSample)+1 
  
  CALL RANDOM_NUMBER(RanNum)
  IF ( (Adsorption%ProbAds(p,q,SurfSide,iSpec).GE.RanNum) .AND. &
       (Adsorption%Coverage(p,q,SurfSide,iSpec).LT.Adsorption%MaxCoverage(SurfSide,iSpec)) ) THEN  
    
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(p,q,SurfSide)
    Adsorption%Coverage(p,q,SurfSide,iSpec) = Adsorption%Coverage(p,q,SurfSide,iSpec) & 
                                                 + Species(iSpec)%MacroParticleFactor/maxPart
    adsindex = 1
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(iSpec)%NumOfAds = Adsorption%AdsorpInfo(iSpec)%NumOfAds + 1
#endif
    ! allocate particle belonging adsorbing side-index and side-subsurface-indexes
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
    ! intersection point with surface
    IntersectionPos(1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
    PartState(PartID,1)  = IntersectionPos(1)
    PartState(PartID,2)  = IntersectionPos(2)
    PartState(PartID,3)  = IntersectionPos(3)
    PartState(PartID,4)  = WallVelo(1)
    PartState(PartID,5)  = WallVelo(2)
    PartState(PartID,6)  = WallVelo(3)
  
    VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
    EtraOld = 0.5 * Species(iSpec)%MassIC * VeloReal**2
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
    IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
      !---- Rotational energy accommodation
      CALL RANDOM_NUMBER(RanNum)
      ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
      ErotNew  = PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
      IntArray(1) = PartStateIntEn(PartID,2)
      IntArray(2) = ErotWall
      IntArray(3) = ErotNew
      PartStateIntEn(PartID,2) = ErotNew
      !---- Vibrational energy accommodation
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        EvibNew = 0.0
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
                 VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
                 VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
        CALL RANDOM_NUMBER(RanNumPoly)
        VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
          CALL RANDOM_NUMBER(RanNumPoly)
          VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
        END DO
        VibQuantNewRPoly(:) = VibQuantsPar(PartID)%Quants(:) + VibACC*(VibQuantWallPoly(:) - VibQuantsPar(PartID)%Quants(:))
        VibQuantNewPoly = INT(VibQuantNewRPoly)
        DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
          CALL RANDOM_NUMBER(RanNum)
          IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
            EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
                      * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
          ELSE
            EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
                      * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
            VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
          END IF
          IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
                      * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(iSpec)%MacroParticleFactor
          IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
                      * SpecDSMC(iSpec)%CharaTVib * Species(iSpec)%MacroParticleFactor
        END DO
      ELSE
        VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) &
                    - DSMC%GammaQuant)
        CALL RANDOM_NUMBER(RanNum)
        VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
        DO WHILE (VibQuantWall.GE.SpecDSMC(iSpec)%MaxVibQuant)
          CALL RANDOM_NUMBER(RanNum)
          VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
        END DO
        VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
        VibQuantNew = INT(VibQuantNewR)
        CALL RANDOM_NUMBER(RanNum)
        IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
          EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
        ELSE
          EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
        END IF
        IntArray(4) = VibQuant + DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
        IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
      END IF
      IntArray(6) = EvibNew
    END IF
    END IF
    !End internal energy accomodation
    
!----  Sampling at walls
    IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
      CALL CalcWallSample(PartID,SurfSide,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap,locBCID)
    END IF
    
  ELSE
    adsindex = 0
  END IF
  
END SUBROUTINE Particle_Wall_Adsorb

SUBROUTINE Calc_PartNum_Wall_Desorb()
!===================================================================================================================================
! calculation of desorbing desorbing particle number when particles deleted at adsorption and inserted at desorption
!===================================================================================================================================
USE MOD_Particle_Vars,          ONLY : nSpecies, Species
USE MOD_DSMC_Vars,              ONLY : Adsorption
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
!    INTEGER                          :: i, PartAds2
   INTEGER                          :: iSurfSide, iSpec, p, q, NPois, WallPartNum
   REAL                             :: PartAds, PartDes, RanNum, Tpois
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(:)%NumOfDes = 0
#endif
  CALL CalcDesorbProb()
  
  DO iSpec = 1,nSpecies
    DO iSurfSide = 1,SurfMesh%nSides
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          IF (Adsorption%Coverage(p,q,iSurfSide,iSpec).GT.0) THEN
            WallPartNum = INT(Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                      * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                      / Species(iSpec)%MacroParticleFactor)
            IF (WallPartNum .GT. 0) THEN
!             IF (DSMC%WallModel.GT.1) THEN
!               PartAds2 = INT(Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
!                       * SurfMesh%SurfaceArea(p,q,iSurfSide) &
!                       / (Species(iSpec)%MacroParticleFactor))
!               Dist_AdsorbateNum = INT (PartAds2 / Dist_MPF)        
!               DO i = 1,Dist_AdsorbateNum
!                 CALL RANDOM_NUMBER(RanNum)
!                 IF (Adsorption%ProbDes(p,q,iSurfSide,iSpec).GT.RanNum) THEN
!                   Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) + 1
!                 END IF
!               END DO
!               Adsorption%SumDesorbPart(:,:,:,:) = Adsorption%SumDesorbPart(:,:,:,:) * Dist_MPF
!             ELSE
              PartAds = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                        * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                        / Species(iSpec)%MacroParticleFactor
              PartDes = PartAds * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
              IF (PartDes.GT.WallPartNum) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = WallPartNum
              ELSE
                CALL RANDOM_NUMBER(RanNum)            
                IF (EXP(-PartDes).LE.TINY(PartDes) &
#if (PP_TimeDiscMethod==42)
                .OR. Adsorption%TPD &
#endif
                ) THEN
                  Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = INT(PartDes + RanNum)
                ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
                  Npois=0
                  Tpois=1.0
                  DO
                    Tpois=RanNum*Tpois
                    IF (Tpois.LT.TINY(Tpois)) THEN
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = INT(PartDes + RanNum)
                      EXIT
                    END IF
                    IF (Tpois.GT.EXP(-PartDes)) THEN
                      Npois=Npois+1
                      CALL RANDOM_NUMBER(RanNum)
                    ELSE
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Npois
                      EXIT
                    END IF
                  END DO
                END IF
              END IF !PartDes.GT.WallPartNum
!             END IF  
            ELSE !not PartAds.GT.0
              Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = 0
            END IF !PartAds.GT.0
#if (PP_TimeDiscMethod==42)
            Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                                  + Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)
#endif
          END IF
        END DO  
      END DO
    END DO
  END DO
END SUBROUTINE Calc_PartNum_Wall_Desorb

! SUBROUTINE Particle_Wall_Desorb(i)
! !===================================================================================================================================
! ! Particle Desorption when particles are remembered (use only if low particle densities, surface coverages 
! ! and high wall temperatures)
! !===================================================================================================================================
! USE MOD_DSMC_Analyze,   ONLY : CalcWallSample
! USE MOD_Particle_Vars
! USE MOD_Mesh_Vars,      ONLY : BC, SideToElem, ElemToSide
! USE MOD_DSMC_Vars,      ONLY : Adsorption, SurfMesh, useDSMC
! USE MOD_DSMC_Vars,      ONLY : PartStateIntEn, SpecDSMC, DSMC
! USE MOD_DSMC_Vars,      ONLY : CollisMode, SampWall
! USE MOD_TimeDisc_Vars,  ONLY : dt, TEnd
! !===================================================================================================================================
!   IMPLICIT NONE
! !===================================================================================================================================
! ! argument list declaration
!    INTEGER, INTENT(IN)              :: i
! ! Local variable declaration
!    INTEGER                          :: iLocSide, SurfSide, globSide
!    INTEGER                          :: Element
!    INTEGER                          :: TriNum
!    REAL                             :: maxPart
!    REAL                             :: TransACC
!    REAL                             :: VibACC
!    REAL                             :: RotACC
!    REAL                             :: WallTemp
!    REAL                             :: WallVelo(1:3)
!    INTEGER                          :: Node1, Node2
!    INTEGER                          :: VibQuant,VibQuantWall,VibQuantNew
!    REAL                             :: VibQuantNewR
!    REAL                             :: nx, ny, nz, nVal
!    REAL                             :: xNod, yNod, zNod, VeloReal, VeloCrad, VeloCx, VeloCy ,VeloCz
!    REAL                             :: EtraOld, EtraNew, RanNum, Cmr, Phi, Fak_D
!    REAL                             :: EtraWall, ErotWall, EvibNew, ErotNew
!    REAL                             :: VelX, VelY, VelZ, VecX, VecY, VecZ
!    REAL                             :: VelXold, VelYold, VelZold
!    REAL                             :: Vector1(1:3), Vector2(1:3)
!    REAL, PARAMETER                  :: PI=3.14159265358979323846
!    REAL                             :: IntersectionPos(1:3)
!    REAL                             :: TransArray(1:6),IntArray(1:6)
! !===================================================================================================================================
!   TransArray(:) = 0.0
!   IntArray(:) = 0.0
!   globSide = PDM%PartAdsorbSideIndx(1,i)
!   SurfSide = SurfMesh%GlobSideToSurfSideMap(globSide)
!     
!   TransACC = 1
!   VibACC = 1
!   RotACC = 1
!   
!   CALL RANDOM_NUMBER(RanNum)
!   IF (Adsorption%AdsorpInfo(PartSpecies(i))%ProbDes(SurfSide).GE.RanNum) THEN 
!  
!     PDM%ParticleAtWall(i)=.FALSE.
!     ! Sample desorping particles
! #if (PP_TimeDiscMethod==42)
!     Adsorption%AdsorpInfo(PartSpecies(i))%NumOfDes(SurfSide) = Adsorption%AdsorpInfo(PartSpecies(i))%NumOfDes(SurfSide) + 1
! #endif
!     
!     maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(SurfSide) !* Adsorption%MaxCoverage(:,iSpec)
!     Adsorption%Coverage(SurfSide,PartSpecies(i)) = Adsorption%Coverage(SurfSide,PartSpecies(i)) & 
!                                                   - Species(PartSpecies(i))%MacroParticleFactor/maxPart
!     WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
!     WallVelo(1:3) = PartBound%WallVelo(1:3,PartBound%MapToPartBC(BC(globSide)))
!     
!     Element = SideToElem(1,globSide)
!     IF (Element.LT.1) THEN
!       Element = SideToElem(2,globSide)
!       iLocSide = SideToElem(4,globSide)
!     ELSE
!       iLocSide = SideToElem(3,globSide)
!     END IF 
!       
!     xNod = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
!     yNod = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
!     zNod = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
! 
!     !---- Calculate normal vector in global coordinates : (points into Element)
!     TriNum = PDM%PartAdsorbSideIndx(2,i)
!     
!     Node1 = TriNum+1     ! normal = cross product of 1-2 and 1-3 for first triangle
!     Node2 = TriNum+2     !          and 1-3 and 1-4 for second triangle
! 
!     Vector1(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - xNod
!     Vector1(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - yNod
!     Vector1(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node1,iLocSide,Element)) - zNod
! 
!     Vector2(1) = GEO%NodeCoords(1,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - xNod
!     Vector2(2) = GEO%NodeCoords(2,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - yNod
!     Vector2(3) = GEO%NodeCoords(3,GEO%ElemSideNodeID(Node2,iLocSide,Element)) - zNod
! 
!     nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2)
!     ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
!     nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)
! 
!     nVal = SQRT(nx*nx + ny*ny + nz*nz)
! 
!     nx = nx/nVal
!     ny = ny/nVal
!     nz = nz/nVal
!     
!     !---- Calculate new velocity vector (Extended Maxwellian Model)
!   !    VeloReal = SQRT(PartState(i,4) * PartState(i,4) + &
!   !              PartState(i,5) * PartState(i,5) + &
!   !              PartState(i,6) * PartState(i,6))
!   !    EtraOld = 0.5 * Species(PartSpecies(i))%MassIC * VeloReal**2
!     VeloReal = 0.0
!     EtraOld = 0.0
!     CALL RANDOM_NUMBER(RanNum)
!     VeloCrad    = SQRT(-LOG(RanNum))
!     CALL RANDOM_NUMBER(RanNum)
!     VeloCz      = SQRT(-LOG(RanNum))
!     Fak_D       = VeloCrad**2 + VeloCz**2
!     EtraWall    = BoltzmannConst * WallTemp * Fak_D
!   !    EtraNew = EtraOld + TransACC * (EtraWall - EtraOld)
!     EtraNew = EtraWall
!     Cmr     = SQRT(2.0 * EtraNew / (Species(PartSpecies(i))%MassIC * Fak_D))
!     CALL RANDOM_NUMBER(RanNum)
!     Phi     = 2 * PI * RanNum
!     VeloCx  = Cmr * VeloCrad * COS(Phi)
!     VeloCy  = Cmr * VeloCrad * SIN(Phi)
!     VeloCz  = Cmr * VeloCz
!     
!     !---- Transformation local distribution -> global coordinates
!     VecX = Vector1(1) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
!     VecY = Vector1(2) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
!     VecZ = Vector1(3) / SQRT( Vector1(1)**2 + Vector1(2)**2 + Vector1(3)**2 )
!     
!     VelX = VecX*VeloCx + (nz*VecY-ny*VecZ)*VeloCy - nx*VeloCz
!     VelY = VecY*VeloCx + (nx*VecZ-nz*VecX)*VeloCy - ny*VeloCz
!     VelZ = VecZ*VeloCx + (ny*VecX-nx*VecY)*VeloCy - nz*VeloCz
! 
!     lastPartPos(i,1) = PartState(i,1)
!     lastPartPos(i,2) = PartState(i,2)
!     lastPartPos(i,3) = PartState(i,3)
!     VelXold = PartState(i,4)
!     VelYold = PartState(i,5)
!     VelZold = PartState(i,6)
!     
!     CALL RANDOM_NUMBER(RanNum)
!     PartState(i,1)   = Partstate(i,1) + RanNum * dt * VelX
!     PartState(i,2)   = Partstate(i,2) + RanNum * dt * VelY
!     PartState(i,3)   = Partstate(i,3) + RanNum * dt * VelZ
!     PartState(i,4) = PartState(i,4) + VelX
!     PartState(i,5) = PartState(i,5) + VelY
!     PartState(i,6) = PartState(i,6) + VelZ
!    
!     TransArray(1) = EtraOld
!     TransArray(2) = EtraWall
!     TransArray(3) = EtraNew
!     TransArray(4) = PartState(i,4)-VelXold
!     TransArray(5) = PartState(i,5)-VelYold
!     TransArray(6) = PartState(i,6)-VelZold
!     
!     !----  Sampling at walls
!     IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)) .OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
!       CALL CalcWallSample(i,SurfSide,TriNum,Transarray,IntArray)
!     END IF
! 
!  ELSE 
!    LastPartPos(i,1) = PartState(i,1)
!    LastPartPos(i,2) = PartState(i,2)
!    LastPartPos(i,3) = PartState(i,3)
!  END IF !End Desorption
! 
! END SUBROUTINE Particle_Wall_Desorb

SUBROUTINE CalcAdsorbProb()
!===================================================================================================================================
! Models for adsorption probability calculation
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_DSMC_Vars,              ONLY : Adsorption
#if (PP_TimeDiscMethod==42)
USE MOD_DSMC_Vars,              ONLY : DSMC
#endif
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
  INTEGER                          :: SurfSide, iSpec, p, q
  REAL                             :: Theta_req, Kfactor, S_0
!=================================================================================================================================== 
  DO iSpec=1,nSpecies
  DO SurfSide=1,SurfMesh%nSides
  DO q = 1,nSurfSample
  DO p = 1,nSurfSample
!===================================================================================================================================
! Kisluik Sticking Model from Kolasinski's Surface Science (book)
    ! enhance later to co-adsorption
    Theta_req = (1.0 - Adsorption%Coverage(p,q,SurfSide,iSpec)/Adsorption%MaxCoverage(SurfSide,iSpec)) &
              **Adsorption%Adsorbexp(SurfSide,iSpec)  
    !----- kann später auf von Wandtemperatur abhängige Werte erweitert werden          
    Kfactor = Adsorption%PrefactorStick(SurfSide,iSpec)
    S_0 = Adsorption%InitStick(SurfSide,iSpec)
    !-----
    IF (Theta_req.EQ.0) THEN
      Adsorption%ProbAds(p,q,SurfSide,iSpec) = 0.
    ELSE
      Adsorption%ProbAds(p,q,SurfSide,iSpec) = S_0 / (1.0 + Kfactor * ( 1.0/Theta_req - 1.0))
    END IF
!=================================================================================================================================== 
! ! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP)  
!     Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)  
!     IF (Theta.GT.0) THEN
!       Q_0A = 26312.
!       D_A  = 59870.
!       n = 2
!       m = Theta * n
!       IF (m.LT.1) THEN
!   !       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A)
!         Q = (Q_0A**2)/(Q_0A/n+D_A)
!       ELSE 
!         sigma = 1/m*(2-1/m)
!   !       Q = (9*Q_0A**2)/(6*Q_0A+16*D_A) *sigma
!         Q = (Q_0A**2)/(Q_0A/n+D_A) * sigma
!       END IF
!       
!       CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct)
!       CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
!       
!       E_des = Q
!       nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
!       rate = nu_des * exp(-E_des/WallTemp)
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = rate * dt
!     ELSE
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = 0.0
!     END IF
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbAds + Adsorption%ProbAds(p,q,SurfSide,iSpec)
    END IF
#endif
  END DO
  END DO
  END DO
  END DO
#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbAds = Adsorption%AdsorpInfo(:)%MeanProbAds / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
#endif
END SUBROUTINE CalcAdsorbProb

SUBROUTINE CalcDesorbProb()
!===================================================================================================================================
! Models for desorption probability calculation
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
  USE MOD_TimeDisc_Vars,          ONLY : dt
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
  USE MOD_DSMC_Vars,              ONLY : DSMC
#endif
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSide, iSpec, globSide, p, q
   REAL                             :: Theta, E_des, nu_des, rate, WallTemp
   REAL                             :: Q_0A, D_A, sigma, Heat, m
   INTEGER                          :: n
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
!===================================================================================================================================
DO SurfSide=1,SurfMesh%nSides
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSide)
! special TPD (temperature programmed desorption) temperature adjustment routine    
#if (PP_TimeDiscMethod==42)
  IF (Adsorption%TPD) THEN
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide))) + (Adsorption%TPD_beta * iter * dt)
    Adsorption%TPD_Temp = Walltemp
  ELSE
    WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
  END IF
#else
  WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
#endif

  DO iSpec = 1,nSpecies 
  DO q = 1,nSurfSample
  DO p = 1,nSurfSample
!===================================================================================================================================
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP)  
    Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)  
    IF (Theta.GT.0) THEN
      Q_0A = 26312.
      D_A  = 59870.
      n = 2
      m = Theta * n
      IF (m.LT.1) THEN
        Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A)
!         Heat = (Q_0A**2)/(Q_0A/n+D_A)
      ELSE 
        sigma = 1/m*(2-1/m)
        Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A) *sigma
!         Heat = (Q_0A**2)/(Q_0A/n+D_A) * sigma
      END IF
      
      CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSide))
      CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
      
      E_des = Heat
      nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
      rate = nu_des * exp(-E_des/WallTemp)
      Adsorption%ProbDes(p,q,SurfSide,iSpec) = rate * dt
    ELSE
      Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.0
    END IF
!===================================================================================================================================
! absolute rate theory / lattice gas model (QCA-quasi-chemical-approximation)
!===================================================================================================================================
! Polanyi-Wigner-eq. from Kolasinski's Surface Science (book) 
!       Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)! / Adsorption%MaxCoverage(SurfSide,iSpec)
!       !----- kann später auf von Wandtemperatur/Translationsenergie abhängige Werte erweitert werden          
!       E_des = Adsorption%DesorbEnergy(SurfSide,iSpec) + Adsorption%Intensification(SurfSide,iSpec) * Theta
!       nu_des = 10**(Adsorption%Nu_a(SurfSide,iSpec) + Adsorption%Nu_b(SurfSide,iSpec) * Theta)/10000
!       !-----
!       rate = nu_des *(Adsorption%DensSurfAtoms(SurfSide)**(Adsorption%Adsorbexp(SurfSide,iSpec)-1)) &
!                     * (Theta**Adsorption%Adsorbexp(SurfSide,iSpec)) * exp(-E_des/WallTemp)
!       IF (Theta.GT.0) THEN
!         Adsorption%ProbDes(p,q,SurfSide,iSpec) = rate * dt
!       ELSE
!         Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.0
!       END IF
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + Adsorption%ProbDes(p,q,SurfSide,iSpec)
    END IF
#endif
  END DO
  END DO
  END DO
END DO
#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbDes = Adsorption%AdsorpInfo(:)%MeanProbDes / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
#endif
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

SUBROUTINE PartitionFuncAct(iSpec, Temp, VarPartitionFuncAct, Surfdensity)
!===================================================================================================================================
! Partitionfunction of activated complex
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,       ONLY : PlanckConst
USE MOD_DSMC_Vars,          ONLY : SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY : BoltzmannConst, Species
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
REAL, INTENT(IN)              :: Surfdensity
!===================================================================================================================================
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncAct
!===================================================================================================================================
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
  Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))/Surfdensity
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
USE MOD_Globals_Vars,       ONLY : PlanckConst
USE MOD_DSMC_Vars,          ONLY : SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY : BoltzmannConst
!===================================================================================================================================
IMPLICIT NONE
!===================================================================================================================================
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!===================================================================================================================================
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncSurf
!===================================================================================================================================
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib
!===================================================================================================================================
  Qtra = 1
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
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

END MODULE MOD_DSMC_SurfModel_Tools
