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
PUBLIC                       :: CalcBackgndPartDesorb
PUBLIC                       :: CalcAdsorbProb
PUBLIC                       :: CalcDesorbProb
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_Update_Wall_Vars()   
!===================================================================================================================================
! Update and sample DSMC-values for adsorption, desorption and reactions on surfaces
!===================================================================================================================================
  USE MOD_PARTICLE_Vars,          ONLY : nSpecies, PDM
  USE MOD_PARTICLE_Vars,          ONLY : KeepWallParticles, PEM, Species
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
!       Adsorption%AdsorpInfo(iSpec)%MeanEAds = 0.
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

    IF (DSMC%WallModel.EQ.2) THEN
      CALL CalcDistNumChange()
      CALL CalcDiffusion()
      CALL CalcSurfDistInteraction()
    ELSE IF (DSMC%WallModel.EQ.3) THEN
!       CALL CalcDistNumChange()
      CALL CalcDiffusion()
    END IF
    
    CALL CalcAdsorbProb()
    IF (KeepWallParticles) CALL CalcDesorbprob()
  END IF

END SUBROUTINE DSMC_Update_Wall_Vars  
    
SUBROUTINE Particle_Wall_Adsorb(PartTrajectory,lengthPartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,adsindex) 
!===================================================================================================================================
! Particle Adsorption after wall collision
!===================================================================================================================================
  USE MOD_DSMC_Analyze,           ONLY : CalcWallSample
  USE MOD_Particle_Vars,          ONLY : WriteMacroValues, KeepWallParticles, PDM
  USE MOD_Particle_Vars,          ONLY : PartState,LastPartPos,Species,BoltzmannConst,PartSpecies,nSpecies
  USE MOD_Mesh_Vars,              ONLY : ElemToSide, BC
  USE MOD_DSMC_Vars,              ONLY : useDSMC, CollisMode, Adsorption
  USE MOD_DSMC_Vars,              ONLY : PartStateIntEn, SpecDSMC, DSMC, PolyatomMolDSMC, VibQuantsPar
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, dXiEQ_SurfSample, Partbound
  USE MOD_TimeDisc_Vars,          ONLY : TEnd, time
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! argument list declaration
  INTEGER,INTENT(INOUT)            :: adsindex
  REAL,INTENT(INOUT)               :: PartTrajectory(1:3), lengthPartTrajectory, alpha
  REAL,INTENT(IN)                  :: xi, eta
  INTEGER,INTENT(IN)               :: PartID, GlobSideID
  LOGICAL,INTENT(IN)               :: IsSpeciesSwap
! LOCAL VARIABLES
  INTEGER                              :: locBCID, VibQuant, VibQuantNew, VibQuantWall
  REAL                                 :: VibQuantNewR
  INTEGER                              :: ElemID
  REAL                                 :: v_help(3)
  REAL                                 :: oneMinus                
  REAL                                 :: VeloReal, RanNum, EtraOld, VeloCrad, Fak_D
  REAL                                 :: EtraWall, EtraNew
  REAL                                 :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
  REAL                                 :: n_loc(1:3), tang1(1:3),tang2(1:3), NewVelo(3)
  REAL                                 :: ErotNew, ErotWall, EVibNew, Phi, Cmr, VeloCx, VeloCy, VeloCz
  REAL                                 :: Xitild,EtaTild
  INTEGER                              :: p,q, SurfSideID
  REAL                                 :: POI_fak, TildPos(3),TildTrajectory(3)
! Polyatomic Molecules
  REAL, ALLOCATABLE                    :: RanNumPoly(:), VibQuantNewRPoly(:)
  INTEGER                              :: iPolyatMole, iDOF
  INTEGER, ALLOCATABLE                 :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
  REAL                                 :: ErotOldPoly(3), ErotNewPoly(3), ErotWallPoly(3), CmrRotPoly(3)
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
          IntArray(5) = IntArray(5) + VibQuantWall * BoltzmannConst &
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
      CALL CalcWallSample(PartID,SurfSide,p,q,Transarray,IntArray,PartTrajectory,lengthPartTrajectory,alpha,IsSpeciesSwap)
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
USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
!    INTEGER                          :: i, PartAds2
   INTEGER                          :: iSurfSide, iSpec, p, q, NPois, WallPartNum, iProbSigma
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
            PartAds = Adsorption%Coverage(p,q,iSurfSide,iSpec) * Adsorption%DensSurfAtoms(iSurfSide) &
                        * SurfMesh%SurfaceArea(p,q,iSurfSide) &
                        / Species(iSpec)%MacroParticleFactor
            
            IF ((DSMC%WallModel.EQ.1) .OR. (DSMC%WallModel.EQ.3)) THEN
              PartDes = PartAds * Adsorption%ProbDes(p,q,iSurfSide,iSpec)
            ELSE IF (DSMC%WallModel.EQ.2) THEN 
              PartDes = 0.
              DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
                PartDes = PartDes + PartAds * Adsorption%ProbSigma(p,q,iSurfSide,iSpec,iProbSigma) &
                                            * Adsorption%ProbSigDes(p,q,iSurfSide,iSpec,iProbSigma) 
              END DO
            END IF
              IF (PartDes.GT.PartAds) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) + INT(PartAds)
              ELSE
                CALL RANDOM_NUMBER(RanNum)            
                IF (EXP(-PartDes).LE.TINY(PartDes) &
#if (PP_TimeDiscMethod==42)
                .OR. Adsorption%TPD &
#endif
                ) THEN
                  Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                  +INT(PartDes + RanNum)
                ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
                  Npois=0
                  Tpois=1.0
                  DO
                    Tpois=RanNum*Tpois
                    IF (Tpois.LT.TINY(Tpois)) THEN
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)&
                      +INT(PartDes + RanNum)
                      EXIT
                    END IF
                    IF (Tpois.GT.EXP(-PartDes)) THEN
                      Npois=Npois+1
                      CALL RANDOM_NUMBER(RanNum)
                    ELSE
                      Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec)+Npois
                      EXIT
                    END IF
                  END DO
                END IF
              END IF !PartDes.GT.WallPartNum
              IF (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).GT.WallPartNum  &
              .OR.Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec).LT.0) THEN
                Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) = WallPartNum
              END IF
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

SUBROUTINE CalcDiffusion()
!===================================================================================================================================
! Model for diffusion calculation
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, globSide, subsurfeta, subsurfxi
   INTEGER                          :: Coord, nSites, nSitesRemain, i, j, AdsorbID
   REAL                             :: WallTemp, Prob_diff, RanNum
   REAL                             :: Q_0A, D_A, Heat_i, Heat_j, Heat_temp, sigma
   INTEGER                          :: bondorder, n_equal_site_Neigh, Indx, Indy, Surfpos, newpos
   INTEGER , ALLOCATABLE            :: free_Neigh_pos(:)
!===================================================================================================================================
D_A  = 59870.

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO iSpec=1,nSpecies
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    Q_0A = Adsorption%HeatOfAdsZero(iSpec)
    
    ALLOCATE ( free_Neigh_pos(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
      n_equal_site_Neigh = 0
      free_Neigh_pos(:) = 0
      
      ! find free Neighbour positions of the same site-coordination
      DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i) .EQ. Coord) THEN
          IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)).EQ.0) ) THEN
            n_equal_site_Neigh = n_equal_site_Neigh + 1
            free_Neigh_pos(n_equal_site_Neigh) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)
          END IF
        END IF
      END DO
      
      ! calculate heat of adsorption for actual site and reduce bond order of bondatoms
      sigma = 0.
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
        
        bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy)
        sigma = sigma + ( (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom) &
              * (1./(bondorder)) * (2.-(1./bondorder)) )
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
      END DO
      Heat_i = (Q_0A**2)/(Q_0A/SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom+D_A) * sigma
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
      
      ! choose Neighbour position with highest heat of adsorption if adsorbate would move there
      Heat_j = 0.
      DO i = 1,n_equal_site_Neigh
        sigma = 0.
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(free_Neigh_pos(i),j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(free_Neigh_pos(i),j) 
          bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
          sigma = sigma + ( (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom) &
                * (1./(bondorder)) * (2.-(1./bondorder)) )
        END DO
        Heat_temp = (Q_0A**2)/(Q_0A/SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom+D_A) * sigma
        IF (Heat_temp .GT. Heat_j) THEN
          Heat_j = Heat_temp
          newpos = free_Neigh_pos(i)
        END IF
      END DO
      
      ! only try to diffuse particle if unoccupied sites available
      IF (n_equal_site_Neigh .GE. 1) THEN
        globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
        WallTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(globSide)))
        Prob_diff = exp(-(Heat_i - Heat_j)/WallTemp) / (1+exp(-(Heat_i - Heat_j)/Walltemp))
        CALL RANDOM_NUMBER(RanNum)
        IF (Prob_diff.GT.RanNum) THEN
        ! move particle to new position and update map
          DO i = 1,nSitesRemain
            IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i).EQ.newpos) THEN
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = newpos
            END IF
          END DO
        ELSE 
          newpos = Surfpos
        END IF
      ELSE
        newpos = Surfpos        
      END IF ! end if (n_equal_site_Neigh >= 1)
      
      ! update surfatom bond order and species map
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(newpos) = iSpec
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(newpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(newpos,j) 
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
      END DO
    END DO
    DEALLOCATE(free_Neigh_pos)
    
  END DO
  
END DO
END DO
END DO

END SUBROUTINE CalcDiffusion

SUBROUTINE CalcBackgndPartDesorb()
!===================================================================================================================================
! Routine for calculation of number of desorbing particles from background surface distribution
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
  USE MOD_TimeDisc_Vars,          ONLY : dt
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
#endif
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, globSide, subsurfeta, subsurfxi
   INTEGER                          :: Coord, i, j, AdsorbID, numSites
   REAL                             :: WallTemp, Prob_diff, RanNum
   REAL                             :: Heat_A, Heat_B, Heat_AB, D_AB, sigma, nu_des, rate, P_actual_des
   REAL , ALLOCATABLE               :: P_des(:), Energy(:)
   INTEGER                          :: bondorder, Indx, Indy, Surfpos
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
   INTEGER                          :: trace, traceNum
   INTEGER , ALLOCATABLE            :: desorbnum(:), adsorbnum(:), nSites(:), nSitesRemain(:), remainNum(:), adsorbates(:)
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
!   Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(:)%NumOfDes = 0
#endif

ALLOCATE (desorbnum(1:nSpecies),&
adsorbnum(1:4),&
nSites(1:4),&
nSitesRemain(1:4),&
remainNum(1:4),&
P_des(1:nSpecies),&
adsorbates(1:nSpecies),&
Energy(1:nSpecies))

#if (PP_TimeDiscMethod==42)
Adsorption%AdsorpInfo(:)%MeanEAds = 0.
#endif

DO SurfSideID = 1,SurfMesh%nSides
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
  
! special TPD (temperature programmed desorption) temperature adjustment part
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

DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample
  desorbnum(:) = 0
  nSites(:) = 0
  nSitesRemain(:) = 0
  remainNum(:) = 0
  P_des(:) = 0.
  Energy(:) = 0.
  
  DO Coord=1,3
    nSites(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    nSites(4) = nSites(4) + nSites(Coord)
    nSitesRemain(4) = nSitesRemain(4) + nSitesRemain(Coord)
  END DO
  adsorbnum(:) = nSites(:) - nSitesRemain(:)
  traceNum = adsorbnum(4)
  adsorbates(:) = 0
  
  IF (adsorbnum(4) .GT. 0) THEN
  DO trace = 1,traceNum
    CALL RANDOM_NUMBER(RanNum)
    AdsorbID = nSitesRemain(4) + 1 + INT((adsorbnum(4)-remainNum(4))*RanNum)
    ! choose right coordination and adjust adsorbate ID
    IF ( AdsorbID.GT.(adsorbnum(1)+adsorbnum(2)+nSitesRemain(4)-remainNum(1)-remainNum(2)) ) THEN
      AdsorbID = AdsorbID - (adsorbnum(1) + adsorbnum(2) + nSitesRemain(1) + nSitesRemain(2) - remainNum(1) - remainNum(2))
      Coord = 3
    ELSE IF ( (AdsorbID.GT.(adsorbnum(1)+nSitesRemain(4)-remainNum(1))) &
             .AND. (AdsorbID.LE.(adsorbnum(1)+adsorbnum(2)+nSitesRemain(4))-remainNum(1)-remainNum(2)) ) THEN
      AdsorbID = AdsorbID - (adsorbnum(1) + nSitesRemain(1) + nSitesRemain(3) - remainNum(1))
      Coord = 2
    ELSE IF (AdsorbID.LE.adsorbnum(1)) THEN
      AdsorbID = AdsorbID - (nSitesRemain(2) + nSitesRemain(3))
      Coord = 1
    ELSE
      CYCLE
    END IF
    
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    
    ! calculate heat of adsorption for actual site and calculate desorption probability
    sigma = 0.
    DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
      Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
      Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
      bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy)
      sigma = sigma + ( (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom) &
            * (1./(bondorder)) * (2.-(1./bondorder)) )
    END DO
    
    Heat_A = Adsorption%HeatOfAdsZero(iSpec) * sigma
    Energy(iSpec) = Energy(iSpec) + Heat_A
    CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID))
    CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
    nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
    rate = nu_des * exp(-Heat_A/WallTemp)
    P_actual_des = rate * dt
    IF (rate*dt.GT.1) THEN
      P_des(iSpec) = P_des(iSpec) + 1
    ELSE
      P_des(iSpec) = P_des(iSpec) + P_actual_des
    END IF
    adsorbates(iSpec) = adsorbates(iSpec) + 1
    
!     ! find Neighbour position with reaction partner
!     DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
!     DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
!       IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species( &
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
!           .EQ.Adsorption%AssocReact(1,ReactNum,iSpec)) ) THEN
!         n_react_Neigh = n_react_Neigh + 1
!         D_AB(n_react_Neigh) = Adsorption%EDissBond((Adsorption%DissNum+ReactNum),iSpec)
!         ReactNum_partner(n_react_Neigh) = ReactNum
!         react_Neigh_pos(n_react_Neigh) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)
!         ! calculate heat of adsorption of product molecule for current site 
!         Heat_AB(n_react_Neigh) = Adsorption%HeatOfAdsZero(Adsorption%AssocReact(2,ReactNum,iSpec)) * sigma
!       END IF
!     END DO
!     END DO
!     ! calculate heat of adsorption for reaction partner site
!     sigma = 0.
!     Heat_j = 0.
!     CALL RANDOM_NUMBER(RanNum)
!     react_Neigh = 1 + INT(n_react_Neigh * RanNum)
!     jSpec = (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(react_Neigh_pos(react_Neigh))
!     DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
!       Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(react_Neigh_pos(react_Neigh),j)
!       Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(react_Neigh_pos(react_Neigh),j)
!       bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy)
!       sigma = sigma + ( (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom) &
!             * (1./(bondorder)) * (2.-(1./bondorder)) )
!     END DO
!     Heat_B = Adsorption%HeatOfAdsZero(jSpec) * sigma
!     ! calculate LH reaction probabilities
!     E_a = 0.5 * ( D_AB(react_Neigh) + Heat_AB(react_Neigh) - (Heat_A + Heat_B) + ((Heat_A * Heat_B)/(Heat_A + Heat_B)))
!     IF (E_a .GT. 0) THEN
!       E_d = (Heat_A + Heat_B) - D_AB(react_Neigh) - E_a
!     ELSE
!       E_d = (Heat_A + Heat_B) - D_AB(react_Neigh)
!     END IF
!     CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID))
!     CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1)
!     CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2)
    
!       IF ( ReactionProb / ( ReactionProb + ReactionProb2) .gt. iRan) THEN
    ! only try to desorb particle if random number higher probability and then update map
    CALL RANDOM_NUMBER(RanNum)
    IF (P_actual_des.GT.RanNum) THEN
      desorbnum(iSpec) = desorbnum(iSpec) + 1
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
      DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
        Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1)
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+1) = Surfpos
      nSitesRemain(Coord) = nSitesRemain(Coord) + 1
      nSitesRemain(4) = nSitesRemain(4) + 1
!       SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
!           SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord))
!       SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)) = Surfpos
!       remainNum(Coord) = remainNum(Coord) + 1
!       remainNum(4) = remainNum(4) + 1
    ELSE
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord))
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)) = Surfpos
      remainNum(Coord) = remainNum(Coord) + 1
      remainNum(4) = remainNum(4) + 1
    END IF
    ! update number of adsorbates
    adsorbnum(:) = nSites(:) - nSitesRemain(:)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain(Coord)
  END DO
  END IF
  
  DO iSpec = 1,nSpecies    
  numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec))
  Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec) = INT((REAL(desorbnum(iSpec)) / REAL(numSites)) &
                                                                  * (Adsorption%DensSurfAtoms(SurfSideID) &
                      * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID) / Species(iSpec)%MacroParticleFactor))
  IF (adsorbates(iSpec).EQ.0) THEN
    Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = 0
    Energy(iSpec) = Adsorption%HeatOfAdsZero(iSpec)
  ELSE
    Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = P_des(iSpec) / REAL(adsorbates(iSpec))
    Energy(iSpec) = Energy(iSpec) / REAL(adsorbates(iSpec))
  END IF
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%NumOfDes = Adsorption%AdsorpInfo(iSpec)%NumOfDes &
                                        + Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec)
  Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
                                            + Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec)
  Adsorption%AdsorpInfo(iSpec)%MeanEAds = Adsorption%AdsorpInfo(iSpec)%MeanEAds + Energy(iSpec)
#endif
  END DO
END DO
END DO
END DO

#if (PP_TimeDiscMethod==42)
IF (.NOT.DSMC%ReservoirRateStatistic) THEN
  Adsorption%AdsorpInfo(:)%MeanProbDes = Adsorption%AdsorpInfo(:)%MeanProbDes / (nSurfSample * nSurfSample * SurfMesh%nSides)
  Adsorption%AdsorpInfo(:)%MeanEAds = Adsorption%AdsorpInfo(:)%MeanEAds / (nSurfSample * nSurfSample * SurfMesh%nSides)
END IF
#endif

DEALLOCATE(desorbnum,adsorbnum,nSites,nSitesRemain,remainNum,P_des,adsorbates,Energy)

END SUBROUTINE CalcBackgndPartDesorb

SUBROUTINE CalcDistNumChange()
!===================================================================================================================================
! updates mapping of adsorbed particles if particles desorbed (sumdesorbpart.GT.0)
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSideID, iSpec, subsurfeta, subsurfxi, Adsorbates
   INTEGER                          :: Coord, nSites, nSitesRemain, iInterAtom, xpos, ypos
   REAL                             :: RanNum
   INTEGER                          :: Surfpos, UsedSiteMapPos, dist, Surfnum, newAdsorbates
   INTEGER , ALLOCATABLE            :: free_Neigh_pos(:)
!===================================================================================================================================

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO iSpec = 1,nSpecies
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    Adsorbates = INT(Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,iSpec) &
                * SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec)))
    IF (nSites-nSitesRemain.GT.Adsorbates) THEN
      newAdsorbates = (nSites-nSitesRemain) - Adsorbates
      ! remove adsorbates randomly from the surface from the correct site and assign surface atom bond order
      dist = 1
      Surfnum = nSites - nSitesRemain
      DO WHILE (dist.LE.newAdsorbates) 
        CALL RANDOM_NUMBER(RanNum)
        Surfpos = 1 + INT(Surfnum * RanNum)
        UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = 0
        ! assign bond order of respective surface atoms in the surfacelattice
        DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
          ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) - 1
        END DO
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+Surfpos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+1) = UsedSiteMapPos
        Surfnum = Surfnum - 1
        nSitesRemain = nSitesRemain + 1
        dist = dist + 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain
    ELSE IF (nSites-nSitesRemain.LT.Adsorbates) THEN
      newAdsorbates = Adsorbates - (nSites-nSitesRemain)
      ! add new Adsorbates to macro adsorbate map
      dist = 1
      Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
      DO WHILE (dist.LE.newAdsorbates) 
        CALL RANDOM_NUMBER(RanNum)
        Surfpos = 1 + INT(Surfnum * RanNum)
        UsedSiteMapPos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(UsedSiteMapPos) = iSpec
        ! assign bond order of respective surface atoms in the surfacelattice
        DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(UsedSiteMapPos,iInterAtom)
          ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(UsedSiteMapPos,iInterAtom)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos) + 1
        END DO
        ! rearrange UsedSiteMap-Surfpos-array
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfpos) = &
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum)
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(Surfnum) = UsedSiteMapPos
        Surfnum = Surfnum - 1
        dist = dist + 1
      END DO
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = Surfnum
    ELSE
      CYCLE
    END IF
  END DO

END DO
END DO
END DO

END SUBROUTINE CalcDistNumChange

SUBROUTINE CalcSurfDistInteraction()
!===================================================================================================================================
! Model for calculating surface distibution and effects of coverage on heat of adsorption
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies
  USE MOD_DSMC_Vars,              ONLY : Adsorption, SurfDistInfo
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!=================================================================================================================================== 
! Local variable declaration
  INTEGER                          :: SurfSideID, subsurfxi, subsurfeta, iSpec, counter
  REAL                             :: sigmasub 
  REAL, ALLOCATABLE                :: sigma(:,:,:,:), ProbSigma(:,:,:,:)
  INTEGER                          :: firstval, secondval, thirdval, fourthval, pos
  INTEGER                          :: Surfpos, nSites, nSitesRemain, xpos, ypos, AdsorbID, iInterAtom, bondorder, Coord
!===================================================================================================================================
ALLOCATE( sigma(1:36*nSpecies,1:36*nSpecies,1:36*nSpecies,1:36*nSpecies),&
          ProbSigma(1:36*nSpecies,1:36*nSpecies,1:36*nSpecies,1:36*nSpecies) )

DO SurfSideID=1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO iSpec=1,nSpecies
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    
    pos = 1
    sigma = 1.
    ProbSigma = 0.
    firstval = 1
    secondval = 1
    thirdval = 1
    fourthval = 1
    counter = 0
      
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)

      sigmasub = 0.
      DO iInterAtom = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
        xpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,iInterAtom)
        ypos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,iInterAtom)
        bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,xpos,ypos)
        sigmasub = sigmasub + (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom &
                  * 1./bondorder * (2.-1./bondorder))
        IF (bondorder.GE.firstval) THEN
          fourthval = thirdval
          thirdval = secondval
          secondval = firstval
          firstval = bondorder
        ELSE IF (bondorder.GE.secondval) THEN
          fourthval = thirdval
          thirdval = secondval
          secondval = bondorder
        ELSE IF (bondorder.GE.thirdval) THEN
          fourthval = thirdval
          thirdval = bondorder
        ELSE
          fourthval = bondorder
        END IF
      END DO
      sigma(firstval,secondval,thirdval,fourthval) = sigmasub
      ProbSigma(firstval,secondval,thirdval,fourthval) = ProbSigma(firstval,secondval,thirdval,fourthval) + 1.
      counter = counter + 1
    END DO
    
    IF (counter.GT.0) THEN
      firstval = 1
      secondval = 1
      thirdval = 1
      fourthval = 1
      DO WHILE ((firstval.LE.4).AND.(secondval.LE.4).AND.(thirdval.LE.4).AND.(fourthval.LE.4))
        DO firstval = 1,4
        DO secondval = 1,firstval
        DO thirdval = 1,secondval
        DO fourthval = 1,thirdval
          Adsorption%Sigma(subsurfxi,subsurfeta,SurfSideID,iSpec,pos) = sigma(firstval,secondval,thirdval,fourthval)
          Adsorption%ProbSigma(subsurfxi,subsurfeta,SurfSideID,iSpec,pos) = ProbSigma(firstval,secondval,thirdval,fourthval)/counter
          pos = pos + 1
        END DO
        END DO
        END DO
        END DO
      END DO
    ELSE
      Adsorption%Sigma(subsurfxi,subsurfeta,SurfSideID,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 1.
      Adsorption%ProbSigma(subsurfxi,subsurfeta,SurfSideID,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 0.
    END IF
    
  END DO! iSpec = 1,nSpecies
    
END DO
END DO
END DO

DEALLOCATE(sigma,ProbSigma)

END SUBROUTINE CalcSurfDistInteraction


SUBROUTINE CalcAdsorbProb()
!===================================================================================================================================
! Models for adsorption probability calculation
!===================================================================================================================================
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
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
  IF (DSMC%WallModel.EQ.1) THEN
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
  ELSE IF (DSMC%WallModel.EQ.2) THEN
!=================================================================================================================================== 
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP) 

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
!       CALL PartitionFuncGas(iSpec, WallTemp, VarPartitionFuncGas)
!       
!       E_des = Q
!       nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
!       rate = nu_des * exp(-E_des/WallTemp)
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = rate * dt
!     ELSE
!       Adsorption%AdsorpInfo(iSpec)%ProbDes(p,q,SurfSide) = 0.0
!     END IF
!===================================================================================================================================
  END IF
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
  USE MOD_Particle_Vars,          ONLY : nSpecies, Species, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
  USE MOD_TimeDisc_Vars,          ONLY : dt
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
#endif
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: SurfSide, iSpec, globSide, p, q
   REAL                             :: Theta, nu_des, rate, WallTemp
   REAL                             :: Q_0A, D_A, m, Heat, E_des, sigma
!    REAL                :: sigma(10)
   INTEGER                          :: n, iProbSigma
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
!===================================================================================================================================
! CALL CalcSurfDistInteraction()
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
  IF (DSMC%WallModel.EQ.1) THEN
!===================================================================================================================================
! Polanyi-Wigner-eq. from Kolasinski's Surface Science (book) 
!===================================================================================================================================
    Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)! / Adsorption%MaxCoverage(SurfSide,iSpec)
    !----- kann später auf von Wandtemperatur/Translationsenergie abhängige Werte erweitert werden          
    E_des = Adsorption%DesorbEnergy(SurfSide,iSpec) + Adsorption%Intensification(SurfSide,iSpec) * Theta
    nu_des = 10**(Adsorption%Nu_a(SurfSide,iSpec) + Adsorption%Nu_b(SurfSide,iSpec) * Theta)!/10000
    !-----
    rate = nu_des &!*(Adsorption%DensSurfAtoms(SurfSide)**(Adsorption%Adsorbexp(SurfSide,iSpec)-1)) &
                  * (Theta**Adsorption%Adsorbexp(SurfSide,iSpec)) * exp(-E_des/WallTemp)
    IF (Theta.GT.0) THEN
      Adsorption%ProbDes(p,q,SurfSide,iSpec) = rate * dt /Theta
    ELSE
      Adsorption%ProbDes(p,q,SurfSide,iSpec) = 0.0
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + Adsorption%ProbDes(p,q,SurfSide,iSpec)
    END IF
#endif
!===================================================================================================================================
  ELSE IF (DSMC%WallModel.EQ.2) THEN
!===================================================================================================================================
! transition state theory(TST) and unity bond index-quadratic exponential potential (UBI-QEP) --> lattice gas model 
    Theta = Adsorption%Coverage(p,q,SurfSide,iSpec)  
    IF (Theta.GT.0) THEN
      Q_0A = Adsorption%HeatOfAdsZero(iSpec)
      D_A  = 59870.
      
      rate = 0.
!       iProbSigma = 1
      DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
        SELECT CASE(Adsorption%Coordination(iSpec))
        CASE(1)
          n = 4
          sigma = Adsorption%Sigma(p,q,SurfSide,iSpec,iProbSigma)
          Heat = (Q_0A**2)/(Q_0A/n+D_A) * sigma
        CASE(2)
          n = 2
!         m = Theta * n
!         IF (m.LT.1) THEN
!           Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A)
! !           Heat = (Q_0A**2)/(Q_0A/n+D_A)
!         ELSE 
!           sigma = 1/m*(2-1/m)
          sigma = Adsorption%Sigma(p,q,SurfSide,iSpec,iProbSigma)
          Heat = (9*Q_0A**2)/(6*Q_0A+16*D_A) * sigma
  !         Heat = (Q_0A**2)/(Q_0A/n+D_A) * sigma
!         END IF
        CASE DEFAULT
        END SELECT
        E_des = Heat
!         rate = exp(-E_des/WallTemp)
!         rate = rate + Adsorption%ProbSigma(p,q,SurfSide,iSpec,iProbSigma) * exp(-E_des/WallTemp)
!       END DO
      
        CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSide))
        CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
        nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
        rate = nu_des * exp(-E_des/WallTemp)
        
        Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) = rate * dt
        IF (Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma).GT.1.) THEN
          Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) = 1.
        END IF
      END DO
    ELSE
      Adsorption%ProbSigDes(p,q,SurfSide,iSpec,(1+(iSpec-1)*36):(36*iSpec)) = 0.0
    END IF
#if (PP_TimeDiscMethod==42)
    IF (.NOT.DSMC%ReservoirRateStatistic) THEN
      DO iProbSigma = (1+(iSpec-1)*36),(36*iSpec)
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = Adsorption%AdsorpInfo(iSpec)%MeanProbDes &
                                                  + Adsorption%ProbSigDes(p,q,SurfSide,iSpec,iProbSigma) &
                                                  * Adsorption%ProbSigma(p,q,SurfSide,iSpec,iProbSigma)
      END DO
    END IF
#endif
!===================================================================================================================================
  END IF ! DSMC%WallModel
! (QCA-quasi-chemical-approximation not considered)
  
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
USE MOD_DSMC_Vars,          ONLY : DSMC, SpecDSMC, PolyatomMolDSMC
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
REAL                          :: Qtra, Qrot, Qvib, Qelec
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
USE MOD_DSMC_Vars,          ONLY : DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY : BoltzmannConst, Species
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

END MODULE MOD_DSMC_SurfModel_Tools