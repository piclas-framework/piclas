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
      Adsorption%AdsorpInfo(iSpec)%Accomodation = 0
      IF (KeepWallParticles) THEN
        Adsorption%AdsorpInfo(iSpec)%NumOfDes = 0
        Adsorption%AdsorpInfo(iSpec)%MeanProbDes = 0.
      END IF
!       Adsorption%AdsorpInfo(iSpec)%MeanEAds = 0.
    END DO
#endif
    
    IF (.NOT.KeepWallParticles) THEN
      ! adjust coverages of all species on surfaces
      DO iSpec = 1,nSpecies
        DO iSurfSide = 1,SurfMesh%nSides
          DO q = 1,nSurfSample
            DO p = 1,nSurfSample
              maxPart = Adsorption%DensSurfAtoms(iSurfSide) * SurfMesh%SurfaceArea(p,q,iSurfSide)
              Adsorption%Coverage(p,q,iSurfSide,iSpec) = Adsorption%Coverage(p,q,iSurfSide,iSpec) &
                  - (Adsorption%SumDesorbPart(p,q,iSurfSide,iSpec) - Adsorption%SumReactPart(p,q,iSurfSide,iSpec)) &
                  * Species(iSpec)%MacroParticleFactor / maxPart
            END DO
          END DO
        END DO
      END DO
      Adsorption%SumDesorbPart(:,:,:,:) = 0
      Adsorption%SumAdsorbPart(:,:,:,:) = 0
      Adsorption%SumReactPart(:,:,:,:) = 0
    END IF

    IF (DSMC%WallModel.EQ.2) THEN
      CALL CalcDistNumChange()
      CALL CalcDiffusion()
      CALL CalcSurfDistInteraction()
    END IF
    
    IF (DSMC%WallModel.EQ.1 .AND. DSMC%WallModel.EQ.2) THEN
      CALL CalcAdsorbProb()
      IF (KeepWallParticles) CALL CalcDesorbprob()
    END IF
  END IF

END SUBROUTINE DSMC_Update_Wall_Vars  
    
SUBROUTINE Particle_Wall_Adsorb(PartTrajectory,alpha,xi,eta,PartID,GlobSideID,IsSpeciesSwap,adsindex,BCSideID) 
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
  USE MOD_Particle_Surfaces_vars, ONLY : SideNormVec,SideType
  USE MOD_Particle_Surfaces,      ONLY : CalcNormAndTangBilinear,CalcNormAndTangBezier
!===================================================================================================================================
   IMPLICIT NONE
!===================================================================================================================================
! argument list declaration
  INTEGER,INTENT(INOUT)            :: adsindex
  REAL,INTENT(INOUT)               :: PartTrajectory(1:3), alpha
  REAL,INTENT(IN)                  :: xi, eta
  INTEGER,INTENT(IN)               :: PartID, GlobSideID
  LOGICAL,INTENT(IN)               :: IsSpeciesSwap
  INTEGER,INTENT(IN),OPTIONAL      :: BCSideID
! LOCAL VARIABLES
  INTEGER                          :: locBCID, VibQuant, VibQuantNew, VibQuantWall
  REAL                             :: VibQuantNewR
  REAL                             :: VeloReal, RanNum, EtraOld
  REAL                             :: EtraWall, EtraNew
  REAL                             :: WallVelo(1:3), WallTemp, TransACC, VibACC, RotACC
  REAL                             :: ErotNew, ErotWall, EVibNew
  REAL                             :: Xitild,EtaTild
  INTEGER                          :: p,q
  REAL                             :: n_loc(1:3), tang1(1:3),tang2(1:3)
  REAL                             :: Adsorption_Prob
  INTEGER                          :: adsorption_case
! Polyatomic Molecules
  REAL, ALLOCATABLE                :: RanNumPoly(:), VibQuantNewRPoly(:)
  INTEGER                          :: iPolyatMole, iDOF
  INTEGER, ALLOCATABLE             :: VibQuantNewPoly(:), VibQuantWallPoly(:), VibQuantTemp(:)
! Local variable declaration            
  REAL                             :: maxPart
  INTEGER                          :: SurfSide, iSpec
  REAL                             :: VelXold, VelYold, VelZold
  REAL, PARAMETER                  :: PI=3.14159265358979323846
  REAL                             :: IntersectionPos(1:3)
  REAL                             :: TransArray(1:6),IntArray(1:6)
  REAL                             :: Norm_velo, Norm_Ec
  INTEGER                          :: outSpec(2)
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
  
  IF (DSMC%WallModel.EQ.3) THEN
    IF(PRESENT(BCSideID))THEN
      SELECT CASE(SideType(BCSideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT)
        n_loc=SideNormVec(1:3,BCSideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,BCSideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,BCSideID)
      END SELECT 
    ELSE
      SELECT CASE(SideType(GlobSideID))
      CASE(PLANAR_RECT,PLANAR_NONRECT)
        n_loc=SideNormVec(1:3,GlobSideID)
      CASE(BILINEAR)
        CALL CalcNormAndTangBilinear(n_loc,tang1,tang2,xi,eta,GlobSideID)
      CASE(CURVED)
        CALL CalcNormAndTangBezier(n_loc,tang1,tang2,xi,eta,GlobSideID)
      END SELECT 
    END IF
    
    Norm_velo = PartState(PartID,4)*n_loc(1) + PartState(PartID,4)*n_loc(2) + PartState(PartID,6)*n_loc(3)
    Norm_Ec = 0.5 * Species(iSpec)%MassIC * Norm_velo**2 + PartStateIntEn(PartID,1) + PartStateIntEn(PartID,2)
    
    CALL CalcBackgndPartAdsorb(p,q,SurfSide,PartID,Norm_Ec,adsorption_case,outSpec)
  ELSE
    adsorption_case = 0
    Adsorption_prob = Adsorption%ProbAds(p,q,SurfSide,iSpec)
    CALL RANDOM_NUMBER(RanNum)
    IF ( (Adsorption_prob.GE.RanNum) .AND. & 
       (Adsorption%Coverage(p,q,SurfSide,iSpec).LT.Adsorption%MaxCoverage(SurfSide,iSpec)) ) THEN
      adsorption_case = 1
    END IF
  END IF
  
  SELECT CASE(adsorption_case)
  !-----------------------------------------------------------------------------------------------------------------------------
  CASE(-1) ! perfect elastic scattering
  !-----------------------------------------------------------------------------------------------------------------------------
    adsindex = -1
  !-----------------------------------------------------------------------------------------------------------------------------
  CASE(1) ! molecular adsorption
  !-----------------------------------------------------------------------------------------------------------------------------    
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(p,q,SurfSide)
    Adsorption%Coverage(p,q,SurfSide,outSpec(1)) = Adsorption%Coverage(p,q,SurfSide,outSpec(1)) & 
                                                 + Species(outSpec(1))%MacroParticleFactor/maxPart
    adsindex = 1
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(outSpec(1))%NumOfAds = Adsorption%AdsorpInfo(outSpec(1))%NumOfAds + 1
#endif
!     ! allocate particle belonging adsorbing side-index and side-subsurface-indexes
!     IF (KeepWallParticles) THEN
!       PDM%PartAdsorbSideIndx(1,PartID) = GlobSideID
!       PDM%PartAdsorbSideIndx(2,PartID) = p
!       PDM%PartAdsorbSideIndx(3,PartID) = q
!     END IF
!   
!     LastPartPos(PartID,1) = PartState(PartID,1)
!     LastPartPos(PartID,2) = PartState(PartID,2)
!     LastPartPos(PartID,3) = PartState(PartID,3)
!     VelXold = PartState(PartID,4)
!     VelYold = PartState(PartID,5)
!     VelZold = PartState(PartID,6)
!     ! intersection point with surface
!     IntersectionPos(1:3) = LastPartPos(PartID,1:3) + PartTrajectory(1:3)*alpha
!     PartState(PartID,1)  = IntersectionPos(1)
!     PartState(PartID,2)  = IntersectionPos(2)
!     PartState(PartID,3)  = IntersectionPos(3)
!     PartState(PartID,4)  = WallVelo(1)
!     PartState(PartID,5)  = WallVelo(2)
!     PartState(PartID,6)  = WallVelo(3)
!   
!     VeloReal = SQRT(VelXold * VelXold + VelYold * VelYold + VelZold * VelZold)
!     EtraOld = 0.5 * Species(iSpec)%MassIC * VeloReal**2
!     EtraWall = 0.0
!     EtraNew = 0.0
!     
!     TransArray(1) = EtraOld
!     TransArray(2) = EtraWall
!     TransArray(3) = EtraNew
!     TransArray(4) = PartState(PartID,4)-VelXold
!     TransArray(5) = PartState(PartID,5)-VelYold
!     TransArray(6) = PartState(PartID,6)-VelZold
!   
!     !---- Internal energy accommodation
!     IF (CollisMode.GT.1) THEN
!     IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
!       !---- Rotational energy accommodation
!       CALL RANDOM_NUMBER(RanNum)
!       ErotWall = - BoltzmannConst * WallTemp * LOG(RanNum)
!       ErotNew  = PartStateIntEn(PartID,2) + RotACC *(ErotWall - PartStateIntEn(PartID,2))
!       IntArray(1) = PartStateIntEn(PartID,2)
!       IntArray(2) = ErotWall
!       IntArray(3) = ErotNew
!       PartStateIntEn(PartID,2) = ErotNew
!       !---- Vibrational energy accommodation
!       IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
!         EvibNew = 0.0
!         iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
!         ALLOCATE(RanNumPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF),VibQuantWallPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                  VibQuantNewRPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), VibQuantNewPoly(PolyatomMolDSMC(iPolyatMole)%VibDOF), &
!                  VibQuantTemp(PolyatomMolDSMC(iPolyatMole)%VibDOF))
!         CALL RANDOM_NUMBER(RanNumPoly)
!         VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
!         DO WHILE (ALL(VibQuantWallPoly.GE.PolyatomMolDSMC(iPolyatMole)%MaxVibQuantDOF))
!           CALL RANDOM_NUMBER(RanNumPoly)
!           VibQuantWallPoly(:) = INT(-LOG(RanNumPoly(:)) * WallTemp / PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(:))
!         END DO
!         VibQuantNewRPoly(:) = VibQuantsPar(PartID)%Quants(:) + VibACC*(VibQuantWallPoly(:) - VibQuantsPar(PartID)%Quants(:))
!         VibQuantNewPoly = INT(VibQuantNewRPoly)
!         DO iDOF = 1, PolyatomMolDSMC(iPolyatMole)%VibDOF
!           CALL RANDOM_NUMBER(RanNum)
!           IF (RanNum.LT.(VibQuantNewRPoly(iDOF) - VibQuantNewPoly(iDOF))) THEN
!             EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant + 1.0d0) &
!                       * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!             VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF) + 1
!           ELSE
!             EvibNew = EvibNew + (VibQuantNewPoly(iDOF) + DSMC%GammaQuant) &
!                       * BoltzmannConst*PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF)
!             VibQuantTemp(iDOF) = VibQuantNewPoly(iDOF)
!           END IF
!           IntArray(4) = IntArray(4) + (VibQuantsPar(PartID)%Quants(iDOF) + DSMC%GammaQuant) * BoltzmannConst &
!                       * PolyatomMolDSMC(iPolyatMole)%CharaTVibDOF(iDOF) * Species(iSpec)%MacroParticleFactor
!           IntArray(5) = IntArray(5) + VibQuantWallPoly(iDOF) * BoltzmannConst &
!                       * SpecDSMC(iSpec)%CharaTVib * Species(iSpec)%MacroParticleFactor
!         END DO
!       ELSE
!         VibQuant     = NINT(PartStateIntEn(PartID,1)/(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) &
!                     - DSMC%GammaQuant)
!         CALL RANDOM_NUMBER(RanNum)
!         VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
!         DO WHILE (VibQuantWall.GE.SpecDSMC(iSpec)%MaxVibQuant)
!           CALL RANDOM_NUMBER(RanNum)
!           VibQuantWall = INT(-LOG(RanNum) * WallTemp / SpecDSMC(iSpec)%CharaTVib)
!         END DO
!         VibQuantNewR = VibQuant + VibACC*(VibQuantWall - VibQuant)
!         VibQuantNew = INT(VibQuantNewR)
!         CALL RANDOM_NUMBER(RanNum)
!         IF (RanNum.LT.(VibQuantNewR - VibQuantNew)) THEN
!           EvibNew = (VibQuantNew + DSMC%GammaQuant + 1.0d0)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
!         ELSE
!           EvibNew = (VibQuantNew + DSMC%GammaQuant)*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
!         END IF
!         IntArray(4) = VibQuant + DSMC%GammaQuant * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
!         IntArray(5) = VibQuantWall * BoltzmannConst * SpecDSMC(iSpec)%CharaTVib
!       END IF
!       IntArray(6) = EvibNew
!     END IF
!     END IF
!     !End internal energy accomodation
!     
! !----  Sampling at walls
!     IF ((DSMC%CalcSurfaceVal.AND.(Time.ge.(1-DSMC%TimeFracSamp)*TEnd)).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroValues)) THEN
!       CALL CalcWallSample(PartID,SurfSide,p,q,Transarray,IntArray,PartTrajectory,alpha,IsSpeciesSwap)
!     END IF
    
  !-----------------------------------------------------------------------------------------------------------------------------
  CASE(2) ! dissociative adsorption
  !-----------------------------------------------------------------------------------------------------------------------------
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(p,q,SurfSide)
    Adsorption%Coverage(p,q,SurfSide,outSpec(1)) = Adsorption%Coverage(p,q,SurfSide,outSpec(1)) & 
                                                 + Species(outSpec(1))%MacroParticleFactor/maxPart
    Adsorption%Coverage(p,q,SurfSide,outSpec(2)) = Adsorption%Coverage(p,q,SurfSide,outSpec(2)) & 
                                                 + Species(outSpec(2))%MacroParticleFactor/maxPart
    adsindex = 1
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(outSpec(1))%NumOfAds = Adsorption%AdsorpInfo(outSpec(1))%NumOfAds + 1
    Adsorption%AdsorpInfo(outSpec(2))%NumOfAds = Adsorption%AdsorpInfo(outSpec(2))%NumOfAds + 1
#endif
  !-----------------------------------------------------------------------------------------------------------------------------
  CASE(3) ! Eley-Rideal reaction
  !-----------------------------------------------------------------------------------------------------------------------------
    maxPart = Adsorption%DensSurfAtoms(SurfSide) * SurfMesh%SurfaceArea(p,q,SurfSide)
    Adsorption%Coverage(p,q,SurfSide,outSpec(1)) = Adsorption%Coverage(p,q,SurfSide,outSpec(1)) & 
                                                 - Species(outSpec(1))%MacroParticleFactor/maxPart
    adsindex = 2
#if (PP_TimeDiscMethod==42)
    Adsorption%AdsorpInfo(outSpec(1))%NumOfDes = Adsorption%AdsorpInfo(outSpec(1))%NumOfDes + 1
#endif
  !-----------------------------------------------------------------------------------------------------------------------------
  CASE DEFAULT ! diffuse reflection
  !-----------------------------------------------------------------------------------------------------------------------------
    adsindex = 0
  END SELECT
  
END SUBROUTINE Particle_Wall_Adsorb

SUBROUTINE Calc_PartNum_Wall_Desorb()
!===================================================================================================================================
! calculation of desorbing particle number when particles deleted at adsorption and inserted at desorption
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
   REAL                             :: Heat_i, Heat_j, Heat_temp, sigma
   INTEGER                          :: bondorder, n_equal_site_Neigh, Indx, Indy, Surfpos, newpos
   INTEGER , ALLOCATABLE            :: free_Neigh_pos(:)
!===================================================================================================================================

DO SurfSideID = 1,SurfMesh%nSides
DO subsurfeta = 1,nSurfSample
DO subsurfxi = 1,nSurfSample

  DO Coord = 1,3
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    
    ALLOCATE ( free_Neigh_pos(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    DO AdsorbID = nSitesRemain+1,nSites,1
      Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
      iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
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
      Heat_i = Adsorption%HeatOfAdsZero(iSpec) * sigma
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
      
      ! choose Neighbour position with highest heat of adsorption if adsorbate would move there
      Heat_j = 0.
      DO i = 1,n_equal_site_Neigh
!       IF (n_equal_site_Neigh .GE. 1) THEN
!         CALL RANDOM_NUMBER(RanNum)
!         i = 1 + INT(n_equal_site_Neigh*RanNum)
        sigma = 0.
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(free_Neigh_pos(i),j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(free_Neigh_pos(i),j) 
          bondorder = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
          sigma = sigma + ( (1./SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom) &
                * (1./(bondorder)) * (2.-(1./bondorder)) )
        END DO
        Heat_temp = Adsorption%HeatOfAdsZero(iSpec) * sigma
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

SUBROUTINE CalcBackgndPartAdsorb(subsurfxi,subsurfeta,SurfSideID,PartID,Norm_Ec,adsorption_case,outSpec)
!===================================================================================================================================
! Particle Adsorption probability calculation for wallmodel 3
!===================================================================================================================================
  USE MOD_Globals_Vars,           ONLY : PlanckConst
  USE MOD_Particle_Vars,          ONLY : PartState, PartSpecies, nSpecies, Species, BoltzmannConst
  USE MOD_Mesh_Vars,              ONLY : BC
  USE MOD_DSMC_Vars,              ONLY : Adsorption, DSMC, SurfDistInfo, SpecDSMC, PartStateIntEn, PolyatomMolDSMC
  USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
  USE MOD_TimeDisc_Vars,          ONLY : dt
  USE MOD_DSMC_Analyze,           ONLY : CalcTVib, CalcTVibPoly
#if (PP_TimeDiscMethod==42)  
  USE MOD_TimeDisc_Vars,          ONLY : iter
#endif
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! argument list declaration
  INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,PartID
  REAL,INTENT(IN)                  :: Norm_Ec
  INTEGER,INTENT(OUT)              :: adsorption_case
  INTEGER,INTENT(OUT)              :: outSpec(2)
! LOCAL VARIABLES
  INTEGER                          :: iSpec, globSide, iPolyatMole
  INTEGER                          :: Coord, Coord2, i, j, AdsorbID, numSites, IDRearrange
  INTEGER                          :: new_adsorbates(2), nSites, nSitesRemain
  REAL                             :: difference, maxPart, new_coverage
  REAL                             :: WallTemp, RanNum
  REAL                             :: Prob_ads
  REAL , ALLOCATABLE               :: P_Eley_Rideal(:), Prob_diss(:)
  INTEGER                          :: Surfpos, ReactNum
  INTEGER , ALLOCATABLE            :: reactadsorbnum(:), adsorbnum(:)
  REAL, PARAMETER                  :: Pi=3.14159265358979323846_8
  INTEGER                          :: jSpec, kSpec, jCoord, kCoord
  REAL                             :: sum_probabilities
  INTEGER , ALLOCATABLE            :: NeighbourID(:,:)
  INTEGER                          :: SiteSpec, Neighpos_j, Neighpos_k, chosen_Neigh_j, chosen_Neigh_k
  INTEGER                          :: n_empty_Neigh(3), n_react_Neigh(3), n_Neigh(3), adsorbates(nSpecies)
  REAL                             :: E_a, c_f, EZeroPoint_Educt, E_col, phi_1, phi_2, Xi_Total, Xi_vib
  REAL                             :: Heat_A, Heat_B, Heat_AB, D_AB, D_A, D_B
  
  REAL                             :: vel_norm, vel_coll, potential_pot, a_const, mu, surfmass, trapping_prob
!===================================================================================================================================
! special TPD (temperature programmed desorption) surface temperature adjustment part
  globSide = Adsorption%SurfSideToGlobSideMap(SurfSideID)
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
  ! calculate number of adsorbates for each species (already on surface)
  adsorbates(:) = 0
  DO Coord = 1,3
  nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
  nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
  DO AdsorbID = 1,nSites-nSitesRemain
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain+AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(iSpec) = adsorbates(iSpec) + 1
  END DO
  END DO
! initialize variables
  Prob_ads = 0.
  ALLOCATE( P_Eley_Rideal(1:Adsorption%ReactNum),&
            Prob_diss(1:Adsorption%ReactNum))
  P_Eley_Rideal(:) = 0.
  Prob_diss(:) = 0.
  
  ! Choose Random surface site with species coordination
  CALL RANDOM_NUMBER(RanNum)
  iSpec = PartSpecies(PartID)
  Coord = Adsorption%Coordination(iSpec)
  AdsorbID = 1 + INT(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)*RanNum)
  Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
  SiteSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate trapping probability
  !---------------------------------------------------------------------------------------------------------------------------------
  surfmass = 1.652E-25
  potential_pot = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.TRUE.)*BoltzmannConst
  vel_norm = - ( 2*(Norm_Ec+potential_pot)/Species(iSpec)%MassIC)**(0.5)
  a_const = (Species(iSpec)%MassIC/(2*BoltzmannConst*WallTemp))**(0.5)
  IF (SiteSpec.EQ.0) THEN
    mu = Species(iSpec)%MassIC / surfmass
  ELSE
    mu = Species(iSpec)%MassIC / (Species(SiteSpec)%MassIC)
  END IF
  vel_coll = 0.5 * ( (1+mu)*(2*(potential_pot/Species(iSpec)%MassIC))**(0.5) + (1-mu)*vel_norm )
  trapping_prob = 0.5 + 0.5*erf(a_const*vel_coll) + ( exp(-(a_const**2)*(vel_coll**2)) / (2*a_const*(vel_norm*PI**0.5)) )
#if (PP_TimeDiscMethod==42)
  Adsorption%AdsorpInfo(iSpec)%Accomodation = Adsorption%AdsorpInfo(iSpec)%Accomodation + trapping_prob
#endif
  CALL RANDOM_NUMBER(RanNum)
  IF (RanNum.GT.trapping_prob) THEN
    outSpec(1) = iSpec
    outSpec(2) = 0
    adsorption_case = -1
    RETURN
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for molecular adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (SiteSpec.EQ.0) THEN
    ! calculation of molecular adsorption probability with TCE
    EZeroPoint_Educt = 0.
    E_a = 0.
    c_f = Adsorption%DensSurfAtoms(SurfSideID) / ( (BoltzmannConst / (2*Pi*Species(iSpec)%MassIC))**0.5 )
    ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
    IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
      IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
        iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
        EZeroPoint_Educt = EZeroPoint_Educt + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
        ! Calculation of the vibrational degree of freedom for the particle 
        IF (PartStateIntEn(PartID,1).GT.PolyatomMolDSMC(iPolyatMole)%EZeroPoint) THEN
          Xi_vib = 2.*(PartStateIntEn(PartID,1)-PolyatomMolDSMC(iPolyatMole)%EZeroPoint) &
                  / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(PartID,1), iSpec))
        ELSE
          Xi_vib = 0.0
        END IF
      ELSE
        EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
        IF((PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).GT.0.0) THEN
!           IF(ChemReac%MeanEVibQua_PerIter(iSpec).GT.0.0) THEN
          Xi_vib = 2.*(PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) &
                  / (BoltzmannConst*CalcTVib(SpecDSMC(iSpec)%CharaTVib, PartStateIntEn(PartID,1), SpecDSMC(iSpec)%MaxVibQuant))
!             Xi_vib = 2.0*ChemReac%MeanEVibQua_PerIter(iSpec) &
!                     * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0)
        ELSE
          Xi_vib = 0.0
        END IF
      END IF
    ELSE
      Xi_vib = 0.0
    END IF
    Xi_Total = Xi_vib + 2 + 3 !Xi_rot + 3
    phi_1 = Adsorption%Ads_Powerfactor(iSpec) - 3./2. + Xi_Total
    phi_2 = 1 - Xi_Total
    Prob_ads = Calc_Beta_Adsorb(iSpec,Xi_Total,c_f) * (Norm_Ec - E_a)**phi_1 * Norm_Ec ** phi_2
  ELSE
    Prob_ads = 0.
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! sort Neighbours to coordinations for search of two random positions for dissociative adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  n_Neigh(:) = 0
  n_empty_Neigh(:) = 0
  ALLOCATE(NeighbourID(1:3,1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))

  DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
  DO Coord2 = 1,3
    IF (Coord2.EQ.SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)) THEN
      IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( & 
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) & 
            .EQ.0) ) THEN
        n_empty_Neigh(Coord2) = n_empty_Neigh(Coord2) + 1
      END IF
      n_Neigh(Coord2) = n_Neigh(Coord2) + 1
      NeighbourID(Coord2,n_Neigh(Coord2)) = i
    END IF
  END DO
  END DO
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for dissociative adsorption
  !---------------------------------------------------------------------------------------------------------------------------------
  IF (SpecDSMC(iSpec)%InterID.EQ.2) THEN
  DO ReactNum = 1,(Adsorption%DissNum)
    jSpec = Adsorption%DissocReact(1,ReactNum,iSpec)
    kSpec = Adsorption%DissocReact(2,ReactNum,iSpec)
    IF ((iSpec.NE.0) .AND. (kSpec.NE.0)) THEN
    jCoord = Adsorption%Coordination(jSpec)
    kCoord = Adsorption%Coordination(kSpec)
    Neighpos_j = 0
    Neighpos_k = 0
    IF ( (jCoord.EQ.Coord) .OR. (kCoord.EQ.Coord) ) THEN
      IF (jCoord.EQ.Coord) THEN
        Neighpos_j = Surfpos
        IF (n_empty_Neigh(kCoord).GT.0) THEN
          ! assign availiable neighbour position for k
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
          Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        END IF
      ELSE
        Neighpos_k = Surfpos
        IF (n_empty_Neigh(jCoord).GT.0) THEN
          ! assign availiable neighbour position for j
          CALL RANDOM_NUMBER(RanNum)
          chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
          Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
          IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
          NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
          NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        END IF
      END IF
    ELSE
      IF ( (jCoord.EQ.kCoord) .AND. (n_empty_Neigh(jCoord).GT.1) ) THEN
        ! assign availiable neighbour positions
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
        NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord)-1)*RanNum)
        Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord)-1)
        NeighbourID(jCoord,n_Neigh(jCoord)-1) = IDRearrange
      ELSE IF ( (jCoord.NE.kCoord) .AND. (n_empty_Neigh(jCoord).GT.0) .AND. (n_empty_Neigh(kCoord).GT.0) ) THEN
        ! assign availiable neighbour positions
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(REAL(n_Neigh(jCoord))*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_j) = NeighbourID(jCoord,n_Neigh(jCoord))
        NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_k = 1 + INT(REAL(n_Neigh(kCoord))*RanNum)
        Neighpos_k = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(kCoord,chosen_Neigh_k))
        IDRearrange = NeighbourID(jCoord,chosen_Neigh_j)
        NeighbourID(jCoord,chosen_Neigh_k) = NeighbourID(kCoord,n_Neigh(kCoord))
        NeighbourID(jCoord,n_Neigh(jCoord)) = IDRearrange
      END IF
    END IF
    IF ( (Neighpos_j.GT.0) .AND. (Neighpos_k.GT.0) ) THEN
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) .AND. &
         (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.0) ) THEN
      ! calculation of activation energy
      Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,Neighpos_j,.TRUE.)
      Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,kSpec,Neighpos_k,.TRUE.)
      Heat_AB = 0. ! direct dissociative adsorption
      D_AB = Adsorption%EDissBond((ReactNum),iSpec)
      D_A = 0. !Adsorption%EDissBond((ReactNum),jSpec)
      D_B = 0. !Adsorption%EDissBond((ReactNum),kSpec)
      E_a = Calc_E_act(Heat_A,Heat_B,Heat_AB,D_AB,D_A,D_B,.TRUE.)
      ! calculation of dissociative adsorption probability with TCE
      EZeroPoint_Educt = 0.
      c_f = Adsorption%DensSurfAtoms(SurfSideID) / ( (BoltzmannConst / (2*Pi*Species(iSpec)%MassIC))**0.5 )
      ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
      IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
        IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
          iPolyatMole = SpecDSMC(iSpec)%SpecToPolyArray
          EZeroPoint_Educt = EZeroPoint_Educt + PolyatomMolDSMC(iPolyatMole)%EZeroPoint
          ! Calculation of the vibrational degree of freedom for the particle 
          IF (PartStateIntEn(PartID,1).GT.PolyatomMolDSMC(iPolyatMole)%EZeroPoint) THEN
            Xi_vib = 2.*(PartStateIntEn(PartID,1)-PolyatomMolDSMC(iPolyatMole)%EZeroPoint) &
                    / (BoltzmannConst*CalcTVibPoly(PartStateIntEn(PartID,1), iSpec))
          ELSE
            Xi_vib = 0.0
          END IF
        ELSE
          EZeroPoint_Educt = EZeroPoint_Educt + DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib
          IF((PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib).GT.0.0) THEN
!           IF(ChemReac%MeanEVibQua_PerIter(iSpec).GT.0.0) THEN
            Xi_vib = 2.*(PartStateIntEn(PartID,1)-DSMC%GammaQuant*BoltzmannConst*SpecDSMC(iSpec)%CharaTVib) &
                    / (BoltzmannConst*CalcTVib(SpecDSMC(iSpec)%CharaTVib, PartStateIntEn(PartID,1), SpecDSMC(iSpec)%MaxVibQuant))
!             Xi_vib = 2.0*ChemReac%MeanEVibQua_PerIter(iSpec) &
!                     * LOG(1.0/ChemReac%MeanEVibQua_PerIter(iSpec) + 1.0)
          ELSE
            Xi_vib = 0.0
          END IF
        END IF
      ELSE
        Xi_vib = 0.0
      END IF
      IF ((Norm_Ec-EZeroPoint_Educt).GT.E_a) THEN
        Xi_Total = Xi_vib + 2 + 3 !Xi_rot + 3
        phi_1 = Adsorption%Ads_Powerfactor(iSpec) - 3./2. + Xi_Total
        phi_2 = 1 - Xi_Total
        Prob_diss(ReactNum) = Calc_Beta_Adsorb(iSpec,Xi_Total,c_f) * (Norm_Ec - E_a)**phi_1 * Norm_Ec ** phi_2
      END IF
    END IF
    END IF
    END IF
  END DO
  END IF
  !---------------------------------------------------------------------------------------------------------------------------------
  ! calculate probability for Eley-Rideal reaction
  !---------------------------------------------------------------------------------------------------------------------------------
  DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
    ! reaction partner
    jSpec = Adsorption%AssocReact(1,ReactNum,iSpec)
    Neighpos_j = 0
    IF (jSpec.EQ.0) CYCLE
    ! reaction results
    kSpec = Adsorption%AssocReact(2,ReactNum,iSpec)
    jCoord = Adsorption%Coordination(jSpec)
    IF (jCoord.EQ.Coord) THEN
      Neighpos_j = Surfpos
    ELSE
      IF ( n_Neigh(jCoord)-n_empty_Neigh(jCoord).GT.0 ) THEN
        CALL RANDOM_NUMBER(RanNum)
        chosen_Neigh_j = 1 + INT(n_Neigh(jCoord)*RanNum)
        Neighpos_j = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(1)%NeighPos(Surfpos,NeighbourID(jCoord,chosen_Neigh_j))
      END IF
    END IF
    IF ( (Neighpos_j.GT.0) ) THEN
    IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(jCoord)%Species(Neighpos_j).EQ.jSpec) ) THEN
      ! calculation of activation energy
      ! calculation of reaction probability with TCE
    END IF
    END IF
  END DO
  !---------------------------------------------------------------------------------------------------------------------------------
  
  sum_probabilities = Prob_ads
  DO ReactNum = 1,Adsorption%ReactNum
    sum_probabilities = sum_probabilities + Prob_diss(ReactNum) + P_Eley_Rideal(ReactNum)
  END DO
  ! choose adsorption case
  adsorption_case = 0
  CALL RANDOM_NUMBER(RanNum)
  IF (sum_probabilities .GT. RanNum) THEN
    ! chose surface reaction case (0=scattering, 1=adsorption, 2=reaction (dissociation), 3=reaction (Eley-Rideal))
    DO ReactNum = 1,(Adsorption%ReactNum)
      CALL RANDOM_NUMBER(RanNum)
      IF ((P_Eley_Rideal(ReactNum)/sum_probabilities).GT.RanNum) THEN
        ! if ER-reaction set output parameters
        adsorption_case = 3
        outSpec(1) = Adsorption%AssocReact(1,ReactNum,iSpec)
        outSpec(2) = Adsorption%AssocReact(2,ReactNum,iSpec)
        ! adjust number of background mapping adsorbates
        maxPart = Adsorption%DensSurfAtoms(SurfSideID) * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID)
        IF ( Adsorption%Coordination(outSpec(1)).EQ.2 ) maxPart = maxPart*2
        new_coverage = Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,outSpec(1)) &
                      - Species(outSpec(1))%MacroParticleFactor / maxPart
        numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(outSpec(1)))
        new_adsorbates(:) = 0
        IF ( (REAL(adsorbates(OutSpec(1)))/REAL(numSites)).GT. new_coverage) THEN
          difference = REAL(adsorbates(OutSpec(1)))/REAL(numSites) - new_coverage
          new_adsorbates(1) = -NINT(difference*REAL(numSites))
        END IF
        IF ( (new_adsorbates(1).LT.0) ) CALL AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,new_adsorbates(1),outSpec(1))
        EXIT
      END IF
      sum_probabilities = sum_probabilities - P_Eley_Rideal(ReactNum)
      CALL RANDOM_NUMBER(RanNum)
      IF ((Prob_diss(ReactNum)/sum_probabilities).GT.RanNum) THEN
        ! if dissocciative adsorption set output parameters
        adsorption_case = 2
        outSpec(1) = Adsorption%DissocReact(1,ReactNum,iSpec)
        outSpec(2) = Adsorption%DissocReact(2,ReactNum,iSpec)
        ! adjust number of background mapping adsorbates 
        new_adsorbates(:) = 0
        DO i = 1,2
          maxPart = Adsorption%DensSurfAtoms(SurfSideID) * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID)
          IF ( Adsorption%Coordination(outSpec(i)).EQ.2 ) maxPart = maxPart*2
          new_coverage = Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,outSpec(i)) &
                        + Species(outSpec(i))%MacroParticleFactor / maxPart
          numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(outSpec(i)))
          IF ( new_coverage .GT. (REAL(adsorbates(OutSpec(i)))/REAL(numSites)) ) THEN
            difference = new_coverage - (REAL(adsorbates(OutSpec(i)))/REAL(numSites))
            new_adsorbates(i) = NINT(difference*REAL(numSites))
          END IF
        END DO
        IF ( (new_adsorbates(1).GT.0) ) CALL AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,new_adsorbates(1),outSpec(1))
        IF ( (new_adsorbates(2).GT.0) ) CALL AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,new_adsorbates(2),outSpec(2))
        EXIT
      END IF
      sum_probabilities = sum_probabilities - Prob_diss(ReactNum)
    END DO
    IF (adsorption_case.EQ.0) THEN
      CALL RANDOM_NUMBER(RanNum)
      IF ((Prob_ads/sum_probabilities).GT.RanNum) THEN
        ! if molecular adsorption set output parameters
        adsorption_case = 1
        outSpec(1) = iSpec
        outSpec(2) = 0
        maxPart = Adsorption%DensSurfAtoms(SurfSideID) * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID)
        IF ( Adsorption%Coordination(outSpec(1)).EQ.2 ) maxPart = maxPart*2
        new_coverage = Adsorption%Coverage(subsurfxi,subsurfeta,SurfSideID,outSpec(1)) &
                      + Species(outSpec(1))%MacroParticleFactor / maxPart
        numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(outSpec(1)))
        new_adsorbates(:) = 0
        IF ( new_coverage .GT. (REAL(adsorbates(OutSpec(1)))/REAL(numSites)) ) THEN
          difference = new_coverage - (REAL(adsorbates(OutSpec(1)))/REAL(numSites))
          new_adsorbates(1) = NINT(difference*REAL(numSites))
        END IF
        IF ( (new_adsorbates(1).GT.0) ) CALL AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,new_adsorbates(1),outSpec(1))
      END IF
    END IF
  END IF
#if (PP_TimeDiscMethod==42)
  IF (.NOT.DSMC%ReservoirRateStatistic) THEN
    Adsorption%AdsorpInfo(iSpec)%MeanProbAds = Adsorption%AdsorpInfo(iSpec)%MeanProbDes + Prob_ads
  END IF
#endif

  DEALLOCATE(P_Eley_Rideal,Prob_diss,NeighbourID)

END SUBROUTINE CalcBackgndPartAdsorb

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
   INTEGER                          :: Coord, Coord2, i, j, AdsorbID, numSites
   REAL                             :: WallTemp, RanNum
   REAL                             :: Heat_A, Heat_B, nu_des, rate, P_actual_des
   REAL , ALLOCATABLE               :: P_des(:), Energy(:)
   INTEGER                          :: Indx, Indy, Surfpos
   REAL                             :: VarPartitionFuncAct, VarPartitionFuncWall
   INTEGER                          :: trace, traceNum
   INTEGER , ALLOCATABLE            :: desorbnum(:), reactdesorbnum(:), adsorbnum(:)
   INTEGER , ALLOCATABLE            :: nSites(:), nSitesRemain(:), remainNum(:), adsorbates(:)
!---------- reaction variables
   INTEGER                          :: react_Neigh, n_react_Neigh, n_empty_Neigh, jSpec, kSpec, ReactNum, PartnerID, LastRemainID
   INTEGER                          :: react_Neigh_pos, surf_react_case, interatom
   REAL                             :: E_a, E_d, E_diff, nu_react, P_diff, nu_diff
   REAL                             :: Heat_AB, D_AB, sum_probabilities
   REAL                             :: VarPartitionFuncWall1, VarPartitionFuncWall2
   INTEGER , ALLOCATABLE            :: NeighbourID(:)
   REAL , ALLOCATABLE               :: P_react(:), P_actual_react(:)
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
!   Adsorption%AdsorpInfo(:)%MeanProbDes = 0.
  Adsorption%AdsorpInfo(:)%NumOfDes = 0
#endif

ALLOCATE (desorbnum(1:nSpecies),&
          reactdesorbnum(1:nSpecies),&
          adsorbnum(1:4),&
          nSites(1:4),&
          nSitesRemain(1:4),&
          remainNum(1:4),&
          P_des(1:nSpecies),&
          adsorbates(1:nSpecies),&
          Energy(1:nSpecies))
ALLOCATE( P_react(1:Adsorption%ReactNum),&
          P_actual_react(1:Adsorption%ReactNum))

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
  reactdesorbnum(:) = 0
  nSites(:) = 0
  nSitesRemain(:) = 0
  remainNum(:) = 0
  P_des(:) = 0.
  P_react(:) = 0.
  Energy(:) = 0.
  
  DO Coord=1,3
    nSites(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain(Coord) = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    nSites(4) = nSites(4) + nSites(Coord)
    nSitesRemain(4) = nSitesRemain(4) + nSitesRemain(Coord)
  END DO
  adsorbnum(:) = nSites(:) - nSitesRemain(:)
  traceNum = adsorbnum(4)
  trace = 1
  adsorbates(:) = 0

  ! calculate number of adsorbates for each species (for analyze)
  DO Coord = 1,3
  DO AdsorbID = 1,nSites(Coord)-nSitesRemain(Coord)
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSitesRemain(Coord)+AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    adsorbates(iSpec) = adsorbates(iSpec) + 1
  END DO
  END DO
  
  IF (adsorbnum(4) .GT. 0) THEN
  DO WHILE (trace.LE.traceNum)
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
    ELSE IF ((AdsorbID.GT.nSitesRemain(4)) .AND. (AdsorbID.LE.(adsorbnum(1)+nSitesRemain(4)-remainNum(1))) ) THEN
      AdsorbID = AdsorbID - (nSitesRemain(2) + nSitesRemain(3))
      Coord = 1
    ELSE
      CYCLE
    END IF
    
    Surfpos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID)
    iSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos)
    
    ! calculate heat of adsorption for actual site
    Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)
!     Energy(iSpec) = Energy(iSpec) + Heat_A

    ! find neighbour positions with reaction partner or empty with same coord for diffusion    
    n_empty_Neigh = 0
    n_react_Neigh = 0
    ALLOCATE(NeighbourID(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours))
    DO i = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nNeighbours
      Coord2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,i)
      IF (Coord2.GT.0) THEN
        IF ( (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
              .EQ.0) ) THEN
          IF (Coord2.EQ.Coord) THEN
            n_empty_Neigh = n_empty_Neigh + 1
            NeighbourID(n_empty_Neigh+n_react_Neigh) = i
          END IF
        ELSE
          DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
            IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species( &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,i)) &
                .EQ.Adsorption%AssocReact(1,ReactNum,iSpec)) THEN
              n_react_Neigh = n_react_Neigh + 1
              NeighbourID(n_empty_Neigh+n_react_Neigh) = i
              EXIT
            END IF
          END DO
        END IF
      END IF
    END DO
    
    IF ((n_empty_Neigh+n_react_Neigh).GT.0) THEN
    ! choose random neighbour site from relevant neighbours
    CALL RANDOM_NUMBER(RanNum)
    react_Neigh = 1 + INT((n_empty_Neigh+n_react_Neigh) * RanNum)
    Coord2 = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighSite(Surfpos,NeighbourID(react_Neigh))
    react_Neigh_pos = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%NeighPos(Surfpos,NeighbourID(react_Neigh))
    DEALLOCATE(NeighbourID)
    jSpec = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species(react_Neigh_pos)
    ! calculate probability for 'reaction' with partner site
      IF (jSpec.EQ.0) THEN !diffusion
        ! calculate heat of adsorption for actual site
        Heat_A = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,.FALSE.)
        ! remove considered particle from bondorder
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
        END DO
        !calculate heat of adsorption for diffusion site
        Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,react_Neigh_pos,.TRUE.)
        ! add considered particle to bondorder again
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
        END DO
        ! calculate diffusion probability
        P_diff = exp(-(Heat_A - Heat_B)/WallTemp) / (1+exp(-(Heat_A - Heat_B)/Walltemp))
  !       interatom = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%nInterAtom
  !       E_diff = ( (interatom-2)/(4*interatom-2))*Heat_A
  !       CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID))
  !       CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
  !       nu_diff = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall) * (3.E-8)**2/4.
  !       rate = nu_diff * exp(-E_diff/WallTemp)
  !       IF (rate*dt.GT.1) THEN
  !         P_diff = 1.
  !       ELSE
  !         P_diff = rate * dt
  !       END IF
        P_actual_react(:) = 0.
      ELSE !reaction
        DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
        IF ( jSpec.EQ.Adsorption%AssocReact(1,ReactNum,iSpec) ) THEN
          kSpec = Adsorption%AssocReact(2,ReactNum,iSpec)
          ! calculate heats of adsorption
  !         Heat_AB = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,kSpec,Surfpos,.FALSE.)
          Heat_AB = 0. ! immediate desorption
          Heat_B = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,jSpec,react_Neigh_pos,.FALSE.)
          D_AB = Adsorption%EDissBond((Adsorption%DissNum+ReactNum),iSpec)
          ! calculate LH reaction probability
          E_d = Calc_E_act(Heat_A,Heat_B,Heat_AB,D_AB,0.,0.,.FALSE.)
  !         E_d = ((Heat_A * Heat_B)/(Heat_A + Heat_B))
          CALL PartitionFuncAct(kSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID))
          CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall1)
          CALL PartitionFuncSurf(jSpec, WallTemp, VarPartitionFuncWall2)
          nu_react = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / (VarPartitionFuncWall1 * VarPartitionFuncWall2))
          rate = nu_react * exp(-E_d/WallTemp)
          P_actual_react(ReactNum) = rate * dt
          IF (rate*dt.GT.1) THEN
            P_react(ReactNum) = P_react(ReactNum) + 1
            P_des(kSpec) = P_des(kSpec) + 1
          ELSE
            P_react(ReactNum) = P_react(ReactNum) + P_actual_react(ReactNum)
            P_des(kSpec) = P_des(kSpec) + P_actual_react(ReactNum)
          END IF
          P_diff = 0.
        END IF
        END DO
      END IF
    ELSE
      DEALLOCATE(NeighbourID)
      P_diff = 0.
      P_actual_react(:) = 0.
      jSpec = 0
    END IF
    
    ! calculate molecular desorption probability
    CALL PartitionFuncAct(iSpec, WallTemp, VarPartitionFuncAct, Adsorption%DensSurfAtoms(SurfSideID))
    CALL PartitionFuncSurf(iSpec, WallTemp, VarPartitionFuncWall)
    nu_des = ((BoltzmannConst*WallTemp)/PlanckConst) * (VarPartitionFuncAct / VarPartitionFuncWall)
    E_d = Calc_E_Act(Heat_A,0.,0.,0.,0.,0.,.FALSE.)
    rate = nu_des * exp(-E_d/WallTemp)
    P_actual_des = rate * dt
    IF (rate*dt.GT.1) THEN
      P_des(iSpec) = P_des(iSpec) + 1
    ELSE
      P_des(iSpec) = P_des(iSpec) + P_actual_des
    END IF
    
    ! initialize choosing reaction
    sum_probabilities = P_actual_des
    IF ((jSpec.GT.0) .AND. ((Adsorption%ReactNum-Adsorption%DissNum).GT.0)) THEN
      DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
        sum_probabilities = sum_probabilities + P_actual_react(ReactNum)
      END DO
    ELSE
      sum_probabilities = sum_probabilities + P_diff
    END IF
    
    CALL RANDOM_NUMBER(RanNum)
    IF (sum_probabilities .GT. RanNum) THEN
      ! chose surface reaction case (1=reaction, 2=diffusion, 3=non-associative desorption)
      IF ((jSpec.GT.0) .AND. ((Adsorption%ReactNum-Adsorption%DissNum).GT.0)) THEN
        DO ReactNum = 1,(Adsorption%ReactNum-Adsorption%DissNum)
          CALL RANDOM_NUMBER(RanNum)
          IF ((P_actual_react(ReactNum)/sum_probabilities).GT.RanNum) THEN
            surf_react_case = 1
            EXIT
          ELSE
            surf_react_case = 3
            sum_probabilities = sum_probabilities - P_actual_react(ReactNum)
          END IF
        END DO
      ELSE
        CALL RANDOM_NUMBER(RanNum)
        IF (P_diff/sum_probabilities.GT.RanNum) THEN
          surf_react_case = 2
        ELSE
          surf_react_case = 3
        END IF
      END IF
      SELECT CASE(surf_react_case)
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(1) !reaction
      !-----------------------------------------------------------------------------------------------------------------------------
        ! update number of desorptions
        desorbnum(Adsorption%AssocReact(2,ReactNum,iSpec)) = desorbnum(Adsorption%AssocReact(2,ReactNum,iSpec)) + 1
        reactdesorbnum(Adsorption%AssocReact(2,ReactNum,iSpec)) = reactdesorbnum(Adsorption%AssocReact(2,ReactNum,iSpec)) + 1
        reactdesorbnum(Adsorption%AssocReact(1,ReactNum,iSpec)) = reactdesorbnum(Adsorption%AssocReact(1,ReactNum,iSpec)) - 1
        reactdesorbnum(iSpec) = reactdesorbnum(iSpec) - 1
!         Energy(Adsorption%AssocReact(2,ReactNum,iSpec)) = Energy(Adsorption%AssocReact(2,ReactNum,iSpec)) + Heat_A
        ! remove adsorbate and update map
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
        ! remove adsorbate reaction partner and update map
        DO PartnerID = nSitesRemain(Coord2)+1,nSites(Coord2)
        IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%UsedSiteMap(PartnerID).EQ.react_Neigh_pos) THEN
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%Species(react_Neigh_pos) = 0
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%nInterAtom
            Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%BondAtomIndx(react_Neigh_pos,j)
            Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%BondAtomIndy(react_Neigh_pos,j)
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(jSpec,Indx,Indy) - 1
          END DO
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%UsedSiteMap(PartnerID) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%UsedSiteMap(nSitesRemain(Coord2)+1)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord2)%UsedSiteMap(nSitesRemain(Coord2)+1) = react_Neigh_pos
          nSitesRemain(Coord2) = nSitesRemain(Coord2) + 1
          nSitesRemain(4) = nSitesRemain(4) + 1
          ! additional increment to trace number
          trace = trace + 1
        END IF
        END DO !PartnerID
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(2) !diffusion
      !-----------------------------------------------------------------------------------------------------------------------------
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(Surfpos) = 0
        SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%Species(react_Neigh_pos) = iSpec
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%nInterAtom
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(Surfpos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(Surfpos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) - 1
          Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndx(react_Neigh_pos,j)
          Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%BondAtomIndy(react_Neigh_pos,j)
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) = &
              SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(iSpec,Indx,Indy) + 1
        END DO
        ! move adsorbate to empty site and update map
        DO i = 1,nSitesRemain(Coord)
          IF (SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i).EQ.react_Neigh_pos) THEN
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(i) = Surfpos
!             SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = react_Neigh_pos
            ! move Surfpos to MapID in remainNum segment
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
                SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
            SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) =react_Neigh_pos
            EXIT
          END IF
        END DO
        remainNum(Coord) = remainNum(Coord) + 1
        remainNum(4) = remainNum(4) + 1
        ! additional increment to trace number
        trace = trace + 1
      !-----------------------------------------------------------------------------------------------------------------------------
      CASE(3) !(desorption)
      !-----------------------------------------------------------------------------------------------------------------------------
        desorbnum(iSpec) = desorbnum(iSpec) + 1
        Energy(iSpec) = Energy(iSpec) + Heat_A
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
      END SELECT
        
    ELSE ! adsorbate remains on surface on same site
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(AdsorbID) = &
          SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord))
      SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coord)%UsedSiteMap(nSites(Coord)-remainNum(Coord)) = Surfpos
      remainNum(Coord) = remainNum(Coord) + 1
      remainNum(4) = remainNum(4) + 1
    END IF
    ! update number of adsorbates
    adsorbnum(:) = nSites(:) - nSitesRemain(:)
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain(Coord)
    trace = trace + 1
  END DO
  END IF
  
  DO iSpec = 1,nSpecies
  ! calculate number of desorbed particles for each species
  numSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Adsorption%Coordination(iSpec))
  Adsorption%SumDesorbPart(subsurfxi,subsurfeta,SurfSideID,iSpec) = INT((REAL(desorbnum(iSpec)) / REAL(numSites)) &
                                                                  * (Adsorption%DensSurfAtoms(SurfSideID) &
                      * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID) / Species(iSpec)%MacroParticleFactor))
  Adsorption%SumReactPart(subsurfxi,subsurfeta,SurfSideID,iSpec) = INT((REAL(reactdesorbnum(iSpec)) / REAL(numSites)) &
                                                                  * (Adsorption%DensSurfAtoms(SurfSideID) &
                      * SurfMesh%SurfaceArea(subsurfxi,subsurfeta,SurfSideID) / Species(iSpec)%MacroParticleFactor))
  ! analyze rate data
  IF (adsorbates(iSpec).EQ.0) THEN !(desorbnum(iSpec).EQ.0) THEN !
    Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = 0
!     Energy(iSpec) = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.FALSE.)
  ELSE
    Adsorption%ProbDes(subsurfxi,subsurfeta,SurfSideID,iSpec) = P_des(iSpec) / REAL(adsorbates(iSpec)) !/REAL(desorbnum(iSpec)) !
!     Energy(iSpec) = Energy(iSpec) /REAL(adsorbates(iSpec))
  END IF
  IF (desorbnum(iSpec).EQ.0) THEN
    Energy(iSpec) = Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,-1,.FALSE.)
  ELSE
    Energy(iSpec) = Energy(iSpec) /REAL(desorbnum(iSpec))
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
DEALLOCATE(P_react,P_actual_react)

END SUBROUTINE CalcBackgndPartDesorb

SUBROUTINE AdjustBackgndAdsNum(subsurfxi,subsurfeta,SurfSideID,adsorbates_num,iSpec)
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
  INTEGER,INTENT(IN)               :: subsurfxi,subsurfeta,SurfSideID,adsorbates_num,iSpec
! LOCAL VARIABLES
  INTEGER                          :: dist
  INTEGER                          :: Coord, Surfnum, Surfpos, UsedSiteMapPos, nSites, nSitesRemain
  INTEGER                          :: iInterAtom, xpos, ypos
  REAL                             :: RanNum
  
!===================================================================================================================================
  IF (adsorbates_num.GT.0) THEN
    ! distribute adsorbates randomly on the surface on the correct site and assign surface atom bond order
    dist = 1
    Coord = Adsorption%Coordination(iSpec)
    Surfnum = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    DO WHILE (dist.LE.adsorbates_num) 
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
  ELSE IF (adsorbates_num.LT.0) THEN
    ! remove adsorbates randomly on the surface on the correct site and assign surface atom bond order
    dist = -1
    Coord = Adsorption%Coordination(iSpec)
    nSites = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%nSites(Coord)
    nSitesRemain = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord)
    Surfnum = nSites - nSitesRemain
    DO WHILE (dist.GE.adsorbates_num) 
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
      dist = dist - 1
    END DO
    SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SitesRemain(Coord) = nSitesRemain
  END IF
    
END SUBROUTINE AdjustBackgndAdsNum

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
  USE MOD_Particle_Vars,          ONLY : nSpecies
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
  USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst
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

SUBROUTINE PartitionFuncGas(iSpec, Temp, VarPartitionFuncGas)
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,       ONLY: PlanckConst
USE MOD_DSMC_Vars,          ONLY: DSMC, SpecDSMC, PolyatomMolDSMC
USE MOD_Particle_Vars,      ONLY: BoltzmannConst, Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)           :: iSpec
REAL, INTENT(IN)              :: Temp
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)             :: VarPartitionFuncGas
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, PARAMETER               :: Pi=3.14159265358979323846_8
INTEGER                       :: iPolyatMole, iDOF
REAL                          :: Qtra, Qrot, Qvib, Qelec
!===================================================================================================================================
  Qtra = (2. * Pi * Species(iSpec)%MassIC * BoltzmannConst * Temp / (PlanckConst**2))**(1.5)
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
!   Qelec = 0.
!   DO iDOF=1, SpecDSMC(iSpec)%NumElecLevels
!     Qelec = Qelec + SpecDSMC(iSpec)%ElectronicState(1,iDOF) * EXP(-SpecDSMC(iSpec)%ElectronicState(2,iDOF) / Temp)
!   END DO
  VarPartitionFuncGas = Qtra * Qrot * Qvib! * Qelec

END SUBROUTINE PartitionFuncGas

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

REAL FUNCTION Calc_Beta_Adsorb(iSpec,Xi_Total,Constant_Adsorb)
!===================================================================================================================================
! Calculates the Beta coefficient for molecular (non-associative) adsorption reaction
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,          ONLY : Adsorption
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)            :: iSpec
  REAL, INTENT(IN)               :: Xi_Total
  REAL, INTENT(IN)               :: Constant_Adsorb
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  IF((Adsorption%Ads_Powerfactor(iSpec) - 0.5 + Xi_Total).GT.0.0) THEN
    Calc_Beta_Adsorb = Adsorption%Ads_Prefactor(iSpec) * Constant_Adsorb &
      * BoltzmannConst**(0.5 - Adsorption%Ads_Powerfactor(iSpec) ) * GAMMA(Xi_Total) &
      / GAMMA(Adsorption%Ads_Powerfactor(iSpec) - 0.5 + Xi_Total) 
  ELSE
    Calc_Beta_Adsorb = 0.0
  END IF

END FUNCTION Calc_Beta_Adsorb

REAL FUNCTION Calc_Adsorb_Heat(subsurfxi,subsurfeta,SurfSideID,iSpec,Surfpos,IsAdsorption)
!===================================================================================================================================
! Calculates the Heat of adsorption for given species and given surface position
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_PARTICLE_Vars,      ONLY : BoltzmannConst, nSpecies
  USE MOD_DSMC_Vars,          ONLY : Adsorption, SurfDistInfo, SpecDSMC
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)            :: subsurfxi, subsurfeta, SurfSideID
  INTEGER, INTENT(IN)            :: iSpec, Surfpos
  LOGICAL, INTENT(IN)            :: IsAdsorption
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                        :: Coordination, i, j, iPolyatMole, Indx, Indy
  REAL , ALLOCATABLE             :: x(:)
  INTEGER , ALLOCATABLE          :: m(:)
  REAL                           :: D_AB, D_AX, D_BX, bondorder
  REAL                           :: Heat_A, Heat_B
  REAL                           :: A, B, sigma
!===================================================================================================================================
  Coordination = Adsorption%Coordination(iSpec)
  ALLOCATE( x(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
  ALLOCATE( m(1:SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom) )
  x(:) = 1. / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  m(:) = 1
  Calc_Adsorb_Heat = 0.
  sigma = 0.
  
  IF (Surfpos.GT.0) THEN
  DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
    Indx = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndx(Surfpos,j)
    Indy = SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%BondAtomIndy(Surfpos,j)
    bondorder = 0.
    DO i = 1,nSpecies
      bondorder = bondorder + REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%SurfAtomBondOrder(i,Indx,Indy))
    END DO
    IF (IsAdsorption) THEN
    ! calculate bond order for heat of adsorption (in case of adsorption treatment)
      m(j) = (bondorder + 1)
    ELSE 
    ! calculate bond order for heat of adsorption (in case of desorption treatment)
      m(j) = bondorder
    END IF
  END DO
  END IF
  DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
    sigma = sigma + (2.*x(j) - x(j)**2) * (2-1/m(j))/m(j)
  END DO
  sigma = sigma / REAL(SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom)
  ! Testing if the adsorption particle is an atom or molecule, if molecule: is it polyatomic?
  ! and calculate right heat of adsorption
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    IF(SpecDSMC(iSpec)%PolyatomicMol) THEN
      D_AB = Adsorption%EDissBond(1,iSpec)
      SELECT CASE(Coordination)
      CASE(1) ! hollow (radical with localized electron like NH2)
        Heat_A = 0.
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
          Heat_A = Heat_A + Adsorption%HeatOfAdsZero(iSpec) * (2.*x(j) - x(j)**2) * (2-1/m(j))/m(j)
        END DO
        Calc_Adsorb_Heat = Heat_A**2/(D_AB*Heat_A)
      CASE(2) ! bridge
        IF (Adsorption%DiCoord(iSpec).EQ.1) THEN ! dicoordination (HCOOH --> M--(HC)O-O(H)--M)
          D_AX = Adsorption%EDissBond(1,iSpec)
          Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/m(1))/m(1)
          Heat_A = Heat_A**2/(D_AX*Heat_A)
          D_BX = Adsorption%EDissBond(1,iSpec)
          Heat_B = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(2) - x(2)**2) * (2-1/m(2))/m(2)
          Heat_B = Heat_B**2/(D_BX*Heat_B)
          Calc_Adsorb_Heat = ( A*B*( A + B ) + D_AB*( A - B )**2 ) / ( A*B + D_AB*( A + B ) )
        ELSE IF (Adsorption%DiCoord(iSpec).EQ.2) THEN ! chelate binding (NO2 --> M--O-N-O--M)
          D_AX = Adsorption%EDissBond(1,iSpec)
          Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/m(1))/m(1)
          Heat_A = Heat_A**2/(D_AX*Heat_A)
          D_BX = Adsorption%EDissBond(1,iSpec)
          Heat_B = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(2) - x(2)**2) * (2-1/m(2))/m(2)
          Heat_B = Heat_B**2/(D_BX*Heat_B)
          A = Heat_A**2 * ( Heat_A + 2*Heat_B ) / ( Heat_A + Heat_B )**2
          B = Heat_B**2 * ( Heat_B + 2*Heat_A ) / ( Heat_A + Heat_B )**2
          Calc_Adsorb_Heat = A + B
        END IF
      CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/(m(1)))/(m(1))
        Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A)
      END SELECT
    ELSE 
      D_AB = Adsorption%EDissBond(1,iSpec)
      SELECT CASE(Coordination)
      CASE(1) ! hollow (radical with localized electron like C-H)
        Heat_A = 0.
        DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
          Heat_A = Heat_A + Adsorption%HeatOfAdsZero(iSpec) * (2.*x(j) - x(j)**2) * (2-1/m(j))/m(j)
        END DO
        Calc_Adsorb_Heat = Heat_A**2/(D_AB*Heat_A)
      CASE(2) ! bridge (closed shell like O2)
        IF (Adsorption%DiCoord(iSpec).EQ.1) THEN
          Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/m(1))/m(1)
          Heat_B = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(2) - x(2)**2) * (2-1/m(2))/m(2)
          A = Heat_A !Heat_A**2 * ( Heat_A + 2*Heat_B ) / ( Heat_A + Heat_B )**2 
          B = Heat_B !Heat_B**2 * ( Heat_B + 2*Heat_A ) / ( Heat_A + Heat_B )**2
!           A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/(m(1)+1))/(m(1)+1)
!           B = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(2) - x(2)**2) * (2-1/(m(2)+1))/(m(2)+1)
          Heat_A = ( A*B*( A + B ) + D_AB*( A - B )**2 ) / ( A*B + D_AB*( A + B ) )
          Calc_Adsorb_Heat = Heat_A! * sigma
        ELSE
          Heat_A = 0.
          DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
            Heat_A = Heat_A + Adsorption%HeatOfAdsZero(iSpec) * (2.*x(j) - x(j)**2) * (2-1/m(j))/m(j)
          END DO
          Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A)
        END IF
      CASE(3) ! on-top (closed shell or open shell with unlocalized electron like CO)
        Heat_A = Adsorption%HeatOfAdsZero(iSpec) * (2.*x(1) - x(1)**2) * (2-1/(m(1)))/(m(1))
        Calc_Adsorb_Heat = Heat_A**2/(D_AB+Heat_A)
      END SELECT
    END IF
  ELSE
    DO j = 1,SurfDistInfo(subsurfxi,subsurfeta,SurfSideID)%AdsMap(Coordination)%nInterAtom
      Calc_Adsorb_Heat = Calc_Adsorb_Heat + Adsorption%HeatOfAdsZero(iSpec) * (2.*x(j) - x(j)**2) * (2-1/m(j))/m(j)
    END DO
  END IF
  
  DEALLOCATE(x,m)

END FUNCTION Calc_Adsorb_Heat

REAL FUNCTION Calc_E_Act(Heat_A,Heat_B,Heat_AB,D_AB,D_A,D_B,IsAdsorption)
!===================================================================================================================================
! Calculates the Heat of adsorption for given species and given surface position
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)               :: Heat_A, Heat_B, Heat_AB, D_AB, D_A, D_B
  LOGICAL, INTENT(IN)            :: IsAdsorption!, IsMolecular
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                           :: Delta_H
!===================================================================================================================================
  Delta_H = ( Heat_AB -Heat_A -Heat_B ) + ( D_AB -D_A -D_B )
  Calc_E_Act = 0.5 * ( Delta_H + (Heat_A*Heat_B / (Heat_A+Heat_B)) )
  IF (Calc_E_Act.LT.0.) Calc_E_Act = 0.
  IF (.NOT.IsAdsorption) THEN
    Calc_E_Act = Calc_E_Act - Delta_H
    IF (Calc_E_Act.LT.0.) Calc_E_Act = 0.
  END IF

END FUNCTION Calc_E_Act

END MODULE MOD_DSMC_SurfModel_Tools
