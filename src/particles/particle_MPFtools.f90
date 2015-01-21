#include "boltzplatz.h"

MODULE MOD_part_MPFtools
!===================================================================================================================================
! CONTAINS THE vMPF part
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE 
  PRIVATE

INTERFACE SplitParticle
  MODULE PROCEDURE SplitParticle
END INTERFACE

INTERFACE MergeParticles
  MODULE PROCEDURE MergeParticles
END INTERFACE

INTERFACE DefinePolyVec
  MODULE PROCEDURE DefinePolyVec
END INTERFACE

INTERFACE DefineSplitVec
  MODULE PROCEDURE DefineSplitVec
END INTERFACE

INTERFACE StartParticleMerge
  MODULE PROCEDURE StartParticleMerge
END INTERFACE

INTERFACE MapToGeo
  MODULE PROCEDURE MapToGeo
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: SplitParticle, MergeParticles, DefinePolyVec, DefineSplitVec, StartParticleMerge, &
            MapToGeo
!===================================================================================================================================

CONTAINS   

SUBROUTINE StartParticleMerge()                                                                
!===================================================================================================================================
! Particle Merge routine
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : doParticleMerge, nSpecies,vMPFMergeParticleTarget,vMPF_SpecNumElem,vMPFSplitParticleTarget
  USE MOD_Mesh_Vars,     ONLY : nElems
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                         :: iElem, iSpec
!===================================================================================================================================
DO iElem = 1, nElems
  DO iSpec= 1, nSpecies
    IF ((vMPFMergeParticleTarget.GT.0).AND.(vMPF_SpecNumElem(iElem,iSpec).GT.vMPFMergeParticleTarget*2)) THEN
      CALL MergeParticles(iElem, vMPFMergeParticleTarget, vMPF_SpecNumElem(iElem,iSpec), iSpec)
    ELSE IF((vMPFSplitParticleTarget.GT.0).AND.(vMPF_SpecNumElem(iElem,iSpec).LT.vMPFSplitParticleTarget/2)) THEN
      CALL MergeParticles(iElem, vMPFSplitParticleTarget, vMPF_SpecNumElem(iElem,iSpec), iSpec)
    END IF
  END DO
END DO
doParticleMerge=.false.
END SUBROUTINE StartParticleMerge
                                                                                        
                                                                                                   
SUBROUTINE SplitParticle(iPart, deltaE)                                                                
!===================================================================================================================================
! Split Particles
!===================================================================================================================================
! MODULES
  USE MOD_Globals,        ONLY : Abort
  USE MOD_Particle_Vars,  ONLY : PDM, PartState, RandomVec, NumRanVec, PartSpecies, PartMPF, PEM, Species  
  USE MOD_DSMC_Vars,      ONLY : useDSMC, CollisMode, PartStateIntEn                                                      
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN)              :: iPart
  REAL, INTENT(IN)                :: deltaE
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                         :: PositionNbr, iVec
  REAL                            :: beta, iRan
  REAL                            :: v_old(3)
!===================================================================================================================================

  v_old(1:3) = PartState(iPart,4:6)

!.... Get free particle index for the new particle produced
  PDM%ParticleVecLength = PDM%ParticleVecLength + 1
  PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1 
  PositionNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
  IF (PositionNbr.EQ.0) THEN
    CALL Abort(&
       __STAMP__,&
      'ERROR in SplitParticle: New Particle Number greater max Part Num!')
  END IF

!Set new particle parameters
  PDM%ParticleInside(PositionNbr) = .true.
  PartSpecies(PositionNbr) = PartSpecies(iPart)
  PartState(PositionNbr,1:3) = PartState(iPart,1:3)
  IF (useDSMC.AND.(CollisMode.GT.1)) THEN
    PartStateIntEn(PositionNbr, 1) = PartStateIntEn(iPart, 1)
    PartStateIntEn(PositionNbr, 2) =   PartStateIntEn(iPart, 2)
  END IF
  PEM%Element(PositionNbr) = PEM%Element(iPart)

!calulating beta = sqrt(deltaE/MPF_old)
  beta = SQRT(2*deltaE/(PartMPF(iPart)*Species(PartSpecies(iPart))%MassIC))

!set new MPFs
  PartMPF(iPart) =  PartMPF(iPart) / 2
  PartMPF(PositionNbr) = PartMPF(iPart)

!set new velocity v1
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec*iRan + 1)
  PartState(iPart,4:6) = v_old(1:3) - beta * RandomVec(iVec, 1:3)
  PartState(PositionNbr,4:6) = v_old(1:3) + beta * RandomVec(iVec, 1:3)

END SUBROUTINE SplitParticle


SUBROUTINE MergeParticles(iElem, NumFinPart, SpecNum, SpecID)                                                                
!===================================================================================================================================
! Merge Particles
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Particle_Vars
  USE Levenberg_Marquardt
  USE MOD_Eval_xyz,            ONLY:eval_xyz_elemcheck
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                      
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES                                                                        
  INTEGER,INTENT(IN)              :: iElem, SpecNum, SpecID
  INTEGER,INTENT(IN)              :: NumFinPart
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                  :: CellTemp
  INTEGER               :: iPart, iLoop, iLoop2, PositionNbr
!===================================================================================================================================

  ALLOCATE(vMPFOldVelo(3,SpecNum))
  ALLOCATE(vMPFOldPos(3,SpecNum), vMPFOldMPF(SpecNum))
  iLoop2=1
  ALLOCATE(PartStateMap(SpecNum,3))
  IF (SpecNum.GE.NumFinPart) THEN
    ALLOCATE(PartStatevMPFSpec(SpecNum))
  ELSE
    ALLOCATE(PartStatevMPFSpec(NumFinPart))
  END IF
  iPart = PEM%pStart(iElem)                        
  DO iLoop = 1, PEM%pNumber(iElem)
    IF(PartSpecies(iPart).EQ.SpecID) THEN
      CALL Eval_XYZ_ElemCheck(PartState(iPart,1:3), PartStateMap(iLoop2,1:3), iElem)
      PartStatevMPFSpec(iLoop2) = iPart
      iLoop2 = iLoop2 + 1
    END IF
    iPart = PEM%pNext(iPart)    
  END DO

  SWRITE(*,*) 'Start Particle Split/Merge'
  IF (NumFinPart.GT.SpecNum) THEN
    DO iLoop = 1 , NumFinPart - SpecNum
      PDM%ParticleVecLength = PDM%ParticleVecLength + 1
      PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1 
      PositionNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
      IF (PositionNbr.EQ.0) THEN
        CALL Abort(&
           __STAMP__,&
          'ERROR in SplitParticle: New Particle Number greater max Part Num!')
      END IF
      PartStatevMPFSpec(SpecNum + iLoop) = PositionNbr
      !Set new particle parameters
      PDM%ParticleInside(PositionNbr) = .true.
      PartSpecies(PositionNbr) = SpecID
      PEM%Element(PositionNbr) = iElem
    END DO
  END IF


  CALL SplitRegion(SpecNum)
  CALL DeleteParticlesMPF(NumFinPart, CellTemp, SpecNum, SpecID)
  CALL SetMPFParticlePosCube(iElem, NumFinPart)
  CALL SetNewvMPF(NumFinPart)
  CALL SetNewVelos(NumFinPart, CellTemp, SpecNum, SpecID)
  SWRITE(*,*) 'Finish Particle Split/Merge'

  DEALLOCATE(PartStateMap, PartStatevMPFSpec, vMPFOldVelo, vMPFOldPos, vMPFOldMPF)

END SUBROUTINE MergeParticles


SUBROUTINE fcn (m, n, x, fvec, fjac, iflag)
!===================================================================================================================================
! do not know
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
  INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
! INPUT VARIABLES
  INTEGER, INTENT(IN)        :: m, n
  REAL (dp), INTENT(IN)      :: x(:)
  REAL (dp), INTENT(IN OUT)  :: fvec(:)
  INTEGER, INTENT(IN OUT)    :: iflag
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL (dp), INTENT(OUT)     :: fjac(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  
  IF (iflag == 1) CALL ssqfcn (m, n, x, fvec)
  IF (iflag == 2) CALL ssqjac (m, n, fjac)
  RETURN
END SUBROUTINE fcn

SUBROUTINE ssqjac (m, n, fjac)
!===================================================================================================================================
! ssqjac
!===================================================================================================================================
! MODULES
  USE Levenberg_Marquardt
  USE MOD_Particle_Vars,    ONLY:vMPFPolyPoint, vMPF_OrderVec
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)     :: m, n
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL (dp), INTENT(OUT)  :: fjac(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: i_Sol, i_DOFIn
!===================================================================================================================================

  DO i_Sol = 1 ,m
    DO i_DOFIn =1, n   
      FJAC(i_Sol,i_DOFIn) = - vMPFPolyPoint(1,i_Sol)**(vMPF_OrderVec(1,i_DOFIn)) &
                            * vMPFPolyPoint(2,i_Sol)**(vMPF_OrderVec(2,i_DOFIn)) &
                            * vMPFPolyPoint(3,i_Sol)**(vMPF_OrderVec(3,i_DOFIn))
    END DO
  END DO

  RETURN
END SUBROUTINE ssqjac

SUBROUTINE ssqfcn (m, n, x, fvec)
!===================================================================================================================================
! ssqfcn?!
!===================================================================================================================================
! MODULES
  USE Levenberg_Marquardt
  USE MOD_Particle_Vars,    ONLY:vMPFPolyPoint, vMPF_OrderVec, vMPFPolySol,vMPFOldMPF
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)     :: m, n
  REAL (dp), INTENT(IN)   :: x(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL (dp), INTENT(OUT)  :: fvec(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                 :: i_Sol, i_DOFIn
!===================================================================================================================================

  DO i_Sol = 1 ,m
    FVEC(i_Sol) = vMPFPolySol(i_Sol) 
    DO i_DOFIn =1, n   
      FVEC(i_Sol) = FVEC(i_Sol) - x(i_DOFIn) *vMPFPolyPoint(1,i_Sol)**(vMPF_OrderVec(1,i_DOFIn)) &
                        *vMPFPolyPoint(2,i_Sol)**(vMPF_OrderVec(2,i_DOFIn)) &
                        *vMPFPolyPoint(3,i_Sol)**(vMPF_OrderVec(3,i_DOFIn))
    END DO
    FVEC(i_Sol)=FVEC(i_Sol)*vMPFOldMPF(i_Sol)
  END DO

  RETURN
  END SUBROUTINE ssqfcn

SUBROUTINE DefinePolyVec(VecOrder)                                                                
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : vMPF_OrderVec                                                      
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: VecOrder
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: DOF_Poly
  INTEGER                         :: x_dir, y_dir, z_dir, sum_dir
!===================================================================================================================================

  DOF_Poly = (VecOrder+3)*(VecOrder+2)*(VecOrder+1)/6
  sum_dir = 1
  ALLOCATE(vMPF_OrderVec(3,DOF_Poly))
  DO x_dir = 0, VecOrder
    DO y_dir = 0, VecOrder
      DO z_dir = 0, VecOrder
        IF ((x_dir+y_dir+z_dir).GT.VecOrder) CYCLE
        vMPF_OrderVec(1, sum_dir) = x_dir
        vMPF_OrderVec(2, sum_dir) = y_dir
        vMPF_OrderVec(3, sum_dir) = z_dir
        sum_dir = sum_dir + 1
      END DO
    END DO
  END DO

END SUBROUTINE DefinePolyVec


SUBROUTINE DefineSplitVec(SplitOrder)                                                                !
!===================================================================================================================================
! blabla
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : vMPF_SplitVec ,vMPF_SplitVecBack                                                     !
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: SplitOrder
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: DOF_Split
  INTEGER                         :: x_dir, y_dir, z_dir, sum_dir
!===================================================================================================================================
  DOF_Split = (SplitOrder+1)**3
  sum_dir = 1
  ALLOCATE(vMPF_SplitVec(3,DOF_Split))
  ALLOCATE(vMPF_SplitVecBack(SplitOrder+1,SplitOrder+1,SplitOrder+1))
  DO x_dir = 0, SplitOrder
    DO y_dir = 0, SplitOrder
      DO z_dir = 0, SplitOrder
        vMPF_SplitVec(1, sum_dir) = x_dir
        vMPF_SplitVec(2, sum_dir) = y_dir
        vMPF_SplitVec(3, sum_dir) = z_dir
        vMPF_SplitVecBack(x_dir+1,y_dir+1,z_dir+1) = sum_dir
        sum_dir = sum_dir + 1
      END DO
    END DO
  END DO

END SUBROUTINE DefineSplitVec


SUBROUTINE SplitRegion(SpecNum)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : vMPFMergeCellSplitOrder,PartStateMap, vMPF_SplitVec, vMPFPolyPoint, &
                                 vMPFPolySol, PartMPF, vMPF_oldMPFSum, vMPF_SplitVecBack, PartStatevMPFSpec
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
  INTEGER, INTENT(IN)   :: SpecNum
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  REAL, ALLOCATABLE     :: RegPartNum(:)
  INTEGER               :: iPart, PolOrder, x_cube, y_cube, z_cube
  REAL                  :: ZoneLenght, RegPartSum
!===================================================================================================================================

  ZoneLenght = 2.0/(vMPFMergeCellSplitOrder+1)
  PolOrder = (vMPFMergeCellSplitOrder+1)**3
  ALLOCATE(RegPartNum(PolOrder))
  RegPartNum = 0                        

  DO iPart = 1, SpecNum
    x_cube = MIN(INT((PartStateMap(iPart,1)+1)/2*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    y_cube = MIN(INT((PartStateMap(iPart,2)+1)/2*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    z_cube = MIN(INT((PartStateMap(iPart,3)+1)/2*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    RegPartNum(vMPF_SplitVecBack(x_cube,y_cube,z_cube)) = RegPartNum(vMPF_SplitVecBack(x_cube,y_cube,z_cube))  &
                          + PartMPF(PartStatevMPFSpec(iPart))
  END DO

  ALLOCATE(vMPFPolyPoint(3,PolOrder))
  ALLOCATE(vMPFPolySol(PolOrder))
  RegPartSum = SUM(RegPartNum)

  vMPF_oldMPFSum = RegPartSum
  vMPFPolySol = RegPartNum/RegPartSum
  vMPFPolyPoint = vMPF_SplitVec*ZoneLenght-1+ZoneLenght/2

END SUBROUTINE SplitRegion


SUBROUTINE DeleteParticlesMPF(FinPartNum, Temp, SpecNum, SpecID)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, vMPF_oldEngSum, vMPF_oldMomSum ,Species, PartMPF, PDM, &
                                  BoltzmannConst, PartStatevMPFSpec, vMPFOldPos, vMPFOldVelo, vMPFOldMPF
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
  INTEGER, INTENT(IN)   :: FinPartNum, SpecNum, SpecID
  REAL, INTENT(OUT)     :: Temp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER               :: iLoop
  REAL                  :: PartV_2(3), PartV2(3), RealPartNum
!===================================================================================================================================

  vMPF_oldMomSum = 0.0
  vMPF_oldEngSum = 0.0 
  PartV_2 = 0.0
  PartV2 = 0.0
  RealPartNum = 0.0
                          
  DO iLoop = 1, SpecNum
    vMPF_oldEngSum = vMPF_oldEngSum+  0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iLoop)) &
              * (PartState(PartStatevMPFSpec(iLoop),4)**2 + PartState(PartStatevMPFSpec(iLoop),5)**2  &
               + PartState(PartStatevMPFSpec(iLoop),6)**2)
    vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) + Species(SpecID)%MassIC &
                    * PartMPF(PartStatevMPFSpec(iLoop)) * PartState(PartStatevMPFSpec(iLoop),4:6)  
    PartV_2 = PartV_2 + PartState(PartStatevMPFSpec(iLoop),4:6) * PartMPF(PartStatevMPFSpec(iLoop))
    PartV2 = PartV2 + PartState(PartStatevMPFSpec(iLoop),4:6)**2 * PartMPF(PartStatevMPFSpec(iLoop))
    RealPartNum = RealPartNum + PartMPF(PartStatevMPFSpec(iLoop))
    vMPFOldVelo(1:3, iLoop) = PartState(PartStatevMPFSpec(iLoop),4:6)
    vMPFOldPos(1:3, iLoop) = PartState(PartStatevMPFSpec(iLoop),1:3)
    vMPFOldMPF(iLoop) = PartMPF(PartStatevMPFSpec(iLoop))
    IF (iLoop.GT.FinPartNum) PDM%ParticleInside(PartStatevMPFSpec(iLoop)) = .false.
  END DO
  PartV_2 = (PartV_2/RealPartNum)**2
  PartV2 = PartV2/RealPartNum
  Temp = Species(SpecID)%MassIC/(3*BoltzmannConst)*(PartV2(1)+PartV2(2)+PartV2(3)-PartV_2(1)-PartV_2(2)-PartV_2(3))

END SUBROUTINE DeleteParticlesMPF


#ifdef DONTCOMPILETHIS
SUBROUTINE SetMPFParticlePos(FinPartNum,x)
!===================================================================================================================================
!
!===================================================================================================================================
! Modules
  USE MOD_Particle_Vars, ONLY : PartState,vMPFMergePolyOrder, vMPF_OrderVec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
  INTEGER, INTENT(IN)   :: FinPartNum
  REAL, INTENT(IN)      :: x(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER               :: iPart, iLoop, iDOF, DOF_LMInput
  REAL                  :: RandVac(3), ProbPos,  iRan
!===================================================================================================================================

  DOF_LMInput = (vMPFMergePolyOrder+1)*(vMPFMergePolyOrder+2)*(vMPFMergePolyOrder+3)/6                 
  DO iLoop = 1, FinPartNum
  ProbPos = 0.0
    CALL RANDOM_NUMBER(RandVac)  
    RandVac = RandVac * 2.0 - 1.0
    DO iDOF =1, DOF_LMInput  
      ProbPos = ProbPos + x(iDOF) *RandVac(1)**(vMPF_OrderVec(1,iDOF)) &
                        *RandVac(2)**(vMPF_OrderVec(2,iDOF)) &
                        *RandVac(3)**(vMPF_OrderVec(3,iDOF))
    END DO
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GE.ProbPos)
      ProbPos = 0.0
      CALL RANDOM_NUMBER(RandVac)  
      RandVac = RandVac * 2.0 - 1.0

      DO iDOF =1, DOF_LMInput  
        ProbPos = ProbPos + x(iDOF) *RandVac(1)**(vMPF_OrderVec(1,iDOF)) &
                          *RandVac(2)**(vMPF_OrderVec(2,iDOF)) &
                          *RandVac(3)**(vMPF_OrderVec(3,iDOF))
      END DO
      CALL RANDOM_NUMBER(iRan)
    END DO
    PartState(iLoop, 1:3) = RandVac
  END DO

END SUBROUTINE SetMPFParticlePos
#endif /*DONTCOMPILETHIS*/


SUBROUTINE SetMPFParticlePosCube(iElem, FinPartNum)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, vMPFMergeCellSplitOrder &
                        , vMPFPolySol, vMPF_SplitVecBack, GEO, PartStatevMPFSpec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
  INTEGER, INTENT(IN)   :: iElem, FinPartNum
!----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
  INTEGER               :: iLoop, x_cube, y_cube, z_cube, iNode, iLoop2
  REAL                  :: RandVac(3), ProbPos,  iRan, P(3,8)
!===================================================================================================================================

  DO iNode = 1,8
    P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
  END DO
  iLoop2 = 1                      

  DO iLoop = 1, FinPartNum
    ProbPos = 0.0
    CALL RANDOM_NUMBER(RandVac)  
    x_cube = MIN(INT(RandVac(1)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    y_cube = MIN(INT(RandVac(2)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    z_cube = MIN(INT(RandVac(3)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
    ProbPos = vMPFPolySol(vMPF_SplitVecBack(x_cube,y_cube,z_cube))
    CALL RANDOM_NUMBER(iRan)
    DO WHILE (iRan.GE.ProbPos)
      ProbPos = 0.0
      CALL RANDOM_NUMBER(RandVac)  
      x_cube = MIN(INT(RandVac(1)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
      y_cube = MIN(INT(RandVac(2)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
      z_cube = MIN(INT(RandVac(3)*(vMPFMergeCellSplitOrder+1)+1),(vMPFMergeCellSplitOrder+1))
      ProbPos = vMPFPolySol(vMPF_SplitVecBack(x_cube,y_cube,z_cube))
      CALL RANDOM_NUMBER(iRan)
    END DO
    RandVac = RandVac * 2.0 - 1.0
    PartState(PartStatevMPFSpec(iLoop), 1:3) = MapToGeo(RandVac, P)
  END DO

END SUBROUTINE SetMPFParticlePosCube


SUBROUTINE SetNewVelos(NewPartNum, Temp, SpecNum, SpecID)                                                                
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars
  USE Levenberg_Marquardt
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
  INTEGER,INTENT(IN)              :: NewPartNum, SpecNum, SpecID
  REAL,INTENT(IN)                 :: Temp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: iPart, iLoop, iDir, DOF_LMInput, info
  INTEGER, ALLOCATABLE            :: iwa(:)
  DOUBLE PRECISION, ALLOCATABLE   :: fjac(:,:), fvec(:), x(:) !!evtl double
  DOUBLE PRECISION                :: tol
  REAL                            :: iRan
!===================================================================================================================================

  tol = SQRT( EPSILON(tol) )
  DOF_LMInput = (vMPFMergePolyOrder+1)*(vMPFMergePolyOrder+2)*(vMPFMergePolyOrder+3)/6
  ALLOCATE(iwa(DOF_LMInput), & 
           x(DOF_LMInput), &
           fjac(SpecNum,DOF_LMInput), &
           fvec(SpecNum))
  DEALLOCATE(vMPFPolyPoint,vMPFPolySol)
  ALLOCATE(vMPFPolyPoint(3,SpecNum))
  ALLOCATE(vMPFPolySol(SpecNum))

  IF (vMPF_velocityDistribution.EQ.'OVDR') THEN
    ALLOCATE(vMPFOldBrownVelo(SpecNum,3))
  END IF

  vMPFPolyPoint = vMPFOldPos
  DO iDir = 1, 3
    vMPFPolySol(:) = vMPFOldVelo(iDir,:)
    x = 1.0
    CALL lmder1 (fcn, SpecNum, DOF_LMInput, x, fvec, fjac, tol, info, iwa)

  !Sample old brownian velo
    IF (vMPF_velocityDistribution.EQ.'OVDR') THEN
      DO iPart = 1, SpecNum
        vMPFOldBrownVelo(iPart,iDir) = vMPFOldVelo(iDir,iPart)
        DO iLoop =1, DOF_LMInput  
          vMPFOldBrownVelo(iPart,iDir)  = vMPFOldBrownVelo(iPart,iDir)  - x(iLoop) &
                            *vMPFPolyPoint(1,iPart)**(vMPF_OrderVec(1,iLoop)) &
                            *vMPFPolyPoint(2,iPart)**(vMPF_OrderVec(2,iLoop)) &
                            *vMPFPolyPoint(3,iPart)**(vMPF_OrderVec(3,iLoop))
        END DO
      END DO
    END IF

    DO iPart = 1, NewPartNum - 1
      PartState(PartStatevMPFSpec(iPart),iDir + 3) = 0.0      
      DO iLoop =1, DOF_LMInput  
        PartState(PartStatevMPFSpec(iPart),iDir + 3) = PartState(PartStatevMPFSpec(iPart),iDir + 3) + x(iLoop) &
                          *PartState(PartStatevMPFSpec(iPart), 1)**(vMPF_OrderVec(1,iLoop)) &
                          *PartState(PartStatevMPFSpec(iPart), 2)**(vMPF_OrderVec(2,iLoop)) &
                          *PartState(PartStatevMPFSpec(iPart), 3)**(vMPF_OrderVec(3,iLoop))
      END DO
    END DO
  END DO


  DO iPart = 1, NewPartNum -1
    vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
            * (PartState(PartStatevMPFSpec(iPart),4)**2 + PartState(PartStatevMPFSpec(iPart),5)**2 &
            + PartState(PartStatevMPFSpec(iPart),6)**2)
    vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
                         * PartState(PartStatevMPFSpec(iPart),4:6)
  END DO
  PartState(PartStatevMPFSpec(NewPartNum),4:6) = vMPF_oldMomSum(1:3) &
                  /(Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)))
  vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
              * (PartState(PartStatevMPFSpec(NewPartNum),4)**2 + PartState(PartStatevMPFSpec(NewPartNum),5)**2 &
              + PartState(PartStatevMPFSpec(NewPartNum),6)**2)
  vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
                         * PartState(PartStatevMPFSpec(NewPartNum),4:6)

  IF (vMPF_velocityDistribution.EQ.'MBDR') THEN
    CALL SetNewTemp_2(Temp, NewPartNum)
  ELSE IF (vMPF_velocityDistribution.EQ.'OVDR') THEN
    CALL SetNewDistrVelo(NewPartNum, 100, SpecNum)
  END IF 

  IF (vMPF_oldEngSum.GT.0) THEN
    CALL RANDOM_NUMBER(iRan)
    iPart = INT(NewPartNum * iRan +1)
    CALL SplitParticle(PartStatevMPFSpec(iPart), vMPF_oldEngSum)
  ELSE IF (vMPF_oldEngSum.LT.0) THEN
    WRITE(*,*) 'Particles could not be merged/split! Continue without merge/split process.'
    DO iPart = 1, SpecNum
      PartState(PartStatevMPFSpec(iPart), 1:3) = vMPFOldPos(1:3,iPart)
      PartState(PartStatevMPFSpec(iPart), 4:6) = vMPFOldVelo(1:3,iPart)
      PartMPF(PartStatevMPFSpec(iPart)) =  vMPFOldMPF(iPart)
      PDM%ParticleInside(PartStatevMPFSpec(iPart)) = .true.
    END DO
  END IF

  DEALLOCATE(iwa, x, fjac, fvec, vMPFPolyPoint, vMPFPolySol)

  IF (vMPF_velocityDistribution.EQ.'OVDR') THEN
    DEALLOCATE(vMPFOldBrownVelo)
  END IF

END SUBROUTINE SetNewVelos


FUNCTION MapToGeo(xi,P)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL,INTENT(IN)          :: xi(3)      ! 
  REAL,INTENT(IN)          :: P(3,8)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                     :: MapToGeo(3)
!===================================================================================================================================

  MapToGeo =0.125 *(P(:,1)*(1-xi(1)) * (1-xi(2)) * (1-xi(3))  &
                + P(:,2)*(1+xi(1)) * (1-xi(2)) * (1-xi(3))  &
                + P(:,3)*(1+xi(1)) * (1+xi(2)) * (1-xi(3))  &
                + P(:,4)*(1-xi(1)) * (1+xi(2)) * (1-xi(3))  &
                + P(:,5)*(1-xi(1)) * (1-xi(2)) * (1+xi(3))  &
                + P(:,6)*(1+xi(1)) * (1-xi(2)) * (1+xi(3))  &
                + P(:,7)*(1+xi(1)) * (1+xi(2)) * (1+xi(3))  &
                + P(:,8)*(1-xi(1)) * (1+xi(2)) * (1+xi(3))) 

END FUNCTION MapToGeo 


#ifdef DONTCOMPILETHIS
SUBROUTINE SetNewTemp(PartIndx, Temp, iPart)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, Species, PartSpecies, vMPF_oldEngSum, vMPF_oldMomSum, PartMPF, BoltzmannConst
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: PartIndx, iPart
  REAL,INTENT(IN)                 :: Temp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  REAL                            :: SumRan, RandVal(2), ran1, ran2
!===================================================================================================================================

  SumRan = 2
  DO WHILE (SumRan .GT. 1)
   CALL RANDOM_NUMBER(RandVal)
   ran1 = 2*RandVal(1) - 1
   ran2 = 2*RandVal(2) - 1
   SumRan = ran1**2 + ran2**2
  END DO
  PartState(PartIndx, 4) = PartState(PartIndx, 4) & 
              + ran1*SQRT(-2*BoltzmannConst*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)
  PartState(PartIndx, 5) = PartState(PartIndx, 5) &
              + ran2*SQRT(-2*BoltzmannConst*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)

  SumRan = 2
  DO WHILE (SumRan .GT. 1)
   CALL RANDOM_NUMBER(RandVal)
   ran1 = 2*RandVal(1) - 1
   ran2 = 2*RandVal(2) - 1
   SumRan = ran1**2 + ran2**2
  END DO
  PartState(PartIndx, 6) =PartState(PartIndx, 6) &
              + ran1*SQRT(-2*BoltzmannConst*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)

  vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(PartSpecies(PartIndx))%MassIC * PartMPF(PartIndx) &
          * (PartState(PartIndx,4)**2 + PartState(PartIndx,5)**2 &
          + PartState(PartIndx,6)**2)
  vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(PartSpecies(PartIndx))%MassIC * PartMPF(PartIndx) &
                       * PartState(PartIndx,4:6)
  IF(vMPF_oldEngSum.lT.0) then
   print*, 'mist: ', iPart
    read*
  end if

END SUBROUTINE SetNewTemp
#endif /*DONTCOMPILETHIS*/


SUBROUTINE SetNewvMPF(FinPartNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartMPF, vMPF_oldMPFSum, PartStatevMPFSpec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: FinPartNum
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: iLoop
  REAL                            :: NewMPF
!===================================================================================================================================

  NewMPF = vMPF_oldMPFSum/FinPartNum
  DO iLoop = 1, FinPartNum  
    PartMPF(PartStatevMPFSpec(iLoop)) = NewMPF
  END DO

END SUBROUTINE SetNewvMPF


SUBROUTINE SetNewTemp_2(Temp, NewPartNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, Species, PartSpecies, vMPF_oldEngSum, vMPF_oldMomSum, PartMPF, PartStatevMPFSpec, &
                                 BoltzmannConst
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: NewPartNum
  REAL,INTENT(IN)                 :: Temp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  REAL                            :: SumRan, RandVal(2), ran1, ran2, v2_sum, v_sum(1:3), maxwellfac, v_merge
  INTEGER                         :: iPart, distnum, SpecID, iLoop, iDir, iPart2
  REAL, ALLOCATABLE               :: PartTemp(:,:)
  REAL                      :: TempPartVelo(NewPartNum,3)
!===================================================================================================================================

  ALLOCATE(PartTemp(NewPartNum,1:3))
  v_sum(1:3) = 0.0
  v2_sum = 0.0
  SpecID = PartSpecies(PartStatevMPFSpec(1))
  iPart = 1
  DO WHILE (iPart .le. NewPartNum) 
    DO distnum = 1, 3
      CALL RANDOM_NUMBER(RandVal)
      ran1 = 2.0*RandVal(1)-1.0
      ran2 = 2.0*RandVal(2)-1.0
      SumRan= ran1**2+ran2**2
      DO WHILE ((SumRan.LE.0).OR.(SumRan.GE.1))
        CALL RANDOM_NUMBER(RandVal)
        ran1 = 2.0*RandVal(1)-1.0
        ran2 = 2.0*RandVal(2)-1.0
        SumRan= ran1**2+ran2**2
      END DO
      PartTemp(iPart,distnum) = ran1*SQRT(-2*LOG(SumRan)/SumRan)
    END DO                    
    v_sum(1:3) = v_sum(1:3) + PartTemp(iPart,1:3)
    v2_sum = v2_sum + PartTemp(iPart,1)**2+PartTemp(iPart,2)**2+PartTemp(iPart,3)**2
    iPart = iPart + 1
  END DO
  v_sum(1:3) = v_sum(1:3) / (NewPartNum) 
  v2_sum = v2_sum / (NewPartNum) 
  maxwellfac = SQRT(3. * BoltzmannConst * Temp/ &              ! velocity of maximum
                 (Species(SpecID)%MassIC*v2_sum))


  DO iPart =1, NewPartNum -1
    vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
            * (PartState(PartStatevMPFSpec(iPart),4)**2 + PartState(PartStatevMPFSpec(iPart),5)**2 &
            + PartState(PartStatevMPFSpec(iPart),6)**2)
    vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) + Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
                         * PartState(PartStatevMPFSpec(iPart),4:6)   

    PartState(PartStatevMPFSpec(iPart),4:6) = PartState(PartStatevMPFSpec(iPart),4:6) &
                        + (PartTemp(iPart,1:3) - v_sum(1:3)) * maxwellfac 
    vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
            * (PartState(PartStatevMPFSpec(iPart),4)**2 + PartState(PartStatevMPFSpec(iPart),5)**2 &
            + PartState(PartStatevMPFSpec(iPart),6)**2)
    vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
                         * PartState(PartStatevMPFSpec(iPart),4:6)  
  END DO

  vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
          * (PartState(PartStatevMPFSpec(NewPartNum),4)**2 + PartState(PartStatevMPFSpec(NewPartNum),5)**2 &
          + PartState(PartStatevMPFSpec(NewPartNum),6)**2)
  vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) + Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
                         * PartState(PartStatevMPFSpec(NewPartNum),4:6)  

  PartState(PartStatevMPFSpec(NewPartNum),4:6) =vMPF_oldMomSum(1:3) &
                        / (Species(SpecID)%MassIC*PartMPF(PartStatevMPFSpec(NewPartNum)) )
  vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
          * (PartState(PartStatevMPFSpec(NewPartNum),4)**2 + PartState(PartStatevMPFSpec(NewPartNum),5)**2 &
          + PartState(PartStatevMPFSpec(NewPartNum),6)**2)

  IF (vMPF_oldEngSum.LT.0.0) THEN
    DO iPart = 1, NewPartNum  
      TempPartVelo(iPart,1:3) = PartState(PartStatevMPFSpec(iPart),4:6)
    END DO

    iLoop = 0
     DO WHILE (vMPF_oldEngSum.LT.0.0)
      CALL RANDOM_NUMBER(ran1)    
      iDir = INT(3*ran1 + 1) + 3
      iPart2 = MAXLOC(TempPartVelo(:,iDir-3),1)
      iPart = MINLOC(TempPartVelo(:,iDir-3),1)
      IF (iPart2.EQ.iPart) CYCLE      
      vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart2)) &
              * (PartState(PartStatevMPFSpec(iPart2),iDir)**2)
      vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
              * (PartState(PartStatevMPFSpec(iPart),iDir)**2)

      CALL RANDOM_NUMBER(ran1) 
      v_merge = (PartState(PartStatevMPFSpec(iPart2),iDir) - PartState(PartStatevMPFSpec(iPart), iDir))
      PartState(PartStatevMPFSpec(iPart2),iDir)=PartState(PartStatevMPFSpec(iPart),iDir) + v_merge*ran1
      PartState(PartStatevMPFSpec(iPart), iDir) =PartState(PartStatevMPFSpec(iPart),iDir) + v_merge*(1.0-ran1)

      vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart2)) &
              * (PartState(PartStatevMPFSpec(iPart2),iDir)**2)
      vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
              * (PartState(PartStatevMPFSpec(iPart),iDir)**2)
      iLoop= iLoop + 1
      TempPartVelo(iPart,iDir-3)=PartState(PartStatevMPFSpec(iPart), iDir)
      TempPartVelo(iPart2,iDir-3)=PartState(PartStatevMPFSpec(iPart2), iDir)
    END DO
    WRITE(*,*)'Loops for energy transformation needed: ', iLoop
  END IF

  DEALLOCATE(PartTemp)

END SUBROUTINE SetNewTemp_2

SUBROUTINE SetNewDistrVelo(NewPartNum, nDist, SpecNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, Species, PartSpecies, vMPF_oldEngSum, vMPF_oldMomSum, &
                     PartMPF, PartStatevMPFSpec, vMPFOldBrownVelo, vMPFOldMPF,vMPF_oldMPFSum
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: NewPartNum, nDist, SpecNum
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  REAL                            :: v_min, v_max, v_width, iRan, iRan2, v_merge, ran1
  INTEGER                         :: iDir, iPart, iBar, SpecID, iLoop, iPart2
  REAL, ALLOCATABLE            :: numDist(:)
  REAL                      :: TempPartVelo(NewPartNum,3)
!===================================================================================================================================

  SpecID = PartSpecies(PartStatevMPFSpec(1))
  ALLOCATE(numDist(nDist))
  DO iDir = 1, 3
    v_min = MINVAL(vMPFOldBrownVelo(:, iDir))
    v_max = MAXVAL(vMPFOldBrownVelo(:, iDir))
    v_width = (v_max - v_min)/nDist
    numDist = 0
    DO iPart = 1, SpecNum
      iBar = MIN(INT((vMPFOldBrownVelo(iPart, iDir)-v_min)/v_width+1), nDist)
      numDist(iBar) = numDist(iBar) + vMPFOldMPF(iPart)
    END DO

    numDist = numDist / vMPF_oldMPFSum

    
    print*, 'Verteilen, Dir: ', iDir
    DO iPart =1, NewPartNum -1
      IF (iDir.EQ.1) THEN 
        vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
              * (PartState(PartStatevMPFSpec(iPart),4)**2 + PartState(PartStatevMPFSpec(iPart),5)**2 &
              + PartState(PartStatevMPFSpec(iPart),6)**2)
        vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) + Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
                         * PartState(PartStatevMPFSpec(iPart),4:6)   
      END IF
      
      CALL RANDOM_NUMBER(iRan)  
      iBar = INT(iRan*nDist + 1)    
      CALL RANDOM_NUMBER(iRan2)
      DO WHILE (iRan2.GE.numDist(iBar))
        CALL RANDOM_NUMBER(iRan)  
        iBar = INT(iRan*nDist + 1)    
        CALL RANDOM_NUMBER(iRan2)
      END DO

      CALL RANDOM_NUMBER(iRan)
      PartState(PartStatevMPFSpec(iPart),iDir+3) = PartState(PartStatevMPFSpec(iPart),iDir+3) &
                          + (v_min + v_width*(iBar-1) + v_width*iRan)
    END DO
  END DO

  DO iPart=1, NewPartNum -1
    vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
            * (PartState(PartStatevMPFSpec(iPart),4)**2 + PartState(PartStatevMPFSpec(iPart),5)**2 &
            + PartState(PartStatevMPFSpec(iPart),6)**2)
    vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
                         * PartState(PartStatevMPFSpec(iPart),4:6)    
  END DO

  vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
          * (PartState(PartStatevMPFSpec(NewPartNum),4)**2 + PartState(PartStatevMPFSpec(NewPartNum),5)**2 &
          + PartState(PartStatevMPFSpec(NewPartNum),6)**2)
  vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) + Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
                         * PartState(PartStatevMPFSpec(NewPartNum),4:6)  

  PartState(PartStatevMPFSpec(NewPartNum),4:6) =vMPF_oldMomSum(1:3) &
                        / (Species(SpecID)%MassIC*PartMPF(PartStatevMPFSpec(NewPartNum)) )
  vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(NewPartNum)) &
          * (PartState(PartStatevMPFSpec(NewPartNum),4)**2 + PartState(PartStatevMPFSpec(NewPartNum),5)**2 &
          + PartState(PartStatevMPFSpec(NewPartNum),6)**2)


  IF (vMPF_oldEngSum.LT.0.0) THEN
    DO iPart = 1, NewPartNum  
      TempPartVelo(iPart,1:3) = PartState(PartStatevMPFSpec(iPart),4:6)
    END DO

    iLoop = 0
     DO WHILE (vMPF_oldEngSum.LT.0.0)
      CALL RANDOM_NUMBER(ran1)    
      iDir = INT(3*ran1 + 1) + 3
      iPart2 = MAXLOC(TempPartVelo(:,iDir-3),1)
      iPart = MINLOC(TempPartVelo(:,iDir-3),1)
      IF (iPart2.EQ.iPart) CYCLE      
      vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart2)) &
              * (PartState(PartStatevMPFSpec(iPart2),iDir)**2)
      vMPF_oldEngSum = vMPF_oldEngSum + 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
              * (PartState(PartStatevMPFSpec(iPart),iDir)**2)

      CALL RANDOM_NUMBER(ran1) 
      v_merge = (PartState(PartStatevMPFSpec(iPart2),iDir) - PartState(PartStatevMPFSpec(iPart), iDir))
      PartState(PartStatevMPFSpec(iPart2),iDir)=PartState(PartStatevMPFSpec(iPart),iDir) + v_merge*ran1
      PartState(PartStatevMPFSpec(iPart), iDir) =PartState(PartStatevMPFSpec(iPart),iDir) + v_merge*(1.0-ran1)

      vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart2)) &
              * (PartState(PartStatevMPFSpec(iPart2),iDir)**2)
      vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(SpecID)%MassIC * PartMPF(PartStatevMPFSpec(iPart)) &
              * (PartState(PartStatevMPFSpec(iPart),iDir)**2)
      iLoop= iLoop + 1
      TempPartVelo(iPart,iDir-3)=PartState(PartStatevMPFSpec(iPart), iDir)
      TempPartVelo(iPart2,iDir-3)=PartState(PartStatevMPFSpec(iPart2), iDir)
    END DO
    WRITE(*,*)'Loops for energy transformation needed: ', iLoop
  END IF

END SUBROUTINE SetNewDistrVelo

END MODULE MOD_part_MPFtools
