MODULE MOD_part_MPFtools                                                                         !
!===================================================================================================================================
! CONTAINS THE vMPF part
!===================================================================================================================================
  IMPLICIT NONE 
  PRIVATE                                                                                          !
!===================================================================================================================================

!===================================================================================================================================
! PUBLIC 
PUBLIC :: SplitParticle, MergeParticles, DefineElemT_inv, DefinePolyVec, DefineSplitVec, StartParticleMerge, GeoCoordToMap &
        , MapToGeo
!===================================================================================================================================

CONTAINS   

SUBROUTINE StartParticleMerge()                                                                
!===================================================================================================================================
! Particle Merge routine
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY :doParticleMerge, nSpecies,vMPFMergeParticleTarget,vMPF_SpecNumElem,vMPFSplitParticleTarget
  USE MOD_Mesh_Vars,     ONLY : nElems
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                      
!----------------------------------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE StartParticleMerge
                                                                                        
                                                                                                   
SUBROUTINE SplitParticle(iPart, deltaE)                                                                
!===================================================================================================================================
! Split Particles
!===================================================================================================================================
  USE MOD_Particle_Vars, ONLY : PDM, PartState, RandomVec, NumRanVec, PartSpecies, PartMPF, PEM, Species  
  USE MOD_DSMC_Vars, ONLY : useDSMC, CollisMode, PartStateIntEn                                                      
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER,INTENT(IN)              :: iPart
  REAL, INTENT(IN)                :: deltaE         
!----------------------------------------------------------------------------------------------------------------------------------
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
    PRINT*, 'New Particle Number greater max Part Num'
    STOP
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
  USE MOD_Particle_Vars, ONLY :  vMPFMergePolyOrder , PDM, PartState, PEM, PartStateMap, vMPFPolyPoint &
                            , vMPF_OrderVec, vMPFPolySol, vMPFMergeCellSplitOrder, vMPF_SplitVec, vMPFOldVelo &
                            , PartStatevMPFSpec, vMPFOldPos, vMPFOldMPF, PartSpecies, PartMPF, Species
  USE MOD_DSMC_Vars, ONLY :  CollisMode                                                     
  USE Levenberg_Marquardt
  USE MOD_PICInterpolation_Vars, ONLY :ElemT_inv
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                      
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        
  INTEGER,INTENT(IN)              :: iElem, SpecNum, SpecID
  INTEGER,INTENT(IN)              :: NumFinPart
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION                                                                       
  REAL                  :: PosMapped(3), CellTemp, test_vor(3)
  INTEGER               :: iPart, iLoop, iLoop2, PositionNbr
!===================================================================================================================================

!test_vor = 0.0

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
    CALL GeoCoordToMap(PartState(iPart,1:3), PartStateMap(iLoop2,1:3), iElem)
    PartStatevMPFSpec(iLoop2) = iPart
    iLoop2 = iLoop2 + 1
!    test_vor(1:3) = test_vor(1:3) + PartState(iPart,4:6)*PartMPF(iPart)*Species(PartSpecies(iPart))%MassIC
  END IF
  iPart = PEM%pNext(iPart)    
END DO

WRITE(*,*) 'Start Particle Split/Merge'
IF (NumFinPart.GT.SpecNum) THEN
  DO iLoop = 1 , NumFinPart - SpecNum
    PDM%ParticleVecLength = PDM%ParticleVecLength + 1
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + 1 
    PositionNbr = PDM%nextFreePosition(PDM%CurrentNextFreePosition)
    IF (PositionNbr.EQ.0) THEN
      PRINT*, 'New Particle Number greater max Part Num'
      STOP
    END IF
    PartStatevMPFSpec(SpecNum + iLoop) = PositionNbr
    !Set new particle parameters
    PDM%ParticleInside(PositionNbr) = .true.
    PartSpecies(PositionNbr) = SpecID
    PEM%Element(PositionNbr) = iElem
  END DO
END IF


CALL SplitRegion(iElem, SpecNum)
CALL DeleteParticlesMPF(iElem, NumFinPart, CellTemp, SpecNum, SpecID)
!CALL SetMPFParticlePos(iElem, NumFinPart, x)
CALL SetMPFParticlePosCube(iElem, NumFinPart,SpecNum)
CALL SetNewvMPF(NumFinPart, iElem)
CALL SetNewVelos(NumFinPart, iElem, CellTemp, SpecNum, SpecID)
WRITE(*,*) 'Finish Particle Split/Merge'

DEALLOCATE(PartStateMap, PartStatevMPFSpec, vMPFOldVelo, vMPFOldPos, vMPFOldMPF)

END SUBROUTINE MergeParticles


SUBROUTINE fcn (m, n, x, fvec, fjac, iflag)
!===================================================================================================================================
! do not know
!===================================================================================================================================
! IMPLICIT VARIABLES HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag
!===================================================================================================================================
IF (iflag == 1) CALL ssqfcn (m, n, x, fvec)
IF (iflag == 2) CALL ssqjac (m, n, x, fjac)
RETURN
END SUBROUTINE fcn


SUBROUTINE ssqjac (m, n, x, fjac)
!===================================================================================================================================
! ssqjqc
!===================================================================================================================================
! MODULES
USE Levenberg_Marquardt
USE MOD_Particle_Vars,    ONLY:vMPFPolyPoint, vMPF_OrderVec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER, INTENT(IN)     :: m, n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: fjac(:,:)
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
!===================================================================================================================================
! MODULES
USE Levenberg_Marquardt
USE MOD_Particle_Vars,    ONLY:vMPFPolyPoint, vMPF_OrderVec, vMPFPolySol,vMPFOldMPF
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER, INTENT(IN)     :: m, n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: fvec(:)
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

SUBROUTINE GeoCoordToMap(x_in,xi_Out,iElem)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions x
! first get xi,eta,zeta from x,y,z...then do tenso product interpolation
! xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
USE MOD_Basis,ONLY: LagrangeInterpolationPolys
USE MOD_Interpolation_Vars,ONLY: wBary,xGP
USE MOD_Particle_Vars,ONLY:GEO
USE MOD_PICInterpolation_Vars, ONLY: ElemT_inv
USE MOD_Eval_xyz, ONLY: Calc_F, Calc_dF_inv 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: iElem                                 ! elem index
REAL,INTENT(IN)     :: x_in(3)                                  ! physical position of particle 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: xi_Out(3)  ! Interpolated Pos
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out,i,j,k, iNode
REAL                :: xi(3)     
REAL              :: P(3,8), F(3), dF_inv(3,3), s(3) 
REAL, PARAMETER   :: EPS=1E-8
REAL              :: T_inv(3,3), DP(3)
!===================================================================================================================================
! --------------------------------------------------
! 1.) Mapping: get xi,eta,zeta value from x,y,z
! --------------------------------------------------
! 1.1.) initial guess from linear part:
DO iNode = 1,8
  P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
END DO
T_inv(1:3,1:3) = ElemT_inv(1:3,1:3,iElem)

! transform also the physical coordinate of the point 
! into the unit element (this is the solution of the 
! linear problem already 
xi = 0.
DP = x_in - P(:,1)
DO i=1,3
  DO j=1,3
    xi(i)= xi(i) + T_inv(i,j) * DP(j) 
  END DO
END DO
xi = xi - (/1.,1.,1./)
! 1.2.) Newton-Method to solve non-linear part
!       If linear elements then F should becom 0 and no 
!       Newton step is required.

F = Calc_F(xi,x_in,P)
DO WHILE(SUM(ABS(F)).GE.EPS) 
  dF_inv = Calc_dF_inv(xi,P)
  s=0.
  DO j = 1,3
    DO k = 1,3
      s(j) = s(j) + dF_inv(j,k) * F(k) 
    END DO ! k
  END DO ! j
  xi = xi - s
  F = Calc_F(xi,x_in,P)
END DO ! i
xi_Out = xi

END SUBROUTINE GeoCoordToMap


SUBROUTINE DefineElemT_inv()                                                                !
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : GEO                                                     !
  USE MOD_PICInterpolation_Vars, ONLY :ElemT_inv,  InterpolationType 
  USE MOD_Mesh_Vars,              ONLY : nElems
  USE MOD_PICInterpolation, ONLY : Calc_inv
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: iElem, iNode
  REAL                      :: P(3,8), T(3,3), T_inv(3,3)
!===================================================================================================================================
IF(TRIM(InterpolationType).NE.'particle_position') THEN
  ALLOCATE(ElemT_inv(1:3,1:3,1:nElems))
  
  DO iElem = 1,nElems
    DO iNode = 1,8
      P(1:3,iNode) = GEO%NodeCoords(1:3,GEO%ElemToNodeID(iNode,iElem))
    END DO
    T(:,1) = 0.5 * (P(:,2)-P(:,1))
    T(:,2) = 0.5 * (P(:,4)-P(:,1))
    T(:,3) = 0.5 * (P(:,5)-P(:,1))
    T_inv = Calc_inv(T)
    ElemT_inv(1:3,1:3,iElem) = T_inv
  END DO
END IF
END SUBROUTINE DefineElemT_inv

SUBROUTINE SplitRegion(iElem, SpecNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : vMPFMergeCellSplitOrder,PartStateMap, PEM, vMPF_SplitVec, vMPFPolyPoint &
                                ,vMPFPolySol, PartMPF, vMPF_oldMPFSum, vMPF_SplitVecBack, PartStatevMPFSpec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!--------------------------------------------------------------------------------------------------!
! ARGUMENT LIST DECLARATION                                                                        !
INTEGER, INTENT(IN)   :: iElem, SpecNum
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
REAL, ALLOCATABLE     :: RegPartNum(:)
INTEGER               :: iPart, iLoop, PolOrder, iPart_indx, x_cube, y_cube, z_cube
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


SUBROUTINE DeleteParticlesMPF(iElem, FinPartNum, Temp, SpecNum, SpecID)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, vMPF_oldEngSum, vMPF_oldMomSum , PEM, Species, PartMPF, PDM &
                          , PartSpecies, BoltzmannConst, PartStatevMPFSpec, vMPFOldPos, vMPFOldVelo, vMPFOldMPF
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
INTEGER, INTENT(IN)   :: iElem, FinPartNum, SpecNum, SpecID
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
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
END SUBROUTINE DeleteParticlesMPF


SUBROUTINE SetMPFParticlePos(iElem, FinPartNum,x)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! Modules
  USE MOD_Particle_Vars, ONLY : PartState, PEM, vMPFMergePolyOrder, vMPF_OrderVec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
INTEGER, INTENT(IN)   :: iElem, FinPartNum
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
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
END SUBROUTINE SetMPFParticlePos


SUBROUTINE SetMPFParticlePosCube(iElem, FinPartNum, SpecNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, PEM, vMPFMergeCellSplitOrder, vMPF_OrderVec &
                        , vMPFPolySol, vMPF_SplitVecBack, GEO, PartStatevMPFSpec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
INTEGER, INTENT(IN)   :: iElem, FinPartNum, SpecNum
!----------------------------------------------------------------------------------------------------------------------------------
! Local variable declaration
INTEGER               :: iPart, iLoop, iDOF, DOF_LMInput, x_cube, y_cube, z_cube, iNode, iLoop2, PositionNbr
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
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
END SUBROUTINE SetMPFParticlePosCube


SUBROUTINE SetNewVelos(NewPartNum, iElem, Temp, SpecNum, SpecID)                                                                
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : vMPFOldVelo, PEM, vMPFOldPos, PartState, vMPFMergePolyOrder &
                        ,vMPFPolyPoint,vMPFPolySol, vMPF_OrderVec, vMPF_oldEngSum, vMPF_oldMomSum & 
                        ,Species, PartSpecies, PartMPF, PartStatevMPFSpec, PDM, vMPFOldMPF, vMPFOldBrownVelo &  
                        , vMPF_velocityDistribution
  USE Levenberg_Marquardt
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: NewPartNum, iElem, SpecNum, SpecID
  REAL,INTENT(IN)                 :: Temp
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: iPart, iLoop, iDir, DOF_LMInput, info, iPart_indx
  INTEGER, ALLOCATABLE            :: iwa(:)
  DOUBLE PRECISION, ALLOCATABLE   :: fjac(:,:), fvec(:), x(:) !!evtl double
  DOUBLE PRECISION                :: tol
  REAL                            :: iRan, test_nach(3), test
  integer       :: k,l,m
!===================================================================================================================================

test_nach = 0.0

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
!--------------------------------------------------------------------------------------------------!
END SUBROUTINE SetNewVelos
!--------------------------------------------------------------------------------------------------!


FUNCTION MapToGeo(xi,P)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: xi(3)      ! 
REAL,INTENT(IN)          :: P(3,8)     ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: MapToGeo(3)  !  
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
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


SUBROUTINE SetNewTemp(PartIndx, Temp, iPart)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, Species, PartSpecies, vMPF_oldEngSum, vMPF_oldMomSum, PartMPF
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
            + ran1*SQRT(-2*1.3806488E-23*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)
PartState(PartIndx, 5) = PartState(PartIndx, 5) &
            + ran2*SQRT(-2*1.3806488E-23*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)

SumRan = 2
DO WHILE (SumRan .GT. 1)
 CALL RANDOM_NUMBER(RandVal)
 ran1 = 2*RandVal(1) - 1
 ran2 = 2*RandVal(2) - 1
 SumRan = ran1**2 + ran2**2
END DO
PartState(PartIndx, 6) =PartState(PartIndx, 6) &
            + ran1*SQRT(-2*1.3806488E-23*Temp/Species(PartSpecies(PartIndx))%MassIC*LOG(SumRan)/SumRan)

vMPF_oldEngSum = vMPF_oldEngSum - 0.5 * Species(PartSpecies(PartIndx))%MassIC * PartMPF(PartIndx) &
        * (PartState(PartIndx,4)**2 + PartState(PartIndx,5)**2 &
        + PartState(PartIndx,6)**2)
vMPF_oldMomSum(1:3) = vMPF_oldMomSum(1:3) - Species(PartSpecies(PartIndx))%MassIC * PartMPF(PartIndx) &
                     * PartState(PartIndx,4:6)
IF(vMPF_oldEngSum.lT.0) then
 print*, 'mist: ', iPart
  read*
end if

!--------------------------------------------------------------------------------------------------!
END SUBROUTINE SetNewTemp

!--------------------------------------------------------------------------------------------------!

SUBROUTINE SetNewvMPF(FinPartNum, iElem)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartMPF, vMPF_oldMPFSum, PEM, PartStatevMPFSpec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE                                                                                    !
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION                                                                        !
  INTEGER,INTENT(IN)              :: FinPartNum, iElem
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                         :: iPart, iLoop
  REAL                            :: NewMPF
!===================================================================================================================================

NewMPF = vMPF_oldMPFSum/FinPartNum
DO iLoop = 1, FinPartNum  
  PartMPF(PartStatevMPFSpec(iLoop)) = NewMPF
END DO

!--------------------------------------------------------------------------------------------------!
END SUBROUTINE SetNewvMPF


SUBROUTINE SetNewTemp_2(Temp, NewPartNum)                                                                !
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Particle_Vars, ONLY : PartState, Species, PartSpecies, vMPF_oldEngSum, vMPF_oldMomSum, PartMPF, PartStatevMPFSpec
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
maxwellfac = SQRT(3. * 1.380650524E-23 * Temp/ &              ! velocity of maximum
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

!--------------------------------------------------------------------------------------------------!
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
  REAL                            :: v_min, v_max, v_width, v_dist, iRan, iRan2, v_merge, ran1, ran2
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
