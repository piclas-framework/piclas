#include "boltzplatz.h"

MODULE MOD_DSMC_ParticlePairing
!===================================================================================================================================
! Module including different particle pairing algorithms
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_pairing_octree
  MODULE PROCEDURE DSMC_pairing_octree
END INTERFACE

INTERFACE DSMC_pairing_statistical
  MODULE PROCEDURE DSMC_pairing_statistical
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_pairing_octree, DSMC_pairing_statistical
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_pairing_octree(iElem)
!===================================================================================================================================
! Pairing subroutine for octree and nearest neighbour, decides whether to create a new octree node or start nearest neighbour search
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : tTreeNode, DSMC
  USE MOD_Particle_Vars,          ONLY : PEM, GEO
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, nPart, iNode
  REAL                          :: ApproxElemMid(1:3)
  TYPE(tTreeNode), POINTER      :: TreeNode
!===================================================================================================================================

  nPart = PEM%pNumber(iElem)

  NULLIFY(TreeNode)
  ApproxElemMid(1:3)= 0.0

  ALLOCATE(TreeNode)
  ALLOCATE(TreeNode%iPartIndx_Node(nPart)) ! List of particles in the cell neccessary for stat pairing
  TreeNode%iPartIndx_Node(1:nPart) = 0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    TreeNode%iPartIndx_Node(iLoop) = iPart
    iPart = PEM%pNext(iPart)    
  END DO

  !Set approx element mid point Sollte man in die ini packen!!!!!
  IF(nPart.GT.DSMC%PartNumOctreeNode) THEN
    DO iNode = 1, 8
      TreeNode%MidPoint(1:3) = TreeNode%MidPoint(1:3) + GEO%NodeCoords(1:3, GEO%ElemToNodeID(iNode,iElem))          
    END DO
    TreeNode%MidPoint(1:3) = TreeNode%MidPoint(1:3) / 8.0
    TreeNode%PNum_Node = nPart
    CALL AddOctreeNode(TreeNode, iElem)
  ELSE  IF (nPart.GT.1) THEN
    CALL FindNearestNeigh(TreeNode%iPartIndx_Node, nPart &
                              , iElem, GEO%Volume(iElem))
  END IF

  DEALLOCATE(TreeNode%iPartIndx_Node) 
  DEALLOCATE(TreeNode)

END SUBROUTINE DSMC_pairing_octree


RECURSIVE SUBROUTINE AddOctreeNode(TreeNode, iElem)
!===================================================================================================================================
! Adds additional octree node/branch
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : tTreeNode, DSMC
  USE MOD_Particle_Vars,          ONLY : PartState
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)                     :: iElem
  TYPE(tTreeNode),INTENT(IN), POINTER     :: TreeNode
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPart, iLoop, iPartIndx
  INTEGER, ALLOCATABLE         :: iPartIndx_ChildNode(:,:)
  REAL                          :: ChildMidPoints(3,8), NodeVolume(8)
  INTEGER                       :: PartNumChildNode(8), PairNumChildNode(8)
!===================================================================================================================================

  ALLOCATE(iPartIndx_ChildNode(8,TreeNode%PNum_Node))
  PartNumChildNode(:) = 0
  DO iLoop = 1, 8
    ChildMidPoints(1:3, iLoop) = TreeNode%MidPoint(1:3)
  END DO

  !         Numbering of the 8 ChildNodes of the Octree
  !          __________
  !         /    /    /|   |z
  !        /  7 / 6  / |   |
  !       /----/----/|6|   |
  !      /_3__/_2__/ |/|   |______Y
  !     |    |    |2 |5|   /
  !     | 3  | 2  |/ |/   /
  !     |----|----|1 /   /x
  !     | 4  | 1  | /
  !     |____|____|/

  ! particle to Octree ChildNode sorting
  DO iPart=1,TreeNode%PNum_Node
    iPartIndx = TreeNode%iPartIndx_Node(iPart)
    IF ((PartState(iPartIndx,1).GE.TreeNode%MidPoint(1)).AND.(PartState(iPartIndx,2).GE.TreeNode%MidPoint(2)) &
        .AND.(PartState(iPartIndx,3).LE.TreeNode%MidPoint(3))) THEN
      ! Locate new adapted mid point for child node
      IF (PartState(iPartIndx,1).GT.ChildMidPoints(1, 1)) ChildMidPoints(1, 1) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).GT.ChildMidPoints(2, 1)) ChildMidPoints(2, 1) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).LT.ChildMidPoints(3, 1)) ChildMidPoints(3, 1) = PartState(iPartIndx,3)
      PartNumChildNode(1) = PartNumChildNode(1) + 1
      iPartIndx_ChildNode(1,PartNumChildNode(1)) = iPartIndx
    ELSE IF((PartState(iPartIndx,1).GE.TreeNode%MidPoint(1)).AND.(PartState(iPartIndx,2).GE.TreeNode%MidPoint(2))) THEN
      IF (PartState(iPartIndx,1).GT.ChildMidPoints(1, 2)) ChildMidPoints(1, 2) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).GT.ChildMidPoints(2, 2)) ChildMidPoints(2, 2) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).GT.ChildMidPoints(3, 2)) ChildMidPoints(3, 2) = PartState(iPartIndx,3)
      PartNumChildNode(2) = PartNumChildNode(2) + 1
      iPartIndx_ChildNode(2,PartNumChildNode(2)) = iPartIndx
    ELSE IF((PartState(iPartIndx,1).GE.TreeNode%MidPoint(1)).AND.(PartState(iPartIndx,3).GE.TreeNode%MidPoint(3))) THEN
      IF (PartState(iPartIndx,1).GT.ChildMidPoints(1, 3)) ChildMidPoints(1, 3) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).LT.ChildMidPoints(2, 3)) ChildMidPoints(2, 3) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).GT.ChildMidPoints(3, 3)) ChildMidPoints(3, 3) = PartState(iPartIndx,3)
      PartNumChildNode(3) = PartNumChildNode(3) + 1
      iPartIndx_ChildNode(3,PartNumChildNode(3)) = iPartIndx
    ELSE IF (PartState(iPartIndx,1).GE.TreeNode%MidPoint(1)) THEN
      IF (PartState(iPartIndx,1).GT.ChildMidPoints(1, 4)) ChildMidPoints(1, 4) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).LT.ChildMidPoints(2, 4)) ChildMidPoints(2, 4) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).LT.ChildMidPoints(3, 4)) ChildMidPoints(3, 4) = PartState(iPartIndx,3)
      PartNumChildNode(4) = PartNumChildNode(4) + 1
      iPartIndx_ChildNode(4,PartNumChildNode(4)) = iPartIndx
    ELSE IF((PartState(iPartIndx,2).GE.TreeNode%MidPoint(2)).AND.(PartState(iPartIndx,3).LE.TreeNode%MidPoint(3))) THEN
      IF (PartState(iPartIndx,1).LT.ChildMidPoints(1, 5)) ChildMidPoints(1, 5) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).GT.ChildMidPoints(2, 5)) ChildMidPoints(2, 5) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).LT.ChildMidPoints(3, 5)) ChildMidPoints(3, 5) = PartState(iPartIndx,3)
      PartNumChildNode(5) = PartNumChildNode(5) + 1
      iPartIndx_ChildNode(5,PartNumChildNode(5)) = iPartIndx
    ELSE IF (PartState(iPartIndx,2).GE.TreeNode%MidPoint(2)) THEN
      IF (PartState(iPartIndx,1).LT.ChildMidPoints(1, 6)) ChildMidPoints(1, 6) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).GT.ChildMidPoints(2, 6)) ChildMidPoints(2, 6) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).GT.ChildMidPoints(3, 6)) ChildMidPoints(3, 6) = PartState(iPartIndx,3)
      PartNumChildNode(6) = PartNumChildNode(6) + 1
      iPartIndx_ChildNode(6,PartNumChildNode(6)) = iPartIndx
    ELSE IF (PartState(iPartIndx,3).GE.TreeNode%MidPoint(3)) THEN
      IF (PartState(iPartIndx,1).LT.ChildMidPoints(1, 7)) ChildMidPoints(1, 7) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).LT.ChildMidPoints(2, 7)) ChildMidPoints(2, 7) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).GT.ChildMidPoints(3, 7)) ChildMidPoints(3, 7) = PartState(iPartIndx,3)
      PartNumChildNode(7) = PartNumChildNode(7) + 1
      iPartIndx_ChildNode(7,PartNumChildNode(7)) = iPartIndx
    ELSE 
      IF (PartState(iPartIndx,1).LT.ChildMidPoints(1, 8)) ChildMidPoints(1, 8) = PartState(iPartIndx,1)
      IF (PartState(iPartIndx,2).LT.ChildMidPoints(2, 8)) ChildMidPoints(2, 8) = PartState(iPartIndx,2)
      IF (PartState(iPartIndx,3).LT.ChildMidPoints(3, 8)) ChildMidPoints(3, 8) = PartState(iPartIndx,3)
      PartNumChildNode(8) = PartNumChildNode(8) + 1
      iPartIndx_ChildNode(8,PartNumChildNode(8)) = iPartIndx
    END IF
  END DO

  ! calculate new pair number of child nodes
  DO iLoop=1, 8
    NodeVolume(iLoop) = ABS(TreeNode%MidPoint(1) - ChildMidPoints(1, iLoop)) &
                      * ABS(TreeNode%MidPoint(2) - ChildMidPoints(2, iLoop)) &
                      * ABS(TreeNode%MidPoint(3) - ChildMidPoints(3, iLoop)) 
    ChildMidPoints(1:3, iLoop) = TreeNode%MidPoint(1:3) - 0.5 * (TreeNode%MidPoint(1:3) - ChildMidPoints(1:3, iLoop))
  END DO

  DO iLoop = 1, 8
    IF(PartNumChildNode(iLoop).GT.DSMC%PartNumOctreeNode) THEN
      NULLIFY(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode)
      ALLOCATE(TreeNode%ChildNode%iPartIndx_Node(PartNumChildNode(iLoop)))
      TreeNode%ChildNode%MidPoint(1:3) = ChildMidPoints(1:3, iLoop)
      TreeNode%ChildNode%iPartIndx_Node(1:PartNumChildNode(iLoop)) = iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop))
      TreeNode%ChildNode%PNum_Node = PartNumChildNode(iLoop)
      TreeNode%ChildNode%PairNum_Node = PairNumChildNode(iLoop)    
      CALL AddOctreeNode(TreeNode%ChildNode, iElem)
      DEALLOCATE(TreeNode%ChildNode%iPartIndx_Node)
      DEALLOCATE(TreeNode%ChildNode)
    ELSE IF (PartNumChildNode(iLoop).GT.1) THEN

      CALL FindNearestNeigh(iPartIndx_ChildNode(iLoop, 1:PartNumChildNode(iLoop)), PartNumChildNode(iLoop) &
                            , iElem, NodeVolume(iLoop))
    END IF
  END DO

END SUBROUTINE AddOctreeNode


SUBROUTINE FindNearestNeigh(iPartIndx_Node, PartNum, iElem, NodeVolume)
!===================================================================================================================================
! Finds nearest neighbour for collision pairing
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : CollInf, tTreeNode, CollisMode, ChemReac, PartStateIntEn, Coll_pData
  USE MOD_DSMC_Vars,              ONLY : DSMC, PairE_vMPF
  USE MOD_Particle_Vars,          ONLY : PartState, nSpecies, PartSpecies, usevMPF, PartMPF
  USE MOD_DSMC_ChemReact,         ONLY : SetMeanVibQua   
  USE MOD_Particle_Analyze_Vars,  ONLY : CalcEkin
  USE MOD_DSMC_CollisionProb,     ONLY : DSMC_prob_calc
  USE MOD_DSMC_Collis,            ONLY : DSMC_perform_collision
  USE MOD_vmpf_collision,         ONLY : DSMC_vmpf_prob
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)                        :: NodeVolume
  INTEGER, INTENT(IN)                     :: PartNum
  INTEGER, INTENT(IN)                     :: iElem
  INTEGER, INTENT(INOUT)                  :: iPartIndx_Node(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: iPair, iPart1, iPart2, iLoop, iPart, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase , PairNum_Node
  REAL                          :: Dist1, Dist2, iRan
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

  PairNum_Node = INT(PartNum/2)
  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
  END IF

  IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) THEN
    ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0
    DO iPart = 1, PartNum
      ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) = &
                            ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) &
                          + PartStateIntEn(iPart,1)
    END DO
  END IF

  DO iPart = 1, PartNum
    CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) = &
              CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx_Node(iPart))) + 1 
  END DO

  IF (CollisMode.EQ.3) THEN
    ChemReac%nPartForRec = PartNum
    ChemReac%nPairForRec = PairNum_Node
  END IF
  ALLOCATE(Coll_pData(PairNum_Node))
  nPart = PartNum

  IF (usevMPF) MPFFac = 1
  DO iPair = 1, PairNum_Node
    CALL RANDOM_NUMBER(iRan) 
    iPart1 = 1 + INT(nPart * iRan)
    Coll_pData(iPair)%iPart_p1 = iPartIndx_Node(iPart1)
    iPartIndx_Node(iPart1) = iPartIndx_Node(nPart)
    nPart = nPart - 1
    iPart2 = 1
    Dist1 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
           - PartState(iPartIndx_Node(iPart2),1))**2 &
           +(PartState(Coll_pData(iPair)%iPart_p1,2) &
           - PartState(iPartIndx_Node(iPart2),2))**2 &
           +(PartState(Coll_pData(iPair)%iPart_p1,3) &
           - PartState(iPartIndx_Node(iPart2),3))**2 
    DO iLoop = 2, nPart
      Dist2 = (PartState(Coll_pData(iPair)%iPart_p1,1) &
             - PartState(iPartIndx_Node(iLoop),1))**2 &
             +(PartState(Coll_pData(iPair)%iPart_p1,2) &
             - PartState(iPartIndx_Node(iLoop),2))**2 &
             +(PartState(Coll_pData(iPair)%iPart_p1,3) &
             - PartState(iPartIndx_Node(iLoop),3))**2 
      IF (Dist2.LT.Dist1) THEN
        iPart2 = iLoop
        Dist1 = Dist2
      END IF
    END DO
    Coll_pData(iPair)%iPart_p2 = iPartIndx_Node(iPart2)
    iPartIndx_Node(iPart2) = iPartIndx_Node(nPart)
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1   
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2
    
    IF (usevMPF) THEN
      TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
      IF (TempMPFFac .GE. MPFFac) THEN
          MPFFac = TempMPFFac
          PairE_vMPF(1) = iPair
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
        ELSE
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
        END IF
      END IF
    END IF
    iCase = CollInf%Coll_Case(cSpec1, cSpec2) 
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2 
    Coll_pData(iPair)%PairType = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
  END DO
  IF ((PairNum_Node.NE.0).AND.(CollisMode.EQ.3).and.(MOD(PartNum, PairNum_Node).NE.0)) THEN
    ChemReac%RecombParticle = iPartIndx_Node(1)
  END IF

  IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) THEN
    CALL SetMeanVibQua()
  END IF

  DO iPair = 1,  PairNum_Node
    IF(.NOT.Coll_pData(iPair)%NeedForRec) THEN  
      IF (usevMPF) THEN            ! calculation of collision prob
        CALL DSMC_vmpf_prob(iElem, iPair, NodeVolume)
      ELSE
        CALL DSMC_prob_calc(iElem, iPair, NodeVolume)
      END IF
      CALL RANDOM_NUMBER(iRan)
      IF (Coll_pData(iPair)%Prob.ge.iRan) THEN
#if (PP_TimeDiscMethod==42)
        IF(CalcEkin.OR.DSMC%ReservoirSimu) THEN
#else
        IF(CalcEkin) THEN
#endif
          DSMC%NumColl(Coll_pData(iPair)%PairType) = DSMC%NumColl(Coll_pData(iPair)%PairType) + 1
          DSMC%NumColl(CollInf%NumCase + 1) = DSMC%NumColl(CollInf%NumCase + 1) + 1
        END IF
        CALL DSMC_perform_collision(iPair, iElem, NodeVolume, PartNum)  ! call collision from octree
      END IF
    END IF
  END DO
  DEALLOCATE(Coll_pData)

END SUBROUTINE FindNearestNeigh


SUBROUTINE DSMC_pairing_statistical(iElem)
!===================================================================================================================================
! Classic statistical pairing method
!===================================================================================================================================
! MODULES
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, CollisMode, PartStateIntEn, ChemReac, PairE_vMPF
  USE MOD_Particle_Vars,          ONLY : PEM, PartSpecies, nSpecies, PartState, usevMPF, PartMPF
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  INTEGER, INTENT(IN)           :: iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER                       :: nPair, iPair, iPart, iLoop, cPart1, cPart2, nPart
  INTEGER                       :: cSpec1, cSpec2, iCase
  INTEGER, ALLOCATABLE         :: iPartIndx(:) ! List of particles in the cell nec for stat pairing
  REAL                          :: iRan
  REAL                          :: TempMPFFac, MPFFac
!===================================================================================================================================

  nPart = PEM%pNumber(iElem)
  nPair = INT(nPart/2)

  IF (CollisMode.EQ.3) THEN
    ChemReac%RecombParticle = 0
  END IF

  CollInf%Coll_SpecPartNum = 0
  CollInf%Coll_CaseNum = 0

  ALLOCATE(Coll_pData(nPair))
  ALLOCATE(iPartIndx(nPart))
  Coll_pData%Ec=0
  iPartIndx = 0
  IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) ChemReac%MeanEVib_PerIter(1:nSpecies) = 0.0

  iPart = PEM%pStart(iElem)                         ! create particle index list for pairing
  DO iLoop = 1, nPart
    iPartIndx(iLoop) = iPart
    CollInf%Coll_SpecPartNum(PartSpecies(iPart)) = CollInf%Coll_SpecPartNum(PartSpecies(iPart)) + 1 
    ! counter for part num of spec per cell
    IF ((CollisMode.EQ.3).AND.ChemReac%MeanEVib_Necc) ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) = &
                          ChemReac%MeanEVib_PerIter(PartSpecies(iPart)) &
                        + PartStateIntEn(iPart,1) !Calculation of mean evib per cell and iter, necessary for disso prob
    iPart = PEM%pNext(iPart)    
  END DO
  IF (CollisMode.EQ.3) THEN
    ChemReac%nPartForRec = nPart
    ChemReac%nPairForRec = nPair
  END IF

  IF (usevMPF) MPFFac = 1

  DO iPair = 1, nPair                               ! statistical pairing
    CALL RANDOM_NUMBER(iRan) 
    cPart1 = 1 + INT(nPart * iRan)                       ! first pair particle
    Coll_pData(iPair)%iPart_p1 = iPartIndx(cPart1)
    iPartIndx(cPart1) = iPartIndx(nPart)
    nPart = nPart - 1
    CALL RANDOM_NUMBER(iRan)
    cPart2 = 1 + INT(nPart * iRan)                       ! second pair particle
    Coll_pData(iPair)%iPart_p2 = iPartIndx(cPart2)
    iPartIndx(cPart2) = iPartIndx(nPart)
    nPart = nPart - 1

    cSpec1 = PartSpecies(Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(Coll_pData(iPair)%iPart_p2) !spec of particle 2

    IF (usevMPF) THEN
      TempMPFFac = PartMPF(Coll_pData(iPair)%iPart_p1) + PartMPF(Coll_pData(iPair)%iPart_p2)
      IF (TempMPFFac .GE. MPFFac) THEN
          MPFFac = TempMPFFac
          PairE_vMPF(1) = iPair
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p2
        ELSE
          PairE_vMPF(2) = Coll_pData(iPair)%iPart_p1
        END IF
      END IF
    END IF

    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)
    Coll_pData(iPair)%CRela2 = (PartState(Coll_pData(iPair)%iPart_p1,4) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,4))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,5) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,5))**2 &
                             + (PartState(Coll_pData(iPair)%iPart_p1,6) &
                             -  PartState(Coll_pData(iPair)%iPart_p2,6))**2
    Coll_pData(iPair)%PairType = iCase
    Coll_pData(iPair)%NeedForRec = .FALSE.
  END DO  
  IF ((nPair.NE.0).AND.(CollisMode.EQ.3).AND.(MOD(nPart, nPair).NE.0)) THEN
    ChemReac%RecombParticle = iPartIndx(1)
  END IF
  
  DEALLOCATE(iPartIndx)
END SUBROUTINE DSMC_pairing_statistical

END MODULE MOD_DSMC_ParticlePairing
