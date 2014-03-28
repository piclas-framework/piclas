#include "boltzplatz.h"

MODULE MOD_LD_Init
!===================================================================================================================================
! Initialisation of LD variables!
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

PUBLIC :: InitLD, CalcDegreeOfFreedom
!===================================================================================================================================

CONTAINS

SUBROUTINE InitLD()
!===================================================================================================================================
! Init of DSMC Vars
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, nSides, SideToElem, ElemToSide, nNodes, Elem_xGP
USE MOD_Particle_Vars,         ONLY : GEO, PDM, Species, PartSpecies, nSpecies
USE nr,                        ONLY : gaussj 
USE MOD_DSMC_Init,             ONLY : InitDSMC
USE MOD_DSMC_Vars,             ONLY : useDSMC, SpecDSMC, CollisMode
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iElem, trinum, iLocSide, iNode, iPart, iInit, iSpec
REAL                    :: IdentityMat(8,8), MatForGaussj(8,8)
REAL                    :: NVecTest
INTEGER,ALLOCATABLE     :: LocalNumberOfElemPerNode(:)
TYPE tElemIndxPerNode 
  INTEGER,ALLOCATABLE   :: NodeCellIndex(:)
END TYPE
TYPE(tElemIndxPerNode), ALLOCATABLE  :: ElemIndxPerNode(:)
CHARACTER(32)           :: hilf
!===================================================================================================================================

  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' LD INIT ...'

  LD_SecantMeth%Guess = GETREAL('LD-InitialGuess','10')
  LD_SecantMeth%MaxIter = GETINT('LD-MaxIterNumForLagVelo','100')
  LD_SecantMeth%Accuracy = GETREAL('LD-AccuracyForLagVelo','0.001')
  LD_RepositionFak = GETREAL('LD-RepositionsFaktor','0.0')
  LD_RelaxationFak = GETREAL('LD-RelaxationsFaktor','0.0')
  LD_CalcDelta_t=GETLOGICAL('LD-CFL-CalcDelta-t','.FALSE.')

  ALLOCATE(TempDens(nElems))

  ALLOCATE(BulkValues(nElems))
  BulkValues(1:nElems)%CellV(1)        = 0.0
  BulkValues(1:nElems)%CellV(2)        = 0.0
  BulkValues(1:nElems)%CellV(3)        = 0.0
  BulkValues(1:nElems)%DegreeOfFreedom = 0.0 
  BulkValues(1:nElems)%Beta            = 0.0  
  BulkValues(1:nElems)%MassDens        = 0.0 
!  BulkValues(1:nElems)%CellType        = 4

  ALLOCATE(BulkValuesOpenBC(nElems))
  BulkValuesOpenBC(1:nElems)%CellV(1)        = 0.0
  BulkValuesOpenBC(1:nElems)%CellV(2)        = 0.0
  BulkValuesOpenBC(1:nElems)%CellV(3)        = 0.0
  BulkValuesOpenBC(1:nElems)%DegreeOfFreedom = 0.0 
  BulkValuesOpenBC(1:nElems)%Beta            = 0.0  
  BulkValuesOpenBC(1:nElems)%MassDens        = 0.0

  ALLOCATE(SurfLagValues(6,nElems,2))
  SurfLagValues(1:6,1:nElems,1:2)%LagVelo    = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(1)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(2)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaM(3)  = 0.0
  SurfLagValues(1:6,1:nElems,1:2)%DeltaE     = 0.0

  ALLOCATE(MeanSurfValues(6,nElems))
  MeanSurfValues(1:6,1:nElems)%MeanLagVelo     = 0.0
  MeanSurfValues(1:6,1:nElems)%MeanBaseD       = 0.0
  MeanSurfValues(1:6,1:nElems)%MeanBaseD2      = 0.0
  MeanSurfValues(1:6,1:nElems)%MeanNormVec(1)  = 0.0
  MeanSurfValues(1:6,1:nElems)%MeanNormVec(2)  = 0.0
  MeanSurfValues(1:6,1:nElems)%MeanNormVec(3)  = 0.0

!!!!  ALLOCATE(NewNodePosIndx(1:3,nNodes))  !!! nur für "Tetra-Methode"

  ALLOCATE(IsDoneLagVelo(nSides))
  IsDoneLagVelo(1:nSides)   = .FALSE.
  DO iElem=1, nElems
    DO iLocSide = 1, 6
      DO trinum=1, 2
        SurfLagValues(iLocSide, iElem, trinum)%Area = CalcTriNumArea(iLocSide, iElem, trinum)
        CALL CalcLagNormVec(iLocSide, iElem, trinum)
        NVecTest = (GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(1,0,0,0,iElem)) &
                 * SurfLagValues(iLocSide, iElem,trinum)%LagNormVec(1) &
                 + (GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(2,0,0,0,iElem)) &
                 * SurfLagValues(iLocSide, iElem,trinum)%LagNormVec(2) &
                 + (GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(3,0,0,0,iElem)) &
                 * SurfLagValues(iLocSide, iElem,trinum)%LagNormVec(3)
        IF (NVecTest.LE.0.0) THEN
          WRITE(*,*)'=============================================='
          WRITE(*,*)' ERROR in Calculation of NormVec for Element:', iElem
          STOP
        END IF    

      END DO
      CALL SetMeanSurfValues(iLocSide, iElem)
      NVecTest = (GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(1,0,0,0,iElem)) &  
               * MeanSurfValues(iLocSide, iElem)%MeanNormVec(1) &
               + (GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(2,0,0,0,iElem)) &
               * MeanSurfValues(iLocSide, iElem)%MeanNormVec(2) &
               + (GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,iElem))-Elem_xGP(3,0,0,0,iElem)) & 
               * MeanSurfValues(iLocSide, iElem)%MeanNormVec(3)
      IF (NVecTest.LE.0.0) THEN
        WRITE(*,*)'=============================================='
        WRITE(*,*)' ERROR in Calculation of NormVec for Element:', iElem
        STOP
      END IF
    END DO
  END DO
  ALLOCATE(PartStateBulkValues(PDM%maxParticleNumber,5))
  ALLOCATE(LD_RHS(PDM%maxParticleNumber,3))
  LD_RHS = 0.0

! Set Particle Bulk Values
  DO iSpec = 1, nSpecies
    IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
      IF (.NOT.((CollisMode.EQ.2).OR.(CollisMode.EQ.3))) THEN
        WRITE(UNIT=hilf,FMT='(I2)') iSpec
        SpecDSMC(iSpec)%CharaTVib  = GETREAL('Part-Species'//TRIM(hilf)//'-CharaTempVib','0.')  
        SpecDSMC(iSpec)%Ediss_eV   = GETREAL('Part-Species'//TRIM(hilf)//'-Ediss_eV','0.')
      END IF
    END IF
  END DO

  DO iPart = 1, PDM%maxParticleNumber
    IF (PDM%ParticleInside(ipart)) THEN
      IF (Species(PartSpecies(iPart))%NumberOfInits.EQ.0) THEN
        PartStateBulkValues(iPart,1) = Species(PartSpecies(iPart))%VeloVecIC(1) * Species(PartSpecies(iPart))%VeloIC
        PartStateBulkValues(iPart,2) = Species(PartSpecies(iPart))%VeloVecIC(2) * Species(PartSpecies(iPart))%VeloIC
        PartStateBulkValues(iPart,3) = Species(PartSpecies(iPart))%VeloVecIC(3) * Species(PartSpecies(iPart))%VeloIC
        PartStateBulkValues(iPart,4) = Species(PartSpecies(iPart))%MWTemperatureIC
        PartStateBulkValues(iPart,5) = CalcDegreeOfFreedom(iPart)
      ELSE
        iInit = PDM%PartInit(iPart)
        PartStateBulkValues(iPart,1) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(1) &
                                     * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
        PartStateBulkValues(iPart,2) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(2) &
                                     * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
        PartStateBulkValues(iPart,3) = Species(PartSpecies(iPart))%Init(iInit)%VeloVecIC(3) &
                                     * Species(PartSpecies(iPart))%Init(iInit)%VeloIC
        PartStateBulkValues(iPart,4) = Species(PartSpecies(iPart))%Init(iInit)%MWTemperatureIC
        PartStateBulkValues(iPart,5) = CalcDegreeOfFreedom(iPart)
      END IF
    END IF
  END DO

  SWRITE(UNIT_stdOut,'(A)')' INIT LD DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitLD

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcDegreeOfFreedom(iPart)

  USE MOD_LD_Vars
  USE MOD_DSMC_Vars,          ONLY : SpecDSMC
  USE MOD_Particle_Vars,      ONLY : Species, PartSpecies, BoltzmannConst
!--------------------------------------------------------------------------------------------------!
! perform chemical init
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE 
! LOCAL VARIABLES
!--------------------------------------------------------------------------------------------------!
  REAL                          :: ZetaRot, ZetaVib, ModTvibToTemp, TvibToTemp, JToEv
  INTEGER                       :: iSpec
!--------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!--------------------------------------------------------------------------------------------------!
  INTEGER, INTENT(IN)           :: iPart
!#ifdef MPI
!#endif
!===================================================================================================
  JToEv = 1.602176565E-19
  iSpec = PartSpecies(iPart)
  IF(SpecDSMC(iSpec)%InterID.EQ.2) THEN
    ZetaRot = 2.0
    TvibToTemp = SpecDSMC(iSpec)%CharaTVib / PartStateBulkValues(iPart,4)
    ModTvibToTemp = SpecDSMC(iSpec)%Ediss_eV * JToEv / (BoltzmannConst * PartStateBulkValues(iPart,4))
    ZetaVib = 2.0 * TvibToTemp / (EXP(TvibToTemp)-1) - 2.0 * ModTvibToTemp / (EXP(ModTvibToTemp)-1)
  ELSE
    ZetaRot = 0.0
    ZetaVib = 0.0
  END IF
  CalcDegreeOfFreedom = 3.0 + ZetaRot + ZetaVib

END FUNCTION CalcDegreeOfFreedom

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

REAL FUNCTION CalcTriNumArea(iLocSide, Element, trinum)
  USE MOD_Particle_Vars,          ONLY : GEO
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER, INTENT(IN)         :: iLocSide, Element, trinum
  INTEGER                     :: Nod2, Nod3
  REAL                        :: xNod1, xNod2, xNod3 
  REAL                        :: yNod1, yNod2, yNod3 
  REAL                        :: zNod1, zNod2, zNod3 
  REAL                        :: Vector1(1:3), Vector2(1:3)

!--------------------------------------------------------------------------------------------------!

!--- Node 1 ---
   xNod1 = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNod1 = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNod1 = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

!--- Node 2 ---
   Nod2 = trinum + 1      ! vector 1-2 for first triangle
                          ! vector 1-3 for second triangle
   xNod2 = GEO%NodeCoords(1,GEO%ElemSideNodeID(Nod2,iLocSide,Element))
   yNod2 = GEO%NodeCoords(2,GEO%ElemSideNodeID(Nod2,iLocSide,Element))
   zNod2 = GEO%NodeCoords(3,GEO%ElemSideNodeID(Nod2,iLocSide,Element))

!--- Node 3 ---
   Nod3 = trinum + 2      ! vector 1-3 for first triangle
                          ! vector 1-4 for second triangle
   xNod3 = GEO%NodeCoords(1,GEO%ElemSideNodeID(Nod3,iLocSide,Element))
   yNod3 = GEO%NodeCoords(2,GEO%ElemSideNodeID(Nod3,iLocSide,Element))
   zNod3 = GEO%NodeCoords(3,GEO%ElemSideNodeID(Nod3,iLocSide,Element))

   Vector1(1) = xNod2 - xNod1
   Vector1(2) = yNod2 - yNod1
   Vector1(3) = zNod2 - zNod1

   Vector2(1) = xNod3 - xNod1
   Vector2(2) = yNod3 - yNod1
   Vector2(3) = zNod3 - zNod1

!--- 2 * Area = |cross product| of vector 1-2 and 1-3 for first triangle or
                                 ! vector 1-3 and 1-4 for second triangle
   CalcTriNumArea = 0.5*SQRT((Vector1(2)*Vector2(3)-Vector1(3)*Vector2(2))**2 &
           + (-Vector1(1)*Vector2(3)+Vector1(3)*Vector2(1))**2 &
           + (Vector1(1)*Vector2(2)-Vector1(2)*Vector2(1))**2)

  RETURN

END FUNCTION CalcTriNumArea
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE CalcLagNormVec(iLocSide, Element, trinum)

  USE MOD_LD_Vars
  USE MOD_Particle_Vars,          ONLY : GEO
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER, INTENT(IN)         :: iLocSide, Element, trinum
  INTEGER                     :: Nod2, Nod3
  REAL                        :: xNod1, xNod2, xNod3 
  REAL                        :: yNod1, yNod2, yNod3 
  REAL                        :: zNod1, zNod2, zNod3 
  REAL                        :: Vector1(1:3), Vector2(1:3)
  REAL                        :: nx, ny, nz, nVal
!--------------------------------------------------------------------------------------------------!
!--- Node 1 ---
   xNod1 = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
   yNod1 = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
   zNod1 = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))

!--- Node 2 ---
   Nod2 = trinum + 1      ! vector 1-2 for first triangle
                          ! vector 1-3 for second triangle
   xNod2 = GEO%NodeCoords(1,GEO%ElemSideNodeID(Nod2,iLocSide,Element))
   yNod2 = GEO%NodeCoords(2,GEO%ElemSideNodeID(Nod2,iLocSide,Element))
   zNod2 = GEO%NodeCoords(3,GEO%ElemSideNodeID(Nod2,iLocSide,Element))

!--- Node 3 ---
   Nod3 = trinum + 2      ! vector 1-3 for first triangle
                          ! vector 1-4 for second triangle
   xNod3 = GEO%NodeCoords(1,GEO%ElemSideNodeID(Nod3,iLocSide,Element))
   yNod3 = GEO%NodeCoords(2,GEO%ElemSideNodeID(Nod3,iLocSide,Element))
   zNod3 = GEO%NodeCoords(3,GEO%ElemSideNodeID(Nod3,iLocSide,Element))

   Vector1(1) = xNod2 - xNod1
   Vector1(2) = yNod2 - yNod1
   Vector1(3) = zNod2 - zNod1

   Vector2(1) = xNod3 - xNod1
   Vector2(2) = yNod3 - yNod1
   Vector2(3) = zNod3 - zNod1

   nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2) ! n is inward normal vector
   ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
   nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)

   nVal = SQRT(nx*nx + ny*ny + nz*nz)

   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(1) = nx/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(2) = ny/nVal
   SurfLagValues(iLocSide, Element,trinum)%LagNormVec(3) = nz/nVal

END SUBROUTINE CalcLagNormVec
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
SUBROUTINE SetMeanSurfValues(iLocSide, Element)

  USE MOD_LD_Vars
  USE MOD_Particle_Vars,          ONLY : GEO
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
  INTEGER, INTENT(IN)         :: iLocSide, Element                                                 !
! Local variable declaration  
  REAL                        :: xNod1, xNod2, xNod3, xNod4
  REAL                        :: yNod1, yNod2, yNod3, yNod4 
  REAL                        :: zNod1, zNod2, zNod3, zNod4 
  REAL                        :: Vector1(3), Vector2(3),BaseVectorS(3) 
  REAL                        :: nx, ny, nz, nVal
!--------------------------------------------------------------------------------------------------!

  !--- Node 1 ---
  xNod1 = GEO%NodeCoords(1,GEO%ElemSideNodeID(1,iLocSide,Element))
  yNod1 = GEO%NodeCoords(2,GEO%ElemSideNodeID(1,iLocSide,Element))
  zNod1 = GEO%NodeCoords(3,GEO%ElemSideNodeID(1,iLocSide,Element))
  !--- Node 2 ---
  xNod2 = GEO%NodeCoords(1,GEO%ElemSideNodeID(2,iLocSide,Element))
  yNod2 = GEO%NodeCoords(2,GEO%ElemSideNodeID(2,iLocSide,Element))
  zNod2 = GEO%NodeCoords(3,GEO%ElemSideNodeID(2,iLocSide,Element))
  !--- Node 3 ---
  xNod3 = GEO%NodeCoords(1,GEO%ElemSideNodeID(3,iLocSide,Element))
  yNod3 = GEO%NodeCoords(2,GEO%ElemSideNodeID(3,iLocSide,Element))
  zNod3 = GEO%NodeCoords(3,GEO%ElemSideNodeID(3,iLocSide,Element))
  !--- Node 4 ---
  xNod4 = GEO%NodeCoords(1,GEO%ElemSideNodeID(4,iLocSide,Element))
  yNod4 = GEO%NodeCoords(2,GEO%ElemSideNodeID(4,iLocSide,Element))
  zNod4 = GEO%NodeCoords(3,GEO%ElemSideNodeID(4,iLocSide,Element))
  Vector1(1) = xNod3 - xNod1
  Vector1(2) = yNod3 - yNod1
  Vector1(3) = zNod3 - zNod1
  Vector2(1) = xNod4 - xNod2
  Vector2(2) = yNod4 - yNod2
  Vector2(3) = zNod4 - zNod2
  nx = Vector1(2) * Vector2(3) - Vector1(3) * Vector2(2) ! n is inward normal vector
  ny = Vector1(3) * Vector2(1) - Vector1(1) * Vector2(3)
  nz = Vector1(1) * Vector2(2) - Vector1(2) * Vector2(1)
  BaseVectorS(1:3) = 0.25 *( &
                   + GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1,iLocSide,Element)) &
                   + GEO%NodeCoords(1:3,GEO%ElemSideNodeID(2,iLocSide,Element)) &
                   + GEO%NodeCoords(1:3,GEO%ElemSideNodeID(3,iLocSide,Element)) &
                   + GEO%NodeCoords(1:3,GEO%ElemSideNodeID(4,iLocSide,Element)) )
  nVal = SQRT(nx*nx + ny*ny + nz*nz)
  MeanSurfValues(iLocSide, Element)%MeanNormVec(1) = nx/nVal
  MeanSurfValues(iLocSide, Element)%MeanNormVec(2) = ny/nVal
  MeanSurfValues(iLocSide, Element)%MeanNormVec(3) = nz/nVal
  MeanSurfValues(iLocSide, Element)%MeanBaseD = MeanSurfValues(iLocSide, Element)%MeanNormVec(1) * BaseVectorS(1) &
                                  + MeanSurfValues(iLocSide, Element)%MeanNormVec(2) * BaseVectorS(2) &
                                  + MeanSurfValues(iLocSide, Element)%MeanNormVec(3) * BaseVectorS(3)

END SUBROUTINE SetMeanSurfValues
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

END MODULE MOD_LD_Init
