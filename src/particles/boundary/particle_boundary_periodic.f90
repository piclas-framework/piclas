MODULE MOD_part_boundary_periodic
#include "boltzplatz.h"
!===================================================================================================================================
!===================================================================================================================================
  IMPLICIT NONE
  PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
  PUBLIC                 :: InitPeriodic
!===================================================================================================================================


CONTAINS

SUBROUTINE InitPeriodic()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_ReadInTools
  USE MOD_Particle_Vars, ONLY : GEO, PartBound, PDM 
  USE MOD_part_MPI_Vars, ONLY : NbrOfCases, casematrix, partShiftVector
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
  INTEGER                :: iVec, ind, ind2
  CHARACTER(32)          :: hilf
!===================================================================================================================================

  GEO%nPeriodicVectors       = GETINT('Part-nPeriodicVectors','0')
  DO iVec = 1, SIZE(PartBound%TargetBoundCond)
    IF((PartBound%TargetBoundCond(iVec).EQ.PartBound%PeriodicBC).AND.(GEO%nPeriodicVectors.EQ.0))THEN
      WRITE(*,*) 'Part-PeriodicVectors need to be assigned in the ini file'
      STOP
    END IF
  END DO
  ALLOCATE(GEO%PeriodicVectors(1:3,1:GEO%nPeriodicVectors))
  DO iVec = 1, GEO%nPeriodicVectors
    WRITE(UNIT=hilf,FMT='(I2)') iVec
    GEO%PeriodicVectors(1:3,iVec) = GETREALARRAY('Part-PeriodicVector'//TRIM(hilf),3,'1.,0.,0.')
  END DO
  IF(GEO%nPeriodicVectors.GE.0)THEN
    CALL getPeriodicVectors()
    CALL MapPeriodicVectorsToSides()
  END IF
  ! build periodic case matrix
  IF (GEO%nPeriodicVectors.GT.0) THEN
    ! build case matrix
    NbrOfCases = 3**GEO%nPeriodicVectors
    ALLOCATE(casematrix(1:NbrOfCases,1:3))
    casematrix(:,:) = 0
    IF (GEO%nPeriodicVectors.EQ.1) THEN
      casematrix(1,1) = 1
      casematrix(3,1) = -1
    END IF
    IF (GEO%nPeriodicVectors.EQ.2) THEN
      casematrix(1:3,1) = 1
      casematrix(7:9,1) = -1
      DO ind = 1,3
        casematrix(ind*3-2,2) = 1
        casematrix(ind*3,2) = -1
      END DO
    END IF
    IF (GEO%nPeriodicVectors.EQ.3) THEN
      casematrix(1:9,1) = 1
      casematrix(19:27,1) = -1
      DO ind = 1,3
        casematrix(ind*9-8:ind*9-6,2) = 1
        casematrix(ind*9-2:ind*9,2) = -1
        DO ind2 = 1,3
          casematrix((ind2*3-2)+(ind-1)*9,3) = 1
          casematrix((ind2*3)+(ind-1)*9,3) = -1
        END DO
      END DO
    END IF
  ELSE
    NbrOfCases = 1
    ALLOCATE(casematrix(1:1,1:3))
    casematrix(:,:) = 0
  END IF
  IF (GEO%nPeriodicVectors.GT.0) THEN
    ALLOCATE(partShiftVector(1:3,1:PDM%maxParticleNumber))
    partShiftVector = 0.
  END IF
END SUBROUTINE InitPeriodic

SUBROUTINE MapPeriodicVectorsToSides()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_ReadInTools
  USE MOD_Particle_Vars,      ONLY : GEO 
  USE MOD_Mesh_Vars,          ONLY : nElems 
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  INTEGER                :: iElem, iLocSide, iPV, maxindex
  REAL, ALLOCATABLE      :: PVNorm(:,:)
  REAL                   :: PVabs, v1(1:3), v2(1:3), normal(1:3), normalAbs, dotprod, maxvalue
  LOGICAL                :: neg
!===================================================================================================================================

  !--- normalize periodic vectors
  ALLOCATE(PVNorm(1:3,1:GEO%nPeriodicVectors))
  DO iPV = 1, GEO%nPeriodicVectors
    PVabs=SQRT(GEO%PeriodicVectors(1,iPV)**2+GEO%PeriodicVectors(2,iPV)**2+GEO%PeriodicVectors(3,iPV)**2)
    PVNorm(1:3,iPV) = GEO%PeriodicVectors(1:3,iPV)/PVabs
  END DO
  !--- go through all elements and all local sides
  DO iElem = 1, nElems
    DO iLocSide = 1,6
      IF (GEO%PeriodicElemSide(iLocSide,iElem).EQ.-1) THEN   ! it is a periodic side    
        ! --- build normal of the side (assumption: nearly planar side, should be ok for periodic)
        v1(1:3)   = GEO%NodeCoords(1:3,GEO%ElemSideNodeID(2,iLocSide,iElem)) - &
                    GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1,iLocSide,iElem))
        v2(1:3)   = GEO%NodeCoords(1:3,GEO%ElemSideNodeID(3,iLocSide,iElem)) - &
                    GEO%NodeCoords(1:3,GEO%ElemSideNodeID(1,iLocSide,iElem))
        normal(1) = v1(2)*v2(3)-v1(3)*v2(2)
        normal(2) = v1(3)*v2(1)-v1(1)*v2(3)
        normal(3) = v1(1)*v2(2)-v1(2)*v2(1)
        normalAbs = SQRT(normal(1)**2+normal(2)**2+normal(3)**2)
        normal(1:3) = normal(1:3)/normalAbs
        ! --- check which periodic vector it fits best
        maxvalue = 0.
        maxindex = 0
        neg = .FALSE.
        DO iPV = 1, GEO%nPeriodicVectors
          dotprod = PVNorm(1,iPV)*normal(1) + PVNorm(2,iPV)*normal(2) + PVNorm(3,iPV)*normal(3)
          IF(abs(dotprod).GE.maxvalue)THEN
            maxvalue = abs(dotprod)
            maxindex = iPV
            IF (dotprod.LT.0) ThEN
              neg = .TRUE.
            ELSE
              neg = .FALSE.
            END IF
          END IF
        END DO
        ! --- mapping
        IF (neg) THEN
          GEO%PeriodicElemSide(iLocSide,iElem) = -maxindex
        ELSE
          GEO%PeriodicElemSide(iLocSide,iElem) = maxindex
        END IF
      END IF
    END DO
  END DO
END SUBROUTINE MapPeriodicVectorsToSides

SUBROUTINE getPeriodicVectors()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Globals,       ONLY : Logging,UNIT_errOut,UNIT_logOut,abort
  USE MOD_Particle_Vars, ONLY : GEO
  USE MOD_PICDepo_Vars
!  USE TypesDef_mod    , ONLY : MESH,CALC,UNIT_LogOut,UNIT_ErrOut
!  USE TypesDef_PIC_mod, ONLY : PIC
!  USE pointerbasis_mod, ONLY : abort
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! ARGUMENT LIST DECLARATION
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE DECLARATION
  LOGICAL                :: directions(1:3)
  INTEGER                :: iPV
  REAL                   :: eps(1:3)
!===================================================================================================================================

  LOGWRITE(*,*)'nPeriodicVectors = ',GEO%nPeriodicVectors
  IF ((GEO%nPeriodicVectors.GT.3).OR.(GEO%nPeriodicVectors.LT.0)) THEN
    CALL abort(__STAMP__,&
                        'nPeriodicVectors must be >= 0 and <= 3!',GEO%nPeriodicVectors,999.)
  ELSE IF (GEO%nPeriodicVectors.EQ.0) THEN
    RETURN
  END IF
  directions(1:3)=.FALSE.
  DO iPV = 1,GEO%nPeriodicVectors
    LOGWRITE(*,*)'PeriodicVectors(1:3),',iPV,')=',GEO%PeriodicVectors(1:3,iPV)
    IF      (GEO%PeriodicVectors(1,iPV).NE.0) THEN
      IF ((GEO%PeriodicVectors(2,iPV).NE.0).OR.(GEO%PeriodicVectors(3,iPV).NE.0)) THEN
        CALL abort(__STAMP__,&
                            'Periodic Vector not in Cartesian direction!',iPV,999.)
      END IF
      IF (.NOT.directions(1)) THEN
        directions(1) = .TRUE.
      ELSE
        CALL abort(__STAMP__,&
                            '2 Periodic Vectors in x-direction!',iPV,999.)
      END IF
    ELSE IF (GEO%PeriodicVectors(2,iPV).NE.0) THEN
      IF ((GEO%PeriodicVectors(1,iPV).NE.0).OR.(GEO%PeriodicVectors(3,iPV).NE.0)) THEN
        CALL abort(__STAMP__,&
                            'Periodic Vector not in Cartesian direction!',iPV,999.)
      END IF
      IF (.NOT.directions(2)) THEN
        directions(2) = .TRUE.
      ELSE
        CALL abort(__STAMP__,&
                            '2 Periodic Vectors in y-direction!',iPV,999.)
      END IF
    ELSE IF (GEO%PeriodicVectors(3,iPV).NE.0) THEN
      IF ((GEO%PeriodicVectors(1,iPV).NE.0).OR.(GEO%PeriodicVectors(2,iPV).NE.0)) THEN
        CALL abort(__STAMP__,&
                            'Periodic Vector not in Cartesian direction!',iPV,999.)
      END IF
      IF (.NOT.directions(3)) THEN
        directions(3) = .TRUE.
      ELSE
        CALL abort(__STAMP__,&
                            '2 Periodic Vectors in z-direction!',iPV,999.)
      END IF
    ELSE
      CALL abort(__STAMP__,'Periodic Vector = 0!',iPV,999.)
    END IF
  END DO
  eps(1)=1.E-9*(GEO%FIBGMDeltas(1))
  eps(2)=1.E-9*(GEO%FIBGMDeltas(2))
  eps(3)=1.E-9*(GEO%FIBGMDeltas(3))
  IF      (ABS(SUM(GEO%PeriodicVectors(1,:))-NINT(SUM(GEO%PeriodicVectors(1,:))/GEO%FIBGMDeltas(1))*GEO%FIBGMDeltas(1)) &
           .GT.eps(1)) THEN
    ERRWRITE(*,*)'SUM(PeriodicVectors(1,:))   =',SUM(GEO%PeriodicVectors(1,:))
    ERRWRITE(*,*)'GEO%FIBGMDeltas(1)          =',GEO%FIBGMDeltas(1)
    ERRWRITE(*,*)'1.E-9*(FIBGMDeltas(1))      =',eps(1)
    ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(1))*D(1))=',ABS(SUM(GEO%PeriodicVectors(1,:))-&
                                                NINT(SUM(GEO%PeriodicVectors(1,:))/GEO%FIBGMDeltas(1))*GEO%FIBGMDeltas(1))
    CALL abort(__STAMP__,&
                        'Periodic Vector in x-direction is not a multiple of FIBGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(1,:))-NINT(SUM(GEO%PeriodicVectors(1,:))/GEO%FIBGMDeltas(1))*GEO%FIBGMDeltas(1)))
  ELSE IF (ABS(SUM(GEO%PeriodicVectors(2,:))-NINT(SUM(GEO%PeriodicVectors(2,:))/GEO%FIBGMDeltas(2))*GEO%FIBGMDeltas(2)) &
           .GT.eps(2)) THEN
    ERRWRITE(*,*)'SUM(PeriodicVectors(2,:))   =',SUM(GEO%PeriodicVectors(2,:))
    ERRWRITE(*,*)'GEO%FIBGMDeltas(2)          =',GEO%FIBGMDeltas(2)
    ERRWRITE(*,*)'1.E-9*(FIBGMDeltas(2))      =',eps(2)
    ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(2))*D(2))=',ABS(SUM(GEO%PeriodicVectors(2,:))-&
                                                NINT(SUM(GEO%PeriodicVectors(2,:))/GEO%FIBGMDeltas(2))*GEO%FIBGMDeltas(2))
    CALL abort(__STAMP__,&
                        'Periodic Vector in y-direction is not a multiple of FIBGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(2,:))-NINT(SUM(GEO%PeriodicVectors(2,:))/GEO%FIBGMDeltas(2))*GEO%FIBGMDeltas(2)))
  ELSE IF (ABS(SUM(GEO%PeriodicVectors(3,:))-NINT(SUM(GEO%PeriodicVectors(3,:))/GEO%FIBGMDeltas(3))*GEO%FIBGMDeltas(3)) &
           .GT.eps(3)) THEN
    ERRWRITE(*,*)'SUM(PeriodicVectors(3,:))   =',SUM(GEO%PeriodicVectors(3,:))
    ERRWRITE(*,*)'GEO%FIBGMDeltas(3)          =',GEO%FIBGMDeltas(3)
    ERRWRITE(*,*)'1.E-9*(FIBGMDeltas(3))      =',eps(3)
    ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(3))*D(3))=',ABS(SUM(GEO%PeriodicVectors(3,:))-&
                                                NINT(SUM(GEO%PeriodicVectors(3,:))/GEO%FIBGMDeltas(3))*GEO%FIBGMDeltas(3))
    CALL abort(__STAMP__,&
                        'Periodic Vector in z-direction is not a multiple of FIBGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(3,:))-NINT(SUM(GEO%PeriodicVectors(3,:))/GEO%FIBGMDeltas(3))*GEO%FIBGMDeltas(3)))
  END IF
  IF((DepositionType.EQ.'cartmesh_volumeweighting').OR.(DepositionType.EQ.'cartmesh_splines'))THEN
    IF (ABS(SUM(GEO%PeriodicVectors(1,:))-NINT(SUM(GEO%PeriodicVectors(1,:))/BGMDeltas(1))*BGMDeltas(1)) &
         .GT.eps(1)) THEN
      ERRWRITE(*,*)'SUM(PeriodicVectors(1,:))   =',SUM(GEO%PeriodicVectors(1,:))
      ERRWRITE(*,*)'BGMDeltas(1)                =',BGMDeltas(1)
      ERRWRITE(*,*)'1.E-9*(BGMDeltas(1))        =',eps(1)
      ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(1))*D(1))=',ABS(SUM(GEO%PeriodicVectors(1,:))-&
           NINT(SUM(GEO%PeriodicVectors(1,:))/BGMDeltas(1))*BGMDeltas(1))
      CALL abort(__STAMP__,&
           'Periodic Vector in x-direction is not a multiple of BGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(1,:))-NINT(SUM(GEO%PeriodicVectors(1,:))/BGMDeltas(1))*BGMDeltas(1)))
    ELSE IF (ABS(SUM(GEO%PeriodicVectors(2,:))-NINT(SUM(GEO%PeriodicVectors(2,:))/BGMDeltas(2))*BGMDeltas(2)) &
         .GT.eps(2)) THEN
      ERRWRITE(*,*)'SUM(PeriodicVectors(2,:))   =',SUM(GEO%PeriodicVectors(2,:))
      ERRWRITE(*,*)'BGMDeltas(2)                =',BGMDeltas(2)
      ERRWRITE(*,*)'1.E-9*(BGMDeltas(2))        =',eps(2)
      ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(2))*D(2))=',ABS(SUM(GEO%PeriodicVectors(2,:))-&
           NINT(SUM(GEO%PeriodicVectors(2,:))/BGMDeltas(2))*BGMDeltas(2))
      CALL abort(__STAMP__,&
           'Periodic Vector in y-direction is not a multiple of BGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(2,:))-NINT(SUM(GEO%PeriodicVectors(2,:))/BGMDeltas(2))*BGMDeltas(2)))
    ELSE IF (ABS(SUM(GEO%PeriodicVectors(3,:))-NINT(SUM(GEO%PeriodicVectors(3,:))/BGMDeltas(3))*BGMDeltas(3)) &
         .GT.eps(3)) THEN
      ERRWRITE(*,*)'SUM(PeriodicVectors(3,:))   =',SUM(GEO%PeriodicVectors(3,:))
      ERRWRITE(*,*)'BGMDeltas(3)                =',BGMDeltas(3)
      ERRWRITE(*,*)'1.E-9*(BGMDeltas(3))      =',eps(3)
      ERRWRITE(*,*)'ABS(SUM-NINT(SUM/D(3))*D(3))=',ABS(SUM(GEO%PeriodicVectors(3,:))-&
           NINT(SUM(GEO%PeriodicVectors(3,:))/BGMDeltas(3))*BGMDeltas(3))
      CALL abort(__STAMP__,&
           'Periodic Vector in z-direction is not a multiple of BGMDeltas!',999,&
           ABS(SUM(GEO%PeriodicVectors(3,:))-NINT(SUM(GEO%PeriodicVectors(3,:))/BGMDeltas(3))*BGMDeltas(3)))
    END IF
  END IF
END SUBROUTINE getPeriodicVectors

END MODULE MOD_part_boundary_periodic
