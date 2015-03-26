#include "boltzplatz.h"

MODULE MOD_Partilce_Periodic_BC
!===================================================================================================================================
! Module initialization of periodic vectors for particle treatment
!===================================================================================================================================
IMPLICIT NONE
PRIVATE

INTERFACE InitPeriodicBC
  MODULE PROCEDURE InitPeriodicBC
END INTERFACE
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC                 :: InitPeriodicBC
!===================================================================================================================================


CONTAINS

SUBROUTINE InitPeriodicBC()
!===================================================================================================================================
! Computes the periodic-displacement vector 
! Both periodic sides have to be planer and parallel!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools,        ONLY:GETINT,GETREALARRAY
USE MOD_Particle_Mesh_Vars, ONLY:GEO
USE MOD_Particle_Vars,      ONLY:PartBound, PDM 
USE MOD_Particle_MPI_Vars,  ONLY:NbrOfCases, casematrix!, partShiftVector
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                :: iVec, ind, ind2
CHARACTER(32)          :: hilf
!===================================================================================================================================


GEO%nPeriodicVectors       = GETINT('Part-nPeriodicVectors','0')
DO iVec = 1, SIZE(PartBound%TargetBoundCond)
  IF((PartBound%TargetBoundCond(iVec).EQ.PartBound%PeriodicBC).AND.(GEO%nPeriodicVectors.EQ.0))THEN
    CALL abort(__STAMP__,&
        'Part-PeriodicVectors need to be assigned in the ini file')
  END IF
END DO

ALLOCATE(GEO%PeriodicVectors(1:3,1:GEO%nPeriodicVectors))
DO iVec = 1, GEO%nPeriodicVectors
  WRITE(UNIT=hilf,FMT='(I2)') iVec
  GEO%PeriodicVectors(1:3,iVec) = GETREALARRAY('Part-PeriodicVector'//TRIM(hilf),3,'1.,0.,0.')
END DO

CALL GetPeriodicVectors()
CALL MapPeriodicVectorsToSides()

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
!IF (GEO%nPeriodicVectors.GT.0) THEN
!  ALLOCATE(partShiftVector(1:3,1:PDM%maxParticleNumber))
!  partShiftVector = 0.
!END IF

END SUBROUTINE InitPeriodicBC

SUBROUTINE MapPeriodicVectorsToSides()
!===================================================================================================================================
! Per direction, two different periodic displacements are possible: +- periodic_vector
! each side gets an periodictype:
! 0 - non-periodic
! 1 - +pv(1)
! 2 - -pv(1)
! ...
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals,                 ONLY:AlmostZero,AlmostEqual,CROSSNORM
!USE MOD_Mesh_Vars,               ONLY:Elems,offsetElem,nSides, ElemToSide
USE MOD_Mesh_Vars,               ONLY:nBCSides,nInnerSides,nMPISides_MINE,nMPISides_YOUR
USE MOD_Mesh_Vars,               ONLY:nSides, ElemToSide,nBCSides,XCL_NGeo,NGeo
USE MOD_Particle_Mesh_Vars,      ONLY:GEO, SidePeriodicType, SidePeriodicDisplacement
USE MOD_Particle_Surfaces_Vars,  ONLY:SideNormVec
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iElem,ilocSide,SideID,flip,iSide
INTEGER             :: iPV,iDisplace,nDisplacement
INTEGER             :: nPeriodicSides,tmpBCSides
INTEGER             :: nLocalSides
REAL                :: xmin,xmax,ymin,ymax,zmin,zmax
REAL                :: v1(3),v2(3),nVec(3)
REAL,ALLOCATABLE    :: normDisplaceVec(:,:)
!===================================================================================================================================

nDisplacement=2*GEO%nPeriodicVectors
!ALLOCATE(SidePeriodicType(1:nSides)                 &
!        ,normDisplaceVec(1:3,nDisplacement)         &
!        ,SidePeriodicDisplacement(1:3,nDisplacement))
!    k

ALLOCATE(normDisplaceVec(1:3,nDisplacement)         &
        ,SidePeriodicDisplacement(1:3,nDisplacement))

!SidePeriodicType=0

iDisplace=1
DO iPV=1,GEO%nPeriodicVectors
  SidePeriodicDisplacement(:,iDisplace)  = GEO%PeriodicVectors(:,iPV)
  normDisplaceVec(:,iDisplace  ) = SidePeriodicDisplacement(:,iDisplace) &
                                 / SQRT(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),SidePeriodicDisplacement(:,iDisplace)))
  SidePeriodicDisplacement(:,iDisplace+1)=-GEO%PeriodicVectors(:,iPV)
  normDisplaceVec(:,iDisplace+1) =-normDisplaceVec(:,iDisplace+1)
  iDisplace=iDisplace+2
END DO ! iPV

nLocalSides = nBCSides+nInnerSides+nMPISides_MINE
!--- Initialize Periodic Side Info
DO iElem=1,PP_nElems
  DO ilocSide=1,6
    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
    flip   = ElemToSide(E2S_FLIP,ilocSide,iElem)
    IF(SideperiodicType(SideID).EQ.-1)THEN
      SELECT CASE(iLocSide)
      CASE(XI_MINUS)
        v1=XCL_NGeo(1:3,0   ,NGeo,0,iElem) - XCL_NGeo(1:3,0   ,0,0,iElem)
        v2=XCL_NGeo(1:3,0   ,0,NGeo,iElem) - XCL_NGeo(1:3,0   ,0,0,iElem)
      CASE(XI_PLUS)
        v1=XCL_NGeo(1:3,NGeo,NGeo,0,iElem) - XCL_NGeo(1:3,NGeo,0,0,iElem)
        v2=XCL_NGeo(1:3,NGeo,0,NGeo,iElem) - XCL_NGeo(1:3,NGeo,0,0,iElem)
      CASE(ETA_MINUS)
        v1=XCL_NGeo(1:3,NGeo,0,0   ,iElem)-XCL_NGeo(1:3,0,0,0,iElem)
        v2=XCL_NGeo(1:3,0   ,0,Ngeo,iElem)-XCL_NGeo(1:3,0,0,0,iElem)
      CASE(ETA_PLUS)
        v1=XCL_NGeo(1:3,NGeo,NGeo,0   ,iElem)-XCL_NGeo(1:3,0,NGeo,0,iElem)
        v2=XCL_NGeo(1:3,0   ,NGeo,NGeo,iElem)-XCL_NGeo(1:3,0,NGeo,0,iElem)
      CASE(ZETA_MINUS)
        v1=XCL_NGeo(1:3,NGeo,0,0   ,iElem)-XCL_NGeo(1:3,0,0,0   ,iElem)
        v2=XCL_NGeo(1:3,0,NGeo,0   ,iElem)-XCL_NGeo(1:3,0,0,0   ,iElem)
      CASE(ZETA_PLUS)
        v1=XCL_NGeo(1:3,NGeo,0,NGeo,iElem)-XCL_NGeo(1:3,0,0,NGeo,iElem)
        v2=XCL_NGeo(1:3,0,NGeo,NGeo,iElem)-XCL_NGeo(1:3,0,0,NGeo,iElem)
      END SELECT
      nVec=CROSSNORM(v1,v2)
      IF(SideID.LE.nLocalSides)THEN
        IF(flip.EQ.0)THEN ! master side
          DO iDisplace=1,nDisplacement
            IF(ALMOSTEQUAL(DOT_PRODUCT(SideNormVec(:,SideID),SidePeriodicDisplacement(:,iDisplace)),-1.0)) &
                SidePeriodicType(SideID)=iDisplace
          END DO ! iDisplace
        END IF ! flip
      ELSE ! mpi Side & no master flip
        DO iDisplace=1,nDisplacement
          IF(ALMOSTEQUAL(DOT_PRODUCT(SideNormVec(:,SideID),SidePeriodicDisplacement(:,iDisplace)),1.0)) &
              SidePeriodicType(SideID)=iDisplace
        END DO ! iDisplace
      END IF ! SideID.GT.nLocalSides
    END IF ! is periodic side
  END DO ! ilocSide
END DO  ! iElem

!DO iElem=1,PP_nElems
!  DO ilocSide=1,6
!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    flip   = ElemToSide(E2S_FLIP,ilocSide,iElem)
!    ! check if periodic side
!!    IF ((Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%BCindex.GT.0)&
!!      .AND.(ASSOCIATED(Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%connection))) THEN
!    IF(SideperiodicType(SideID).EQ.-1)THEN
!      ! method with flip
!      IF(flip.EQ.0)THEN ! master side
!        DO iDisplace=1,nDisplacement
!          IF(ALMOSTEQUAL(DOT_PRODUCT(SideNormVec(:,SideID),SidePeriodicDisplacement(:,iDisplace)),-1.0)) &
!              SidePeriodicType(SideID)=iDisplace
!        END DO ! iDisplace
!      ELSE ! slave side
!        DO iDisplace=1,nDisplacement
!          IF(ALMOSTEQUAL(DOT_PRODUCT(SideNormVec(:,SideID),SidePeriodicDisplacement(:,iDisplace)),1.0)) &
!              SidePeriodicType(SideID)=iDisplace
!        END DO ! iDisplace
!      END IF ! flip
!    END IF
!  END DO
!END DO


!!--- Initialize Periodic Side Info
!DO iElem=1,PP_nElems
!  DO ilocSide=1,6
!    SideID = ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
!    flip   = ElemToSide(E2S_FLIP,ilocSide,iElem)
!    ! check if periodic side
!    IF ((Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%BCindex.GT.0)&
!      .AND.(ASSOCIATED(Elems(iElem+offsetElem)%ep%Side(iLocSide)%sp%connection))) THEN
! method with geometry data
!      ! find which periodic side, remember, this is a cartesian grid
!      IF(ALMOSTEQUAL(ABS(DOT_PRODUCT(SideNormVec(:,SideID),(/1.0,0.0,0.0/))),1.0))THEN
!        IF(ALMOSTEQUAL(BezierControlPoints3D(1,0,0,SideID),GEO%xminglob))THEN
!          ! requires positive mapping
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/1.0,0.0,0.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        ELSE ! has to be xmax side
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/-1.0,0.0,0.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        END IF ! AlmostEqual in x
!      END IF
!      !  y-direction
!      IF(ALMOSTEQUAL(ABS(DOT_PRODUCT(SideNormVec(:,SideID),(/0.0,1.0,0.0/))),1.0))THEN
!        IF(ALMOSTEQUAL(BezierControlPoints3D(2,0,0,SideID),GEO%yminglob))THEN
!          ! requires positive mapping
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/0.0,1.0,0.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        ELSE ! has to be ymax side
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/0.0,-1.0,0.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        END IF ! AlmostEqual in y
!      END IF
!      !  z-direction
!      IF(ALMOSTEQUAL(ABS(DOT_PRODUCT(SideNormVec(:,SideID),(/0.0,0.0,1.0/))),1.0))THEN
!        IF(ALMOSTEQUAL(BezierControlPoints3D(3,0,0,SideID),GEO%zminglob))THEN
!          ! requires positive mapping
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/0.0,0.0,1.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        ELSE ! has to be ymax side
!          DO iDisplace=1,nDisplacement
!            IF(DOT_PRODUCT(SidePeriodicDisplacement(:,iDisplace),(/0.0,0.0,-1.0/)).GT.0.) SidePeriodicType(SideID)=iDisplace
!          END DO ! iDisplace
!        END IF ! AlmostEqual in x
!      END IF
!      IF(SidePeriodicType(SideID).EQ.0)THEN
!         CALL abort(__STAMP__,&
!           ' Error with periodic sides for particles. Found no displacement vector for side:', SideID)
!      END IF
!    END IF
!  END DO
!END DO

SDEALLOCATE(normDisplaceVec)

END SUBROUTINE MapPeriodicVectorsToSides

SUBROUTINE GetPeriodicVectors()
!===================================================================================================================================
! Check the periodic vectors for consistency
! For particles, each periodic vector has to stastisfy following conditions
! 1) only a cartesian displacement/ periodicity is supported, e.g. periodicity in x,y,z
! 2) Mesh has to fit into the FIBGM, therefore, the discplacement is a multiple of the FIBGM-delta
! 3) Additionally for PIC with Volume or BSpline weighting/deposition
!    Periodic displacement has to be multiple of BGMdeltas of deposition method
!===================================================================================================================================
! MODULES
USE MOD_Globals,            ONLY: Logging,UNIT_errOut,UNIT_logOut,abort
USE MOD_Particle_Mesh_Vars, ONLY: GEO
USE MOD_PICDepo_Vars
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
LOGICAL                :: directions(1:3)
INTEGER                :: iPV
REAL                   :: eps(1:3)
!===================================================================================================================================

LOGWRITE(*,*)'nPeriodicVectors = ',GEO%nPeriodicVectors
IF ((GEO%nPeriodicVectors.GT.3).OR.(GEO%nPeriodicVectors.LT.0)) THEN
  CALL abort(__STAMP__,&
                      'nPeriodicVectors must be >= 0 and <= 3!',GEO%nPeriodicVectors,999.)
END IF

IF(GEO%nPeriodicVectors.EQ.0) RETURN

! check if all periodic vectors are cartesian
directions(1:3)=.FALSE.
DO iPV = 1,GEO%nPeriodicVectors
  LOGWRITE(*,*)'PeriodicVectors(1:3),',iPV,')=',GEO%PeriodicVectors(1:3,iPV)
  IF(GEO%PeriodicVectors(1,iPV).NE.0) THEN
    IF((GEO%PeriodicVectors(2,iPV).NE.0).OR.(GEO%PeriodicVectors(3,iPV).NE.0)) THEN
      CALL abort(__STAMP__,&
                          'Periodic Vector not in Cartesian direction!',iPV)
    END IF
    IF (.NOT.directions(1)) THEN
      directions(1) = .TRUE.
    ELSE
      CALL abort(__STAMP__,&
                          '2 Periodic Vectors in x-direction!',iPV)
    END IF
  ELSE IF (GEO%PeriodicVectors(2,iPV).NE.0) THEN
    IF ((GEO%PeriodicVectors(1,iPV).NE.0).OR.(GEO%PeriodicVectors(3,iPV).NE.0)) THEN
      CALL abort(__STAMP__,&
                          'Periodic Vector not in Cartesian direction!',iPV)
    END IF
    IF (.NOT.directions(2)) THEN
      directions(2) = .TRUE.
    ELSE
      CALL abort(__STAMP__,&
                          '2 Periodic Vectors in y-direction!',iPV)
    END IF
  ELSE IF (GEO%PeriodicVectors(3,iPV).NE.0) THEN
    IF ((GEO%PeriodicVectors(1,iPV).NE.0).OR.(GEO%PeriodicVectors(2,iPV).NE.0)) THEN
      CALL abort(__STAMP__,&
                          'Periodic Vector not in Cartesian direction!',iPV)
    END IF
    IF (.NOT.directions(3)) THEN
      directions(3) = .TRUE.
    ELSE
      CALL abort(__STAMP__,&
                          '2 Periodic Vectors in z-direction!',iPV)
    END IF
  ELSE
    CALL abort(__STAMP__,'Periodic Vector = 0!',iPV)
  END IF
END DO

! check if periodic vector is multiple of FIBGM-deltas
! some tolerance
eps(1)=1.E-9*(GEO%FIBGMDeltas(1))
eps(2)=1.E-9*(GEO%FIBGMDeltas(2))
eps(3)=1.E-9*(GEO%FIBGMDeltas(3))
IF(ABS(SUM(GEO%PeriodicVectors(1,:))-NINT(SUM(GEO%PeriodicVectors(1,:))/GEO%FIBGMDeltas(1))*GEO%FIBGMDeltas(1)) &
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

! check if periodic vector is multiple of BGM-Delta. This BGM is for the deposition with volume or spline weighting 
! functions
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
END SUBROUTINE GetPeriodicVectors

END MODULE MOD_Partilce_Periodic_BC
