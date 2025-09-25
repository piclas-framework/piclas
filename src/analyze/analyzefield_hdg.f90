!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_AnalyzeField_HDG
!===================================================================================================================================
! Module containing subroutines and functions that are required only for Hybrid Discontinuous Galerkin (HDG) methods
!===================================================================================================================================
! MODULES
#if USE_HDG
USE MOD_Globals, ONLY:UNIT_stdOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: CalculateAverageElectricPotential
PUBLIC :: GetAverageElectricPotentialPlane
PUBLIC :: FinalizeAverageElectricPotential
PUBLIC :: CalculateElectricDisplacementCurrentSurface
#if defined(PARTICLES)
PUBLIC :: CalculateElectricPotentialAndFieldBoundaryVDL
#endif /*defined(PARTICLES)*/
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
SUBROUTINE CalculateAverageElectricPotential()
!===================================================================================================================================
! Calculation of the average electric potential with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is
! missing ... ups
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems, N_SurfMesh, offSetElem
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_Analyze_Vars       ,ONLY: isAverageElecPotSide,AverageElectricPotential,AverageElectricPotentialFaces
USE MOD_Interpolation_Vars ,ONLY: N_Inter,NMax
USE MOD_DG_Vars            ,ONLY: U_N,U_Surf_N,N_DG_Mapping
#if USE_MPI
USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem, SideID,ilocSide,Nloc
INTEGER          :: p,q,l
REAL             :: Uface(PP_nVar,0:NMax,0:NMax)
LOGICAL,PARAMETER:: Prolong=.TRUE.
REAL             :: AverageElectricPotentialProc
REAL             :: area_loc,integral_loc
!===================================================================================================================================

AverageElectricPotentialProc = 0.

DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  Do ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isAverageElecPotSide(SideID)) CYCLE
      IF(Prolong)THEN
#if (PP_NodeType==1) /* for Gauss-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,q,p)=U_N(iElem)%U(:,0,p,q)*N_Inter(Nloc)%L_Minus(0)
              DO l=1,Nloc
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U_N(iElem)%U(:,l,p,q)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,p,q)=U_N(iElem)%U(:,p,0,q)*N_Inter(Nloc)%L_Minus(0)
              DO l=1,Nloc
                Uface(:,p,q)=Uface(:,p,q)+U_N(iElem)%U(:,p,l,q)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_MINUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,q,p)=U_N(iElem)%U(:,p,q,0)*N_Inter(Nloc)%L_Minus(0)
              DO l=1,Nloc
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U_N(iElem)%U(:,p,q,l)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! qfirst stuff
        CASE(XI_PLUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,p,q)=U_N(iElem)%U(:,0,p,q)*N_Inter(Nloc)%L_Plus(0)
              DO l=1,Nloc
                Uface(:,p,q)=Uface(:,p,q)+U_N(iElem)%U(:,l,p,q)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_PLUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,Nloc-p,q)=U_N(iElem)%U(:,p,0,q)*N_Inter(Nloc)%L_Plus(0)
              DO l=1,Nloc
                ! switch to right hand system
                Uface(:,Nloc-p,q)=Uface(:,Nloc-p,q)+U_N(iElem)%U(:,p,l,q)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,p,q)=U_N(iElem)%U(:,p,q,0)*N_Inter(Nloc)%L_Plus(0)
              DO l=1,Nloc
                Uface(:,p,q)=Uface(:,p,q)+U_N(iElem)%U(:,p,q,l)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        END SELECT
#else /* for Gauss-Lobatto-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,q,p)=U_N(iElem)%U(:,0,p,q)
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          Uface(:,0:Nloc,0:Nloc)=U_N(iElem)%U(:,:,0,:)
        CASE(ZETA_MINUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,q,p)=U_N(iElem)%U(:,p,q,0)
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
          Uface(:,0:Nloc,0:Nloc)=U_N(iElem)%U(:,Nloc,:,:)
        CASE(ETA_PLUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,Nloc-p,q)=U_N(iElem)%U(:,p,Nloc,q)
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,Nloc
            DO p=0,Nloc
              Uface(:,p,q)=U_N(iElem)%U(:,p,q,Nloc)
            END DO ! p
          END DO ! q
        END SELECT
#endif
        ELSE ! no prolonge to face
          Uface=U_Surf_N(SideID)%U_master(:,:,:)
        END IF ! Prolong

        ! multiplied by surface element and Gauss Points
        area_loc     = SUM(N_SurfMesh(SideID)%SurfElem(:,:) * N_Inter(Nloc)%wGPSurf(:,:))
        integral_loc = SUM(Uface(1,0:Nloc,0:Nloc) * N_SurfMesh(SideID)%SurfElem(:,:) * N_Inter(Nloc)%wGPSurf(:,:))
        AverageElectricPotentialProc = AverageElectricPotentialProc + integral_loc/area_loc
    END IF ! flip =0
  END DO ! iSides
END DO ! iElems
!AverageElectricPotentialProc = AverageElectricPotentialProc / (1e-4 * 1.28e-2)

#if USE_MPI
  CALL MPI_ALLREDUCE(AverageElectricPotentialProc , AverageElectricPotential , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
  AverageElectricPotential = AverageElectricPotentialProc
#endif /*USE_MPI*/

! Get global average value
AverageElectricPotential = AverageElectricPotential / AverageElectricPotentialFaces

END SUBROUTINE CalculateAverageElectricPotential


SUBROUTINE GetAverageElectricPotentialPlane()
!===================================================================================================================================
!> Initializes Poynting vector integral variables and check every side: set "isPoyntingIntSide(SideID) = .TRUE." if a side coincides
!> with a defined Poynting vector integral plane.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars    ,ONLY: nSides,nElems
USE MOD_Mesh_Vars    ,ONLY: ElemToSide,N_SurfMesh
USE MOD_Analyze_Vars ,ONLY: isAverageElecPotSide,AverageElectricPotentialCoordErr,PosAverageElectricPotential
USE MOD_ReadInTools  ,ONLY: GETINT,GETREAL
USE MOD_Globals      ,ONLY: abort
#if USE_MPI
USE MOD_Globals
#endif
USE MOD_Analyze_Vars ,ONLY: AverageElectricPotential,AverageElectricPotentialFaces
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iSide, SideID, Nloc
INTEGER             :: nAverageElecPotSides
REAL                :: diff
INTEGER             :: p,q
!===================================================================================================================================

AverageElectricPotentialCoordErr = GETREAL('AvgPotential-Plane-Tolerance')
SWRITE(UNIT_stdOut,'(A)') ' GET PLANE TO CALCULATE AVERAGE ELECTRIC POTENTIAL ...'

! Initialize number of Poynting plane sides zero and set all sides to false
ALLOCATE(isAverageElecPotSide(1:nSides))
isAverageElecPotSide = .FALSE.

! Counter
AverageElectricPotential=0. ! initialize
nAverageElecPotSides = 0

! Loop over all elements
DO iElem=1,nElems
  ! Loop over all local sides
  DO iSide=1,6
    IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
      ASSOCIATE( AverageElecPotNormalDir1 => 2 ,&
                 AverageElecPotNormalDir2 => 3 ,&
                 AverageElecPotMainDir    => 1 )
        ! First search only planes with normal vector parallel to direction of "MainDir"
        IF((     N_SurfMesh(SideID)%NormVec(AverageElecPotNormalDir1,0,0)  < AverageElectricPotentialCoordErr) .AND. &
           (     N_SurfMesh(SideID)%NormVec(AverageElecPotNormalDir2,0,0)  < AverageElectricPotentialCoordErr) .AND. &
           ( ABS(N_SurfMesh(SideID)%NormVec(AverageElecPotMainDir   ,0,0)) > AverageElectricPotentialCoordErr))THEN
        ! Loop over all Points on Face
        Nloc = N_SurfMesh(SideID)%NSide
          DO q=0,Nloc
            DO p=0,Nloc
              diff = ABS(N_SurfMesh(SideID)%Face_xGP(AverageElecPotMainDir,p,q) - PosAverageElectricPotential)
              IF (diff < AverageElectricPotentialCoordErr) THEN
                IF (.NOT.isAverageElecPotSide(SideID)) THEN
                  nAverageElecPotSides = nAverageElecPotSides +1
                  isAverageElecPotSide(SideID) = .TRUE.
                END IF
              END IF ! diff < eps
            END DO !p
          END DO !q
        END IF
      END ASSOCIATE
    END IF ! flip = 0 master side
  END DO ! iSides
END DO !iElem=1,nElems

#if USE_MPI
CALL MPI_ALLREDUCE(nAverageElecPotSides , AverageElectricPotentialFaces , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_PICLAS , IERROR)
#else
AverageElectricPotentialFaces=nAverageElecPotSides
#endif /*USE_MPI*/
SWRITE(UNIT_stdOut,'(A,I10,A)') ' A total of',AverageElectricPotentialFaces,&
    ' surfaces for the average electric potential calculation are found.'
SWRITE(UNIT_stdOut,'(A)') ' ... AVERAGE ELECTRIC POTENTIAL INITIALIZATION DONE.'
#if USE_MPI
IF(MPIRoot)THEN
#endif /*USE_MPI*/
  IF(AverageElectricPotentialFaces.EQ.0)THEN
    SWRITE(UNIT_stdOut,*) 'ERROR with: PosAverageElectricPotential = ',PosAverageElectricPotential
    CALL abort(__STAMP__&
    ,'Found zero faces for averaging the electric potential. Please make sure \nthat the x-coordinate coincides with element'//&
    ' interfaces. Planes cutting through elements in currently not implemented.')
  END IF ! AverageElectricPotentialFaces.EQ.0
#if USE_MPI
END IF ! MPIRoot
#endif /*USE_MPI*/

END SUBROUTINE GetAverageElectricPotentialPlane


SUBROUTINE FinalizeAverageElectricPotential()
!===================================================================================================================================
! Finalize Poynting Integral
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars ,ONLY:isAverageElecPotSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! DEALLOCATE ALL
SDEALLOCATE(isAverageElecPotSide)
END SUBROUTINE FinalizeAverageElectricPotential


#if defined(PARTICLES)
SUBROUTINE CalculateElectricPotentialAndFieldBoundaryVDL()
!===================================================================================================================================
!> Calculation of the average electric potential with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is
!> missing ... ups
!>
!> 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!      interpolate the vector field E = (/Ex, Ey, Ez/) to the boundary face
!> 2.) Apply the normal vector: Uface(1,:,:)=DOT_PRODUCT(Uface(1:3,:,:),NormVec(1:3,:,:,SideID))*NormVec(1:3,:,:,SideID)
!> 3.) Apply the E-field correction factor
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals
#endif
USE MOD_Globals                ,ONLY: VECNORM
USE MOD_Mesh_Vars              ,ONLY: N_SurfMesh,SideToElem,nBCSides,N_SurfMesh,offSetElem,BC
USE MOD_DG_Vars                ,ONLY: U_N,N_DG_Mapping
USE MOD_Particle_Boundary_Vars ,ONLY: N_SurfVDL,PartBound,ElementThicknessVDLPerSide
USE MOD_ProlongToFace          ,ONLY: ProlongToFace_Side
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ElemID,SideID,ilocSide,ilocSide2,Nloc,iPartBound,p,q
REAL,ALLOCATABLE :: Uface(:,:,:),Eface(:,:,:)
REAL             :: MinOfElem,MinOfFace,MaxOfElem,MaxOfFace,Edir
!===================================================================================================================================
! 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!     interpolate the vector field E = (/Ex, Ey, Ez/) to the boundary face
DO SideID=1,nBCSides
  ! Get the local element index
  ElemID = SideToElem(S2E_ELEM_ID,SideID)
  ! Get local polynomial degree of the element
  Nloc   = N_DG_Mapping(2,ElemID+offSetElem)
  ! Get particle boundary index
  iPartBound = PartBound%MapToPartBC(BC(SideID))

  ! Skip sides that are not a VDL boundary (these sides are still in the list of sides)
  IF(PartBound%ThicknessVDL(iPartBound).GT.0.0)THEN

    ! Allocate Uface and Eface depending on the local polynomial degree
    ALLOCATE(Uface(  1,0:Nloc,0:Nloc))
    ALLOCATE(Eface(1:3,0:Nloc,0:Nloc))
    ! Get local side index
    ilocSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    ! Prolong-to-face depending on orientation in reference element
    CALL ProlongToFace_Side(3, Nloc, ilocSide, 0, U_N(ElemID)%E(:,:,:,:), Eface)
    ! Get opposing face for U because the Dirichlet value is present on the actual boundary face, which is already known
    SELECT CASE(ilocSide)
    CASE(XI_MINUS)
      ilocSide2 = XI_PLUS
    CASE(XI_PLUS)
      ilocSide2 = XI_MINUS
    CASE(ETA_MINUS)
      ilocSide2 = ETA_PLUS
    CASE(ETA_PLUS)
      ilocSide2 = ETA_MINUS
    CASE(ZETA_MINUS)
      ilocSide2 = ZETA_PLUS
    CASE(ZETA_PLUS)
      ilocSide2 = ZETA_MINUS
    END SELECT
    ! Defaulting to master sides with flip=0 in the following call. However, this only changes the sorting for Gauss-Lobatto
    CALL ProlongToFace_Side(1, Nloc, ilocSide2, 0, U_N(ElemID)%U(1,:,:,:), Uface)

    ! 2.) Apply the normal vector to get the normal electric field
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE( Ecorr => N_SurfVDL(SideID)%U(2:4,p,q), normal => N_SurfMesh(SideID)%NormVec(1:3,p,q), E => Eface(1:3,p,q))
          ! E_normal =  <E,normal>*normal
          Ecorr = DOT_PRODUCT(E,normal)*normal
        END ASSOCIATE
      END DO ! p
    END DO ! q

    ! Calculate the corrected E-field
    N_SurfVDL(SideID)%U(2:4,:,:) = N_SurfVDL(SideID)%U(2:4,:,:) * (ElementThicknessVDLPerSide(SideID)/PartBound%ThicknessVDL(iPartBound))

    ! Get Phi_F
    MinOfElem = MINVAL(U_N(ElemID)%U(1,:,:,:))
    MinOfFace = MINVAL(Uface)
    MaxOfElem = MAXVAL(U_N(ElemID)%U(1,:,:,:))
    MaxOfFace = MAXVAL(Uface)
    DO q=0,Nloc
      DO p=0,Nloc
        ASSOCIATE( E => N_SurfVDL(SideID)%U(2:4,p,q), normal => N_SurfMesh(SideID)%NormVec(1:3,p,q) )
          ! Normal vector points outwards on BC sides, hence, invert it
          Edir = DOT_PRODUCT(E,-normal)
          ! Check in which direction the E-field points
          IF(Edir.GT.0.0)THEN
            ! Get minimum of Phi as Phi_F (Phi_Max)
            N_SurfVDL(SideID)%U(5:5,p,q) = MIN(MinOfElem,MinOfFace)
            ! Reconstruct Phi_F from E
            N_SurfVDL(SideID)%U(1:1,p,q) = -VECNORM(E)*PartBound%ThicknessVDL(iPartBound)
            ! Reconstruct E from Phi_Max via E = Phi/d
            N_SurfVDL(SideID)%U(6:8,p,q) = -N_SurfVDL(SideID)%U(5,p,q)/PartBound%ThicknessVDL(iPartBound)*normal
          ELSE
            ! Get maximum of Phi as Phi_F (Phi_Max)
            N_SurfVDL(SideID)%U(5:5,p,q) = MAX(MaxOfElem,MaxOfFace)
            ! Reconstruct Phi_F from E
            N_SurfVDL(SideID)%U(1:1,p,q) =  VECNORM(E)*PartBound%ThicknessVDL(iPartBound)
            ! Reconstruct E from Phi_Max via E = Phi/d
            N_SurfVDL(SideID)%U(6:8,p,q) =  N_SurfVDL(SideID)%U(5,p,q)/PartBound%ThicknessVDL(iPartBound)*normal
          END IF ! Edir.GT.0
        END ASSOCIATE
      END DO ! p
    END DO ! q

    DEALLOCATE(Uface)
    DEALLOCATE(Eface)

  END IF ! PartBound%ThicknessVDL(iPartBound).GT.0.0
END DO ! SideID=1,nBCSides

END SUBROUTINE CalculateElectricPotentialAndFieldBoundaryVDL
#endif /*defined(PARTICLES)*/


SUBROUTINE CalculateElectricDisplacementCurrentSurface()
!===================================================================================================================================
!> Calculation of the average electric potential with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is
!> missing ... ups
!>
!> 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!      interpolate the vector field Dt = (/Dtx, Dty, Dtz/) to the boundary face
!> 2.) Apply the normal vector: Uface(1,:,:)=DOT_PRODUCT(Uface(1:3,:,:),NormVec(1:3,:,:,SideID))
!      Store result of dot product in first array index
!> 3.) Get BC index and EDC index and the mapping of the SideID boundary to the EDC boundary ID and store the integrated current
!> 4.) Communicate the integrated current values on each boundary to the MPI root process (the root outputs the values to .csv)
!===================================================================================================================================
! MODULES
#if USE_MPI
USE MOD_Globals
#endif
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh,SideToElem,nBCSides,N_SurfMesh,BC, offSetElem
USE MOD_Analyze_Vars       ,ONLY: EDC
USE MOD_Interpolation_Vars ,ONLY: N_Inter
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping
USE MOD_ProlongToFace      ,ONLY: ProlongToFace_Side
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ElemID,SideID,ilocSide,Nloc
REAL,ALLOCATABLE :: Eface(:,:,:)
INTEGER          :: iBC,iEDCBC
!REAL             :: SIP(0:PP_N,0:PP_N)
!REAL             :: AverageElectricPotentialProc
!REAL             :: area_loc,integral_loc
!===================================================================================================================================
! Nullify
EDC%Current = 0.

! 1.) Loop over all processor-local BC sides and therein find the local side ID which corresponds to the reference element and
!     interpolate the vector field Dt = (/Dtx, Dty, Dtz/) to the boundary face
DO SideID=1,nBCSides
  ! Get the local element index
  ElemID   = SideToElem(S2E_ELEM_ID,SideID)
  ! Get local polynomial degree of the element
  Nloc = N_DG_Mapping(2,ElemID+offSetElem)
  ! Allocate Eface depending on the local polynomial degree
  ALLOCATE(Eface(1:3,0:Nloc,0:Nloc))
  ! Get local side index
  ilocSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
  ! Prolong-to-face depending on orientation in reference element
  CALL ProlongToFace_Side(3, Nloc, ilocSide, 0, U_N(ElemID)%Dt, Eface)

  ! 2.) Apply the normal vector: Eface(1,:,:)=DOT_PRODUCT(Eface(1:3,:,:),NormVec(1:3,:,:,SideID))
  !     Store result of dot product in first array index
  Eface(1,:,:) =   Eface(1,:,:) * N_SurfMesh(SideID)%NormVec(1,:,:) &
                 + Eface(2,:,:) * N_SurfMesh(SideID)%NormVec(2,:,:) &
                 + Eface(3,:,:) * N_SurfMesh(SideID)%NormVec(3,:,:)

  ! 3.) Get BC index and EDC index and the mapping of the SideID boundary to the EDC boundary ID and store the integrated current
  iBC    = BC(SideID)
  iEDCBC = EDC%BCIDToEDCBCID(iBC)
  EDC%Current(iEDCBC) = EDC%Current(iEDCBC) + SUM(Eface(1,:,:) * N_SurfMesh(SideID)%SurfElem(:,:) * N_Inter(Nloc)%wGPSurf(:,:))

  DEALLOCATE(Eface)
END DO ! SideID=1,nBCSides

#if USE_MPI
! 4.) Communicate the integrated current values on each boundary to the MPI root process (the root outputs the values to .csv)
DO iEDCBC = 1, EDC%NBoundaries
  IF(EDC%COMM(iEDCBC)%UNICATOR.NE.MPI_COMM_NULL)THEN
    ASSOCIATE( Current => EDC%Current(iEDCBC), COMM => EDC%COMM(iEDCBC)%UNICATOR)
      IF(MPIroot)THEN
        CALL MPI_REDUCE(MPI_IN_PLACE , Current , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , COMM , IERROR)
      ELSE
        CALL MPI_REDUCE(Current      , 0       , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , 0 , COMM , IERROR)
      END IF ! myLeaderGroupRank.EQ.0
    END ASSOCIATE
  END IF ! EDC%COMM(iEDCBC)%UNICATOR.NE.MPI_COMM_NULL
END DO ! iEDCBC = 1, EDC%NBoundaries
#endif /*USE_MPI*/

END SUBROUTINE CalculateElectricDisplacementCurrentSurface
#endif /*USE_HDG*/

END MODULE MOD_AnalyzeField_HDG