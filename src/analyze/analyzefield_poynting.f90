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

MODULE MOD_AnalyzeField_Poynting
!===================================================================================================================================
! Contains the Poynting Vector Integral part for the power analysis of the field vector
!===================================================================================================================================
USE MOD_Globals, ONLY:UNIT_stdout
USE MOD_PreProc
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
#if !(USE_FV) || (USE_HDG)
#if (PP_nVar>=6)
PUBLIC:: CalcPoyntingIntegral,GetPoyntingIntPlane,FinalizePoyntingInt
#endif /*(PP_nVar>=6)*/
#endif /*no FV alone - !(USE_FV) || (USE_HDG)*/
!===================================================================================================================================

CONTAINS

#if !(USE_FV) || (USE_HDG)
#if (PP_nVar>=6)
SUBROUTINE CalcPoyntingIntegral(PoyntingIntegral)
!===================================================================================================================================
! Calculation of Poynting Integral with its own Prolong to face // check if Gauss-Lobatto or Gauss Points is used is missing ... ups
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems, N_SurfMesh, offSetElem
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_Analyze_Vars       ,ONLY: nPoyntingIntPlanes,isPoyntingIntSide,SideIDToPoyntingSide,PoyntingMainDir,Poynting
USE MOD_DG_Vars            ,ONLY: U_N,N_DG_Mapping
USE MOD_Globals_Vars       ,ONLY: smu0
USE MOD_Dielectric_Vars    ,ONLY: isDielectricFace,PoyntingUseMuR_Inv,DoDielectric
USE MOD_Globals
#if (PP_NodeType==1) /* for Gauss-points*/
USE MOD_Interpolation_Vars ,ONLY: N_inter
#endif /*(PP_NodeType==1)*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: PoyntingIntegral(1:nPoyntingIntPlanes)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,SideID,ilocSide,iPoyntingSide,Nloc
INTEGER          :: p,q
#if (PP_NodeType==1) /* for Gauss-points*/
INTEGER          :: l
#endif /*(PP_NodeType==1)*/
#if USE_MPI
REAL             :: SumSabs(nPoyntingIntPlanes)
#endif
!===================================================================================================================================
PoyntingIntegral = 0.

iPoyntingSide = 0
DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  Do ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isPoyntingIntSide(SideID)) CYCLE ! Skip other sides
      ! calculate Poynting vector
      iPoyntingSide = iPoyntingSide + 1
      Poynting(iPoyntingSide)%S = 0.
#if (PP_NodeType==1) /* for Gauss-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,q,p)=U_N(iElem)%U(:,0,p,q)*N_Inter(Nloc)%L_Minus(0)
            DO l=1,Nloc
                ! switch to right hand system
              Poynting(iPoyntingSide)%Uface(:,q,p)=Poynting(iPoyntingSide)%Uface(:,q,p)+U_N(iElem)%U(:,l,p,q)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,p,q)=U_N(iElem)%U(:,p,0,q)*N_Inter(Nloc)%L_Minus(0)
            DO l=1,Nloc
              Poynting(iPoyntingSide)%Uface(:,p,q)=Poynting(iPoyntingSide)%Uface(:,p,q)+U_N(iElem)%U(:,p,l,q)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_MINUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,q,p)=U_N(iElem)%U(:,p,q,0)*N_Inter(Nloc)%L_Minus(0)
            DO l=1,Nloc
                ! switch to right hand system
              Poynting(iPoyntingSide)%Uface(:,q,p)=Poynting(iPoyntingSide)%Uface(:,q,p)+U_N(iElem)%U(:,p,q,l)*N_Inter(Nloc)%L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! qfirst stuff
        CASE(XI_PLUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,p,q)=U_N(iElem)%U(:,0,p,q)*N_Inter(Nloc)%L_Plus(0)
            DO l=1,Nloc
              Poynting(iPoyntingSide)%Uface(:,p,q)=Poynting(iPoyntingSide)%Uface(:,p,q)+U_N(iElem)%U(:,l,p,q)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_PLUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,Nloc-p,q)=U_N(iElem)%U(:,p,0,q)*N_Inter(Nloc)%L_Plus(0)
            DO l=1,Nloc
                ! switch to right hand system
              Poynting(iPoyntingSide)%Uface(:,Nloc-p,q)=Poynting(iPoyntingSide)%Uface(:,Nloc-p,q)+U_N(iElem)%U(:,p,l,q)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,p,q)=U_N(iElem)%U(:,p,q,0)*N_Inter(Nloc)%L_Plus(0)
            DO l=1,Nloc
              Poynting(iPoyntingSide)%Uface(:,p,q)=Poynting(iPoyntingSide)%Uface(:,p,q)+U_N(iElem)%U(:,p,q,l)*N_Inter(Nloc)%L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        END SELECT
#else /* for Gauss-Lobatto-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,q,p)=U_N(iElem)%U(:,0,p,q)
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
        Poynting(iPoyntingSide)%Uface(:,:,:)=U_N(iElem)%U(:,:,0,:)
        CASE(ZETA_MINUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,q,p)=U_N(iElem)%U(:,p,q,0)
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
        Poynting(iPoyntingSide)%Uface(:,:,:)=U_N(iElem)%U(:,Nloc,:,:)
        CASE(ETA_PLUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,Nloc-p,q)=U_N(iElem)%U(:,p,Nloc,q)
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
        DO q=0,Nloc
          DO p=0,Nloc
            Poynting(iPoyntingSide)%Uface(:,p,q)=U_N(iElem)%U(:,p,q,Nloc)
            END DO ! p
          END DO ! q
        END SELECT
#endif
        ! check if dielectric regions are involved
        IF(DoDielectric)THEN
          IF(PoyntingUseMuR_Inv.AND.isDielectricFace(SideID))THEN
          CALL abort(__STAMP__,'CALL PoyntingVectorDielectric() not implemented')
          !CALL PoyntingVectorDielectric(Uface(:,:,:),S(:,:,:,iPoyntingSide),Dielectric_MuR_Master_inv(0:Nloc,0:Nloc,SideID))
          ELSE
          CALL PoyntingVector(Nloc,Poynting(iPoyntingSide)%Uface(:,:,:),Poynting(iPoyntingSide)%S(:,:,:))
          END IF
        ELSE
        CALL PoyntingVector(Nloc,Poynting(iPoyntingSide)%Uface(:,:,:),Poynting(iPoyntingSide)%S(:,:,:))
        END IF

      ASSOCIATE( NormVec     => N_SurfMesh(SideID)%NormVec(:,0,0)       ,& ! Use first vector entry and assume planar face
                 SIP         => Poynting(iPoyntingSide)%SIP(:,:)    )
                 !S           => Poynting(iPoyntingSide)%S(:,:,:)        )! ,&
                 !SurfElemwGP => Poynting(iPoyntingSide)%SurfElemwGP(:,:) )
        IF(NormVec(PoyntingMainDir).GT.0.0)THEN
          SIP(:,:) =   Poynting(iPoyntingSide)%S(1,:,:) * NormVec(1) &
                     + Poynting(iPoyntingSide)%S(2,:,:) * NormVec(2) &
                     + Poynting(iPoyntingSide)%S(3,:,:) * NormVec(3)
        ELSE ! NormVec(PoyntingMainDir).LT.0
          SIP(:,:) = - Poynting(iPoyntingSide)%S(1,:,:) * NormVec(1) &
                     - Poynting(iPoyntingSide)%S(2,:,:) * NormVec(2) &
                     - Poynting(iPoyntingSide)%S(3,:,:) * NormVec(3)
        END IF ! NormVec(PoyntingMainDir,:,:,iPoyntingSide)
      END ASSOCIATE

        ! multiplied by surface element and  Gauss Points
      Poynting(iPoyntingSide)%SIP(:,:) = Poynting(iPoyntingSide)%SIP(:,:) * Poynting(iPoyntingSide)%SurfElemwGP(:,:)

        ! total flux through each plane
      PoyntingIntegral(SideIDToPoyntingSide(SideID)) = PoyntingIntegral(SideIDToPoyntingSide(SideID)) + smu0*SUM(Poynting(iPoyntingSide)%SIP(:,:))
    END IF ! flip =0
  END DO ! iSides
END DO ! iElems

#if USE_MPI
  CALL MPI_REDUCE(PoyntingIntegral(:) , SumSabs(:) , nPoyntingIntPlanes , MPI_DOUBLE_PRECISION ,MPI_SUM, 0, MPI_COMM_PICLAS,IERROR)
  PoyntingIntegral(:) = SumSabs(:)
#endif /*USE_MPI*/

END SUBROUTINE CalcPoyntingIntegral


PPURE SUBROUTINE PoyntingVector(Nloc,U,S)
!===================================================================================================================================
!> Calculate the Poynting Vector on a certain face for vacuum properties
!>
!> ATTENTION: permeability is not applied here due to performance gain
!> Definition: S = E x H = 1/mu0 * ( E x H )
!> Here      : S = E x B (i.e. mu0 is applied later)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: Nloc
REAL,INTENT(IN)       :: U(PP_nVar,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: S(1:3,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
!===================================================================================================================================

! calculate the Poynting vector at each node, additionally the abs of the Poynting vector only based on E
DO p = 0,Nloc
  DO q = 0,Nloc
    S(1,p,q)  =  U(2,p,q)*U(6,p,q) - U(3,p,q)*U(5,p,q)
    S(2,p,q)  = -U(1,p,q)*U(6,p,q) + U(3,p,q)*U(4,p,q)
    S(3,p,q)  =  U(1,p,q)*U(5,p,q) - U(2,p,q)*U(4,p,q)
  END DO ! q - Nloc
END DO  ! p - Nloc

END SUBROUTINE PoyntingVector


! PPURE SUBROUTINE PoyntingVectorDielectric(Uface_in,Sloc,mu_r_inv)
! !===================================================================================================================================
! !> Calculate the Poynting Vector on a certain face for dielectric properties (consider mu_r here, but not mu0)
! !>
! !> ATTENTION: permeability is not applied here due to performance gain
! !> Definition: S = E x H = 1/(mu_r*mu_0) * ( E x H )
! !> Here      : S = 1/mu_r * E x B (i.e. mu0 is applied later)
! !===================================================================================================================================
! ! MODULES
! ! IMPLICIT VARIABLE HANDLING
! IMPLICIT NONE
! !----------------------------------------------------------------------------------------------------------------------------------
! ! INPUT VARIABLES
! REAL,INTENT(IN)       :: Uface_in(PP_nVar,0:PP_N,0:PP_N)
! REAL,INTENT(IN)       :: mu_r_inv(0:PP_N,0:PP_N)         ! 1/mu_r for every face DOF (may vary on face depending on position)
! !                                                        ! (isotropic property for permittivity)
! !----------------------------------------------------------------------------------------------------------------------------------
! ! OUTPUT VARIABLES
! REAL,INTENT(OUT)      :: Sloc(1:3,0:PP_N,0:PP_N)
! !----------------------------------------------------------------------------------------------------------------------------------
! ! LOCAL VARIABLES
! INTEGER               :: p,q
! !===================================================================================================================================
!
! ! calculate the Poynting vector at each node, additionally the abs of the Poynting vector only based on E
! DO p = 0,PP_N
!   DO q = 0,PP_N
!     Sloc(1,p,q)  = (  Uface_in(2,p,q)*Uface_in(6,p,q) - Uface_in(3,p,q)*Uface_in(5,p,q) ) * mu_r_inv(p,q)
!     Sloc(2,p,q)  = ( -Uface_in(1,p,q)*Uface_in(6,p,q) + Uface_in(3,p,q)*Uface_in(4,p,q) ) * mu_r_inv(p,q)
!     Sloc(3,p,q)  = (  Uface_in(1,p,q)*Uface_in(5,p,q) - Uface_in(2,p,q)*Uface_in(4,p,q) ) * mu_r_inv(p,q)
!   END DO ! q - PP_N
! END DO  ! p - PP_N
!
! END SUBROUTINE PoyntingVectorDielectric


SUBROUTINE GetPoyntingIntPlane()
!===================================================================================================================================
!> Initializes Poynting vector integral variables and check every side: set "isPoyntingIntSide(SideID) = .TRUE." if a side coincides
!> with a defined Poynting vector integral plane.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars       ,ONLY: nSides,nElems, offSetElem
USE MOD_Mesh_Vars       ,ONLY: ElemToSide,N_SurfMesh
USE MOD_Analyze_Vars    ,ONLY: PoyntingIntCoordErr,nPoyntingIntPlanes,PosPoyntingInt!,S,STEM
USE MOD_Analyze_Vars    ,ONLY: isPoyntingIntSide,SideIDToPoyntingSide,PoyntingMainDir
USE MOD_ReadInTools     ,ONLY: GETINT,GETREAL
USE MOD_Dielectric_Vars ,ONLY: DoDielectric,nDielectricElems,DielectricVol,ElemToDielectric,isDielectricInterFace
USE MOD_Dielectric_Vars ,ONLY: isDielectricFace,PoyntingUseMuR_Inv
USE MOD_Globals         ,ONLY: abort
#if USE_MPI
USE MOD_Globals
#else
USE MOD_Globals         ,ONLY: CollectiveStop
#endif
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iSide, iPlane, SideID, Nloc
INTEGER,ALLOCATABLE :: nFaces(:)
REAL                :: diff
INTEGER             :: p,q
CHARACTER(LEN=32)   :: index_plane
INTEGER,ALLOCATABLE :: sumFaces(:)
INTEGER             :: sumAllfaces
LOGICAL             :: CheckDielectricSides
INTEGER             :: PoyntingNormalDir1,PoyntingNormalDir2
INTEGER             :: nPoyntingIntSides    !< Sides for the calculation of the Poynting vector integral
!===================================================================================================================================

LBWRITE(UNIT_stdOut,'(A)') ' GET PLANES TO CALCULATE POYNTING VECTOR INTEGRAL ...'

! Initialize number of Poynting plane sides zero and set all sides to false
nPoyntingIntSides=0
ALLOCATE(isPoyntingIntSide(1:nSides))
isPoyntingIntSide = .FALSE.

! Get the number of Poynting planes and coordinates
nPoyntingIntPlanes = GETINT('PoyntingVecInt-Planes')
PoyntingMainDir = GETINT('PoyntingMainDir') ! default "3" is z-direction
SELECT CASE (PoyntingMainDir)
  CASE (1) ! poynting vector integral in x-direction
    PoyntingNormalDir1=2
    PoyntingNormalDir2=3
  CASE (2) ! poynting vector integral in y-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=3
  CASE (3) ! poynting vector integral in z-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=2
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,'Poynting vector itnegral currently only in x,y,z!')
END SELECT
ALLOCATE(PosPoyntingInt(nPoyntingIntPlanes))
ALLOCATE(SideIDToPoyntingSide(nSides))
ALLOCATE(nFaces(nPoyntingIntPlanes))
SideIDToPoyntingSide = -1
nFaces(:) = 0

! Get z-coordinates and factors for every Poynting plane
DO iPlane=1,nPoyntingIntPlanes
 WRITE(UNIT=index_plane,FMT='(I2.2)') iPlane
 SELECT CASE (PoyntingMainDir)
    CASE (1)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-x-coord')
    CASE (2)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-y-coord')
    CASE (3)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-z-coord')
  END SELECT
END DO
PoyntingIntCoordErr=GETREAL('Plane-Tolerance')

! Dielectric Sides:
! 1.) check if a dielectric region (only permeability, NOT permittivity is important) coincides with a Poynting vector
!     integral plane. Dielectric interfaces with mu_r .NE. 1.0 cannot compute a Poynting vector because of the jump in material
!     parameter of mu_r
CheckDielectricSides=.FALSE.
IF(DoDielectric)THEN
  DO iElem = 1, nDielectricElems
    IF(ANY(ABS(DielectricVol(iElem)%DielectricMu(:,:,:)-1.0).GT.0.0))THEN
    CheckDielectricSides=.TRUE.
      EXIT
  END IF
  END DO
END IF

! 2.) for dielectric sides (NOT interface sides between dielectric and some other region), determine mu_r on face for Poynting vector
PoyntingUseMuR_Inv=.FALSE.

! Loop over all planes
DO iPlane = 1, nPoyntingIntPlanes
  ! Loop over all elements
  DO iElem=1,nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ! Loop over all local sides
    DO iSide=1,6
      IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
        SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
        ! First search only planes with normal vector parallel to direction of "MainDir"
        IF((     N_SurfMesh(SideID)%NormVec(PoyntingNormalDir1,0,0)  < PoyntingIntCoordErr) .AND. &
           (     N_SurfMesh(SideID)%NormVec(PoyntingNormalDir2,0,0)  < PoyntingIntCoordErr) .AND. &
           ( ABS(N_SurfMesh(SideID)%NormVec(PoyntingMainDir   ,0,0)) > PoyntingIntCoordErr))THEN
        ! Loop over all Points on Face
          DO q=0,Nloc
            DO p=0,Nloc
              diff = ABS(N_SurfMesh(SideID)%Face_xGP(PoyntingMainDir,p,q) - PosPoyntingInt(iPlane))
              IF (diff < PoyntingIntCoordErr) THEN
                IF (.NOT.isPoyntingIntSide(SideID)) THEN
                  nPoyntingIntSides = nPoyntingIntSides +1
                  SideIDToPoyntingSide(SideID) = iPlane
                  isPoyntingIntSide(SideID) = .TRUE.
                  nFaces(iPlane) = nFaces(iPlane) + 1

                  ! Dielectric sides
                  IF(CheckDielectricSides)THEN
                    ! 1.) Check for illegal sides in dielectrics: mu_r != 1.0 on dielectric interface
                    IF(isDielectricInterFace(SideID))THEN
                      IF(ANY(ABS(DielectricVol(ElemToDielectric(iElem))%DielectricMu(:,:,:)-1.0).GT.0.0))THEN
                        ! If the Poynting vector integral SideID additionally is a dielectric interface between a dielectric region
                        ! with a permittivity and vacuum, then mu_r might be unequal to 1.0 on the interface and the calculation of
                        ! the Poynting vector is not implemented for this case
                        IPWRITE(UNIT_stdOut,*) " "
                        IPWRITE(UNIT_stdOut,*) "Found illegal Poyting plane side. SideID= ",SideID,&
                            " z-coordinate= ",PosPoyntingInt(iPlane)
                        CALL abort(__STAMP__&
                            ,'GetPoyntingIntPlane: Found SideID for Poynting vector integral which is attached to an element'//&
                            ' within which the dielectric permittivity mu_r is not euqal to 1.0 everywhere. The value could be'//&
                            ' unequal to 1.0 on the interface and this is not implemented. TODO: determine mu_r on interface,'//&
                            ' communicate it via MPI (do not forget Mortar sides) and calculate the Poynting vector on that'//&
                            ' interface via some method.')
                      END IF
                    END IF

                    ! 2.) Check for legal sides in dielectrics: mu_r != 1.0 within dielectric region
                    IF(isDielectricFace(SideID))THEN
                      !IPWRITE(UNIT_stdOut,*) "found dielectric face: ",SideID,"z= ",PosPoyntingInt(iPlane)
                      PoyntingUseMuR_Inv=.TRUE.
                    END IF
                  END IF

                END IF
              END IF ! diff < eps
            END DO !p
          END DO !q
        END IF ! n parallel gyrotron axis
      END IF ! flip = 0 master side
    END DO ! iSides
  END DO !iElem=1,nElems
END DO ! iPlanes

! Dielectric sides:
#if USE_MPI
! Send info to ALL MPI ranks:
! TODO: If 1/mu_r is never needed on master AND slave procs, this routine can be adjusted so that only master procs determine the
! prolonged values of mu_r and no MPI information has to be sent. The master side cannot currently be outside of the dielectric
! region (e.g. in vacuum) because that is not allowed. If this would be allowed that MPI rank would need the information of the
! prolonged dielectric material properties from the slave side
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,PoyntingUseMuR_Inv,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_PICLAS,iError)
#endif
! Determine mu_r on faces within a dielectric region for calculating the Poynting vector and communicate the
! prolonged values via MPI
#if (PP_nVar>=6)
IF(PoyntingUseMuR_Inv) CALL SetDielectricFaceProfileForPoynting()
#endif /*(PP_nVar>=6)*/

ALLOCATE(sumFaces(nPoyntingIntPlanes))
#if USE_MPI
sumFaces=0
sumAllFaces=0
  CALL MPI_REDUCE(nFaces , sumFaces , nPoyntingIntPlanes , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  !nFaces(:) = sumFaces(:)
  CALL MPI_REDUCE(nPoyntingIntSides , sumAllFaces , 1 , MPI_INTEGER, MPI_SUM,0, MPI_COMM_PICLAS, IERROR)
  !nPoyntingIntSides = sumAllFaces
#else
sumFaces=nFaces
sumAllFaces=nPoyntingIntSides
#endif /*USE_MPI*/

DO iPlane= 1, nPoyntingIntPlanes
  LBWRITE(UNIT_stdOut,'(A,I2,A,I10,A)') 'Processed plane no.: ',iPlane,'. Found ',sumFaces(iPlane),' surfaces.'
END DO
LBWRITE(UNIT_stdOut,'(A,I10,A)') 'A total of',sumAllFaces, &
                        ' surfaces for the poynting vector integral calculation are found.'

!ALLOCATE(S    (1:3,0:PP_N,0:PP_N,1:nPoyntingIntSides) , &
         !STEM     (0:PP_N,0:PP_N,1:nPoyntingIntSides)  )
! Allocate the Poynting vector containers S, SIP and Uface
CALL AllocatePoyntingVector()

LBWRITE(UNIT_stdOut,'(A)') ' ... POYNTING VECTOR INTEGRAL INITIALIZATION DONE.'

END SUBROUTINE GetPoyntingIntPlane


!===================================================================================================================================
!> Allocate the Poynting vector containers S, SIP and Uface on the element-local polynomial degree Nloc
!===================================================================================================================================
SUBROUTINE AllocatePoyntingVector()
! MODULES
USE MOD_Mesh_Vars          ,ONLY: nElems, N_SurfMesh, offSetElem
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_Analyze_Vars       ,ONLY: Poynting,isPoyntingIntSide
USE MOD_Interpolation_Vars ,ONLY: NMax
USE MOD_DG_Vars            ,ONLY: DG_Elems_master,DG_Elems_slave,N_DG_Mapping
USE MOD_Interpolation_Vars ,ONLY: N_Inter,PREF_VDM
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iElem,SideID,ilocSide,iPoyntingSide,Nloc,NSideMax,NbrOfPoyntingFace
REAL,DIMENSION(0:NMax,0:NMax) :: SurfElem,tmp
!===================================================================================================================================

! Loop over all elements
NbrOfPoyntingFace = 0
DO iElem = 1, nElems
  ! Loop over the 6 local sides
  DO ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem).EQ.0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isPoyntingIntSide(SideID)) CYCLE ! Skip other sides
      NbrOfPoyntingFace = NbrOfPoyntingFace + 1
    END IF ! ElemToSide(E2S_FLIP,ilocSide,iElem).EQ.0
  END DO ! ilocSide = 1, 6
END DO ! iElem = 1, nElems
ALLOCATE(Poynting(1:NbrOfPoyntingFace))

! Exit is not faces are found: Can this happen?
IF(NbrOfPoyntingFace.EQ.0) RETURN

iPoyntingSide = 0
! Loop over all elements
DO iElem = 1, nElems
  ! Get element-local polynomial degree
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ! Loop over the 6 local sides
  DO ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem).EQ.0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isPoyntingIntSide(SideID)) CYCLE ! Skip other sides

      ! Get index of Poynting side
      iPoyntingSide = iPoyntingSide + 1

      ! Allocate Poynting vector
      ALLOCATE(Poynting(iPoyntingSide)%Uface(PP_nVar,0:Nloc,0:Nloc))
      ALLOCATE(Poynting(iPoyntingSide)%S(1:3,0:Nloc,0:Nloc))
      ALLOCATE(Poynting(iPoyntingSide)%SIP(0:Nloc,0:Nloc))
      ALLOCATE(Poynting(iPoyntingSide)%SurfElemwGP(0:Nloc,0:Nloc))

      ! Map surface elem to Nloc, multiply with wGPSurf and store in SurfElemwGP
      ! SurfElem is built on N = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
      NSideMax = MAX(DG_Elems_master(SideID),DG_Elems_slave(SideID))
      ! TODO NSideMin - What to do here?
      ! Map SurfElem from NSideMax to NSideMin
      IF(Nloc.EQ.NSideMax)THEN ! N is equal
        SurfElem(0:Nloc,0:Nloc) = N_SurfMesh(SideID)%SurfElem(0:NSideMax,0:NSideMax)
      ELSEIF(Nloc.GT.NSideMax)THEN ! N increases: Simply interpolate the lower polynomial degree solution
        CALL ChangeBasis2D(1, NSideMax, Nloc, PREF_VDM(NSideMax, Nloc)%Vdm, N_SurfMesh(SideID)%SurfElem(0:NSideMax,0:NSideMax), SurfElem(0:Nloc,0:Nloc))
      ELSE ! N reduces: This requires an intermediate modal basis
        ! Switch to Legendre basis
        CALL ChangeBasis2D(1, NSideMax, NSideMax, N_Inter(NSideMax)%sVdm_Leg, N_SurfMesh(SideID)%SurfElem(0:NSideMax,0:NSideMax), tmp(0:NSideMax,0:NSideMax))
        ! Switch back to nodal basis but cut-off the higher-order DOFs
        CALL ChangeBasis2D(1, Nloc, Nloc, N_Inter(Nloc)%Vdm_Leg, tmp(0:Nloc,0:Nloc), SurfElem(0:Nloc,0:Nloc))
      END IF ! Nloc.EQ.NSideMax
      ! Calculate SurfElem * N_Inter
      Poynting(iPoyntingSide)%SurfElemwGP(:,:) = SurfElem(0:Nloc,0:Nloc) * N_Inter(Nloc)%wGPSurf(:,:)

    END IF ! ElemToSide(E2S_FLIP,ilocSide,iElem).EQ.0
  END DO ! ilocSide = 1, 6
END DO ! iElem = 1, nElems

END SUBROUTINE AllocatePoyntingVector


SUBROUTINE FinalizePoyntingInt()
!===================================================================================================================================
! Finalize Poynting Integral
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars ,ONLY:PosPoyntingInt,isPoyntingIntSide,SideIDToPoyntingSide,Poynting
!USE MOD_Analyze_Vars ,ONLY:PosPoyntingInt,S,STEM,isPoyntingIntSide,SideIDToPoyntingSide
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
SDEALLOCATE(isPoyntingIntSide)
SDEALLOCATE(PosPoyntingInt)
SDEALLOCATE(SideIDToPoyntingSide)
SDEALLOCATE(Poynting)

END SUBROUTINE FinalizePoyntingInt


SUBROUTINE SetDielectricFaceProfileForPoynting()
!===================================================================================================================================
!> THIS ROUTINE IS ONLY CALLED FOR THE POYNTING VECTOR INTEGRAL CALCULATION ON INITIALIZATION
!>
!> Set the dielectric factor 1./MuR for each face DOF in the array "Dielectric_MuR_Master_inv" (needed for S = E X H calculation).
!> Only the array "Dielectric_MuR_Master_inv" is used in the Riemann solver, as only the master calculates the flux array
!> (maybe slave information is used in the future)
!>
!> Note:
!> for MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!> threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!> is not allowed)
!> re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!> information)
!> This could be overcome by using template subroutines .t90 (see FlexiOS)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: Dielectric_MuR_Master_inv,Dielectric_MuR_Slave_inv
USE MOD_Dielectric_Vars ,ONLY: isDielectricElem,ElemToDielectric,DielectricVol
USE MOD_Mesh_Vars       ,ONLY: nSides
USE MOD_ProlongToFace   ,ONLY: ProlongToFace
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI             ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
!USE MOD_FillMortar      ,ONLY: U_Mortar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE,Dielectric_dummy_Master2S
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,1:nSides)           :: Dielectric_dummy_Master
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,1:nSides)           :: Dielectric_dummy_Slave
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems) :: Dielectric_dummy_elem
#if USE_MPI
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Master2
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides)                 :: Dielectric_dummy_Slave2
INTEGER                                                  :: I,J,iSide
#endif /*USE_MPI*/
INTEGER                                                  :: iElem
!===================================================================================================================================
! General workflow:
! 1.  Initialize dummy arrays for Elem/Face
! 2.  Fill dummy element values for non-Dielectric sides
! 3.  Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar)
! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!     is not allowed)
!     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!     information)
! 5.  Send/Receive MPI data
! 6.  Allocate the actually needed arrays containing the dielectric material information on the sides
! 7.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)

! 1.  Initialize dummy arrays for Elem/Face
Dielectric_dummy_elem    = -1.
Dielectric_dummy_Master  = -1.
Dielectric_dummy_Slave   = -1.

! 2.  Fill dummy element values for non-Dielectric sides
DO iElem=1,PP_nElems
  IF(isDielectricElem(iElem))THEN
    ! set only the first dimension to 1./MuR (the rest are dummies)
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=1.0 / DielectricVol(ElemToDielectric(iElem))%DielectricMu(0:PP_N,0:PP_N,0:PP_N)
  ELSE
    Dielectric_dummy_elem(1,0:PP_N,0:PP_N,0:PP_N,(iElem))=1.0
  END IF
END DO

!3.   Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar)
CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.FALSE.)
CALL abort(__STAMP__,'not implemented')
!CALL U_Mortar(Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.FALSE.)
#if USE_MPI
  CALL ProlongToFace(Dielectric_dummy_elem,Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.TRUE.)
  !CALL U_Mortar(Dielectric_dummy_Master,Dielectric_dummy_Slave,doMPISides=.TRUE.)

  ! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
  !     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
  !     is not allowed)
  !     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
  !     information)
  Dielectric_dummy_Master2 = 0.
  Dielectric_dummy_Slave2  = 0.
  DO I=0,PP_N
    DO J=0,PP_N
      DO iSide=1,nSides
        Dielectric_dummy_Master2(1,I,J,iSide)=Dielectric_dummy_Master(1,I,J,iSide)
        Dielectric_dummy_Slave2 (1,I,J,iSide)=Dielectric_dummy_Slave( 1,I,J,iSide)
      END DO
    END DO
  END DO

  ! 5.  Send Slave Dielectric info (real array with dimension (N+1)*(N+1)) to Master procs
  CALL StartReceiveMPIData(1,Dielectric_dummy_Slave2 ,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
  CALL StartSendMPIData(   1,Dielectric_dummy_Slave2 ,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR

  ! Send Master Dielectric info (real array with dimension (N+1)*(N+1)) to Slave procs
  CALL StartReceiveMPIData(1,Dielectric_dummy_Master2,1,nSides ,RecRequest_U ,SendID=1) ! Receive YOUR
  CALL StartSendMPIData(   1,Dielectric_dummy_Master2,1,nSides,SendRequest_U ,SendID=1) ! Send MINE

  CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE - receive YOUR
  CALL FinishExchangeMPIData(SendRequest_U, RecRequest_U ,SendID=1) !Send YOUR - receive MINE
#endif /*USE_MPI*/

! 6.  Allocate the actually needed arrays containing the dielectric material information on the sides
ALLOCATE(Dielectric_MuR_Master_inv(0:PP_N,0:PP_N,1:nSides))
ALLOCATE(Dielectric_MuR_Slave_inv( 0:PP_N,0:PP_N,1:nSides))


! 7.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
#if USE_MPI
  Dielectric_MuR_Master_inv=Dielectric_dummy_Master2(1,0:PP_N,0:PP_N,1:nSides)
  Dielectric_MuR_Slave_inv =Dielectric_dummy_Slave2( 1,0:PP_N,0:PP_N,1:nSides)
#else
  Dielectric_MuR_Master_inv=Dielectric_dummy_Master(1,0:PP_N,0:PP_N,1:nSides)
  Dielectric_MuR_Slave_inv =Dielectric_dummy_Slave( 1,0:PP_N,0:PP_N,1:nSides)
#endif /*USE_MPI*/

! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)
IF(MINVAL(Dielectric_MuR_Master_inv).LE.0.0)THEN
  CALL abort(&
  __STAMP__&
  ,'Dielectric material values for Riemann solver not correctly determined. MINVAL(Dielectric_MuR_Master_inv)=',&
  RealInfoOpt=MINVAL(Dielectric_MuR_Master_inv))
END IF
END SUBROUTINE SetDielectricFaceProfileForPoynting
#endif /*PP_nVar>=6*/
#endif /*no FV alone - !(USE_FV) || (USE_HDG)*/

END MODULE MOD_AnalyzeField_Poynting