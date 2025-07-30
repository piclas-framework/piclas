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

MODULE MOD_GetBoundaryFlux
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700)) && !(USE_HDG)
!===================================================================================================================================
! Contains FillBoundary (which depends on the considered equation)
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
PUBLIC::GetBoundaryFlux
PUBLIC:: InitBC,FinalizeBC
!===================================================================================================================================

CONTAINS


SUBROUTINE InitBC()
!===================================================================================================================================
! Initialize boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
   CALL abort(__STAMP__,"InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
  ! should not be required || example for MaxBCState
  !IF((locType.NE.22).AND.locType.NE.3) MaxBCState = MAX(MaxBCState,locState)
  !IF((locType.EQ.4).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No temperature (refstate) defined for BC_TYPE',locType)
  !IF((locType.EQ.23).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow Mach number in refstate (x,Ma,x,x,x) defined for BC_TYPE',locType)
  !IF((locType.EQ.24).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow pressure in refstate defined for BC_TYPE',locType)
  !IF((locType.EQ.25).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No outflow pressure in refstate defined for BC_TYPE',locType)
  !IF((locType.EQ.27).AND.(locState.LT.1))&
  !  CALL abort(__STAMP__,&
  !             'No inflow refstate (dens,v1,v2,v3,pressure) in refstate defined for BC_TYPE',locType)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

! Sanity check for BCs
!IF(MaxBCState.GT.nRefState)&
!  CALL abort(__STAMP__,&
!    'ERROR: Boundary RefState not defined! (MaxBCState,nRefState):',MaxBCState,REAL(nRefState))


! ! Initialize State File Boundary condition
! DO i=1,nBCs
!   locType =BoundaryType(i,BC_TYPE)
!   IF(locType.EQ.20)THEN
!     ! Allocate buffer array to store temp data for all BC sides
!     ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_N,nBCSides))
!     BCData=0.
!     CALL ReadBCFlow(BCStateFile)
!     !CALL abort(__STAMP__,&
!     !     'no BC defined in maxwell/getboundaryflux.f90!')
!     EXIT
!   END IF
! END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC


SUBROUTINE GetBoundaryFlux(SideID, t,tDeriv, Nloc, Flux, U_Minus, NormVec, TangVec1, TangVec2, BCFace_xGP)
!===================================================================================================================================
! Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
! BCType: 1...periodic, 2...exact BC
! Attention 1: this is only a tensor of local values U_Face and has to be stored into the right U_Left or U_Right in
!              SUBROUTINE CalcSurfInt
! Attention 2: U_FacePeriodic is only needed in the case of periodic boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Globals         ,ONLY: Abort,CROSS
USE MOD_PreProc
USE MOD_Riemann         ,ONLY: RiemannVacuum,RiemannPML
USE MOD_Riemann         ,ONLY: RiemannDielectric
USE MOD_Equation        ,ONLY: ExactFunc
USE MOD_Globals_Vars    ,ONLY: c,c_inv
USE MOD_Mesh_Vars       ,ONLY: BoundaryType,BC
USE MOD_PML_Vars        ,ONLY: PMLnVar, DoPML
USE MOD_Interfaces_Vars ,ONLY: InterfaceRiemann
USE MOD_Dielectric_vars ,ONLY: DielectricSurf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID                             !< ID of current side
REAL,INTENT(IN)                      :: t                                  !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: tDeriv
INTEGER,INTENT(IN)                   :: Nloc                               !< polynomial degree
REAL,INTENT(IN)                      :: U_Minus(     PP_nVar,0:Nloc,0:Nloc)
REAL,INTENT(IN)                      :: NormVec(           3,0:Nloc,0:Nloc)
REAL,INTENT(IN),OPTIONAL             :: TangVec1(          3,0:Nloc,0:Nloc)
REAL,INTENT(IN),OPTIONAL             :: TangVec2(          3,0:Nloc,0:Nloc)
REAL,INTENT(IN)                      :: BCFace_xGP(        3,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                     :: Flux(PP_nVar+PMLnVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: p,q
INTEGER                              :: BCType,BCState
REAL                                 :: n_loc(3),resul(PP_nVar),epsBC
REAL                                 :: U_Face_loc(PP_nVar,0:Nloc,0:Nloc)
!===================================================================================================================================

BCType =BoundaryType(BC(SideID),BC_TYPE)
BCState=BoundaryType(BC(SideID),BC_STATE)

SELECT CASE(BCType)
CASE(1) !Periodic already filled!

CASE(2) ! exact BC = Dirichlet BC !!
  DO q=0,Nloc
    DO p=0,Nloc
      CALL ExactFunc(BCState,t,tDeriv,BCFace_xGP(:,p,q),U_Face_loc(:,p,q))
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(3) ! 1st order absorbing BC
        ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  epsBC=1e-10

  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero.
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,Nloc
    DO p=0,Nloc
      IF (SUM(abs(U_Minus(4:6,p,q))).GT.epsBC)THEN
        U_Face_loc(7,p,q) = - U_Minus(7,p,q) - c*(DOT_PRODUCT(U_Minus(4:6,p,q),normVec(1:3,p,q)))
        U_Face_loc(8,p,q) = - U_Minus(8,p,q) - c_inv*(DOT_PRODUCT(U_Minus(1:3,p,q),normVec(1:3,p,q)))
      END IF ! sum(abs(B)) > epsBC
    END DO ! p
  END DO ! q

  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(4) ! perfectly conducting surface (MunzOmnesSchneider 2000, pp. 97-98)
  ! Determine the exact BC state
  DO q=0,Nloc
    DO p=0,Nloc
      resul=U_Minus(:,p,q)
      n_loc=normVec(:,p,q)
      U_Face_loc(1:3,p,q) = -resul(1:3) + 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
      U_Face_loc(4:6,p,q) =  resul(4:6) - 2*(DOT_PRODUCT(resul(4:6),n_loc))*n_loc
      U_Face_loc(  7,p,q) =  resul(  7)
      U_Face_loc(  8,p,q) = -resul(  8)
    END DO ! p
  END DO ! q
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(5) ! 1st order absorbing BC
      ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  ! A problem of the absorbing BC arises if E or B is close to zero.
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero.
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,Nloc
    DO p=0,Nloc
      U_Face_loc(7,p,q) = - U_Minus(7,p,q) - c*(DOT_PRODUCT(U_Minus(4:6,p,q),normVec(1:3,p,q)))
      U_Face_loc(8,p,q) = - U_Minus(8,p,q) - c_inv*(DOT_PRODUCT(U_Minus(1:3,p,q),normVec(1:3,p,q)))
    END DO ! p
  END DO ! q
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(6) ! 1st order absorbing BC + fix for low B field
        ! Silver-Mueller BC - Munz et al. 2000 / Computer Physics Communication 130, 83-117
  U_Face_loc=0.
  ! A problem of the absorbing BC arises if E or B is close to zero.
  ! Example: electro(dynamic or static) dominated problem and B is approximately zero, than the Silver-Mueller BC requires
  !          that E cross n is zero which is enforced through the div. cleaning of E
  DO q=0,Nloc
    DO p=0,Nloc
      IF (DOT_PRODUCT(U_Minus(4:6,p,q),U_Minus(4:6,p,q))*c*10.&
      .GT.DOT_PRODUCT(U_Minus(1:3,p,q),U_Minus(1:3,p,q)))THEN
        U_Face_loc(7,p,q) = - U_Minus(7,p,q) - c*(DOT_PRODUCT(U_Minus(4:6,p,q),normVec(1:3,p,q)))
        U_Face_loc(8,p,q) = - U_Minus(8,p,q) - c_inv*(DOT_PRODUCT(U_Minus(1:3,p,q),normVec(1:3,p,q)))
      END IF ! sum(abs(B)) > epsBC
    END DO ! p
  END DO ! q
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(10) ! symmetry BC (perfect MAGNETIC conductor, PMC)
  ! Determine the exact BC state
  DO q=0,Nloc
    DO p=0,Nloc
      resul=U_Minus(:,p,q)
      n_loc=normVec(:,p,q)
      U_Face_loc(1:3,p,q) =  resul(1:3) - 2*(DOT_PRODUCT(resul(1:3),n_loc))*n_loc
      U_Face_loc(4:6,p,q) = -resul(4:6) + 2*(DOT_PRODUCT(resul(4:6),n_loc))*n_loc
      U_Face_loc(  7,p,q) = -resul(  7)
      U_Face_loc(  8,p,q) =  resul(  8)
    END DO ! p
  END DO ! q
  ! Dirichlet means that we use the gradients from inside the grid cell
  !CALL RiemannVacuum(Nloc, Flux(:,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),normal(:,:,:))
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
    CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE(20) ! exact BC = Dirichlet BC !!
  ! SPECIAL BC: BCState uses readin state
  CALL abort(__STAMP__,'no state file BC implemented')
  !U_Face_loc=BCData(:,:,:)
  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC)
   CALL RiemannDielectric(Nloc, Flux(1:8,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),&
                           NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
    IF(DoPML) Flux(9:32,:,:) = 0.
  CASE(RIEMANN_PML)
    CALL RiemannPML(Nloc, Flux(1:PP_nVar+PMLnVar,:,:),U_Minus(:,:,:),U_Face_loc(:,:,:),NormVec(:,:,:))
  CASE DEFAULT
    CALL RiemannVacuum(Nloc, Flux(1:PP_nVar,:,:),U_Minus(:,:,:),U_Face_loc(  :,:,:), NormVec(:,:,:))
    IF(DoPML) Flux(9:32,:,:) = 0.
  END SELECT

CASE DEFAULT ! unknown BCType
  CALL abort(__STAMP__,'no BC defined in maxwell/getboundaryflux.f90!')
END SELECT ! BCType

IF(1.EQ.2)THEN
  epsBC=TangVec1(1,1,1)
  epsBC=TangVec2(1,1,1)
END IF

END SUBROUTINE GetBoundaryFlux


SUBROUTINE FinalizeBC()
!===================================================================================================================================
! Initialize boundary conditions
!===================================================================================================================================
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC


#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700)) && !(USE_HDG)*/
END MODULE MOD_GetBoundaryFlux