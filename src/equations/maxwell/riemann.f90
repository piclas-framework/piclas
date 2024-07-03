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

MODULE MOD_Riemann
!===================================================================================================================================
! Contains routines to compute the riemann (Advection, Diffusion) for a given Face
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

INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE RiemannVacuum
  MODULE PROCEDURE RiemannVacuum
END INTERFACE

INTERFACE RiemannPML                ! is called in src/dg/fillflux.f90 (additional 24 auxiliary variables)
  MODULE PROCEDURE RiemannPML
END INTERFACE

INTERFACE RiemannDielectric         ! is called in src/dg/fillflux.f90 (inner dielectric media RP solver: dielectric <-> dielectric)
  MODULE PROCEDURE RiemannDielectric
END INTERFACE

INTERFACE RiemannDielectricInterFace! is called in src/dg/fillflux.f90 (inter dielectric media RP solver: physical <-> dielectric)
  MODULE PROCEDURE RiemannDielectricInterFace
END INTERFACE

INTERFACE RiemannDielectricInterFace2! is called in src/dg/fillflux.f90 (inter dielectric media RP solver: physical <-> dielectric)
  MODULE PROCEDURE RiemannDielectricInterFace2
END INTERFACE

INTERFACE ExactFlux
  MODULE PROCEDURE ExactFlux
END INTERFACE

PUBLIC::Riemann,RiemannVacuum,RiemannPML,Exactflux,RiemannDielectric,RiemannDielectricInterFace,RiemannDielectricInterFace2
!===================================================================================================================================

CONTAINS

SUBROUTINE Riemann(Nloc,Flux_Master,Flux_Slave,U_Master,U_Slave,NormVec,SideID)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Dielectric_vars    ,ONLY: DielectricSurf
USE MOD_Globals            ,ONLY: abort,UNIT_StdOut
USE MOD_PML_vars           ,ONLY: PMLnVar
USE MOD_Interfaces_Vars    ,ONLY: InterfaceRiemann
USE MOD_DG_Vars            ,ONLY: DG_Elems_master
USE MOD_Interpolation_Vars ,ONLY: PREF_VDM
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
#if USE_MPI
USE MOD_Globals            ,ONLY: myrank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: Nloc
INTEGER,INTENT(IN)   :: SideID
REAL,INTENT(IN)      :: NormVec(PP_nVar,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: U_master(PP_nVar,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: U_slave (PP_nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Flux_Master(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
REAL,INTENT(INOUT)   :: Flux_Slave(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE      :: DieLoc(:,:,:), DieTmp(:,:,:)
INTEGER               :: NDie
!===================================================================================================================================
NDie = DG_Elems_master(SideID)
! Check every face and set the correct identifier for selecting the corresponding Riemann solver
! possible connections are (Master <-> Slave direction is important):
SELECT CASE(InterfaceRiemann(SideID))
  !   - vacuum     <-> vacuum       : RIEMANN_VACUUM            = 0
  !   - PML        <-> vacuum       : RIEMANN_PML               = 1
  !   - PML        <-> PML          : RIEMANN_PML               = 1
  !   - dielectric <-> dielectric   : RIEMANN_DIELECTRIC        = 2
  !   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC    = 3 ! for conservative fluxes (one flux)
  !   - vacuum      -> dielectric   : RIEMANN_VAC2DIELECTRIC    = 4 ! for conservative fluxes (one flux)
  !   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC_NC = 5 ! for non-conservative fluxes (two fluxes)
  !   - vacuum      -> dielectric   : RIEMANN_VAC2DIELECTRIC_NC = 6 ! for non-conservative fluxes (two fluxes)
CASE(RIEMANN_VACUUM)
  ! standard flux
  CALL RiemannVacuum(Nloc,Flux_Master(1:8,:,:),U_Master( :,:,:),U_Slave(  :,:,:),NormVec(:,:,:))
CASE(RIEMANN_PML)
  ! RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations (flux-splitting!)
  CALL RiemannPML(Nloc,Flux_Master(1:32,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:))
CASE(RIEMANN_DIELECTRIC)
  ! dielectric region <-> dielectric region
  ALLOCATE(DieLoc(1,0:Nloc,0:Nloc), DieTmp(1,0:NDie,0:NDie))
  IF (NDie.LT.Nloc) THEN
    DieTmp(1,0:NDie,0:NDie) = DielectricSurf(SideID)%Dielectric_Master(0:NDie,0:NDie)
    CALL ChangeBasis2D(1, NDie, Nloc, PREF_VDM(NDie,Nloc)%Vdm, DieTmp(1:1,0:NDie,0:NDie), DieLoc(1:1,0:Nloc,0:Nloc))
  ELSE
    DieLoc(1,0:Nloc,0:Nloc) = DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc)
  END IF
  CALL RiemannDielectric(Nloc,Flux_Master(1:8,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:),DieLoc(1,0:Nloc,0:Nloc))
  SDEALLOCATE(DieLoc)
  SDEALLOCATE(DieTmp)
CASE(RIEMANN_DIELECTRIC2VAC)
  ALLOCATE(DieLoc(1,0:Nloc,0:Nloc), DieTmp(1,0:NDie,0:NDie))
  IF (NDie.LT.Nloc) THEN
    DieTmp(1,0:NDie,0:NDie) = DielectricSurf(SideID)%Dielectric_Master(0:NDie,0:NDie)
    CALL ChangeBasis2D(1, NDie, Nloc, PREF_VDM(NDie,Nloc)%Vdm, DieTmp(1:1,0:NDie,0:NDie), DieLoc(1:1,0:Nloc,0:Nloc))
  ELSE
    DieLoc(1,0:Nloc,0:Nloc) = DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc)
  END IF
  ! master is DIELECTRIC and slave PHYSICAL: A+(Eps0,Mu0) and A-(EpsR,MuR)
  CALL RiemannDielectricInterFace2(Nloc,Flux_Master(1:8,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:),DieLoc(1,0:Nloc,0:Nloc))
  SDEALLOCATE(DieLoc)
  SDEALLOCATE(DieTmp)
CASE(RIEMANN_VAC2DIELECTRIC)
  ALLOCATE(DieLoc(1,0:Nloc,0:Nloc), DieTmp(1,0:NDie,0:NDie))
  IF (NDie.LT.Nloc) THEN
    DieTmp(1,0:NDie,0:NDie) = DielectricSurf(SideID)%Dielectric_Master(0:NDie,0:NDie)
    CALL ChangeBasis2D(1, NDie, Nloc, PREF_VDM(NDie,Nloc)%Vdm, DieTmp(1:1,0:NDie,0:NDie), DieLoc(1:1,0:Nloc,0:Nloc))
  ELSE
    DieLoc(1,0:Nloc,0:Nloc) = DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc)
  END IF
  ! master is PHYSICAL and slave DIELECTRIC: A+(EpsR,MuR) and A-(Eps0,Mu0)
  CALL RiemannDielectricInterFace(Nloc,Flux_Master(1:8,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:),DieLoc(1,0:Nloc,0:Nloc))
  SDEALLOCATE(DieLoc)
  SDEALLOCATE(DieTmp)
CASE(RIEMANN_DIELECTRIC2VAC_NC)  ! use non-conserving fluxes (two different fluxes for master and slave side)
  ALLOCATE(DieLoc(1,0:Nloc,0:Nloc), DieTmp(1,0:NDie,0:NDie))
  IF (NDie.LT.Nloc) THEN
    DieTmp(1,0:NDie,0:NDie) = DielectricSurf(SideID)%Dielectric_Master(0:NDie,0:NDie)
    CALL ChangeBasis2D(1, NDie, Nloc, PREF_VDM(NDie,Nloc)%Vdm, DieTmp(1:1,0:NDie,0:NDie), DieLoc(1:1,0:Nloc,0:Nloc))
  ELSE
    DieLoc(1,0:Nloc,0:Nloc) = DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc)
  END IF
  ! 1.) dielectric master side
  CALL RiemannDielectric(Nloc,Flux_Master(1:8,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:),DieLoc(1,0:Nloc,0:Nloc))
  ! 2.) vacuum slave side
  CALL RiemannVacuum(Nloc,Flux_Slave(1:8,:,:),U_Master( :,:,:),U_Slave(  :,:,:),NormVec(:,:,:))
  SDEALLOCATE(DieLoc)
  SDEALLOCATE(DieTmp)
CASE(RIEMANN_VAC2DIELECTRIC_NC) ! use non-conserving fluxes (two different fluxes for master and slave side)
  ALLOCATE(DieLoc(1,0:Nloc,0:Nloc), DieTmp(1,0:NDie,0:NDie))
  IF (NDie.LT.Nloc) THEN
    DieTmp(1,0:NDie,0:NDie) = DielectricSurf(SideID)%Dielectric_Master(0:NDie,0:NDie)
    CALL ChangeBasis2D(1, NDie, Nloc, PREF_VDM(NDie,Nloc)%Vdm, DieTmp(1:1,0:NDie,0:NDie), DieLoc(1:1,0:Nloc,0:Nloc))
  ELSE
    DieLoc(1,0:Nloc,0:Nloc) = DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc)
  END IF
  ! 1.) dielectric slave side
  CALL RiemannDielectric(Nloc,Flux_Slave(1:8,:,:),U_Master(:,:,:),U_Slave(:,:,:),NormVec(:,:,:),DieLoc(1,0:Nloc,0:Nloc))
  ! 2.) vacuum master side
  CALL RiemannVacuum(Nloc,Flux_Master(1:8,:,:),U_Master( :,:,:),U_Slave(  :,:,:),NormVec(:,:,:))
  SDEALLOCATE(DieLoc)
  SDEALLOCATE(DieTmp)
CASE DEFAULT
  IPWRITE(UNIT_StdOut,*) "SideID                   = ", SideID
  IPWRITE(UNIT_StdOut,*) "InterfaceRiemann(SideID) = ", InterfaceRiemann(SideID)
  CALL abort(__STAMP__,'Unknown interface type for Riemann solver (vacuum, dielectric, PML ...)')
END SELECT

END SUBROUTINE Riemann


SUBROUTINE RiemannVacuum(Nloc,F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: eta_c,c_corr,c_corr_c,c_corr_c2,CentralFlux
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),A_p(8,8),A_n(8,8), A(8,8)
INTEGER                                          :: Count_1,Count_2
!REAL                                             :: D(3,3)                  ! auxiliary matrices used
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!===================================================================================================================================

IF(CentralFlux)THEN
  ! Gauss point i,j
  DO Count_2=0,Nloc
    DO Count_1=0,Nloc
      n_loc(:)=nv(:,Count_1,Count_2)
      A=0.
      ! check again 7 and 8
      A(5,1)=n_loc(3)
      A(6,1)=-n_loc(2)
      A(8,1)=c_corr*n_loc(1)
      A(4,2)=-n_loc(3)
      A(6,2)=n_loc(1)
      A(8,2)=c_corr*n_loc(2)
      A(4,3)=n_loc(2)
      A(5,3)=-n_loc(1)
      A(8,3)=c_corr*n_loc(3)

      A(2,4)=       -c2*n_loc(3)
      A(3,4)=        c2*n_loc(2)
      A(7,4)= c_corr_c2*n_loc(1)
      A(1,5)=        c2*n_loc(3)
      A(3,5)=       -c2*n_loc(1)
      A(7,5)= c_corr_c2*n_loc(2)
      A(1,6)=       -c2*n_loc(2)
      A(2,6)=        c2*n_loc(1)
      A(7,6)= c_corr_c2*n_loc(3)

      A(4,7)= c_corr*n_loc(1)
      A(5,7)= c_corr*n_loc(2)
      A(6,7)= c_corr*n_loc(3)

      A(1,8)= c_corr_c2*n_loc(1)
      A(2,8)= c_corr_c2*n_loc(2)
      A(3,8)= c_corr_c2*n_loc(3)

      F(:,Count_1,Count_2)=0.5*(MATMUL(A,U_R(:,Count_1,Count_2))+MATMUL(A,U_L(:,Count_1,Count_2)))
    END DO
  END DO
ELSE
  ! Gauss point i,j
  DO Count_2=0,Nloc
    DO Count_1=0,Nloc
      n_loc(:)=nv(:,Count_1,Count_2)

  !--- for original version see below (easier to understand)

      A_p(7,1:3)=0.
      A_p(1:3,7)=0.
      A_p(8,4:7)=0.
      A_p(4:7,8)=0.

      !D-Teilmatrix: Since chi and gamma is equal we
      ! consider D(chi,gamma) = D(gamma,chi)
      ! ATTENTION: if chi .ne. gamma this have to be changed.
      ! Then we need D_1 and D_2 (see commented section below)
      A_p(1,1) = c + n_loc(1)*n_loc(1)*eta_c   !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
      A_p(1,2) = n_loc(1)*n_loc(2)*eta_c       !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(1,3) = n_loc(1)*n_loc(3)*eta_c       !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(2,1) = A_p(1,2)                      !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(2,2) = c + n_loc(2)*n_loc(2)*eta_c   !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
      A_p(2,3) = n_loc(2)*n_loc(3)*eta_c       !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_p(3,1) = A_p(1,3)                      !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(3,2) = A_p(2,3)                      !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_p(3,3) = c+n_loc(3)*n_loc(3)*eta_c     !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c
      ! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_p(1,4:6)= (/0.,c2*n_loc(3),-c2*n_loc(2)/)
      A_p(2,4:6)= (/-c2*n_loc(3),0.,c2*n_loc(1)/)
      A_p(3,4:6)= (/c2*n_loc(2),-c2*n_loc(1),0./)
      A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      !composition of the Matrix
      !positive A-Matrx
      A_p(4:6,4:6)=A_p(1:3,1:3)
      !negative A-Matrix
      A_n(1:3,1:3)=-A_p(1:3,1:3)
      A_n(1:3,4:6)= A_p(1:3,4:6)   ! c*c*E(:,:)
      A_n(4:6,1:3)= A_p(4:6,1:3)
      A_n(4:6,4:6)=-A_p(1:3,1:3)

     ! !positive A-Matrix-Divergence-Correction-Term
      A_p(1,8) = c_corr_c2*n_loc(1)
      A_p(2,8) = c_corr_c2*n_loc(2)
      A_p(3,8) = c_corr_c2*n_loc(3)
      A_p(4,7) = c_corr*n_loc(1)
      A_p(5,7) = c_corr*n_loc(2)
      A_p(6,7) = c_corr*n_loc(3)
      A_p(7,4) = c_corr_c2*n_loc(1)
      A_p(7,5) = c_corr_c2*n_loc(2)
      A_p(7,6) = c_corr_c2*n_loc(3)
      A_p(7,7) = c_corr_c
      A_p(8,1) = c_corr*n_loc(1)
      A_p(8,2) = c_corr*n_loc(2)
      A_p(8,3) = c_corr*n_loc(3)
      A_p(8,8) = c_corr_c
      !negative A-Matrix-Divergence-Correction-Term
      A_n(1:7,8) = A_p(1:7,8) !c_corr*c*c*n(1)
      A_n(1:6,7) = A_p(1:6,7) !c_corr*n(1)
      A_n(7,1:6)=  A_p(7,1:6)
      A_n(7,7)  = -A_p(7,7)
      A_n(8,1:7)= A_p(8,1:7)
      A_n(8,8)  =-A_p(8,8)
      ! Warum 0.5 -> Antwort im Taube/Dumbser-Paper. Im Munz/Schneider Paper fehlt das 0.5 lustigerweise.

  !--- Original Version:

  !    A_p=0.
  !    A_n=0.
  !
  !    !D-Teilmatrix: Since chi and gamma is equal we
  !    ! consider D(chi,gamma) = D(gamma,chi)
  !    ! ATTENTION: if chi .ne. gamma this have to be changed.
  !    ! Then we need D_1 and D_2
  !    D(1,1) = c + n_loc(1)*n_loc(1)*eta_c   !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
  !    D(1,2) = n_loc(1)*n_loc(2)*eta_c            !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
  !    D(1,3) = n_loc(1)*n_loc(3)*eta_c            !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
  !    D(2,1) = D(1,2)                          !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
  !    D(2,2) = c + n_loc(2)*n_loc(2)*eta_c   !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
  !    D(2,3) = n_loc(2)*n_loc(3)*eta_c            !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
  !    D(3,1) = D(1,3)                          !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
  !    D(3,2) = D(2,3)                          !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c
  !    D(3,3) = c+n_loc(3)*n_loc(3)*eta_c     !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c
  !    ! epsilon-Teilmatrix
  !    !E_trans=transpose(E)
  !    E(1,:)= (/0.,n_loc(3),-n_loc(2)/)
  !    E(2,:)= (/-n_loc(3),0.,n_loc(1)/)
  !    E(3,:)= (/n_loc(2),-n_loc(1),0./)
  !    E_trans(1,:)= (/0.,-n_loc(3),n_loc(2)/)
  !    E_trans(2,:)= (/n_loc(3),0.,-n_loc(1)/)
  !    E_trans(3,:)= (/-n_loc(2),n_loc(1),0./)
  !    !composition of the Matrix
  !    !positive A-Matrx
  !    A_p(1:3,1:3)=D(:,:)
  !    A_p(1:3,4:6)=c2*E(:,:) ! c*c*E(:,:)
  !    A_p(4:6,1:3)=E_trans(:,:)
  !    A_p(4:6,4:6)=D(:,:)
  !    !negative A-Matrix
  !    A_n(1:3,1:3)=-D(:,:)
  !    A_n(1:3,4:6)= A_p(1:3,4:6)   ! c*c*E(:,:)
  !    A_n(4:6,1:3)= E_trans(:,:)
  !    A_n(4:6,4:6)=-D(:,:)
  !
  !   ! !positive A-Matrix-Divergence-Correction-Term
  !    A_p(1,8) = c_corr_c2*n_loc(1)
  !    A_p(2,8) = c_corr_c2*n_loc(2)
  !    A_p(3,8) = c_corr_c2*n_loc(3)
  !    A_p(4,7) = c_corr*n_loc(1)
  !    A_p(5,7) = c_corr*n_loc(2)
  !    A_p(6,7) = c_corr*n_loc(3)
  !    A_p(7,4) = c_corr_c2*n_loc(1)
  !    A_p(7,5) = c_corr_c2*n_loc(2)
  !    A_p(7,6) = c_corr_c2*n_loc(3)
  !    A_p(7,7) = c_corr_c
  !    A_p(8,1) = c_corr*n_loc(1)
  !    A_p(8,2) = c_corr*n_loc(2)
  !    A_p(8,3) = c_corr*n_loc(3)
  !    A_p(8,8) = c_corr_c
  !    !negative A-Matrix-Divergence-Correction-Term
  !    A_n(1,8) = A_p(1,8) !c_corr*c*c*n(1)
  !    A_n(2,8) = A_p(2,8) !c_corr*c*c*n(2)
  !    A_n(3,8) = A_p(3,8) !c_corr*c*c*n(3)
  !    A_n(4,7) = A_p(4,7) !c_corr*n(1)
  !    A_n(5,7) = A_p(5,7) !c_corr*n(2)
  !    A_n(6,7) = A_p(6,7) !c_corr*n(3)
  !    !A_n(7,:)=(/0.,0.,0.,c_corr*c*c*n(1),c_corr*c*c*n(2),c_corr*c*c*n(3),-c_corr*c,0./)
  !    A_n(7,1:6)=  A_p(7,1:6)
  !    A_n(7,7)  = -A_p(7,7)
  !    !A_n(8,:)=(/c_corr*n(1),c_corr*n(2),c_corr*n(3),0.,0.,0.,0.,-c_corr*c/)
  !    A_n(8,1:7)= A_p(8,1:7)
  !    A_n(8,8)  =-A_p(8,8)
  !    ! Warum 0.5 -> Antwort im Taube/Dumbser-Paper. Im Munz/Schneider Paper fehlt das 0.5 lustigerweise.

      F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
    END DO
  END DO
END IF

END SUBROUTINE RiemannVacuum


SUBROUTINE RiemannPML(Nloc,F,U_L,U_R,nv)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: c_corr,c_corr_c,c_corr_c2
USE MOD_PML_vars      ,ONLY: PMLnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc) !,t1(3,0:Nloc,0:Nloc),t2(3,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar+PMLnVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3)
INTEGER                                          :: Count_1,Count_2
!REAL                                             :: D(3,3)                  ! auxiliary matrices used
!REAL                                             :: E(3,3), E_trans(3,3)    ! auxiliary matrices used
!===================================================================================================================================
! Gauss point i,j

DO Count_2=0,Nloc
  DO Count_1=0,Nloc
    n_loc(:)=nv(:,Count_1,Count_2)

    ! A^-*U_R + A^+*U_L
    !===============================================================================================================================
    !P1, P2, P3 and Ex
    F(9,Count_1,Count_2)  = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c_corr_c*n_loc(1)*n_loc(1),-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(1)/) )+&
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c_corr_c*n_loc(1)*n_loc(1), c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(1)/) ))

    F(10,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c*n_loc(2)*n_loc(2), c*n_loc(1)*n_loc(2),0.,0.,0.,-c2*n_loc(2), c*n_loc(2)*n_loc(3),0./) )+ &
                  DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c*n_loc(2)*n_loc(2),-c*n_loc(1)*n_loc(2),0.,0.,0.,-c2*n_loc(2),-c*n_loc(2)*n_loc(3),0./)))

    F(11,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/-c*n_loc(3)*n_loc(3),0., c*n_loc(1)*n_loc(3),0.,c2*n_loc(3),0.,-c*n_loc(2)*n_loc(3),0./) )+ &
                  DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/ c*n_loc(3)*n_loc(3),0.,-c*n_loc(1)*n_loc(3),0.,c2*n_loc(3),0., c*n_loc(2)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(1,Count_1,Count_2)  = F(9,Count_1,Count_2) + F(10,Count_1,Count_2) + F(11,Count_1,Count_2)
    !!==============================================================================================================================
    !Q1, Q2, Q3 and Ey
    F(12,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/ c*n_loc(1)*n_loc(2),-c*n_loc(1)*n_loc(1),0.,0.,0.,c2*n_loc(1),-c*n_loc(1)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/-c*n_loc(1)*n_loc(2), c*n_loc(1)*n_loc(1),0.,0.,0.,c2*n_loc(1), c*n_loc(1)*n_loc(3),0./) ))

    F(13,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
(/-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
(/ c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(2)*n_loc(2), c_corr_c*n_loc(2)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(2)/)))

    F(14,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
          (/0.,-c*n_loc(3)*n_loc(3), c*n_loc(2)*n_loc(3),-c2*n_loc(3),0.,0., c*n_loc(1)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
          (/0., c*n_loc(3)*n_loc(3),-c*n_loc(2)*n_loc(3),-c2*n_loc(3),0.,0.,-c*n_loc(1)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(2,Count_1,Count_2)  = F(12,Count_1,Count_2) + F(13,Count_1,Count_2) + F(14,Count_1,Count_2)
    !!==============================================================================================================================
    !R1, R2, R3 and Ez
    F(15,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/ c*n_loc(1)*n_loc(3),0.,-c*n_loc(1)*n_loc(1),0.,-c2*n_loc(1),0., c*n_loc(1)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/-c*n_loc(1)*n_loc(3),0., c*n_loc(1)*n_loc(1),0.,-c2*n_loc(1),0.,-c*n_loc(1)*n_loc(2),0./) ))

    F(16,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0., c*n_loc(2)*n_loc(3),-c*n_loc(2)*n_loc(2),c2*n_loc(2),0.,0.,-c*n_loc(1)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,-c*n_loc(2)*n_loc(3), c*n_loc(2)*n_loc(2),c2*n_loc(2),0.,0., c*n_loc(1)*n_loc(2),0./)))

    F(17,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
(/-c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(3)/))+&
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
(/ c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),0.,0.,0.,0.,c_corr_c2*n_loc(3)/)))
    !!==============================================================================================================================
    F(3,Count_1,Count_2)  = F(15,Count_1,Count_2) + F(16,Count_1,Count_2) + F(17,Count_1,Count_2)
    !!==============================================================================================================================
    !L1, L2, L3 and Bx
    F(18,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(1),-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(1)*n_loc(3),c_corr*n_loc(1),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(1), c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(1)*n_loc(3),c_corr*n_loc(1),0./) ))

    F(19,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0., n_loc(2),-c*n_loc(2)*n_loc(2), c*n_loc(1)*n_loc(2),0.,0.,-c*n_loc(2)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0., n_loc(2), c*n_loc(2)*n_loc(2),-c*n_loc(1)*n_loc(2),0.,0., c*n_loc(2)*n_loc(3)/)))

    F(20,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,-n_loc(3),0.,-c*n_loc(3)*n_loc(3),0., c*n_loc(1)*n_loc(3),0., c*n_loc(2)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,-n_loc(3),0., c*n_loc(3)*n_loc(3),0.,-c*n_loc(1)*n_loc(3),0.,-c*n_loc(2)*n_loc(3)/)))
    !!==============================================================================================================================
    F(4,Count_1,Count_2)  = F(18,Count_1,Count_2) + F(19,Count_1,Count_2) + F(20,Count_1,Count_2)
    !!==============================================================================================================================
    !!M1, M2, M3 and By
    F(21,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,-n_loc(1), c*n_loc(1)*n_loc(2),-c*n_loc(1)*n_loc(1),0.,0., c*n_loc(1)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,-n_loc(1),-c*n_loc(1)*n_loc(2), c*n_loc(1)*n_loc(1),0.,0.,-c*n_loc(1)*n_loc(3)/) ))

    F(22,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(2),-c_corr_c*n_loc(2)*n_loc(3),c_corr*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(2), c_corr_c*n_loc(2)*n_loc(2), c_corr_c*n_loc(2)*n_loc(3),c_corr*n_loc(2),0./)))

    F(23,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/n_loc(3),0.,0.,0.,-c*n_loc(3)*n_loc(3), c*n_loc(2)*n_loc(3),0.,-c*n_loc(1)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/n_loc(3),0.,0.,0., c*n_loc(3)*n_loc(3),-c*n_loc(2)*n_loc(3),0., c*n_loc(1)*n_loc(3)/)))
    !!==============================================================================================================================
    F(5,Count_1,Count_2)  = F(21,Count_1,Count_2) + F(22,Count_1,Count_2) + F(23,Count_1,Count_2)
    !!==============================================================================================================================
    !N1, N2, N3 and Bz
    F(24,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,n_loc(1),0., c*n_loc(1)*n_loc(3),0.,-c*n_loc(1)*n_loc(1),0.,-c*n_loc(1)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,n_loc(1),0.,-c*n_loc(1)*n_loc(3),0., c*n_loc(1)*n_loc(1),0., c*n_loc(1)*n_loc(2)/) ))

    F(25,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/-n_loc(2),0.,0.,0., c*n_loc(2)*n_loc(3),-c*n_loc(2)*n_loc(2),0., c*n_loc(1)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/-n_loc(2),0.,0.,0.,-c*n_loc(2)*n_loc(3), c*n_loc(2)*n_loc(2),0.,-c*n_loc(1)*n_loc(2)/)))

    F(26,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),c_corr*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,0., c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),c_corr*n_loc(3),0./)))
    !!==============================================================================================================================
    F(6,Count_1,Count_2)  = F(24,Count_1,Count_2) + F(25,Count_1,Count_2) + F(26,Count_1,Count_2)
    !!==============================================================================================================================
    !S1, S2, S3 and PhiB
    F(27,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/0., c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(2),c_corr_c2*n_loc(1),0.,0.,-c_corr_c*n_loc(1)*n_loc(1),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/0.,-c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(1)*n_loc(2),c_corr_c2*n_loc(1),0.,0., c_corr_c*n_loc(1)*n_loc(1),0./) ))

    F(28,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/-c_corr_c*n_loc(2)*n_loc(3),0., c_corr_c*n_loc(1)*n_loc(2),0.,c_corr_c2*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(2),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/ c_corr_c*n_loc(2)*n_loc(3),0.,-c_corr_c*n_loc(1)*n_loc(2),0.,c_corr_c2*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(2),0./)))

    F(29,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
    (/ c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,c_corr_c2*n_loc(3),-c_corr_c*n_loc(3)*n_loc(3),0./) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
    (/-c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(1)*n_loc(3),0.,0.,0.,c_corr_c2*n_loc(3), c_corr_c*n_loc(3)*n_loc(3),0./)))
    !!==============================================================================================================================
    F(7,Count_1,Count_2)  = F(27,Count_1,Count_2) + F(28,Count_1,Count_2) + F(29,Count_1,Count_2)
    !!==============================================================================================================================
    !T1, T2, T3 and PhiE
    F(30,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/c_corr*n_loc(1),0.,0.,0.,-c_corr_c*n_loc(1)*n_loc(3), c_corr_c*n_loc(1)*n_loc(2),0.,-c_corr_c*n_loc(1)*n_loc(1)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/c_corr*n_loc(1),0.,0.,0., c_corr_c*n_loc(1)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(2),0., c_corr_c*n_loc(1)*n_loc(1)/) ))

    F(31,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,c_corr*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(3),0.,-c_corr_c*n_loc(1)*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(2)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,c_corr*n_loc(2),0.,-c_corr_c*n_loc(2)*n_loc(3),0., c_corr_c*n_loc(1)*n_loc(2),0., c_corr_c*n_loc(2)*n_loc(2)/)))

    F(32,Count_1,Count_2) = 0.5 * (DOT_PRODUCT(U_R(:,Count_1,Count_2), &
           (/0.,0.,c_corr*n_loc(3),-c_corr_c*n_loc(2)*n_loc(3), c_corr_c*n_loc(1)*n_loc(3),0.,0.,-c_corr_c*n_loc(3)*n_loc(3)/) )+ &
                                   DOT_PRODUCT(U_L(:,Count_1,Count_2), &
           (/0.,0.,c_corr*n_loc(3), c_corr_c*n_loc(2)*n_loc(3),-c_corr_c*n_loc(1)*n_loc(3),0.,0., c_corr_c*n_loc(3)*n_loc(3)/)))
    !!==============================================================================================================================
    !F(8) passt, durch ausgabe getestet
    F(8,Count_1,Count_2)  = F(30,Count_1,Count_2) + F(31,Count_1,Count_2) + F(32,Count_1,Count_2)
    !!==============================================================================================================================
  END DO
END DO

END SUBROUTINE RiemannPML


SUBROUTINE ExactFlux(SideID,t,tDeriv,Nloc,Flux_Master,Flux_Slave,U_Master, U_slave,NormVec,Face_xGP,SurfElem)
!===================================================================================================================================
! Routine to add an exact function to a Riemann-Problem Face
! used at PML interfaces to emit a wave
! The ExactFlux is a non-conservative flux and is only emitted in ONE direction
! mapping
!               |
!       Master  |  Slave
!               |
!        nvec ----->
!               |
!               |
! Flux_Master   | Flux_Slave
! CAUTION: This routine has to multiply the new flux with the SurfElem
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
!USE MOD_DG_Vars,         ONLY:U_Master_loc,U_Slave_loc,Flux_loc
USE MOD_Equation_Vars,   ONLY:IniExactFunc,ExactFluxDir
USE MOD_Equation,        ONLY:ExactFunc
USE MOD_PML_vars,        ONLY:PMLnVar
USE MOD_Interfaces_Vars, ONLY:InterfaceRiemann
USE MOD_Dielectric_vars, ONLY:DielectricSurf
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: SideID
INTEGER,INTENT(IN)    :: Nloc
REAL,INTENT(IN)       :: t
INTEGER,INTENT(IN)    :: tDeriv
REAL,INTENT(IN)       :: NormVec( 1:3      ,0:Nloc,0:Nloc)
REAL,INTENT(IN)       :: Face_xGP(1:3      ,0:Nloc,0:Nloc)
REAL,INTENT(IN)       :: U_master(1:PP_nVar,0:Nloc,0:Nloc)
REAL,INTENT(IN)       :: U_slave (1:PP_nVar,0:Nloc,0:Nloc)
REAL,INTENT(IN)       :: SurfElem(          0:Nloc,0:Nloc)
! Allocate arrays to hold the face flux to reduce memory churn
REAL                  :: U_Master_loc(1:PP_nVar        ,0:Nloc,0:Nloc)
REAL                  :: U_Slave_loc (1:PP_nVar        ,0:Nloc,0:Nloc)
REAL                  :: Flux_loc(    1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: Flux_Master(1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
REAL,INTENT(INOUT)    :: Flux_Slave (1:PP_nVar+PMLnVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
LOGICAL               :: UseMaster
REAL                  :: U_loc(1:PP_nVar)
!===================================================================================================================================

UseMaster=.TRUE.
! emission over plane, hence, first entry decides orientation of  plane
IF(NormVec(ExactFluxDir,0,0).GT.0)THEN
  UseMaster=.FALSE.
ELSE IF(NormVec(ExactFluxDir,0,0).LT.0)THEN
  UseMaster=.TRUE.
ELSE
  CALL abort(__STAMP__,'ExactFlux has encountered an error with the normal vector.')
END IF

U_Slave_loc  = U_Slave
U_Master_loc = U_Master

DO q=0,Nloc
  DO p=0,Nloc
    ! the second state is always zero and already computed
    ! caution, rotation
    CALL ExactFunc(IniExactFunc,t,tDeriv,Face_xGP(:,p,q),U_loc(:))
    IF(UseMaster)THEN
      U_Slave_loc(:,p,q)  = U_Slave_loc(:,p,q)  + U_loc
    ELSE
      U_Master_loc(:,p,q) = U_Master_loc(:,p,q) + U_loc
    END IF
  END DO ! p
END DO ! q

Flux_loc=0.

! check interface type for Riemann solver selection
SELECT CASE(InterfaceRiemann(SideID))
CASE(RIEMANN_VACUUM) ! standard flux
  CALL RiemannVacuum(Nloc,Flux_loc(1:8,:,:),U_Master_loc( :,:,:),U_Slave_loc(  :,:,:),NormVec(:,:,:))
CASE(RIEMANN_PML) ! RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations (flux-splitting!)
  CALL RiemannPML(Nloc,Flux_loc(1:32,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:), NormVec(:,:,:))
CASE(RIEMANN_DIELECTRIC) ! dielectric region <-> dielectric region
  CALL RiemannDielectric(Nloc,Flux_loc(1:8,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:),&
                         NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
CASE(RIEMANN_DIELECTRIC2VAC) ! master is DIELECTRIC and slave PHYSICAL: A+(Eps0,Mu0) and A-(EpsR,MuR)
  CALL RiemannDielectricInterFace2(Nloc,Flux_loc(1:8,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:),&
                                              NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
CASE(RIEMANN_VAC2DIELECTRIC) ! master is PHYSICAL and slave DIELECTRIC: A+(EpsR,MuR) and A-(Eps0,Mu0)
  CALL RiemannDielectricInterFace(Nloc,Flux_loc(1:8,:,:),U_Master_loc(:,:,:),U_Slave_loc(:,:,:),&
                                              NormVec(:,:,:),DielectricSurf(SideID)%Dielectric_Master(0:Nloc,0:Nloc))
CASE DEFAULT
  CALL abort(__STAMP__,'Unknown interface type for Riemann solver (vacuum, dielectric, PML ...)')
END SELECT

IF(UseMaster)THEN
  DO q=0,Nloc
    DO p=0,Nloc
      Flux_Master(:,p,q) = Flux_loc(:,p,q)*SurfElem(p,q)
    END DO ! p
  END DO ! q
ELSE
  DO q=0,Nloc
    DO p=0,Nloc
      Flux_Slave(:,p,q)  = Flux_loc(:,p,q)*SurfElem(p,q)
    END DO ! p
  END DO ! q
END IF

END SUBROUTINE ExactFlux


SUBROUTINE RiemannDielectric(Nloc,F,U_L,U_R,nv,Dielectric_Master)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c
USE MOD_Equation_Vars ,ONLY: CentralFlux,c_corr,c_corr_c
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)
REAL,INTENT(IN)                                  :: Dielectric_Master(0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),A_p(8,8),A_n(8,8)
INTEGER                                          :: Count_1,Count_2
REAL                                             :: eta_c_dielectric,c_dielectric,c2_dielectric
REAL                                             :: c_corr_c_dielectric,c_corr_c2_dielectric
!===================================================================================================================================

IF(CentralFlux)THEN
  CALL abort(__STAMP__,'central flux for dielectric media not implemented!')
ELSE
  ! Gauss point i,j
  DO Count_2=0,Nloc
    DO Count_1=0,Nloc
      n_loc(:)=nv(:,Count_1,Count_2)

      ! set dielectric values
      c_dielectric         = c*Dielectric_Master(Count_1,Count_2)                       !            c/sqrt(EpsR*MuR)
      c_corr_c_dielectric  = c_corr_c*Dielectric_Master(Count_1,Count_2)                !        chi*c/sqrt(EpsR*MuR)
      eta_c_dielectric     = c_corr_c_dielectric - c_dielectric                         ! (chi - 1.)*c/sqrt(EpsR*MuR)
      c2_dielectric        = c_dielectric*c_dielectric                                  !             c**2/(EpsR*MuR)
      c_corr_c2_dielectric = c_corr * c2_dielectric                                     !         chi*c**2/(EpsR*MuR)

      A_p = -1234567.
      A_n = -1234567.

      A_p(7,1:3)=0.
      A_p(1:3,7)=0.
      A_p(8,4:7)=0.
      A_p(4:7,8)=0.

      !D-Teilmatrix: Since chi and gamma is equal we
      ! consider D(chi,gamma) = D(gamma,chi)
      ! ATTENTION: if chi .ne. gamma this have to be changed.
      ! Then we need D_1 and D_2 (see commented section below)
      A_p(1,1) = c_dielectric + n_loc(1)*n_loc(1)*eta_c_dielectric        !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(chi-1.))*c/sqrt(EpsR*MuR)
      A_p(1,2) = n_loc(1)*n_loc(2)*eta_c_dielectric                       !  D(1,2)=n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(1,3) = n_loc(1)*n_loc(3)*eta_c_dielectric                       !  D(1,3)=n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(2,1) = A_p(1,2)                                                 !  D(2,1)=n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(2,2) = c_dielectric + n_loc(2)*n_loc(2)*eta_c_dielectric        !  D(2,2)=(1./sqrt(EpsR*MuR)+n_loc(2)*n_loc(2)*(chi-1.))*c
      A_p(2,3) = n_loc(2)*n_loc(3)*eta_c_dielectric                       !  D(2,3)=n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,1) = A_p(1,3)                                                 !  D(3,1)=n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,2) = A_p(2,3)                                                 !  D(3,2)=n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,3) = c_dielectric+n_loc(3)*n_loc(3)*eta_c_dielectric          !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c/sqrt(EpsR*MuR)
      ! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_p(1,4:6)= (/0.,c2_dielectric*n_loc(3),-c2_dielectric*n_loc(2)/)
      A_p(2,4:6)= (/-c2_dielectric*n_loc(3),0.,c2_dielectric*n_loc(1)/)
      A_p(3,4:6)= (/c2_dielectric*n_loc(2),-c2_dielectric*n_loc(1),0./)
      A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      !composition of the Matrix
      !positive A-Matrx
      A_p(4:6,4:6)=A_p(1:3,1:3)
      !negative A-Matrix
      A_n(1:3,1:3)=-A_p(1:3,1:3)
      A_n(1:3,4:6)= A_p(1:3,4:6)   ! c*c*E(:,:)
      A_n(4:6,1:3)= A_p(4:6,1:3)
      A_n(4:6,4:6)=-A_p(1:3,1:3)

     ! !positive A-Matrix-Divergence-Correction-Term
      A_p(1,8) = c_corr_c2_dielectric*n_loc(1)
      A_p(2,8) = c_corr_c2_dielectric*n_loc(2)
      A_p(3,8) = c_corr_c2_dielectric*n_loc(3)
      A_p(4,7) = c_corr*n_loc(1)
      A_p(5,7) = c_corr*n_loc(2)
      A_p(6,7) = c_corr*n_loc(3)
      A_p(7,4) = c_corr_c2_dielectric*n_loc(1)
      A_p(7,5) = c_corr_c2_dielectric*n_loc(2)
      A_p(7,6) = c_corr_c2_dielectric*n_loc(3)
      A_p(7,7) = c_corr_c_dielectric
      A_p(8,1) = c_corr*n_loc(1)
      A_p(8,2) = c_corr*n_loc(2)
      A_p(8,3) = c_corr*n_loc(3)
      A_p(8,8) = c_corr_c_dielectric
      !negative A-Matrix-Divergence-Correction-Term
      A_n(1:7,8) =  A_p(1:7,8) !c_corr*c*c*n(1)
      A_n(1:6,7) =  A_p(1:6,7) !c_corr*n(1)
      A_n(7,1:6) =  A_p(7,1:6)
      A_n(7,7)   = -A_p(7,7)
      A_n(8,1:7) =  A_p(8,1:7)
      A_n(8,8)   = -A_p(8,8)

      IF(MINVAL(ABS(A_n+1234567))==0)THEN
        CALL abort(__STAMP__,'A_n RiemannDielectric')
      END IF
      IF(MINVAL(ABS(A_p+1234567))==0)THEN
        CALL abort(__STAMP__,'A_p RiemannDielectric')
      END IF
      F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
    END DO
  END DO
END IF

END SUBROUTINE RiemannDielectric


SUBROUTINE RiemannDielectricInterFace(Nloc,F,U_L,U_R,nv,Dielectric_Master)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
! master is PHYSICAL and slave DIELECTRIC: A+(EpsR,MuR) and A-(Eps0,Mu0)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: CentralFlux,eta_c,c_corr,c_corr_c,c_corr_c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)
REAL,INTENT(IN)                                  :: Dielectric_Master(0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),A_p(8,8),A_n(8,8)
INTEGER                                          :: Count_1,Count_2
REAL                                             :: eta_c_dielectric,c_dielectric,c2_dielectric
REAL                                             :: c_corr_c_dielectric,c_corr_c2_dielectric
!===================================================================================================================================

IF(CentralFlux)THEN
  CALL abort(__STAMP__,'central flux for dielectric media not implemented!')
ELSE
  ! Gauss point i,j
  DO Count_2=0,Nloc
    DO Count_1=0,Nloc
      n_loc(:)=nv(:,Count_1,Count_2)

      ! set dielectric values
      c_dielectric         = c*Dielectric_Master(Count_1,Count_2)                       !            c/sqrt(EpsR*MuR)
      c_corr_c_dielectric  = c_corr_c*Dielectric_Master(Count_1,Count_2)                !        chi*c/sqrt(EpsR*MuR)
      eta_c_dielectric     = c_corr_c_dielectric - c_dielectric                         ! (chi - 1.)*c/sqrt(EpsR*MuR)
      c2_dielectric        = c_dielectric*c_dielectric                                  !             c**2/(EpsR*MuR)
      c_corr_c2_dielectric = c_corr * c2_dielectric                                     !         chi*c**2/(EpsR*MuR)

      A_n(7,1:3)=0.
      A_n(1:3,7)=0.
      A_n(8,4:7)=0.
      A_n(4:7,8)=0.

      A_p(7,1:3)=0.
      A_p(1:3,7)=0.
      A_p(8,4:7)=0.
      A_p(4:7,8)=0.

      !D-Teilmatrix: Since chi and gamma is equal we
      ! consider D(chi,gamma) = D(gamma,chi)
      ! ATTENTION: if chi .ne. gamma this have to be changed.
      ! Then we need D_1 and D_2 (see commented section below)
      A_p(1,1) = c_dielectric + n_loc(1)*n_loc(1)*eta_c_dielectric        !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(chi-1.))*c/sqrt(EpsR*MuR)
      A_p(1,2) = n_loc(1)*n_loc(2)*eta_c_dielectric                       !  D(1,2)=n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(1,3) = n_loc(1)*n_loc(3)*eta_c_dielectric                       !  D(1,3)=n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(2,1) = A_p(1,2)                                                 !  D(2,1)=n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(2,2) = c_dielectric + n_loc(2)*n_loc(2)*eta_c_dielectric        !  D(2,2)=(1./sqrt(EpsR*MuR)+n_loc(2)*n_loc(2)*(chi-1.))*c
      A_p(2,3) = n_loc(2)*n_loc(3)*eta_c_dielectric                       !  D(2,3)=n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,1) = A_p(1,3)                                                 !  D(3,1)=n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,2) = A_p(2,3)                                                 !  D(3,2)=n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)
      A_p(3,3) = c_dielectric+n_loc(3)*n_loc(3)*eta_c_dielectric          !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c/sqrt(EpsR*MuR)

! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_p(1,4:6)= (/0.,c2_dielectric*n_loc(3),-c2_dielectric*n_loc(2)/)
      A_p(2,4:6)= (/-c2_dielectric*n_loc(3),0.,c2_dielectric*n_loc(1)/)
      A_p(3,4:6)= (/c2_dielectric*n_loc(2),-c2_dielectric*n_loc(1),0./)
      A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      !composition of the Matrix
      !positive A-Matrx
      A_p(4:6,4:6)=A_p(1:3,1:3)

      !negative A-Matrix
      A_n(1,1) = -( c + n_loc(1)*n_loc(1)*eta_c )      !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(chi-1.))*c
      A_n(1,2) = -( n_loc(1)*n_loc(2)*eta_c     )      !  D(1,2)=n_loc(1)*n_loc(2)*(chi-1)*c
      A_n(1,3) = -( n_loc(1)*n_loc(3)*eta_c     )      !  D(1,3)=n_loc(1)*n_loc(3)*(chi-1)*c
      A_n(2,1) = -( A_n(1,2)                    )      !  D(2,1)=n_loc(1)*n_loc(2)*(chi-1)*c
      A_n(2,2) = -( c + n_loc(2)*n_loc(2)*eta_c )      !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(chi-1.))*c
      A_n(2,3) = -( n_loc(2)*n_loc(3)*eta_c     )      !  D(2,3)=n_loc(2)*n_loc(3)*(chi-1)*c
      A_n(3,1) = -( A_n(1,3)                    )      !  D(3,1)=n_loc(1)*n_loc(3)*(chi-1)*c
      A_n(3,2) = -( A_n(2,3)                    )      !  D(3,2)=n_loc(2)*n_loc(3)*(chi-1)*c
      A_n(3,3) = -( c+n_loc(3)*n_loc(3)*eta_c   )      !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c

      A_n(4:6,4:6)=A_n(1:3,1:3)

      A_n(1,4:6)= (/0.,c2*n_loc(3),-c2*n_loc(2)/)
      A_n(2,4:6)= (/-c2*n_loc(3),0.,c2*n_loc(1)/)
      A_n(3,4:6)= (/c2*n_loc(2),-c2*n_loc(1),0./)
      A_n(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_n(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_n(6,1:3)= (/-n_loc(2),n_loc(1),0./)

      !positive A-Matrix-Divergence-Correction-Term
      A_p(1,8) = c_corr_c2_dielectric*n_loc(1)
      A_p(2,8) = c_corr_c2_dielectric*n_loc(2)
      A_p(3,8) = c_corr_c2_dielectric*n_loc(3)
      A_p(4,7) = c_corr*n_loc(1)
      A_p(5,7) = c_corr*n_loc(2)
      A_p(6,7) = c_corr*n_loc(3)
      A_p(7,4) = c_corr_c2_dielectric*n_loc(1)
      A_p(7,5) = c_corr_c2_dielectric*n_loc(2)
      A_p(7,6) = c_corr_c2_dielectric*n_loc(3)
      A_p(7,7) = c_corr_c_dielectric
      A_p(8,1) = c_corr*n_loc(1)
      A_p(8,2) = c_corr*n_loc(2)
      A_p(8,3) = c_corr*n_loc(3)
      A_p(8,8) = c_corr_c_dielectric

      !negative A-Matrix-Divergence-Correction-Term
      A_n(1,8) = c_corr_c2*n_loc(1)
      A_n(2,8) = c_corr_c2*n_loc(2)
      A_n(3,8) = c_corr_c2*n_loc(3)
      A_n(4,7) = c_corr*n_loc(1)
      A_n(5,7) = c_corr*n_loc(2)
      A_n(6,7) = c_corr*n_loc(3)
      A_n(7,4) = c_corr_c2*n_loc(1)
      A_n(7,5) = c_corr_c2*n_loc(2)
      A_n(7,6) = c_corr_c2*n_loc(3)
      A_n(7,7) = -c_corr_c
      A_n(8,1) = c_corr*n_loc(1)
      A_n(8,2) = c_corr*n_loc(2)
      A_n(8,3) = c_corr*n_loc(3)
      A_n(8,8) = -c_corr_c

      F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
    END DO
  END DO
END IF

END SUBROUTINE RiemannDielectricInterFace


SUBROUTINE RiemannDielectricInterFace2(Nloc,F,U_L,U_R,nv,Dielectric_Master)
!===================================================================================================================================
! Computes the numerical flux
! Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
! master is DIELECTRIC and slave PHYSICAL: A+(Eps0,Mu0) and A-(EpsR,MuR)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Globals_Vars  ,ONLY: c,c2
USE MOD_Equation_Vars ,ONLY: CentralFlux,eta_c,c_corr,c_corr_c,c_corr_c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                               :: Nloc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc),INTENT(IN) :: U_L,U_R
REAL,INTENT(IN)                                  :: nv(3,0:Nloc,0:Nloc)
REAL,INTENT(IN)                                  :: Dielectric_Master(0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                                 :: F(PP_nVar,0:Nloc,0:Nloc)
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                             :: n_loc(3),A_p(8,8),A_n(8,8)
INTEGER                                          :: Count_1,Count_2
REAL                                             :: eta_c_dielectric,c_dielectric,c2_dielectric
REAL                                             :: c_corr_c_dielectric,c_corr_c2_dielectric
!===================================================================================================================================

IF(CentralFlux)THEN
  CALL abort(__STAMP__,'central flux for dielectric media not implemented!')
ELSE
  ! Gauss point i,j
  DO Count_2=0,Nloc
    DO Count_1=0,Nloc
      n_loc(:)=nv(:,Count_1,Count_2)

      ! set dielectric values
      c_dielectric         = c*Dielectric_Master(Count_1,Count_2)                       !            c/sqrt(EpsR*MuR)
      c_corr_c_dielectric  = c_corr_c*Dielectric_Master(Count_1,Count_2)                !        chi*c/sqrt(EpsR*MuR)
      eta_c_dielectric     = c_corr_c_dielectric - c_dielectric                         ! (chi - 1.)*c/sqrt(EpsR*MuR)
      c2_dielectric        = c_dielectric*c_dielectric                                  !             c**2/(EpsR*MuR)
      c_corr_c2_dielectric = c_corr * c2_dielectric                                     !         chi*c**2/(EpsR*MuR)

      A_n(7,1:3)=0.
      A_n(1:3,7)=0.
      A_n(8,4:7)=0.
      A_n(4:7,8)=0.

      A_p(7,1:3)=0.
      A_p(1:3,7)=0.
      A_p(8,4:7)=0.
      A_p(4:7,8)=0.

      !D-Teilmatrix: Since chi and gamma is equal we
      ! consider D(chi,gamma) = D(gamma,chi)
      ! ATTENTION: if chi .ne. gamma this have to be changed.
      ! Then we need D_1 and D_2 (see commented section below)
      A_p(1,1) = c + n_loc(1)*n_loc(1)*eta_c   !  D(1,1)=(1.+n_loc(1)*n_loc(1)*(eta-1.))*c
      A_p(1,2) = n_loc(1)*n_loc(2)*eta_c       !  D(1,2)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(1,3) = n_loc(1)*n_loc(3)*eta_c       !  D(1,3)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(2,1) = A_p(1,2)                      !  D(2,1)=n_loc(1)*n_loc(2)*(eta-1)*c
      A_p(2,2) = c + n_loc(2)*n_loc(2)*eta_c   !  D(2,2)=(1.+n_loc(2)*n_loc(2)*(eta-1.))*c
      A_p(2,3) = n_loc(2)*n_loc(3)*eta_c       !  D(2,3)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_p(3,1) = A_p(1,3)                      !  D(3,1)=n_loc(1)*n_loc(3)*(eta-1)*c
      A_p(3,2) = A_p(2,3)                      !  D(3,2)=n_loc(2)*n_loc(3)*(eta-1)*c
      A_p(3,3) = c+n_loc(3)*n_loc(3)*eta_c     !  D(3,3)=(1.+n_loc(3)*n_loc(3)*(mu-1.))*c

      ! epsilon-Teilmatrix
      !E_trans=transpose(E)
      A_p(1,4:6)= (/0.,c2*n_loc(3),-c2*n_loc(2)/)
      A_p(2,4:6)= (/-c2*n_loc(3),0.,c2*n_loc(1)/)
      A_p(3,4:6)= (/c2*n_loc(2),-c2*n_loc(1),0./)
      A_p(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_p(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_p(6,1:3)= (/-n_loc(2),n_loc(1),0./)

      !composition of the Matrix
      !positive A-Matrx
      A_p(4:6,4:6)=A_p(1:3,1:3)

      !negative A-Matrix
      A_n(1,1)= -( c_dielectric + n_loc(1)*n_loc(1)*eta_c_dielectric )!D(1,1)=-1.*((1.+n_loc(1)*n_loc(1)*(chi-1.))*c/sqrt(EpsR*MuR))
      A_n(1,2)= -( n_loc(1)*n_loc(2)*eta_c_dielectric                )!D(1,2)=-1.*(n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(1,3)= -( n_loc(1)*n_loc(3)*eta_c_dielectric                )!D(1,3)=-1.*(n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(2,1)= -( A_n(1,2)                                          )!D(2,1)=-1.*(n_loc(1)*n_loc(2)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(2,2)= -( c_dielectric + n_loc(2)*n_loc(2)*eta_c_dielectric )!D(2,2)=-1.*((1./sqrt(EpsR*MuR)+n_loc(2)*n_loc(2)*(chi-1.))*c)
      A_n(2,3)= -( n_loc(2)*n_loc(3)*eta_c_dielectric                )!D(2,3)=-1.*(n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(3,1)= -( A_n(1,3)                                          )!D(3,1)=-1.*(n_loc(1)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(3,2)= -( A_n(2,3)                                          )!D(3,2)=-1.*(n_loc(2)*n_loc(3)*(chi-1.)*c/sqrt(EpsR*MuR)     )
      A_n(3,3)= -( c_dielectric+n_loc(3)*n_loc(3)*eta_c_dielectric   )!D(3,3)=-1.*((1.+n_loc(3)*n_loc(3)*(mu-1.))*c/sqrt(EpsR*MuR) )

      A_n(1,4:6)= (/0.,c2_dielectric*n_loc(3),-c2_dielectric*n_loc(2)/)
      A_n(2,4:6)= (/-c2_dielectric*n_loc(3),0.,c2_dielectric*n_loc(1)/)
      A_n(3,4:6)= (/c2_dielectric*n_loc(2),-c2_dielectric*n_loc(1),0./)
      A_n(4,1:3)= (/0.,-n_loc(3),n_loc(2)/)
      A_n(5,1:3)= (/n_loc(3),0.,-n_loc(1)/)
      A_n(6,1:3)= (/-n_loc(2),n_loc(1),0./)
      A_n(4:6,4:6)=A_n(1:3,1:3)


     ! !positive A-Matrix-Divergence-Correction-Term
      A_p(1,8) = c_corr_c2*n_loc(1)
      A_p(2,8) = c_corr_c2*n_loc(2)
      A_p(3,8) = c_corr_c2*n_loc(3)
      A_p(4,7) = c_corr*n_loc(1)
      A_p(5,7) = c_corr*n_loc(2)
      A_p(6,7) = c_corr*n_loc(3)
      A_p(7,4) = c_corr_c2*n_loc(1)
      A_p(7,5) = c_corr_c2*n_loc(2)
      A_p(7,6) = c_corr_c2*n_loc(3)
      A_p(7,7) = c_corr_c
      A_p(8,1) = c_corr*n_loc(1)
      A_p(8,2) = c_corr*n_loc(2)
      A_p(8,3) = c_corr*n_loc(3)
      A_p(8,8) = c_corr_c

      !negative A-Matrix-Divergence-Correction-Term
      A_n(1,8) = c_corr_c2_dielectric*n_loc(1)
      A_n(2,8) = c_corr_c2_dielectric*n_loc(2)
      A_n(3,8) = c_corr_c2_dielectric*n_loc(3)
      A_n(4,7) = c_corr*n_loc(1)
      A_n(5,7) = c_corr*n_loc(2)
      A_n(6,7) = c_corr*n_loc(3)
      A_n(7,4) = c_corr_c2_dielectric*n_loc(1)
      A_n(7,5) = c_corr_c2_dielectric*n_loc(2)
      A_n(7,6) = c_corr_c2_dielectric*n_loc(3)
      A_n(7,7) = -c_corr_c_dielectric
      A_n(8,1) = c_corr*n_loc(1)
      A_n(8,2) = c_corr*n_loc(2)
      A_n(8,3) = c_corr*n_loc(3)
      A_n(8,8) = -c_corr_c_dielectric

      F(:,Count_1,Count_2)=0.5*(MATMUL(A_n,U_R(:,Count_1,Count_2))+MATMUL(A_p,U_L(:,Count_1,Count_2)))
    END DO
  END DO
END IF

END SUBROUTINE RiemannDielectricInterFace2


END MODULE MOD_Riemann
