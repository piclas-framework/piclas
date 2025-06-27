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

MODULE MOD_FillFlux
!===================================================================================================================================
! Fills the inner, periodic and bc fluxes
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
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700)) && !(USE_HDG)
INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!===================================================================================================================================

CONTAINS

SUBROUTINE FillFlux(t,tDeriv,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_GLobals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: nBCSides,N_SurfMesh
USE MOD_DG_Vars            ,ONLY: U_Surf_N,DG_Elems_slave,DG_Elems_master
USE MOD_GetBoundaryFlux    ,ONLY: GetBoundaryFlux
USE MOD_Mesh_Vars          ,ONLY: firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_PML_vars           ,ONLY: PMLnVar
USE MOD_Equation_Vars      ,ONLY: DoExactFlux,isExactFluxInterFace
#ifdef maxwell
USE MOD_Riemann            ,ONLY: ExactFlux
#endif /*maxwell*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           ! time
INTEGER,INTENT(IN) :: tDeriv      ! deriv
!REAL,INTENT(IN)    :: U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides)
!REAL,INTENT(IN)    :: U_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(OUT)   :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!REAL,INTENT(OUT)   :: Flux_Slave(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID,N_master,N_slave,N_max
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
! Set the side range according to MPI or no MPI
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides (where the local proc is master)
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID       = firstMPISide_MINE
  lastSideID        = lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID       = firstBCSide    ! include BCs for master sides
  lastSideID        = lastInnerSide
END IF

! =============================
! Workflow:
!
!  1.  compute flux for non-BC sides
!  1.1) advective flux
!  ( 1.2) viscous flux )
!  ( 1.3) add up viscous flux to Flux_master )
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem
!  4.  copy flux from Flux_master to Flux_slave
!  ( 5.  convert FV flux to DG flux at mixed interfaces )
!  6. Exact flux determination (inner BC)
!==============================

! 1. compute flux for non-BC sides: Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID_wo_BC, lastSideID

  ! Get polynomial degrees of master/slave sides
  N_master = DG_Elems_master(SideID)
  N_slave  = DG_Elems_slave (SideID)
  N_max    = MAX(N_master,N_slave)

  CALL GetSurfaceFlux(SideID, N_master, N_slave, N_max&
                    , U_Surf_N(SideID)%Flux_Master(1:PP_nVar+PMLnVar , 0:N_master , 0:N_master) &
                    , U_Surf_N(SideID)%Flux_Slave( 1:PP_nVar+PMLnVar , 0:N_slave  , 0:N_slave ) &
                    , U_Surf_N(SideID)%U_Master(   1:PP_nVar         , 0:N_master , 0:N_master) &
                    , U_Surf_N(SideID)%U_Slave(    1:PP_nVar         , 0:N_slave  , 0:N_slave ) &
                    , N_SurfMesh(SideID)%NormVec(  1:3               , 0:N_max    , 0:N_max   ) &
                    , N_SurfMesh(SideID)%SurfElem(                     0:N_max    , 0:N_max   ) )
END DO ! SideID

! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  DO SideID=1,nBCSides
    ! BC sides are always master sides
    N_master = DG_Elems_master(SideID)
    CALL GetBoundaryFlux(SideID,t,tDeriv,N_master&
                                 ,U_Surf_N(SideID)%Flux_Master(1:PP_nVar+PMLnVar , 0:N_master , 0:N_master) &
                                 ,U_Surf_N(SideID)%U_Master   (1:PP_nVar         , 0:N_master , 0:N_master) &
                                 ,N_SurfMesh(SideID)%NormVec  (1:3               , 0:N_master , 0:N_master) &
                                 ,N_SurfMesh(SideID)%TangVec1 (1:3               , 0:N_master , 0:N_master) &
                                 ,N_SurfMesh(SideID)%TangVec2 (1:3               , 0:N_master , 0:N_master) &
                                 ,N_SurfMesh(SideID)%Face_XGP (1:3               , 0:N_master , 0:N_master) )
    DO q=0,N_master; DO p=0,N_master
      U_Surf_N(SideID)%Flux_Master(:,p,q)=U_Surf_N(SideID)%Flux_Master(:,p,q)*N_SurfMesh(SideID)%SurfElem(p,q)
    END DO; END DO
  END DO
END IF

!  6. Exact flux determination (inner BC)
IF(DoExactFlux) THEN
  DO SideID=firstSideID,lastSideID
    IF (isExactFluxInterFace(SideID))THEN ! CAUTION: Multiplication with SurfElem is done in ExactFlux
      ! Get polynomial degrees of master/slave sides
      N_master = DG_Elems_master(SideID)
      N_slave  = DG_Elems_slave (SideID)
      N_max    = MAX(N_master,N_slave)
      CALL ExactFlux(SideID,t,tDeriv, N_master, N_slave, N_max&
                    , U_Surf_N(SideID)%Flux_Master(1:PP_nVar+PMLnVar , 0:N_master , 0:N_master)&
                    , U_Surf_N(SideID)%Flux_Slave( 1:PP_nVar+PMLnVar , 0:N_slave  , 0:N_slave )&
                    , U_Surf_N(SideID)%U_Master(   1:PP_nVar         , 0:N_master , 0:N_master)&
                    , U_Surf_N(SideID)%U_Slave(    1:PP_nVar         , 0:N_slave  , 0:N_slave )&
                    , N_SurfMesh(SideID)%NormVec(  1:3               , 0:N_max    , 0:N_max   )&
                    , N_SurfMesh(SideID)%Face_xGP( 1:3               , 0:N_max    , 0:N_max   )&
                    , N_SurfMesh(SideID)%SurfElem(                     0:N_max    , 0:N_max   ))
    END IF ! isExactFluxFace(SideID)
  END DO ! SideID
END IF

END SUBROUTINE FillFlux


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE GetSurfaceFlux(SideID,N_master,N_slave,N_max,Flux_Master,Flux_Slave,U_Master,U_Slave,NormVec,SurfElem)
! MODULES
USE MOD_Riemann            ,ONLY: Riemann
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
USE MOD_Interpolation_Vars ,ONLY: PREF_VDM,N_Inter
USE MOD_Interfaces_Vars    ,ONLY: InterfaceRiemann
USE MOD_PML_vars           ,ONLY: PMLnVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID
INTEGER,INTENT(IN) :: N_master
INTEGER,INTENT(IN) :: N_slave
INTEGER,INTENT(IN) :: N_max
REAL,INTENT(INOUT) :: Flux_Master(1:PP_nVar+PMLnVar,0:N_master,0:N_master)
REAL,INTENT(INOUT) :: Flux_Slave( 1:PP_nVar+PMLnVar,0:N_slave ,0:N_slave)
REAL,INTENT(INOUT) :: U_Master(   1:PP_nVar        ,0:N_master,0:N_master)
REAL,INTENT(INOUT) :: U_Slave(    1:PP_nVar        ,0:N_slave ,0:N_slave)
REAL,INTENT(INOUT) :: NormVec(    1:3              ,0:N_max   ,0:N_max)
REAL,INTENT(INOUT) :: SurfElem(                     0:N_max   ,0:N_max)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE   :: Uloc(:,:,:),Fluxloc(:,:,:),Fluxdie(:,:,:)
INTEGER            :: p,q
!===================================================================================================================================

IF(N_master.EQ.N_slave) THEN ! both sides have the same polynomial degree, nothing to be done
  CALL Riemann(N_master, Flux_Master, Flux_Slave, U_Master, U_Slave, NormVec, SideID)

  DO q=0,N_master; DO p=0,N_master
    Flux_Master(:,p,q)=Flux_Master(:,p,q)*SurfElem(p,q)
  END DO; END DO

  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC,RIEMANN_VAC2DIELECTRIC_NC)
    ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! slaves sides have already been calculated
    DO q=0,N_slave; DO p=0,N_slave
      Flux_Slave(:,p,q)=Flux_Slave(:,p,q)*SurfElem(p,q)
    END DO; END DO
  CASE DEFAULT
    ! 4. copy flux from master side to slave side: DO not change sign
    Flux_Slave(:,:,:) = Flux_Master(:,:,:)
  END SELECT

ELSEIF(N_master.GT.N_slave) THEN
  ALLOCATE(Uloc(   PP_nVar        ,0:N_master,0:N_master))
  ALLOCATE(Fluxloc(PP_nVar+PMLnVar,0:N_master,0:N_master))
  CALL ChangeBasis2D(PP_nVar, N_slave, N_master, PREF_VDM(N_slave,N_master)%Vdm, &
           U_Slave(1:PP_nVar , 0:N_slave  , 0:N_slave)   , &
              Uloc(1:PP_nVar , 0:N_master , 0:N_master))
  CALL Riemann(N_master, Flux_Master, Fluxloc, U_Master, Uloc, NormVec, SideID)

  DO q=0,N_master; DO p=0,N_master
    Flux_Master(:,p,q)=Flux_Master(:,p,q)*SurfElem(p,q)
  END DO; END DO

  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC,RIEMANN_VAC2DIELECTRIC_NC)
    ALLOCATE(Fluxdie(PP_nVar+PMLnVar,0:N_master,0:N_master))
    ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! slaves sides have already been calculated
    DO q=0,N_master; DO p=0,N_master
      Fluxloc(:,p,q)=Fluxloc(:,p,q)*SurfElem(p,q)
    END DO; END DO

    !transform the slave side to the same degree as the master: switch to Legendre basis
    CALL ChangeBasis2D(PP_nVar+PMLnVar, N_master, N_master,N_Inter(N_master)%sVdm_Leg,Fluxloc(1:PP_nVar+PMLnVar,0:N_master,0:N_master), Fluxdie(1:PP_nVar+PMLnVar,0:N_master,0:N_master))
    !Fluxdie(:, N_slave+1:N_master,         0:N_master) = 0.0 ! set unnecessary modes to zero
    !Fluxdie(:,         0:N_master, N_slave+1:N_master) = 0.0 ! set unnecessary modes to zero
    ! switch back to nodal basis
    CALL ChangeBasis2D(PP_nVar+PMLnVar, N_slave , N_slave, N_Inter(N_slave)%Vdm_Leg  ,Fluxdie(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave),Flux_Slave(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave))
  CASE DEFAULT
    ! 4. copy flux from master side to slave side: DO not change sign
    !Flux_Slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)

    !transform the slave side to the same degree as the master: switch to Legendre basis
    CALL ChangeBasis2D(PP_nVar+PMLnVar, N_master,N_master,N_Inter(N_master)%sVdm_Leg,Flux_Master(1:PP_nVar+PMLnVar,0:N_master,0:N_master), Fluxloc(1:PP_nVar+PMLnVar,0:N_master,0:N_master))
    !Fluxdie(:, N_slave+1:N_master,         0:N_master) = 0.0 ! set unnecessary modes to zero
    !Fluxdie(:,         0:N_master, N_slave+1:N_master) = 0.0 ! set unnecessary modes to zero
    ! switch back to nodal basis
    CALL ChangeBasis2D(PP_nVar+PMLnVar, N_slave ,N_slave ,N_Inter(N_slave)%Vdm_Leg  ,    Fluxloc(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave),Flux_Slave(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave))
  END SELECT
ELSE ! N_slave > N_master
  ALLOCATE(Uloc(   PP_nVar        ,0:N_slave,0:N_slave))
  ALLOCATE(Fluxloc(PP_nVar+PMLnVar,0:N_slave,0:N_slave))
  CALL ChangeBasis2D(PP_nVar, N_master, N_slave, PREF_VDM(N_master,N_slave)%Vdm, &
          U_Master(1:PP_nVar , 0:N_master , 0:N_master) , &
              Uloc(1:PP_nVar , 0:N_slave  , 0:N_slave))
  CALL Riemann(N_slave, Fluxloc, Flux_Slave, Uloc, U_Slave, NormVec, SideID)

  DO q=0,N_slave; DO p=0,N_slave
    Fluxloc(:,p,q)=Fluxloc(:,p,q)*SurfElem(p,q)
  END DO; END DO

  ALLOCATE(Fluxdie(PP_nVar+PMLnVar,0:N_slave,0:N_slave))
  !transform the slave side to the same degree as the master: switch to Legendre basis
  CALL ChangeBasis2D(PP_nVar+PMLnVar, N_slave,N_slave, N_Inter(N_slave)%sVdm_Leg,Fluxloc(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave),      Fluxdie(1:PP_nVar+PMLnVar,0:N_slave,0:N_slave))
  !Fluxdie(:, N_slave+1:N_slave,         0:N_slave) = 0.0 ! set unnecessary modes to zero
  !Fluxdie(:,         0:N_slave, N_slave+1:N_slave) = 0.0 ! set unnecessary modes to zero
  ! switch back to nodal basis
  CALL ChangeBasis2D(PP_nVar+PMLnVar,N_master,N_master,N_Inter(N_master)%Vdm_Leg,Fluxdie(1:PP_nVar+PMLnVar,0:N_master,0:N_master),Flux_Master(1:PP_nVar+PMLnVar,0:N_master,0:N_master))

  SELECT CASE(InterfaceRiemann(SideID))
  CASE(RIEMANN_DIELECTRIC2VAC_NC,RIEMANN_VAC2DIELECTRIC_NC)
    ! use non-conserving fluxes (two different fluxes for master and slave side)
    ! slaves sides have already been calculated
    DO q=0,N_slave; DO p=0,N_slave
      Flux_Slave(:,p,q)=Flux_Slave(:,p,q)*SurfElem(p,q)
    END DO; END DO
  CASE DEFAULT
    ! 4. copy flux from master side to slave side: DO not change sign
    Flux_Slave(:,:,:) = Fluxloc(:,:,:)
  END SELECT
END IF
SDEALLOCATE(Uloc)
SDEALLOCATE(Fluxloc)
SDEALLOCATE(Fluxdie)

END SUBROUTINE GetSurfaceFlux

#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400) || (PP_TimeDiscMethod==700)) && !(USE_HDG)*/

END MODULE MOD_FillFlux