#include "boltzplatz.h"

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
#ifndef PP_HDG
INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE

PUBLIC::FillFlux
!===================================================================================================================================

CONTAINS

SUBROUTINE FillFlux(t,tDeriv,Flux_Master,Flux_Slave,U_master,U_slave,doMPISides)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_GLobals
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY:nSides,nBCSides
USE MOD_Riemann,         ONLY:Riemann,RiemannPML
USE MOD_Riemann,         ONLY:RiemannDielectric,RiemannDielectricInterFace,RiemannDielectricInterFace2
USE MOD_Mesh_Vars,       ONLY:NormVec,TangVec1, tangVec2, SurfElem,Face_xGP
USE MOD_GetBoundaryFlux, ONLY:GetBoundaryFlux
USE MOD_Mesh_Vars,       ONLY:firstMPISide_MINE,lastMPISide_MINE,firstInnerSide,firstBCSide,lastInnerSide
USE MOD_PML_vars,        ONLY:isPMLFace,DoPML,isPMLFace,isPMLInterFace,PMLnVar
USE MOD_Dielectric_vars, ONLY:isDielectricFace,DoDielectric,isDielectricFace,isDielectricInterFace,isDielectricElem
USE MOD_Mesh_Vars,       ONLY: SideToElem
#ifdef maxwell
USE MOD_Riemann,         ONLY:ExactFlux
#endif /*maxwell*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  != .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
REAL,INTENT(IN)    :: t           ! time
INTEGER,INTENT(IN) :: tDeriv      ! deriv
REAL,INTENT(IN)    :: U_master(PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(IN)    :: U_slave (PP_nVar,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Flux_Master(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
REAL,INTENT(OUT)   :: Flux_Slave(1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID,ElemID
!===================================================================================================================================

! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver
IF(doMPISides)THEN 
  ! fill only flux for MINE MPISides
  ! fill only flux for MINE MPISides (where the local proc is master) 
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
  lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
  lastSideID = lastInnerSide
END IF

! Compute fluxes on PP_N, no additional interpolation required
DO SideID=firstSideID,lastSideID
  ! 0.) Sanity: It is forbidden to connect a PML to a dielectric region!
  IF(DoPML.AND.DoDielectric)THEN
    IF(isPMLFace(SideID).AND.isDielectricFace(SideID))THEN
      CALL abort(&
      __STAMP__&
      ,'It is forbidden to connect a PML to a dielectric region!')
    END IF
  END IF

  ! 1.) Check Perfectly Matched Layer
  IF(DoPML) THEN
    IF (isPMLFace(SideID))THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
                                 !     (flux-splitting!)
      CALL RiemannPML(Flux_Master(1:32,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
      CYCLE
    !ELSE ! 2.) no PML, standard flux
      !CALL Riemann(Flux_Master(1:8,:,:,SideID), U_Master(:,:,:,SideID),  U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
    END IF
  !ELSE ! 3.) no PML, standard flux
    !CALL Riemann(Flux_Master(:,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
    !CYCLE
  END IF ! DoPML

  ! 2.) Check Dielectric Medium
  IF(DoDielectric) THEN
    IF (isDielectricFace(SideID))THEN ! 1.) RiemannDielectric
      IF(isDielectricInterFace(SideID))THEN     ! a) physical <-> dielectric region: in Riemann: A+(EpsR,MuR) and A-(Eps0,Mu0)
        !ElemID = SideToElem(S2E_ELEM_ID,SideID) ! get master element ID for checking if it is in a physical or dielectric region
        !IF(isDielectricElem(ElemID))THEN        ! check if master is DIELECTRIC and slave PHYSICAL
        !  ! CAUTION: switch U_Master and U_Slave
        !  CALL RiemannDielectricInterFace(Flux_Master(1:8,:,:,SideID),U_Slave(:,:,:,SideID),U_Master(:,:,:,SideID),&
        !                                                                                                      NormVec(:,:,:,SideID))
        !                                                                                 !NormVec(:,:,:,SideID),SwapNormal=.TRUE.))
        !ELSE ! check if master is PHYSICAL and slave DIELECTRIC
        !  CALL RiemannDielectricInterFace(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
        !                                                                                                      NormVec(:,:,:,SideID))
        !END IF
        ElemID = SideToElem(S2E_ELEM_ID,SideID) ! get master element ID for checking if it is in a physical or dielectric region
        IF(isDielectricElem(ElemID))THEN        ! check if master is DIELECTRIC and slave PHYSICAL
          ! A+(Eps0,Mu0) and A-(EpsR,MuR)
          CALL RiemannDielectricInterFace2(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                                                                                                              NormVec(:,:,:,SideID))
        ELSE ! check if master is PHYSICAL and slave DIELECTRIC
          ! A+(EpsR,MuR) and A-(Eps0,Mu0)
          CALL RiemannDielectricInterFace(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),&
                                                                                                              NormVec(:,:,:,SideID))
        END IF
      ELSE ! b) dielectric region <-> dielectric region
        CALL RiemannDielectric(Flux_Master(1:8,:,:,SideID),U_Master(:,:,:,SideID),U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
      END IF
    ELSE ! 2.) no Dielectric, standard flux
      CALL Riemann(Flux_Master(1:8,:,:,SideID), U_Master(:,:,:,SideID),  U_Slave(:,:,:,SideID),NormVec(:,:,:,SideID))
    END IF
  ELSE ! 3.) no Dielectric, standard flux
    CALL Riemann(Flux_Master(:,:,:,SideID),U_Master( :,:,:,SideID),U_Slave(  :,:,:,SideID),NormVec(:,:,:,SideID))
  END IF ! DoDielectric
END DO ! SideID
  
IF(.NOT.doMPISides)THEN
  CALL GetBoundaryFlux(t,tDeriv,Flux_Master    (1:PP_nVar+PMLnVar,0:PP_N,0:PP_N,1:nBCSides) &
                               ,U_master       (1:PP_nVar        ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,NormVec        (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec1       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,TangVec2       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) &
                               ,Face_XGP       (1:3              ,0:PP_N,0:PP_N,1:nBCSides) )
END IF

! Apply surface element size
DO SideID=firstSideID,lastSideID
  DO q=0,PP_N; DO p=0,PP_N
    Flux_Master(:,p,q,SideID)=Flux_Master(:,p,q,SideID)*SurfElem(p,q,SideID)
  END DO; END DO
END DO

! copy flux from Master side to slave side, DO not change sign
Flux_slave(:,:,:,firstSideID:lastSideID) = Flux_master(:,:,:,firstSideID:lastSideID)

#ifdef maxwell
IF(DoPML) THEN
  DO SideID=firstSideID,lastSideID
    IF (isPMLFace(SideID))THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations 
      ! CAUTION: Multiplication with SurfElem is done in ExactFlux
      IF(isPMLInterFace(SideID)) CALL ExactFlux(t,tDeriv                                        &
                                               , Flux_Master(1:PP_nVar+PMLnVar,:,:,SideID)      &
                                               , Flux_Slave(1:PP_nVar+PMLnVar,:,:,SideID)       &
                                               , U_Master(:,:,:,SideID)                         &
                                               , U_Slave(:,:,:,SideID)                          &
                                               , NormVec(:,:,:,SideID)                          &
                                               , Face_xGP(1:3,:,:,SideID)                       &
                                               , SurfElem(:,:,SideID)                           )
    END IF
  END DO ! SideID
END IF                                           
#endif /*maxwell*/


END SUBROUTINE FillFlux
#endif

END MODULE MOD_FillFlux
