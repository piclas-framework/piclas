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
MODULE MOD_VolInt
!===================================================================================================================================
! Containes the different DG volume integrals
! Computes the volume integral contribution based on U and updates Ut
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
#if !(USE_HDG)
INTERFACE VolInt
  MODULE PROCEDURE VolInt_weakForm
END INTERFACE


PUBLIC::VolInt
#endif
!===================================================================================================================================


CONTAINS



#if !(USE_HDG)
SUBROUTINE VolInt_weakForm(Ut,dofirstElems)
!===================================================================================================================================
! Computes the volume integral of the weak DG form a la Kopriva
! Attention 1: 1/J(i,j,k) is not yet accounted for
! Attention 2: ut is initialized and is updated with the volume flux derivatives
!===================================================================================================================================
! MODULES
USE MOD_DG_Vars,           ONLY:D_hat
USE MOD_Mesh_Vars,         ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_PML_Vars,          ONLY: DoPML,ElemToPML,isPMLElem,U2t
USE MOD_Dielectric_Vars,   ONLY: DoDielectric,isDielectricElem
USE MOD_PreProc
USE MOD_Flux,ONLY:EvalFlux3D,EvalFlux3DDielectric                      ! computes volume fluxes in local coordinates
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
LOGICAL,INTENT(IN)                                  :: dofirstElems
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: f,g,h                ! volume fluxes at all Gauss points
REAL,DIMENSION(PP_nVar)                           :: fTilde,gTilde,hTilde ! auxiliary variables needed to store the fluxes at one GP
INTEGER                                           :: i,j,k,iElem
INTEGER                                           :: l                    ! row index for matrix vector product
INTEGER                                           :: firstElemID, lastElemID
!===================================================================================================================================

IF(dofirstElems)THEN
  firstElemID = 1
  lastElemID  = PP_nElems/2+1
ELSE ! second half of elements
  firstElemID = PP_nElems/2+2
  lastElemID  = PP_nElems
END IF

DO iElem=firstElemID,lastElemID
!DO iElem=1,PP_nElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  IF(DoDielectric)THEN
    IF(isDielectricElem(iElem)) THEN ! 1.) PML version - PML element
      CALL EvalFlux3DDielectric(iElem,f,g,h)
    ELSE
      CALL EvalFlux3D(iElem,f,g,h)
    END IF
  ELSE
    CALL EvalFlux3D(iElem,f,g,h)
  END IF

  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        fTilde=f(:,i,j,k)
        gTilde=g(:,i,j,k)
        hTilde=h(:,i,j,k)
        ! Compute the transformed fluxes with the metric terms
        ! Attention 1: we store the transformed fluxes in f,g,h again
        f(:,i,j,k) = fTilde(:)*Metrics_fTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_fTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_fTilde(3,i,j,k,iElem)
        g(:,i,j,k) = fTilde(:)*Metrics_gTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_gTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_gTilde(3,i,j,k,iElem)
        h(:,i,j,k) = fTilde(:)*Metrics_hTilde(1,i,j,k,iElem) + &
                     gTilde(:)*Metrics_hTilde(2,i,j,k,iElem) + &
                     hTilde(:)*Metrics_hTilde(3,i,j,k,iElem)
      END DO ! i
    END DO ! j
  END DO ! k


  IF(DoPML)THEN
    IF(isPMLElem(iElem)) THEN ! 1.) PML version - PML element
      DO l=0,PP_N
        DO k=0,PP_N
          DO j=0,PP_N
            DO i=0,PP_N
              ! Update the time derivative with the spatial derivatives of the transformed fluxes
              Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*f(:,l,j,k) + &
                                                      D_hat(j,l)*g(:,i,l,k) + &
                                                      D_hat(k,l)*h(:,i,j,l)
              ! Update the time derivatives of the auxiliary variables of the CFS-PML region
              U2t(1 : 3,i,j,k,ElemToPML(iElem)) = U2t(1 : 3,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(1,l,j,k) , D_hat(j,l)*g(1,i,l,k) , D_hat(k,l)*h(1,i,j,l) /)
              U2t(4 : 6,i,j,k,ElemToPML(iElem)) = U2t(4 : 6,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(2,l,j,k) , D_hat(j,l)*g(2,i,l,k) , D_hat(k,l)*h(2,i,j,l) /)
              U2t(7 : 9,i,j,k,ElemToPML(iElem)) = U2t(7 : 9,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(3,l,j,k) , D_hat(j,l)*g(3,i,l,k) , D_hat(k,l)*h(3,i,j,l) /)
              U2t(10:12,i,j,k,ElemToPML(iElem)) = U2t(10:12,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(4,l,j,k) , D_hat(j,l)*g(4,i,l,k) , D_hat(k,l)*h(4,i,j,l) /)
              U2t(13:15,i,j,k,ElemToPML(iElem)) = U2t(13:15,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(5,l,j,k) , D_hat(j,l)*g(5,i,l,k) , D_hat(k,l)*h(5,i,j,l) /)
              U2t(16:18,i,j,k,ElemToPML(iElem)) = U2t(16:18,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(6,l,j,k) , D_hat(j,l)*g(6,i,l,k) , D_hat(k,l)*h(6,i,j,l) /)
              !Phi_B
              U2t(19:21,i,j,k,ElemToPML(iElem)) = U2t(19:21,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(7,l,j,k) , D_hat(j,l)*g(7,i,l,k) , D_hat(k,l)*h(7,i,j,l) /)
              !Phi_E
              U2t(22:24,i,j,k,ElemToPML(iElem)) = U2t(22:24,i,j,k,ElemToPML(iElem)) + &
                                              (/ D_hat(i,l)*f(8,l,j,k) , D_hat(j,l)*g(8,i,l,k) , D_hat(k,l)*h(8,i,j,l) /)
              !print *,"FPt(:,i,j,k,ElemToPML(iElem))",FPt(:,i,j,k,ElemToPML(iElem))
            END DO !i
          END DO ! j
        END DO ! k
      END DO ! l
    ELSE ! 2.) PML version - physical element
      DO l=0,PP_N
        DO k=0,PP_N
          DO j=0,PP_N
            DO i=0,PP_N
              ! Update the time derivative with the spatial derivatives of the transformed fluxes
              Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*f(:,l,j,k) + &
                                                      D_hat(j,l)*g(:,i,l,k) + &
                                                      D_hat(k,l)*h(:,i,j,l)
            END DO !i
          END DO ! j
        END DO ! k
      END DO ! l
    END IF
  ELSE ! 3.) physical element
    DO l=0,PP_N
      DO k=0,PP_N
        DO j=0,PP_N
          DO i=0,PP_N
            ! Update the time derivative with the spatial derivatives of the transformed fluxes
            Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*f(:,l,j,k) + &
                                                    D_hat(j,l)*g(:,i,l,k) + &
                                                    D_hat(k,l)*h(:,i,j,l)
          END DO !i
        END DO ! j
      END DO ! k
    END DO ! l
  END IF ! isPMLElem(iElem)

  !CALL VolInt_Metrics(f,g,h,Metrics_fTilde(:,:,:,:,iElem),&
  !                          Metrics_gTilde(:,:,:,:,iElem),&
  !                          Metrics_hTilde(:,:,:,:,iElem))
  !DO k=0,PP_N
  !  DO j=0,PP_N
  !    DO i=0,PP_N
  !      Ut(:,i,j,k,iElem) = D_Hat_T(0,i)*f(:,0,j,k) + &
  !                          D_Hat_T(0,j)*g(:,i,0,k) + &
  !                          D_Hat_T(0,k)*h(:,i,j,0)
  !      DO l=1,PP_N
  !        ! Update the time derivative with the spatial derivatives of the transformed fluxes
  !        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
  !                                                D_Hat_T(l,j)*g(:,i,l,k) + &
  !                                                D_Hat_T(l,k)*h(:,i,j,l)
  !      END DO ! l
  !    END DO !i
  !  END DO ! j
  !END DO ! k
END DO ! iElem
END SUBROUTINE VolInt_weakForm
#endif

#ifdef donotcompilethis
SUBROUTINE VolInt_Metrics(f,g,h,Mf,Mg,Mh)
!===================================================================================================================================
! Compute the tranformed states for all conservative variables
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:nTotal_vol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(3,nTotal_Vol),INTENT(IN)          :: Mf,Mg,Mh             ! Metrics
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,nTotal_vol),INTENT(INOUT) :: f,g,h                ! volume fluxes at all Gauss points
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                        :: i
REAL,DIMENSION(PP_nVar)                        :: fTilde,gTilde,hTilde ! auxiliary variables needed to store the fluxes at one GP
!===================================================================================================================================
DO i=1,nTotal_Vol
  fTilde=f(:,i)
  gTilde=g(:,i)
  hTilde=h(:,i)
  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i) + &
           hTilde*Mf(3,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i) + &
           hTilde*Mg(3,i)
  h(:,i) = fTilde*Mh(1,i) + &
           gTilde*Mh(2,i) + &
           hTilde*Mh(3,i)
END DO ! i
END SUBROUTINE VolInt_Metrics
#endif /* donotcompilethis */

END MODULE MOD_VolInt
