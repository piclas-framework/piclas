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
SUBROUTINE VolInt_weakForm(dofirstElems)
!===================================================================================================================================
! Computes the volume integral of the weak DG form a la Kopriva
! Attention 1: 1/J(i,j,k) is not yet accounted for
! Attention 2: ut is initialized and is updated with the volume flux derivatives
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY: N_DG,DGB_N,U_N
USE MOD_Mesh_Vars          ,ONLY: N_VolMesh
USE MOD_Interpolation_Vars ,ONLY: Nmax
USE MOD_PML_Vars           ,ONLY: DoPML,ElemToPML,isPMLElem,U2t
USE MOD_Dielectric_Vars    ,ONLY: DoDielectric,isDielectricElem
USE MOD_Flux               ,ONLY: EvalFlux3D,EvalFlux3DDielectric              ! computes volume fluxes in local coordinates
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(IN)                                  :: dofirstElems
! Adds volume contribution to time derivative Ut contained in MOD_DG_Vars (=aufschmutzen!)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:NMax,0:NMax,0:NMax)      :: f,g,h                ! volume fluxes at all Gauss points
REAL,DIMENSION(PP_nVar)                           :: fTilde,gTilde,hTilde ! auxiliary variables needed to store the fluxes at one GP
INTEGER                                           :: i,j,k,iElem
INTEGER                                           :: l                    ! row index for matrix vector product
INTEGER                                           :: firstElemID, lastElemID
INTEGER                                           :: Nloc
!===================================================================================================================================

IF(dofirstElems)THEN
  firstElemID = 1
  lastElemID  = PP_nElems/2+1
ELSE ! second half of elements
  firstElemID = PP_nElems/2+2
  lastElemID  = PP_nElems
END IF

DO iElem=firstElemID,lastElemID
  Nloc = N_DG(iElem)
!DO iElem=1,PP_nElems
  ! Cut out the local DG solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components
  !ASSOCIATE( f => f(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc)  ,&
             !g => g(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc)  ,&
             !h => h(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc)  )
    IF(DoDielectric)THEN
      IF(isDielectricElem(iElem)) THEN ! 1.) PML version - PML element
        CALL EvalFlux3DDielectric(Nloc,iElem,f(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),g(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),h(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc))
      ELSE
        CALL EvalFlux3D(Nloc,iElem,f(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),g(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),h(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc))
      END IF
    ELSE
      CALL EvalFlux3D(Nloc,iElem,f(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),g(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc),h(1:PP_nVar,0:Nloc,0:Nloc,0:Nloc))
      !Fortran runtime error: Index '0' of dimension 2 of array 'f' outside of expected range (1:2)

    END IF
  !END ASSOCIATE

  DO k=0,NLoc
    DO j=0,NLoc
      DO i=0,NLoc
        fTilde(1:PP_nVar)=f(1:PP_nVar,i,j,k)
        gTilde(1:PP_nVar)=g(1:PP_nVar,i,j,k)
        hTilde(1:PP_nVar)=h(1:PP_nVar,i,j,k)
        ! Compute the transformed fluxes with the metric terms
        ! Attention 1: we store the transformed fluxes in f,g,h again
        f(1:PP_nVar,i,j,k) = fTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_fTilde(1,i,j,k) + &
                             gTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_fTilde(2,i,j,k) + &
                             hTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_fTilde(3,i,j,k)
        g(1:PP_nVar,i,j,k) = fTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_gTilde(1,i,j,k) + &
                             gTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_gTilde(2,i,j,k) + &
                             hTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_gTilde(3,i,j,k)
        h(1:PP_nVar,i,j,k) = fTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_hTilde(1,i,j,k) + &
                             gTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_hTilde(2,i,j,k) + &
                             hTilde(1:PP_nVar)*N_VolMesh(iElem)%Metrics_hTilde(3,i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k


  IF(DoPML)THEN
    IF(isPMLElem(iElem)) THEN ! 1.) PML version - PML element
      DO l=0,Nloc
        DO k=0,Nloc
          DO j=0,Nloc
            DO i=0,Nloc
              ! Update the time derivative with the spatial derivatives of the transformed fluxes
              U_N(iElem)%Ut(:,i,j,k) = U_N(iElem)%Ut(:,i,j,k) + DGB_N(Nloc)%D_hat(i,l)*f(:,l,j,k) + &
                                                                DGB_N(Nloc)%D_hat(j,l)*g(:,i,l,k) + &
                                                                DGB_N(Nloc)%D_hat(k,l)*h(:,i,j,l)
              ! Update the time derivatives of the auxiliary variables of the CFS-PML region
              U2t(1 : 3,i,j,k,ElemToPML(iElem)) = U2t(1 : 3,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(1,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(1,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(1,i,j,l) /)
              U2t(4 : 6,i,j,k,ElemToPML(iElem)) = U2t(4 : 6,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(2,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(2,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(2,i,j,l) /)
              U2t(7 : 9,i,j,k,ElemToPML(iElem)) = U2t(7 : 9,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(3,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(3,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(3,i,j,l) /)
              U2t(10:12,i,j,k,ElemToPML(iElem)) = U2t(10:12,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(4,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(4,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(4,i,j,l) /)
              U2t(13:15,i,j,k,ElemToPML(iElem)) = U2t(13:15,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(5,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(5,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(5,i,j,l) /)
              U2t(16:18,i,j,k,ElemToPML(iElem)) = U2t(16:18,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(6,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(6,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(6,i,j,l) /)
              !Phi_B
              U2t(19:21,i,j,k,ElemToPML(iElem)) = U2t(19:21,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(7,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(7,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(7,i,j,l) /)
              !Phi_E
              U2t(22:24,i,j,k,ElemToPML(iElem)) = U2t(22:24,i,j,k,ElemToPML(iElem)) + &
                                              (/ DGB_N(Nloc)%D_hat(i,l)*f(8,l,j,k) , DGB_N(Nloc)%D_hat(j,l)*g(8,i,l,k) , DGB_N(Nloc)%D_hat(k,l)*h(8,i,j,l) /)
              !print *,"FPt(:,i,j,k,ElemToPML(iElem))",FPt(:,i,j,k,ElemToPML(iElem))
            END DO !i
          END DO ! j
        END DO ! k
      END DO ! l
    ELSE ! 2.) PML version - physical element
      DO l=0,Nloc
        DO k=0,Nloc
          DO j=0,Nloc
            DO i=0,Nloc
              ! Update the time derivative with the spatial derivatives of the transformed fluxes
              U_N(iElem)%Ut(:,i,j,k) = U_N(iElem)%Ut(:,i,j,k) + DGB_N(Nloc)%D_hat(i,l)*f(:,l,j,k) + &
                                                                DGB_N(Nloc)%D_hat(j,l)*g(:,i,l,k) + &
                                                                DGB_N(Nloc)%D_hat(k,l)*h(:,i,j,l)
            END DO !i
          END DO ! j
        END DO ! k
      END DO ! l
    END IF
  ELSE ! 3.) physical element
    DO l=0,Nloc
      DO k=0,Nloc
        DO j=0,Nloc
          DO i=0,Nloc
            ! Update the time derivative with the spatial derivatives of the transformed fluxes
            U_N(iElem)%Ut(:,i,j,k) = U_N(iElem)%Ut(:,i,j,k) + DGB_N(Nloc)%D_hat(i,l)*f(:,l,j,k) + &
                                                              DGB_N(Nloc)%D_hat(j,l)*g(:,i,l,k) + &
                                                              DGB_N(Nloc)%D_hat(k,l)*h(:,i,j,l)
          END DO !i
        END DO ! j
      END DO ! k
    END DO ! l
  END IF ! isPMLElem(iElem)

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
