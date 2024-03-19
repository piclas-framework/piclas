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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif


!===================================================================================================================================
! Module for the HDG element matrices
!===================================================================================================================================
MODULE MOD_Elem_Mat
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_HDG
INTERFACE Elem_Mat
  MODULE PROCEDURE Elem_Mat
END INTERFACE

!INTERFACE BuildPrecond
!  MODULE PROCEDURE BuildPrecond
!END INTERFACE

!INTERFACE PostProcessGradient
!  MODULE PROCEDURE PostProcessGradient
!END INTERFACE

PUBLIC :: Elem_Mat
PUBLIC :: BuildPrecond
PUBLIC :: PostProcessGradientHDG
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
SUBROUTINE Elem_Mat(td_iter)
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDG_Vars
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_Equation_Vars      ,ONLY: chi
#if defined(IMPA) || defined(ROS)
USE MOD_LinearSolver_Vars  ,ONLY: DoPrintConvInfo
#else
USE MOD_TimeDisc_Vars      ,ONLY: IterDisplayStep,DoDisplayIter
#endif
USE MOD_Interpolation_Vars ,ONLY: N_Inter,NMax
USE MOD_Mesh_Vars          ,ONLY: N_VolMesh,offSetElem
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh
USE MOD_Mesh_Vars          ,ONLY: N_Mesh
#ifdef VDM_ANALYTICAL
USE MOD_Mathtools          ,ONLY: INVERSE_LU
#else
USE MOD_Basis              ,ONLY: getSPDInverse
#endif
#if defined(PARTICLES)
USE MOD_HDG_Vars           ,ONLY: UseBRElectronFluid
#endif /*defined(PARTICLES)*/
#if USE_PETSC
USE PETSc
USE MOD_Mesh_Vars          ,ONLY: SideToElem, nSides
USE MOD_Mesh_Vars          ,ONLY: BoundaryType,nSides,BC
#endif
USE MOD_Mesh_Vars          ,ONLY: ElemToSide
USE MOD_DG_Vars            ,ONLY: DG_Elems_slave,DG_Elems_master
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: td_iter
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: l,p,q,g1,g2,g3,Nloc,NSideMin
INTEGER              :: i,j,iElem, i_m,i_p,j_m,j_p
INTEGER              :: iDir,jDir
INTEGER              :: iLocSide, jLocSide
INTEGER              :: SideID(6),Flip(6),iSide
REAL                 :: TauS(2,3),Fdiag_i
REAL                 :: Dhat(nGP_vol(Nmax),nGP_vol(Nmax))
REAL                 :: Ktilde(3,3)
REAL                 :: Stmp1(nGP_vol(Nmax),nGP_face(Nmax)), Stmp2(nGP_face(Nmax),nGP_face(Nmax))
INTEGER              :: idx(3),jdx(3),gdx(3)
REAL                 :: time0, time
#if USE_PETSC
PetscErrorCode       :: ierr
INTEGER              :: iSideID,jSideID
INTEGER              :: ElemID, BCsideID
INTEGER              :: iBCSide,locBCSideID
INTEGER              :: iPETScGlobal, jPETScGlobal
INTEGER              :: locSideID
REAL                 :: intMat(nGP_face, nGP_face)
INTEGER              :: BCState
#endif
!===================================================================================================================================

#if defined(IMPA) || defined(ROS)
IF(DoPrintConvInfo)THEN
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(*,*)'HDG ELEM_MAT: Pre-compute HDG local element matrices...'
  time0=PICLASTIME()
END IF
#else
IF(DoDisplayIter)THEN
  IF(HDGDisplayConvergence.AND.(MOD(td_iter,IterDisplayStep).EQ.0)) THEN
    time0=PICLASTIME()
    SWRITE(UNIT_stdOut,'(132("-"))')
    SWRITE(*,*)'HDG ELEM_MAT: Pre-compute HDG local element matrices...'
  END IF
END IF
#endif


DO iElem=1,PP_nElems
  HDG_Vol_N(iElem)%Ehat = 0.0
  HDG_Vol_N(iElem)%Smat = 0.0
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  SideID(:)=ElemToSide(E2S_SIDE_ID,:,iElem)
  Flip(:)  =ElemToSide(E2S_FLIP,:,iElem)

  ! Loop over the Gauss points with indexes (g1,g2,g3); for each
  ! point, compute all the i,j contributions in the local matrices.
  Dhat = 0.0
  DO g3=0,Nloc
    DO g2=0,Nloc
      DO g1=0,Nloc
        ASSOCIATE( JaCon1 => N_VolMesh(iElem)%Metrics_fTilde(:,g1,g2,g3) , &
                   JaCon2 => N_VolMesh(iElem)%Metrics_gTilde(:,g1,g2,g3) , &
                   JaCon3 => N_VolMesh(iElem)%Metrics_hTilde(:,g1,g2,g3) , &
                   chi    =>               chi(iElem)%tens(:,:,g1,g2,g3) )
          Ktilde(1,1) =  SUM( JaCon1 * MATMUL( chi , JaCon1 ))
          Ktilde(2,2) =  SUM( JaCon2 * MATMUL( chi , JaCon2 ))
          Ktilde(3,3) =  SUM( JaCon3 * MATMUL( chi , JaCon3 ))
          Ktilde(2,1) =  SUM( JaCon2 * MATMUL( chi , JaCon1 ))
          Ktilde(3,1) =  SUM( JaCon3 * MATMUL( chi , JaCon1 ))
          Ktilde(3,2) =  SUM( JaCon3 * MATMUL( chi , JaCon2 ))
        END ASSOCIATE
        Ktilde(1,2)=Ktilde(2,1)
        Ktilde(1,3)=Ktilde(3,1)
        Ktilde(2,3)=Ktilde(3,2)
        !scale with omega_ijk/J
        Ktilde=(N_VolMesh(iElem)%sJ(g1,g2,g3)*N_Inter(Nloc)%wGP_vol(index_3to1(g1,g2,g3,Nloc)) )*Ktilde

        ! scaled tau: omega* tau * SurfElem
        DO iLocSide=1,6
          ASSOCIATE( p => N_Mesh(Nloc)%VolToSideA(1,g1,g2,g3,Flip(iLocSide),iLocSide), &
                     q => N_Mesh(Nloc)%VolToSideA(2,g1,g2,g3,Flip(iLocSide),iLocSide), &
                     l1=> N_Mesh(Nloc)%VolToSideIJKA(1,g1,g2,g3,Flip(iLocSide),iLocSide), &
                     l2=> N_Mesh(Nloc)%VolToSideIJKA(2,g1,g2,g3,Flip(iLocSide),iLocSide)  )
             iSide = SideID(iLocSide)
             NSideMin = N_SurfMesh(iSide)%NSideMin
             IF(Nloc.EQ.NSideMin)THEN
               Taus(pm(iLocSide),SideDir(iLocSide))=N_Inter(Nloc)%wGP(l1)*N_Inter(Nloc)%wGP(l2)*N_SurfMesh(iSide)%SurfElemMin(p,q)
             ELSE  
               Taus(pm(iLocSide),SideDir(iLocSide))=N_Inter(Nloc)%wGP(l1)*N_Inter(Nloc)%wGP(l2)*N_SurfMesh(iSide)%SurfElem(p,q)
             END IF ! Nloc.EQ.NSideMin
           END ASSOCIATE
         END DO !iLocSide

        Taus=Taus*Tau(ielem)

        !---------------------------------------------------------------
        ! Dhat = D - B A^{-1} B^T

#if defined(PARTICLES)
        !  D  volume contribution for nonlinear stuff
        IF (UseBRElectronFluid.AND.(HDGNonLinSolver.EQ.1)) THEN
          j = index_3to1(g1,g2,g3,Nloc)
          Dhat(j,j) = Dhat(j,j) - HDG_Vol_N(iElem)%JwGP_vol( j)*NonlinVolumeFac(j,iElem)
        END IF
#endif /*defined(PARTICLES)*/
        !  D  surface contribution

        gdx=(/g1,g2,g3/)

        j = index_3to1(g1,g2,g3,Nloc)
        DO iDir=1,3
          idx = gdx
          DO l=0,Nloc
            idx(iDir) = l
            i = index_3to1(idx(1),idx(2),idx(3),Nloc)
            Dhat(i,j) = Dhat(i,j) - (TauS(1,iDir)*N_Inter(Nloc)%LL_minus(l,gdx(iDir))+TauS(2,iDir)*N_Inter(Nloc)%LL_plus(l,gdx(iDir)))
          END DO !l
        END DO !iDir


        ! [- B A^{-1} B^T]  contribution
        DO jDir=1,3
          jdx = gdx
          DO q=0,Nloc
            jdx(jDir)=q
            j = index_3to1(jdx(1),jdx(2),jdx(3),Nloc)
            DO iDir=1,3
              idx = gdx
              DO p=0,Nloc
                idx(iDir) = p
                i = index_3to1(idx(1),idx(2),idx(3),Nloc)
                Dhat(i,j) = Dhat(i,j) - Ktilde(iDir,jDir)*N_Inter(Nloc)%Domega(p,gdx(iDir))*N_Inter(Nloc)%Domega(q,gdx(jDir))
              END DO !p
            END DO !iDir
          END DO !q
        END DO !jDir
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Ehat = E - B A^{-1} C^T
        DO iDir=1,3
          ASSOCIATE(mLocSide=>dirPm2iSide(1,iDir), &
                    pLocSide=>dirPm2iSide(2,iDir) )
          ! X direction
          i_m = sindex_3to1(g1,g2,g3,mLocSide,Nloc) ! index on the side
          i_p = sindex_3to1(g1,g2,g3,pLocSide,Nloc)
          ASSOCIATE( Ehat_m => HDG_Vol_N(iElem)%Ehat(i_m,:,mLocSide) , &
                     Ehat_p => HDG_Vol_N(iElem)%Ehat(i_p,:,pLocSide) )
            !  E  contribution
            j = index_3to1(g1,g2,g3, Nloc)
            Ehat_m(j) = Ehat_m(j) + TauS(1,iDir)*N_Inter(Nloc)%wGP(gdx(iDir))*(-N_Inter(Nloc)%Lomega_m(gdx(iDir)))
            Ehat_p(j) = Ehat_p(j) + TauS(2,iDir)*N_Inter(Nloc)%wGP(gdx(iDir))*( N_Inter(Nloc)%Lomega_p(gdx(iDir)))
            !  [- B A^{-1} C^T]  contribution
            DO q=0,Nloc
              j = index_3to1( q,g2,g3, Nloc)
              Ehat_m(j) = Ehat_m(j) + Ktilde(1,iDir)*N_Inter(Nloc)%Domega(q,g1)*N_Inter(Nloc)%Lomega_m(gdx(iDir))
              Ehat_p(j) = Ehat_p(j) + Ktilde(1,iDir)*N_Inter(Nloc)%Domega(q,g1)*N_Inter(Nloc)%Lomega_p(gdx(iDir))
              j = index_3to1(g1, q,g3, Nloc)
              Ehat_m(j) = Ehat_m(j) + Ktilde(2,iDir)*N_Inter(Nloc)%Domega(q,g2)*N_Inter(Nloc)%Lomega_m(gdx(iDir))
              Ehat_p(j) = Ehat_p(j) + Ktilde(2,iDir)*N_Inter(Nloc)%Domega(q,g2)*N_Inter(Nloc)%Lomega_p(gdx(iDir))
              j = index_3to1(g1,g2, q, Nloc)
              Ehat_m(j) = Ehat_m(j) + Ktilde(3,iDir)*N_Inter(Nloc)%Domega(q,g3)*N_Inter(Nloc)%Lomega_m(gdx(iDir))
              Ehat_p(j) = Ehat_p(j) + Ktilde(3,iDir)*N_Inter(Nloc)%Domega(q,g3)*N_Inter(Nloc)%Lomega_p(gdx(iDir))
            END DO !q
          END ASSOCIATE
          END ASSOCIATE
        END DO !iDir

        !---------------------------------------------------------------

        !---------------------------------------------------------------
        ! Smat:  C A(-1) C^T  contribution

        ASSOCIATE( SmatK => HDG_Vol_N(iElem)%Smat(:,:,:,:) )
        DO jDir=1,3
          ! TODO: it would be better to have another index to loop
          ! over PLUS and MINUS.
          ASSOCIATE( jLS_m => dirPm2iSide(1,jDir) , &
                     jLS_p => dirPm2iSide(2,jDir) )
          j_m = sindex_3to1(g1,g2,g3,jLS_m,Nloc)
          j_p = sindex_3to1(g1,g2,g3,jLS_p,Nloc)

          DO iDir=1,3
            ASSOCIATE( iLS_m => dirPm2iSide(1,iDir) , &
                       iLS_p => dirPm2iSide(2,iDir) )
            i_m = sindex_3to1(g1,g2,g3,iLS_m,Nloc)
            i_p = sindex_3to1(g1,g2,g3,iLS_p,Nloc)

            SmatK(i_m,j_m,iLS_m,jLS_m) = SmatK(i_m,j_m,iLS_m,jLS_m) &
              + Ktilde(iDir,jDir) * N_Inter(Nloc)%Lomega_m(gdx(iDir)) * N_Inter(Nloc)%Lomega_m(gdx(jDir))

            SmatK(i_p,j_m,iLS_p,jLS_m) = SmatK(i_p,j_m,iLS_p,jLS_m) &
              + Ktilde(iDir,jDir) * N_Inter(Nloc)%Lomega_p(gdx(iDir)) * N_Inter(Nloc)%Lomega_m(gdx(jDir))

            SmatK(i_m,j_p,iLS_m,jLS_p) = SmatK(i_m,j_p,iLS_m,jLS_p) &
              + Ktilde(iDir,jDir) * N_Inter(Nloc)%Lomega_m(gdx(iDir)) * N_Inter(Nloc)%Lomega_p(gdx(jDir))

            SmatK(i_p,j_p,iLS_p,jLS_p) = SmatK(i_p,j_p,iLS_p,jLS_p) &
              + Ktilde(iDir,jDir) * N_Inter(Nloc)%Lomega_p(gdx(iDir)) * N_Inter(Nloc)%Lomega_p(gdx(jDir))

            END ASSOCIATE
          END DO !iDir
          END ASSOCIATE
        END DO !jDir
        END ASSOCIATE

        !---------------------------------------------------------------

      END DO !g1
    END DO !g2
  END DO !g3

! Invert Dhat
#ifdef VDM_ANALYTICAL
  ! Computes InvDhat via analytical expression (only works for Lagrange polynomials, hence the "analytical"
  ! pre-processor flag) when Lapack fails
  ! For Bezier (Bernstein basis) polynomial: use INVERSE_LU function
  HDG_Vol_N(iElem)%InvDhat(:,:)=INVERSE_LU(Dhat(1:nGP_vol(Nloc),1:nGP_vol(Nloc)))
#else
  HDG_Vol_N(iElem)%InvDhat(:,:)=-getSPDInverse(nGP_vol(Nloc),-Dhat(1:nGP_vol(Nloc),1:nGP_vol(Nloc)))
#endif /*VDM_ANALYTICAL*/
  ! Compute for each side pair  Ehat Dhat^{-1} Ehat^T
  DO jLocSide=1,6
    !Stmp1 = TRANSPOSE( MATMUL( Ehat(:,:,jLocSide,iElem) , InvDhat(:,:,iElem) ) )
    CALL DSYMM('L','U',nGP_vol(Nloc),nGP_face(Nloc),1., &
                HDG_Vol_N(iElem)%InvDhat(:,:),nGP_vol(Nloc), &
                TRANSPOSE( HDG_Vol_N(iElem)%Ehat(:,:,jLocSide)),nGP_vol(Nloc),0., &
                Stmp1(1:nGP_vol(Nloc),1:nGP_face(Nloc)),nGP_vol(Nloc))
    ! diagonal term
    !Stmp2 = MATMUL( Ehat(:,:,jLocSide,iElem) , Stmp1 )
    CALL DGEMM('N','N',nGP_face(Nloc),nGP_face(Nloc),nGP_vol(Nloc),1., &
                        HDG_Vol_N(iElem)%Ehat(:,:,jLocSide), nGP_face(Nloc), &
                        Stmp1(1:nGP_vol(Nloc),1:nGP_face(Nloc)),nGP_vol(Nloc),0.,&
                        Stmp2(1:nGP_face(Nloc),1:nGP_face(Nloc)),nGP_face(Nloc))
    HDG_Vol_N(iElem)%Smat(:,:,jLocSide,jLocSide) = HDG_Vol_N(iElem)%Smat(:,:,jLocSide,jLocSide) + Stmp2(1:nGP_face(Nloc),1:nGP_face(Nloc))
    !standard diagonal side mass matrix Fdiag =-Tau(elem)*wGP_pq*surfelem_pq
    ! then combined with to Smat  = Smat - F
    DO q=0,Nloc; DO p=0,Nloc
      i=q*(Nloc+1)+p+1
      iSide = SideID(jLocSide)
      NSideMin = N_SurfMesh(iSide)%NSideMin
      IF(Nloc.EQ.NSideMin)THEN
        Fdiag_i = - Tau(ielem)*N_SurfMesh(SideID(jLocSide))%SurfElemMin(p,q)*N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)
      ELSE
        Fdiag_i = - Tau(ielem)*N_SurfMesh(SideID(jLocSide))%SurfElem(p,q)   *N_Inter(Nloc)%wGP(p)*N_Inter(Nloc)%wGP(q)
      END IF ! Nloc.EQ.NSideMin
      HDG_Vol_N(iElem)%Smat(i,i,jLocSide,jLocSide) = HDG_Vol_N(iElem)%Smat(i,i,jLocSide,jLocSide) -Fdiag_i
    END DO; END DO !p,q

    ! off-diagonal terms
    DO iLocSide=jLocSide+1,6
      !Stmp2 = MATMUL( Ehat(:,:,iLocSide,iElem) , Stmp1 )
      CALL DGEMM('N','N',nGP_face(Nloc),nGP_face(Nloc),nGP_vol(Nloc),1., &
                          HDG_Vol_N(iElem)%Ehat(:,:,iLocSide), nGP_face(Nloc), &
                          Stmp1(1:nGP_vol(Nloc) ,1:nGP_face(Nloc)),nGP_vol(Nloc),0.,&
                          Stmp2(1:nGP_face(Nloc),1:nGP_face(Nloc)),nGP_face(Nloc))
      ! Using the fact that Smat is symmetric
      HDG_Vol_N(iElem)%Smat(:,:,iLocSide,jLocSide) = HDG_Vol_N(iElem)%Smat(:,:,iLocSide,jLocSide) + Stmp2(1:nGP_face(Nloc),1:nGP_face(Nloc))
      HDG_Vol_N(iElem)%Smat(:,:,jLocSide,iLocSide) = HDG_Vol_N(iElem)%Smat(:,:,jLocSide,iLocSide) + TRANSPOSE(Stmp2(1:nGP_face(Nloc),1:nGP_face(Nloc)))
    END DO !iLocSide
  END DO !jLocSide

END DO !iElem

#if USE_PETSC
! Fill Smat Petsc, TODO do this without filling Smat

! Change Smat for all small mortar sides to account for the interpolation from big to small side
DO iSide=1,nSides
  IF (SmallMortarInfo(iSide).NE.0) THEN
    locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,iSide)
    iElem    = SideToElem(S2E_NB_ELEM_ID,iSide)
    IF (iElem.LT.0) CYCLE
  ELSE
    CYCLE
  END IF
  intMat = IntMatMortar(:,:,SmallMortarType(2,iSide),SmallMortarType(1,iSide))
  DO iLocSide=1,6
    Smat(:,:,iLocSide,locSideID,iElem) = MATMUL(Smat(:,:,iLocSide,locSideID,iElem),intMat)
    Smat(:,:,locSideID,iLocSide,iElem) = MATMUL(TRANSPOSE(intMat),Smat(:,:,locSideID,iLocSide,iElem))
  END DO
END DO

! Fill Dirichlet BC Smat
DO iBCSide=1,nDirichletBCSides
  BCSideID=DirichletBC(iBCSide)
  locBCSideID = SideToElem(S2E_LOC_SIDE_ID,BCSideID)
  ElemID    = SideToElem(S2E_ELEM_ID,BCSideID)
  DO iLocSide=1,6
    Smat_BC(:,:,iLocSide,iBCSide) = Smat(:,:,iLocSide,locBCSideID,ElemID)
  END DO
END DO
! Fill Smat for PETSc with remaining DOFs
DO iElem=1,PP_nElems
  DO iLocSide=1,6
    iSideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    iPETScGlobal=PETScGlobal(iSideID)
    IF (iPETScGlobal.EQ.-1) CYCLE
    DO jLocSide=1,6
      jSideID=ElemToSide(E2S_SIDE_ID,jLocSide,iElem)
      jPETScGlobal=PETScGlobal(jSideID)
      IF (iPETScGlobal.GT.jPETScGlobal) CYCLE
      IF(SetZeroPotentialDOF.AND.(iPETScGlobal.EQ.0)) THEN
        ! The first DOF is set to constant 0 -> lambda_{1,1} = 0
        Smat(:,1,jLocSide,iLocSide,iElem) = 0
        IF(jPETScGlobal.EQ.iPETScGlobal) Smat(1,1,jLocSide,iLocSide,iElem) = 1
      END IF
      PetscCallA(MatSetValuesBlocked(Smat_petsc,1,iPETScGlobal,1,jPETScGlobal,Smat(:,:,jLocSide,iLocSide,iElem),ADD_VALUES,ierr))
    END DO
  END DO
END DO
! Set Conductor matrix
DO BCsideID=1,nConductorBCsides
  jSideID=ConductorBC(BCsideID)
  iElem=SideToElem(S2E_ELEM_ID,jSideID)
  jLocSide=SideToElem(S2E_LOC_SIDE_ID,jSideID)

  BCState = BoundaryType(BC(jSideID),BC_STATE)
  jPETScGlobal=nPETScUniqueSidesGlobal-FPC%nUniqueFPCBounds+FPC%Group(BCState,2)-1
  DO iLocSide=1,6
    iSideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    iPETScGlobal=PETScGlobal(iSideID)
    DO j=2,nGP_face; DO i=1,nGP_face ! Sum up all columns
      Smat(1,i,jLocSide,iLocSide,iElem) = Smat(1,i,jLocSide,iLocSide,iElem) + Smat(j,i,jLocSide,iLocSide,iElem)
      Smat(j,i,jLocSide,iLocSide,iElem) = 0.
    END DO; END DO
    IF(MaskedSide(iSideID).EQ.2) THEN
      DO i=2,nGP_face ! Sum up all rows
        Smat(1,1,jLocSide,iLocSide,iElem) = Smat(1,1,jLocSide,iLocSide,iElem) + Smat(1,i,jLocSide,iLocSide,iElem)
        Smat(1,i,jLocSide,iLocSide,iElem) = 0.
        Smat(i,i,jLocSide,iLocSide,iElem) = 1. ! Add diagonal entries for unused DOFs
      END DO
      iPETScGlobal=nPETScUniqueSidesGlobal-FPC%nUniqueFPCBounds+FPC%Group(BCState,2)-1
    ELSEIF(iPETScGlobal.EQ.-1) THEN
      CYCLE
    END IF
    PetscCallA(MatSetValuesBlocked(Smat_petsc,1,iPETScGlobal,1,jPETScGlobal,Smat(:,:,jLocSide,iLocSide,iElem),ADD_VALUES,ierr))
  END DO
END DO
PetscCallA(MatAssemblyBegin(Smat_petsc,MAT_FINAL_ASSEMBLY,ierr))
PetscCallA(MatAssemblyEnd(Smat_petsc,MAT_FINAL_ASSEMBLY,ierr))
#endif


#if defined(IMPA) || defined(ROS)
IF(DoPrintConvInfo)THEN
  time=PICLASTIME()
  SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' HDG ELEME_MAT DONE! [',Time-time0,' sec ]'
  SWRITE(UNIT_stdOut,'(132("-"))')
END IF
#else
IF(DoDisplayIter)THEN
  IF(HDGDisplayConvergence.AND.(MOD(td_iter,IterDisplayStep).EQ.0)) THEN
    time=PICLASTIME()
    SWRITE(UNIT_stdOut,'(A,F14.2,A)') ' HDG ELEME_MAT DONE! [',Time-time0,' sec ]'
    SWRITE(UNIT_stdOut,'(132("-"))')
  END IF
END IF
#endif

CONTAINS

 PPURE FUNCTION index_3to1(i1,i2,i3,Nloc) RESULT(i)
  INTEGER, INTENT(IN) :: i1, i2, i3, Nloc
  INTEGER :: i
   i = i3*(Nloc+1)**2 + i2*(Nloc+1) + i1 + 1
 END FUNCTION index_3to1

 PPURE FUNCTION sindex_3to1(i1,i2,i3,iLocSide,Nloc) RESULT(i)
  INTEGER, INTENT(IN) :: i1, i2, i3, iLocSide, Nloc
  INTEGER :: i
  !local variables
  INTEGER :: p, q

   p = N_Mesh(Nloc)%VolToSideA(1,i1,i2,i3,Flip(iLocSide),iLocSide)
   q = N_Mesh(Nloc)%VolToSideA(2,i1,i2,i3,Flip(iLocSide),iLocSide)

   i = q*(Nloc+1) + p + 1

 END FUNCTION sindex_3to1

END SUBROUTINE Elem_Mat


SUBROUTINE BuildPrecond()
!===================================================================================================================================
! Build a block-diagonal preconditioner for the lambda system
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI            ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData,Mask_MPIsides
#endif /*USE_MPI*/
#if USE_PETSC
USE PETSc
#else
USE MOD_Mesh_Vars      ,ONLY: nSides,SideToElem,nMPIsides_YOUR,N_SurfMesh, offSetElem
USE MOD_FillMortar_HDG ,ONLY: SmallToBigMortarPrecond_HDG
USE MOD_DG_Vars        ,ONLY: DG_Elems_master, DG_Elems_slave,N_DG_Mapping
#endif
USE MOD_Interpolation_Vars ,ONLY: Nmax
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_PETSC
PetscErrorCode   :: ierr
PC               :: pc
PetscInt         :: lens(nPETScUniqueSides)
#else
INTEGER          :: ElemID, locSideID, SideID, igf
INTEGER          :: lapack_info
#endif
INTEGER          :: NSideMin
!===================================================================================================================================

#if USE_PETSC
PetscCallA(KSPGetPC(ksp,pc,ierr))
SELECT CASE(PrecondType)
CASE(0)
  PetscCallA(PCSetType(pc,PCNONE,ierr))
CASE(1)
  PetscCallA(PCSetType(pc,PCJACOBI,ierr))
CASE(2)
  PetscCallA(PCHYPRESetType(pc,PCILU,ierr))
CASE(3)
  PetscCallA(PCHYPRESetType(pc,PCSPAI,ierr))
CASE(4)
  lens=nGP_Face
  PetscCallA(PCSetType(pc,PCBJACOBI,ierr))
  PetscCallA(PCBJacobiSetLocalBlocks(pc,nPETScUniqueSides,lens,ierr))
  PetscCallA(KSPSetUp(ksp,ierr))
case(10)
  PetscCallA(PCSetType(pc,PCCHOLESKY,ierr))
case(11)
  PetscCallA(PCSetType(pc,PCLU,ierr))
END SELECT
#else
SELECT CASE(PrecondType)
CASE(0)
! do nothing
CASE(1)
  DO SideID=1,nSides
    NSideMin = N_SurfMesh(SideID)%NSideMin
    IF(.NOT.ALLOCATED(HDG_Surf_N(SideID)%Precond)) ALLOCATE(HDG_Surf_N(SideID)%Precond(nGP_face(NSideMin),nGP_face(NSideMin)))
    HDG_Surf_N(SideID)%Precond = 0.
    !master element
    locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
    IF(locSideID.NE.-1)THEN
      ElemID    = SideToElem(S2E_ELEM_ID,SideID)
      HDG_Surf_N(SideID)%Precond(:,:) = HDG_Surf_N(SideID)%Precond(:,:)+HDG_Vol_N(ElemID)%Smat(:,:,locSideID,locSideID)
    END IF !locSideID.NE.-1
    ! neighbour element
    locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    IF(locSideID.NE.-1)THEN
      ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
      HDG_Surf_N(SideID)%Precond(:,:) = HDG_Surf_N(SideID)%Precond(:,:)+HDG_Vol_N(ElemID)%Smat(:,:,locSideID,locSideID)
    END IF !locSideID.NE.-1
  END DO ! SideID=1,nSides
#if USE_MPI
  CALL Mask_MPIsides('Precond')
#endif /*USE_MPI*/
  CALL SmallToBigMortarPrecond_HDG(PrecondType) !assemble big side
  DO SideID=1,nSides-nMPIsides_YOUR
    IF(MaskedSide(SideID).GT.0)CYCLE
    NSideMin = N_SurfMesh(SideID)%NSideMin
    ! do choleski and store into Precond
    CALL DPOTRF('U',nGP_face(NSideMin),HDG_Surf_N(SideID)%Precond(:,:),nGP_face(NSideMin),lapack_info)
    IF (lapack_info .NE. 0) CALL abort(__STAMP__,'MATRIX INVERSION FAILED for PrecondType=1!')
  END DO ! SideID=1,nSides

CASE(2)
  DO SideID=1,nSides
    NSideMin = N_SurfMesh(SideID)%NSideMin
    IF(.NOT.ALLOCATED(HDG_Surf_N(SideID)%InvPrecondDiag)) ALLOCATE(HDG_Surf_N(SideID)%InvPrecondDiag(nGP_face(NSideMin)))
    HDG_Surf_N(SideID)%InvPrecondDiag = 0.
    !master element
    locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
    IF(locSideID.NE.-1)THEN
      ElemID    = SideToElem(S2E_ELEM_ID,SideID)
      IF(NSideMin.NE.NMax) CALL abort(__STAMP__,'not implemented for different polynomial degrees')
      DO igf = 1, nGP_face(NSideMin)
        HDG_Surf_N(SideID)%InvPrecondDiag(igf) = HDG_Surf_N(SideID)%InvPrecondDiag(igf)+ &
                              HDG_Vol_N(ElemID)%Smat(igf,igf,locSideID,locSideID)
      END DO ! igf
    END IF !locSideID.NE.-1
    ! neighbour element
    locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    IF(locSideID.NE.-1)THEN
      ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
      DO igf = 1, nGP_face(NSideMin)
        HDG_Surf_N(SideID)%InvPrecondDiag(igf) = HDG_Surf_N(SideID)%InvPrecondDiag(igf)+ &
                              HDG_Vol_N(ElemID)%Smat(igf,igf,locSideID,locSideID)
      END DO ! igf
    END IF !locSideID.NE.-1
  END DO ! SideID=1,nSides
#if USE_MPI
  CALL Mask_MPIsides('InvPrecondDiag')
#endif /*USE_MPI*/
  CALL SmallToBigMortarPrecond_HDG(PrecondType) !assemble big side
  !inverse of the preconditioner matrix
  DO SideID=1,nSides-nMPIsides_YOUR
    IF(MaskedSide(SideID).GT.0)CYCLE
    IF (MAXVAL(ABS(HDG_Surf_N(SideID)%InvPrecondDiag(:))).GT.1.0e-12) THEN
      HDG_Surf_N(SideID)%InvPrecondDiag(:)=1./HDG_Surf_N(SideID)%InvPrecondDiag(:)
    ELSE
      STOP 'DIAGONAL MATRIX ENTRIES <1.0e-12,  INVERSION FAILED!'
    END IF
  END DO !1,nSides-nMPIsides_YOUR
END SELECT
#endif
END SUBROUTINE BuildPrecond


SUBROUTINE PostProcessGradientHDG()
!===================================================================================================================================
! Build a block-diagonal preconditioner for the lambda system
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_HDG_Vars
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,N_VolMesh,N_Mesh,nSides,N_SurfMesh, offSetElem
USE MOD_DG_Vars            ,ONLY: DG_Elems_master,DG_Elems_slave,N_DG_Mapping,U_N
USE MOD_Interpolation_Vars ,ONLY: N_Inter,PREF_VDM,NMax
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem,Nloc,idx_m,idx_p,SideIDm,SideIDp,NSideMin,NSideMax,iSide
INTEGER          :: SideID(6),Flip(6)
INTEGER          :: q,g1,g2,g3,gdx(3),jdx(3),jDir
REAL             :: aCon(3,3),q_loc
REAL,DIMENSION(1,0:NMax,0:NMax) :: tmp2,tmp3
!===================================================================================================================================
DO iSide = 1, nSides
  NSideMin = N_SurfMesh(iSide)%NSideMin
  NSideMax = MAX(DG_Elems_master(iSide),DG_Elems_slave(iSide))
  IF(NSideMin.EQ.NSideMax)THEN
    HDG_Surf_N(iSide)%lambdaMax(1,:) = HDG_Surf_N(iSide)%lambda(1,:)
  ELSE
    ! From low to high
    tmp2(1:1,0:NSideMin,0:NSideMin) = RESHAPE(HDG_Surf_N(iSide)%lambda(1,:),(/1,NSideMin+1,NSideMin+1/))
    CALL ChangeBasis2D(1, NSideMin, NSideMax, PREF_VDM(NSideMin,NSideMax)%Vdm , tmp2(1,0:NSideMin,0:NSideMin), tmp3(1,0:NSideMax,0:NSideMax))
    HDG_Surf_N(iSide)%lambdaMax(1,:) = RESHAPE(tmp3(1,0:NSideMax,0:NSideMax),(/nGP_face(NSideMax)/))
  END IF ! NSideMin.EQ.NSideMax
END DO ! iSide = 1, nSides

DO iElem=1,PP_nElems
  U_N(iElem)%E = 0.
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  SideID(:)=ElemToSide(E2S_SIDE_ID,:,iElem)
  Flip(:)  =ElemToSide(E2S_FLIP,:,iElem)

  ! Loop over the Gauss points with indexes (g1,g2,g3); for each
  ! point, compute all the i,j contributions in the local matrices.
  DO g3=0,Nloc
    DO g2=0,Nloc
      DO g1=0,Nloc
        ! IF q is the flux, compute  K a^J
        !ASSOCIATE( JaCon1 => Metrics_fTilde(:,g1,g2,g3,iElem) , &
        !           JaCon2 => Metrics_gTilde(:,g1,g2,g3,iElem) , &
        !           JaCon3 => Metrics_hTilde(:,g1,g2,g3,iElem) , &
        !           chi    =>      chitens(:,:,g1,g2,g3,iElem) )
        !  aCon(:,1)=sJ(g1,g2,g3,iElem)*MATMUL(chi,JaCon1(:))
        !  aCon(:,2)=sJ(g1,g2,g3,iElem)*MATMUL(chi,JaCon2(:))
        !  aCon(:,3)=sJ(g1,g2,g3,iElem)*MATMUL(chi,JaCon3(:))
        !END ASSOCIATE
        ! If q is the gradient gradu_I=-K^{-1} q
        aCon(:,1)=N_VolMesh(iElem)%sJ(g1,g2,g3)*N_VolMesh(iElem)%Metrics_fTilde(:,g1,g2,g3)
        aCon(:,2)=N_VolMesh(iElem)%sJ(g1,g2,g3)*N_VolMesh(iElem)%Metrics_gTilde(:,g1,g2,g3)
        aCon(:,3)=N_VolMesh(iElem)%sJ(g1,g2,g3)*N_VolMesh(iElem)%Metrics_hTilde(:,g1,g2,g3)

        gdx=(/g1,g2,g3/)
        !---------------------------------------------------------------
        ! q =- K^{-1} ( -A^{-1} B^T *u - A^{-1}C^T *lambda )
        DO jDir=1,3
          q_loc=0.
          jdx = gdx
          DO q=0,Nloc
            jdx(jDir)=q
            q_loc=q_loc + N_Inter(Nloc)%Domega(q,gdx(jDir))*U_N(iElem)%U(1,jdx(1),jdx(2),jdx(3))
          END DO !q
          ASSOCIATE(mLocSide=>dirPm2iSide(1,jDir), &
                    pLocSide=>dirPm2iSide(2,jDir) )

            ! X direction
            ASSOCIATE(p_m => N_Mesh(Nloc)%VolToSideA(1,g1,g2,g3,Flip(mLocSide),mLocSide), &
                      q_m => N_Mesh(Nloc)%VolToSideA(2,g1,g2,g3,Flip(mLocSide),mLocSide), &
                      p_p => N_Mesh(Nloc)%VolToSideA(1,g1,g2,g3,Flip(pLocSide),pLocSide), &
                      q_p => N_Mesh(Nloc)%VolToSideA(2,g1,g2,g3,Flip(pLocSide),pLocSide), &
                      mNSideMin => N_SurfMesh(SideID(mLocSide))%NSideMin, &
                      pNSideMin => N_SurfMesh(SideID(pLocSide))%NSideMin  )
              idx_m = q_m*(Nloc+1)+p_m+1
              idx_p = q_p*(Nloc+1)+p_p+1

              IF(Nloc.EQ.mNSideMin)THEN
                q_loc = q_loc - N_Inter(Nloc)%Lomega_m(gdx(jDir))*HDG_Surf_N(SideID(mLocSide))%lambda(1,idx_m)
              ELSE
                q_loc = q_loc - N_Inter(Nloc)%Lomega_m(gdx(jDir))*HDG_Surf_N(SideID(mLocSide))%lambdaMax(1,idx_m)
              END IF ! Nloc.EQ.mNSideMin

              IF(Nloc.EQ.pNSideMin)THEN
                q_loc = q_loc - N_Inter(Nloc)%Lomega_p(gdx(jDir))*HDG_Surf_N(SideID(pLocSide))%lambda(1,idx_p)
              ELSE
                q_loc = q_loc - N_Inter(Nloc)%Lomega_p(gdx(jDir))*HDG_Surf_N(SideID(pLocSide))%lambdaMax(1,idx_p)
              END IF ! Nloc.EQ.mNSideMin
            END ASSOCIATE
          END ASSOCIATE
          U_N(iElem)%E(:,g1,g2,g3)=U_N(iElem)%E(:,g1,g2,g3)+aCon(:,jDir)*q_loc
        END DO !jDir

        !---------------------------------------------------------------
      END DO !g1
    END DO !g2
  END DO !g3

END DO !iElem

END SUBROUTINE PostProcessGradientHDG
#endif /*USE_HDG*/
END MODULE MOD_Elem_Mat
