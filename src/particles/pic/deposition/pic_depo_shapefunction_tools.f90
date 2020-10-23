!==================================================================================================================================
! Copyright (c) 2010 - 2018 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
!
! This file is part of PICLas (gitlab.com/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "piclas.h"

MODULE MOD_PICDepo_Shapefunction_Tools
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
INTERFACE calcSfSource
  MODULE PROCEDURE calcSfSource
END INTERFACE

PUBLIC:: calcSfSource
!===================================================================================================================================

CONTAINS

SUBROUTINE calcSfSource(SourceSize_in,ChargeMPF,PartPos,PartID,PartVelo)
!============================================================================================================================
! deposit charges on DOFs via shapefunction including periodic displacements and mirroring
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars,           ONLY:DepositionType
USE MOD_Globals
!USE MOD_Particle_Mesh_Vars,     ONLY:casematrix,NbrOfCases
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: SourceSize_in
REAL, INTENT(IN)                 :: ChargeMPF,PartPos(3)
INTEGER, INTENT(IN)              :: PartID
REAL, INTENT(IN), OPTIONAL       :: PartVelo(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!#if ((USE_HDG) && (PP_nVar==1))
!yes, PartVelo and SourceSize_in are not used, but the subroutine-call and -head would be ugly with the preproc-flags...
!INTEGER, PARAMETER               :: SourceSize=1
!REAL                             :: Fac(4:4), Fac2(4:4)
!#else
INTEGER                          :: SourceSize
REAL                             :: Fac(4-SourceSize_in+1:4), Fac2(4-SourceSize_in+1:4)
!#endif
!INTEGER                          :: iCase, ind
!REAL                             :: ShiftedPart(1:3), caseShiftedPart(1:3), n_loc(1:3)
!----------------------------------------------------------------------------------------------------------------------------------
!#if !((USE_HDG) && (PP_nVar==1))
SourceSize=SourceSize_in
!#endif
IF (SourceSize.EQ.1) THEN
  Fac2= ChargeMPF
!#if !((USE_HDG) && (PP_nVar==1))
ELSE IF (SourceSize.EQ.4) THEN
  Fac2(1:3) = PartVelo*ChargeMPF
  Fac2(4)= ChargeMPF
!#endif
ELSE
  CALL abort(&
__STAMP__ &
,'SourceSize has to be either 1 or 4!',SourceSize)
END IF

!  DO iCase = 1, NbrOfCases
!    DO ind = 1,3
!      ShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
!        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
!    END DO
    Fac = Fac2
  SELECT CASE(TRIM(DepositionType))
  CASE('shape_function')
    CALL depoChargeOnDOFs_sf(PartPos,SourceSize,Fac)
  CASE('shape_function_cc')
    CALL depoChargeOnDOFs_sfChargeCon(PartPos,SourceSize,Fac)
  CASE('shape_function_adaptive')
    CALL depoChargeOnDOFs_sfAdaptive(PartPos,SourceSize,Fac,PartID)
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
        'Unknown ShapeFunction Method!')
  END SELECT

!  END DO ! iCase (periodicity)

END SUBROUTINE calcSfSource


SUBROUTINE depoChargeOnDOFs_sfNew(Position,SourceSize,Fac)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars,           ONLY:PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf
USE MOD_Mesh_Vars,              ONLY:nElems,offSetElem
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,ElemBaryNgeo, FIBGM_offsetElem, FIBGM_nElems, FIBGM_Element, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo
USE MOD_Preproc
USE MOD_Mesh_Tools,             ONLY: GetCNElemID
#if USE_MPI
USE MOD_PICDepo_Vars,           ONLY:PartSource_Shared_Win
USE MOD_MPI_Shared_Vars,        ONLY:nComputeNodeTotalElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
!#if ((USE_HDG) && (PP_nVar==1))
!REAL, INTENT(IN)                 :: Fac(4:4)
!#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
!#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: globElemID, CNElemID
REAL                             :: radius2, S, S1
REAL                             :: PartSourceLoc(4-SourceSize+1:4,0:PP_N,0:PP_N,0:PP_N)
INTEGER                          :: PartSourceSize, PartSourceSizeTarget, Request
INTEGER                          :: expo
#if USE_MPI
LOGICAL                          :: chargedone(1:nComputeNodeTotalElems)
#else
LOGICAL                          :: chargedone(1:nElems)
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
PartSourceSize =  SourceSize*(PP_N+1)**3
#if USE_MPI
PartSourceSizeTarget = 4*(PP_N+1)**3*nComputeNodeTotalElems
#else
PartSourceSizeTarget = 4*(PP_N+1)**3*nElems
#endif /*USE_MPI*/
chargedone(:) = .FALSE.
!-- determine which background mesh cells (and interpolation points within) need to be considered
kmax = CEILING((Position(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
kmax = MIN(kmax,GEO%FIBGMimax)
kmin = FLOOR((Position(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
kmin = MAX(kmin,GEO%FIBGMimin)
lmax = CEILING((Position(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
lmax = MIN(lmax,GEO%FIBGMjmax)
lmin = FLOOR((Position(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
lmin = MAX(lmin,GEO%FIBGMjmin)
mmax = CEILING((Position(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
mmax = MIN(mmax,GEO%FIBGMkmax)
mmin = FLOOR((Position(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
mmin = MAX(mmin,GEO%FIBGMkmin)
DO kk = kmin,kmax
  DO ll = lmin, lmax
    DO mm = mmin, mmax
      !--- go through all mapped elements not done yet
      DO ppp = 1,FIBGM_nElems(kk,ll,mm)
        globElemID = FIBGM_Element(FIBGM_offsetElem(kk,ll,mm)+ppp)
        CNElemID = GetCNElemID(globElemID)
        IF (chargedone(CNElemID)) CYCLE
        IF (VECNORM(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        ! loadbalance for halo region?
        IF (((globElemID-offSetElem).GE.1).AND.(globElemID-offSetElem).LE.nElems) &
          nDeposPerElem(globElemID-offSetElem)=nDeposPerElem(globElemID-offSetElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        PartSourceLoc = 0.0
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SUM((Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))**2.)
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO
            IF (SourceSize.EQ.1) THEN
              PartSourceLoc(4,k,l,m) = PartSourceLoc(4,k,l,m) + Fac(4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
            ELSE IF (SourceSize.EQ.4) THEN
              PartSourceLoc(1:4,k,l,m) = PartSourceLoc(1:4,k,l,m) + Fac(1:4) * S1
!#endif
            END IF
          END IF
        END DO; END DO; END DO
        chargedone(CNElemID) = .TRUE.
#if USE_MPI
!        CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,0,MPI_INFO_NULL,PartSource_Shared_Win, IERROR)
        CALL MPI_RGet_accumulate(PartSourceLoc(4-SourceSize+1:4,:,:,:),            PartSourceSize       ,MPI_DOUBLE_PRECISION, &
                                 PartSource(   4-SourceSize+1  ,0,0,0,globElemID), PartSourceSizeTarget, MPI_DOUBLE_PRECISION, 0, &
            INT((4-SourceSize+1)*(PP_N+1)**3*(globElemID-1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND, &
            PartSourceSize, MPI_DOUBLE_PRECISION, MPI_SUM, PartSource_Shared_Win,Request, IERROR)
!        CALL MPI_WAIT(Request, MPI_STATUS_IGNORE, IERROR)
!         CALL MPI_WIN_UNLOCK(0,PartSource_Shared_Win, IERROR)
!        PartSource(4-SourceSize+1:4,:,:,:,globElemID) = PartSource(4-SourceSize+1:4,:,:,:,globElemID) &
!            + PartSourceLoc(4-SourceSize+1:4,:,:,:)
!        CALL MPI_Win_flush(0,PartSource_Shared_Win, IERROR)
#else
        PartSource(4-SourceSize+1:4,:,:,:,globElemID) = PartSource(4-SourceSize+1:4,:,:,:,globElemID) &
            + PartSourceLoc(4-SourceSize+1:4,:,:,:)
#endif
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

END SUBROUTINE depoChargeOnDOFs_sfNew


SUBROUTINE depoChargeOnDOFs_sf(Position,SourceSize,Fac)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars,           ONLY:r_sf,r2_sf,r2_sf_inv,alpha_sf,PartSource
USE MOD_Mesh_Vars,              ONLY:nElems,offSetElem
USE MOD_Particle_Mesh_Vars,     ONLY:GEO,ElemBaryNgeo,FIBGM_offsetElem,FIBGM_nElems,FIBGM_Element,Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo
USE MOD_Preproc
USE MOD_Mesh_Tools,             ONLY:GetCNElemID
#if USE_MPI
USE MOD_MPI_Shared_Vars,        ONLY:nComputeNodeTotalElems
USE MOD_PICDepo_Vars,           ONLY:PartSourceProc
USE MOD_PICDepo_Vars,           ONLY:SendElemShapeID
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
!#if ((USE_HDG) && (PP_nVar==1))
!REAL, INTENT(IN)                 :: Fac(4:4)
!#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
!#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: globElemID, CNElemID
REAL                             :: radius2, S, S1
INTEGER                          :: expo, nUsedElems
#if USE_MPI
LOGICAL                          :: chargedone(1:nComputeNodeTotalElems)
!INTEGER                          :: usedElems(   nComputeNodeTotalElems)
#else
LOGICAL                          :: chargedone(1:nElems)
!INTEGER                          :: usedElems(   nElems)
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
chargedone(:) = .FALSE.
nUsedElems = 0
!-- determine which background mesh cells (and interpolation points within) need to be considered
kmax = CEILING((Position(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
kmax = MIN(kmax,GEO%FIBGMimax)
kmin = FLOOR((Position(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
kmin = MAX(kmin,GEO%FIBGMimin)
lmax = CEILING((Position(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
lmax = MIN(lmax,GEO%FIBGMjmax)
lmin = FLOOR((Position(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
lmin = MAX(lmin,GEO%FIBGMjmin)
mmax = CEILING((Position(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
mmax = MIN(mmax,GEO%FIBGMkmax)
mmin = FLOOR((Position(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
mmin = MAX(mmin,GEO%FIBGMkmin)
DO kk = kmin,kmax
  DO ll = lmin, lmax
    DO mm = mmin, mmax
      !--- go through all mapped elements not done yet
      DO ppp = 1,FIBGM_nElems(kk,ll,mm)
        globElemID = FIBGM_Element(FIBGM_offsetElem(kk,ll,mm)+ppp)
        CNElemID = GetCNElemID(globElemID)
        IF (chargedone(CNElemID)) CYCLE
        IF (VECNORM(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        IF (((globElemID-offSetElem).GE.1).AND.(globElemID-offSetElem).LE.nElems) &
          nDeposPerElem(globElemID-offSetElem)=nDeposPerElem(globElemID-offSetElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SUM((Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))**2.)
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
!            nUsedElems = nUsedElems + 1
!            usedElems(nUsedElems) = CNElemID
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO

#if USE_MPI
            IF (((globElemID-offSetElem).GE.1).AND.(globElemID-offSetElem).LE.nElems) THEN
#endif /*USE_MPI*/
              IF (SourceSize.EQ.1) THEN
                PartSource(4,k,l,m, CNElemID) = PartSource(4,k,l,m, CNElemID) + Fac(4) * S1
  !#if !((USE_HDG) && (PP_nVar==1))
              ELSE IF (SourceSize.EQ.4) THEN
                PartSource(1:4,k,l,m, CNElemID) = PartSource(1:4,k,l,m, CNElemID) + Fac(1:4) * S1
  !#endif
              END IF
#if USE_MPI
            ELSE
              IF (SourceSize.EQ.1) THEN
                PartSourceProc(4,k,l,m, SendElemShapeID(CNElemID)) =  &
                    PartSourceProc(4,k,l,m, SendElemShapeID(CNElemID)) + Fac(4) * S1
  !#if !((USE_HDG) && (PP_nVar==1))
              ELSE IF (SourceSize.EQ.4) THEN
                PartSourceProc(1:4,k,l,m, SendElemShapeID(CNElemID)) = &
                    PartSourceProc(1:4,k,l,m, SendElemShapeID(CNElemID)) + Fac(1:4) * S1
  !#endif
              END IF
            END IF
#endif /*USE_MPI*/
          END IF
        END DO; END DO; END DO
        chargedone(CNElemID) = .TRUE.
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk
END SUBROUTINE depoChargeOnDOFs_sf


SUBROUTINE depoChargeOnDOFs_sfChargeCon(Position,SourceSize,Fac)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_PICDepo_Vars,           ONLY:r_sf, r2_sf, r2_sf_inv,alpha_sf,PartSource,w_sf
USE MOD_Mesh_Vars,              ONLY:nElems, offSetElem
USE MOD_Particle_Mesh_Vars,     ONLY:GEO, ElemBaryNgeo, FIBGM_offsetElem, FIBGM_nElems, FIBGM_Element, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo, ElemsJ
USE MOD_Preproc
USE MOD_Mesh_Tools,             ONLY:GetCNElemID
USE MOD_Interpolation_Vars,     ONLY:wGP
#if USE_MPI
USE MOD_MPI_Shared_Vars,        ONLY:nComputeNodeTotalElems
USE MOD_PICDepo_Vars,           ONLY:SendElemShapeID,PartSourceProc
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
!#if ((USE_HDG) && (PP_nVar==1))
!REAL, INTENT(IN)                 :: Fac(4:4)
!#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
!#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: firstElem,elemDone
INTEGER                          :: k, l, m
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: globElemID, CNElemID
INTEGER                          :: expo, nUsedElems, localElem
REAL                             :: radius2, S, S1
REAL                             :: totalCharge, alpha
REAL                             :: PartSourcetmp(1:4,0:PP_N,0:PP_N,0:PP_N)
#if USE_MPI
LOGICAL                          :: chargedone(1:nComputeNodeTotalElems)
#else
LOGICAL                          :: chargedone(1:nElems)
#endif /*USE_MPI*/
TYPE SPElem
  REAL, ALLOCATABLE     :: PartSourceLoc(:,:,:,:)
  INTEGER               :: globElemID
  TYPE (SPElem), POINTER :: next => null()
END TYPE
TYPE (SPElem), POINTER :: first => null()
TYPE (SPElem), POINTER :: element
!----------------------------------------------------------------------------------------------------------------------------------
chargedone(:) = .FALSE.
firstElem = .TRUE.

ALLOCATE(first)
ALLOCATE(first%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
nUsedElems = 0
totalCharge = 0.0
!-- determine which background mesh cells (and interpolation points within) need to be considered
kmax = CEILING((Position(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
kmax = MIN(kmax,GEO%FIBGMimax)
kmin = FLOOR((Position(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
kmin = MAX(kmin,GEO%FIBGMimin)
lmax = CEILING((Position(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
lmax = MIN(lmax,GEO%FIBGMjmax)
lmin = FLOOR((Position(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
lmin = MAX(lmin,GEO%FIBGMjmin)
mmax = CEILING((Position(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
mmax = MIN(mmax,GEO%FIBGMkmax)
mmin = FLOOR((Position(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
mmin = MAX(mmin,GEO%FIBGMkmin)
DO kk = kmin,kmax
  DO ll = lmin, lmax
    DO mm = mmin, mmax
      !--- go through all mapped elements not done yet
      DO ppp = 1,FIBGM_nElems(kk,ll,mm)
        globElemID = FIBGM_Element(FIBGM_offsetElem(kk,ll,mm)+ppp)
        elemDone = .FALSE.
        CNElemID = GetCNElemID(globElemID)
        localElem = globElemID-offSetElem
        IF (chargedone(CNElemID)) CYCLE
        IF (VECNORM(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        IF ((localElem.GE.1).AND.localElem.LE.nElems) nDeposPerElem(localElem)=nDeposPerElem(localElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SUM((Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))**2.)
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
            IF (.NOT.elemDone) THEN
              PartSourcetmp = 0.0
              nUsedElems = nUsedElems + 1
              elemDone = .TRUE.
            END IF
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO
            IF (SourceSize.EQ.1) THEN
              PartSourcetmp(4,k,l,m) = Fac(4) * S1
            ELSE
              PartSourcetmp(1:4,k,l,m) = Fac(1:4) * S1
            END IF
            totalCharge = totalCharge  + wGP(k)*wGP(l)*wGP(m)*PartSourcetmp(4,k,l,m)/ElemsJ(k,l,m,CNElemID)
          END IF
        END DO; END DO; END DO

        IF (elemDone) THEN
          IF (firstElem) THEN
            first%PartSourceLoc(:,:,:,:) = PartSourcetmp(:,:,:,:)
            first%globElemID = globElemID
            firstElem = .FALSE.
          ELSE
            ALLOCATE(element)
            ALLOCATE(element%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
            element%next => first%next
            first%next => element
            element%PartSourceLoc(:,:,:,:) = PartSourcetmp(:,:,:,:)
            element%globElemID = globElemID
          END IF
        END IF
        chargedone(CNElemID) = .TRUE.
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

element => first
firstElem = .TRUE.
IF (nUsedElems.GT.0) THEN
  alpha = (Fac(4)/w_sf) / totalCharge
  DO ppp=1, nUsedElems
    globElemID = element%globElemID
    localElem = globElemID-offSetElem
    CNElemID = GetCNElemID(globElemID)
#if USE_MPI
    IF (((localElem).GE.1).AND.(localElem).LE.nElems) THEN
#endif /*USE_MPI*/
      IF (SourceSize.EQ.1) THEN
        PartSource(4,:,:,:, CNElemID) = PartSource(4,:,:,:, CNElemID) + alpha*element%PartSourceLoc(4,:,:,:)
      ELSE IF (SourceSize.EQ.4) THEN
        PartSource(1:4,:,:,:, CNElemID) = PartSource(1:4,:,:,:, CNElemID) + alpha*element%PartSourceLoc(1:4,:,:,:)
      END IF
#if USE_MPI
    ELSE
      IF (SourceSize.EQ.1) THEN
        PartSourceProc(4,:,:,:, SendElemShapeID(CNElemID)) =  &
            PartSourceProc(4,:,:,:, SendElemShapeID(CNElemID)) + alpha * element%PartSourceLoc(4,:,:,:)
      ELSE IF (SourceSize.EQ.4) THEN
        PartSourceProc(1:4,:,:,:, SendElemShapeID(CNElemID)) = &
            PartSourceProc(1:4,:,:,:, SendElemShapeID(CNElemID)) + alpha * element%PartSourceLoc(1:4,:,:,:)
      END IF
    END IF
#endif /*USE_MPI*/
    first => first%next
    DEALLOCATE(element%PartSourceLoc)
    DEALLOCATE(element)
    element => first
  END DO
END IF

END SUBROUTINE depoChargeOnDOFs_sfChargeCon


SUBROUTINE depoChargeOnDOFs_sfAdaptive(Position,SourceSize,Fac,PartIdx)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_PICDepo_Vars,           ONLY:alpha_sf,PartSource,w_sf,SFElemr2_Shared
USE MOD_Mesh_Vars,              ONLY:nElems, offSetElem
USE MOD_Particle_Mesh_Vars,     ONLY:ElemBaryNgeo, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo, ElemsJ, ElemToElemMapping,ElemToElemInfo
USE MOD_Preproc
USE MOD_Mesh_Tools,             ONLY:GetCNElemID, GetGlobalElemID
USE MOD_Interpolation_Vars,     ONLY:wGP
USE MOD_Particle_Vars,          ONLY:PEM
#if USE_MPI
USE MOD_MPI_Shared_Vars,        ONLY:nComputeNodeTotalElems
USE MOD_PICDepo_Vars,           ONLY:SendElemShapeID,PartSourceProc
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,       ONLY:nDeposPerElem
#endif  /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
INTEGER, INTENT(IN)              :: SourceSize
INTEGER, INTENT(IN)              :: PartIdx
!#if ((USE_HDG) && (PP_nVar==1))
!REAL, INTENT(IN)                 :: Fac(4:4)
!#else
REAL, INTENT(IN)                 :: Fac(4-SourceSize+1:4)
!#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: firstElem,elemDone
INTEGER                          :: k, l, m
INTEGER                          :: ppp
INTEGER                          :: globElemID, CNElemID, OrigCNElemID, OrigElem
INTEGER                          :: expo, nUsedElems, localElem
REAL                             :: radius2, S, S1
REAL                             :: totalCharge, alpha
REAL                             :: PartSourcetmp(1:4,0:PP_N,0:PP_N,0:PP_N)
#if USE_MPI
LOGICAL                          :: chargedone(1:nComputeNodeTotalElems)
#else
LOGICAL                          :: chargedone(1:nElems)
#endif /*USE_MPI*/
TYPE SPElem
  REAL, ALLOCATABLE     :: PartSourceLoc(:,:,:,:)
  INTEGER               :: globElemID
  TYPE (SPElem), POINTER :: next => null()
END TYPE
TYPE (SPElem), POINTER :: first => null()
TYPE (SPElem), POINTER :: element
!----------------------------------------------------------------------------------------------------------------------------------
chargedone(:) = .FALSE.
firstElem = .TRUE.

ALLOCATE(first)
ALLOCATE(first%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
nUsedElems = 0
totalCharge = 0.0
!-- determine which background mesh cells (and interpolation points within) need to be considered
!--- go through all mapped elements not done yet
OrigElem = PEM%GlobalElemID(PartIdx)
OrigCNElemID = GetCNElemID(OrigElem)
DO ppp = 0,ElemToElemMapping(2,OrigCNElemID)  
  IF (ppp.EQ.0) THEN
    globElemID = OrigElem
  ELSE    
    globElemID = GetGlobalElemID(ElemToElemInfo(ElemToElemMapping(1,OrigCNElemID)+ppp))
  END IF 
  elemDone = .FALSE.
  CNElemID = GetCNElemID(globElemID)
  localElem = globElemID-offSetElem
  IF (chargedone(CNElemID)) CYCLE
  IF (VECNORM(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(SFElemr2_Shared(1,CNElemID)+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
  IF ((localElem.GE.1).AND.localElem.LE.nElems) nDeposPerElem(localElem)=nDeposPerElem(localElem)+1
#endif /*USE_LOADBALANCE*/
    !--- go through all gauss points
  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    !-- calculate distance between gauss and particle
    radius2 = SUM((Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))**2.)
    !-- calculate charge and current density at ip point using a shape function
    !-- currently only one shapefunction available, more to follow (including structure change)
    IF (radius2 .LE. SFElemr2_Shared(2,CNElemID)) THEN
      IF (.NOT.elemDone) THEN
        PartSourcetmp = 0.0
        nUsedElems = nUsedElems + 1
        elemDone = .TRUE.
      END IF
      S = 1. - radius2/SFElemr2_Shared(2,CNElemID)
      S1 = S*S
      DO expo = 3, alpha_sf
        S1 = S*S1
      END DO
      IF (SourceSize.EQ.1) THEN
        PartSourcetmp(4,k,l,m) = Fac(4) * S1
      ELSE
        PartSourcetmp(1:4,k,l,m) = Fac(1:4) * S1
      END IF
      totalCharge = totalCharge  + wGP(k)*wGP(l)*wGP(m)*PartSourcetmp(4,k,l,m)/ElemsJ(k,l,m,CNElemID)
    END IF
  END DO; END DO; END DO

  IF (elemDone) THEN
    IF (firstElem) THEN
      first%PartSourceLoc(:,:,:,:) = PartSourcetmp(:,:,:,:)
      first%globElemID = globElemID
      firstElem = .FALSE.
    ELSE
      ALLOCATE(element)
      ALLOCATE(element%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
      element%next => first%next
      first%next => element
      element%PartSourceLoc(:,:,:,:) = PartSourcetmp(:,:,:,:)
      element%globElemID = globElemID
    END IF
  END IF
  chargedone(CNElemID) = .TRUE.
END DO ! ppp

element => first
firstElem = .TRUE.
IF (nUsedElems.GT.0) THEN
  alpha = (Fac(4)/w_sf) / totalCharge
  DO ppp=1, nUsedElems
    globElemID = element%globElemID
    localElem = globElemID-offSetElem
    CNElemID = GetCNElemID(globElemID)
#if USE_MPI
    IF (((localElem).GE.1).AND.(localElem).LE.nElems) THEN
#endif /*USE_MPI*/
      IF (SourceSize.EQ.1) THEN
        PartSource(4,:,:,:, CNElemID) = PartSource(4,:,:,:, CNElemID) + alpha*element%PartSourceLoc(4,:,:,:)
      ELSE IF (SourceSize.EQ.4) THEN
        PartSource(1:4,:,:,:, CNElemID) = PartSource(1:4,:,:,:, CNElemID) + alpha*element%PartSourceLoc(1:4,:,:,:)
      END IF
#if USE_MPI
    ELSE
      IF (SourceSize.EQ.1) THEN
        PartSourceProc(4,:,:,:, SendElemShapeID(CNElemID)) =  &
            PartSourceProc(4,:,:,:, SendElemShapeID(CNElemID)) + alpha * element%PartSourceLoc(4,:,:,:)
      ELSE IF (SourceSize.EQ.4) THEN
        PartSourceProc(1:4,:,:,:, SendElemShapeID(CNElemID)) = &
            PartSourceProc(1:4,:,:,:, SendElemShapeID(CNElemID)) + alpha * element%PartSourceLoc(1:4,:,:,:)
      END IF
    END IF
#endif /*USE_MPI*/
    first => first%next
    DEALLOCATE(element%PartSourceLoc)
    DEALLOCATE(element)
    element => first
  END DO
END IF

END SUBROUTINE depoChargeOnDOFs_sfAdaptive


END MODULE MOD_PICDepo_Shapefunction_Tools
