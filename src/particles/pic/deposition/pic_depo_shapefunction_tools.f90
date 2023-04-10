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

INTERFACE SFNorm
  MODULE PROCEDURE SFNorm
END INTERFACE

PUBLIC:: calcSfSource
PUBLIC:: SFNorm
PUBLIC:: InitShapeFunctionDimensionalty
!===================================================================================================================================

CONTAINS

SUBROUTINE calcSfSource(SourceSize_in,ChargeMPF,PartPos,PartID,PartVelo)
!============================================================================================================================
! deposit charges on DOFs via shapefunction including periodic displacements and mirroring
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: DepositionType,dimFactorSF,sfDepo3D,r_sf,r2_sf,r2_sf_inv,SFElemr2_Shared,w_sf
USE MOD_PICDepo_Vars       ,ONLY: totalChargePeriodicSF,SFAdaptiveSmoothing
USE MOD_Particle_Mesh_Vars ,ONLY: NbrOfPeriodicSFCases
USE MOD_Particle_Vars      ,ONLY: PEM
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
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
!REAL                             :: Fac(4:4)
!#else
INTEGER           :: SourceSize
REAL              :: Fac(5-SourceSize_in:4)
!#endif
INTEGER           :: iCase
REAL              :: PartPosShifted(1:3)
INTEGER           :: OrigElem,OrigCNElemID
REAL              :: r_sf_tmp,r2_sf_tmp,r2_sf_inv_tmp
!----------------------------------------------------------------------------------------------------------------------------------
!#if !((USE_HDG) && (PP_nVar==1))
SourceSize=SourceSize_in
!#endif
IF (SourceSize.EQ.1) THEN
  Fac= ChargeMPF
!#if !((USE_HDG) && (PP_nVar==1))
ELSE IF (SourceSize.EQ.4) THEN
  Fac(1:3) = PartVelo*ChargeMPF
  Fac(4)= ChargeMPF
!#endif
ELSE
  CALL abort(__STAMP__,'SourceSize has to be either 1 or 4!',SourceSize)
END IF

! Check if periodic sides are present, otherwise simply use PartPos for deposition
IF(NbrOfPeriodicSFCases.GT.1)THEN

  IF(TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
    ! Get radius/radius squared/inverse radius squared for each CN element when using adaptive SF
    OrigElem     = PEM%GlobalElemID(PartID)
    OrigCNElemID = GetCNElemID(OrigElem)
    r_sf_tmp     = SFElemr2_Shared(1,OrigCNElemID)
    r2_sf_tmp    = SFElemr2_Shared(2,OrigCNElemID)
    r2_sf_inv_tmp= 1./SFElemr2_Shared(2,OrigCNElemID)
  ELSE
    r_sf_tmp     = r_sf
    r2_sf_tmp    = r2_sf
    r2_sf_inv_tmp= r2_sf_inv
  END IF ! TRIM(DepositionType).EQ.'shape_function_adaptive'

  IF(TRIM(DepositionType).EQ.'shape_function')THEN
    ! Incorporate the coefficient that considers the volume integral of the kernel
    Fac = Fac*w_sf
  ELSE
    ! Nullify for each particle
    totalChargePeriodicSF = 0.

    ! Virtually deposit each particle and add up the charge contribution
    DO iCase = 1, NbrOfPeriodicSFCases
      PartPosShifted(1:3) = GetPartPosShifted(iCase,PartPos(1:3))
      CALL calcTotalChargePeriodic_cc(PartPosShifted , Fac(4) , totalChargePeriodicSF , r_sf_tmp , r2_sf_tmp , r2_sf_inv_tmp)
    END DO

    ! Check if the charge is to be distributed over a line (1D) or area (2D)
    IF(.NOT.sfDepo3D)THEN
      totalChargePeriodicSF = totalChargePeriodicSF / dimFactorSF
    END IF

    ! Add scaling factor for charge conservation
    Fac(1:4) = Fac(1:4) * Fac(4)/totalChargePeriodicSF
  END IF ! TRIM(DepositionType).NE.'shape_function'

  ! Call adaptive periodic shape function routine with corrected charge contribution stored in Fac
  DO iCase = 1 , NbrOfPeriodicSFCases
    PartPosShifted(1:3) = GetPartPosShifted(iCase,PartPos(1:3))
    CALL depoChargeOnDOFsSF(PartPosShifted , SourceSize , Fac , r_sf_tmp , r2_sf_tmp , r2_sf_inv_tmp)
  END DO

ELSE
  ! Select deposition type
  SELECT CASE(TRIM(DepositionType))
  CASE('shape_function')
    ! Consider the integration factor w_sf for standard (uncorrected) deposition method
    CALL depoChargeOnDOFsSF(          PartPos , SourceSize , Fac*w_sf , r_sf , r2_sf , r2_sf_inv)
  CASE('shape_function_cc')
    CALL depoChargeOnDOFsSFChargeCon( PartPos , SourceSize , Fac      , r_sf , r2_sf , r2_sf_inv)
  CASE('shape_function_adaptive')
    IF(SFAdaptiveSmoothing)THEN
      ! Get radius/radius squared/inverse radius squared for each CN element when using adaptive SF
      OrigElem     = PEM%GlobalElemID(PartID)
      OrigCNElemID = GetCNElemID(OrigElem)
      r_sf_tmp     = SFElemr2_Shared(1,OrigCNElemID)
      r2_sf_tmp    = SFElemr2_Shared(2,OrigCNElemID)
      r2_sf_inv_tmp= 1./SFElemr2_Shared(2,OrigCNElemID)
      CALL depoChargeOnDOFsSFChargeCon( PartPos , SourceSize , Fac      , r_sf_tmp , r2_sf_tmp , r2_sf_inv_tmp)
    ELSE
      CALL depoChargeOnDOFsSFAdaptive(  PartPos , SourceSize , Fac      , PartID )
    END IF ! SFAdaptiveSmoothing
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
        'Unknown ShapeFunction Method!')
  END SELECT ! DepositionType
END IF ! TRIM(DepositionType).EQ.'shape_function_cc'.AND.(NbrOfPeriodicSFCases.GT.1)


END SUBROUTINE calcSfSource


#ifdef WIP
SUBROUTINE depoChargeOnDOFsSF_RGetAccumulate(Position,SourceSize,Fac)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf,ChargeSFDone
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,ElemBaryNgeo, FIBGM_offsetElem, FIBGM_nElems, FIBGM_Element, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemRadiusNGeo
USE MOD_Preproc
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: PartSource_Shared_Win
USE MOD_MPI_Shared_Vars    ,ONLY: nComputeNodeTotalElems
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
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
INTEGER                          :: expo,I
!----------------------------------------------------------------------------------------------------------------------------------
I=5-SourceSize
PartSourceSize =  SourceSize*(PP_N+1)**3
#if USE_MPI
PartSourceSizeTarget = 4*(PP_N+1)**3*nComputeNodeTotalElems
#else
PartSourceSizeTarget = 4*(PP_N+1)**3*nElems
#endif /*USE_MPI*/
ChargeSFDone(:) = .FALSE.
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
        IF (ChargeSFDone(CNElemID)) CYCLE
        IF (SFNorm(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        ! loadbalance for halo region?
        IF (((globElemID-offSetElem).GE.1).AND.(globElemID-offSetElem).LE.nElems) &
          nDeposPerElem(globElemID-offSetElem)=nDeposPerElem(globElemID-offSetElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        PartSourceLoc = 0.0
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SFRadius2(Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO
            PartSourceLoc(I:4,k,l,m) = PartSourceLoc(I:4,k,l,m) + Fac(I:4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
!#endif
          END IF
        END DO; END DO; END DO
        ChargeSFDone(CNElemID) = .TRUE.
#if USE_MPI
!        CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,0,MPI_INFO_NULL,PartSource_Shared_Win, IERROR)
        CALL MPI_RGet_accumulate(PartSourceLoc(4-SourceSize+1:4,:,:,:),            PartSourceSize       ,MPI_DOUBLE_PRECISION, &
                                 PartSource(   4-SourceSize+1  ,0,0,0,globElemID), PartSourceSizeTarget, MPI_DOUBLE_PRECISION, 0, &
            INT((4-SourceSize+1)*(PP_N+1)**3*(globElemID-1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND, &
            PartSourceSize, MPI_DOUBLE_PRECISION, MPI_SUM, PartSource_Shared_Win,Request, IERROR)
!        CALL MPI_WAIT(Request, MPI_STATUS_IGNORE, IERROR)
!         CALL MPI_WIN_UNLOCK_ALL(0,PartSource_Shared_Win, IERROR)
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

END SUBROUTINE depoChargeOnDOFsSF_RGetAccumulate
#endif /*WIP*/


SUBROUTINE depoChargeOnDOFsSF(Position, SourceSize, Fac, r_sf, r2_sf, r2_sf_inv)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: alpha_sf,ChargeSFDone
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,ElemBaryNgeo,FIBGM_offsetElem,FIBGM_nElems,FIBGM_Element,Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemRadiusNGeo
USE MOD_Preproc
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: nElems,offSetElem
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
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
REAL, INTENT(IN)                 :: r_sf       !< shape function radius
REAL, INTENT(IN)                 :: r2_sf      !< shape function radius squared
REAL, INTENT(IN)                 :: r2_sf_inv  !< inverse of shape function radius squared
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: globElemID, CNElemID
REAL                             :: radius2, S, S1
INTEGER                          :: expo, nUsedElems,I
!----------------------------------------------------------------------------------------------------------------------------------
I=5-SourceSize
ChargeSFDone(:) = .FALSE.
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
        IF (ChargeSFDone(CNElemID)) CYCLE
        IF (SFNorm(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        IF (((globElemID-offSetElem).GE.1).AND.(globElemID-offSetElem).LE.nElems) &
          nDeposPerElem(globElemID-offSetElem)=nDeposPerElem(globElemID-offSetElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SFRadius2(Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))
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
            CALL UpdatePartSource(I,k,l,m,globElemID,S1*Fac(I:4))
          END IF
        END DO; END DO; END DO
        ChargeSFDone(CNElemID) = .TRUE.
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk
END SUBROUTINE depoChargeOnDOFsSF


SUBROUTINE calcTotalChargePeriodic_cc(Position, Fac, totalCharge, r_sf, r2_sf, r2_sf_inv)
!============================================================================================================================
! sum up the charge contribution deposited by a particle and all of its periodically shifted positions for a globally fixed
! shape function radius
!============================================================================================================================
! use MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: alpha_sf,ChargeSFDone
USE MOD_Mesh_Vars          ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: GEO, ElemBaryNgeo, FIBGM_offsetElem, FIBGM_nElems, FIBGM_Element, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemRadiusNGeo, ElemsJ
USE MOD_Preproc
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Interpolation_Vars ,ONLY: wGP
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
#endif  /*USE_LOADBALANCE*/
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)                 :: Position(3)
REAL, INTENT(IN)                 :: Fac
REAL, INTENT(INOUT)              :: totalCharge
REAL, INTENT(IN)                 :: r_sf       !< shape function radius
REAL, INTENT(IN)                 :: r2_sf      !< shape function radius squared
REAL, INTENT(IN)                 :: r2_sf_inv  !< inverse of shape function radius squared
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: globElemID, CNElemID
INTEGER                          :: expo, localElem
REAL                             :: radius2, S, S1
!----------------------------------------------------------------------------------------------------------------------------------
ChargeSFDone(:) = .FALSE.

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
        localElem = globElemID-offSetElem
        IF (ChargeSFDone(CNElemID)) CYCLE
        IF (SFNorm(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        IF ((localElem.GE.1).AND.localElem.LE.nElems) nDeposPerElem(localElem)=nDeposPerElem(localElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          radius2 = SFRadius2(Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO
            totalCharge = totalCharge  + wGP(k)*wGP(l)*wGP(m)*Fac*S1/ElemsJ(k,l,m,CNElemID)
          END IF
        END DO; END DO; END DO

        ChargeSFDone(CNElemID) = .TRUE.
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

END SUBROUTINE calcTotalChargePeriodic_cc


SUBROUTINE depoChargeOnDOFsSFChargeCon(Position,SourceSize,Fac,r_sf, r2_sf, r2_sf_inv)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: alpha_sf,ChargeSFDone,PartSourceTmp
USE MOD_Mesh_Vars          ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: GEO, ElemBaryNgeo, FIBGM_offsetElem, FIBGM_nElems, FIBGM_Element, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemRadiusNGeo, ElemsJ
USE MOD_Preproc
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Interpolation_Vars ,ONLY: wGP
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
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
REAL, INTENT(IN)                 :: r_sf       !< shape function radius
REAL, INTENT(IN)                 :: r2_sf      !< shape function radius squared
REAL, INTENT(IN)                 :: r2_sf_inv  !< inverse of shape function radius squared
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
TYPE SPElem
  REAL, ALLOCATABLE     :: PartSourceLoc(:,:,:,:)
  INTEGER               :: globElemID
  TYPE (SPElem), POINTER :: next => null()
END TYPE
TYPE (SPElem), POINTER :: first => null()
TYPE (SPElem), POINTER :: element
INTEGER           :: I
!----------------------------------------------------------------------------------------------------------------------------------
I=5-SourceSize
ChargeSFDone(:) = .FALSE.
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
        IF (ChargeSFDone(CNElemID)) CYCLE
        !IF (VECNORM(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
        IF (SFNorm(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(r_sf+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
        IF ((localElem.GE.1).AND.localElem.LE.nElems) nDeposPerElem(localElem)=nDeposPerElem(localElem)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
        DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
          !-- calculate distance between gauss and particle
          !radius2 = SUM((Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))**2.)
          radius2 = SFRadius2(Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))
          !-- calculate charge and current density at ip point using a shape function
          !-- currently only one shapefunction available, more to follow (including structure change)
          IF (radius2 .LE. r2_sf) THEN
            IF (.NOT.elemDone) THEN
              PartSourceTmp = 0.0
              nUsedElems = nUsedElems + 1
              elemDone = .TRUE.
            END IF
            S = 1. - r2_sf_inv * radius2
            S1 = S*S
            DO expo = 3, alpha_sf
              S1 = S*S1
            END DO
            IF (SourceSize.EQ.1) THEN
              PartSourceTmp(4,k,l,m) = Fac(4) * S1
            ELSE
              PartSourceTmp(1:4,k,l,m) = Fac(1:4) * S1
            END IF
            totalCharge = totalCharge  + wGP(k)*wGP(l)*wGP(m)*PartSourceTmp(4,k,l,m)/ElemsJ(k,l,m,CNElemID)
          END IF
        END DO; END DO; END DO

        IF (elemDone) THEN
          IF (firstElem) THEN
            first%PartSourceLoc(:,:,:,:) = PartSourceTmp(:,:,:,:)
            first%globElemID = globElemID
            firstElem = .FALSE.
          ELSE
            ALLOCATE(element)
            ALLOCATE(element%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
            element%next => first%next
            first%next => element
            element%PartSourceLoc(:,:,:,:) = PartSourceTmp(:,:,:,:)
            element%globElemID = globElemID
          END IF
        END IF
        ChargeSFDone(CNElemID) = .TRUE.
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

element => first
firstElem = .TRUE.
IF (nUsedElems.GT.0) THEN
  alpha = Fac(4) / totalCharge
  DO ppp=1, nUsedElems
    globElemID = element%globElemID
    DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
      CALL UpdatePartSource(I,k,l,m,globElemID,alpha*element%PartSourceLoc(I:4,k,l,m))
    END DO;END DO;END DO;
    first => first%next
    DEALLOCATE(element%PartSourceLoc)
    DEALLOCATE(element)
    element => first
  END DO
END IF

IF(ASSOCIATED(first)) THEN
  SDEALLOCATE(first%PartSourceLoc)
  DEALLOCATE(first)
END IF

END SUBROUTINE depoChargeOnDOFsSFChargeCon


SUBROUTINE depoChargeOnDOFsSFAdaptive(Position,SourceSize,Fac,PartIdx)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: alpha_sf,SFElemr2_Shared,ChargeSFDone,sfDepo3D,dimFactorSF
USE MOD_Mesh_Vars          ,ONLY: offSetElem
USE MOD_Particle_Mesh_Vars ,ONLY: ElemBaryNgeo, Elem_xGP_Shared
USE MOD_Particle_Mesh_Vars ,ONLY: ElemRadiusNGeo, ElemsJ, ElemToElemMapping,ElemToElemInfo
USE MOD_Preproc
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID, GetGlobalElemID
USE MOD_Interpolation_Vars ,ONLY: wGP
USE MOD_Particle_Vars      ,ONLY: PEM
#if USE_LOADBALANCE
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
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
REAL                             :: PartSourceTmp(1:4,0:PP_N,0:PP_N,0:PP_N)
TYPE SPElem
  REAL, ALLOCATABLE     :: PartSourceLoc(:,:,:,:)
  INTEGER               :: globElemID
  TYPE (SPElem), POINTER :: next => null()
END TYPE
TYPE (SPElem), POINTER :: first => null()
TYPE (SPElem), POINTER :: element
INTEGER           :: I
!----------------------------------------------------------------------------------------------------------------------------------
I=5-SourceSize
ChargeSFDone(:) = .FALSE.
firstElem = .TRUE.

ALLOCATE(first)
ALLOCATE(first%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
nUsedElems = 0
totalCharge = 0.0
!-- determine which background mesh cells (and interpolation points within) need to be considered
!--- go through all mapped elements not done yet
OrigElem = PEM%GlobalElemID(PartIdx)
OrigCNElemID = GetCNElemID(OrigElem)
globElemID = OrigElem
DO ppp = 0,ElemToElemMapping(2,OrigCNElemID)
  IF (ppp.GT.0) globElemID = GetGlobalElemID(ElemToElemInfo(ElemToElemMapping(1,OrigCNElemID)+ppp))
  elemDone = .FALSE.
  CNElemID = GetCNElemID(globElemID)
  localElem = globElemID-offSetElem
  IF (ChargeSFDone(CNElemID)) CYCLE
  IF (SFNorm(Position(1:3)-ElemBaryNgeo(1:3,CNElemID)).GT.(SFElemr2_Shared(1,OrigCNElemID)+ElemRadiusNGeo(CNElemID))) CYCLE
#if USE_LOADBALANCE
  IF ((localElem.GE.1).AND.localElem.LE.nElems) nDeposPerElem(localElem)=nDeposPerElem(localElem)+1
#endif /*USE_LOADBALANCE*/
    !--- go through all gauss points
  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    !-- calculate distance between gauss and particle
    radius2 = SFRadius2(Position(1:3) - Elem_xGP_Shared(1:3,k,l,m,globElemID))
    !-- calculate charge and current density at ip point using a shape function
    !-- currently only one shapefunction available, more to follow (including structure change)
    IF (radius2 .LE. SFElemr2_Shared(2,OrigCNElemID)) THEN
      IF (.NOT.elemDone) THEN
        PartSourceTmp = 0.0
        nUsedElems = nUsedElems + 1
        elemDone = .TRUE.
      END IF
      S = 1. - radius2/SFElemr2_Shared(2,OrigCNElemID)
      S1 = S*S
      DO expo = 3, alpha_sf
        S1 = S*S1
      END DO
      IF (SourceSize.EQ.1) THEN
        PartSourceTmp(4,k,l,m) = Fac(4) * S1
      ELSE
        PartSourceTmp(1:4,k,l,m) = Fac(1:4) * S1
      END IF
      totalCharge = totalCharge  + wGP(k)*wGP(l)*wGP(m)*PartSourceTmp(4,k,l,m)/ElemsJ(k,l,m,CNElemID)
    END IF
  END DO; END DO; END DO

  IF (elemDone) THEN
    IF (firstElem) THEN
      first%PartSourceLoc(:,:,:,:) = PartSourceTmp(:,:,:,:)
      first%globElemID = globElemID
      firstElem = .FALSE.
    ELSE
      ALLOCATE(element)
      ALLOCATE(element%PartSourceLoc(1:4,0:PP_N,0:PP_N,0:PP_N))
      element%next => first%next
      first%next => element
      element%PartSourceLoc(:,:,:,:) = PartSourceTmp(:,:,:,:)
      element%globElemID = globElemID
    END IF
  END IF
  ChargeSFDone(CNElemID) = .TRUE.
END DO ! ppp

! Check if the charge is to be distributed over a line (1D) or area (2D)
IF(.NOT.sfDepo3D)THEN
  totalCharge = totalCharge / dimFactorSF
END IF

element => first
firstElem = .TRUE.
IF (nUsedElems.GT.0) THEN
  alpha = Fac(4) / totalCharge
  DO ppp=1, nUsedElems
    globElemID = element%globElemID
    DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
      CALL UpdatePartSource(I,k,l,m,globElemID,alpha*element%PartSourceLoc(I:4,k,l,m))
    END DO;END DO;END DO;
    first => first%next
    DEALLOCATE(element%PartSourceLoc)
    DEALLOCATE(element)
    element => first
  END DO
END IF

IF(ASSOCIATED(first)) THEN
  SDEALLOCATE(first%PartSourceLoc)
  DEALLOCATE(first)
END IF

END SUBROUTINE depoChargeOnDOFsSFAdaptive


SUBROUTINE UpdatePartSource(dim1,k,l,m,globElemID,Source)
!============================================================================================================================
! Update PartSource (if the element where the deposition occurs is mine) or PartSourceProc (if the element where the deposition
! occurs is on another processor)
!============================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_PICDepo_Vars       ,ONLY: PartSource
USE MOD_Mesh_Vars          ,ONLY: offsetElem
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
#if USE_MPI
USE MOD_PICDepo_Vars       ,ONLY: SendElemShapeID, CNRankToSendRank, ShapeMapping !PartSourceGlob
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared
!USE MOD_MPI_Shared_Vars, ONLY: myComputeNodeRank, ComputeNodeRootRank
!USE MOD_Particle_Mesh_Vars   ,ONLY: nComputeNodeElems, ElemInfo_Shared
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
INTEGER, INTENT(IN) :: dim1
INTEGER, INTENT(IN) :: k,l,m
INTEGER, INTENT(IN) :: globElemID
REAL, INTENT(IN)    :: Source(dim1:4)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: localElem, CNElemID
#if USE_MPI
INTEGER           :: ExRankID
#endif /*USE_MPI*/
!===================================================================================================================================
localElem = globElemID-offSetElem
CNElemID = GetCNElemID(globElemID)
#if USE_MPI
IF (ElemOnMyProc(globElemID)) THEN
#endif /*USE_MPI*/
  PartSource(dim1:4,k,l,m,localElem) = PartSource(dim1:4,k,l,m,localElem) + Source(dim1:4)
!#if !((USE_HDG) && (PP_nVar==1))
!#endif
#if USE_MPI
ELSE
!  IF (myComputeNodeRank.EQ.0) THEN  
!    PartSourceGlob(dim1:4,k,l,m,CNElemID) = PartSourceGlob(dim1:4,k,l,m,CNElemID) + Source(dim1:4)
!  ELSE
    ASSOCIATE( ShapeID => SendElemShapeID(CNElemID))
      IF(ShapeID.EQ.-1)THEN
        IPWRITE(UNIT_StdOut,*) "CNElemID   =", CNElemID
        IPWRITE(UNIT_StdOut,*) "globElemID =", globElemID
        CALL abort(__STAMP__,'SendElemShapeID(CNElemID)=-1 and therefore not correctly mapped. Increase Particles-HaloEpsVelo!')
      END IF
!      IF (CNElemID.GT.nComputeNodeElems) THEN
!        ExRankID = CNRankToSendRank(0)
!      ELSE        
!      ExRankID = CNRankToSendRank(ElemInfo_Shared(ELEM_RANK,globElemID)-ComputeNodeRootRank)
      ExRankID = CNRankToSendRank(ElemInfo_Shared(ELEM_RANK,globElemID))
!      END IF
  !    PartSourceProc(dim1:4,k,l,m,ShapeID) = PartSourceProc(dim1:4,k,l,m,ShapeID) + Source(dim1:4)
      ShapeMapping(ExRankID)%SendBuffer(dim1:4,k,l,m,ShapeID) = ShapeMapping(ExRankID)%SendBuffer(dim1:4,k,l,m,ShapeID) + Source(dim1:4)
    END ASSOCIATE
!  END IF
!#if !((USE_HDG) && (PP_nVar==1))
!#endif
END IF
#endif /*USE_MPI*/

END SUBROUTINE UpdatePartSource


PPURE LOGICAL FUNCTION ElemOnMyProc(globElemID)
!============================================================================================================================
! Check if the global element ID is an element on my processor (i.e. local element index is between 1 and nElems)
!============================================================================================================================
USE MOD_Mesh_Vars ,ONLY: nElems
USE MOD_Mesh_Vars ,ONLY: offsetElem
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: globElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: localElem
!===================================================================================================================================
localElem = globElemID-offsetElem

IF ((localElem.GE.1).AND.(localElem.LE.nElems)) THEN
  ElemOnMyProc = .TRUE.
ELSE
  ElemOnMyProc = .FALSE.
END IF
END FUNCTION ElemOnMyProc


PPURE REAL FUNCTION SFNorm(v1)
!============================================================================================================================
! Return the shape function norm by calculating the corresponding 1D, 2D or 3D distance from the input vector 'v1'
!============================================================================================================================
USE MOD_Globals      ,ONLY: VECNORM
USE MOD_PICDepo_Vars ,ONLY: dim_sf,dim_sf_dir,dim_sf_dir1,dim_sf_dir2
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: v1(1:3) !< Input vector for which the norm is calculated
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Check in which dimension the norm is to be calculated
SELECT CASE (dim_sf)
CASE (1)
  SFNorm = ABS(v1(dim_sf_dir))
CASE (2)
  SFNorm = SQRT(v1(dim_sf_dir1)**2+v1(dim_sf_dir2)**2)
CASE (3)
  SFNorm = VECNORM(v1)
CASE DEFAULT
  SFNorm = 0.
END SELECT

END FUNCTION SFNorm


PPURE REAL FUNCTION SFRadius2(v1)
!============================================================================================================================
! Return the squared distance by calculating the corresponding 1D, 2D or 3D value from the input vector 'v1'
!============================================================================================================================
USE MOD_Globals      ,ONLY: VECNORM
USE MOD_PICDepo_Vars ,ONLY: dim_sf,dim_sf_dir,dim_sf_dir1,dim_sf_dir2
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: v1(1:3) !< Input vector for which the norm is calculated
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! Check in which dimension the norm is to be calculated
SELECT CASE (dim_sf)
CASE (1)
  SFRadius2 = v1(dim_sf_dir)**2
CASE (2)
  SFRadius2 = v1(dim_sf_dir1)**2 + v1(dim_sf_dir2)**2
CASE (3)
  SFRadius2 = SUM(v1(1:3)**2)
CASE DEFAULT
  SFRadius2 = 0.
END SELECT

END FUNCTION SFRadius2


PPURE FUNCTION GetPartPosShifted(iCase,PartPos)
!============================================================================================================================
! Return the shape function norm by calculating the corresponding 1D, 2D or 3D distance from the input vector 'vec'
!============================================================================================================================
USE MOD_Globals            ,ONLY: VECNORM
USE MOD_PICDepo_Vars       ,ONLY: dim_sf,dim_sf_dir1,dim_sf_dir2,dim_sf_dir,dim_periodic_vec1,dim_periodic_vec2
USE MOD_Particle_Mesh_Vars ,ONLY: PeriodicSFCaseMatrix
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: iCase        !< Periodic case
REAL, INTENT(IN)    :: PartPos(1:3) !< Particle position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL    :: GetPartPosShifted(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: I
!===================================================================================================================================
! Select dimension of shape function
SELECT CASE (dim_sf)
CASE (1)! 1D shape function

  ! NbrOfPeriodicSFCases is either 1 (no periodic vectors: 0) or 3=2+1 (one periodic vector: -1 0 +1) cases
  GetPartPosShifted = (/0., 0., 0./)
  GetPartPosShifted(dim_sf_dir) = PartPos(dim_sf_dir) + PeriodicSFCaseMatrix(iCase,1)*GEO%PeriodicVectors(dim_sf_dir,dim_sf_dir)

CASE (2)! 2D shape function

  ! NbrOfPeriodicSFCases is either 1 (no periodic vectors: 0) or 9=8+1 (two periodic vectors) cases
  ! Constant deposition direction
  GetPartPosShifted(dim_sf_dir) = PartPos(dim_sf_dir)

  ! 1st periodic vector
  ! 1st deposition direction
  GetPartPosShifted(dim_sf_dir1) = PartPos(dim_sf_dir1)&
                                 + PeriodicSFCaseMatrix(iCase,1)*GEO%PeriodicVectors(dim_sf_dir1,dim_periodic_vec1)
  ! 2nd deposition direction
  GetPartPosShifted(dim_sf_dir2) = PartPos(dim_sf_dir2)&
                                 + PeriodicSFCaseMatrix(iCase,1)*GEO%PeriodicVectors(dim_sf_dir2,dim_periodic_vec1)

  ! 2nd periodic vector (if available)
  IF(dim_periodic_vec2.GT.0)THEN
    ! 1st deposition direction
    GetPartPosShifted(dim_sf_dir1) = GetPartPosShifted(dim_sf_dir1)&
                                   + PeriodicSFCaseMatrix(iCase,2)*GEO%PeriodicVectors(dim_sf_dir1,dim_periodic_vec2)
    ! 2nd deposition direction
    GetPartPosShifted(dim_sf_dir2) = GetPartPosShifted(dim_sf_dir2)&
                                   + PeriodicSFCaseMatrix(iCase,2)*GEO%PeriodicVectors(dim_sf_dir2,dim_periodic_vec2)
  END IF ! dim_periodic_vec2.GT.0

CASE DEFAULT!CASE(3) Standard 3D shape function

  ! NbrOfPeriodicSFCases is either 1 (no periodic vectors: 0) or 27=26+1 (three periodic vectors) cases
  DO I = 1,3
    GetPartPosShifted(I) = PartPos(I)&
                         + PeriodicSFCaseMatrix(iCase,1)*GEO%PeriodicVectors(I,1)&
                         + PeriodicSFCaseMatrix(iCase,2)*GEO%PeriodicVectors(I,2)&
                         + PeriodicSFCaseMatrix(iCase,3)*GEO%PeriodicVectors(I,3)
  END DO
END SELECT
END FUNCTION GetPartPosShifted


!===================================================================================================================================
!> Set dimension (1D, 2D or 3D) of shape function and calculate the corresponding line, area of volume to which the charge is to be
!> deposited. Output the average number of DOF that are captured by the shape function deposition kernel
!===================================================================================================================================
SUBROUTINE InitShapeFunctionDimensionalty()
! MODULES
USE MOD_Preproc
USE MOD_Globals            ,ONLY: UNIT_stdOut,MPIRoot,abort
USE MOD_PICDepo_Vars       ,ONLY: dim_sf,BetaFac,w_sf,r_sf,r2_sf,r2_sf_inv,alpha_sf,dim_sf_dir,sfDepo3D,dim_sf_dir1,dim_sf_dir2
USE MOD_PICDepo_Vars       ,ONLY: DepositionType,dimFactorSF
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,MeshVolume
USE MOD_ReadInTools        ,ONLY: PrintOption
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Mesh_Vars          ,ONLY: nGlobalElems
USE MOD_PICDepo_Tools      ,ONLY: beta
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)             :: hilf_geo
CHARACTER(1)              :: hilf_dim
INTEGER                   :: nTotalDOF
REAL                      :: VolumeShapeFunction
REAL                      :: r_sf_loc,r2_sf_loc  !> temporary variables
!===================================================================================================================================
! 0. Set global radius squared and the inverse of that
IF(TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
  r_sf_loc  = 1.
  r2_sf_loc = 1.
ELSE
  r_sf_loc  = r_sf
  r2_sf     = r_sf * r_sf  ! Radius squared
  r2_sf_loc = r2_sf
  r2_sf_inv = 1./r2_sf ! Inverse of radius squared
END IF

! 1. Initialize auxiliary variables
hilf_geo='volume'
WRITE(UNIT=hilf_dim,FMT='(I0)') dim_sf

! 2. Set the scaling factor for the shape function depending on 1D, 2D or 3D shape function and how the charge is to be distributed
SELECT CASE (dim_sf)

  CASE (1) ! --- 1D shape function -------------------------------------------------------------------------------------------------
    ! Set perpendicular directions
    IF(dim_sf_dir.EQ.1)THEN ! Shape function deposits charge in x-direction
      dimFactorSF = (GEO%ymaxglob-GEO%yminglob)*(GEO%zmaxglob-GEO%zminglob)
    ELSE IF (dim_sf_dir.EQ.2)THEN ! Shape function deposits charge in y-direction
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)*(GEO%zmaxglob-GEO%zminglob)
    ELSE IF (dim_sf_dir.EQ.3)THEN ! Shape function deposits charge in z-direction
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)*(GEO%ymaxglob-GEO%yminglob)
    END IF

    ! Set prefix factor
    IF(sfDepo3D)THEN ! Distribute the charge over the volume (3D)
      w_sf = GAMMA(REAL(alpha_sf)+1.5)/(SQRT(PI)*r_sf_loc*GAMMA(REAL(alpha_sf+1))*dimFactorSF)
    ELSE ! Distribute the charge over the line (1D)
      w_sf = GAMMA(REAL(alpha_sf)+1.5)/(SQRT(PI)*r_sf_loc*GAMMA(REAL(alpha_sf+1)))
      hilf_geo='line'
    END IF

    ! Set shape function length (3D volume)
    VolumeShapeFunction=2*r_sf_loc*dimFactorSF
    ! Calculate number of 1D DOF (assume second and third direction with 1 cell layer and area given by dimFactorSF)
    nTotalDOF=nGlobalElems*(PP_N+1)

  CASE (2) ! --- 2D shape function -------------------------------------------------------------------------------------------------
    ! Set perpendicular direction
    IF(dim_sf_dir.EQ.1)THEN ! Shape function deposits charge in y-z-direction (const. in x)
      dimFactorSF = (GEO%xmaxglob-GEO%xminglob)
    ELSE IF (dim_sf_dir.EQ.2)THEN ! Shape function deposits charge in x-z-direction (const. in y)
      dimFactorSF = (GEO%ymaxglob-GEO%yminglob)
    ELSE IF (dim_sf_dir.EQ.3)THEN! Shape function deposits charge in x-y-direction (const. in z)
      dimFactorSF = (GEO%zmaxglob-GEO%zminglob)
    END IF

    ! Set prefix factor
    IF(sfDepo3D)THEN ! Distribute the charge over the volume (3D)
      w_sf = (REAL(alpha_sf)+1.0)/(PI*r2_sf_loc*dimFactorSF)
    ELSE ! Distribute the charge over the area (2D)
      w_sf = (REAL(alpha_sf)+1.0)/(PI*r2_sf_loc)
      hilf_geo='area'
    END IF

    ! set the two perpendicular directions used for deposition
    dim_sf_dir1 = MERGE(1,2,dim_sf_dir.EQ.2)
    dim_sf_dir2 = MERGE(1,MERGE(3,3,dim_sf_dir.EQ.2),dim_sf_dir.EQ.3)
    LBWRITE(UNIT_stdOut,'(A,I0,A,I0,A,I0,A)') ' Shape function 2D with const. distribution in dir ',dim_sf_dir,&
        ' and variable distrbution in ',dim_sf_dir1,' and ',dim_sf_dir2,' (1: x, 2: y and 3: z)'

    ! Set shape function length (3D volume)
    VolumeShapeFunction=PI*(r_sf_loc**2)*dimFactorSF
    ! Calculate number of 2D DOF (assume third direction with 1 cell layer and width dimFactorSF)
    nTotalDOF=nGlobalElems*(PP_N+1)**2

  CASE (3) ! --- 3D shape function -------------------------------------------------------------------------------------------------
    ! Set prefix factor (not for shape_function_adaptive)
    BetaFac = beta(1.5, REAL(alpha_sf) + 1.)
    w_sf = 1./(2. * BetaFac * REAL(alpha_sf) + 2 * BetaFac) * (REAL(alpha_sf) + 1.)/(PI*(r_sf_loc**3))

  CASE DEFAULT
    CALL abort(__STAMP__,'Shape function dimension must be 1, 2 or 3')
END SELECT

LBWRITE(UNIT_stdOut,'(A)') ' The complete charge is '//TRIM(hilf_geo)//' distributed (via '//TRIM(hilf_dim)//'D shape function)'

IF(.NOT.sfDepo3D)THEN
  LBWRITE(UNIT_stdOut,'(A)') ' Note that the integral of the charge density over the mesh volume is larger than the complete charge'
  LBWRITE(UNIT_stdOut,'(A)') ' because the charge is spread out over either a line (1D shape function) or an area (2D shape function)!'
END IF

! 3. Output info regarding charge distribution and points per shape function resolution
IF(.NOT.TRIM(DepositionType).EQ.'shape_function_adaptive')THEN
  ASSOCIATE(nTotalDOFin3D             => nGlobalElems*(PP_N+1)**3 ,&
            VolumeShapeFunctionSphere => 4./3.*PI*r_sf**3         )

    ! Output shape function volume
    IF(dim_sf.EQ.1)THEN
      CALL PrintOption('Shape function volume (corresponding to a cuboid in 3D)'  , 'CALCUL.', RealOpt=VolumeShapeFunction)
    ELSEIF(dim_sf.EQ.2)THEN
      CALL PrintOption('Shape function volume (corresponding to a cylinder in 3D)', 'CALCUL.', RealOpt=VolumeShapeFunction)
    ELSE
      VolumeShapeFunction = VolumeShapeFunctionSphere
      CALL PrintOption('Shape function volume (corresponding to a sphere in 3D)'  , 'CALCUL.', RealOpt=VolumeShapeFunction)
    END IF

    ! Sanity check: Shape function volume is not allowed to be larger than the complete mesh simulation domain
    IF(MPIRoot)THEN
      IF(VolumeShapeFunction.GT.MeshVolume)THEN
        CALL PrintOption('Mesh volume ('//TRIM(hilf_dim)//')', 'CALCUL.' , RealOpt=MeshVolume)
        WRITE(UNIT_stdOut,'(A)') ' Maybe wrong perpendicular direction (PIC-shapefunction-direction)?'
        CALL abort(__STAMP__,'ShapeFunctionVolume > MeshVolume ('//TRIM(hilf_dim)//' shape function)')
      END IF
    END IF

    ! Display 1D or 2D deposition info
    IF(dim_sf.NE.3)THEN
      CALL PrintOption('Average DOFs in Shape-Function '//TRIM(hilf_geo)//' ('//TRIM(hilf_dim)//')' , 'CALCUL.' , RealOpt=&
          REAL(nTotalDOF)*VolumeShapeFunction/MeshVolume)
    END IF ! dim_sf.NE.3

    CALL PrintOption('Average DOFs in Shape-Function (corresponding 3D sphere)' , 'CALCUL.' , RealOpt=&
        REAL(nTotalDOFin3D)*VolumeShapeFunctionSphere/MeshVolume)
  END ASSOCIATE
END IF ! .NOT.TRIM(DepositionType).EQ.'shape_function_adaptive'

END SUBROUTINE InitShapeFunctionDimensionalty


END MODULE MOD_PICDepo_Shapefunction_Tools
