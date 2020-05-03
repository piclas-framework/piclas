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

INTERFACE DepoSFParticleLocally
  MODULE PROCEDURE DepoSFParticleLocally
END INTERFACE

PUBLIC:: calcSfSource,DepoSFParticleLocally
!===================================================================================================================================

CONTAINS

SUBROUTINE calcSfSource(SourceSize_in,ChargeMPF,Vec1,Vec2,Vec3,PartPos,PartIdx,PartVelo,const_opt)
!============================================================================================================================
! deposit charges on DOFs via shapefunction including periodic displacements and mirroring with SFdepoFixes
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars,           ONLY:r_sf,DepositionType
USE MOD_PICDepo_Vars,           ONLY:NbrOfSFdepoFixes,SFdepoFixesGeo,SFdepoFixesBounds,SFdepoFixesChargeMult
USE MOD_PICDepo_Vars,           ONLY:SFdepoFixesPartOfLink,SFdepoFixesEps,NbrOfSFdepoFixLinks,SFdepoFixLinks
USE MOD_Globals
USE MOD_Particle_Mesh_Vars,     ONLY:casematrix,NbrOfCases
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: SourceSize_in,PartIdx
REAL, INTENT(IN)                 :: ChargeMPF,PartPos(3),Vec1(3),Vec2(3),Vec3(3)
REAL, INTENT(IN), OPTIONAL       :: PartVelo(3)
LOGICAL, INTENT(IN), OPTIONAL    :: const_opt
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
INTEGER                          :: iCase, ind
REAL                             :: ShiftedPart(1:3), caseShiftedPart(1:3), n_loc(1:3)
INTEGER                          :: iSFfix, LinkLoopEnd(2), iSFfixLink, iTwin, iLinkRecursive, SFfixIdx, SFfixIdx2
LOGICAL                          :: DoCycle, DoNotDeposit
REAL                             :: SFfixDistance, SFfixDistance2
LOGICAL , ALLOCATABLE            :: SFdepoFixDone(:)
LOGICAL                          :: const
!----------------------------------------------------------------------------------------------------------------------------------
IF (PRESENT(const_opt)) THEN
  const=const_opt
ELSE
  const=.FALSE.
END IF
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
IF (NbrOfSFdepoFixes.EQ.0) THEN
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      ShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    Fac = Fac2
    SELECT CASE(TRIM(DepositionType))
    CASE('shape_function')
      CALL depoChargeOnDOFs_sf(ShiftedPart,SourceSize,Fac,const)
    CASE('shape_function_simple')
      CALL depoChargeOnDOFs_sf_simple(ShiftedPart,SourceSize,Fac,const)
    END SELECT
  END DO ! iCase (periodicity)
ELSE ! NbrOfSFdepoFixes.NE.0
  ALLOCATE(SFdepoFixDone(0:NbrOfSFdepoFixes))
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      caseShiftedPart(ind) = PartPos(ind) + casematrix(iCase,1)*Vec1(ind) + &
        casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    SFdepoFixDone=.FALSE.
    DO iSFfix=0,NbrOfSFdepoFixes
      IF (SFdepoFixesPartOfLink(iSFfix)) CYCLE !this SFfix will be already covered as part of a FixLink
      IF (iSFfix.EQ.0) THEN
        LinkLoopEnd(1)=NbrOfSFdepoFixLinks !non-mirrored position: consider FixLinks
      ELSE
        LinkLoopEnd(1)=0
      END IF
      DO iSFfixLink=0,LinkLoopEnd(1)
        DoCycle=.FALSE.
        !-- strategy for SFdepoFixLink:
        !-- 1.: create 2 identities by iTwin for iLinkRecursive=0 (non-mirror and mirror at SFdepoFixLinks(*,1))
        !-- 2.: mirror recursevely further for iLinkRecursive>0 (alternating at SFdepoFixLinks(*,2) and SFdepoFixLinks(*,1))
        DO iTwin=1,2
          IF (iSFfixLink.EQ.0 .AND. iTwin.EQ.2) EXIT !no SFdepoFixLink
          Fac = Fac2
          ShiftedPart=caseShiftedPart
          IF (iSFfixLink.EQ.0) THEN
            LinkLoopEnd(2)=0 !no SFdepoFixLink
          ELSE
            LinkLoopEnd(2)=ABS(SFdepoFixLinks(iSFfixLink,3))-1 !might be negative as flag for special case
          END IF
          DO iLinkRecursive=0,LinkLoopEnd(2)
            IF (iLinkRecursive.EQ.0) THEN !(first!)
              IF (iTwin.EQ.1) THEN
                SFfixIdx = iSFfix
              ELSE
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            ELSE
              IF (MOD(iLinkRecursive,2).NE.0) THEN !uneven (second and later: 1, 3, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,2)
              ELSE !even (third and later: 2, 4, ...)
                SFfixIdx =SFdepoFixLinks(iSFfixLink,1)
              END IF
            END IF
            DoNotDeposit=.FALSE.
            IF (iLinkRecursive.EQ.0 .OR. (iLinkRecursive.EQ.1 .AND. iTwin.EQ.1)) THEN !single or one-time-mirrored identity
              IF (SFdepoFixDone(SFfixIdx)) DoNotDeposit=.TRUE. !do not deposit this charge (but save position for further mirroring)
              SFdepoFixDone(SFfixIdx) = .TRUE.
            ELSE IF (iSFfixLink.GT.0) THEN
              IF (SFdepoFixLinks(iSFfixLink,3).LT.0) EXIT !skip double mirroring for 120 deg case
            END IF
            IF (SFfixIdx.GT.0) THEN
              IF ( SFdepoFixesBounds(SFfixIdx,1,1).GT.ShiftedPart(1) .OR. ShiftedPart(1).GT.SFdepoFixesBounds(SFfixIdx,2,1) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,2).GT.ShiftedPart(2) .OR. ShiftedPart(2).GT.SFdepoFixesBounds(SFfixIdx,2,2) .OR. &
                SFdepoFixesBounds(SFfixIdx,1,3).GT.ShiftedPart(3) .OR. ShiftedPart(3).GT.SFdepoFixesBounds(SFfixIdx,2,3) ) THEN
                DoCycle=.TRUE.
                EXIT !do not shift this particle (-> go to next single SFfixIdx -> iSFfixLink-loop)
              END IF
              SFfixDistance = SFdepoFixesGeo(SFfixIdx,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx,1,1)) &
                + SFdepoFixesGeo(SFfixIdx,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx,1,2)) &
                + SFdepoFixesGeo(SFfixIdx,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx,1,3))
              SFfixIdx2=0 !init
              IF (SFfixDistance .GT. SFdepoFixesEps) THEN
                IPWRITE(*,'(I4,A,3(x,E12.5))')' original case-pos:',caseShiftedPart
                IPWRITE(*,'(I4,A,3(x,E12.5))')' current pos:',ShiftedPart
                IPWRITE(*,'(I4,7(A,I0))') &
                  ' iCase: ',iCase,', iSFfix:',iSFfix,', iSFfixLink:',iSFfixLink,', iTwin:',iTwin,', iLinkRec:',iLinkRecursive,&
                  ', SFfixIdx:',SFfixIdx,', PartIdx:',PartIdx
                CALL abort(&
                  __STAMP__ &
                  ,'Particle is outside of SF-Fix-Plane! (For Layer-/Resample-Parts: try -UseFixBounds)',SFfixIdx,SFfixDistance)
              ELSE IF ( (SFfixDistance.LT.-r_sf) ) THEN !.OR. (SFfixDistance.GT.0.) ) THEN !SFfixDistance>0 are particle within eps
                DoNotDeposit=.TRUE. !too far inside so that mirrored SF would not reach any DOF
              ELSE IF (iSFfixLink.GT.0) THEN !check which is the other SFfixIdx of current link
                IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,1)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,2)
                ELSE IF (SFfixIdx.EQ.SFdepoFixLinks(iSFfixLink,2)) THEN
                  SFfixIdx2=SFdepoFixLinks(iSFfixLink,1)
                ELSE
                  CALL abort(&
                    __STAMP__ &
                    ,'Something is wrong with iSFfixLink',iSFfixLink)
                END IF
              END IF
              ShiftedPart(1:3) = ShiftedPart(1:3) - 2.*SFfixDistance*SFdepoFixesGeo(SFfixIdx,2,1:3)
              Fac = Fac * SFdepoFixesChargeMult(SFfixIdx)
!#if !((USE_HDG) && (PP_nVar==1))
              IF (SourceSize.EQ.4) THEN
                ! change velocity
                n_loc = SFdepoFixesGeo(SFfixIdx,2,1:3)
                Fac(1:3) = Fac2(1:3) -2.*DOT_PRODUCT(Fac2(1:3),n_loc)*n_loc
              END IF
!#endif
              IF (SFfixIdx2.NE.0) THEN !check if new position would not reach a dof because of the other plane
                SFfixDistance2 = SFdepoFixesGeo(SFfixIdx2,2,1)*(ShiftedPart(1)-SFdepoFixesGeo(SFfixIdx2,1,1)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,2)*(ShiftedPart(2)-SFdepoFixesGeo(SFfixIdx2,1,2)) &
                  + SFdepoFixesGeo(SFfixIdx2,2,3)*(ShiftedPart(3)-SFdepoFixesGeo(SFfixIdx2,1,3))
                IF (SFfixDistance2 .GT. r_sf) DoNotDeposit=.TRUE. !too far outside of plane
              END IF
            END IF
            IF (DoNotDeposit) CYCLE !(-> do not deposit but save position for possible further recursive mirroring)
            !------------- actual deposition:
            SELECT CASE(TRIM(DepositionType))
            CASE('shape_function')
              CALL depoChargeOnDOFs_sf(ShiftedPart,SourceSize,Fac,const)
            CASE('shape_function_simple')
              CALL depoChargeOnDOFs_sf_simple(ShiftedPart,SourceSize,Fac,const)
            END SELECT
          END DO ! iLinkRecursive
          IF (DoCycle) EXIT
        END DO ! iTwin
      END DO ! iSFfixLink
    END DO ! iSFfix
  END DO ! iCase (periodicity)
END IF !NbrOfSFdepoFixes

END SUBROUTINE calcSfSource


SUBROUTINE depoChargeOnDOFs_sf(Position,SourceSize,Fac,const)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars,           ONLY:PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf, ElemDepo_xGP, PartSourceConst
USE MOD_Mesh_Vars,              ONLY:nElems
USE MOD_Particle_Mesh_Vars,     ONLY:GEO
USE MOD_PreProc,                ONLY:PP_N
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
LOGICAL, INTENT(IN)              :: const
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
LOGICAL                          :: chargedone(1:nElems)
INTEGER                          :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER                          :: kk, ll, mm, ppp
INTEGER                          :: ElemID
REAL                             :: radius2, S, S1
REAL                             :: dx,dy,dz
INTEGER                          :: expo
!----------------------------------------------------------------------------------------------------------------------------------

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
      DO ppp = 1,GEO%FIBGM(kk,ll,mm)%nElem
        ElemID = GEO%FIBGM(kk,ll,mm)%Element(ppp)
        IF(ElemID.GT.nElems) CYCLE
        IF (.NOT.chargedone(ElemID)) THEN
#if USE_LOADBALANCE
          nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*USE_LOADBALANCE*/
          !--- go through all gauss points
          DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
            !-- calculate distance between gauss and particle
            dX = ABS(Position(1) - ElemDepo_xGP(1,k,l,m,ElemID))
            IF(dX.GT.r_sf) CYCLE
            dY = ABS(Position(2) - ElemDepo_xGP(2,k,l,m,ElemID))
            IF(dY.GT.r_sf) CYCLE
            dZ = ABS(Position(3) - ElemDepo_xGP(3,k,l,m,ElemID))
            IF(dZ.GT.r_sf) CYCLE
            radius2 = dX*dX+dY*dY+dZ*dZ
            !-- calculate charge and current density at ip point using a shape function
            !-- currently only one shapefunction available, more to follow (including structure change)
            IF (radius2 .LT. r2_sf) THEN
              S = 1. - r2_sf_inv * radius2
              S1 = S*S
              DO expo = 3, alpha_sf
                S1 = S*S1
              END DO
              IF (const) THEN
                IF (SourceSize.EQ.1) THEN
                  PartSourceConst(4,k,l,m,ElemID) = PartSourceConst(4,k,l,m,ElemID) + Fac(4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
                ELSE IF (SourceSize.EQ.4) THEN
                  PartSourceConst(1:4,k,l,m,ElemID) = PartSourceConst(1:4,k,l,m,ElemID) + Fac(1:4) * S1
!#endif
                END IF
              ELSE !.NOT.const
                IF (SourceSize.EQ.1) THEN
                  PartSource(4,k,l,m,ElemID) = PartSource(4,k,l,m,ElemID) + Fac(4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
                ELSE IF (SourceSize.EQ.4) THEN
                  PartSource(1:4,k,l,m,ElemID) = PartSource(1:4,k,l,m,ElemID) + Fac(1:4) * S1
!#endif
                END IF
              END IF !const.
            END IF
          END DO; END DO; END DO
          chargedone(ElemID) = .TRUE.
        END IF
      END DO ! ppp
    END DO ! mm
  END DO ! ll
END DO ! kk

END SUBROUTINE depoChargeOnDOFs_sf


SUBROUTINE depoChargeOnDOFs_sf_simple(Position,SourceSize,Fac,const)
!============================================================================================================================
! actual deposition of single charge on DOFs via shapefunction_simple (i.e. loop through all elems instead of part-dependency: efficient for small elem-nbr!)
!============================================================================================================================
! use MODULES
USE MOD_Mesh_Vars,              ONLY:ElemBaryNGeo
USE MOD_PICDepo_Vars,           ONLY:PartSource, r_sf, r2_sf, r2_sf_inv, alpha_sf, ElemDepo_xGP, ElemRadius2_sf, PartSourceConst
USE MOD_Particle_Mesh_Vars,     ONLY:ElemRadiusNGeo
USE MOD_PreProc,                ONLY:PP_N, PP_nElems
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
LOGICAL, INTENT(IN)              :: const
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: k, l, m
INTEGER                          :: ElemID
REAL                             :: radius2, S, S1
REAL                             :: dx,dy,dz
INTEGER                          :: expo
!----------------------------------------------------------------------------------------------------------------------------------

DO ElemID=1,PP_nElems
  dX = ABS(Position(1) - ElemBaryNgeo(1,ElemID))
  IF(dX.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  dY = ABS(Position(2) - ElemBaryNgeo(2,ElemID))
  IF(dY.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  dZ = ABS(Position(3) - ElemBaryNgeo(3,ElemID))
  IF(dZ.GT.r_sf+ElemRadiusNGeo(ElemID)) CYCLE
  radius2 = dX*dX+dY*dY+dZ*dZ
  IF(radius2.GT.ElemRadius2_sf(ElemID)) CYCLE
#if USE_LOADBALANCE
  nDeposPerElem(ElemID)=nDeposPerElem(ElemID)+1
#endif /*USE_LOADBALANCE*/
  DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    !-- calculate distance between gauss and particle
    dX = ABS(Position(1) - ElemDepo_xGP(1,k,l,m,ElemID))
    IF(dX.GT.r_sf) CYCLE
    dY = ABS(Position(2) - ElemDepo_xGP(2,k,l,m,ElemID))
    IF(dY.GT.r_sf) CYCLE
    dZ = ABS(Position(3) - ElemDepo_xGP(3,k,l,m,ElemID))
    IF(dZ.GT.r_sf) CYCLE
    radius2 = dX*dX+dY*dY+dZ*dZ
    !-- calculate charge and current density at ip point using a shape function
    !-- currently only one shapefunction available, more to follow (including structure change)
    IF (radius2 .GT. r2_sf) CYCLE
    S = 1. - r2_sf_inv * radius2
    S1 = S*S
    DO expo = 3, alpha_sf
      S1 = S*S1
    END DO
    IF (const) THEN
      IF (SourceSize.EQ.1) THEN
        PartSourceConst(4,k,l,m,ElemID) = PartSourceConst(4,k,l,m,ElemID) + Fac(4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
      ELSE IF (SourceSize.EQ.4) THEN
        PartSourceConst(1:4,k,l,m,ElemID) = PartSourceConst(1:4,k,l,m,ElemID) + Fac(1:4) * S1
!#endif
      END IF
    ELSE !.NOT.const
      IF (SourceSize.EQ.1) THEN
        PartSource(4,k,l,m,ElemID) = PartSource(4,k,l,m,ElemID) + Fac(4) * S1
!#if !((USE_HDG) && (PP_nVar==1))
      ELSE IF (SourceSize.EQ.4) THEN
        PartSource(1:4,k,l,m,ElemID) = PartSource(1:4,k,l,m,ElemID) + Fac(1:4) * S1
!#endif
      END IF
    END IF !const
  END DO; END DO; END DO
END DO !ElemID=1,PP_nElems

END SUBROUTINE depoChargeOnDOFs_sf_simple


!==================================================================================================================================
!> Check whether a particle is inside of a local deposition element, where instead of the shape function, a local deposition method
!> is used.
!==================================================================================================================================
SUBROUTINE DepoSFParticleLocally(DepoLoc,ElemID,PartID)
! MODULES                                                                                                                          !
USE MOD_PreProc
USE MOD_PICDepo_Vars           ,ONLY: DoSFLocalDepoAtBounds,CellVolWeight_Volumes,cellvolweightfac,PartSource
USE MOD_Particle_Vars          ,ONLY: PEM
USE MOD_Particle_Mesh_Vars     ,ONLY: IsLocalDepositionBCElem
USE MOD_Particle_Vars          ,ONLY: PartState,PartPosRef,PartSpecies,Species,PartMPF,usevMPF
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: ElemID  !< Element ID
INTEGER,INTENT(IN)  :: PartID  !< Particle ID
LOGICAL,INTENT(OUT) :: DepoLoc !< Returns true when particle is deposited locally, else returns false
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: BGMSourceCellVol_loc(0:1,0:1,0:1,1:4)
REAL                :: Charge, TSource(1:4), alpha1, alpha2, alpha3, temppartpos(1:3)
INTEGER             :: k,l,m
!===================================================================================================================================

! =============================
! Workflow:
!
!  1.  Check if local deposition is used. If not: return
!  2.  Check Check if particle is inside of a local deposition element and if the current element is the same
!  3.  Local deposition via cell vol weight method
!==============================

! 1.  Check if local deposition is used. If not: return
IF(.NOT.DoSFLocalDepoAtBounds)THEN
  DepoLoc=.FALSE.
  RETURN
END IF

! 2.  Check Check if particle is inside of a local deposition element and if the current element is the same
IF(IsLocalDepositionBCElem(PEM%Element(PartID)))THEN ! Particle element is a local deposition element
  DepoLoc=.TRUE.
  ! Check if particle is NOT in current element to prevent deposition in neighboring elements.
  IF(ElemID.NE.PEM%Element(PartID)) RETURN
ELSE ! Particle element is NOT a local deposition element: perform deposition via shape function
  DepoLoc=.FALSE.
  RETURN
END IF

! 3.  Local deposition via cell vol weight method
BGMSourceCellVol_loc = 0.0
IF (usevMPF) THEN
  Charge= Species(PartSpecies(PartID))%ChargeIC * PartMPF(PartID)
ELSE
  Charge= Species(PartSpecies(PartID))%ChargeIC * Species(PartSpecies(PartID))%MacroParticleFactor
END IF ! usevMPF
IF(DoRefMapping)THEN
  TempPartPos(1:3)=PartPosRef(1:3,PartID)
ELSE
  CALL GetPositionInRefElem(PartState(1:3,PartID),TempPartPos,ElemID,ForceMode=.TRUE.)
END IF
TSource(:) = 0.0
!#if (PP_nVar==8)
TSource(1) = PartState(4,PartID)*Charge
TSource(2) = PartState(5,PartID)*Charge
TSource(3) = PartState(6,PartID)*Charge
!#endif
TSource(4) = Charge
alpha1=(TempPartPos(1)+1.0)/2.0
alpha2=(TempPartPos(2)+1.0)/2.0
alpha3=(TempPartPos(3)+1.0)/2.0
BGMSourceCellVol_loc(0,0,0,1:4) = BGMSourceCellVol_loc(0,0,0,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
BGMSourceCellVol_loc(0,0,1,1:4) = BGMSourceCellVol_loc(0,0,1,1:4) + (TSource(1:4)*(1-alpha1)*(1-alpha2)*(alpha3))
BGMSourceCellVol_loc(0,1,0,1:4) = BGMSourceCellVol_loc(0,1,0,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(1-alpha3))
BGMSourceCellVol_loc(0,1,1,1:4) = BGMSourceCellVol_loc(0,1,1,1:4) + (TSource(1:4)*(1-alpha1)*(alpha2)*(alpha3))
BGMSourceCellVol_loc(1,0,0,1:4) = BGMSourceCellVol_loc(1,0,0,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(1-alpha3))
BGMSourceCellVol_loc(1,0,1,1:4) = BGMSourceCellVol_loc(1,0,1,1:4) + (TSource(1:4)*(alpha1)*(1-alpha2)*(alpha3))
BGMSourceCellVol_loc(1,1,0,1:4) = BGMSourceCellVol_loc(1,1,0,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(1-alpha3))
BGMSourceCellVol_loc(1,1,1,1:4) = BGMSourceCellVol_loc(1,1,1,1:4) + (TSource(1:4)*(alpha1)*(alpha2)*(alpha3))

BGMSourceCellVol_loc(0,0,0,:) = BGMSourceCellVol_loc(0,0,0,1:4)/CellVolWeight_Volumes(0,0,0,ElemID)
BGMSourceCellVol_loc(0,0,1,:) = BGMSourceCellVol_loc(0,0,1,1:4)/CellVolWeight_Volumes(0,0,1,ElemID)
BGMSourceCellVol_loc(0,1,0,:) = BGMSourceCellVol_loc(0,1,0,1:4)/CellVolWeight_Volumes(0,1,0,ElemID)
BGMSourceCellVol_loc(0,1,1,:) = BGMSourceCellVol_loc(0,1,1,1:4)/CellVolWeight_Volumes(0,1,1,ElemID)
BGMSourceCellVol_loc(1,0,0,:) = BGMSourceCellVol_loc(1,0,0,1:4)/CellVolWeight_Volumes(1,0,0,ElemID)
BGMSourceCellVol_loc(1,0,1,:) = BGMSourceCellVol_loc(1,0,1,1:4)/CellVolWeight_Volumes(1,0,1,ElemID)
BGMSourceCellVol_loc(1,1,0,:) = BGMSourceCellVol_loc(1,1,0,1:4)/CellVolWeight_Volumes(1,1,0,ElemID)
BGMSourceCellVol_loc(1,1,1,:) = BGMSourceCellVol_loc(1,1,1,1:4)/CellVolWeight_Volumes(1,1,1,ElemID)

DO k = 0, PP_N
  DO l = 0, PP_N
    DO m = 0, PP_N
      alpha1 = CellVolWeightFac(k)
      alpha2 = CellVolWeightFac(l)
      alpha3 = CellVolWeightFac(m)
      PartSource(1:4,k,l,m,ElemID) =PartSource(1:4,k,l,m,ElemID)     + &
          BGMSourceCellVol_loc(0,0,0,1:4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
          BGMSourceCellVol_loc(0,0,1,1:4) * (1-alpha1) * (1-alpha2) *   (alpha3) + &
          BGMSourceCellVol_loc(0,1,0,1:4) * (1-alpha1) *   (alpha2) * (1-alpha3) + &
          BGMSourceCellVol_loc(0,1,1,1:4) * (1-alpha1) *   (alpha2) *   (alpha3) + &
          BGMSourceCellVol_loc(1,0,0,1:4) *   (alpha1) * (1-alpha2) * (1-alpha3) + &
          BGMSourceCellVol_loc(1,0,1,1:4) *   (alpha1) * (1-alpha2) *   (alpha3) + &
          BGMSourceCellVol_loc(1,1,0,1:4) *   (alpha1) *   (alpha2) * (1-alpha3) + &
          BGMSourceCellVol_loc(1,1,1,1:4) *   (alpha1) *   (alpha2) *   (alpha3)
    END DO !m
  END DO !l
END DO !k

END SUBROUTINE DepoSFParticleLocally


END MODULE MOD_PICDepo_Shapefunction_Tools
