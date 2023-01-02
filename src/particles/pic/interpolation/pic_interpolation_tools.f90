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

MODULE  MOD_PICInterpolation_tools
!===================================================================================================================================
!
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE GetExternalFieldAtParticle
  MODULE PROCEDURE GetExternalFieldAtParticle
END INTERFACE

INTERFACE GetInterpolatedFieldPartPos
  MODULE PROCEDURE GetInterpolatedFieldPartPos
END INTERFACE

INTERFACE GetEMField
  MODULE PROCEDURE GetEMField
END INTERFACE

INTERFACE InterpolateVariableExternalField1D
  MODULE PROCEDURE InterpolateVariableExternalField1D
END INTERFACE

PUBLIC :: GetExternalFieldAtParticle
PUBLIC :: GetInterpolatedFieldPartPos
PUBLIC :: GetEMField
PUBLIC :: InterpolateVariableExternalField1D
!===================================================================================================================================

CONTAINS


PPURE FUNCTION GetExternalFieldAtParticle(pos)
!===================================================================================================================================
! Get the external field (analytic, variable, etc.) for the particle at position pos
! 4 Methods can be used:
!   0. External field from analytic function (only for convergence tests and compiled with CODE_ANALYZE=ON)
!   1. External E field from user-supplied vector (const.) and
!      B field from CSV file (only Bz) that is interpolated to the particle z-coordinate
!   2. External field from CSV file (only Bz) that is interpolated to the particle z-coordinate
!   3. External E and B field from user-supplied vector (const.)
!===================================================================================================================================
! MODULES
USE MOD_PICInterpolation_Vars ,ONLY: externalField,useVariableExternalField,useAlgebraicExternalField,VariableExternalFieldDim
#ifdef CODE_ANALYZE
USE MOD_PICInterpolation_Vars ,ONLY: DoInterpolationAnalytic
#endif /*CODE_ANALYZE*/
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN) :: pos(3) ! position x,y,z
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GetExternalFieldAtParticle(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetExternalFieldAtParticle=0.
#ifdef CODE_ANALYZE
! 0. External field from analytic function (only for convergence tests)
IF(DoInterpolationAnalytic)THEN ! use analytic/algebraic functions for the field interpolation
  GetExternalFieldAtParticle(1:6) = GetAnalyticFieldAtParticle(pos)
ELSE ! use variable or fixed external field
#endif /*CODE_ANALYZE*/
!#if (PP_nVar==8))
  IF(useVariableExternalField)THEN
    ! 1. External E field from user-supplied vector (const.) and
    GetExternalFieldAtParticle(1:6) = externalField(1:6)
    ! Select spatial dimension for interpolation
    SELECT CASE(VariableExternalFieldDim)
    CASE(1)
      ! B field from CSV file (only Bz) that is interpolated to the particle z-coordinate (linear)
      GetExternalFieldAtParticle(6) = InterpolateVariableExternalField1D(pos(3))
    CASE(2)
      ! B field from .h5 file (only Bz and Br) that is interpolated to the particle z- and r-coordinate (bilinear)
      GetExternalFieldAtParticle(4:6) = InterpolateVariableExternalField2D(pos(1:3))
    CASE(3)
      ! B field from .h5 file (Bx, By and Bz) that is interpolated to the particle x-, y- and z-coordinate (trilinear)
      GetExternalFieldAtParticle(4:6) = InterpolateVariableExternalField3D(pos(1:3))
    END SELECT
  ELSEIF(useAlgebraicExternalField)THEN
    ! 2. External E and B field from algebraic expression that is interpolated to the particle position
    GetExternalFieldAtParticle(1:6) = InterpolateAlgebraicExternalField(pos)
  ELSE
    ! 3. External E and B field from user-supplied vector (const.)
    GetExternalFieldAtParticle(1:6) = externalField(1:6)
  END IF
!#endif /*(PP_nVar==8))*/
#ifdef CODE_ANALYZE
END IF
#endif /*CODE_ANALYZE*/

END FUNCTION GetExternalFieldAtParticle


#ifdef CODE_ANALYZE
PPURE FUNCTION GetAnalyticFieldAtParticle(PartPos)
!===================================================================================================================================
! Calculate the electro-(magnetic) field at the particle's position form an analytic solution
!===================================================================================================================================
! MODULES
USE MOD_PICInterpolation_Vars ,ONLY: AnalyticInterpolationType,AnalyticInterpolationSubType
USE MOD_PICInterpolation_Vars ,ONLY: AnalyticInterpolationPhase,AnalyticInterpolationGamma,AnalyticInterpolationE
USE MOD_Globals_Vars          ,ONLY: c
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: PartPos(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GetAnalyticFieldAtParticle(1:6)
REAL :: v_perp ! Perpendicular velocity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
GetAnalyticFieldAtParticle(1:6) = 0.
SELECT CASE(AnalyticInterpolationType)
CASE(0) ! 0: const. magnetostatic field: B = B_z = (/ 0 , 0 , 1 T /) = const.

  ! 0: non-relativistic, 1: relativistic
  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0) ! 0: non-relativistic
    GetAnalyticFieldAtParticle(6) = 1.0
  CASE(1) ! 1: relativistic
    ASSOCIATE( gamma1 => AnalyticInterpolationGamma ,& ! Lorentz factor
               m      => 1.0                        ,& ! [kg] particle mass
               q      => 1.0                        ,& ! [C] particle charge
               phi    => AnalyticInterpolationPhase )  ! [rad] phase shift
      !-- get Lorentz factor gamma1(n)
      v_perp = c*SQRT(1.0 - 1/(gamma1**2))
      !-- Set const. magnetic field [T]
      GetAnalyticFieldAtParticle(6) = gamma1*v_perp
    END ASSOCIATE
  END SELECT

CASE(1) ! magnetostatic field: B = B_z = B_0 * EXP(x/l)

  ASSOCIATE( B_0 => 1.0, l => 1.0  )
    GetAnalyticFieldAtParticle(6) = B_0 * EXP(PartPos(1) / l)
  END ASSOCIATE

CASE(2)

  ! const. electromagnetic field: B = B_z = (/ 0 , 0 , (x^2+y^2)^0.5 /) = const.
  !                                  E = 1e-2/(x^2+y^2)^(3/2) * (/ x , y , 0. /)
  ! Example from Paper by H. Qin: Why is Boris algorithm so good? (2013)
  ! http://dx.doi.org/10.1063/1.4818428
  ASSOCIATE( x => PartPos(1), y => PartPos(2) )
    ! Ex and Ey
    GetAnalyticFieldAtParticle(1) = 1.0e-2 * (x**2+y**2)**(-1.5) * x
    GetAnalyticFieldAtParticle(2) = 1.0e-2 * (x**2+y**2)**(-1.5) * y
    ! Bz
    GetAnalyticFieldAtParticle(6) = SQRT(x**2+y**2)
  END ASSOCIATE

CASE(3) ! 3: const. electric field: E = E_x = (/ 1 V/m , 0 , 0 /) = const.

  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0,1) ! 0: non-relativistic, 1: relativistic
    GetAnalyticFieldAtParticle(1) = AnalyticInterpolationE
  END SELECT

CASE(4) ! 4: const. electric field: E = E_x = (/ x V/m , 0 , 0 /) = const.

  SELECT CASE(AnalyticInterpolationSubType)
  CASE(0,1) ! 0: non-relativistic, 1: relativistic
    GetAnalyticFieldAtParticle(1) = AnalyticInterpolationE
  END SELECT

END SELECT
END FUNCTION GetAnalyticFieldAtParticle
#endif /*CODE_ANALYZE*/


FUNCTION GetInterpolatedFieldPartPos(ElemID,PartID)
!===================================================================================================================================
! Evaluate the electro-(magnetic) field using the reference position and return the field
!===================================================================================================================================
! MODULES
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: PartPosRef,PDM,PartState,PEM
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
USE MOD_Particle_Vars          ,ONLY: DoSurfaceFlux
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID !< Global element ID
INTEGER,INTENT(IN) :: PartID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GetInterpolatedFieldPartPos(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: PartPosRef_loc(1:3)
LOGICAL                      :: SucRefPos
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
LOGICAL                      :: NotMappedSurfFluxParts
#else
LOGICAL,PARAMETER            :: NotMappedSurfFluxParts=.FALSE.
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/
!===================================================================================================================================

! Check Surface Flux Particles
#if (PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)
NotMappedSurfFluxParts=DoSurfaceFlux !Surfaceflux particles inserted before interpolation and tracking. Field at wall is needed!
#endif /*(PP_TimeDiscMethod>=500) && (PP_TimeDiscMethod<=509)*/

SucRefPos = .TRUE. ! Initialize for all methods

! Check if reference position is required
IF(NotMappedSurfFluxParts .AND.(TrackingMethod.EQ.REFMAPPING))THEN
  IF(PDM%dtFracPush(PartID)) CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef_loc(1:3),ElemID)
ELSEIF(TrackingMethod.NE.REFMAPPING)THEN
  CALL GetPositionInRefElem(PartState(1:3,PartID),PartPosRef_loc(1:3),ElemID, isSuccessful = SucRefPos)
ELSE
  PartPosRef_loc(1:3) = PartPosRef(1:3,PartID)
END IF

! Interpolate the field and return the vector
IF ((.NOT.SucRefPos).AND.(TRIM(DepositionType).EQ.'cell_volweight_mean')) THEN
  GetInterpolatedFieldPartPos(1:6) =  GetEMFieldDW(PEM%LocalElemID(PartID),PartState(1:3,PartID))
ELSE
  GetInterpolatedFieldPartPos(1:6) =  GetEMField(PEM%LocalElemID(PartID),PartPosRef_loc(1:3))
END IF
END FUNCTION GetInterpolatedFieldPartPos


PPURE FUNCTION GetEMField(ElemID,PartPosRef_loc)
!===================================================================================================================================
! Evaluate the electro-(magnetic) field using the reference position and return the field
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Eval_xyz      ,ONLY: EvaluateFieldAtRefPos
#if ! (USE_HDG)
USE MOD_DG_Vars       ,ONLY: U
#endif
#ifdef PP_POIS
USE MOD_Equation_Vars ,ONLY: E
#endif
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars ,ONLY: B
#else
USE MOD_Equation_Vars ,ONLY: B,E
#endif /*PP_nVar==1*/
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID !< Local element ID
REAL,INTENT(IN)    :: PartPosRef_loc(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GetEMField(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if defined PP_POIS || (USE_HDG && PP_nVar==4)
REAL :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
#endif /*(PP_POIS||USE_HDG)*/
!===================================================================================================================================
GetEMField(1:6)=0.
!--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),6,PP_N,HelperU,6,GetEMField(1:6),ElemID)
#else
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),6,PP_N,U(1:6,:,:,:,ElemID),6,GetEMField(1:6),ElemID)
#endif
#else
#ifdef PP_POIS
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),3,PP_N,E(1:3,:,:,:,ElemID),3,GetEMField(1:3),ElemID)
#elif USE_HDG
#if PP_nVar==1
#if (PP_TimeDiscMethod==507) || (PP_TimeDiscMethod==508)
! Boris or HC: consider B-Field, e.g., from SuperB
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),3,PP_N,E(1:3,:,:,:,ElemID),6,GetEMField(1:6),ElemID)
#else
! Consider only electric fields
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),3,PP_N,E(1:3,:,:,:,ElemID),3,GetEMField(1:3),ElemID)
#endif
#elif PP_nVar==3
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),3,PP_N,B(1:3,:,:,:,ElemID),3,GetEMField(4:6),ElemID)
#else
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),6,PP_N,HelperU,6,GetEMField(1:6),ElemID)
#endif
#else
CALL EvaluateFieldAtRefPos(PartPosRef_loc(1:3),3,PP_N,U(1:3,:,:,:,ElemID),3,GetEMField(1:3),ElemID)
#endif
#endif
END FUNCTION GetEMField


PPURE FUNCTION GetEMFieldDW(ElemID, PartPos_loc)
!===================================================================================================================================
! Evaluate the electro-(magnetic) field using the reference position and return the field
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars             ,ONLY: Elem_xGP
USE MOD_PICInterpolation_Vars ,ONLY: useBGField
USE MOD_Interpolation_Vars    ,ONLY: BGField,BGType,BGDataSize
USE MOD_Globals
USE MOD_PreProc
#if ! (USE_HDG)
USE MOD_DG_Vars       ,ONLY: U
#endif
#ifdef PP_POIS
USE MOD_Equation_Vars ,ONLY: E
#endif
#if USE_HDG
#if PP_nVar==1
USE MOD_Equation_Vars ,ONLY: E
#elif PP_nVar==3
USE MOD_Equation_Vars ,ONLY: B
#else
USE MOD_Equation_Vars ,ONLY: B,E
#endif /*PP_nVar==1*/
#endif /*USE_HDG*/
!----------------------------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: ElemID !< Local element ID
REAL,INTENT(IN)    :: PartPos_loc(1:3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: GetEMFieldDW(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: HelperU(1:6,0:PP_N,0:PP_N,0:PP_N)
REAL    :: PartDistDepo(0:PP_N,0:PP_N,0:PP_N), DistSum
INTEGER :: k,l,m,ind1,ind2
REAL    :: norm
!===================================================================================================================================
GetEMFieldDW(1:6)=0.
!--- evaluate at Particle position
#if (PP_nVar==8)
#ifdef PP_POIS
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
HelperU(4:6,:,:,:) = U(4:6,:,:,:,ElemID)
#else
HelperU(1:6,:,:,:) = U(1:6,:,:,:,ElemID)
#endif
#else
#ifdef PP_POIS
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
#elif USE_HDG
#if PP_nVar==1
#if (PP_TimeDiscMethod==507) || (PP_TimeDiscMethod==508)
! Boris or HC: consider B-Field, e.g., from SuperB
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
#else
! Consider only electric fields
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
#endif
#elif PP_nVar==3
HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
#else
HelperU(1:3,:,:,:) = E(1:3,:,:,:,ElemID)
HelperU(4:6,:,:,:) = B(1:3,:,:,:,ElemID)
#endif
#else
HelperU(1:3,:,:,:) = U(1:3,:,:,:,ElemID)
#endif
#endif

DistSum = 0.0
DO k = 0, PP_N; DO l=0, PP_N; DO m=0, PP_N
  norm = VECNORM(Elem_xGP(1:3,k,l,m, ElemID)-PartPos_loc(1:3))
  IF(norm.GT.0.)THEN
    PartDistDepo(k,l,m) = 1./norm
  ELSE
    PartDistDepo(:,:,:) = 0.
    PartDistDepo(k,l,m) = 1.
    DistSum = 1.
    EXIT
  END IF ! norm.GT.0.
  DistSum = DistSum + PartDistDepo(k,l,m) 
END DO; END DO; END DO

GetEMFieldDW = 0.0
DO k = 0, PP_N; DO l=0, PP_N; DO m=0, PP_N
  GetEMFieldDW(1:6) = GetEMFieldDW(1:6) + PartDistDepo(k,l,m)/DistSum*HelperU(1:6,k,l,m)
END DO; END DO; END DO

! Check whether magnetic background field is activated (superB)
IF(useBGField)THEN
  ! Check BG type and set dimensions
  SELECT CASE(BGType)
  CASE(1) ! Ex,Ey,Ez
    ind1 = 1
    ind2 = 3
  CASE(2) ! Bx,By,Bz
    ind1 = 4
    ind2 = 6
  CASE(3) ! Ex,Ey,Ez,Bx,By,Bz
    ind1 = 1
    ind2 = 6
  END SELECT
  ! Add contribution of the magnetic field
  DO k = 0, PP_N; DO l=0, PP_N; DO m=0, PP_N
    GetEMFieldDW(ind1:ind2) = GetEMFieldDW(ind1:ind2) + PartDistDepo(k,l,m)/DistSum*BGField(1:BGDataSize,k,l,m,ElemID)
  END DO; END DO; END DO
END IF ! useBGField

END FUNCTION GetEMFieldDW


PPURE FUNCTION InterpolateVariableExternalField1D(Pos)
!===================================================================================================================================
!> Interpolates the variable external field to the z-position
!> NO z-values smaller than VariableExternalField(1,1) are allowed!
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_PICInterpolation_Vars   ,ONLY:DeltaExternalField,nIntPoints,VariableExternalField
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos                               !< particle z-position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField1D  !< Bz (magnetic field in z-direction)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos                              !< index in array (equidistant subdivision assumed)
!===================================================================================================================================
iPos = INT((Pos-VariableExternalField(1,1))/DeltaExternalField(1)) + 1
IF(iPos.GE.nIntPoints)THEN ! particle outside of range (greater -> use constant value)
  InterpolateVariableExternalField1D = VariableExternalField(2,nIntPoints)
ELSEIF(iPos.LT.1)THEN ! particle outside of range (lower -> use constant value)
  InterpolateVariableExternalField1D = VariableExternalField(2,1)
ELSE ! Linear Interpolation between iPos and iPos+1 B point
  InterpolateVariableExternalField1D = (VariableExternalField(2,iPos+1) - VariableExternalField(2,iPos)) & !  dy
                                   / (VariableExternalField(1,iPos+1) - VariableExternalField(1,iPos)) & ! /dx
                             * (Pos - VariableExternalField(1,iPos) ) + VariableExternalField(2,iPos)    ! *(z - z_i) + z_i
END IF
END FUNCTION InterpolateVariableExternalField1D


PPURE FUNCTION InterpolateVariableExternalField2D(Pos)
!===================================================================================================================================
!> Interpolates the variable external field to the r- and z-position via bilinear interpolation
!>
!>   1.2 --------- 2.2
!>    |             |
!>  r |             |
!>    |             |
!>   1.1 --------- 2.1
!>           z
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars ,ONLY: DeltaExternalField,VariableExternalField
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldN
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldMin,VariableExternalFieldMax
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos(1:3)                                 !< particle z-position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField2D(1:3)  !< Bz (magnetic field in z-direction)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos,jPos                                !< index in array (equidistant subdivision assumed)
REAL                     :: r,delta,f(1:2),mat(2,2),dx(1:2),dy(1:2),vec(1:2)
INTEGER                  :: idx1,idx2,idx3,idx4,i
!===================================================================================================================================
ASSOCIATE(&
      x => Pos(1) ,&
      y => Pos(2) ,&
      z => Pos(3)  &
      )
  r = SQRT(x*x + y*y)

  IF(r.GT.VariableExternalFieldMax(1))THEN
    InterpolateVariableExternalField2D = 0.
  ELSEIF(r.LT.VariableExternalFieldMin(1))THEN
    InterpolateVariableExternalField2D = 0.
  ELSEIF(z.GT.VariableExternalFieldMax(2))THEN
    InterpolateVariableExternalField2D = 0.
  ELSEIF(z.LT.VariableExternalFieldMin(2))THEN
    InterpolateVariableExternalField2D = 0.
  ELSE

    ! Get index in r and z
    iPos = INT((r-VariableExternalField(1,1))/DeltaExternalField(1)) + 1 ! dr = DeltaExternalField(1)
    jPos = INT((z-VariableExternalField(2,1))/DeltaExternalField(2)) + 1 ! dz = DeltaExternalField(2)

    ! Catch problem when r or z are exactly at the upper boundary and INT() does not round to the lower integer (do not add +1 in
    ! this case)
    iPos = MIN(iPos, VariableExternalFieldN(2) - 1 )
    jPos = MIN(jPos, VariableExternalFieldN(1) - 1 )


    ! Shift all points by Nz = EmissionDistributionNum(1)
    ASSOCIATE( Nz => VariableExternalFieldN(1) )
      ! 1.1
      idx1 = (iPos-1)*Nz + jPos
      ! 2.1
      idx2 = (iPos-1)*Nz + jPos + 1
      ! 1.2
      idx3 = iPos*Nz + jPos
      ! 2.2
      idx4 = iPos*Nz + jPos + 1
    END ASSOCIATE

    ! Interpolate
    delta = DeltaExternalField(1)*DeltaExternalField(2)
    delta = 1./delta

    dx(1) = VariableExternalField(1,idx4)-r
    dx(2) = DeltaExternalField(1) - dx(1)

    dy(1) = VariableExternalField(2,idx4)-z
    dy(2) = DeltaExternalField(2) - dy(1)

    DO i = 1, 2
      mat(1,1) = VariableExternalField(2+i,idx1)
      mat(2,1) = VariableExternalField(2+i,idx2)
      mat(1,2) = VariableExternalField(2+i,idx3)
      mat(2,2) = VariableExternalField(2+i,idx4)

      vec(1) = dx(1)
      vec(2) = dx(2)
      f(i) = delta * DOT_PRODUCT( vec, MATMUL(mat, (/dy(1),dy(2)/) ) )
    END DO ! i = 1, 2

    ! Transform from Br, Bz to Bx, By, Bz
    r=1./r
    InterpolateVariableExternalField2D(1) = f(1)*x*r
    InterpolateVariableExternalField2D(2) = f(1)*y*r
    InterpolateVariableExternalField2D(3) = f(2)
  END IF ! r.GT.VariableExternalFieldMax(1)
END ASSOCIATE

END FUNCTION InterpolateVariableExternalField2D


PPURE FUNCTION InterpolateVariableExternalField3D(Pos)
!===================================================================================================================================
!> Interpolates the variable external field to the x-, y- and z-position via trilinear interpolation
!>
!>        1.2.2 ---------- 2.2.2
!>         /|               /|
!>        / |              / |
!>       /  |             /  |
!>      /   |            /   |
!>   1.2.1 ---------- 2.2.1  |
!>     |    |           |    |
!>   y |  1.1.2 --------|- 2.1.2
!>     |   /            |   /
!>     |  /             |  /
!>     | /              | / z
!>     |/               |/
!>   1.1.1 ---------- 2.1.1
!>            x
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PICInterpolation_Vars ,ONLY: DeltaExternalField,VariableExternalField
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldN
USE MOD_PICInterpolation_Vars ,ONLY: VariableExternalFieldMin,VariableExternalFieldMax
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)          :: Pos(1:3)                                 !< particle z-position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: InterpolateVariableExternalField3D(1:3)  !< Magnetic field B
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iPos,jPos,kPos                           !< index in array (equidistant subdivision assumed)
REAL,dimension(3)        :: c0,c1,c00,c01,c10,c11
INTEGER                  :: idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,Nxy
REAL                     :: xd,yd,zd
!===================================================================================================================================
ASSOCIATE(&
      x  => Pos(1) ,&
      y  => Pos(2) ,&
      z  => Pos(3) ,&
      Nx => VariableExternalFieldN(1)  ,&
      Ny => VariableExternalFieldN(2)  ,&
      Nz => VariableExternalFieldN(3)   &
      )

  ! Magnetic field outside of interpolation domain results in B=0
  IF(x.GT.VariableExternalFieldMax(1))THEN
    InterpolateVariableExternalField3D = 0.
  ELSEIF(x.LT.VariableExternalFieldMin(1))THEN
    InterpolateVariableExternalField3D = 0.
  ELSEIF(y.GT.VariableExternalFieldMax(2))THEN
    InterpolateVariableExternalField3D = 0.
  ELSEIF(y.LT.VariableExternalFieldMin(2))THEN
    InterpolateVariableExternalField3D = 0.
  ELSEIF(z.GT.VariableExternalFieldMax(3))THEN
    InterpolateVariableExternalField3D = 0.
  ELSEIF(z.LT.VariableExternalFieldMin(3))THEN
    InterpolateVariableExternalField3D = 0.
  ELSE

    ! Get index in x, y and z
    iPos = INT((x-VariableExternalField(1,1))/DeltaExternalField(1)) ! 0 to Nx-1
    jPos = INT((y-VariableExternalField(2,1))/DeltaExternalField(2)) ! 0 to Ny-1
    kPos = INT((z-VariableExternalField(3,1))/DeltaExternalField(3)) ! 0 to Nz-1
    ! Catch problem when coordinates are exactly at the upper boundary and INT() does not round to the lower integer
    ! e.g. when x.EQ.VariableExternalFieldMax(1) or y.EQ.VariableExternalFieldMax(2) or z.EQ.VariableExternalFieldMax(3)
    IF(iPos.EQ.Nx) iPos = Nx-1
    IF(jPos.EQ.Ny) jPos = Ny-1
    IF(kPos.EQ.Nz) kPos = Nz-1
    Nxy  = Nx*Ny

    ! Get corner node indices
    ! 1.1.1
    idx1 = iPos + jPos*Ny + kPos*Nxy + 1
    ! 2.1.1
    idx2 = idx1 + 1
    ! 1.2.1
    idx3 = idx1 + Ny ! iPos + (jPos+1)*Ny + kPos*Nxy +1
    ! 2.2.1
    idx4 = idx3 + 1

    ! 1.1.2
    idx5 = idx1 + Nxy
    ! 2.1.2
    idx6 = idx2 + Nxy
    ! 1.2.2
    idx7 = idx3 + Nxy
    ! 2.2.2
    idx8 = idx4 + Nxy

    ! Deltas
    xd = (x-VariableExternalField(1,idx1))/DeltaExternalField(1)
    yd = (y-VariableExternalField(2,idx1))/DeltaExternalField(2)
    zd = (z-VariableExternalField(3,idx1))/DeltaExternalField(3)

    ! Interpolate in x
    c00(1:3) = VariableExternalField(4:6,idx1)*(1.0-xd) + VariableExternalField(4:6,idx2)*xd
    c01(1:3) = VariableExternalField(4:6,idx3)*(1.0-xd) + VariableExternalField(4:6,idx4)*xd
    c10(1:3) = VariableExternalField(4:6,idx5)*(1.0-xd) + VariableExternalField(4:6,idx6)*xd
    c11(1:3) = VariableExternalField(4:6,idx7)*(1.0-xd) + VariableExternalField(4:6,idx8)*xd

    ! Interpolate in y: Note that c01 and c10 are switched
    c0(1:3) = c00(1:3)*(1.0-yd) + c01(1:3)*yd
    c1(1:3) = c10(1:3)*(1.0-yd) + c11(1:3)*yd

    ! Interpolate in z
    InterpolateVariableExternalField3D(1:3) = c0(1:3)*(1.0-zd) + c1(1:3)*zd
  END IF ! r.GT.VariableExternalFieldMax(1)
END ASSOCIATE

END FUNCTION InterpolateVariableExternalField3D


PPURE FUNCTION InterpolateAlgebraicExternalField(Pos)
!===================================================================================================================================
!> Interpolates the variable external field to the z-position
!> NO z-values smaller than VariableExternalField(1,1) are allowed!
!===================================================================================================================================
! MODULES
!USE MOD_Globals
USE MOD_PICInterpolation_Vars ,ONLY: externalField,AlgebraicExternalField,AlgebraicExternalFieldDelta
#if USE_HDG
USE MOD_Analyze_Vars          ,ONLY: AverageElectricPotential
#endif /*USE_HDG*/
USE MOD_Particle_Mesh_Vars    ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: Pos(1:3) !< particle position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: InterpolateAlgebraicExternalField(1:6)  !< E and B field at particle position
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: r !< radius factor
!===================================================================================================================================
ASSOCIATE(&
      x => Pos(1) ,&
      y => Pos(2) ,&
      z => Pos(3)  &
      )
  SELECT CASE (AlgebraicExternalField)
  CASE (1) ! Charoy 2019
#if USE_HDG
    ! Set Ex = Ue/xe
    ASSOCIATE( &
          Ue => AverageElectricPotential ,&
          xe => 2.4e-2                   )
      InterpolateAlgebraicExternalField(1) = Ue/xe
    END ASSOCIATE
#else
    InterpolateAlgebraicExternalField(1) = externalField(1) ! default
#endif /*USE_HDG*/

    ! Set Ey, Ez, Bx and By
    InterpolateAlgebraicExternalField(2:5) = externalField(2:5)

    ! Calc Bz
    ! Original formula
    !ASSOCIATE(&
    !      x     => Pos(1)   ,&
    !      Lx    => 2.50E-2  ,&
    !      xBmax => 0.750E-2 ,&
    !      B0    => 6.00E-3  ,&
    !      BLx   => 1.00E-3  ,&
    !      Bmax  => 10.0E-3  ,&
    !      sigma => 0.625E-2  &
    !      )
    ASSOCIATE( xBmax => 0.750E-2 )
      IF(Pos(1).LT.xBmax)THEN
        ! Original formula
        !ASSOCIATE(a1 => (Bmax-B0)/(1.0-EXP(-0.5*(xBmax/sigma)**2))                              ,&
        !          b1 => (B0 - Bmax*EXP(-0.5*(xBmax/sigma)**2))/(1.0-EXP(-0.5*(xBmax/sigma)**2))  )
        !  InterpolateAlgebraicExternalField(6) = a1 * EXP(-0.5*((x-xBmax)/sigma)**2) + b1
        !END ASSOCIATE
        InterpolateAlgebraicExternalField(6) = 7.7935072222899814E-003 * EXP(-12800.0*(x-xBmax)**2) + 2.2064927777100192E-003
      ELSE
        ! Original formula
        !ASSOCIATE(a2 => (Bmax-BLx)/(1.0-EXP(-0.5*((Lx-xBmax)/sigma)**2))                                  ,&
        !          b2 => (BLx - Bmax*EXP(-0.5*((Lx-xBmax)/sigma)**2))/(1.0-EXP(-0.5*((Lx-xBmax)/sigma)**2))  )
        !  InterpolateAlgebraicExternalField(6) = a2 * EXP(-0.5*((x-xBmax)/sigma)**2) + b2
        !END ASSOCIATE
        InterpolateAlgebraicExternalField(6) = 9.1821845944997683E-003 * EXP(-12800.0*(x-xBmax)**2) + 8.1781540550023306E-004
      END IF ! Pos(1).LT.xBmax
    END ASSOCIATE
    !END ASSOCIATE
  CASE (2) ! Liu 2010: 2D case
    ! Set Ex, Ey, Ez, Bx and Bz
    InterpolateAlgebraicExternalField(1:4) = externalField(1:4)
    InterpolateAlgebraicExternalField(6)   = externalField(6)

    ! Calculate By(x)=Br(x)=Br,0*(x/L)^delta
    ASSOCIATE( L => (GEO%xmaxglob-GEO%xminglob), Br0 => 400.0e-4, d => AlgebraicExternalFieldDelta ) ! Gauss to Tesla is *1e-4
      InterpolateAlgebraicExternalField(5) = Br0*((x/L)**d)
    END ASSOCIATE
  CASE (3) ! Liu 2010: 3D case
    ! Set Ex, Ey, Ez and Bz
    InterpolateAlgebraicExternalField(1:3) = externalField(1:3)
    InterpolateAlgebraicExternalField(6)   = externalField(6)

    ! Calculate By(z)=Br(z)=Br,0*(z/L)^delta
    ASSOCIATE( L => (GEO%zmaxglob-GEO%zminglob), Br0 => 400.0e-4, d => AlgebraicExternalFieldDelta ) ! Gauss to Tesla is *1e-4
      r=Br0*((z/L)**d)/SQRT(x**2+y**2)
      InterpolateAlgebraicExternalField(4) = x*r
      InterpolateAlgebraicExternalField(5) = y*r
    END ASSOCIATE
  END SELECT
END ASSOCIATE
END FUNCTION InterpolateAlgebraicExternalField


END MODULE MOD_PICInterpolation_tools
