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

MODULE MOD_PICDepo_Method
!===================================================================================================================================
! Module containing the different deposition methods (NGP, linear (inter-cell) weighting, shape function
!===================================================================================================================================
IMPLICIT NONE
PRIVATE

INTERFACE PartRHS
  PROCEDURE DepositionMethod
END INTERFACE

!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC :: DepositionMethod
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE DepositionMethodInterface(DoInnerParts,doPartInExists,doParticle_In)
    USE MOD_Particle_Vars ,ONLY: PDM
    LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
    LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
    LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
  END SUBROUTINE
END INTERFACE

PROCEDURE(DepositionMethodInterface),POINTER :: DepositionMethod    !< pointer defining the standard inner Riemann solver

INTEGER,PARAMETER      :: PRM_DEPO_SF   = 0  ! shape_function
INTEGER,PARAMETER      :: PRM_DEPO_SF1D = 1  ! shape_function_1d
INTEGER,PARAMETER      :: PRM_DEPO_SF2D = 2  ! shape_function_2d
INTEGER,PARAMETER      :: PRM_DEPO_SFC  = 3  ! shape_function_cylindrical
INTEGER,PARAMETER      :: PRM_DEPO_SFS  = 4  ! shape_function_spherical
INTEGER,PARAMETER      :: PRM_DEPO_SFSi = 5  ! shape_function_spherical
INTEGER,PARAMETER      :: PRM_DEPO_CVW  = 6  ! cell_volweight
INTEGER,PARAMETER      :: PRM_DEPO_NGP  = 7  ! nearest_gausspoint
INTEGER,PARAMETER      :: PRM_DEPO_DD   = 8  ! delta_distri
INTEGER,PARAMETER      :: PRM_DEPO_MVW  = 9  ! cartmesh_volumeweighting
INTEGER,PARAMETER      :: PRM_DEPO_MS   = 10 ! cartmesh_splines
INTEGER,PARAMETER      :: PRM_DEPO_NBC  = 11 ! nearest_blurrycenter
INTEGER,PARAMETER      :: PRM_DEPO_CVWM = 12 ! cell_volweight_mean

INTERFACE InitDepositionMethod
  MODULE PROCEDURE InitDepositionMethod
END INTERFACE

PUBLIC :: InitDepositionMethod
!==================================================================================================================================

PUBLIC :: DefineParametersDepositionMethod
CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%CreateLogicalOption('PIC-DoDeposition', 'Switch deposition of charge (and current density) on/off', '.TRUE.')

CALL prms%CreateIntFromStringOption('PIC-Deposition-Type', "Type/Method used in the deposition step: \n"//&
                                                             '1.1)  shape_function\n'                   //&
                                                             '1.2)  shape_function_1d\n'                //&
                                                             '1.3)  shape_function_2d\n'                //&
                                                             '1.4)  shape_function_cylindrical\n'       //&
                                                             '1.5)  shape_function_spherical\n'         //&
                                                             '1.6)  shape_function_simple\n'            //&
                                                             '      1.1) to 1.6) require\n'             //&
                                                             '        PIC-shapefunction-radius\n'       //&
                                                             '        PIC-shapefunction-alpha\n'        //&
                                                             '      1.2) and 1.3) require\n'            //&
                                                             '        PIC-shapefunction1d-direction\n'  //&
                                                             '      1.4) and 1.5) require\n'            //&
                                                             '        PIC-shapefunction-radius0\n'      //&
                                                             '        PIC-shapefunction-scale\n'        //&
                                                             '2.)   cell_volweight\n'                   //&
                                                             '3.)   nearest_gausspoint\n'               //&
                                                             '4.)   delta_distri\n'                     //&
                                                             '      requires PIC-DeltaType\n'           //&
                                                             '               PIC-DeltaType-N\n'         //&
                                                             '5.1)  cartmesh_volumeweighting\n'         //&
                                                             '5.2)  cartmesh_splines\n'                 //&
                                                             '      requires PIC-BGMdeltas\n'           //&
                                                             '               PIC-FactorBGM\n'           //&
                                                             '6.)   nearest-blurrycenter\n'             //&
                                                             '7.)   cell_volweight_mean'                &
                                                       ,"cell_volweight")

CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function'             , PRM_DEPO_SF)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_1d'          , PRM_DEPO_SF1D)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_2d'          , PRM_DEPO_SF2D)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_cylindrical' , PRM_DEPO_SFC)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_spherical'   , PRM_DEPO_SFS)
CALL addStrListEntry('PIC-Deposition-Type' , 'shape_function_simple'      , PRM_DEPO_SFSi)
CALL addStrListEntry('PIC-Deposition-Type' , 'cell_volweight'             , PRM_DEPO_CVW)
CALL addStrListEntry('PIC-Deposition-Type' , 'nearest_gausspoint'         , PRM_DEPO_NGP)
CALL addStrListEntry('PIC-Deposition-Type' , 'delta_distri'               , PRM_DEPO_DD)
CALL addStrListEntry('PIC-Deposition-Type' , 'cartmesh_volumeweighting'   , PRM_DEPO_MVW)
CALL addStrListEntry('PIC-Deposition-Type' , 'cartmesh_splines'           , PRM_DEPO_MS)
CALL addStrListEntry('PIC-Deposition-Type' , 'nearest_blurrycenter'       , PRM_DEPO_NBC)
CALL addStrListEntry('PIC-Deposition-Type' , 'cell_volweight_mean'        , PRM_DEPO_CVWM)
END SUBROUTINE DefineParametersDepositionMethod


!==================================================================================================================================!
!> Initialize particle RHS functions
!==================================================================================================================================!
SUBROUTINE InitDepositionMethod()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools            ,ONLY: GETINTFROMSTR
USE MOD_PICDepo_Vars           ,ONLY: DepositionType
USE MOD_PICInterpolation_Vars  ,ONLY: InterpolationType
USE MOD_Particle_Tracking_Vars ,ONLY: TriaTracking
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: DepositionType_loc
!==================================================================================================================================
DepositionType_loc = GETINTFROMSTR('PIC-Deposition-Type')
! check for interpolation type incompatibilities (cannot be done at interpolation_init
! because DepositionType_loc is not known yet)
IF((TRIM(InterpolationType).EQ.'nearest_gausspoint').AND. &
   (DepositionType_loc.NE.PRM_DEPO_NGP)) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in pic_depo.f90: Interpolation type nearest_gausspoint only allowed with same deposition type!')
END IF
IF((DepositionType_loc.EQ.PRM_DEPO_CVWM).AND. &
   (.NOT.(TriaTracking))) THEN
  CALL abort(&
  __STAMP__&
  ,'ERROR in pic_depo.f90: PIC-Deposition-Type = cell_volweight_mean only allowed with TriaTracking!')
END IF


SELECT CASE(DepositionType_loc)
Case(PRM_DEPO_SF) ! shape_function
  DepositionMethod => DepositionMethod_SF
  DepositionType   = 'shape_function'
Case(PRM_DEPO_SFSi) ! shape_function_simple
  DepositionType   = 'shape_function_simple'
  DepositionMethod => DepositionMethod_SF
Case(PRM_DEPO_SF1D) ! shape_function_1d
  DepositionType   = 'shape_function_1d'
  DepositionMethod => DepositionMethod_SF1D
Case(PRM_DEPO_SF2D) ! shape_function_2d
  DepositionType   = 'shape_function_2d'
  DepositionMethod => DepositionMethod_SF2D
Case(PRM_DEPO_SFC) ! shape_function_cylindrical
  DepositionType   = 'shape_function_cylindrical'
  DepositionMethod => DepositionMethod_SFCS
Case(PRM_DEPO_SFS) ! shape_function_spherical
  DepositionType   = 'shape_function_spherical'
  DepositionMethod => DepositionMethod_SFCS
Case(PRM_DEPO_CVW) ! cell_volweight
  DepositionType   = 'cell_volweight'
  DepositionMethod => DepositionMethod_CVW
Case(PRM_DEPO_NGP) ! nearest_gausspoint
  DepositionType   = 'nearest_gausspoint'
  DepositionMethod => DepositionMethod_NGP
Case(PRM_DEPO_DD) ! delta_distri
  DepositionType   = 'delta_distri'
  DepositionMethod => DepositionMethod_DD
Case(PRM_DEPO_MVW) ! cartmesh_volumeweighting
  DepositionType   = 'cartmesh_volumeweighting'
  DepositionMethod => DepositionMethod_MVW
Case(PRM_DEPO_MS) ! cartmesh_splines
  DepositionType   = 'cartmesh_splines'
  DepositionMethod => DepositionMethod_MS
Case(PRM_DEPO_NBC) ! nearest_blurrycenter
  DepositionType   = 'nearest_blurrycenter'
  DepositionMethod => DepositionMethod_NBC
Case(PRM_DEPO_CVWM) ! cell_volweight_mean
  DepositionType   = 'cell_volweight_mean'
  DepositionMethod => DepositionMethod_CVWM
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
      'Unknown DepositionMethod!' ,IntInfo=DepositionType_loc)
      !'Unknown DepositionMethod! ['//TRIM(DepositionType)//']')
END SELECT

! Suppress compiler warning
RETURN
CALL DepositionMethod_SF(.FALSE.)
CALL DepositionMethod_SF1D(.FALSE.)
CALL DepositionMethod_SF2D(.FALSE.)
CALL DepositionMethod_SFCS(.FALSE.)
CALL DepositionMethod_CVW(.FALSE.)
CALL DepositionMethod_NGP(.FALSE.)
CALL DepositionMethod_DD(.FALSE.)
CALL DepositionMethod_MVW(.FALSE.)
CALL DepositionMethod_MS(.FALSE.)
CALL DepositionMethod_NBC(.FALSE.)
CALL DepositionMethod_CVWM(.FALSE.)

END SUBROUTINE InitDepositionMethod


SUBROUTINE DepositionMethod_NGP(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! 'non-relativistic'
! Particle Right-Hand-Side: Non-relativistic push
! Former FUNCTION NON_RELATIVISTIC_PUSH
!===================================================================================================================================
! MODULES
USE MOD_PreProc               ,ONLY: PP_N,PP_nElems
USE MOD_Particle_Vars         ,ONLY: Species, PartSpecies,PDM,PEM,usevMPF,PartPosRef,PartMPF
USE MOD_PICDepo_Vars          ,ONLY: PartSource,gaussborder
!#if (PP_nVar==8)
USE MOD_Particle_Vars         ,ONLY: PartState
!#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers    ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers    ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Interpolation_Vars    ,ONLY: wGP
USE MOD_Mesh_Vars             ,ONLY: sJ
USE MOD_Part_Tools            ,ONLY: isDepositParticle
USE MOD_PICInterpolation_Vars ,ONLY: InterpolationType
USE MOD_Eval_xyz              ,ONLY: GetPositionInRefElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: prefac
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
LOGICAL            :: SAVE_GAUSS
INTEGER            :: a,b, ii
INTEGER            :: k, l, m, iElem, iPart
INTEGER            :: firstPart,LastPart
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
SAVE_GAUSS = .FALSE.
IF(TRIM(InterpolationType).EQ.'nearest_gausspoint') SAVE_GAUSS = .TRUE.
IF(MOD(PP_N,2).EQ.0) THEN
  a = PP_N/2
  b = a
ELSE
  a = (PP_N+1)/2
  b = a-1
END IF
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem=1,PP_nElems
  DO iPart=firstPart,LastPart
    ! TODO: Info why and under which conditions the following 'CYCLE' is called
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! Don't deposit neutral particles!
    IF(.NOT.isDepositParticle(iPart)) CYCLE
    IF(PEM%Element(iPart).EQ.iElem)THEN
      IF (usevMPF) THEN
        prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
      ELSE
        prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
      END IF ! usevMPF
      ! Map Particle to -1|1 space (re-used in interpolation)
      !IF(.NOT.DoRefMapping)THEN
      !  CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),iElem,iPart)
      !END IF
      ! Find out which gausspoint is closest and add up charges and currents
      !! x-direction
      IF(.NOT.SAVE_GAUSS) THEN
        k = a
        DO ii = 0,b-1
          IF(ABS(PartPosRef(1,iPart)).GE.GaussBorder(PP_N-ii))THEN
            k = PP_N-ii
            EXIT
          END IF
        END DO
        k = NINT((PP_N+SIGN(2.0*k-PP_N,PartPosRef(1,iPart)))/2)
        !! y-direction
        l = a
        DO ii = 0,b-1
          IF(ABS(PartPosRef(2,iPart)).GE.GaussBorder(PP_N-ii))THEN
            l = PP_N-ii
            EXIT
          END IF
        END DO
        l = NINT((PP_N+SIGN(2.0*l-PP_N,PartPosRef(2,iPart)))/2)
        !! z-direction
        m = a
        DO ii = 0,b-1
          IF(ABS(PartPosRef(3,iPart)).GE.GaussBorder(PP_N-ii))THEN
            m = PP_N-ii
            EXIT
          END IF
        END DO
        m = NINT((PP_N+SIGN(2.0*m-PP_N,PartPosRef(3,iPart)))/2)
      END IF
!#if (PP_nVar==8)
      PartSource(1:3,k,l,m,iElem) = PartSource(1:3,k,l,m,iElem) + PartState(4:6,iPart) * prefac
!#endif
      PartSource( 4 ,k,l,m,iElem) = PartSource( 4 ,k,l,m,iElem) + prefac
      !IF (SAVE_GAUSS) THEN
      !  PartPosGauss(iPart,1) = k
      !  PartPosGauss(iPart,2) = l
      !  PartPosGauss(iPart,3) = m
      !END IF
    END IF ! Element .EQ. iElem
  END DO ! iPart
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iElem=1,PP_nElems
IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
    ! get densities by dividing by gauss volume
!#if (PP_nVar==8)
      PartSource(1:4,k,l,m,iElem) = PartSource(1:4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#else
!        PartSource(4,k,l,m,iElem) = PartSource(4,k,l,m,iElem) * sJ(k,l,m,iElem)/(wGP(k)*wGP(l)*wGP(m))
!#endif
    END DO; END DO; END DO
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! iElem=1,PP_nElems
END IF


!SELECT CASE(TRIM(DepositionType))
!CASE('nearest_blurrycenter')
!CASE('cell_volweight')
!CASE('cell_volweight_mean','cell_volweight_mean2')
!CASE('epanechnikov')
!CASE('shape_function','shape_function_simple')
!CASE('shape_function_1d')
!CASE('shape_function_2d')
!CASE('shape_function_cylindrical','shape_function_spherical')
!CASE('delta_distri')
!CASE('nearest_gausspoint')
!CASE('cartmesh_volumeweighting')
!CASE('cartmesh_splines')
!CASE DEFAULT
  !CALL abort(&
  !__STAMP__&
  !,'Unknown DepositionType in pic_depo.f90')
!END SELECT



END SUBROUTINE DepositionMethod_NGP


SUBROUTINE DepositionMethod_NBC(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! 'nearest_blurrycenter'
! Deposits the complete particle charge at the center of the cell -> cell-constant deposition
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_nElems
USE MOD_Particle_Vars      ,ONLY: Species,PartSpecies,PDM,PEM,usevMPF,PartMPF
!#if (PP_nVar==8)
USE MOD_Particle_Vars      ,ONLY: PartState
!#endif
USE MOD_PICDepo_Vars       ,ONLY: PartSource
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Part_Tools         ,ONLY: isDepositParticle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE  :: ElemSource(:,:)
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER :: SourceDim=1
LOGICAL, PARAMETER :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL            :: doCalculateCurrentDensity
INTEGER            :: SourceDim
#endif
INTEGER            :: iElem,iPart
INTEGER            :: firstPart,LastPart
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
ALLOCATE(ElemSource(SourceDim:4,1:nElems))
ElemSource=0.0
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem=1,PP_nElems
  DO iPart=firstPart,lastPart
    ! TODO: Info why and under which conditions the following 'CYCLE' is called
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! Don't deposit neutral particles!
    IF(.NOT.isDepositParticle(iPart)) CYCLE
    IF(PEM%Element(iPart).EQ.iElem)THEN
      IF(usevMPF)THEN
        IF(doCalculateCurrentDensity)THEN
          ElemSource(1:3,iElem) = ElemSource(1:3,iElem)+ PartState(4:6,iPart)*Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
        END IF
        ElemSource(4,iElem) = ElemSource(4,iElem)  + Species(PartSpecies(iPart))%ChargeIC* PartMPF(iPart)
      ELSE
        IF(doCalculateCurrentDensity)THEN
          ElemSource(1:3,iElem) = ElemSource(1:3,iElem)+ &
              PartState(4:6,iPart)* Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
        END IF
        ElemSource(4,iElem) = ElemSource(4,iElem) + &
            Species(PartSpecies(iPart))%ChargeIC* Species(PartSpecies(iPart))%MacroParticleFactor
      END IF ! usevMPF
    END IF ! Element(iPart).EQ.iElem
  END DO ! iPart
  IF(doCalculateCurrentDensity)THEN
    PartSource(1,:,:,:,iElem) = PartSource(1,:,:,:,iElem)+ElemSource(1,iElem)
    PartSource(2,:,:,:,iElem) = PartSource(2,:,:,:,iElem)+ElemSource(2,iElem)
    PartSource(3,:,:,:,iElem) = PartSource(3,:,:,:,iElem)+ElemSource(3,iElem)
  END IF
  PartSource(4,:,:,:,iElem) = PartSource(4,:,:,:,iElem)+ElemSource(4,iElem)
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iElem=1,PP_nElems

IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    PartSource(SourceDim:4,:,:,:,iElem) = PartSource(SourceDim:4,:,:,:,iElem) / GEO%Volume(iElem)
  END DO ! iElem=1,PP_nElems
#if USE_LOADBALANCE
  CALL LBElemPauseTime_avg(tLBStart) ! Average over the number of elems
#endif /*USE_LOADBALANCE*/
END IF ! .NOT. DoInnerParts
END SUBROUTINE DepositionMethod_NBC




SUBROUTINE DepositionMethod_CVW(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! 'nearest_blurrycenter'
! Deposits the complete particle charge at the center of the cell -> cell-constant deposition
!===================================================================================================================================
! MODULES
USE MOD_PreProc                ,ONLY: PP_N
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies,PDM,PEM,PartPosRef,usevMPF,PartMPF
USE MOD_Particle_Vars          ,ONLY: PartState
USE MOD_PICDepo_Vars           ,ONLY: PartSource,CellVolWeight_Volumes
USE MOD_Part_Tools             ,ONLY: isDepositParticle
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers     ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_PICDepo_Vars           ,ONLY: PartSource,CellVolWeightFac
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE  :: BGMSourceCellVol(:,:,:,:,:)
REAL               :: Charge, TSource(1:4)
REAL               :: alpha1, alpha2, alpha3, TempPartPos(1:3)
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER :: SourceDim=1
LOGICAL, PARAMETER :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL            :: doCalculateCurrentDensity
INTEGER            :: SourceDim
#endif
INTEGER            :: kk, ll, mm
INTEGER            :: firstPart,LastPart,iPart,iElem
!===================================================================================================================================
! Return here for 2nd Deposition() call as it is not required for this deposition method, 
! because the MPI communication is done here directly
IF(.NOT.DoInnerParts) RETURN

ALLOCATE(BGMSourceCellVol(SourceDim:4,0:1,0:1,0:1,1:nElems))
BGMSourceCellVol(:,:,:,:,:) = 0.0
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iPart = firstPart, lastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  IF (usevMPF) THEN
    Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
  ELSE
    Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF ! usevMPF
  iElem = PEM%Element(iPart)
  IF(DoRefMapping)THEN
    TempPartPos(1:3)=PartPosRef(1:3,iPart)
  ELSE
    CALL GetPositionInRefElem(PartState(1:3,iPart),TempPartPos,iElem,ForceMode=.TRUE.)
  END IF
  IF(doCalculateCurrentDensity)THEN
    TSource(1:3) = PartState(4:6,iPart)*Charge
  ELSE
    TSource(1:3) = 0.0
  END IF
  TSource(4) = Charge
  alpha1=(TempPartPos(1)+1.0)/2.0
  alpha2=(TempPartPos(2)+1.0)/2.0
  alpha3=(TempPartPos(3)+1.0)/2.0
  BGMSourceCellVol(:,0,0,0,iElem) = BGMSourceCellVol(:,0,0,0,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
  BGMSourceCellVol(:,0,0,1,iElem) = BGMSourceCellVol(:,0,0,1,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(alpha3))
  BGMSourceCellVol(:,0,1,0,iElem) = BGMSourceCellVol(:,0,1,0,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(alpha2)*(1-alpha3))
  BGMSourceCellVol(:,0,1,1,iElem) = BGMSourceCellVol(:,0,1,1,iElem) + (TSource(SourceDim:4)*(1-alpha1)*(alpha2)*(alpha3))
  BGMSourceCellVol(:,1,0,0,iElem) = BGMSourceCellVol(:,1,0,0,iElem) + (TSource(SourceDim:4)*(alpha1)*(1-alpha2)*(1-alpha3))
  BGMSourceCellVol(:,1,0,1,iElem) = BGMSourceCellVol(:,1,0,1,iElem) + (TSource(SourceDim:4)*(alpha1)*(1-alpha2)*(alpha3))
  BGMSourceCellVol(:,1,1,0,iElem) = BGMSourceCellVol(:,1,1,0,iElem) + (TSource(SourceDim:4)*(alpha1)*(alpha2)*(1-alpha3))
  BGMSourceCellVol(:,1,1,1,iElem) = BGMSourceCellVol(:,1,1,1,iElem) + (TSource(SourceDim:4)*(alpha1)*(alpha2)*(alpha3))
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to iElem
#endif /*USE_LOADBALANCE*/
END DO

DO iElem=1, nElems
  BGMSourceCellVol(:,0,0,0,iElem) = BGMSourceCellVol(:,0,0,0,iElem)/CellVolWeight_Volumes(0,0,0,iElem)
  BGMSourceCellVol(:,0,0,1,iElem) = BGMSourceCellVol(:,0,0,1,iElem)/CellVolWeight_Volumes(0,0,1,iElem)
  BGMSourceCellVol(:,0,1,0,iElem) = BGMSourceCellVol(:,0,1,0,iElem)/CellVolWeight_Volumes(0,1,0,iElem)
  BGMSourceCellVol(:,0,1,1,iElem) = BGMSourceCellVol(:,0,1,1,iElem)/CellVolWeight_Volumes(0,1,1,iElem)
  BGMSourceCellVol(:,1,0,0,iElem) = BGMSourceCellVol(:,1,0,0,iElem)/CellVolWeight_Volumes(1,0,0,iElem)
  BGMSourceCellVol(:,1,0,1,iElem) = BGMSourceCellVol(:,1,0,1,iElem)/CellVolWeight_Volumes(1,0,1,iElem)
  BGMSourceCellVol(:,1,1,0,iElem) = BGMSourceCellVol(:,1,1,0,iElem)/CellVolWeight_Volumes(1,1,0,iElem)
  BGMSourceCellVol(:,1,1,1,iElem) = BGMSourceCellVol(:,1,1,1,iElem)/CellVolWeight_Volumes(1,1,1,iElem)
END DO

DO iElem = 1, nElems
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
        alpha1 = CellVolWeightFac(kk)
        alpha2 = CellVolWeightFac(ll)
        alpha3 = CellVolWeightFac(mm)
        PartSource(SourceDim:4,kk,ll,mm,iElem) =PartSource(SourceDim:4,kk,ll,mm,iElem) + &
            BGMSourceCellVol(:,0,0,0,iElem) * (1-alpha1) * (1-alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,0,0,1,iElem) * (1-alpha1) * (1-alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,0,1,0,iElem) * (1-alpha1) *   (alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,0,1,1,iElem) * (1-alpha1) *   (alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,1,0,0,iElem) *   (alpha1) * (1-alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,1,0,1,iElem) *   (alpha1) * (1-alpha2) *   (alpha3)    + &
            BGMSourceCellVol(:,1,1,0,iElem) *   (alpha1) *   (alpha2) * (1-alpha3)    + &
            BGMSourceCellVol(:,1,1,1,iElem) *   (alpha1) *   (alpha2) *   (alpha3)
      END DO ! mm
    END DO ! ll
  END DO ! kk
END DO ! iElem
#if USE_LOADBALANCE
CALL LBElemSplitTime_avg(tLBStart) ! Average over the number of elems (and Start again)
#endif /*USE_LOADBALANCE*/
DEALLOCATE(BGMSourceCellVol)
END SUBROUTINE DepositionMethod_CVW



SUBROUTINE DepositionMethod_CVWM(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits the complete particle charge at the center of the cell -> cell-constant deposition
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N
USE MOD_Particle_Vars      ,ONLY: Species, PartSpecies,PDM,PEM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_PICDepo_Vars       ,ONLY: PartSource,CellVolWeightFac,NodeSourceExtTmp,NodeSourceExt,CellLocNodes_Volumes,DepositionType
#if USE_MPI
USE MOD_Particle_MPI       ,ONLY: AddHaloNodeData
#endif  /*USE_MPI*/
USE MOD_Mesh_Vars          ,ONLY: nElems, nNodes
USE MOD_Eval_xyz           ,ONLY: GetPositionInRefElem
USE MOD_Dielectric_Vars    ,ONLY: DoDielectricSurfaceCharge
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime,LBElemSplitTime,LBElemPauseTime_avg
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL, ALLOCATABLE  :: NodeSource(:,:), tempNodeSource(:,:)
REAL               :: Charge, TSource(1:4)
REAL               :: alpha1, alpha2, alpha3, TempPartPos(1:3)
INTEGER            :: NodeID(1:8)
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER :: SourceDim=1
LOGICAL, PARAMETER :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL            :: doCalculateCurrentDensity
INTEGER            :: SourceDim
#endif
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER            :: kk, ll, mm, iPart,iElem
!===================================================================================================================================
! Return here for 2nd Deposition() call as it is not required for this deposition method, 
! because the MPI communication is done here directly
IF(.NOT.DoInnerParts) RETURN

! Allocate NodeSource array and deallocate at the end of this procedure
ALLOCATE(NodeSource(SourceDim:4,1:nNodes))
NodeSource = 0.0

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iPart=1,PDM%ParticleVecLength
  IF (PDM%ParticleInside(iPart)) THEN
    IF (usevMPF) THEN
      Charge = Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)
    ELSE
      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
    END IF
    iElem = PEM%Element(iPart)
    CALL GetPositionInRefElem(PartState(1:3,iPart),TempPartPos(1:3),iElem,ForceMode=.TRUE.)
    TSource(:) = 0.0
    IF(doCalculateCurrentDensity)THEN
      TSource(1:3) = PartState(4:6,iPart)*Charge
    END IF
    TSource(4) = Charge

    alpha1=0.5*(TempPartPos(1)+1.0)
    alpha2=0.5*(TempPartPos(2)+1.0)
    alpha3=0.5*(TempPartPos(3)+1.0)
    NodeID=GEO%ElemToNodeID(1:8,iElem)
    NodeSource(:,NodeID(1)) = NodeSource(:,NodeID(1))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
    NodeSource(:,NodeID(2)) = NodeSource(:,NodeID(2))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*(1-alpha3))
    NodeSource(:,NodeID(3)) = NodeSource(:,NodeID(3))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*(1-alpha3))
    NodeSource(:,NodeID(4)) = NodeSource(:,NodeID(4))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*(1-alpha3))
    NodeSource(:,NodeID(5)) = NodeSource(:,NodeID(5))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*  (alpha3))
    NodeSource(:,NodeID(6)) = NodeSource(:,NodeID(6))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*  (alpha3))
    NodeSource(:,NodeID(7)) = NodeSource(:,NodeID(7))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*  (alpha3))
    NodeSource(:,NodeID(8)) = NodeSource(:,NodeID(8))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*  (alpha3))
#if USE_LOADBALANCE
   CALL LBElemSplitTime(iElem,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to iElem
#endif /*USE_LOADBALANCE*/
  END IF
END DO


! Node MPI communication
#if USE_MPI
IF(doCalculateCurrentDensity)THEN
  CALL AddHaloNodeData(NodeSource(1,:))
  CALL AddHaloNodeData(NodeSource(2,:))
  CALL AddHaloNodeData(NodeSource(3,:))
END IF
CALL AddHaloNodeData(NodeSource(4,:))

! Communicate dielectric surface charges stored in NodeSourceExtTmp
IF(DoDielectricSurfaceCharge)THEN
  CALL AddHaloNodeData(NodeSourceExtTmp)
END IF ! DoDielectricSurfaceCharge
#endif /*USE_MPI*/

IF(DoDielectricSurfaceCharge)THEN
  ! Update external node source containing dielectric surface charges and nullify
  NodeSourceExt    = NodeSourceExt + NodeSourceExtTmp
  NodeSourceExtTmp = 0.

  ! Add external node source (e.g. surface charging)
  NodeSource(4,:) = NodeSource(4,:) + NodeSourceExt
END IF ! DoDielectricSurfaceCharge



! Currently also "Nodes" are included in time measurement that is averaged across all elements. Can this be improved?
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem=1, nNodes
  NodeSource(SourceDim:4,iElem) = NodeSource(SourceDim:4,iElem)/CellLocNodes_Volumes(iElem)
END DO

IF (TRIM(DepositionType).EQ.'cell_volweight_mean2') THEN
  ALLOCATE(tempNodeSource(SourceDim:4,1:nNodes))
  tempNodeSource = 0.0
  DO iElem=1, nNodes
    tempNodeSource(SourceDim:4,iElem) = NodeSource(SourceDim:4,iElem)
    DO kk =1, GEO%NeighNodesOnNode(iElem)
      tempNodeSource(SourceDim:4,iElem) = tempNodeSource(SourceDim:4,iElem) + NodeSource(SourceDim:4,GEO%NodeToNeighNode(iElem)%ElemID(kk))
    END DO
    tempNodeSource(SourceDim:4,iElem) = tempNodeSource(SourceDim:4,iElem) / (GEO%NeighNodesOnNode(iElem) + 1.0)
  END DO
  NodeSource = tempNodeSource
END IF

! Interpolate node source values to volume polynomial
DO iElem = 1, nElems
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
        alpha1 = CellVolWeightFac(kk)
        alpha2 = CellVolWeightFac(ll)
        alpha3 = CellVolWeightFac(mm)
        NodeID=GEO%ElemToNodeID(1:8,iElem)
        Partsource(SourceDim:4,kk,ll,mm,iElem) = &
             NodeSource(SourceDim:4,NodeID(1)) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(2)) * (alpha1) * (1-alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(3)) * (alpha1) * (alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(4)) * (1-alpha1) * (alpha2) * (1-alpha3) + &
             NodeSource(SourceDim:4,NodeID(5)) * (1-alpha1) * (1-alpha2) * (alpha3) + &
             NodeSource(SourceDim:4,NodeID(6)) * (alpha1) * (1-alpha2) * (alpha3) + &
             NodeSource(SourceDim:4,NodeID(7)) * (alpha1) * (alpha2) * (alpha3) + &
             NodeSource(SourceDim:4,NodeID(8)) * (1-alpha1) * (alpha2) * (alpha3)
      END DO !mm
    END DO !ll
  END DO !kk
END DO !iEle
#if USE_LOADBALANCE
CALL LBElemPauseTime_avg(tLBStart) ! Average over the number of elems
#endif /*USE_LOADBALANCE*/
DEALLOCATE(NodeSource)

! Suppress compiler warning
RETURN
IF(doPartInExists.AND.doParticle_In(1))kk=0
END SUBROUTINE DepositionMethod_CVWM


SUBROUTINE DepositionMethod_SF(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc                     ,ONLY: PP_N,PP_nElems
USE MOD_globals                     ,ONLY: abort
USE MOD_Particle_Vars               ,ONLY: Species, PartSpecies,PDM,PartMPF,usevMPF
USE MOD_Particle_Vars               ,ONLY: PartState
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO
USE MOD_PICDepo_Vars                ,ONLY: PartSource,SFdepoLayersGeo,SFdepoLayersBaseVector,LastAnalyzeSurfCollis,PartSourceConst
USE MOD_PICDepo_Vars                ,ONLY: RelaxFac,sfdepofixesgeo,SFdepoLayersSpace,sfdepolayersbounds,SFdepoLayersUseFixBounds
USE MOD_PICDepo_Vars                ,ONLY: sfdepolayersradius,sfdepolayerspartnum,PartSourceOld,w_sf,SFdepoLayersMPF,SFdepoLayersSpec
USE MOD_PICDepo_Vars                ,ONLY: Vdm_EquiN_GaussN,SFResampleAnalyzeSurfCollis,SFdepoLayersAlreadyDone,RelaxDeposition,r_SF
USE MOD_PICDepo_Vars                ,ONLY: PartSourceConstExists,NbrOfSFdepoLayers,NbrOfSFdepoFixes,DoSFEqui
USE MOD_PICDepo_Vars                ,ONLY: ConstantSFdepoLayers
USE MOD_PICDepo_Shapefunction_Tools ,ONLY: calcSfSource
USE MOD_Particle_MPI_Vars           ,ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM,ExtPartMPF,NbrOfextParticles
USE MOD_TimeDisc_Vars               ,ONLY: dtWeight
USE MOD_Part_Tools                  ,ONLY: isDepositParticle
USE MOD_Mesh_Vars                   ,ONLY: nElems
USE MOD_ChangeBasis                 ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
REAL               :: RandVal, RandVal2(2), layerPartPos(3), PartRadius, FractPush(3), SFfixDistance
INTEGER            :: iElem, iPart, iPart2, iSFfix
INTEGER            :: iLayer, layerParts
INTEGER            :: firstPart,LastPart
INTEGER            :: kk, ll, mm
LOGICAL            :: DoCycle
!===================================================================================================================================
!-- "normal" particles
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
Vec1(1:3) = 0.
Vec2(1:3) = 0.
Vec3(1:3) = 0.
IF (GEO%nPeriodicVectors.EQ.1) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
END IF
IF (GEO%nPeriodicVectors.EQ.2) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
END IF
IF (GEO%nPeriodicVectors.EQ.3) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
END IF
IF (usevMPF) THEN
  DO iPart=firstPart,LastPart
    ! TODO: Info why and under which conditions the following 'CYCLE' is called
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! Don't deposit neutral particles!
    IF(.NOT.isDepositParticle(iPart)) CYCLE
    CALL calcSfSource(4,Species(PartSpecies(iPart))%ChargeIC*PartMPF(iPart)*w_sf &
      ,Vec1,Vec2,Vec3,PartState(1:3,iPart),iPart,PartVelo=PartState(4:6,iPart))
  END DO ! iPart
ELSE
  DO iPart=firstPart,LastPart
    ! TODO: Info why and under which conditions the following 'CYCLE' is called
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    ! Don't deposit neutral particles!
    IF(.NOT.isDepositParticle(iPart)) CYCLE
    CALL calcSfSource(4,Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor*w_sf &
      ,Vec1,Vec2,Vec3,PartState(1:3,iPart),iPart,PartVelo=PartState(4:6,iPart))
  END DO ! iPart
END IF ! usevMPF
IF(.NOT.DoInnerParts)THEN

  !-- layer particles (only once, i.e., during call with .NOT.DoInnerParts)
  DO iLayer=1,NbrOfSFdepoLayers
    IF (SFdepoLayersAlreadyDone) EXIT
    CALL RANDOM_NUMBER(RandVal)
    IF (SFdepoLayersPartNum(iLayer).GT.0.) THEN
      layerParts=INT(SFdepoLayersPartNum(iLayer)+RandVal)
    ELSE
      layerParts=0
    END IF
    DO iPart=1,layerParts
      SELECT CASE (TRIM(SFdepoLayersSpace(iLayer)))
      CASE('cuboid')
        CALL RANDOM_NUMBER(RandVal2)
        layerPartPos = SFdepoLayersGeo(iLayer,1,:) + &
          (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
      CASE('cylinder')
        PartRadius = SFdepoLayersRadius(iLayer) + 1
        DO WHILE(PartRadius.GT.SFdepoLayersRadius(iLayer))
          CALL RANDOM_NUMBER(RandVal2)
          RandVal2 = RandVal2 * 2. - 1.
          layerPartPos = SFdepoLayersGeo(iLayer,1,:) + SFdepoLayersRadius(iLayer) * &
            (RandVal2(1)*SFdepoLayersBaseVector(iLayer,1,:) + RandVal2(2)*SFdepoLayersBaseVector(iLayer,2,:))
          PartRadius = SQRT( (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) * &
            (layerPartPos(1)-SFdepoLayersGeo(iLayer,1,1)) + &
            (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) * &
            (layerPartPos(2)-SFdepoLayersGeo(iLayer,1,2)) + &
            (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) * &
            (layerPartPos(3)-SFdepoLayersGeo(iLayer,1,3)) )
        END DO
      CASE DEFAULT
        CALL abort(__STAMP__, &
          ' Wrong Space for SFdepoLayer: only cuboid and cylinder implemented!')
      END SELECT
      CALL RANDOM_NUMBER(RandVal)
      layerPartPos = layerPartPos + RandVal*SFdepoLayersGeo(iLayer,2,:)
      IF ( SFdepoLayersUseFixBounds(iLayer) ) THEN
        DoCycle=.FALSE.
        DO iSFfix=1,NbrOfSFdepoFixes
          SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
            + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
            + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
          IF (SFfixDistance .GT. 0.) THEN !outside of plane
            DoCycle=.TRUE.
            EXIT
          END IF
        END DO
        IF (DoCycle) CYCLE
      ELSE IF ( SFdepoLayersBounds(iLayer,1,1).GT.layerPartPos(1) .OR. layerPartPos(1).GT.SFdepoLayersBounds(iLayer,2,1) .OR. &
        SFdepoLayersBounds(iLayer,1,2).GT.layerPartPos(2) .OR. layerPartPos(2).GT.SFdepoLayersBounds(iLayer,2,2) .OR. &
        SFdepoLayersBounds(iLayer,1,3).GT.layerPartPos(3) .OR. layerPartPos(3).GT.SFdepoLayersBounds(iLayer,2,3) ) THEN
        CYCLE !outside of bounds
      END IF
      CALL calcSfSource(1 &
        ,Species(SFdepoLayersSpec(iLayer))%ChargeIC*SFdepoLayersMPF(iLayer)*w_sf &
        ,Vec1,Vec2,Vec3,layerPartPos,-iPart,const_opt=ConstantSFdepoLayers)
    END DO ! iPart
    IF (iLayer.EQ.NbrOfSFdepoLayers .AND. ConstantSFdepoLayers) THEN
      SFdepoLayersAlreadyDone=.TRUE.
    END IF
  END DO ! iLayer=1,NbrOfSFdepoLayers

  !--SFResampleAnalyzeSurfCollis
  IF (SFResampleAnalyzeSurfCollis) THEN
    iPart=0
    DO iPart2=1,LastAnalyzeSurfCollis%PartNumberDepo
      !get random (equal!) position between [1,PartNumberSamp]
      CALL RANDOM_NUMBER(RandVal)
      iPart=MIN(1+INT(RandVal*REAL(LastAnalyzeSurfCollis%PartNumberSamp)),LastAnalyzeSurfCollis%PartNumberSamp)
      !perform surfaceflux-like push into sf-layer outside of mesh
      CALL RANDOM_NUMBER(RandVal)
      FractPush = RandVal*LastAnalyzeSurfCollis%pushTimeStep*LastAnalyzeSurfCollis%WallState(4:6,iPart)
      IF ( DOT_PRODUCT(LastAnalyzeSurfCollis%NormVecOfWall,FractPush).LE.r_SF  ) THEN
        layerPartPos = LastAnalyzeSurfCollis%WallState(1:3,iPart) + FractPush
        IF ( LastAnalyzeSurfCollis%UseFixBounds ) THEN
          DoCycle=.FALSE.
          DO iSFfix=1,NbrOfSFdepoFixes
            SFfixDistance = SFdepoFixesGeo(iSFfix,2,1)*(layerPartPos(1)-SFdepoFixesGeo(iSFfix,1,1)) &
              + SFdepoFixesGeo(iSFfix,2,2)*(layerPartPos(2)-SFdepoFixesGeo(iSFfix,1,2)) &
              + SFdepoFixesGeo(iSFfix,2,3)*(layerPartPos(3)-SFdepoFixesGeo(iSFfix,1,3))
            IF (SFfixDistance .GT. 0.) THEN !outside of plane
              DoCycle=.TRUE.
              EXIT
            END IF
          END DO
          IF (DoCycle) CYCLE
        ELSE IF ( LastAnalyzeSurfCollis%Bounds(1,1).GT.layerPartPos(1) .OR. &
          layerPartPos(1).GT.LastAnalyzeSurfCollis%Bounds(2,1) .OR. &
          LastAnalyzeSurfCollis%Bounds(1,2).GT.layerPartPos(2) .OR. &
          layerPartPos(2).GT.LastAnalyzeSurfCollis%Bounds(2,2) .OR. &
          LastAnalyzeSurfCollis%Bounds(1,3).GT.layerPartPos(3) .OR. &
          layerPartPos(3).GT.LastAnalyzeSurfCollis%Bounds(2,3) ) THEN
          CYCLE !outside of bounds
        END IF
      ELSE
        CYCLE !outside of r_SF
      END IF
      CALL calcSfSource(4 &
        ,Species(LastAnalyzeSurfCollis%Species(iPart))%ChargeIC &
        *Species(LastAnalyzeSurfCollis%Species(iPart))%MacroParticleFactor*w_sf &
        ,Vec1,Vec2,Vec3,layerPartPos,-iPart2,PartVelo=LastAnalyzeSurfCollis%WallState(4:6,iPart))
    END DO ! iPart2
  END IF !SFResampleAnalyzeSurfCollis

  !-- external particles
#if USE_MPI
  IF (usevMPF) THEN
    DO iPart=1,NbrOfextParticles  !external Particles
      CALL calcSfSource(4,Species(ExtPartSpecies(iPart))%ChargeIC*ExtPartMPF(iPart)*w_sf &
        ,Vec1,Vec2,Vec3,ExtPartState(1:3,iPart),-iPart,PartVelo=ExtPartState(4:6,iPart))
    END DO
  ELSE
    DO iPart=1,NbrOfextParticles  !external Particles
      CALL calcSfSource(4 &
        ,Species(ExtPartSpecies(iPart))%ChargeIC*Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf &
        ,Vec1,Vec2,Vec3,ExtPartState(1:3,iPart),-iPart,PartVelo=ExtPartState(4:6,iPart))
    END DO
  END IF ! usevMPF
  ! deallocate external state
  SDEALLOCATE(ExtPartState)
  SDEALLOCATE(ExtPartSpecies)
  SDEALLOCATE(ExtPartToFIBGM)
  SDEALLOCATE(ExtPartMPF)
  NbrOfExtParticles=0
#endif /*USE_MPI*/

  !-- add const. PartSource and relaxation (only once, i.e., during call with .NOT.DoInnerParts)
  IF (PartSourceConstExists) THEN
    DO iElem = 1, nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
            PartSource(1:4,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem) + PartSourceConst(1:4,mm,ll,kk,iElem)
            IF (RelaxDeposition) THEN
#if ((USE_HDG) && (PP_nVar==1))
              PartSource(4,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem) * RelaxFac*dtWeight &
                                           + PartSourceOld(1,1,mm,ll,kk,iElem) * (1.0-RelaxFac*dtWeight)
              PartSourceOld(1,1,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem)
#else
              PartSource(1:4,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem) * RelaxFac*dtWeight &
                                             + PartSourceOld(1:4,1,mm,ll,kk,iElem) * (1.0-RelaxFac*dtWeight)
              PartSourceOld(1:4,1,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem)
#endif
            END IF
          END DO !mm
        END DO !ll
      END DO !kk
    END DO !iElem
  ELSE IF (RelaxDeposition) THEN
    DO iElem = 1, nElems
      DO kk = 0, PP_N
        DO ll = 0, PP_N
          DO mm = 0, PP_N
#if ((USE_HDG) && (PP_nVar==1))
            PartSource(4,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem) * RelaxFac*dtWeight &
                                         + PartSourceOld(1,1,mm,ll,kk,iElem) * (1.0-RelaxFac*dtWeight)
            PartSourceOld(1,1,mm,ll,kk,iElem) = PartSource(4,mm,ll,kk,iElem)
#else
            PartSource(1:4,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem) * RelaxFac*dtWeight &
                                           + PartSourceOld(1:4,1,mm,ll,kk,iElem) * (1.0-RelaxFac*dtWeight)
            PartSourceOld(1:4,1,mm,ll,kk,iElem) = PartSource(1:4,mm,ll,kk,iElem)
#endif
          END DO !mm
        END DO !ll
      END DO !kk
    END DO !iElem
  END IF !PartSourceConstExists
END IF !.NOT.DoInnerParts

IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
  ! map PartSource from Equististant points on Gauss-Points
  DO iElem=1,PP_nElems
    CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
  END DO ! iElem=1,PP_nElems
END IF

END SUBROUTINE DepositionMethod_SF


SUBROUTINE DepositionMethod_SF1D(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N,PP_nElems
USE MOD_Particle_Vars      ,ONLY: Species, PartSpecies,PDM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
USE MOD_PICDepo_Vars       ,ONLY: PartSource,w_sf,Vdm_EquiN_GaussN,sf1d_dir,r_sf,r2_sf_inv,r2_sf,ElemDepo_xGP,DoSFEqui,alpha_sf
USE MOD_Part_Tools         ,ONLY: isDepositParticle
#if USE_MPI
USE MOD_Particle_MPI_Vars  ,ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM,ExtPartMPF,NbrOfextParticles
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Particle_Mesh_Vars ,ONLY: NbrOfCases,casematrix
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
REAL               :: ShiftedPart(1:3)
LOGICAL            :: chargedone(1:nElems)!, SAVE_GAUSS
REAL               :: radius2, S, S1, Fac(1:4)!, Fac2(4)
INTEGER            :: firstPart,LastPart
INTEGER            :: kk, ll, mm, ppp
INTEGER            :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER            :: k,l,m,iPart,ind,iCase,iElem,expo,ElemID
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
Vec1(1:3) = 0.
Vec2(1:3) = 0.
Vec3(1:3) = 0.
IF (GEO%nPeriodicVectors.EQ.1) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
END IF
IF (GEO%nPeriodicVectors.EQ.2) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
END IF
IF (GEO%nPeriodicVectors.EQ.3) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
END IF
DO iPart=firstPart,LastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  ! Set charge pre-factor
  IF (usevMPF) THEN
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
  ELSE
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
  END IF ! usevMPF
  !IF(fac(4).GT.0.) print*,'charge pos'
  Fac(1:3) = PartState(4:6,iPart)*Fac(4)
  !-- determine which background mesh cells (and interpolation points within) need to be considered
  chargedone(:) = .FALSE.
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      ShiftedPart(ind) = PartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
           casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    IF(sf1d_dir.EQ.1)THEN
      kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = GEO%FIBGMjmax
      lmin = GEO%FIBGMjmin
      mmax = GEO%FIBGMkmax
      mmin = GEO%FIBGMkmin
    ELSEIF(sf1d_dir.EQ.2)THEN
      kmax = GEO%FIBGMimax
      kmin = GEO%FIBGMimin
      lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = GEO%FIBGMkmax
      mmin = GEO%FIBGMkmin
    ELSE
      kmax = GEO%FIBGMimax
      kmin = GEO%FIBGMimin
      lmax = GEO%FIBGMjmax
      lmin = GEO%FIBGMjmin
      mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
    END IF
    !-- go through all these cells
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
                radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                        * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                !-- calculate charge and current density at ip point using a shape function
                !-- currently only one shapefunction available, more to follow (including structure change)
                IF (radius2 .LT. r2_sf) THEN
                  S = 1. - r2_sf_inv * radius2
                !radius2=GaussDistance(k,l,m)
                !IF (radius2 .LT. 1.0) THEN
                !  S = 1 -  radius2
                  S1 = S*S
                  DO expo = 3, alpha_sf
                    S1 = S*S1
                  END DO
                  PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                  PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                END IF
              END DO; END DO; END DO
              chargedone(ElemID) = .TRUE.
            END IF
          END DO ! ppp
        END DO ! mm
      END DO ! ll
    END DO ! kk
  END DO ! iCase (periodicity)
END DO ! i
#if USE_MPI
IF(.NOT.DoInnerParts)THEN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (NbrOfextParticles .GT. 0) THEN
    IF (GEO%nPeriodicVectors.EQ.1) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    END IF
    IF (GEO%nPeriodicVectors.EQ.2) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    END IF
    IF (GEO%nPeriodicVectors.EQ.3) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
    END IF
  END IF

  DO iPart=1,NbrOfextParticles  !external Particles
    ! Set charge pre-factor
    IF (usevMPF) THEN
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
    ELSE
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
    END IF ! usevMPF
    Fac(1:3) = ExtPartState(4:6,iPart)*Fac(4)
    chargedone(:) = .FALSE.
    !-- determine which background mesh cells (and interpolation points within) need to be considered
    DO iCase = 1, NbrOfCases
      DO ind = 1,3
        ShiftedPart(ind) = ExtPartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
             casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
      END DO
      IF(sf1d_dir.EQ.1)THEN
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        lmax = GEO%FIBGMjmax
        lmin = GEO%FIBGMjmin
        mmax = GEO%FIBGMkmax
        mmin = GEO%FIBGMkmin
      ELSEIF(sf1d_dir.EQ.2)THEN
        kmax = GEO%FIBGMimax
        kmin = GEO%FIBGMimin
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        mmax = GEO%FIBGMkmax
        mmin = GEO%FIBGMkmin
      ELSE
        kmax = GEO%FIBGMimax
        kmin = GEO%FIBGMimin
        lmax = GEO%FIBGMjmax
        lmin = GEO%FIBGMjmin
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
      END IF
      !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain)
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
                    radius2 = (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID)) &
                            * (ShiftedPart(sf1d_dir) - ElemDepo_xGP(sf1d_dir,k,l,m,ElemID))
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. r2_sf) THEN
                      S = 1. - r2_sf_inv * radius2
                      !IF(S.LT.0.) print*,'dist neg '
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END DO
  END DO
  ! deallocate external state
  SDEALLOCATE(ExtPartState)
  SDEALLOCATE(ExtPartSpecies)
  SDEALLOCATE(ExtPartToFIBGM)
  SDEALLOCATE(ExtPartMPF)
  NbrOfExtParticles=0
END IF
#endif /*USE_MPI*/

IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
  ! map source from Equististant points on Gauss-Points
  DO iElem=1,PP_nElems
    CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
  END DO ! iElem=1,PP_nElems
END IF

END SUBROUTINE DepositionMethod_SF1D


SUBROUTINE DepositionMethod_SF2D(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc                     ,ONLY: PP_N,PP_nElems
USE MOD_Particle_Vars               ,ONLY: Species, PartSpecies,PDM,usevMPF,PartMPF
USE MOD_Particle_Vars               ,ONLY: PartState
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO
USE MOD_PICDepo_Vars                ,ONLY: PartSource,w_sf,Vdm_EquiN_GaussN
USE MOD_PICDepo_Vars                ,ONLY: sf1d_dir,r_sf,r2_sf_inv,r2_sf,ElemDepo_xGP,DoSFEqui,alpha_sf,w_sf
USE MOD_PICDepo_Shapefunction_Tools ,ONLY: DepoSFParticleLocally
#if USE_MPI
USE MOD_Particle_MPI_Vars           ,ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM,ExtPartMPF,NbrOfextParticles
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars            ,ONLY: nDeposPerElem
USE MOD_LoadBalance_Timers          ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers          ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars          ,ONLY: GEO,casematrix, NbrOfCases
USE MOD_Mesh_Vars                   ,ONLY: nElems
USE MOD_Part_Tools                  ,ONLY: isDepositParticle
USE MOD_ChangeBasis                 ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
LOGICAL            :: chargedone(1:nElems)!, SAVE_GAUSS
REAL               :: radius2, S, S1, Fac(1:4)!, Fac2(4)
INTEGER            :: ElemID, iCase, ind
REAL               :: ShiftedPart(1:3)
INTEGER            :: firstPart,LastPart
INTEGER            :: kk, ll, mm, ppp
INTEGER            :: kmin, kmax, lmin, lmax, mmin, mmax
INTEGER            :: k,l,m,I,J,iPart,iElem
INTEGER            :: expo
REAL               :: dx,dy
LOGICAL            :: DepoLoc
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
Vec1(1:3) = 0.
Vec2(1:3) = 0.
Vec3(1:3) = 0.
IF (GEO%nPeriodicVectors.EQ.1) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
END IF
IF (GEO%nPeriodicVectors.EQ.2) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
END IF
IF (GEO%nPeriodicVectors.EQ.3) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
END IF
DO iPart=firstPart,LastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  ! Set charge pre-factor
  IF (usevMPF) THEN
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf
  ELSE
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf
  END IF ! usevMPF
  !IF(fac(4).GT.0.) print*,'charge pos'
  Fac(1:3) = PartState(4:6,iPart)*Fac(4)
  !-- determine which background mesh cells (and interpolation points within) need to be considered
  chargedone(:) = .FALSE.
  DO iCase = 1, NbrOfCases
    DO ind = 1,3
      ShiftedPart(ind) = PartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
           casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    IF(sf1d_dir.EQ.1)THEN
      ! x
      kmax = GEO%FIBGMimax
      kmin = GEO%FIBGMimin
      ! y
      lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      ! z
      mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      ! Set main directions
      I = 2
      J = 3
    ELSEIF(sf1d_dir.EQ.2)THEN
      ! x
      kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      ! y
      lmax = GEO%FIBGMjmax
      lmin = GEO%FIBGMjmin
      ! z
      mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      ! Set main directions
      I = 1
      J = 3
    ELSE
      ! x
      kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      ! y
      lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      ! z
      mmax = GEO%FIBGMkmax
      mmin = GEO%FIBGMkmin
      ! Set main directions
      I = 1
      J = 2
    END IF
    !-- go through all these cells
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
              ! Check whether the SF particle has to be locally deposited (set DepoLoc=T/F)
              CALL DepoSFParticleLocally(DepoLoc,ElemID,iPart)

              ! Shape function deposition
              IF(.NOT.DepoLoc)THEN
                !--- go through all gauss points
                DO m=0,PP_N; DO l=0,PP_N; DO k=0,PP_N
                  !-- calculate distance between gauss and particle
                  dX = ABS(ShiftedPart(I) - ElemDepo_xGP(I,k,l,m,ElemID))
                  IF(dX.GT.r_sf) CYCLE
                  dY = ABS(ShiftedPart(J) - ElemDepo_xGP(J,k,l,m,ElemID))
                  IF(dY.GT.r_sf) CYCLE
                  radius2 = dX*dX+dY*dY
                  !-- calculate charge and current density at ip point using a shape function
                  !-- currently only one shapefunction available, more to follow (including structure change)
                  IF (radius2 .LT. r2_sf) THEN
                    S = 1. - r2_sf_inv * radius2
                  !radius2=GaussDistance(k,l,m)
                  !IF (radius2 .LT. 1.0) THEN
                  !  S = 1 -  radius2
                    S1 = S*S
                    DO expo = 3, alpha_sf
                      S1 = S*S1
                    END DO
                    PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                    PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                  END IF
                END DO; END DO; END DO
              END IF ! DepoLoc


              chargedone(ElemID) = .TRUE.
            END IF
          END DO ! ppp
        END DO ! mm
      END DO ! ll
    END DO ! kk
  END DO ! iCase (periodicity)
END DO ! i
#if USE_MPI
IF(.NOT.DoInnerParts)THEN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (NbrOfextParticles .GT. 0) THEN
    IF (GEO%nPeriodicVectors.EQ.1) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    END IF
    IF (GEO%nPeriodicVectors.EQ.2) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    END IF
    IF (GEO%nPeriodicVectors.EQ.3) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
    END IF
  END IF

  DO iPart=1,NbrOfextParticles  !external Particles
    ! Set charge pre-factor
    IF (usevMPF) THEN
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf
    ELSE
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf
    END IF ! usevMPF
    Fac(1:3) = ExtPartState(4:6,iPart)*Fac(4)
    chargedone(:) = .FALSE.
    !-- determine which background mesh cells (and interpolation points within) need to be considered
    DO iCase = 1, NbrOfCases
      DO ind = 1,3
        ShiftedPart(ind) = ExtPartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
             casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
      END DO
      IF(sf1d_dir.EQ.1)THEN
        ! x
        kmax = GEO%FIBGMimax
        kmin = GEO%FIBGMimin
        ! y
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        ! z
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        ! Set main directions
        I = 2
        J = 3
      ELSEIF(sf1d_dir.EQ.2)THEN
        ! x
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        ! y
        lmax = GEO%FIBGMjmax
        lmin = GEO%FIBGMjmin
        ! z
        mmax = CEILING((ShiftedPart(3)+r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
        mmax = MIN(mmax,GEO%FIBGMkmax)
        mmin = FLOOR((ShiftedPart(3)-r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
        mmin = MAX(mmin,GEO%FIBGMkmin)
        ! Set main directions
        I = 1
        J = 3
      ELSE
        ! x
        kmax = CEILING((ShiftedPart(1)+r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
        kmax = MIN(kmax,GEO%FIBGMimax)
        kmin = FLOOR((ShiftedPart(1)-r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
        kmin = MAX(kmin,GEO%FIBGMimin)
        ! y
        lmax = CEILING((ShiftedPart(2)+r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
        lmax = MIN(lmax,GEO%FIBGMjmax)
        lmin = FLOOR((ShiftedPart(2)-r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
        lmin = MAX(lmin,GEO%FIBGMjmin)
        ! z
        mmax = GEO%FIBGMkmax
        mmin = GEO%FIBGMkmin
        ! Set main directions
        I = 1
        J = 2
      END IF
      !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain)
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
                    dX = ABS(ShiftedPart(I) - ElemDepo_xGP(I,k,l,m,ElemID))
                    IF(dX.GT.r_sf) CYCLE
                    dY = ABS(ShiftedPart(J) - ElemDepo_xGP(J,k,l,m,ElemID))
                    IF(dY.GT.r_sf) CYCLE
                    radius2 = dX*dX+dY*dY
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. r2_sf) THEN
                      S = 1. - r2_sf_inv * radius2
                      !IF(S.LT.0.) print*,'dist neg '
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END DO
  END DO
  ! deallocate external state
  SDEALLOCATE(ExtPartState)
  SDEALLOCATE(ExtPartSpecies)
  SDEALLOCATE(ExtPartToFIBGM)
  SDEALLOCATE(ExtPartMPF)
  NbrOfExtParticles=0
END IF
#endif /*USE_MPI*/

IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
  ! map source from Equististant points on Gauss-Points
  DO iElem=1,PP_nElems
    CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
  END DO ! iElem=1,PP_nElems
END IF

END SUBROUTINE DepositionMethod_SF2D


SUBROUTINE DepositionMethod_SFCS(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N,PP_nElems
USE MOD_Globals_Vars       ,ONLY: PI
USE MOD_Particle_Vars      ,ONLY: Species, PartSpecies,PDM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_PICDepo_Vars       ,ONLY: PartSource,Vdm_EquiN_GaussN,SfRadiusInt,r_sf_scale,r_sf0
USE MOD_PICDepo_Vars       ,ONLY: ElemDepo_xGP,DoSFEqui,alpha_sf,w_sf
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: nDeposPerElem
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime,LBElemPauseTime,LBElemSplitTime,LBElemPauseTime_avg
USE MOD_LoadBalance_Timers ,ONLY: LBElemSplitTime_avg
#endif /*USE_LOADBALANCE*/
USE MOD_Particle_Mesh_Vars ,ONLY: GEO,casematrix, NbrOfCases
#if USE_MPI
USE MOD_Particle_MPI_Vars  ,ONLY: ExtPartState,ExtPartSpecies,ExtPartToFIBGM,ExtPartMPF,NbrOfextParticles
#endif /*USE_MPI*/
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Part_Tools         ,ONLY: isDepositParticle
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Vec1(1:3), Vec2(1:3), Vec3(1:3)
REAL               :: radius2, S, S1, Fac(1:4)!, Fac2(4)
LOGICAL            :: chargedone(1:nElems)!, SAVE_GAUSS
INTEGER            :: ElemID, iCase, ind
REAL               :: ShiftedPart(1:3)
INTEGER            :: firstPart,LastPart
INTEGER            :: kmin, kmax, lmin, lmax, mmin, mmax
REAL               :: dx,dy,dz
INTEGER            :: ppp
INTEGER            :: k,l,m,iPart,iElem
REAL               :: local_r_sf, local_r2_sf, local_r2_sf_inv
INTEGER            :: expo
INTEGER            :: kk,ll,mm
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
Vec1(1:3) = 0.
Vec2(1:3) = 0.
Vec3(1:3) = 0.
IF (GEO%nPeriodicVectors.EQ.1) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
END IF
IF (GEO%nPeriodicVectors.EQ.2) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
END IF
IF (GEO%nPeriodicVectors.EQ.3) THEN
  Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
  Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
  Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
END IF
DO iPart=firstPart,LastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  ! compute local radius
  local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(1:SfRadiusInt,iPart),PartState(1:SfRadiusInt,iPart)))
  local_r2_sf=local_r_sf*local_r_sf
  local_r2_sf_inv=1./local_r2_sf
  IF (usevMPF) THEN
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
  ELSE
    Fac(4)= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
  END IF ! usevMPF
  !IF(fac(4).GT.0.) print*,'charge pos'
  Fac(1:3) = PartState(4:6,iPart)*Fac(4)
  !-- determine which background mesh cells (and interpolation points within) need to be considered
  DO iCase = 1, NbrOfCases
    chargedone(:) = .FALSE.
    DO ind = 1,3
      ShiftedPart(ind) = PartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
           casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
    END DO
    kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
    kmax = MIN(kmax,GEO%FIBGMimax)
    kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
    kmin = MAX(kmin,GEO%FIBGMimin)
    lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
    lmax = MIN(lmax,GEO%FIBGMjmax)
    lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
    lmin = MAX(lmin,GEO%FIBGMjmin)
    mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
    mmax = MIN(mmax,GEO%FIBGMkmax)
    mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
    mmin = MAX(mmin,GEO%FIBGMkmin)
    !-- go through all these cells
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
                dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                IF(dX.GT.local_r_sf) CYCLE
                dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                IF(dY.GT.local_r_sf) CYCLE
                dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                IF(dZ.GT.local_r_sf) CYCLE
                radius2 = dX*dX+dY*dY+dZ*dZ
                !-- calculate charge and current density at ip point using a shape function
                !-- currently only one shapefunction available, more to follow (including structure change)
                IF (radius2 .LT. local_r2_sf) THEN
                  S = 1. - local_r2_sf_inv * radius2
                !radius2=GaussDistance(k,l,m)
                !IF (radius2 .LT. 1.0) THEN
                !  S = 1 -  radius2
                  S1 = S*S
                  DO expo = 3, alpha_sf
                    S1 = S*S1
                  END DO
                  PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                  PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                END IF
              END DO; END DO; END DO
              chargedone(ElemID) = .TRUE.
            END IF
          END DO ! ppp
        END DO ! mm
      END DO ! ll
    END DO ! kk
  END DO ! iCase (periodicity)
END DO ! i
#if USE_MPI
IF(.NOT.DoInnerParts)THEN
  Vec1(1:3) = 0.
  Vec2(1:3) = 0.
  Vec3(1:3) = 0.
  IF (NbrOfextParticles .GT. 0) THEN
    IF (GEO%nPeriodicVectors.EQ.1) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
    END IF
    IF (GEO%nPeriodicVectors.EQ.2) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
    END IF
    IF (GEO%nPeriodicVectors.EQ.3) THEN
      Vec1(1:3) = GEO%PeriodicVectors(1:3,1)
      Vec2(1:3) = GEO%PeriodicVectors(1:3,2)
      Vec3(1:3) = GEO%PeriodicVectors(1:3,3)
    END IF
  END IF

  DO iPart=1,NbrOfextParticles  !external Particles
    local_r_sf= r_sf0 * (1.0 + r_sf_scale*DOT_PRODUCT(PartState(1:SfRadiusInt,iPart),PartState(1:SfRadiusInt,iPart)))
    local_r2_sf=local_r_sf*local_r_sf
    local_r2_sf_inv=1./local_r2_sf
    IF (usevMPF) THEN
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC * ExtPartMPF(iPart)*w_sf/(PI*(local_r_sf**3))
    ELSE
      Fac(4)= Species(ExtPartSpecies(iPart))%ChargeIC &
            * Species(ExtPartSpecies(iPart))%MacroParticleFactor*w_sf/(PI*(local_r_sf**3))
    END IF ! usevMPF
    Fac(1:3) = ExtPartState(4:6,iPart)*Fac(4)
    !-- determine which background mesh cells (and interpolation points within) need to be considered
    DO iCase = 1, NbrOfCases
      chargedone(:) = .FALSE.
      DO ind = 1,3
        ShiftedPart(ind) = ExtPartState(ind,iPart) + casematrix(iCase,1)*Vec1(ind) + &
             casematrix(iCase,2)*Vec2(ind) + casematrix(iCase,3)*Vec3(ind)
      END DO
      kmax = CEILING((ShiftedPart(1)+local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1))
      kmax = MIN(kmax,GEO%FIBGMimax)
      kmin = FLOOR((ShiftedPart(1)-local_r_sf-GEO%xminglob)/GEO%FIBGMdeltas(1)+1)
      kmin = MAX(kmin,GEO%FIBGMimin)
      lmax = CEILING((ShiftedPart(2)+local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2))
      lmax = MIN(lmax,GEO%FIBGMjmax)
      lmin = FLOOR((ShiftedPart(2)-local_r_sf-GEO%yminglob)/GEO%FIBGMdeltas(2)+1)
      lmin = MAX(lmin,GEO%FIBGMjmin)
      mmax = CEILING((ShiftedPart(3)+local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3))
      mmax = MIN(mmax,GEO%FIBGMkmax)
      mmin = FLOOR((ShiftedPart(3)-local_r_sf-GEO%zminglob)/GEO%FIBGMdeltas(3)+1)
      mmin = MAX(mmin,GEO%FIBGMkmin)
      !-- go through all these cells (should go through non if periodic and shiftedpart not in my domain
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
                    dX = ABS(ShiftedPart(1) - ElemDepo_xGP(1,k,l,m,ElemID))
                    IF(dX.GT.local_r_sf) CYCLE
                    dY = ABS(ShiftedPart(2) - ElemDepo_xGP(2,k,l,m,ElemID))
                    IF(dY.GT.local_r_sf) CYCLE
                    dZ = ABS(ShiftedPart(3) - ElemDepo_xGP(3,k,l,m,ElemID))
                    IF(dZ.GT.local_r_sf) CYCLE
                    radius2 = dX*dX+dY*dY+dZ*dZ
                    !-- calculate charge and current density at ip point using a shape function
                    !-- currently only one shapefunction available, more to follow (including structure change)
                    IF (radius2 .LT. local_r2_sf) THEN
                      S = 1. - local_r2_sf_inv * radius2
                      !IF(S.LT.0.) print*,'dist neg '
                    !radius2=GaussDistance(k,l,m)
                    !IF (radius2 .LT. 1.0) THEN
                    !  S = 1 -  radius2
                      S1 = S*S
                      DO expo = 3, alpha_sf
                        S1 = S*S1
                      END DO
                      PartSource(1:3,k,l,m,ElemID) = PartSource(1:3,k,l,m,ElemID) + Fac(1:3) * S1
                      PartSource( 4 ,k,l,m,ElemID) = PartSource( 4 ,k,l,m,ElemID) + Fac(4) * S1
                    END IF
                END DO; END DO; END DO
                chargedone(ElemID) = .TRUE.
              END IF
            END DO ! ppp
          END DO ! mm
        END DO ! ll
      END DO ! kk
    END DO
  END DO
  ! deallocate external state
  SDEALLOCATE(ExtPartState)
  SDEALLOCATE(ExtPartSpecies)
  SDEALLOCATE(ExtPartToFIBGM)
  SDEALLOCATE(ExtPartMPF)
  NbrOfExtParticles=0
END IF
#endif /*USE_MPI*/
IF( .NOT.DoInnerParts .AND. DoSFEqui) THEN
  ! map PartSource from Equististant points on Gauss-Points
  DO iElem=1,PP_nElems
    CALL ChangeBasis3D(4,PP_N,PP_N,Vdm_EquiN_GaussN,PartSource(:,:,:,:,iElem),PartSource(:,:,:,:,iElem))
  END DO ! iElem=1,PP_nElems
END IF

END SUBROUTINE DepositionMethod_SFCS


SUBROUTINE DepositionMethod_DD(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc                ,ONLY: PP_nElems,PP_N
USE MOD_Particle_Vars          ,ONLY: Species, PartSpecies,PDM,PEM,usevMPF,PartMPF
USE MOD_Particle_Vars          ,ONLY: PartState,PartPosRef
USE MOD_PICDepo_Vars           ,ONLY: PartSource,DeltaType,XiNDepo,wBaryNDepo,Vdm_NDepo_GaussN,NKnots,Knots,NDepoChooseK,NDepo
USE MOD_PICDepo_Vars           ,ONLY: DDMassInv,DeltaDistriChangeBasis
USE MOD_PICDepo_Tools          ,ONLY: DeBoorRef,DeBoor
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Part_Tools             ,ONLY: isDepositParticle
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping
USE MOD_Basis                  ,ONLY: LagrangeInterpolationPolys,BernSteinPolynomial
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers     ,ONLY: LBStartTime,LBSplitTime,LBPauseTime,LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3,0:NDepo) :: L_xi
INTEGER                   :: firstPart,LastPart
#if USE_LOADBALANCE
REAL                      :: tLBStart
#endif /*USE_LOADBALANCE*/
REAL                      :: DeltaIntCoeff,prefac
INTEGER                   :: iElem,i,j,k,iPart
!===================================================================================================================================
! TODO: Info why and under which conditions the following 'RETURN' is called
IF((DoInnerParts).AND.(LastPart.LT.firstPart)) RETURN
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem=1,PP_nElems
  DO iPart=firstPart,LastPart
    ! TODO: Info why and under which conditions the following 'CYCLE' is called
    IF(doPartInExists)THEN
      IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
    ELSE
      IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
    END IF
    IF(PEM%Element(iPart).EQ.iElem)THEN
      ! Don't deposit neutral particles!
      IF(.NOT.isDepositParticle(iPart)) CYCLE
      ! Set pre-factor
      IF (usevMPF) THEN
        prefac= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
      ELSE
        prefac= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
      END IF ! usevMPF
      ! Map Particle to -1|1 space (re-used in interpolation)
      IF(.NOT.DoRefMapping)THEN
        CALL GetPositionInRefElem(PartState(1:3,iPart),PartPosRef(1:3,iPart),iElem)
      END IF
      ! get value of test function at particle position
      SELECT CASE(DeltaType)
      CASE(1)
        ! xi   -direction
        CALL LagrangeInterpolationPolys(PartPosRef(1,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(1,:))
        ! eta  -direction
        CALL LagrangeInterpolationPolys(PartPosRef(2,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(2,:))
        ! zeta -direction
        CALL LagrangeInterpolationPolys(PartPosRef(3,iPart),NDepo,XiNDepo,wBaryNDepo,L_xi(3,:))
      CASE(2)
        DO i=0,NDepo
          ! xi   -direction
           CALL BernsteinPolynomial(NDepo,i,PartPosRef(1,iPart),L_xi(1,i),NDepoChooseK)
          ! eta  -direction
           CALL BernsteinPolynomial(NDepo,i,PartPosRef(2,iPart),L_xi(2,i),NDepoChooseK)
          ! zeta  -direction
           CALL BernsteinPolynomial(NDepo,i,PartPosRef(3,iPart),L_xi(3,i),NDepoChooseK)
        END DO ! i
      CASE(3)
        ! xi - direction
        CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(1,iPart),L_xi(1,:))
        ! eta - direction
        CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(2,iPart),L_xi(2,:))
        ! zeta - direction
        CALL DeBoorRef(NDepo,NKnots,Knots,PartPosRef(3,iPart),L_xi(3,:))
      END SELECT
      DO k=0,NDepo
        DO j=0,NDepo
          DO i=0,NDepo
       !     print*,'i,j,k,L',i,j,k,L_xi(1,i)* L_xi(2,j)* L_xi(3,k)
            DeltaIntCoeff = L_xi(1,i)* L_xi(2,j)* L_xi(3,k)*prefac
            PartSource(1:3,i,j,k,iElem) = PartSource(1:3,i,j,k,iElem) + DeltaIntCoeff*PartState(4:6,iPart)
            PartSource( 4 ,i,j,k,iElem) = PartSource( 4 ,i,j,k,iElem) + DeltaIntCoeff
          END DO ! i
        END DO ! j
      END DO ! k
    END IF ! Particle in Element
  END DO ! ParticleVecLength
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO ! iElem
IF(.NOT.DoInnerParts)THEN
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  DO iElem=1,PP_nElems
    DO k=0,NDepo
      DO j=0,NDepo
        DO i=0,NDepo
          PartSource( : ,i,j,k,iElem) = PartSource( : ,i,j,k,iElem) * DDMassInv(i,j,k,iElem)
        END DO ! i
      END DO ! j
    END DO ! k
    IF(DeltaDistriChangeBasis)THEN
      CALL ChangeBasis3D(4,NDepo,PP_N,Vdm_NDepo_GaussN,PartSource(1:4,0:NDepo,0:NDepo,0:NDepo,iElem)&
                                                      ,PartSource(1:4,0:PP_N ,0:PP_N ,0:PP_N, iElem))
    END IF
#if USE_LOADBALANCE
    CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
  END DO ! loop over all elems
END IF ! DoInnerParts
END SUBROUTINE DepositionMethod_DD


SUBROUTINE DepositionMethod_MVW(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N
USE MOD_Particle_Vars      ,ONLY: Species, PartSpecies,PDM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_PICDepo_Vars       ,ONLY: PartSource,bgmdeltas,BGMSource,GaussBGMFactor,GaussBGMIndex,BGMVolume
#if USE_MPI
USE MOD_PICDepo_MPI        ,ONLY: MPISourceExchangeBGM
#else /*NOT USE_MPI*/
USE MOD_PICDepo_MPI        ,ONLY: PeriodicSourceExchange
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
#endif /*USE_MPI*/
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Part_Tools         ,ONLY: isDepositParticle
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime,LBElemSplitTime
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_LOADBALANCE
REAL                        :: tLBStart
#endif /*USE_LOADBALANCE*/
REAL                        :: Charge, TSource(1:4)
INTEGER                     :: kk,ll,mm,k,l,m
INTEGER                     :: firstPart,LastPart,iPart,iElem,i
REAL                        :: alpha1,alpha2,alpha3
!===================================================================================================================================
! Step 1: Deposition of all particles onto background mesh -> densities
! IF(DoInnerParts) BGMSource=0.0 ! not possible due to periodic stuff --> two communications
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
BGMSource(:,:,:,:) = 0.0
DO iPart = firstPart, lastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
  IF (usevMPF) THEN
    Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
  ELSE
    Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF ! usevMPF
  !Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
  k = FLOOR(PartState(1,iPart)/BGMdeltas(1))
  l = FLOOR(PartState(2,iPart)/BGMdeltas(2))
  m = FLOOR(PartState(3,iPart)/BGMdeltas(3))
  alpha1 = (PartState(1,iPart) / BGMdeltas(1)) - k
  alpha2 = (PartState(2,iPart) / BGMdeltas(2)) - l
  alpha3 = (PartState(3,iPart) / BGMdeltas(3)) - m
  TSource(:) = 0.0
!#if (PP_nVar==8)
  TSource(1) = PartState(4,iPart)*Charge
  TSource(2) = PartState(5,iPart)*Charge
  TSource(3) = PartState(6,iPart)*Charge
!#endif
  TSource(4) = Charge

  BGMSource(k,l,m,1:4)       = BGMSource(k,l,m,1:4)       + (TSource * (1.-alpha1)*(1.-alpha2)*(1.-alpha3))
  BGMSource(k,l,m+1,1:4)     = BGMSource(k,l,m+1,1:4)     + (TSource * (1.-alpha1)*(1.-alpha2)*(alpha3))
  BGMSource(k,l+1,m,1:4)     = BGMSource(k,l+1,m,1:4)     + (TSource * (1.-alpha1)*(alpha2)*(1.-alpha3))
  BGMSource(k,l+1,m+1,1:4)   = BGMSource(k,l+1,m+1,1:4)   + (TSource * (1.-alpha1)*(alpha2)*(alpha3))
  BGMSource(k+1,l,m,1:4)     = BGMSource(k+1,l,m,1:4)     + (TSource * (alpha1)*(1.-alpha2)*(1.-alpha3))
  BGMSource(k+1,l,m+1,1:4)   = BGMSource(k+1,l,m+1,1:4)   + (TSource * (alpha1)*(1.-alpha2)*(alpha3))
  BGMSource(k+1,l+1,m,1:4)   = BGMSource(k+1,l+1,m,1:4)   + (TSource * (alpha1)*(alpha2)*(1.-alpha3))
  BGMSource(k+1,l+1,m+1,1:4) = BGMSource(k+1,l+1,m+1,1:4) + (TSource * (alpha1)*(alpha2)*(alpha3))
END DO
BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_LOADBALANCE
CALL LBPauseTime(LB_CARTMESHDEPO,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
! should be treated in this way, unforunately, we would negelct the periodic stuff
!IF(.NOT.DoInnerParts)
CALL MPISourceExchangeBGM()
#else
IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif

! Step 2: Interpolation of densities onto grid
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem = 1, nElems
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
       k = GaussBGMIndex(1,kk,ll,mm,iElem)
       l = GaussBGMIndex(2,kk,ll,mm,iElem)
       m = GaussBGMIndex(3,kk,ll,mm,iElem)
       alpha1 = GaussBGMFactor(1,kk,ll,mm,iElem)
       alpha2 = GaussBGMFactor(2,kk,ll,mm,iElem)
       alpha3 = GaussBGMFactor(3,kk,ll,mm,iElem)
!#if (PP_nVar==8)
       DO i = 1,3
         PartSource(i,kk,ll,mm,iElem) = PartSource(i,kk,ll,mm,iElem)            + &
              BGMSource(k,l,m,i) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k,l,m+1,i) * (1-alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k,l+1,m,i) * (1-alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k,l+1,m+1,i) * (1-alpha1) * (alpha2) * (alpha3) + &
              BGMSource(k+1,l,m,i) * (alpha1) * (1-alpha2) * (1-alpha3) + &
              BGMSource(k+1,l,m+1,i) * (alpha1) * (1-alpha2) * (alpha3) + &
              BGMSource(k+1,l+1,m,i) * (alpha1) * (alpha2) * (1-alpha3) + &
              BGMSource(k+1,l+1,m+1,i) * (alpha1) * (alpha2) * (alpha3)
       END DO
!#endif
         PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem)          + &
            BGMSource(k,l,m,4) * (1-alpha1) * (1-alpha2) * (1-alpha3) + &
            BGMSource(k,l,m+1,4) * (1-alpha1) * (1-alpha2) * (alpha3) + &
            BGMSource(k,l+1,m,4) * (1-alpha1) * (alpha2) * (1-alpha3) + &
            BGMSource(k,l+1,m+1,4) * (1-alpha1) * (alpha2) * (alpha3) + &
            BGMSource(k+1,l,m,4) * (alpha1) * (1-alpha2) * (1-alpha3) + &
            BGMSource(k+1,l,m+1,4) * (alpha1) * (1-alpha2) * (alpha3) + &
            BGMSource(k+1,l+1,m,4) * (alpha1) * (alpha2) * (1-alpha3) + &
            BGMSource(k+1,l+1,m+1,4) * (alpha1) * (alpha2) * (alpha3)
     END DO !mm
   END DO !ll
 END DO !kk
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO !iElem
!DEALLOCATE(BGMSource)

! Suppress compiler warning
RETURN
IF(DoInnerParts)k=0
END SUBROUTINE DepositionMethod_MVW


SUBROUTINE DepositionMethod_MS(DoInnerParts,doPartInExists,doParticle_In)
!===================================================================================================================================
! ''
! Deposits 
!===================================================================================================================================
! MODULES
USE MOD_PreProc            ,ONLY: PP_N
USE MOD_Particle_Vars      ,ONLY: Species, PartSpecies,PDM,usevMPF,PartMPF
USE MOD_Particle_Vars      ,ONLY: PartState
USE MOD_PICDepo_Vars       ,ONLY: PartSource,BGMdeltas,GaussBGMIndex,PartSource,GPWeight,BGMSource,BGMVolume
#if USE_MPI
USE MOD_PICDepo_MPI        ,ONLY: MPISourceExchangeBGM
#else /*NOT USE_MPI*/
USE MOD_PICDepo_MPI        ,ONLY: PeriodicSourceExchange
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
#endif /*USE_MPI*/
USE MOD_Part_Tools         ,ONLY: isDepositParticle
USE MOD_Mesh_Vars          ,ONLY: nElems
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime,LBElemSplitTime
#endif /*USE_LOADBALANCE*/
USE MOD_PICDepo_Tools      ,ONLY: DeBoor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)          :: DoInnerParts ! TRUE: do cell-local particles, FALSE: do external particles from other (MPI) cells
LOGICAL,INTENT(IN),OPTIONAL :: doParticle_In(1:PDM%ParticleVecLength) ! TODO: definition of this variable
LOGICAL,INTENT(IN),OPTIONAL :: doPartInExists
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: auxiliary(0:3),weight(1:3,0:3), locweight
INTEGER            :: PosInd(3),r,ss,t,u,v,w, dir, weightrun
#if USE_LOADBALANCE
REAL               :: tLBStart
#endif /*USE_LOADBALANCE*/
INTEGER            :: k,l,m,kk,ll,mm,iPart,iElem
INTEGER            :: firstPart,LastPart
REAL               :: Charge
!===================================================================================================================================
! Step 1: Deposition of all particles onto background mesh -> densities
!ALLOCATE(BGMSource(BGMminX:BGMmaxX,BGMminY:BGMmaxY,BGMminZ:BGMmaxZ,1:4))
! IF(DoInnerParts) BGMSource=0. not possible due to periodic stuff
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
BGMSource(:,:,:,:) = 0.0
DO iPart = firstPart, lastPart
  ! TODO: Info why and under which conditions the following 'CYCLE' is called
  IF(doPartInExists)THEN
    IF (.NOT.(PDM%ParticleInside(iPart).AND.doParticle_In(iPart))) CYCLE
  ELSE
    IF (.NOT.PDM%ParticleInside(iPart)) CYCLE
  END IF
  ! Don't deposit neutral particles!
  IF(.NOT.isDepositParticle(iPart)) CYCLE
!      Charge = Species(PartSpecies(iPart))%ChargeIC*Species(PartSpecies(iPart))%MacroParticleFactor
  IF (usevMPF) THEN
    Charge= Species(PartSpecies(iPart))%ChargeIC * PartMPF(iPart)
  ELSE
    Charge= Species(PartSpecies(iPart))%ChargeIC * Species(PartSpecies(iPart))%MacroParticleFactor
  END IF ! usevMPF
  PosInd(1) = FLOOR(PartState(1,iPart)/BGMdeltas(1))
  PosInd(2) = FLOOR(PartState(2,iPart)/BGMdeltas(2))
  PosInd(3) = FLOOR(PartState(3,iPart)/BGMdeltas(3))
  !print*,'posind(1:3),charge',posInd,charge
  DO dir = 1,3               ! x,y,z direction
    DO weightrun = 0,3
      DO mm = 0, 3
        IF (mm.EQ.weightrun) then
          auxiliary(mm) = 1.0
        ELSE
          auxiliary(mm) = 0.0
        END IF
      END DO
      CALL DeBoor(PosInd(dir),auxiliary,PartState(dir,iPart),weight(dir,weightrun),dir)
    END DO
  END DO
  DO k = PosInd(1)-1, PosInd(1)+2
    kk = abs(k - PosInd(1) - 2)
    DO l = PosInd(2)-1, PosInd(2)+2
      ll = abs(l - PosInd(2) - 2)
      DO m = PosInd(3)-1, PosInd(3)+2
        mm = abs(m - PosInd(3) - 2)
        locweight = weight(1,kk)*weight(2,ll)*weight(3,mm)*charge
!#if (PP_nVar==8)
        BGMSource(k,l,m,1) = BGMSource(k,l,m,1) + PartState(4,iPart)* locweight
        BGMSource(k,l,m,2) = BGMSource(k,l,m,2) + PartState(5,iPart)* locweight
        BGMSource(k,l,m,3) = BGMSource(k,l,m,3) + PartState(6,iPart)* locweight
!#endif
        BGMSource(k,l,m,4) = BGMSource(k,l,m,4) + locweight
     !   print*,'BMGSOURCE4',BGMSOURCE(k,l,m,4)
      END DO
    END DO
  END DO
END DO
BGMSource(:,:,:,:) = BGMSource(:,:,:,:) / BGMVolume

#if USE_LOADBALANCE
CALL LBPauseTime(LB_CARTMESHDEPO,tLBStart)
#endif /*USE_LOADBALANCE*/
#if USE_MPI
!IF(.NOT.DoInnerParts)THEN has to be communicated each time :(
CALL MPISourceExchangeBGM()
#else
IF (GEO%nPeriodicVectors.GT.0) CALL PeriodicSourceExchange()
#endif

! Step 2: Interpolation of densities onto grid
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
DO iElem = 1, nElems
  DO kk = 0, PP_N
    DO ll = 0, PP_N
      DO mm = 0, PP_N
        k = GaussBGMIndex(1,kk,ll,mm,iElem)
        l = GaussBGMIndex(2,kk,ll,mm,iElem)
        m = GaussBGMIndex(3,kk,ll,mm,iElem)
        DO r = k-1,k+2
          u = r-k+2
          DO ss = l-1,l+2
            v = ss-l+2
            DO t = m-1,m+2
              w = t-m+2
!#if (PP_nVar==8)
              PartSource(1:4,kk,ll,mm,iElem) = PartSource(1:4,kk,ll,mm,iElem) + BGMSource(r,ss,t,1:4) * GPWeight(iElem,kk,ll,mm,u,v,w)
              !DO i = 1,3
              !  PartSource(i,kk,ll,mm,iElem) = PartSource(i,kk,ll,mm,iElem) + BGMSource(r,ss,t,i) * GPWeight(iElem,kk,ll,mm,u,v,w)
              !END DO
!#endif
              !PartSource(4,kk,ll,mm,iElem) = PartSource(4,kk,ll,mm,iElem) + BGMSource(r,ss,t,4) * GPWeight(iElem,kk,ll,mm,u,v,w)
            END DO !t
          END DO !s
        END DO !r
      END DO !mm
    END DO !ll
  END DO !kk
#if USE_LOADBALANCE
  CALL LBElemSplitTime(iElem,tLBStart)
#endif /*USE_LOADBALANCE*/
END DO !iElem
!DEALLOCATE(BGMSource)

! Suppress compiler warning
RETURN
IF(DoInnerParts)k=0
END SUBROUTINE DepositionMethod_MS


END MODULE MOD_PICDepo_Method
