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


MODULE MOD_PML
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
INTERFACE InitPML
  MODULE PROCEDURE InitPML
END INTERFACE
INTERFACE FinalizePML
  MODULE PROCEDURE FinalizePML
END INTERFACE
INTERFACE PMLTimeRamping
  MODULE PROCEDURE PMLTimeRamping
END INTERFACE
INTERFACE CalcPMLSource
  MODULE PROCEDURE CalcPMLSource
END INTERFACE
INTERFACE PMLTimeDerivative
  MODULE PROCEDURE PMLTimeDerivative
END INTERFACE

PUBLIC::InitPML,FinalizePML,PMLTimeRamping,CalcPMLSource,PMLTimeDerivative
!===================================================================================================================================
PUBLIC::DefineParametersPML
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

CONTAINS

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!==================================================================================================================================
!> Define parameters for PML
!==================================================================================================================================
SUBROUTINE DefineParametersPML()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Perfectly Matched Layer (PML) Region")

CALL prms%CreateLogicalOption(  'DoPML'            , 'Use a PML region' , '.FALSE.')
CALL prms%CreateRealOption(     'PMLzeta0'         , 'Damping parameter maximum' , '0.')
CALL prms%CreateRealOption(     'PMLalpha0'        , 'Complex frequency shift parameter' , '0.')
CALL prms%CreateIntOption(      'PMLzetaShape'     , 'Parameter for selecting the shape of the damping function'&
    //'0=constant, 1=linear, 2=sinusoidal, 3=polynomial of degree 4', '0')
CALL prms%CreateRealOption(     'PMLRampLength'    , 'Parameter for defining the relative length of the ramping region'//&
    ' in the PML. Ranges from 0.0 to 1.0' , '1.')
CALL prms%CreateIntOption(      'PMLspread'        , 'Use damping parameter spreading, which distributes the values to all DOF'&
    , '0')
CALL prms%CreateIntOption(      'PMLwriteFields'   , 'Write the PML region information to .h5 file (CURRENTLY DEPRICATED)'&
    , '0')
CALL prms%CreateLogicalOption(  'PMLzetaNorm'      , 'Use damping parameter normalization (CURRENTLY DEPRICATED)' , '.FALSE.')

CALL prms%CreateLogicalOption(  'DoPMLTimeRamp'    , 'Use damping parameter ramping over time' , '.FALSE.')
CALL prms%CreateRealOption(     'PMLTimeRamptStart', 'Time to start the damping parameter ramping over time' , '-1.')
CALL prms%CreateRealOption(     'PMLTimeRamptEnd'  , 'Time to stop the damping parameter ramping over time' , '-1.')

CALL prms%CreateRealArrayOption('xyzPhysicalMinMax'     , '[xmin, xmax, ymin, ymax, zmin, zmax] vector for defining a'&
    //' PML region by giving the bounding box coordinates of the PHYSICAL region', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealArrayOption('xyzPMLzetaShapeOrigin' , 'The origin for defining a rising of declining slope'&
    //' of the PML region damping parameters', '0.0 , 0.0 , 0.0')
CALL prms%CreateRealArrayOption('xyzPMLMinMax'          , '[xmin, xmax, ymin, ymax, zmin, zmax] vector for defining a '&
    //' PML region by giving the bounding box coordinates of the PML region', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')

END SUBROUTINE DefineParametersPML

SUBROUTINE InitPML()
!===================================================================================================================================
! Initialize perfectly matched layer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_PML_Vars
USE MOD_HDF5_output       ,ONLY: GatheredWriteArray,WriteAttributeToHDF5,WriteHDF5Header
USE MOD_HDF5_Output_Fields,ONLY: WritePMLzetaGlobalToHDF5
USE MOD_Interfaces        ,ONLY: FindInterfacesInRegion,FindElementInRegion,CountAndCreateMappings,DisplayRanges,SelectMinMaxRegion
USE MOD_IO_HDF5           ,ONLY: AddToElemData,ElementOut
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars  ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT PML...'
!===================================================================================================================================
! Readin
!===================================================================================================================================
DoPML                  = GETLOGICAL('DoPML','.FALSE.')
IF(.NOT.DoPML) THEN
  LBWRITE(UNIT_stdOut,'(A)') ' PML region deactivated. '
  PMLnVar=0
  nPMLElems=0
  RETURN
ELSE
#if PP_nVar == 1 || PP_nVar == 4
  CALL abort(__STAMP__,'Equation system does not support a PML!')
#endif
  PMLnVar=24
END IF

! Get PML parameters
PMLzeta0               = GETREAL('PMLzeta0','0.')
PMLalpha0              = GETREAL('PMLalpha0','0.')
PMLzetaShape           = GETINT('PMLzetaShape','0')
PMLRampLength          = GETREAL('PMLRampLength','1.')
PMLspread              = GETINT('PMLspread','0')
PMLwriteFields         = GETINT('PMLwriteFields','0')
PMLzetaNorm            = GETLOGICAL('PMLzetaNorm','.FALSE.')

! Temporal PML ramping
DoPMLTimeRamp          = GETLOGICAL('DoPMLTimeRamp','.FALSE.')
PMLTimeRamptStart      = GETREAL('PMLTimeRamptStart','-1.')
PMLTimeRamptEnd        = GETREAL('PMLTimeRamptEnd','-1.')
PMLsDeltaT             = 0.0 ! init
PMLTimeRampCoeff       = 0.0 ! init
IF(ANY((/PMLTimeRamptStart,PMLTimeRamptEnd/).LT.0.0))THEN
  PMLTimeRamptStart    = 0.0
  PMLTimeRamptEnd      = 0.0
  DoPMLTimeRamp        = .FALSE.
ELSE
  IF(ALMOSTEQUALRELATIVE(PMLTimeRamptStart,PMLTimeRamptEnd,1E-3))THEN
    LBWRITE(UNIT_stdOut,'(A)') ' WARNING: PML time ramp uses the same times for tStart and tEnd. Relative difference is < 1E-3'
    PMLsDeltaT         = 1e12 ! set no a very high value
  ELSE
    IF(PMLTimeRamptStart.GT.PMLTimeRamptEnd)THEN
      CALL abort(__STAMP__,' PMLTimeRamptStart must be smaller than PMLTimeRamptEnd.')
    END IF
    PMLsDeltaT         = 1/(PMLTimeRamptEnd-PMLTimeRamptStart)
    PMLTimeRampCoeff   = -PMLTimeRamptStart * PMLsDeltaT
  END IF
END IF
PMLTimeRamp            = 1.0 ! init

! determine PML elements
xyzPhysicalMinMax(1:6)     = GETREALARRAY('xyzPhysicalMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
xyzPMLzetaShapeOrigin(1:3) = GETREALARRAY('xyzPMLzetaShapeOrigin',3,'0.0,0.0,0.0')
xyzPMLMinMax(1:6)          = GETREALARRAY('xyzPMLMinMax',6,'0.0,0.0,0.0,0.0,0.0,0.0')
! use xyzPhysicalMinMax before xyzPMLMinMax:
! 1.) check for xyzPhysicalMinMax
! 2.) check for xyzPMLMinMax
CALL SelectMinMaxRegion('PML',usePMLMinMax,&
                        'xyzPhysicalMinMax',xyzPhysicalMinMax,&
                        'xyzPMLMinMax',xyzPMLMinMax)

! display ranges of PML region depending on usePMLMinMax
!CALL DisplayRanges('usePMLMinMax',usePMLMinMax,'xyzPMLMinMax',xyzPMLMinMax(1:6),'xyzPhysicalMinMax',xyzPhysicalMinMax(1:6))

! find all elements in the PML region
IF(usePMLMinMax)THEN ! find all elements located inside of 'xyzPMLMinMax'
  CALL FindElementInRegion(isPMLElem,xyzPMLMinMax,&
                           ElementIsInside=.TRUE. ,DoRadius=.FALSE.,Radius=-1.,DisplayInfo=.TRUE.)
ELSE ! find all elements located outside of 'xyzPhysicalMinMax'
  CALL FindElementInRegion(isPMLElem,xyzPhysicalMinMax,&
                           ElementIsInside=.FALSE.,DoRadius=.FALSE.,Radius=-1.,DisplayInfo=.TRUE.)
END IF

! find all faces in the PML region
CALL FindInterfacesInRegion(isPMLFace,isPMLInterFace,isPMLElem,info_opt='find all faces in the PML region')

! Get number of PML Elems, Faces and Interfaces. Create Mappngs PML <-> physical region
CALL CountAndCreateMappings('PML',&
                            isPMLElem,isPMLFace,isPMLInterFace,&
                            nPMLElems,nPMLFaces, nPMLInterFaces,&
                            ElemToPML,PMLToElem,&
                            FaceToPML,PMLToFace,&
                            FaceToPMLInter,PMLInterToFace,&
                            DisplayInfo=.TRUE.)

! nPMLElems is determined, now allocate the PML field correnspondingly
!ALLOCATE(U2       (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))  =>  U_N(iElem)%U2
!ALLOCATE(U2t      (1:PMLnVar,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))  =>  U_N(iElem)%U2t
!U2 =0.  =>  U_N(iElem)%U2
!U2t=0.  =>  U_N(iElem)%U2t

CALL AddToElemData(ElementOut,'PMLElem',LogArray=isPMLElem(:))

! Set the damping profile function zeta=f(x,y) in the PML region
CALL SetPMLdampingProfile()

! create a HDF5 file containing the PMLzetaGlobal field
CALL WritePMLzetaGlobalToHDF5()

PMLInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT PML DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitPML



PPURE SUBROUTINE PMLTimeRamping(t,RampingFactor)
!===================================================================================================================================
! set the scaling factor which ramps the damping factor over time
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_PML_Vars,      ONLY: PMLTimeRamptStart,PMLTimeRamptEnd,PMLsDeltaT,PMLTimeRampCoeff,DoPMLTimeRamp
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: RampingFactor
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER             :: i,j,k,iPMLElem,m
!===================================================================================================================================
IF(.NOT.DoPMLTimeRamp) RETURN

! Ramp if t in [0,1]
IF(t.LT.PMLTimeRamptStart)THEN
  RampingFactor = 0.0                             ! set PMLTimeRamp to 0.0
ELSEIF(t.GT.PMLTimeRamptEnd)THEN
  RampingFactor = 1.0                             ! set PMLTimeRamp to 1.0
ELSE
  RampingFactor = PMLsDeltaT*t + PMLTimeRampCoeff ! set PMLTimeRamp to [0,1]
END IF
END SUBROUTINE PMLTimeRamping


SUBROUTINE CalcPMLSource()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals  ,ONLY: abort
USE MOD_DG_Vars  ,ONLY: U_N
USE MOD_PML_Vars ,ONLY: nPMLElems,PMLToElem,PML,PMLTimeRamp
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iPMLElem,m,iElem
!===================================================================================================================================
! sources for the standard variables
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  iElem = PMLToElem(iPMLElem)
  DO m=1,8
    U_N(iElem)%Ut(m,i,j,k) = U_N(iElem)%Ut(m,i,j,k) - PMLTimeRamp*(&
        PML(iPMLElem)%zeta(1,i,j,k)*U_N(iElem)%U2(m*3-2,i,j,k) +&   ! = 1,4,7,10,13,16,19,22
        PML(iPMLElem)%zeta(2,i,j,k)*U_N(iElem)%U2(m*3-1,i,j,k) +&   ! = 2,5,8,11,12,17,20,23
        PML(iPMLElem)%zeta(3,i,j,k)*U_N(iElem)%U2(m*3  ,i,j,k) )    ! = 3,6,9,12,15,18,21,24
  END DO
END DO; END DO; END DO !nPMLElems,k,j,i
END DO
END SUBROUTINE CalcPMLSource


SUBROUTINE PMLTimeDerivative()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY: U_N
USE MOD_PML_Vars      ,ONLY: nPMLElems,PMLToElem,PMLnVar
USE MOD_PML_Vars      ,ONLY: PML,PMLTimeRamp
USE MOD_Mesh_Vars     ,ONLY: N_VolMesh
USE MOD_Equation_Vars ,ONLY: fDamping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iPMLElem,iPMLVar,iElem
!===================================================================================================================================
! We have to take the inverse of the Jacobians into account
! the '-' sign is due to the movement of the term to the right-hand-side of the equation
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  iElem = PMLToElem(iPMLElem)
  DO iPMLVar=1,PMLnVar
    U_N(iElem)%U2t(iPMLVar,i,j,k) = - U_N(iElem)%U2t(iPMLVar,i,j,k) * N_VolMesh(iElem)%sJ(i,j,k)
  END DO
END DO; END DO; END DO; END DO !nPMLElems,k,j,i


! Add Source Terms
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  iElem = PMLToElem(iPMLElem)
  U_N(iElem)%U2t(1 : 3,i,j,k) = U_N(iElem)%U2t(1 : 3,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(1 : 3,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(4 : 6,i,j,k) = U_N(iElem)%U2t(4 : 6,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(4 : 6,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(7 : 9,i,j,k) = U_N(iElem)%U2t(7 : 9,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(7 : 9,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(10:12,i,j,k) = U_N(iElem)%U2t(10:12,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(10:12,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(13:15,i,j,k) = U_N(iElem)%U2t(13:15,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(13:15,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(16:18,i,j,k) = U_N(iElem)%U2t(16:18,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(16:18,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(19:21,i,j,k) = U_N(iElem)%U2t(19:21,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(19:21,i,j,k) * PMLTimeRamp
  U_N(iElem)%U2t(22:24,i,j,k) = U_N(iElem)%U2t(22:24,i,j,k) - PML(iPMLElem)%zetaEff(1:3,i,j,k) * U_N(iElem)%U2(22:24,i,j,k) * PMLTimeRamp
END DO; END DO; END DO; END DO !nPMLElems,k,j,i


! 1.) DEBUGPML: apply the damping factor also to PML source terms
! copied from: U(7:8,i,j,k,iElem) = U(7:8,i,j,k,iElem) * fDamping
!U2 = U2 * fDamping

! 2.) DEBUGPML: apply the damping factor only to PML variables for Phi_E and Phi_B
!               to prevent charge-related instabilities (accumulation of divergence compensation over time)

DO iPMLElem=1,nPMLElems
  iElem = PMLToElem(iPMLElem)
  ASSOCIATE( U2 => U_N(iElem)%U2(:,:,:,:) )
    U2(19:24,:,:,:) = fDamping* U2(19:24,:,:,:)
  END ASSOCIATE
END DO

END SUBROUTINE PMLTimeDerivative


SUBROUTINE SetPMLdampingProfile()
!===================================================================================================================================
! Determine the local PML damping factor in x,y and z-direction using a constant/linear/polynomial/... function
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh      ,ONLY: GetMeshMinMaxBoundaries
USE MOD_Mesh_Vars ,ONLY: N_VolMesh,xyzMinMax,offSetElem
USE MOD_PML_Vars  ,ONLY: PML,usePMLMinMax,xyzPMLzetaShapeOrigin,xyzPMLMinMax
USE MOD_PML_Vars  ,ONLY: nPMLElems,PMLToElem
USE MOD_PML_Vars  ,ONLY: PMLzeta0,PMLalpha0,xyzPhysicalMinMax,PMLzetaShape
USE MOD_DG_Vars   ,ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iPMLElem,iElem,Nloc
REAL                :: XiN(3)
REAL                :: fFuncType
INTEGER             :: iDir,PMLDir
REAL                :: xMin,xMax
LOGICAL             :: SetValue
!===================================================================================================================================
!ALLOCATE(PMLRamp          (0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
!ALLOCATE(PMLzeta      (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
!ALLOCATE(PMLzetaEff   (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
!ALLOCATE(PMLalpha     (1:3,0:PP_N,0:PP_N,0:PP_N,1:nPMLElems))
ALLOCATE(PML(1:nPMLElems))
DO iPMLElem=1,nPMLElems
  iElem = PMLToElem(iPMLElem)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(PML(iPMLElem)%zeta(1:3,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(PML(iPMLElem)%zetaEff(1:3,0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(PML(iPMLElem)%alpha(1:3,0:Nloc,0:Nloc,0:Nloc))
  PML(iPMLElem)%zeta = 0.
  PML(iPMLElem)%zetaEff = 0.
  PML(iPMLElem)%alpha = PMLalpha0 ! currently only constant a alpha distribution in the PML region is used
END DO
!PMLzeta=0.
!PMLRamp=1. ! goes from 1 to 0
!PMLzetaEff=0.
!PMLalpha=PMLalpha0 ! currently only constant a alpha distribution in the PML region is used

! get xyzMinMax
CALL GetMeshMinMaxBoundaries()

!determine PMLzeta values for each interpolation point according to ramping function (const., linear, sinusoidal, polynomial)
IF(usePMLMinMax)THEN ! use xyPMLMinMax -> define the PML region
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! CURRENTLY ONLY SUPPORTS ONE DIRECTION, EITHER X- Y- or Z-DIRECTION
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! the PML ramp is oriented via "xyzPMLzetaShapeOrigin" which is used as the origin for defining a rising of declining slope
  ! of the PML ramping profile: below are two examples with a linear profile where the PML region is defined by the same values
  ! the only differens is that the origin "xyzPMLzetaShapeOrigin" located at "P(x_PML)" in the two examples is moved in the domain
  ! to a larger x-value (in this example)
  !
  !                     /                        \                                              !
  !   example 1        /                          \       example 2                             !
  !                   /                            \                                            !
  !                  /                              \                                           !
  !     origin      /  PML                      PML  \    origin                                !
  !                /  ramp                      ramp  \                                         !
  !    P(x_PML)   /                                    \    P(x_PML)                            !
  !              /                                      \                                       !
  !  ___________/                                        \___________                           !
  !                                                                                             !
  !  ---------------------->                     ---------------------->                        !
  !      x_PML              x                                 x_PML      x                      !
  !                                                                                             !
  ! --------------------------------------------------------------------------------------------------------------------------------
  DO iDir=1,3 !1=x, 2=y, 3=z
    IF((xyzPMLzetaShapeOrigin(iDir).GT.xyzPMLMinMax(2*iDir-1)).AND.(xyzPMLzetaShapeOrigin(iDir).LT.xyzPMLMinMax(2*iDir)))THEN
      IF(iDir.LT.3) CYCLE ! if all directions are true, the the point must be indisde the region
    ELSE
      PMLDir=iDir
      EXIT ! if one direction is outside, the point must be outside of the region
    END IF
    SWRITE(UNIT_stdOut,'(ES25.14E3,ES25.14E3,ES25.14E3)') xyzPMLzetaShapeOrigin(1),xyzPMLzetaShapeOrigin(2),xyzPMLzetaShapeOrigin(3)
    CALL abort(__STAMP__,'The origin reference point "xyzPMLzetaShapeOrigin" cannot lie within the PML region defined by "xyzPMLMinMax"')
  END DO

  ! set new values for minimum and maximum to the domain boundary values
  xyzPMLMinMax(2*PMLDir-1) = MAX(xyzPMLMinMax(2*PMLDir-1),xyzMinMax(2*PMLDir-1)) ! minimum
  xyzPMLMinMax(2*PMLDir  ) = MIN(xyzPMLMinMax(2*PMLDir  ),xyzMinMax(2*PMLDir  )) ! maximum
  SWRITE(UNIT_stdOut,'(A,I2)') 'Setting xyzPMLMinMax to <=xyzMinMax for iDir=',PMLDir
  DO iPMLElem=1,nPMLElems; DO k=1,PP_N+1; DO j=1,PP_N+1; DO i=1,PP_N+1
    iElem = PMLToElem(iPMLElem)
    SetValue=.FALSE.
    ASSOCIATE( Elem_xGP => N_VolMesh(iElem)%Elem_xGP(:,:,:,:) )
      IF((Elem_xGP(PMLDir,i,j,k).GE.xyzPMLMinMax(2*PMLDir-1)).AND.&
         (Elem_xGP(PMLDir,i,j,k).LE.xyzPMLMinMax(2*PMLDir)))THEN ! Point is in [PMLDir]-direction region
        xMin = xyzPMLMinMax(2*PMLDir-1)-xyzPMLzetaShapeOrigin(PMLDir) ! Min of region defined for PML region
        xMax = xyzPMLMinMax(2*PMLDir  )-xyzPMLzetaShapeOrigin(PMLDir) ! Max of region defined for PML region
        XiN(PMLDir) = ( Elem_xGP(PMLDir,i,j,k) - xyzPMLzetaShapeOrigin(PMLDir)-MIN(xMin,xMax) )/&
                      ( MAX(xMin,xMax)                                        -MIN(xMin,xMax) )
        !PMLzeta(PMLDir,i,j,k,iPMLElem) = PMLzeta0*fFuncType(&
                                       !( Elem_xGP(PMLDir,i,j,k) - xyzPMLzetaShapeOrigin(PMLDir)-MIN(xMin,xMax) )/&
                                       !( MAX(xMin,xMax)                                        -MIN(xMin,xMax) ),&
                                       !PMLzetashape)
        SetValue=.TRUE.
      END IF
    END ASSOCIATE
    IF(SetValue) PML(iPMLElem)%zeta(iDir,i-1,j-1,k-1) = PMLzeta0*fFuncType(XiN(iDir),PMLzetaShape)
  END DO; END DO; END DO; END DO !iPMLElem,k,j,i
! ----------------------------------------------------------------------------------------------------------------------------------
ELSE ! use xyzPhysicalMinMax -> define the physical region
  DO iPMLElem=1,nPMLElems; DO k=1,PP_N+1; DO j=1,PP_N+1; DO i=1,PP_N+1
    iElem = PMLToElem(iPMLElem)
    ASSOCIATE( Elem_xGP => N_VolMesh(iElem)%Elem_xGP(:,:,:,:) )
      DO iDir=1,3 !1=x, 2=y, 3=z
        IF          (Elem_xGP(iDir,i,j,k) .LT.   xyzPhysicalMinMax(2*iDir-1)) THEN ! region is in lower part
          XiN(iDir) = (ABS(Elem_xGP(iDir,i,j,k)) - ABS(xyzPhysicalMinMax(2*iDir-1)))/&   ! of the domain
                      (ABS(xyzMinMax(2*iDir-1))  - ABS(xyzPhysicalMinMax(2*iDir-1)))
          !PMLzeta(iDir,i,j,k,iPMLElem)   = PMLzeta0*fFuncType(XiN,PMLzetaShape)
        ELSEIF      (Elem_xGP(iDir,i,j,k) .GT.   xyzPhysicalMinMax(2*iDir)) THEN ! region is in upper part
          XiN(iDir) = (ABS(Elem_xGP(iDir,i,j,k)) - ABS(xyzPhysicalMinMax(2*iDir)))/&   ! of the domain
                      (ABS(xyzMinMax(2*iDir))    - ABS(xyzPhysicalMinMax(2*iDir)))
          !PMLzeta(iDir,i,j,k,iPMLElem)   = PMLzeta0*fFuncType(XiN,PMLzetaShape)
        END IF
      END DO
    END ASSOCIATE
    DO iDir=1,3 !1=x, 2=y, 3=z
      PML(iPMLElem)%zeta(iDir,i-1,j-1,k-1)   = PMLzeta0*fFuncType(XiN(iDir),PMLzetaShape)
    END DO
  END DO; END DO; END DO; END DO !iElem,k,j,i

END IF ! usePMLMinMax
! ----------------------------------------------------------------------------------------------------------------------------------
! CFS-PML formulation: calculate zeta eff using the complex frequency shift PMLalpha
DO iPMLElem=1,nPMLElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  PML(iPMLElem)%zetaEff(:,i,j,k) = ( PML(iPMLElem)%alpha(:,i,j,k)+PML(iPMLElem)%zeta(:,i,j,k) )
END DO; END DO; END DO; END DO !iPMLElem,k,j,i
!DEALLOCATE(PMLalpha)
DO iPMLElem = 1, nPMLElems
  DEALLOCATE(PML(iPMLElem)%alpha)
END DO ! iPMLElem = 1, nPMLElems

END SUBROUTINE SetPMLdampingProfile


SUBROUTINE FinalizePML()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_PML_Vars,            ONLY: PML
USE MOD_PML_Vars,            ONLY: ElemToPML,PMLToElem,DoPML,isPMLElem,isPMLFace,PMLToFace,FaceToPML
USE MOD_PML_Vars,            ONLY: PMLInterToFace,FaceToPMLInter,isPMLInterFace
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!RETURN
IF(.NOT.DoPML) RETURN
SDEALLOCATE(PML)
SDEALLOCATE(PMLToElem)
SDEALLOCATE(ElemToPML)
SDEALLOCATE(PMLToFace)
SDEALLOCATE(FaceToPML)
!SDEALLOCATE(PMLRamp)
SDEALLOCATE(isPMLElem)
SDEALLOCATE(isPMLFace)
SDEALLOCATE(isPMLInterFace)
SDEALLOCATE(FaceToPMLInter)
SDEALLOCATE(PMLInterToFace)
!SDEALLOCATE(PMLZetaEff)
END SUBROUTINE FinalizePML
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

END MODULE MOD_PML


#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!-----------------------------------------------------------------------------------------------------------------------------------
! local SUBROUTINES and FUNCTIONS
!-----------------------------------------------------------------------------------------------------------------------------------

!===================================================================================================================================
! switch between different types of ramping functions for the calculation of the local zeta damping value field
!===================================================================================================================================
REAL FUNCTION fFuncType(x,PMLzetaShape)
! MODULES
USE MOD_Globals,       ONLY: abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,    INTENT(IN) :: x
INTEGER, INTENT(IN) :: PMLzetaShape ! linear, polynomial, const., sinusoidal ramping function
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: fLinear,fSinus,fPolynomial
!===================================================================================================================================
fFuncType=0. ! Initialize
SELECT CASE (PMLzetaShape)
CASE(0) !Constant Distribution of the Damping Coefficient
  fFuncType=1.
CASE(1) ! Linear Distribution of the Damping Coefficient
  fFuncType=fLinear(x)
CASE(2) ! Sinusoidal  Distribution of the Damping Coefficient
  fFuncType=fSinus(x)
CASE(3) ! polynomial
  fFuncType=fPolynomial(x)
CASE DEFAULT
  CALL abort(__STAMP__,'Shape function for damping coefficient in PML region not specified!')
END SELECT ! PMLzetaShape

END FUNCTION fFuncType


!===================================================================================================================================
! Evaluates a linear function
!===================================================================================================================================
REAL FUNCTION fLinear(x)
! MODULES
USE MOD_PML_Vars,            ONLY: PMLRampLength
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.PMLRampLength) THEN
  x_temp = x/PMLRampLength
  fLinear = x_temp
ELSE
  fLinear = 1.
END IF
END FUNCTION fLinear


!===================================================================================================================================
! Evaluates a sin function
!===================================================================================================================================
REAL FUNCTION fSinus(x)
! MODULES
USE MOD_PML_Vars,            ONLY: PMLRampLength
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.PMLRampLength) THEN
  x_temp = x/PMLRampLength
  fSinus = x_temp-SIN(2*ACOS(-1.)*x_temp)/(2*ACOS(-1.))
ELSE
  fSinus = 1.
END IF
END FUNCTION fSinus


!===================================================================================================================================
! Evaluates the polynomial -3x^4 + 4x^3
!===================================================================================================================================
REAL FUNCTION fPolynomial(x)
! MODULES
USE MOD_PML_Vars,            ONLY: PMLRampLength
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: x
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: x_temp ![0,1] -> [0,1] sinusodial distribution
!===================================================================================================================================
IF (x.LE.PMLRampLength) THEN
  x_temp = x/PMLRampLength
  fPolynomial = -3*x_temp**4+4*x_temp**3
ELSE
  fPolynomial = 1.
END IF
END FUNCTION fPolynomial
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
