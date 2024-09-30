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

MODULE MOD_Interfaces
!===================================================================================================================================
!> Contains the routines to
!> - identify interfaces between special regions, e.g., PML <-> physical region
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitInterfaces
  MODULE PROCEDURE InitInterfaces
END INTERFACE
INTERFACE FindInterfacesInRegion
  MODULE PROCEDURE FindInterfacesInRegion
END INTERFACE
INTERFACE FindElementInRegion
  MODULE PROCEDURE FindElementInRegion
END INTERFACE
INTERFACE CountAndCreateMappings
  MODULE PROCEDURE CountAndCreateMappings
END INTERFACE
INTERFACE FinalizeInterfaces
  MODULE PROCEDURE FinalizeInterfaces
END INTERFACE
INTERFACE DisplayRanges
  MODULE PROCEDURE DisplayRanges
END INTERFACE
INTERFACE SelectMinMaxRegion
  MODULE PROCEDURE SelectMinMaxRegion
END INTERFACE

PUBLIC::InitInterfaces
PUBLIC::FindElementInRegion
PUBLIC::FindInterfacesInRegion
PUBLIC::CountAndCreateMappings
PUBLIC::FinalizeInterfaces
PUBLIC::DisplayRanges
PUBLIC::SelectMinMaxRegion
!===================================================================================================================================
CONTAINS


SUBROUTINE InitInterfaces
!===================================================================================================================================
!> Check every face and set the correct identifier for selecting the corresponding Riemann solver
!> possible connections are (Master <-> Slave direction is important):
!>   - vacuum     <-> vacuum       : RIEMANN_VACUUM            = 0
!>   - PML        <-> vacuum       : RIEMANN_PML               = 1
!>   - PML        <-> PML          : RIEMANN_PML               = 1
!>   - dielectric <-> dielectric   : RIEMANN_DIELECTRIC        = 2
!>   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC    = 3 ! for conservative fluxes (one flux)
!>   - vacuum      -> dielectric   : RIEMANN_VAC2DIELECTRIC    = 4 ! for conservative fluxes (one flux)
!>   - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC_NC = 5 ! for non-conservative fluxes (two fluxes)
!>   - vacuum      -> dielectric   : RIEMANN_VAC2DIELECTRIC_NC = 6 ! for non-conservative fluxes (two fluxes)
!===================================================================================================================================
! MODULES
USE MOD_globals
USE MOD_Mesh_Vars        ,ONLY: nSides,Face_xGP,NGeo,MortarType
#if ! (USE_HDG)
USE MOD_PML_vars         ,ONLY: DoPML,isPMLFace
#endif /*NOT HDG*/
USE MOD_Dielectric_vars  ,ONLY: DoDielectric,isDielectricFace,isDielectricInterFace,isDielectricElem,DielectricFluxNonConserving
USE MOD_Interfaces_Vars  ,ONLY: InterfaceRiemann,InterfacesInitIsDone
USE MOD_Globals          ,ONLY: abort,UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: mpiroot
#endif
USE MOD_Mesh_Vars        ,ONLY: SideToElem
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,ElemID
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT INTERFACES...'
ALLOCATE(InterfaceRiemann(1:nSides))
DO SideID=1,nSides
  InterfaceRiemann(SideID)=-2 ! set default to invalid number: check later
  ! 0.) Sanity: It is forbidden to connect a PML to a dielectric region because it is not implemented!
#if !(USE_HDG) /*pure Maxwell simulations*/
  IF(DoPML.AND.DoDielectric)THEN
    IF(isPMLFace(SideID).AND.isDielectricFace(SideID)) CALL abort(__STAMP__,&
        'It is forbidden to connect a PML to a dielectric region! (Not implemented)')
  END IF

  ! 1.) Check Perfectly Matched Layer
  ! - PML <-> vacuum              : RIEMANN_PML            = 1
  ! - PML <-> PML                 : RIEMANN_PML            = 1
  IF(DoPML) THEN
    IF (isPMLFace(SideID))THEN ! 1.) RiemannPML additionally calculates the 24 fluxes needed for the auxiliary equations
                                 !     (flux-splitting!)
      InterfaceRiemann(SideID)=RIEMANN_PML
      CYCLE ! don't check the following if the flux has already been calculated here -> continue with next side
    END IF
  END IF ! DoPML
#endif /*NOT HDG*/

  ! 2.) Check Dielectric Medium
  ! c), d) - vacuum     <-> vacuum       : RIEMANN_VACUUM         = 0
  ! b)     - dielectric <-> dielectric   : RIEMANN_DIELECTRIC     = 2
  ! a1)    - dielectric  -> vacuum       : RIEMANN_DIELECTRIC2VAC = 3 or 5 (when using non-conservative fluxes)
  ! a2)    - vacuum      -> dielectric   : RIEMANN_VAC2DIELECTRIC = 4 or 6 (when using non-conservative fluxes)
  ! am1)   - vacuum      -> dielectric mortar : RIEMANN_VAC2DIELECTRIC = 4 or 6 (when using non-conservative fluxes)
  ! am2)   - dielectric  -> vacuum mortar     : RIEMANN_DIELECTRIC2VAC = 3 or 5 (when using non-conservative fluxes)
  IF(DoDielectric) THEN
    IF (isDielectricFace(SideID))THEN ! 1.) RiemannDielectric
      IF(isDielectricInterFace(SideID))THEN
        ! a) physical <-> dielectric region: for Riemann solver, select A+ and A- as functions of f(Eps0,Mu0) or f(EpsR,MuR)
        ElemID = SideToElem(S2E_ELEM_ID,SideID) ! get master element ID for checking if it is in a physical or dielectric region
        IF(ElemID.EQ.-1) THEN
          IF(MortarType(1,SideID).EQ.0) THEN
            ! small mortar slave sides have no corresponding master element
            IF(SideToElem(S2E_NB_ELEM_ID,SideID).GT.0) THEN
              IF(isDielectricElem(SideToElem(S2E_NB_ELEM_ID,SideID)))THEN
                ! am1) big elem is PHYSICAL and small slave DIELECTRIC
                IF(DielectricFluxNonConserving)THEN
                  InterfaceRiemann(SideID)=RIEMANN_VAC2DIELECTRIC_NC ! use two different Riemann solvers
                ELSE
                  InterfaceRiemann(SideID)=RIEMANN_VAC2DIELECTRIC ! A+(EpsR,MuR) and A-(Eps0,Mu0)
                END IF
              ELSE
                ! am2) big elem is DIELECTRIC and small slave PHYSICAL
                IF(DielectricFluxNonConserving)THEN ! use one flux (conserving) or two fluxes (non-conserving) at the interface
                  InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC2VAC_NC ! use two different Riemann solvers
                ELSE
                  InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC2VAC ! A+(Eps0,Mu0) and A-(EpsR,MuR)
                END IF
              END IF
            ELSE
              InterfaceRiemann(SideID)=-1
            END IF
          ELSE
            InterfaceRiemann(SideID)=-1
          END IF
        ELSE
          IF(isDielectricElem(ElemID))THEN
            ! a1) master is DIELECTRIC and slave PHYSICAL
            IF(DielectricFluxNonConserving)THEN ! use one flux (conserving) or two fluxes (non-conserving) at the interface
              InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC2VAC_NC ! use two different Riemann solvers
            ELSE
              InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC2VAC ! A+(Eps0,Mu0) and A-(EpsR,MuR)
            END IF
          ELSE
            ! a2) master is PHYSICAL and slave DIELECTRIC
            IF(DielectricFluxNonConserving)THEN
              InterfaceRiemann(SideID)=RIEMANN_VAC2DIELECTRIC_NC ! use two different Riemann solvers
            ELSE
              InterfaceRiemann(SideID)=RIEMANN_VAC2DIELECTRIC ! A+(EpsR,MuR) and A-(Eps0,Mu0)
            END IF
          END IF
        END IF
      ELSE
        ! b) dielectric region <-> dielectric region
        InterfaceRiemann(SideID)=RIEMANN_DIELECTRIC
      END IF
    ELSE
      ! c) no Dielectric, standard flux
      InterfaceRiemann(SideID)=RIEMANN_VACUUM
    END IF ! IF(isDielectricFace(SideID))
  ELSE
    ! d) no Dielectric, standard flux
    InterfaceRiemann(SideID)=RIEMANN_VACUUM
  END IF ! DoDielectric
END DO ! SideID


! Check if all sides have correctly been set
DO SideID=1,nSides
  IF(InterfaceRiemann(SideID).EQ.-2)THEN ! check if the default value remains unchanged
#if !(USE_HDG) /*pure Maxwell simulations*/
    IPWRITE(UNIT_StdOut,*) "DoPML                          = ", DoPML
#endif /*NOT HDG*/
    IPWRITE(UNIT_StdOut,*) "DoDielectric                   = ", DoDielectric
    IPWRITE(UNIT_StdOut,*) "SideID                         = ", SideID
    IPWRITE(UNIT_StdOut,*) "MortarType(1,SideID)           = ", MortarType(1,SideID)
    IPWRITE(UNIT_StdOut,*) "InterfaceRiemann(SideID)       = ", InterfaceRiemann(SideID)
    IPWRITE(UNIT_StdOut,*) "SideToElem(S2E_ELEM_ID,SideID) = ", SideToElem(S2E_ELEM_ID,SideID)
    IPWRITE(UNIT_StdOut,*) "Face_xGP(1:3,0,0,SideID)       = ", Face_xGP(1:3,0,0,SideID)
    IPWRITE(UNIT_StdOut,*) "Face_xGP(1:3 , 0    , NGeo , SideID) =" , Face_xGP(1:3 , 0    , NGeo , SideID)
    IPWRITE(UNIT_StdOut,*) "Face_xGP(1:3 , NGeo , 0    , SideID) =" , Face_xGP(1:3 , NGeo , 0    , SideID)
    IPWRITE(UNIT_StdOut,*) "Face_xGP(1:3 , NGeo , NGeo , SideID) =" , Face_xGP(1:3 , NGeo , NGeo , SideID)
    CALL abort(__STAMP__,'Interface for Riemann solver not correctly determined (vacuum, dielectric, PML)')
  END IF
END DO ! SideID

InterfacesInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT INTERFACES DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitInterfaces


SUBROUTINE FindElementInRegion(isElem,region,ElementIsInside,DoRadius,Radius,DisplayInfo,GeometryName,GeometryAxis)
!===================================================================================================================================
!> Determine whether an element resides within or outside of a special region (e.g. PML or dielectric region)
!> Additionally, a radius can be supplied for determining if an element belongs to a special region or not
!> Note: As soon as only one DOF is not inside/outside of the region, the complete element is excluded
!> Method 1.) check DOF by using a bounding box
!> Method 2.) Additionally check radius (e.g. when creating dielectric regions in form of a half sphere)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals         ,ONLY: abort,UNIT_stdOut
#if USE_MPI
USE MOD_Globals         ,ONLY: mpiroot
#endif /*USE_MPI*/
USE MOD_Mesh_Vars       ,ONLY: Elem_xGP
USE MOD_Dielectric_Vars ,ONLY: DielectricRadiusValueB
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                     :: ElementIsInside ! Check whether and element DOF in inside/outside of a region or radius
REAL,INTENT(IN)                        :: region(1:6)     ! MIN/MAX for x,y,z of bounding box region
LOGICAL,INTENT(IN)                     :: DoRadius        ! Check if DOF is inside/outside of radius
REAL,INTENT(IN)                        :: Radius          ! Check if DOF is inside/outside of radius
LOGICAL,INTENT(IN),OPTIONAL            :: DisplayInfo     ! Output to stdOut with region size info
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: GeometryName    ! Name of special geometry with user-defined coordinates
INTEGER,INTENT(IN),OPTIONAL            :: GeometryAxis    ! Spatial direction for circle/cylinder definition
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT):: isElem(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,m, dim_1, dim_2
REAL                :: r
REAL                :: rInterpolated
!===================================================================================================================================
! Display region information (how elements are checked + min-max region extensions)
IF(PRESENT(DisplayInfo))THEN
  IF(DisplayInfo)THEN
    IF(ElementIsInside)THEN ! Display information regarding the orientation of the element-region-search
      LBWRITE(UNIT_stdOut,'(A)')'  Checking for elements INSIDE region'
    ELSE
      LBWRITE(UNIT_stdOut,'(A)')'  Checking for elements OUTSIDE region'
    END IF
  CALL DisplayMinMax(region) ! Display table with min-max info of region
  END IF
END IF

! set logical vector for each element with default .FALSE.
IF(ALLOCATED(isElem))THEN
  CALL abort(__STAMP__,'Attempting to allocate already allocated LOGICAL,ALLOCATABLE,INTENT(INOUT) variable in FindElementInRegion')
END IF ! ALLOCATED(isElem)
ALLOCATE(isElem(1:PP_nElems))
isElem=.FALSE.

IF(ElementIsInside)THEN ! Display information regarding the orientation of the element-region-search
  isElem(:)=.TRUE.
ELSE
  isElem(:)=.FALSE.
END IF

! ----------------------------------------------------------------------------------------------------------------------------------
! 1.) use standard bounding box region
! ----------------------------------------------------------------------------------------------------------------------------------
! all DOF in an element must be inside the region, if one DOF is outside, the element is excluded
DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  DO m=1,3 ! m=x,y,z
    IF ( (Elem_xGP(m,i,j,k,iElem) .LT. region(2*m-1)) .OR. & ! 1,3,5
         (Elem_xGP(m,i,j,k,iElem) .GT. region(2*m)) ) THEN   ! 2,4,6 ! element is outside
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
    END IF
  END DO
END DO; END DO; END DO; END DO !iElem,k,j,i

! ----------------------------------------------------------------------------------------------------------------------------------
! 2.) Additionally check a radius (e.g. half sphere regions)
! ----------------------------------------------------------------------------------------------------------------------------------
! if option 'DoRadius' is applied, elements are double-checked if they are within a certain radius
IF(DoRadius.AND.Radius.GT.0.0)THEN
  DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    r = SQRT(Elem_xGP(1,i,j,k,iElem)**2+&
             Elem_xGP(2,i,j,k,iElem)**2+&
             Elem_xGP(3,i,j,k,iElem)**2  )
    ! check if r is larger than the supplied value .AND.
    ! if r is not almost equal to the radius
    IF(r.GT.Radius.AND.(.NOT.ALMOSTEQUALRELATIVE(r,Radius,1e-3)))THEN
      IF(isElem(iElem).EQV.ElementIsInside)THEN ! only check elements that were not EXCLUDED in 1.) and invert them
        isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
      END IF
    END IF
  END DO; END DO; END DO; END DO !iElem,k,j,i
END IF

! ----------------------------------------------------------------------------------------------------------------------------------
! 3.) Additionally check special geometries (e.g. Lenses defined by rotationally symmetry geometries)
! ----------------------------------------------------------------------------------------------------------------------------------
IF(PRESENT(GeometryName))THEN
  ! Set the geometrical coordinates (e.g. Axial symmetric with r(x) dependency)
  CALL SetGeometry(GeometryName)

  ! Inquire if DOFs/Elems are within/outside of a region
  SELECT CASE(TRIM(GeometryName))
  CASE('FH_lens')
    ! Loop every element and compare the DOF position
    DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ! x-axis symmetric geometry: get interpolated radius of lens geometry -> r_interpolated(x)
      CALL InterpolateGeometry(Elem_xGP(1,i,j,k,iElem),dim_x=1,dim_y=2,x_OUT=rInterpolated) ! Scale radius

      ! Calculate 2D radius for y-z-plane for comparison with interpolated lens radius
      r = SQRT(&
          Elem_xGP(2,i,j,k,iElem)**2+&
          Elem_xGP(3,i,j,k,iElem)**2  )

      ! Check if r is larger than the interpolated radius of the geometry .AND.
      ! if r is not almost equal to the radius invert the "isElem" logical value
      IF(r.GT.rInterpolated.AND.(.NOT.ALMOSTEQUALRELATIVE(r,rInterpolated,1e-3)))THEN
        IF(isElem(iElem).EQV.ElementIsInside)THEN ! only check elements that were not EXCLUDED in 1.) and invert them
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
        END IF
      END IF
    END DO; END DO; END DO; END DO !iElem,k,j,i
  CASE('FishEyeLens')
    ! Nothing to do, because the geometry is set by using the sphere's radius in step 2.)
  CASE('DielectricResonatorAntenna') ! Radius is checked, but only in x-y (not z)
    DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IF(isElem(iElem).EQV.ElementIsInside)THEN ! only check elements that were not EXCLUDED in 1.) and invert them

        ! Calculate 2D radius for x-y-plane
        r = SQRT(Elem_xGP(1,i,j,k,iElem)**2 + Elem_xGP(2,i,j,k,iElem)**2)

        ! Only perform check for elements in z = Elem_xGP(3) > 0
        IF(Elem_xGP(3,i,j,k,iElem).GT.0.0)THEN
          ! Check if r is larger than the supplied value .AND. if r is not almost equal to the radius
          IF(r.GT.Radius.AND.(.NOT.ALMOSTEQUALRELATIVE(r,Radius,1e-3)))THEN
              isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
          END IF ! r.GT.Radius.AND.(.NOT.ALMOSTEQUALRELATIVE(r,Radius,1e-3))
        END IF ! Elem_xGP(3,i,j,k,iElem).GT.0.0
      END IF ! isElem(iElem).EQV.ElementIsInside
    END DO; END DO; END DO; END DO !iElem,k,j,i
  CASE('Circle') ! Radius is checked
    SELECT CASE(GeometryAxis)
      CASE(1) ! x-axis
        dim_1 = 2
        dim_2 = 3
      CASE(2) ! y-axis
        dim_1 = 1
        dim_2 = 3
      CASE(3) ! z-axis
        dim_1 = 1
        dim_2 = 2
      CASE DEFAULT
        SWRITE(UNIT_stdOut,'(A)') ' '
        CALL abort(__STAMP__,'Error in CALL FindElementInRegion(GeometryName): GeometryAxis is wrong!')
    END SELECT
    DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IF(isElem(iElem).EQV.ElementIsInside)THEN ! only check elements that were not EXCLUDED in 1.)

        ! Calculate 2D radius for x-y-plane
        r = SQRT(Elem_xGP(dim_1,i,j,k,iElem)**2 + Elem_xGP(dim_2,i,j,k,iElem)**2)

        ! Check if r is larger than the supplied value .AND. if r is not almost equal to the radius
        IF(r.GT.Radius.AND.(.NOT.ALMOSTEQUALRELATIVE(r,Radius,1e-3)))THEN
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
        END IF

        ! For dielectric regions, check (optional) 2nd radius and exclude regions within the radius
        ! Check if r is smaller than the radius DielectricRadiusValueB .AND. if r is not almost equal to the radius DielectricRadiusValueB
        IF(DielectricRadiusValueB.GT.0.0)THEN
          IF(r.LT.DielectricRadiusValueB.AND.(.NOT.ALMOSTEQUALRELATIVE(r,DielectricRadiusValueB,1e-3)))THEN
            isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
          END IF ! r.LT.DielectricRadiusValueB
        END IF ! DielectricRadiusValueB.GT.0.0
      END IF ! isElem(iElem).EQV.ElementIsInside
    END DO; END DO; END DO; END DO !iElem,k,j,i
  CASE('HollowCircle') ! Inner (r.LT.DielectricRadiusValueB) and outer radius (Radius) are checked: region between is excluded
    SELECT CASE(GeometryAxis)
      CASE(1) ! x-axis
        dim_1 = 2
        dim_2 = 3
      CASE(2) ! y-axis
        dim_1 = 1
        dim_2 = 3
      CASE(3) ! z-axis
        dim_1 = 1
        dim_2 = 2
      CASE DEFAULT
        SWRITE(UNIT_stdOut,'(A)') ' '
        CALL abort(__STAMP__,'Error in CALL FindElementInRegion(GeometryName): GeometryAxis is wrong!')
    END SELECT

    DO iElem=1,PP_nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IF(isElem(iElem).EQV.ElementIsInside)THEN ! only check elements that were not EXCLUDED in 1.)

        ! Calculate 2D radius for x-y-plane
        r = SQRT(Elem_xGP(dim_1,i,j,k,iElem)**2 + Elem_xGP(dim_2,i,j,k,iElem)**2)

        ! Check if r is smaller than the supplied value .AND. if r is not almost equal to the radius
        IF(r.LT.Radius.AND.(.NOT.ALMOSTEQUALRELATIVE(r,Radius,1e-3)))THEN
          isElem(iElem) = .NOT.ElementIsInside ! EXCLUDE elements outside the region
        END IF

        ! For dielectric regions, check (optional) 2nd radius and INCLUDE regions within the radius
        ! Check if r is smaller than the radius DielectricRadiusValueB .AND. if r is not almost equal to the radius DielectricRadiusValueB
        IF(DielectricRadiusValueB.GT.0.0)THEN
          IF( (r.LT.DielectricRadiusValueB).AND.&
              (.NOT.ALMOSTEQUALRELATIVE(r,DielectricRadiusValueB,1e-3)))THEN
            isElem(iElem) = ElementIsInside ! INCLUDE elements smaller than the radius again
          END IF ! r.LT.DielectricRadiusValueB
        END IF ! DielectricRadiusValueB.GT.0.0
      END IF ! isElem(iElem).EQV.ElementIsInside
    END DO; END DO; END DO; END DO !iElem,k,j,i
  CASE('default')
    ! Nothing to do, because the geometry is set by using the box coordinates
  CASE DEFAULT
    SWRITE(UNIT_stdOut,'(A)') ' '
    SWRITE(UNIT_stdOut,'(A)') ' TRIM(GeometryName)='//TRIM(GeometryName)
    CALL abort(__STAMP__,&
        'Error in CALL FindElementInRegion(GeometryName): GeometryName is not defined! Even dummy geometries must be defined.')
  END SELECT
END IF

END SUBROUTINE  FindElementInRegion


SUBROUTINE FindInterfacesInRegion(isFace,isInterFace,isElem,info_opt)
!===================================================================================================================================
!> Check if a face is in a special region (e.g. Dielectric) and/or connects a special region (e.g. Dielectric) to the physical
!> region. This is used, e.g., for dielectric or PML regions
!> indentifies the following connections and stores them in "isFace_combined"
!> use numbering: 1: Master side is special (e.g. dielectric)
!>                2: Slave  side is special (e.g. dielectric)
!>                3: both sides are special (e.g. dielectric) sides
!>                0: normal face in physical region (no special region involved)
!> ToDo: adjust (reduce) the dimension of "isFace_combined", which currently corresponds to the dimensions of the MPI exchange
!> routines -> need for template fortran files
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars              ,ONLY: NGeo
USE MOD_Mesh_Vars              ,ONLY: nSides,nBCSides
USE MOD_Particle_Surfaces_Vars ,ONLY: BezierControlPoints3D
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                    ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif
!USE MOD_Mesh_Vars             ,ONLY: nBCSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                   :: isElem(1:PP_nElems) ! True/False element: special region
CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: info_opt            ! Optional information regarding the type of faces/interfaces to be found
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: isFace(:)      ! True/False face: special region <-> special region
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: isInterFace(:) ! True/False face: special region <-> physical region (or vice versa)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1,0:PP_N,0:PP_N,1:nSides) :: isFace_Slave,isFace_Master,isFace_combined ! the dimension is only used because of
                                                                                       ! the prolong to face routine and MPI logic
INTEGER                                  :: iSide ! Side iterator
CHARACTER(LEN=255)                       :: info  ! Output info on failure
!===================================================================================================================================
! General workflow:
! 1.  initialize Master, Slave and combined side array (it is a dummy array for which only a scalar value is communicated)
! 2.  prolong elem data 'isElem' (Integer data for true/false to side data (also handles mortar interfaces)
! 3.  MPI: communicate slave sides to master
! 4.  comminucate the values to the slave sides (currently done but not used anywhere)
! 5.  calculate combinded value 'isFace_combined' which determines the type of the interface on the master side, where the
!     information is later used when fluxes are determined
! 6.  loop over all sides and use the calculated value 'isFace_combined' to determine 'isFace' and 'interFace'

! 1.  initialize Master, Slave and combined side array (it is a dummy array for which only a scalar value is communicated)
ALLOCATE(isFace(1:nSides))
ALLOCATE(isInterFace(1:nSides))
isFace=.FALSE.
isInterFace=.FALSE.
! For MPI sides send the info to all other procs
isFace_Slave    = -3.
isFace_Master   = -3.
isFace_combined = -3.


! 2.  prolong elem data 'isElem' (Integer data for true/false to side data (also handles mortar interfaces)
CALL ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.FALSE.) ! Includes Mortar sides
!      !DEBUGGING
!      DO iSide = 1, nSides
!        IF(iSide.LE.nBCSides)THEN
!          WRITE (*,*) "iSide, Master, Slave =", iSide, isFace_Master(1,0,0,iSide), isFace_Slave(1,0,0,iSide)," (BC)"
!        ELSE
!          WRITE (*,*) "iSide, Master, Slave =", iSide, isFace_Master(1,0,0,iSide), isFace_Slave(1,0,0,iSide)
!        END IF ! iSide.LE.nBCSides
!        IF(isFace_Master(1,0,0,iSide).LT.0..OR.isFace_Slave(1,0,0,iSide).LT.0.)THEN
!          EXIT
!        END IF ! isFace_Master(1,0,0,iSide).LT.0..OR.isFace_Slave(1,0,0,iSide).LT.0.
!      END DO ! iSide = 1, nSides
#if USE_MPI
CALL ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides=.TRUE.)  ! Includes Mortar sides

! 3.  MPI: communicate slave sides to master
!         send Slave special region info (real with [0=no special region] or [1=special region] as (N+1)*(N+1) array) to Master proc
CALL StartReceiveMPIData(1,isFace_Slave,1,nSides ,RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIData(   1,isFace_Slave,1,nSides,SendRequest_U2,SendID=2) ! Send YOUR
CALL FinishExchangeMPIData(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE -receive YOUR
#endif /*USE_MPI*/


! 4.  Communicate the values to the slave sides
!         communicate information to slave sides
CALL Flux_Mortar_SideInfo(isFace_Master,isFace_Slave,doMPISides=.FALSE.)

#if USE_MPI
CALL Flux_Mortar_SideInfo(isFace_Master,isFace_Slave,doMPISides=.TRUE.)
! send Master special region info (real with [0=no special region] or [1=special region] as (N+1)*(N+1) array) to Slave procs
CALL StartReceiveMPIData(1,isFace_Master,1,nSides ,RecRequest_U,SendID=1) ! Receive YOUR
CALL StartSendMPIData(   1,isFace_Master,1,nSides,SendRequest_U,SendID=1) ! Send MINE
CALL FinishExchangeMPIData(SendRequest_U ,RecRequest_U ,SendID=1) !Send YOUR -receive MINE
!#else ! is this required for MPI=OFF?
!isFace_Slave=isFace_Master
#endif /*USE_MPI*/


! 5.  Calculate combined value 'isFace_combined' which determines the type of the interface on the master side, where the
!     information is later used when fluxes are determined
!         add isFace_Master to isFace_Slave and send
! Build four-states-array for the 4 different combinations phy/phy(0), spec/phy(1), phy/spec(2) and spec/spec(3) a face can be.
!isFace_combined=2*isFace_Master+isFace_Slave
isFace_combined=2*isFace_Slave+isFace_Master

!      !DEBUGGING
!      DO iSide = 1, nSides
!        IF(iSide.LE.nBCSides)THEN
!          WRITE (*,*) "iSide, isFace_combined =", iSide, isFace_combined(1,0,0,iSide)," (BC)"
!        ELSE
!          WRITE (*,*) "iSide, isFace_combined =", iSide, isFace_combined(1,0,0,iSide)
!        END IF ! iSide.LE.nBCSides
!        IF(isFace_combined(1,0,0,iSide).LT.0.)THEN
!          EXIT
!        END IF ! isFace_Master(1,0,0,iSide).LT.0..OR.isFace_Slave(1,0,0,iSide).LT.0.
!      END DO ! iSide = 1, nSides

! use numbering:  2*isFace_Slave+isFace_Master  = 1: Master side is special (e.g. dielectric)
!                                                 2: Slave  side is special (e.g. dielectric)
!                                                 3: both sides are special (e.g. dielectric) sides
!                                                 0: normal face in physical region (no special region involved)


! 6.  loop over all sides and use the calculated value 'isFace_combined' to determine 'isFace' and 'interFace'
!         set 'isFace' for sides that have at least one special region on either side
!         set 'interFace' for sides that are between two different region, e.g., between a PML and the physical domain

DO iSide=1,nSides
  IF(isFace_combined(1,0,0,iSide).GT.0.)THEN
    isFace(iSide)=.TRUE. ! mixed or pure special region face: when my side is not special but neighbor is special
    IF((NINT(isFace_combined(1,0,0,iSide)).EQ.1).OR.&
       (NINT(isFace_combined(1,0,0,iSide)).EQ.2))THEN
        isInterFace(iSide)=.TRUE. ! set all mixed sides as InterFaces, exclude BCs later on
    END IF
  ! The following sanity check is currently deactivated (numbers cannot be <-6.)
  ELSEIF(isFace_combined(1,0,0,iSide).LT.-100.)THEN ! set back to 0. when fixed
    IF(PRESENT(info_opt))THEN
      info=": "//TRIM(info_opt)
    ELSE
      info=''
    END IF ! PRESENT(info_opt)
    IPWRITE (UNIT_stdOut,*) " ERROR in FindInterfacesInRegion()"//TRIM(info)
    IPWRITE (UNIT_stdOut,*) " iSide =", iSide
    IPWRITE (UNIT_stdOut,*) " X: BezierControlPoints3D(1,0:NGeo,0:NGeo,iSide) =", BezierControlPoints3D(1,0:NGeo,0:NGeo,iSide)
    IPWRITE (UNIT_stdOut,*) " Y: BezierControlPoints3D(2,0:NGeo,0:NGeo,iSide) =", BezierControlPoints3D(2,0:NGeo,0:NGeo,iSide)
    IPWRITE (UNIT_stdOut,*) " Z: BezierControlPoints3D(3,0:NGeo,0:NGeo,iSide) =", BezierControlPoints3D(3,0:NGeo,0:NGeo,iSide)
    IPWRITE (UNIT_stdOut,*) " Corner nodes are:"
    IPWRITE (UNIT_stdOut,*) " X: BezierControlPoints3D(1,0/NGeo,0/NGeo,iSide) =", BezierControlPoints3D(1,0,0,iSide),&
                                                                                  BezierControlPoints3D(1,0,NGeo,iSide),&
                                                                                  BezierControlPoints3D(1,NGeo,0,iSide),&
                                                                                  BezierControlPoints3D(1,NGeo,NGeo,iSide)
    IPWRITE (UNIT_stdOut,*) " X: BezierControlPoints3D(2,0/NGeo,0/NGeo,iSide) =", BezierControlPoints3D(2,0,0,iSide),&
                                                                                  BezierControlPoints3D(2,0,NGeo,iSide),&
                                                                                  BezierControlPoints3D(2,NGeo,0,iSide),&
                                                                                  BezierControlPoints3D(2,NGeo,NGeo,iSide)
    IPWRITE (UNIT_stdOut,*) " X: BezierControlPoints3D(3,0/NGeo,0/NGeo,iSide) =", BezierControlPoints3D(3,0,0,iSide),&
                                                                                  BezierControlPoints3D(3,0,NGeo,iSide),&
                                                                                  BezierControlPoints3D(3,NGeo,0,iSide),&
                                                                                  BezierControlPoints3D(3,NGeo,NGeo,iSide)
    CALL abort(__STAMP__,'isFace_combined(1,0,0,iSide).LT.0. -> ',RealInfoOpt=isFace_combined(1,0,0,iSide))
  END IF
END DO
isInterFace(1:nBCSides)=.FALSE. ! BC sides cannot be interfaces!

END SUBROUTINE FindInterfacesInRegion


SUBROUTINE ProlongToFace_ElementInfo(isElem,isFace_Master,isFace_Slave,doMPISides)
!===================================================================================================================================
! Map the element info (isElem) to the surface array (for later MPI communication)
! 1.) Non-Mortar sides (stored in master/slave arrays)
! 2.) Mortar sides (map large side info to multiple smaller sides (stored in master/slave array))
!    ( Copy isFace information from big mortar sides to the small sides. Compare to U_mortar subroutine. )
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides
USE MOD_Mesh_Vars,          ONLY: nBCSides
USE MOD_Mesh_Vars,          ONLY: MortarType,MortarInfo
USE MOD_Mesh_Vars,          ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,          ONLY: firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
LOGICAL,INTENT(IN)              :: isElem(1:PP_nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)              :: isFace_Master(1,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT)              :: isFace_Slave( 1,0:PP_N,0:PP_N,1:nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,ElemID(2),SideID,flip(2),LocSideID(2)
INTEGER                         :: MortarSideID,locSide
INTEGER                         :: iMortar,nMortars
INTEGER                         :: firstMortarSideID,lastMortarSideID
!===================================================================================================================================
! 1.) Non-Mortar sides (also loops over mortar sides, but they are over-written in 2.) )
IF(.NOT.doMPISides)THEN
  ! loop over all sides (independent if doMPISides=F)
  DO SideID=1,nSides

    ! master side, flip=0
    ElemID(1)    = SideToElem(S2E_ELEM_ID,SideID)
    locSideID(1) = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip(1)=0 ! <<<<<<<<<<< THIS was not set! WHY? SELECT CASE(Flip(i)) produces random integer because memory is not set correctly

    ! neighbor side !ElemID,locSideID and flip =-1 if not existing
    ElemID(2)    = SideToElem(S2E_NB_ELEM_ID,SideID)
    locSideID(2) = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip(2)      = SideToElem(S2E_FLIP,SideID)

    ! loop first master then slave side
    DO i=1,2
      IF(ElemID(i).NE.-1)THEN ! exclude BC and Mortar sides
        SELECT CASE(Flip(i))
          CASE(0)   ! master side
            isFace_Master(:,:,:,SideID)=MERGE(1.,0.,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
          CASE(1:4) ! slave side
            isFace_Slave( :,:,:,SideID)=MERGE(1.,0.,isElem(ElemID(i))) ! if isElem(ElemID(i))=.TRUE. -> 1, else 0
        END SELECT
      END IF
    END DO !i=1,2, masterside & slave side
  END DO !SideID
  isFace_Slave(:,:,:,1:nBCSides)=isFace_Master(:,:,:,1:nBCSides)
END IF


! 2.) Mortar sides (Compare to U_mortar subroutine.)
! Map the solution values from the large side of the mortar interface (which is always stored in master array) to the smaller
! mortar sides (either slave or master)
! get 1st and last SideID depeding on doMPISides=T/F
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

! loop over all mortar sides
DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1) ! get number of mortar sides
  locSide  = MortarType(2,MortarSideID)                 ! get loc number

  ! loop over all local mortar sides
  DO iMortar=1,nMortars

    ! get SideID and flip for mapping the array elements
    SideID   = MortarInfo(MI_SIDEID,iMortar,locSide)
    flip(1)  = MortarInfo(MI_FLIP,iMortar,locSide)

    ! map according to master/slave side
    SELECT CASE(flip(1))
      CASE(0)   ! master side
        isFace_Master(:,:,:,SideID)=isFace_Master(:,:,:,MortarSideID) ! should be taken from master (large mortar side is master)
      CASE(1:4) ! slave side
        isFace_Slave( :,:,:,SideID)=isFace_Master(:,:,:,MortarSideID) ! should be taken from master (large mortar side is master)
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID

END SUBROUTINE ProlongToFace_ElementInfo



SUBROUTINE Flux_Mortar_SideInfo(isFace_Master,isFace_Slave,doMPISides)
!===================================================================================================================================
!> Map the flux values from the small mortar sides (either slave or master) to the larger side (which is always stored in master
!> array)
!===================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)      :: isFace_Master(1,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(INOUT)      :: isFace_Slave( 1,0:PP_N,0:PP_N,1:nSides)
LOGICAL,INTENT(IN)      :: doMPISides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iMortar,nMortars
INTEGER                 :: firstMortarSideID,lastMortarSideID
INTEGER                 :: MortarSideID,SideID,iSide,flip
!===================================================================================================================================
! get 1st and last SideID depending on doMPISides=T/F
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

! loop over all mortar sides
DO MortarSideID=firstMortarSideID,lastMortarSideID
  nMortars = MERGE(4,2,MortarType(1,MortarSideID).EQ.1) ! get number of mortar sides
  iSide    = MortarType(2,MortarSideID)                 ! get loc number

  ! loop over all local mortar sides
  DO iMortar=1,nMortars

    ! get SideID and flip for mapping the array elements
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)

    ! map according to master/slave side
    SELECT CASE(flip)
    CASE(0)   ! master side
      isFace_Master(:,:,:,MortarSideID)=isFace_Master(:,:,:,SideID) ! should be written to Master (large mortar side is master)
    CASE(1:4) ! slave sides (should only occur for MPI)
      isFace_Master(:,:,:,MortarSideID)=isFace_Slave(:,:,:,SideID) ! should be written to Master (large mortar side is master)
    END SELECT
  END DO !iMortar
END DO !MortarSideID

END SUBROUTINE Flux_Mortar_SideInfo


SUBROUTINE CountAndCreateMappings(TypeName,&
                                  isElem,isFace,isInterFace,&
                                  nElems,nFaces, nInterFaces,&
                                  ElemToX,XToElem,&
                                  FaceToX,XToFace,&
                                  FaceToXInter,XInterToFace,&
                                  DisplayInfo)
!===================================================================================================================================
!> 1.) Count the number of Elements, Faces and Interfaces of the PML/Dielectric/... region
!> 2.) Create mappings from general element to PML/Dielectric/... element and vice versa
!>                                  face    to PML/Dielectric/... face or interface and vice vesa
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars        ,ONLY: nSides,nGlobalElems
#if USE_MPI
USE MOD_Mesh_Vars        ,ONLY: ElemToSide
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                :: isElem(:),isFace(:),isInterFace(:)
CHARACTER(LEN=*),INTENT(IN)       :: TypeName
LOGICAL,INTENT(IN),OPTIONAL       :: DisplayInfo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)             :: nFaces,nInterFaces,nElems
INTEGER,ALLOCATABLE,INTENT(INOUT) :: ElemToX(:),XToElem(:),FaceToX(:),XToFace(:),FaceToXInter(:),XInterToFace(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: iElem,iSide,nGlobalSpecialElems,nGlobalFaces,nGlobalInterFaces
INTEGER                           :: iXElem,iXFace,iXInterFace,sumGlobalFaces,sumGlobalInterFaces
#if USE_MPI
INTEGER                           :: SideID,nMasterfaces,nMasterInterFaces
#endif
!===================================================================================================================================
! Get number of Elems
nFaces = 0
nInterFaces = 0
nElems = 0
DO iSide=1,nSides
  IF(isFace(iSide))THEN
    nFaces=nFaces+1
  END IF
END DO ! iSide
DO iSide=1,nSides
  IF(isInterFace(iSide))THEN
    nInterFaces=nInterFaces+1
  END IF
END DO ! iSide
DO iElem=1,PP_nElems
  IF(isElem(iElem))THEN
    nElems=nElems+1
  END IF
END DO ! iElem

!===================================================================================================================================
! display face number infos
!===================================================================================================================================
IF(PRESENT(DisplayInfo))THEN
  IF(DisplayInfo)THEN
#if USE_MPI
    nMasterFaces      = 0
    nMasterInterFaces = 0
    DO iElem=1,nElems ! loop over all local elems
      DO iSide=1,6    ! loop over all local sides
        IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
          SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
          IF(isFace(SideID))THEN
            nMasterfaces=nMasterfaces+1
          END IF
          IF(isInterFace(SideID))THEN
            nMasterInterFaces=nMasterInterFaces+1
          END IF
        END IF
      END DO
    END DO
    sumGlobalFaces      = 0
    sumGlobalInterFaces = 0
    CALL MPI_REDUCE(nElems           ,nGlobalSpecialElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,iError)
    CALL MPI_REDUCE(nMasterfaces     ,nGlobalFaces       ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
    CALL MPI_REDUCE(nMasterInterFaces,nGlobalInterfaces  ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_PICLAS,IERROR)
#else
    nGlobalSpecialElems = nElems
    sumGlobalFaces      = nFaces
    sumGlobalInterFaces = nInterFaces
#endif /*USE_MPI*/
    LBWRITE(UNIT_stdOut,'(A,I10,A,I10,A,F6.2,A)')&
    '  Found [',nGlobalSpecialElems,'] nGlobal'//TRIM(TypeName)//'-Elems      inside of '//TRIM(TypeName)//'-region of ['&
    ,nGlobalElems,'] elems in complete domain [',REAL(nGlobalSpecialElems)/REAL(nGlobalElems)*100.,' %]'
    LBWRITE(UNIT_stdOut,'(A,I10,A)')&
    '  Found [',nGlobalFaces       ,'] nGlobal'//TRIM(TypeName)//'-Faces      inside of '//TRIM(TypeName)//'-region.'
    LBWRITE(UNIT_stdOut,'(A,I10,A)')&
    '  Found [',nGlobalInterfaces  ,'] nGlobal'//TRIM(TypeName)//'-InterFaces inside of '//TRIM(TypeName)//'-region.'
  END IF
END IF


!===================================================================================================================================
! create  mappings: element <-> pml-element
!                      face <-> pml-face
!                      face <-> interface
!===================================================================================================================================
ALLOCATE(ElemToX(PP_nElems))
ALLOCATE(XToElem(nElems))
ALLOCATE(FaceToX(nSides))
ALLOCATE(XToFace(nFaces))
ALLOCATE(FaceToXInter(nSides))
ALLOCATE(XInterToFace(nInterFaces))
ElemToX=0
XToElem=0
FaceToX=0
XToFace=0
FaceToXInter=0
XInterToFace=0
! Create array with mapping
iXElem=0
DO iElem=1,PP_nElems
  IF(isElem(iElem))THEN
    iXElem=iXElem+1
    ElemToX(iElem) = iXElem
    XToElem(iXElem) = iElem
  END IF
END DO
iXFace=0
DO iSide=1,nSides
  IF(isFace(iSide))THEN
    iXFace=iXFace+1
    FaceToX(iSide) = iXFace
    XToFace(iXFace) = iSide
  END IF
END DO
iXInterFace=0
DO iSide=1,nSides
  IF(isInterFace(iSide))THEN
    iXInterFace=iXInterFace+1
    FaceToXInter(iSide) = iXInterFace
    XInterToFace(iXInterFace) = iSide
  END IF
END DO

END SUBROUTINE CountAndCreateMappings


SUBROUTINE DisplayRanges(useMinMax_Name,useMinMax,xyzMinMax_name,xyzMinMax,PhysicalMinMax_Name,xyzPhysicalMinMax)
!===================================================================================================================================
! display the ranges of the used min-max-regions for, e.g., PMLs, dielectric regions, etc.
! usually a, e.g., PML/dielectric region is specified or the inverse region, i.e., the physical region is specified
!===================================================================================================================================
! MODULES
USE MOD_Globals          ,ONLY: UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: MPIRoot
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)       :: useMinMax_Name,xyzMinMax_name,PhysicalMinMax_Name
REAL,INTENT(IN)                   :: xyzMinMax(1:6),xyzPhysicalMinMax(1:6)
LOGICAL,INTENT(IN)                :: useMinMax
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! display ranges of special region depending on useMinMax
!LBWRITE(UNIT_stdOut,'(A,L1,A1)')'  '//TRIM(useMinMax_Name)//'=[',useMinMax,']'
IF(useMinMax)THEN
  LBWRITE(UNIT_stdOut,'(A,L1,A1,A)')'  '//TRIM(useMinMax_Name)//'=[',useMinMax,']',': Ranges '//TRIM(xyzMinMax_name)//'(1:6) are'
  CALL DisplayMinMax(xyzMinMax)
ELSE
  LBWRITE(UNIT_stdOut,'(A,L1,A1,A)')'  '//TRIM(useMinMax_Name)//'=[',useMinMax,']',': Ranges '//TRIM(PhysicalMinMax_Name)//'(1:6) are'
  CALL DisplayMinMax(xyzPhysicalMinMax)
END IF
END SUBROUTINE DisplayRanges


SUBROUTINE DisplayMinMax(MinMax)
!===================================================================================================================================
! Display the ranges of a x-y-z min-max region in the vector MinMax(xmin,xmax,ymin,ymax,zmin,zmax)
!===================================================================================================================================
! MODULES
USE MOD_Globals          ,ONLY: UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: MPIRoot
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: MinMax(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: I
!===================================================================================================================================
LBWRITE(UNIT_stdOut,'(A)') '       [        x-dir         ] [        y-dir         ] [         z-dir        ]'
LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MIN'
DO I=1,3
  LBWRITE(UNIT_stdOut,WRITEFORMAT,ADVANCE='NO')  MinMax(2*I-1)
END DO
LBWRITE(UNIT_stdOut,'(A)') ''
LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  MAX'
DO I=1,3
  LBWRITE(UNIT_stdOut,WRITEFORMAT,ADVANCE='NO')  MinMax(2*I)
END DO
LBWRITE(UNIT_stdOut,'(A)') ''
END SUBROUTINE DisplayMinMax


SUBROUTINE SelectMinMaxRegion(TypeName,useMinMax,region1_name,region1,region2_name,region2)
!===================================================================================================================================
! check whether a MinMax region was defined by the user.
!
!===================================================================================================================================
! MODULES
USE MOD_Globals          ,ONLY: UNIT_stdOut
#if USE_MPI
USE MOD_Globals          ,ONLY: MPIRoot
#endif
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)       :: TypeName,region2_name,region1_name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(INOUT)             :: useMinMax
REAL,INTENT(INOUT)                :: region2(1:6),region1(1:6)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(ALMOSTEQUAL(MAXVAL(region1),MINVAL(region1)))THEN ! if still the initialized values
  region1(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
  IF(ALMOSTEQUAL(MAXVAL(region2),MINVAL(region2)))THEN ! if still the initialized values
    region2(1:6)=(/-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.),-HUGE(1.),HUGE(1.)/)
    useMinMax=.FALSE. ! ! region1 and region2 are undefined -> use HUGE for both
    LBWRITE(UNIT_stdOut,'(A)')'  no '//TRIM(TypeName)//' region supplied, setting '//TRIM(region1_name)//'(1:6): Setting [+-HUGE]'
    LBWRITE(UNIT_stdOut,'(A)')'  no '//TRIM(TypeName)//' region supplied, setting '//TRIM(region2_name)//'(1:6): Setting [+-HUGE]'
  ELSE
    LBWRITE(UNIT_stdOut,'(A)')'  '//TRIM(TypeName)//' region supplied via '//TRIM(region2_name)//'(1:6)'
    useMinMax=.TRUE. ! region1 is undefined but region2 is not
  END IF
ELSE
  LBWRITE(UNIT_stdOut,'(A)')'  '//TRIM(TypeName)//' region supplied via '//TRIM(region1_name)//'(1:6)'
END IF
END SUBROUTINE SelectMinMaxRegion


SUBROUTINE FinalizeInterfaces()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_Interfaces_Vars, ONLY:InterfaceRiemann,InterfacesInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
InterfacesInitIsDone = .FALSE.
SDEALLOCATE(InterfaceRiemann)
END SUBROUTINE FinalizeInterfaces


SUBROUTINE SetGeometry(GeometryName)
!===================================================================================================================================
!> Set the coordinates of the pre-defined geometry
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interfaces_Vars  ,ONLY: GeometryIsSet,Geometry,GeometryMin,GeometryMax,GeometryNPoints
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                :: temp_array(:)            !< temporary array
REAL                            :: array_shift
INTEGER                         :: I
INTEGER                         :: dim_2
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*)  ,INTENT(IN)   :: GeometryName             !< name of the pre-defined geometry
!===================================================================================================================================
IF(GeometryIsSet)RETURN

LBWRITE(UNIT_stdOut,'(A)') 'Selecting geometry: ['//TRIM(GeometryName)//']'
SELECT CASE(TRIM(GeometryName))
CASE('FH_lens')
  array_shift=0.0 !-0.038812
  !array_shift=347.6000
  GeometryNPoints=106 !18
  dim_2=2
  ALLOCATE(Geometry(1:GeometryNPoints,1:dim_2)) ! 385 radial- and axis-coordinates for gyrotron tube radius along the axis
  ALLOCATE(temp_array(1:GeometryNPoints*dim_2))
  temp_array=(/&            ! this array will be re-shaped into [385,2] = [z,r]
    !-0.038812,0.000531 ,&
    !-0.038415,0.010012 ,&
    !-0.037116,0.018273 ,&
    !-0.035497,0.028064 ,&
    !-0.032961,0.038164 ,&
    !-0.028610,0.046741 ,&
    !-0.022700,0.058992 ,&
    !-0.018341,0.068486 ,&
    !-0.012759,0.078290 ,&
    !-0.005938,0.090239 ,&
    ! 0.001479,0.100000 ,&
    ! 0.004412,0.088743 ,&
    ! 0.005838,0.078656 ,&
    ! 0.007845,0.065513 ,&
    ! 0.010171,0.053900 ,&
    ! 0.011571,0.041061 ,&
    ! 0.012047,0.027301 ,&
    ! 0.012527,0.013847 /)
    0.0, 0.0, &
0.014816715, 1.99670444, &
0.059267674, 3.993106501, &
0.133331455, 5.988367882, &
0.237014201, 7.982623035, &
0.370259774, 9.974796787, &
0.5330791, 11.965063522, &
0.72539951, 13.952507031, &
0.947200545, 15.936944013, &
1.198438389, 17.917952527, &
1.479015859, 19.894791526, &
1.78894558, 21.867539088, &
2.128062215, 23.835147592, &
2.496399114, 25.797812611, &
2.893776762, 27.754568018, &
3.320176754, 29.705334883, &
3.775487266, 31.649585499, &
4.259557482, 33.586690716, &
4.772390912, 35.51665657, &
5.313721801, 37.438487845, &
5.883606031, 39.352383004, &
6.481751181, 41.257355615, &
7.108157782, 43.153413717, &
7.76262568, 45.039940565, &
8.444974215, 46.916426019, &
9.155178358, 48.782794565, &
9.892887299, 50.638129252, &
10.658174646, 52.482614178, &
11.450634535, 54.315269161, &
12.270297113, 56.136167729, &
13.116858274, 57.944626529, &
13.9901344, 59.740261244, &
14.890047927, 61.522906284, &
15.816183031, 63.291742958, &
16.768621298, 65.046925451, &
17.746849304, 66.787507079, &
18.750932064, 68.513605862, &
19.780449875, 70.224494821, &
20.852188411, 71.975995606, &
21.923926639, 73.727495887, &
22.995664868, 75.478996169, &
24.067403096, 77.23049645, &
25.139141325, 78.981996731, &
26.210879553, 80.733497013, &
27.282617782, 82.484997294, &
28.35435601, 84.236497575, &
29.426094235, 85.987997851, &
30.497832455, 87.739498119, &
31.569570676, 89.490998388, &
32.641308897, 91.242498657, &
33.713047117, 92.993998925, &
34.784785338, 94.745499194, &
35.856523559, 96.496999463, &
36.928261779, 98.248499731, &
38.0, 100.0, &
38.498243083, 98.080926044, &
38.987216136, 96.15937442, &
39.466858851, 94.235473726, &
39.937113175, 92.30935426, &
40.39796875, 90.380955996, &
40.849409285, 88.450239259, &
41.291378608, 86.517337199, &
41.723823144, 84.582383722, &
42.146737013, 82.645304551, &
42.56010212, 80.706076171, &
42.96386685, 78.764835972, &
43.357982733, 76.82171986, &
43.742448366, 74.87663694, &
44.117243934, 72.929583453, &
44.482322763, 70.980700664, &
44.837642126, 69.030124197, &
45.183204585, 67.077746152, &
45.518989127, 65.123584812, &
45.844954386, 63.167784907, &
46.161064061, 61.210478167, &
46.467324017, 59.251539872, &
46.763712601, 57.291013219, &
47.050194074, 55.329045939, &
47.326738936, 53.365762761, &
47.593355114, 51.401026683, &
47.850021726, 49.43490383, &
48.096708937, 47.467544431, &
48.333394428, 45.499061214, &
48.560087212, 43.529309193, &
48.776768186, 41.558377479, &
48.983413646, 39.586418261, &
49.180008392, 37.613528456, &
49.366561602, 35.639561577, &
49.543057295, 33.66462635, &
49.709478068, 31.688876367, &
49.865815696, 29.712386535, &
50.012078877, 27.735016381, &
50.148255745, 25.756892209, &
50.27433531, 23.778168427, &
50.390315494, 21.798896834, &
50.49620446, 19.818948929, &
50.591995523, 17.838462951, &
50.677684152, 15.857593528, &
50.7532738, 13.876359719, &
50.818772099, 11.894656063, &
50.874178088, 9.912630568, &
50.919493685, 7.930437493, &
50.954726488, 5.948063962, &
50.979884484, 3.965431569, &
50.99497289, 1.982692481, &
51.0, 0.0/)
  !Geometry=RESHAPE(temp_array, (/385, 2/))
  Geometry=TRANSPOSE(RESHAPE(temp_array, (/dim_2,GeometryNPoints/)))
  Geometry=Geometry/1000.
  Geometry(:,1)=Geometry(:,1)-array_shift

  LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Geometry-MIN="
  ALLOCATE(GeometryMin(1:dim_2))
  DO I=1,dim_2
    GeometryMin(I)=MINVAL(Geometry(:,I))
    LBWRITE(UNIT_stdOut,'(F24.12)',ADVANCE='NO') GeometryMin(I)
  END DO
  LBWRITE(UNIT_stdOut,'(A)') ' '

  ALLOCATE(GeometryMax(1:dim_2))
  LBWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') "Geometry-MAX="
  DO I=1,dim_2
    GeometryMax(I)=MAXVAL(Geometry(:,I))
    LBWRITE(UNIT_stdOut,'(F24.12)',ADVANCE='NO') GeometryMax(I)
  END DO
  LBWRITE(UNIT_stdOut,'(A)') ' '

  !!!!               ! use X points for averaged gradient
  !!!!               PMLGradientEntry=0
  !!!!               PMLGradientExit =0
  !!!!               ! average taper radius for entry (lower z-coordinate)
  !!!!               DO I=1,PMLNGradientPoints(1)-1
  !!!!               PMLGradientEntry = PMLGradientEntry + (Geometry(I+1,2)-Geometry(I,2))/(Geometry(I+1,1)-Geometry(I,1))
  !!!!               END DO
  !!!!               ! average taper radius for exit (upper z-coordinate)
  !!!!               DO I=1,PMLNGradientPoints(2)-1
  !!!!               PMLGradientExit  = PMLGradientExit  + (Geometry(GeometryNPoints-I+1,2)-Geometry(GeometryNPoints-I,2))/&
  !!!!                                                     (Geometry(GeometryNPoints-I+1,1)-Geometry(GeometryNPoints-I,1))
  !!!!               END DO
  !!!!               PMLGradientEntry = PMLGradientEntry/(PMLNGradientPoints(1)-1)
  !!!!               PMLGradientExit  = PMLGradientExit /(PMLNGradientPoints(2)-1)

CASE('FishEyeLens')
  ! Nothing to do, because the geometry is set by using the spheres radius in 2.)
CASE('DielectricResonatorAntenna') ! Radius only in x-y (not z)
  ! nothing to set, because rotationally symmetric, as defined by a radius in x-y
CASE('Circle') ! Radius only in x-y (not z)
  ! nothing to set, because rotationally symmetric, as defined by one or two radii in x-y: DielectricRadiusValue and DielectricRadiusValueB (optional)
CASE('HollowCircle') ! Radius only in x-y (not z)
  ! nothing to set, because rotationally symmetric, as defined by one or two radii in x-y: DielectricRadiusValue and DielectricRadiusValueB (optional)
CASE('default')
  ! Nothing to do, because the geometry is set by using the box coordinates
CASE DEFAULT
  SWRITE(UNIT_stdOut,'(A)') ' '
  SWRITE(UNIT_stdOut,'(A)') ' TRIM(GeometryName)='//TRIM(GeometryName)
  CALL abort(__STAMP__,'Error in SetGeometry(GeometryName): GeometryName is not defined! Even dummy geometries must be defined.')
END SELECT

GeometryIsSet=.TRUE.


END SUBROUTINE SetGeometry


SUBROUTINE InterpolateGeometry(x_IN,dim_x,dim_y,x_OUT)
!===================================================================================================================================
!> Set the coordinates of the pre-defined geometry
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Interfaces_Vars, ONLY:Geometry,GeometryMin,GeometryMax,GeometryNPoints
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,ALLOCATABLE                :: temp_array(:)            !< temporary array
INTEGER                         :: location
REAL                            :: x1,x2,y1,y2,m
INTEGER                         :: location_Upper,location_Lower
!INTEGER                         :: GeometryNPoints,dim_2
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,INTENT(IN)                  :: x_IN
INTEGER,INTENT(IN)               :: dim_x,dim_y
! OUTPUT VARIABLES
!CHARACTER(LEN=*)  ,INTENT(IN)   :: GeometryName             !< name of the pre-defined geometry
REAL,INTENT(OUT)                 :: x_OUT ! radius
!===================================================================================================================================

IF(     x_IN .LE. GeometryMin(dim_x) )THEN
  !x_OUT=GeometryMin(dim_y)
  x_OUT=HUGE(1.0)
  RETURN
ELSEIF( x_IN .GE. GeometryMax(dim_x) )THEN
  !x_OUT=GeometryMax(dim_y)
  x_OUT=HUGE(1.0)
  RETURN
ELSE ! interpolate
  location = MINLOC(ABS(Geometry(:,dim_x)-x_IN),1) !Determines the location of the element in the array with min value
  IF(x_IN.LE.Geometry(location,dim_x))THEN ! location is upper of x_IN
    IF(location.LE.GeometryNPoints)THEN
      location_Upper = location
      location_Lower = location-1
    ELSE
      CALL abort(__STAMP__,'InterpolateGeometry Error when location is upper of x_IN: interpolated outside of array dimension?!')
    END IF
  ELSEIF(x_IN.GE.Geometry(location,dim_x))THEN ! location is lower of x_IN
    IF(location.GE.1)THEN
      location_Upper = location+1
      location_Lower = location
    ELSE
      CALL abort(__STAMP__,'InterpolateGeometry Error when location is lower of x_IN: interpolated outside of array dimension?!')
    END IF
  ELSE
    CALL abort(__STAMP__,'InterpolateGeometry Error location: interpolated outside of array dimension?!')
  END IF

  ! sanity check
  IF((location.LT.1).OR.(location.GT.GeometryNPoints))THEN
    CALL abort(__STAMP__,'InterpolateGeometry Error: interpolated outside of array dimension?!')
  END IF
  ! do the actual interpolation
  x1 = Geometry(location_Lower,dim_x)
  x2 = Geometry(location_Upper,dim_x)
  y1 = Geometry(location_Lower,dim_y)
  y2 = Geometry(location_Upper,dim_y)
  ! calc gradient
  m=(y2-y1)/(x2-x1)
END IF


! linear interpolation
x_OUT = y1 + m*(x_IN-x1)




END SUBROUTINE InterpolateGeometry

END MODULE MOD_Interfaces
