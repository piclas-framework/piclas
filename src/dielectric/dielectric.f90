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


MODULE MOD_Dielectric
!===================================================================================================================================
! Dielectric material handling in Maxwell's (maxwell dielectric) or Poisson's (HDG dielectric) equations
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
INTERFACE InitDielectric
  MODULE PROCEDURE InitDielectric
END INTERFACE
INTERFACE FinalizeDielectric
  MODULE PROCEDURE FinalizeDielectric
END INTERFACE

PUBLIC::InitDielectric,FinalizeDielectric
PUBLIC::DefineParametersDielectric
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

CONTAINS

#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!==================================================================================================================================
!> Define parameters for surfaces (particle-sides)
!==================================================================================================================================
SUBROUTINE DefineParametersDielectric()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Dielectric Region")

CALL prms%CreateLogicalOption(  'DoDielectric'                 , 'Use dielectric regions with EpsR and MuR' , '.FALSE.')
CALL prms%CreateLogicalOption(  'DielectricFluxNonConserving'  , 'Use non-conservative fluxes at dielectric interfaces between a dielectric region and vacuum' , '.FALSE.')
CALL prms%CreateRealOption(     'DielectricEpsR'               , 'Relative permittivity' , '1.')
CALL prms%CreateRealOption(     'DielectricMuR'                , 'Relative permeability' , '1.')
CALL prms%CreateLogicalOption(  'DielectricNoParticles'        , 'Do not insert/emit particles into dielectric regions' , '.TRUE.')
CALL prms%CreateStringOption(   'DielectricTestCase'           , 'Specific test cases: "FishEyeLens", "FH_lens", "Circle", "HollowCircle"' , 'default')
CALL prms%CreateRealOption(     'DielectricRmax'               , 'Radius parameter for functions' , '1.')
CALL prms%CreateLogicalOption(  'DielectricCheckRadius'        , 'Use additional parameter "DielectricRadiusValue" for checking if a DOF is within a dielectric region' ,'.FALSE.')
CALL prms%CreateRealOption(     'DielectricRadiusValue'        , 'Additional parameter radius for checking if a DOF is within a dielectric region' , '-1.')
CALL prms%CreateIntOption(      'DielectricAxis'               , 'Additional parameter spatial direction (cylinder) if a DOF is within a dielectric region (Default = z-axis)' , '3')
CALL prms%CreateRealOption(     'DielectricRadiusValueB'       , '2nd radius for cutting out circular areas within a dielectric region' , '-1.')
CALL prms%CreateRealArrayOption('xyzPhysicalMinMaxDielectric'  , '[xmin, xmax, ymin, ymax, zmin, zmax] vector for defining a dielectric region by giving the bounding box coordinates of the PHYSICAL region', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealArrayOption('xyzDielectricMinMax'          , '[xmin, xmax, ymin, ymax, zmin, zmax] vector for defining a dielectric region by giving the bounding box coordinates of the DIELECTRIC region', '0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0')
CALL prms%CreateRealOption(     'Dielectric_E_0'               , 'Electric field strength parameter for functions' , '1.')
CALL prms%CreateIntOption(      'DielectricNbrOfZones'         , 'Number of zones (from hopr) when defining the dielectric elements via zones' , '0')
CALL prms%CreateIntArrayOption( 'DielectricZoneID'             , 'ID for each dielectric zone (only when zones are used for defining the dielectric elements)', no=0)
CALL prms%CreateRealArrayOption('DielectricZoneEpsR'           , 'EpsR for each zone (only when zones are used for defining the dielectric elements)', no=0)
CALL prms%CreateRealArrayOption('DielectricZoneMuR'            , 'MuR for each zone (only when zones are used for defining the dielectric elements)', no=0)

END SUBROUTINE DefineParametersDielectric

SUBROUTINE InitDielectric()
!===================================================================================================================================
!  Initialize perfectly matched layer
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_Dielectric_Vars
USE MOD_HDF5_Output_Fields,ONLY: WriteDielectricGlobalToHDF5
USE MOD_Globals_Vars      ,ONLY: c
USE MOD_Interfaces        ,ONLY: FindInterfacesInRegion,FindElementInRegion,CountAndCreateMappings,DisplayRanges,SelectMinMaxRegion
USE MOD_Mesh_Vars         ,ONLY: nElems,ElemInfo,offsetElem
#if ! (USE_HDG)
USE MOD_Equation_Vars     ,ONLY: c_corr
#endif /*if not USE_HDG*/
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
INTEGER           :: iElem,iZone
!===================================================================================================================================
LBWRITE(UNIT_StdOut,'(132("-"))')
LBWRITE(UNIT_stdOut,'(A)') ' INIT Dielectric...'
!===================================================================================================================================
! Readin
!===================================================================================================================================
DoDielectric                     = GETLOGICAL('DoDielectric','.FALSE.')
IF(.NOT.DoDielectric) THEN
  LBWRITE(UNIT_stdOut,'(A)') ' Dielectric region deactivated. '
  nDielectricElems=0
  RETURN
END IF
DielectricNoParticles            = GETLOGICAL('DielectricNoParticles')
DielectricFluxNonConserving      = GETLOGICAL('DielectricFluxNonConserving')
DielectricEpsR                   = GETREAL('DielectricEpsR')
DielectricMuR                    = GETREAL('DielectricMuR')
DielectricTestCase               = GETSTR('DielectricTestCase')
DielectricRmax                   = GETREAL('DielectricRmax')
IF((DielectricEpsR.LE.0.0).OR.(DielectricMuR.LE.0.0))THEN
  CALL abort(__STAMP__,'Dielectric: MuR or EpsR cannot be negative or zero.')
END IF
DielectricEpsR_inv               = 1./(DielectricEpsR)                   ! 1./EpsR
!DielectricConstant_inv           = 1./(DielectricEpsR*DielectricMuR)    ! 1./(EpsR*MuR)
DielectricConstant_RootInv       = 1./sqrt(DielectricEpsR*DielectricMuR) ! 1./sqrt(EpsR*MuR)
#if !(USE_HDG)
eta_c_dielectric                 = (c_corr-DielectricConstant_RootInv)*c ! ( chi - 1./sqrt(EpsR*MuR) ) * c
#endif /*if USE_HDG*/
c_dielectric                     = c*DielectricConstant_RootInv          ! c/sqrt(EpsR*MuR)
c2_dielectric                    = c*c/(DielectricEpsR*DielectricMuR)    ! c**2/(EpsR*MuR)
DielectricCheckRadius            = GETLOGICAL('DielectricCheckRadius')
DielectricRadiusValue            = GETREAL('DielectricRadiusValue')
DielectricCircleAxis             = GETINT('DielectricAxis')
IF(DielectricRadiusValue.LE.0.0) DielectricCheckRadius=.FALSE.
DielectricRadiusValueB           = GETREAL('DielectricRadiusValueB')
! determine Dielectric elements
useDielectricMinMax = .FALSE. ! default
xyzPhysicalMinMaxDielectric(1:6) = GETREALARRAY('xyzPhysicalMinMaxDielectric',6)
xyzDielectricMinMax(1:6)         = GETREALARRAY('xyzDielectricMinMax',6)

! Use ZONES (from hopr)
DielectricNbrOfZones = GETINT('DielectricNbrOfZones')

! Check whether to use a user-defined region (physical or dielectric) OR set the dielectric values by ZONES
IF(DielectricNbrOfZones.GT.0)THEN
  ALLOCATE(DielectricZoneID(1:DielectricNbrOfZones))
  ALLOCATE(DielectricZoneEpsR(1:DielectricNbrOfZones))
  ALLOCATE(DielectricZoneMuR(1:DielectricNbrOfZones))

  DielectricZoneID   = GETINTARRAY( 'DielectricZoneID'   , DielectricNbrOfZones)
  DielectricZoneEpsR = GETREALARRAY('DielectricZoneEpsR' , DielectricNbrOfZones)
  DielectricZoneMuR  = GETREALARRAY('DielectricZoneMuR'  , DielectricNbrOfZones)

  ! Set true/false flag for each element
  ALLOCATE(isDielectricElem(1:PP_nElems))
  isDielectricElem = .FALSE.

  DO iELem = 1, nElems
    DO iZone = 1, DielectricNbrOfZones
      IF(ElemInfo(ELEM_ZONE,offsetElem+iElem).EQ.DielectricZoneID(iZone))THEN
        isDielectricElem(iElem) = .TRUE.
      END IF ! ElemInfo(ELEM_ZONE,offsetElem+iElem).EQ.DielectricZoneID(iZone)
    END DO ! iZone = 1, DielectricNbrOfZones
  END DO ! iELem = 1, nElems


ELSE
  ! use xyzPhysicalMinMaxDielectric before xyzDielectricMinMax:
  ! 1.) check for xyzPhysicalMinMaxDielectric
  ! 2.) check for xyzDielectricMinMax
  CALL SelectMinMaxRegion('Dielectric'                  , useDielectricMinMax         , &
                          'xyzPhysicalMinMaxDielectric' , xyzPhysicalMinMaxDielectric , &
                          'xyzDielectricMinMax'         , xyzDielectricMinMax)

  ! display ranges of Dielectric region depending on useDielectricMinMax
  !CALL DisplayRanges('useDielectricMinMax',useDielectricMinMax,&
                     !'xyzDielectricMinMax',xyzDielectricMinMax(1:6),&
             !'xyzPhysicalMinMaxDielectric',xyzPhysicalMinMaxDielectric(1:6))

  ! find all elements in the Dielectric region
  IF(useDielectricMinMax)THEN ! find all elements located inside of 'xyzMinMax'
    CALL FindElementInRegion(isDielectricElem, xyzDielectricMinMax,&
                             ElementIsInside = .TRUE.,&
                             DoRadius        = DielectricCheckRadius,&
                             Radius          = DielectricRadiusValue,&
                             DisplayInfo     = .TRUE.,&
                             GeometryName    = DielectricTestCase,&
                             GeometryAxis    = DielectricCircleAxis)
  ELSE ! find all elements located outside of 'xyzPhysicalMinMaxDielectric'
    CALL FindElementInRegion(isDielectricElem, xyzPhysicalMinMaxDielectric,&
                             ElementIsInside = .FALSE.,&
                             DoRadius        = DielectricCheckRadius,&
                             Radius          = DielectricRadiusValue,&
                             DisplayInfo     = .TRUE.,&
                             GeometryName    = DielectricTestCase,&
                             GeometryAxis    = DielectricCircleAxis)
  END IF
END IF ! DielectricNbrOfZones.GT.0



! find all faces in the Dielectric region
CALL FindInterfacesInRegion(isDielectricFace,isDielectricInterFace,isDielectricElem,info_opt='find all faces in the Dielectric region')

! Get number of Dielectric Elems, Faces and Interfaces. Create Mappngs Dielectric <-> physical region
CALL CountAndCreateMappings('Dielectric',&
                            isDielectricElem      , isDielectricFace      , isDielectricInterFace       , &
                            nDielectricElems      , nDielectricFaces      , nDielectricInterFaces       , &
                            ElemToDielectric      , DielectricToElem      , & ! these two are allocated
                            FaceToDielectric      , DielectricToFace      , & ! these two are allocated
                            FaceToDielectricInter , DielectricInterToFace , & ! these two are allocated
                            DisplayInfo=.TRUE.)

! Set the dielectric profile function EpsR,MuR=f(x,y,z) in the Dielectric region: only for Maxwell + HDG
! for HDG the volume field is not used, only for output in .h5 file for checking if the region is used correctly
! because in HDG only a constant profile is implemented
CALL SetDielectricVolumeProfile()

#if !(USE_HDG)
  ! Determine dielectric Values on faces and communicate them: only for Maxwell
  CALL SetDielectricFaceProfile()
#else /*if USE_HDG*/
  ! Set HDG diffusion tensor 'chitens' on faces
  CALL SetDielectricFaceProfile_HDG()
  !IF(ANY(IniExactFunc.EQ.(/200,300/)))THEN ! for dielectric sphere/slab case
    ! set dielectric ratio e_io = eps_inner/eps_outer for dielectric sphere depending on wheter
    ! the dielectric reagion is inside the sphere or outside: currently one reagion is assumed vacuum
    IF(useDielectricMinMax)THEN ! dielectric elements are assumed to be located inside of 'xyzMinMax'
      DielectricRatio=DielectricEpsR
    ELSE ! dielectric elements outside of sphere, hence, the inverse value is taken
      DielectricRatio=DielectricEpsR_inv
    END IF
    ! get the axial electric field strength in x-direction of the dielectric sphere setup
    Dielectric_E_0 = GETREAL('Dielectric_E_0')
  !END IF
#endif /*USE_HDG*/

! create a HDF5 file containing the DielectriczetaGlobal field: only for Maxwell
CALL WriteDielectricGlobalToHDF5()

DielectricInitIsDone=.TRUE.
LBWRITE(UNIT_stdOut,'(A)')' INIT Dielectric DONE!'
LBWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDielectric


SUBROUTINE SetDielectricVolumeProfile()
!===================================================================================================================================
! Determine the local Dielectric damping factor in x,y and z-direction using a constant/linear/polynomial/... function
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,            ONLY: N_VolMesh,ElemInfo,offsetElem
USE MOD_Dielectric_Vars,      ONLY: DielectricVol,DielectricNbrOfZones,DielectricZoneID,DielectricZoneEpsR,DielectricZoneMuR
USE MOD_Dielectric_Vars,      ONLY: nDielectricElems,DielectricToElem
USE MOD_Dielectric_Vars,      ONLY: DielectricRmax,DielectricEpsR,DielectricMuR,DielectricTestCase
USE MOD_DG_Vars,              ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iDielectricElem,Nloc,iElem,iZone
REAL                :: r
#if USE_HDG
REAL                :: DielectricEpsMax,DielectricEpsMin
REAL                :: DielectricMuMax,DielectricMuMin
#endif /*USE_HDG*/
!===================================================================================================================================
! Check if there are dielectric elements
IF(nDielectricElems.LT.1) RETURN

! Allocate field variables
ALLOCATE(DielectricVol(1:nDielectricElems))
DO iDielectricElem = 1, nDielectricElems
  iElem = DielectricToElem(iDielectricElem)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(DielectricVol(iDielectricElem)%DielectricEps(0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(DielectricVol(iDielectricElem)%DielectricMu(0:Nloc,0:Nloc,0:Nloc))
  ALLOCATE(DielectricVol(iDielectricElem)%DielectricConstant_inv(0:Nloc,0:Nloc,0:Nloc))
END DO ! iDielectricElem = 1, nDielectricElems

! Check whether to use  dielectric values by ZONES or not
IF(DielectricNbrOfZones.GT.0)THEN

  ! Loop over all dielectric elements
  DO iDielectricElem = 1, nDielectricElems
    ! Loop over all dielectric zones
    DO iZone = 1, DielectricNbrOfZones
    ! Check to which zone the element belongs
      IF(ElemInfo(ELEM_ZONE,offsetElem+DielectricToElem(iDielectricElem)).EQ.DielectricZoneID(iZone))THEN
        DielectricVol(iDielectricElem)%DielectricEps(:,:,:) = DielectricZoneEpsR(iZone)
        DielectricVol(iDielectricElem)%DielectricMu( :,:,:) = DielectricZoneMuR(iZone)
      END IF ! ElemInfo(ELEM_ZONE,offsetElem+DielectricToElem(iDielectricElem)).EQ.DielectricZoneID(iZone)
    END DO ! iZone = 1, DielectricNbrOfZones
  END DO ! iDielectricElem = 1, nDielectricElems

ELSE

  ! Fish eye lens: half sphere filled with gradually changing dielectric medium
  IF(TRIM(DielectricTestCase).EQ.'FishEyeLens')THEN
    ! use function with radial dependence: EpsR=n0^2 / (1 + (r/r_max)^2)^2
    DO iDielectricElem=1,nDielectricElems
      iElem = DielectricToElem(iDielectricElem)
      Nloc  = N_DG_Mapping(2,iElem+offSetElem)
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        r = SQRT(N_VolMesh(iElem)%Elem_xGP(1,i,j,k)**2+&
                 N_VolMesh(iElem)%Elem_xGP(2,i,j,k)**2+&
                 N_VolMesh(iElem)%Elem_xGP(3,i,j,k)**2  )
        DielectricVol(iDielectricElem)%DielectricEps(i,j,k) = 4./((1+(r/DielectricRmax)**2)**2)
      END DO; END DO; END DO
      DielectricVol(iDielectricElem)%DielectricMu(0:Nloc,0:Nloc,0:Nloc) = DielectricMuR
    END DO !iDielectricElem,k,j,i

  ELSE ! simply set values const.
    DO iDielectricElem=1,nDielectricElems
      iElem = DielectricToElem(iDielectricElem)
      Nloc  = N_DG_Mapping(2,iElem+offSetElem)
      DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
        DielectricVol(iDielectricElem)%DielectricEps(i,j,k) = DielectricEpsR
        DielectricVol(iDielectricElem)%DielectricMu( i,j,k) = DielectricMuR
      END DO; END DO; END DO
    END DO !iDielectricElem
  END IF

END IF ! DielectricNbrOfZones.GT.0

! invert the product of EpsR and MuR
DO iDielectricElem=1,nDielectricElems
  iElem = DielectricToElem(iDielectricElem)
  Nloc  = N_DG_Mapping(2,iElem+offSetElem)
  DielectricVol(iDielectricElem)%DielectricConstant_inv(0:Nloc,0:Nloc,0:Nloc) = 1./& ! 1./(EpsR*MuR)
          (DielectricVol(iDielectricElem)%DielectricEps(0:Nloc,0:Nloc,0:Nloc)*&
           DielectricVol(iDielectricElem)%DielectricMu( 0:Nloc,0:Nloc,0:Nloc))
END DO !iDielectricElem

! check if MPI local values differ for HDG only (variable dielectric values are not implemented)
#if USE_HDG
DielectricEpsMax = -HUGE(1.)
DielectricEpsMin = HUGE(1.)
DO iDielectricElem=1,nDielectricElems
  DielectricEpsMax = MAX(DielectricEpsMax,MAXVAL(DielectricVol(iDielectricElem)%DielectricEps(:,:,:)))
  DielectricEpsMin = MIN(DielectricEpsMin,MINVAL(DielectricVol(iDielectricElem)%DielectricEps(:,:,:)))
END DO !iDielectricElem
IF(.NOT.ALMOSTEQUALRELATIVE(DielectricEpsMax,DielectricEpsMin,1e-8))THEN
  IF(nDielectricElems.GT.0)THEN
    CALL abort(__STAMP__&
        ,'Dielectric material values in HDG solver cannot be spatially variable because this feature is not implemented! Delta Eps_R=',&
    RealInfoOpt=DielectricEpsMax-DielectricEpsMin)
  END IF
END IF
DielectricMuMax = -HUGE(1.)
DielectricMuMin = HUGE(1.)
DO iDielectricElem=1,nDielectricElems
  DielectricMuMax = MAX(DielectricMuMax,MAXVAL(DielectricVol(iDielectricElem)%DielectricMu(:,:,:)))
  DielectricMuMin = MIN(DielectricMuMin,MINVAL(DielectricVol(iDielectricElem)%DielectricMu(:,:,:)))
END DO !iDielectricElem
IF(.NOT.ALMOSTEQUALRELATIVE(DielectricMuMax,DielectricMuMin,1e-8))THEN
  IF(nDielectricElems.GT.0)THEN
    CALL abort(__STAMP__&
        ,'Dielectric material values in HDG solver cannot be spatially variable because this feature is not implemented! Delta Mu_R=',&
    RealInfoOpt=DielectricMuMax-DielectricMuMin)
  END IF
END IF
#endif /*USE_HDG*/
END SUBROUTINE SetDielectricVolumeProfile


#if !(USE_HDG)
SUBROUTINE SetDielectricFaceProfile()
!===================================================================================================================================
!> Set the dielectric factor 1./SQRT(EpsR*MuR) for each face DOF in the array "Dielectric_Master".
!> Only the array "Dielectric_Master" is used in the Riemann solver, as only the master calculates the flux array
!> (maybe slave information is used in the future)
!>
!> Note:
!> for MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!> threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!> is not allowed)
!> re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!> information)
!> This could be overcome by using template subroutines .t90 (see FlexiOS)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: isDielectricElem,ElemToDielectric, DielectricSurf, DielectricVol, DielectricVolDummy
USE MOD_Mesh_Vars       ,ONLY: nSides, SideToElem, nElems, offSetElem
USE MOD_DG_Vars         ,ONLY: DG_Elems_master, DG_Elems_slave, N_DG_Mapping
USE MOD_ProlongToFace   ,ONLY: ProlongToFace_TypeBased
USE MOD_FillMortar      ,ONLY: U_Mortar
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI             ,ONLY: StartReceiveMPIDataType,StartSendMPIDataTypeDielectric,FinishExchangeMPIDataTypeDielectric
#endif
!USE MOD_FillMortar      ,ONLY: U_Mortar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if USE_MPI
INTEGER           :: p,q
#endif /*USE_MPI*/
INTEGER           :: nMaster, nSlave, locSideID, iElem, flip, iSide, Nloc
REAL              :: dummy,MinSlave,MinMaster
!===================================================================================================================================
! General workflow:
! 1.  Initialize dummy arrays for Elem/Face
! 2.  Fill dummy element values for non-Dielectric sides
! 3.  Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar)
! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!     is not allowed)
!     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!     information)
! 5.  Send/Receive MPI data
! 6.  Allocate the actually needed arrays containing the dielectric material information on the sides
! 7.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)

! 1.  Initialize dummy arrays for Elem/Face
ALLOCATE(DielectricSurf(nSides))
DO iSide =1, nSides
  nMaster = DG_Elems_master(iSide)
  nSlave = DG_Elems_slave(iSide)
  ALLOCATE(DielectricSurf(iSide)%Dielectric_Master(        0:nMaster,0:nMaster))
  DielectricSurf(iSide)%Dielectric_Master = 0.
  ALLOCATE(DielectricSurf(iSide)%Dielectric_Slave(         0:nSlave ,0:nSlave))
  DielectricSurf(iSide)%Dielectric_Slave = 0.
  ALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Master(1,0:nMaster,0:nMaster))
  DielectricSurf(iSide)%Dielectric_dummy_Master = 0.
  ALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Slave( 1,0:nSlave ,0:nSlave))
  DielectricSurf(iSide)%Dielectric_dummy_Slave = 0.
#if USE_MPI
  ALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Master2(1,0:nMaster,0:nMaster))
  DielectricSurf(iSide)%Dielectric_dummy_Master2 = 0.
  ALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Slave2( 1,0:nSlave ,0:nSlave))
  DielectricSurf(iSide)%Dielectric_dummy_Slave2 = 0.
#endif /*USE_MPI*/
END DO

ALLOCATE(DielectricVolDummy(1:nElems))
DO iElem = 1, nElems
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  ALLOCATE(DielectricVolDummy(iElem)%U(1,0:Nloc,0:Nloc,0:Nloc))
  DielectricVolDummy(iElem)%U = 0.
  ! 2.  Fill dummy element values for non-Dielectric sides
  IF (isDielectricElem(iElem))THEN
    DielectricVolDummy(iElem)%U(1,0:Nloc,0:Nloc,0:Nloc) = &
        SQRT(DielectricVol(ElemToDielectric(iElem))%DielectricConstant_inv(0:Nloc,0:Nloc,0:Nloc))
  ELSE
    DielectricVolDummy(iElem)%U(1,0:Nloc,0:Nloc,0:Nloc) = -2.0
  END IF
END DO ! iElem = 1, nElems

! 3. Map dummy element values to face arrays (prolong to face needs data of dimension PP_nVar) to
!    DielectricSurf(:)%Dielectric_dummy_Master, DielectricSurf(:)%Dielectric_dummy_Slave
CALL ProlongToFace_TypeBased(doDielectricSides=.TRUE., doMPISides=.FALSE.)
CALL U_Mortar(doDielectricSides=.TRUE., doMPISides=.FALSE.)

#if USE_MPI
CALL ProlongToFace_TypeBased(doDielectricSides=.TRUE., doMPISides=.TRUE.)
CALL U_Mortar(doDielectricSides=.TRUE., doMPISides=.TRUE.)

! 4.  For MPI communication, the data on the faces has to be stored in an array which is completely sent to the corresponding MPI
!     threads (one cannot simply send parts of an array using, e.g., "2:5" for an allocated array of dimension "1:5" because this
!     is not allowed)
!     re-map data from dimension PP_nVar (due to prolong to face routine) to 1 (only one dimension is needed to transfer the
!     information)
DO iSide=1,nSides
  !WRITE (*,*) "MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave)   =", MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave)
  DO p=0,PP_N
    DO q=0,PP_N
      DielectricSurf(iSide)%Dielectric_dummy_Master2(1,p,q) = DielectricSurf(iSide)%Dielectric_dummy_Master(1,p,q)
      DielectricSurf(iSide)%Dielectric_dummy_Slave2 (1,p,q) = DielectricSurf(iSide)%Dielectric_dummy_Slave( 1,p,q)
    END DO
  END DO
  !WRITE (*,*) "MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master2),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave2) =", MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master2),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave2)
END DO
  !write(*,*) ""
  !IF(myrank.eq.0) read*; CALL MPI_BARRIER(MPI_COMM_WORLD,iError)

! 5.  Send Slave Dielectric info (real array with dimension (N+1)*(N+1)) to Master procs
!     Dielectric_dummy_Slave2
CALL StartReceiveMPIDataType(       RecRequest_U2,SendID=2) ! Receive MINE
CALL StartSendMPIDataTypeDielectric(SendRequest_U2,SendID=2) ! Send YOUR

! Send Master Dielectric info (real array with dimension (N+1)*(N+1)) to Slave procs
! Dielectric_dummy_Master2
CALL StartReceiveMPIDataType(       RecRequest_U ,SendID=1) ! Receive YOUR
CALL StartSendMPIDataTypeDielectric(SendRequest_U ,SendID=1) ! Send MINE

CALL FinishExchangeMPIDataTypeDielectric(SendRequest_U2,RecRequest_U2,SendID=2) !Send MINE - receive YOUR
CALL FinishExchangeMPIDataTypeDielectric(SendRequest_U, RecRequest_U ,SendID=1) !Send YOUR - receive MINE
#endif /*USE_MPI*/

DEALLOCATE(DielectricVolDummy)

! 6.  With MPI, use dummy array which was used for sending the MPI data
!     or with single execution, directly use prolonged data on face
DO iSide =1, nSides
  nMaster = DG_Elems_master(iSide)
  nSlave  = DG_Elems_slave(iSide)
  !WRITE (*,*) "MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master2),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave2) =", MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Master2),MINVAL(DielectricSurf(iSide)%Dielectric_dummy_Slave2)
#if USE_MPI
  DielectricSurf(iSide)%Dielectric_Master = DielectricSurf(iSide)%Dielectric_dummy_Master2(1,0:nMaster,0:nMaster)
  DielectricSurf(iSide)%Dielectric_Slave  = DielectricSurf(iSide)%Dielectric_dummy_Slave2( 1,0:nSlave ,0:nSlave )
#else
  DielectricSurf(iSide)%Dielectric_Master = DielectricSurf(iSide)%Dielectric_dummy_Master(1,0:nMaster,0:nMaster)
  DielectricSurf(iSide)%Dielectric_Slave  = DielectricSurf(iSide)%Dielectric_dummy_Slave( 1,0:nSlave ,0:nSlave )
#endif /*USE_MPI*/
  !WRITE (*,*) "MINVAL(DielectricSurf(iSide)%Dielectric_Master),MINVAL(DielectricSurf(iSide)%Dielectric_Slave)               =", MINVAL(DielectricSurf(iSide)%Dielectric_Master),MINVAL(DielectricSurf(iSide)%Dielectric_Slave)
  !write(*,*) ""
END DO ! iSide =1, nSides
!IF(myrank.eq.0) read*; CALL MPI_BARRIER(MPI_COMM_WORLD,iError)

  ! 7. Copy slave side to master side if the dielectric region is on the slave side as the master will calculate the flux for the
  !    master and the slave side and it requires the factor 1./SQRT(EpsR*MuR) for the wave travelling into the dielectric region
DO iSide = 1, nSides
  MinSlave  = MINVAL(DielectricSurf(iSide)%Dielectric_Slave(:,:))
  MinMaster = MINVAL(DielectricSurf(iSide)%Dielectric_Master(:,:))
  IF((MinMaster.LT.0.0).AND.(MinSlave.LT.0.0))THEN
    DielectricSurf(iSide)%Dielectric_Master(:,:) = 1.0
  ELSEIF(MinMaster.LT.0.0)THEN
    DielectricSurf(iSide)%Dielectric_Master(:,:) = DielectricSurf(iSide)%Dielectric_Slave(:,:)
  END IF ! (MinMaster.LT.0.0).AND.(MinSlave.LT.0.0)
END DO ! iSide = 1, nSides

! 8.  Check if the default value remains unchanged (negative material constants are not allowed until now)
DO iSide =1, nSides
  dummy = MINVAL(DielectricSurf(iSide)%Dielectric_Master)
  IF(dummy.LT.0.0)THEN
    IPWRITE(UNIT_StdOut,*) "DielectricSurf(iSide)%Dielectric_Master =", DielectricSurf(iSide)%Dielectric_Master
    CALL abort(__STAMP__,'Dielectric material values for Riemann solver not correctly determined. MINVAL(Dielectric_Master)=',&
        RealInfoOpt=dummy)
  END IF
END DO ! iSide =1, nSides

DO iSide =1, nSides
  DEALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Slave)
  DEALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Master)
#if USE_MPI
  DEALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Slave2)
  DEALLOCATE(DielectricSurf(iSide)%Dielectric_dummy_Master2)
#endif /*USE_MPI*/
END DO ! iSide =1, nSides
!IF(myrank.eq.0) read*; CALL MPI_BARRIER(MPI_COMM_WORLD,iError)

END SUBROUTINE SetDielectricFaceProfile
#endif /* not USE_HDG*/


#if USE_HDG
SUBROUTINE SetDielectricFaceProfile_HDG()
!===================================================================================================================================
! set the dielectric factor EpsR for each face DOF in the array "chitens" (constant. on the diagonal matrix)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars ,ONLY: isDielectricElem,DielectricEpsR
USE MOD_Equation_Vars   ,ONLY: chi
USE MOD_Mesh_Vars       ,ONLY: nInnerSides, offSetElem
USE MOD_Mesh_Vars       ,ONLY: ElemToSide
USE MOD_DG_Vars         ,ONLY: N_DG_Mapping
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER :: i,j,k,iElem,Nloc
INTEGER :: p,q,flip,locSideID,SideID
REAL    :: Invdummy(3,3)
!===================================================================================================================================
DO iElem=1,PP_nElems
  ! cycle the loop if no dielectric element is connected to the side
  IF(.NOT.isDielectricElem(iElem)) CYCLE

  !compute field on Gauss-Lobatto points (continuous!)
  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  DO k=0,Nloc ; DO j=0,Nloc ; DO i=0,Nloc
    !CALL CalcChiTens(Elem_xGP(:,i,j,k,iElem),chitens(:,:,i,j,k,iElem),chitensInv(:,:,i,j,k,iElem),DielectricEpsR)
    CALL CalcChiTens(chi(iElem)%tens(   :,:,i,j,k),&
                     chi(iElem)%tensInv(:,:,i,j,k),DielectricEpsR)
  END DO; END DO; END DO !i,j,k

  !DO locSideID=1,6
  !  flip=ElemToSide(E2S_FLIP,LocSideID,iElem)
  !  SideID=ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
  !  IF(.NOT.((flip.NE.0).AND.(SideID.LE.nInnerSides)))THEN
  !    DO q=0,PP_N; DO p=0,PP_N
  !      !CALL CalcChiTens(Face_xGP(:,p,q),chitens_face(:,:,p,q,SideID),Invdummy(:,:),DielectricEpsR)
  !      CALL CalcChiTens(chitens_face(:,:,p,q,SideID),Invdummy(:,:),DielectricEpsR)
  !    END DO; END DO !p, q
  !  END IF
  !END DO !locSideID
END DO
END SUBROUTINE SetDielectricFaceProfile_HDG


SUBROUTINE CalcChiTens(chitens,chitensInv,DielectricEpsR)
!===================================================================================================================================
! calculate diffusion tensor, diffusion coefficient chi1/chi0 along B vector field plus isotropic diffusion 1.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mathtools ,ONLY: INVERSE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: DielectricEpsR
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: chitens(3,3)
REAL,INTENT(OUT),OPTIONAL       :: chitensInv(3,3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! default
chitens=0.
chitens(1,1)=1.
chitens(2,2)=1.
chitens(3,3)=1.

! set diffusion tensor: currently only constant distribution
chitens(1,1)=DielectricEpsR
chitens(2,2)=DielectricEpsR
chitens(3,3)=DielectricEpsR

! inverse of diffusion 3x3 tensor on each gausspoint
chitensInv(:,:)=INVERSE(chitens(:,:))

END SUBROUTINE calcChiTens
#endif /*USE_HDG*/



SUBROUTINE FinalizeDielectric()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Dielectric_Vars
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
IF(.NOT.DoDielectric) RETURN
SDEALLOCATE(DielectricVol)
SDEALLOCATE(DielectricSurf)
SDEALLOCATE(DielectricToElem)
SDEALLOCATE(ElemToDielectric)
SDEALLOCATE(DielectricToFace)
SDEALLOCATE(FaceToDielectricInter)
SDEALLOCATE(DielectricInterToFace)
SDEALLOCATE(FaceToDielectric)
SDEALLOCATE(isDielectricElem)
SDEALLOCATE(isDielectricFace)
SDEALLOCATE(isDielectricInterFace)
SDEALLOCATE(DielectricZoneID)
SDEALLOCATE(DielectricZoneEpsR)
SDEALLOCATE(DielectricZoneMuR)
END SUBROUTINE FinalizeDielectric
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

END MODULE MOD_Dielectric
