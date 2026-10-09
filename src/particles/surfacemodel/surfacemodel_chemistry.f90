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

MODULE MOD_SurfaceModel_Chemistry
!===================================================================================================================================
!> Module for the initialization of the surface chemistry
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
PUBLIC :: DefineParametersSurfaceChemistry, InitializeVariablesSurfaceChemistry, InitSurfaceModelChemistry
PUBLIC :: SurfaceModelChemistry, SurfaceModelEventProbability, SurfChemCoverage
#if USE_MPI
PUBLIC :: ExchangeSurfChemCoverage
#endif
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for the catalysis
!==================================================================================================================================
SUBROUTINE DefineParametersSurfaceChemistry()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection("Surface Chemistry")
CALL prms%CreateIntOption(      'Surface-NumOfReactions','Number of chemical Surface reactions')
CALL prms%CreateIntOption(      'Surface-Species','Bulk species of the boundary')
CALL prms%CreateStringOption(   'Surface-Reaction[$]-Type',  &
                                'No default, options are:/n'//&
                                ' A (adsorption), D (desorption), ER (Eley-Rideal), LH (Langmuir-Hinshelwood), LHD/n'//&
                                ' P (probability-based)', numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Surface-Reaction[$]-SurfName', 'none' ,numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'OverwriteCatParameters', 'Flag to set catalytic parameters manually', '.FALSE.')
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$] (Reactant1, Reactant2)', '0 , 0' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Products'                                                    &
                              ,'Products of Reaction[$] (Product1, Product2, Product3).\n'                      //&
                               'For adsorption (Type=A) the slots have a fixed meaning:\n'                      //&
                               '  0,0,0 : the impacting species itself is adsorbed\n'                           //&
                               '  X,0,0 : X is adsorbed\n'                                                      //&
                               '  X,Y,0 : X is adsorbed, Y is released into the gas phase (dissociative ads.)\n'//&
                               'For all other types the entries are the reaction products without slot meaning.'  &
                              ,'0 , 0, 0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-ReactHeat', &
                                    'Heat flux to or from the surface due to the reaction [K]', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-HeatScaling', &
                                    'Linear dependence of the heat flux on the coverage', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EnergyAccommodation'                                         &
                              ,'Fraction of the reaction/adsorption enthalpy accommodated by the surface. '       &
                              ,'1.', numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Reaction[$]-Inhibition','Inhibition/Coadsorption behaviour due to other reactions', &
                                '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Inhibitors'  &
                                ,'Inhibitors of Reaction[$] (Species1, Species2, Species3)', '0 , 0, 0' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Reaction[$]-Promotion','Promotion/Coadsorption behaviour due to other reactions', &
                                '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Promotors'  &
                                ,'Promotors of Reaction[$] (Species1, Species2, Species3)', '0 , 0, 0' &
                                , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-StickingCoefficient','Ratio of adsorbed to impinging particles on a\n' //&
                                    'reactive surface, Langmuir or Kisluik model', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-DissOrder'                                                   &
                              ,'Number of surface sites involved in Reaction[$].\n'                             //&
                               'Type A: number of adjacent free sites required for the adsorption.\n'           //&
                               'Type D: desorption order (2 = associative desorption).\n'                       //&
                               'Default: 2 for dissociative adsorption, 1 otherwise.'                             &
                              ,'-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EqConstant',  &
                                    'Equilibrium constant between the adsorption and desorption (K), Langmuir: K=1', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-LateralInteraction', &
                                    'Interaction between neighbouring particles (W), Edes = E0 + W*Coverage', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Ca', &
                                    'Desorption prefactor: nu = 10^(Ca + Coverage*Cb)', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Cb', &
                                    'Desorption prefactor: nu = 10^(Ca + Coverage*Cb)', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Prefactor', &
                                    'Arrhenius prefactor for the reaction/desorption', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Energy', &
                                    'Arrhenius energy for the reaction/desorption [K]', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Reaction[$]-CoverageDependence', '', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-MinimumTotalCoverage', &
                                    'Lower coverage bound for a reaction to take place', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-MaximumTotalCoverage', &
                                    'Upper coverage bound for a reaction to take place', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Species[$]-MinimumCoverage', &
                                    'Lower coverage bound of  species X for a catalytic reaction', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Species[$]-MaximumCoverage', &
                                    'Upper coverage bound of species X for a catalytic reaction', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Reaction[$]-TemperatureDependence', '', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-MinimumTemperature', &
                                    'Lower temperature bound for a reaction to take place', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-MaximumTemperature', &
                                    'Upper temperature bound for a reaction to take place', '10000.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Diffusion', 'Diffusion along the surface', '.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-TotalDiffusion', 'Diffusion along all possible surface', '.FALSE.')
CALL prms%CreateIntOption(      'Surface-Reaction[$]-NumOfBoundaries', 'Num of boundaries for surface reaction.', &
                                    numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Boundaries', 'Array of boundary indices of surface reaction.', &
                                    numberedmulti=.TRUE., no=0)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EventProbability', &
                                    'Event probability for the simple probability-based surface chemistry (Type = P)',numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-ProductAccommodation', &
                                    'Reaction-specific translation thermal accommodation of the product species (Type = P), default is to ' //&
                                    'utilize the surface-specific accommodation coefficient (TransACC)', '-1.',numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersSurfaceChemistry


SUBROUTINE InitializeVariablesSurfaceChemistry()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies, SpeciesDatabase,VarTimeStep
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, nPartBound
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, SurfChemReac, DoChemSurface
USE MOD_io_hdf5
USE MOD_HDF5_input              ,ONLY: ReadAttribute, DatasetExists, AttributeExists
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars        ,ONLY: PerformLoadBalance
#endif /*USE_LOADBALANCE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)     :: dsetname
INTEGER(HID_T)        :: file_id_specdb                       !< File identifier
LOGICAL               :: DataSetFound
LOGICAL               :: Attr_Exists
CHARACTER(LEN=32)     :: hilf, hilfBC, hilfSpec
INTEGER               :: iReac, iReac2, iPartBound, iVal, err
INTEGER               :: ReadInNumOfReact
INTEGER               :: iSpec, SpecID, ReactionPathPerSpecies(nSpecies)
REAL                  :: ReacProbTest(nPartBound)
LOGICAL,ALLOCATABLE   :: ReadFromParameterFile(:)
REAL,PARAMETER        :: DissOrderUnset = -1.
!===================================================================================================================================

IF(SurfChem%NumOfReact.LE.0) RETURN

ReadInNumOfReact = SurfChem%NumOfReact
LBWRITE(*,*) '| Number of considered reaction paths on Surfaces: ', SurfChem%NumOfReact

IF (VarTimeStep%UseSpeciesSpecific) THEN
  CALL abort(__STAMP__,'ERROR: Species-specific time steps are not implemented for the surface chemistry!')
END IF
! -----------------------------------------------------------------------------------------------------------------------------------
! 0.) Allocation and initialization
! -----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(SurfChemReac(ReadInNumOfReact))
ALLOCATE(ReadFromParameterFile(ReadInNumOfReact))
ReadFromParameterFile = .TRUE.

! Surface map
ALLOCATE(SurfChem%BoundIsChemSurf(nPartBound))
SurfChem%BoundIsChemSurf = .FALSE.
ALLOCATE(SurfChem%PSMap(nPartBound))
SurfChem%CatBoundNum = 0
DO iPartBound = 1, nPartBound
  ALLOCATE(SurfChem%PSMap(iPartBound)%PureSurfReac(ReadInNumOfReact))
  SurfChem%PSMap(iPartBound)%PureSurfReac = .FALSE.
END DO

SurfChemReac(:)%CatName = 'empty'

! Probability based surface chemistry model
ALLOCATE(SurfChem%EventProbInfo(nSpecies))
SurfChem%EventProbInfo(:)%NumOfReactionPaths = 0
ReactionPathPerSpecies = 0

! -----------------------------------------------------------------------------------------------------------------------------------
! 1.) Reaction definition: reactants, products, type, boundaries and the boundary parameters
! -----------------------------------------------------------------------------------------------------------------------------------
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac

  SurfChemReac(iReac)%Reactants(:) = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
  SurfChemReac(iReac)%Products(:)  = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',3,'0,0,0')
  SurfChemReac(iReac)%ReactType    = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))
  SurfChemReac(iReac)%NumOfBounds  = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries')
  SurfChemReac(iReac)%TempDep      = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-TemperatureDependence','.FALSE.')
  SurfChemReac(iReac)%CovDep       = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-CoverageDependence','.FALSE.')

  ! Species-specific coverage window
  IF (SurfChemReac(iReac)%CovDep) THEN
    ALLOCATE(SurfChemReac(iReac)%MinCov(nSpecies))
    ALLOCATE(SurfChemReac(iReac)%MaxCov(nSpecies))
    DO iSpec = 1, nSpecies
      WRITE(UNIT=hilfSpec,FMT='(I0)') iSpec
      SurfChemReac(iReac)%MinCov(iSpec) = &
        GETREAL('Surface-Reaction'//TRIM(hilf)//'-Species'//TRIM(hilfSpec)//'-MinimumCoverage','0.')
      SurfChemReac(iReac)%MaxCov(iSpec) = &
        GETREAL('Surface-Reaction'//TRIM(hilf)//'-Species'//TRIM(hilfSpec)//'-MaximumCoverage','1.')
    END DO
  END IF

  IF (SurfChemReac(iReac)%NumOfBounds.EQ.0) THEN
    CALL abort(__STAMP__,'ERROR: At least one boundary must be defined for each surface reaction!',IntInfoOpt=iReac)
  END IF
  SurfChemReac(iReac)%Boundaries = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Boundaries',SurfChemReac(iReac)%NumOfBounds)

  ! Sanity check whether the boundary indices are within the range of nPartBound
  IF (ANY(SurfChemReac(iReac)%Boundaries.LT.1).OR.ANY(SurfChemReac(iReac)%Boundaries.GT.nPartBound)) THEN
    CALL abort(__STAMP__,'ERROR: Boundary index out of range for surface reaction: ',IntInfoOpt=iReac)
  END IF

  ! --- Define the surface model -----------------------------------------------------------------------------
  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

  CASE('P')
    ! Simple probability based surface model
    IF (SurfChemReac(iReac)%Reactants(2).NE.0) THEN
      CALL abort(__STAMP__,'ERROR: Probability based model only supports one reactant!',IntInfoOpt=iReac)
    END IF
    DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
      iPartBound = SurfChemReac(iReac)%Boundaries(iVal)
      IF ((PartBound%SurfaceModel(iPartBound).GT.0).AND.(PartBound%SurfaceModel(iPartBound).NE.2)) THEN
        CALL abort(__STAMP__,'ERROR: The surface model has already been set for boundary index: ',IntInfoOpt=iPartBound)
      END IF
      PartBound%SurfaceModel(iPartBound) = 2
    END DO
    SpecID = SurfChemReac(iReac)%Reactants(1)
    SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths = SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths + 1

  CASE('A','D','LH','LHD','ER')
    SurfChemReac(iReac)%CatName = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-SurfName'))

    ! Assign the surface model. The boundary PARAMETERS are read in section 1b, once per boundary, because
    ! they are properties of the boundary and not of the reaction.
    DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
      iPartBound = SurfChemReac(iReac)%Boundaries(iVal)
      IF (PartBound%SurfaceModel(iPartBound).EQ.20) CYCLE
      IF (PartBound%SurfaceModel(iPartBound).GT.0) THEN
        CALL abort(__STAMP__,'ERROR: The surface model has already been set for boundary index: ',IntInfoOpt=iPartBound)
      END IF
      PartBound%SurfaceModel(iPartBound) = 20
      PartBound%Reactive(iPartBound)     = .TRUE.
    END DO

    ! Flag the boundaries and build the boundary-to-reaction map
    DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
      iPartBound = SurfChemReac(iReac)%Boundaries(iVal)
      SurfChem%BoundIsChemSurf(iPartBound)           = .TRUE.
      SurfChem%PSMap(iPartBound)%PureSurfReac(iReac) = .TRUE.
    END DO
    DoChemSurface = .TRUE.

    ! Coadsorbing species, either in the form of inhibitors or promotors
    IF (TRIM(SurfChemReac(iReac)%ReactType).EQ.'A') THEN
      SurfChemReac(iReac)%Inhibition    = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-Inhibition','.FALSE.')
      SurfChemReac(iReac)%Inhibitors(:) = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Inhibitors',3,'0,0,0')
      SurfChemReac(iReac)%Promotion     = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-Promotion','.FALSE.')
      SurfChemReac(iReac)%Promotors(:)  = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Promotors',3,'0,0,0')
    END IF

  CASE DEFAULT
    SWRITE(*,*) ' Reaction Type does not exist: ', TRIM(SurfChemReac(iReac)%ReactType)
    CALL abort(__STAMP__,' ERROR: Surface Reaction Type does not exist!',IntInfoOpt=iReac)
  END SELECT
END DO

! -----------------------------------------------------------------------------------------------------------------------------------
! 1b.) Boundary parameters of the catalytic surfaces
! -----------------------------------------------------------------------------------------------------------------------------------
DO iPartBound = 1, nPartBound
  IF (.NOT.SurfChem%BoundIsChemSurf(iPartBound)) CYCLE

  ! Number of boundaries carrying a catalytic surface reaction
  SurfChem%CatBoundNum = SurfChem%CatBoundNum + 1

  WRITE(UNIT=hilfBC,FMT='(I0)') iPartBound

  PartBound%LatticeVec(iPartBound)     = GETREAL('Part-Boundary'//TRIM(hilfBC)//'-LatticeVector')
  PartBound%MolPerUnitCell(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilfBC)//'-NbrOfMol-UnitCell')

  DO iSpec = 1, nSpecies
    WRITE(UNIT=hilfSpec,FMT='(I0)') iSpec
    PartBound%CoverageIni(iPartBound,iSpec) = &
      GETREAL('Part-Boundary'//TRIM(hilfBC)//'-Species'//TRIM(hilfSpec)//'-Coverage')
    PartBound%MaxCoverage(iPartBound,iSpec) = &
      GETREAL('Part-Boundary'//TRIM(hilfBC)//'-Species'//TRIM(hilfSpec)//'-MaxCoverage')
    IF (PartBound%CoverageIni(iPartBound,iSpec).GT.PartBound%MaxCoverage(iPartBound,iSpec)) THEN
      CALL abort(__STAMP__,'ERROR: Surface coverage can not be larger than the maximum value',IntInfoOpt=iPartBound)
    END IF
  END DO

  PartBound%TotalCoverage(iPartBound)    = SUM(PartBound%CoverageIni(iPartBound,:))
  PartBound%MaxTotalCoverage(iPartBound) = GETREAL('Part-Boundary'//TRIM(hilfBC)//'-MaxTotalCoverage')

  ! Check whether the maximum of the coverage is already reached
  IF (PartBound%TotalCoverage(iPartBound).GT.PartBound%MaxTotalCoverage(iPartBound)) THEN
    CALL abort(__STAMP__,'ERROR: Maximum surface coverage reached.',IntInfoOpt=iPartBound)
  END IF
END DO

! -----------------------------------------------------------------------------------------------------------------------------------
! 2.) Probability based model: allocate the species specific containers with the number of reaction paths
! -----------------------------------------------------------------------------------------------------------------------------------
DO iSpec = 1, nSpecies
  IF (SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths.GT.0) THEN
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ReactionIndex(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ReactionIndex = 0
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ReactionProb(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ReactionProb = 0.
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ProdTransACC(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ProdTransACC = -1.
  END IF
END DO

! -----------------------------------------------------------------------------------------------------------------------------------
! 3.) Global switches
! -----------------------------------------------------------------------------------------------------------------------------------
! Bulk species involved in the reactions
SurfChem%SurfSpecies = GETINT('Surface-Species','0')
IF ((SurfChem%SurfSpecies.LT.0).OR.(SurfChem%SurfSpecies.GT.nSpecies)) THEN
  CALL abort(__STAMP__,'ERROR: Surface-Species must be 0 or a valid species index!',IntInfoOpt=SurfChem%SurfSpecies)
END IF

! Diffusion
SurfChem%Diffusion    = GETLOGICAL('Surface-Diffusion','.FALSE.')
SurfChem%TotDiffusion = GETLOGICAL('Surface-TotalDiffusion','.FALSE.')

! Species database
SurfChem%OverwriteCatParameters = GETLOGICAL('OverwriteCatParameters','.FALSE.')
IF (SpeciesDatabase.EQ.'none') SurfChem%OverwriteCatParameters = .TRUE.

! -----------------------------------------------------------------------------------------------------------------------------------
! 4.) Rate parameters from the species database
!     Only the setup-independent NUMBERS are taken from the database. Reactants, products, type, boundaries and
!     the reaction name always come from the parameter file, since species indices are setup-specific.
! -----------------------------------------------------------------------------------------------------------------------------------
IF ((.NOT.SurfChem%OverwriteCatParameters).AND.(SurfChem%CatBoundNum.GT.0)) THEN
  CALL H5OPEN_F(err)
  CALL H5FOPEN_F(TRIM(SpeciesDatabase),H5F_ACC_RDONLY_F,file_id_specdb,err)

  DO iReac = 1, ReadInNumOfReact
    IF (TRIM(SurfChemReac(iReac)%ReactType).EQ.'P') CYCLE

    dsetname = TRIM('/Surface-Chemistry/'//TRIM(SurfChemReac(iReac)%CatName))
    CALL DatasetExists(file_id_specdb,TRIM(dsetname),DataSetFound)

    IF (.NOT.DataSetFound) THEN
      SWRITE(*,*) 'WARNING: DataSet not found: ['//TRIM(dsetname)//'] ['//TRIM(SpeciesDatabase)//']'
      SWRITE(*,*) '         Falling back to the parameter file for surface reaction ', iReac
      CYCLE
    END IF

    ReadFromParameterFile(iReac) = .FALSE.

    ! The database may define a different reaction type than the parameter file
    CALL ReadAttribute(file_id_specdb,'Type',1,DatasetName=dsetname,StrScalar=SurfChemReac(iReac)%ReactType)

    CALL GetDBReal('ReactHeat'           ,SurfChemReac(iReac)%EReact           ,    0.)
    CALL GetDBReal('HeatScaling'         ,SurfChemReac(iReac)%EScale           ,    0.)
    CALL GetDBReal('EnergyAccommodation' ,SurfChemReac(iReac)%HeatAccommodation,    1.)
    CALL GetDBReal('MinimumTotalCoverage',SurfChemReac(iReac)%MinCovTotal      ,    0.)
    CALL GetDBReal('MaximumTotalCoverage',SurfChemReac(iReac)%MaxCovTotal      ,    1.)
    CALL GetDBReal('MinimumTemperature'  ,SurfChemReac(iReac)%MinTemp          ,    0.)
    CALL GetDBReal('MaximumTemperature'  ,SurfChemReac(iReac)%MaxTemp          ,10000.)

    SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

    CASE('A')
      CALL GetDBReal('StickingCoefficient',SurfChemReac(iReac)%S_initial ,1.)
      CALL GetDBReal('EqConstant'         ,SurfChemReac(iReac)%EqConstant,1.)
      CALL GetDBReal('DissOrder'          ,SurfChemReac(iReac)%DissOrder ,DissOrderUnset)
      ! The adsorption products are derived from Products(:) in section 6.

    CASE('D')
      CALL GetDBReal('LateralInteraction',SurfChemReac(iReac)%W_interact,0.)
      CALL GetDBReal('Ca'                ,SurfChemReac(iReac)%C_a       ,0.)
      CALL GetDBReal('Cb'                ,SurfChemReac(iReac)%C_b       ,0.)
      ! Prefactor 0. switches PureSurfChemistry to the Ca/Cb correlation
      CALL GetDBReal('Prefactor'         ,SurfChemReac(iReac)%Prefactor ,0.)
      CALL GetDBReal('Energy'            ,SurfChemReac(iReac)%E_initial ,0.)
      CALL GetDBReal('DissOrder'         ,SurfChemReac(iReac)%DissOrder ,DissOrderUnset)

    CASE('LH','LHD','ER')
      CALL GetDBReal('Energy'   ,SurfChemReac(iReac)%ArrheniusEnergy,0.)
      CALL GetDBReal('Prefactor',SurfChemReac(iReac)%Prefactor      ,1.)

    CASE DEFAULT
      SWRITE(*,*) ' Reaction Type does not exist: ', TRIM(SurfChemReac(iReac)%ReactType)
      CALL abort(__STAMP__,'ERROR: Surface Reaction Type from the species database does not exist!',IntInfoOpt=iReac)
    END SELECT
  END DO ! iReac

  CALL H5FCLOSE_F(file_id_specdb,err)
  CALL H5CLOSE_F(err)
END IF

! -----------------------------------------------------------------------------------------------------------------------------------
! 5.) Rate parameters from the parameter file, for every reaction the database did not cover
! -----------------------------------------------------------------------------------------------------------------------------------
DO iReac = 1, ReadInNumOfReact
  IF (.NOT.ReadFromParameterFile(iReac)) CYCLE
  WRITE(UNIT=hilf,FMT='(I0)') iReac

  SurfChemReac(iReac)%EReact            = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ReactHeat','0.')
  SurfChemReac(iReac)%EScale            = GETREAL('Surface-Reaction'//TRIM(hilf)//'-HeatScaling','0.')
  SurfChemReac(iReac)%HeatAccommodation = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EnergyAccommodation','1.')
  SurfChemReac(iReac)%MinCovTotal       = GETREAL('Surface-Reaction'//TRIM(hilf)//'-MinimumTotalCoverage','0.')
  SurfChemReac(iReac)%MaxCovTotal       = GETREAL('Surface-Reaction'//TRIM(hilf)//'-MaximumTotalCoverage','1.')
  SurfChemReac(iReac)%MinTemp           = GETREAL('Surface-Reaction'//TRIM(hilf)//'-MinimumTemperature','0.')
  SurfChemReac(iReac)%MaxTemp           = GETREAL('Surface-Reaction'//TRIM(hilf)//'-MaximumTemperature','10000.')

  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

  CASE('A')
    SurfChemReac(iReac)%S_initial  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-StickingCoefficient','1.')
    SurfChemReac(iReac)%EqConstant = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EqConstant','1.')
    ! negative sentinel, the type-dependent default is applied in section 6
    SurfChemReac(iReac)%DissOrder  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','-1.')
    ! The adsorption products are derived from Products(:) in section 6

  CASE('D')
    SurfChemReac(iReac)%W_interact = GETREAL('Surface-Reaction'//TRIM(hilf)//'-LateralInteraction','0.')
    SurfChemReac(iReac)%C_a        = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Ca','0.')
    SurfChemReac(iReac)%C_b        = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Cb','0.')
    ! Prefactor 0. switches PureSurfChemistry to the Ca/Cb correlation
    SurfChemReac(iReac)%Prefactor  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','0.')
    SurfChemReac(iReac)%E_initial  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
    SurfChemReac(iReac)%DissOrder  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','-1.')

  CASE('LH','LHD','ER')
    SurfChemReac(iReac)%ArrheniusEnergy = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
    SurfChemReac(iReac)%Prefactor       = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

  CASE('P')
    SpecID = SurfChemReac(iReac)%Reactants(1)
    ReactionPathPerSpecies(SpecID) = ReactionPathPerSpecies(SpecID) + 1
    SurfChem%EventProbInfo(SpecID)%ReactionIndex(ReactionPathPerSpecies(SpecID)) = iReac
    SurfChem%EventProbInfo(SpecID)%ReactionProb(ReactionPathPerSpecies(SpecID))  = &
      GETREAL('Surface-Reaction'//TRIM(hilf)//'-EventProbability')
    SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID))  = &
      GETREAL('Surface-Reaction'//TRIM(hilf)//'-ProductAccommodation')

    ! Sanity checks of the reaction-specific accommodation coefficient
    IF (SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).NE.-1.) THEN
      ! Must be a fraction
      IF ((SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).LT.0.).OR. &
          (SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).GT.1.)) THEN
        CALL abort(__STAMP__,'Reaction-specific thermal accommodation must be between 0 and 1 for reaction: ',IntInfoOpt=iReac)
      END IF
      ! A non-zero accommodation requires a wall temperature
      IF (SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).GT.0.) THEN
        DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
          iPartBound = SurfChemReac(iReac)%Boundaries(iVal)
          IF (PartBound%WallTemp(iPartBound).EQ.0.) THEN
            CALL abort(__STAMP__,'Reaction-specific thermal accommodation requires a wall temperature for boundary '//&
                       TRIM(PartBound%SourceBoundName(iPartBound))//' used for reaction: ',IntInfoOpt=iReac)
          END IF
        END DO
      END IF
    END IF
  END SELECT
END DO

! -----------------------------------------------------------------------------------------------------------------------------------
! 6.) Derivation of the dependent quantities and validation of the input
!     Runs after BOTH read-in paths, because the species database may overwrite SurfChemReac%ReactType and only
!     here is the final reaction type known.
! -----------------------------------------------------------------------------------------------------------------------------------
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac

  ! --- Checks that apply to every reaction type -------------------------------------------------------------
  IF (SurfChemReac(iReac)%Reactants(1).EQ.0) THEN
    CALL abort(__STAMP__,'ERROR Surface chemistry: Reactant1 must be defined for reaction: '//TRIM(hilf))
  END IF

  !-- The accommodation coefficient is a fraction of the released enthalpy
  IF ((SurfChemReac(iReac)%HeatAccommodation.LT.0.).OR.(SurfChemReac(iReac)%HeatAccommodation.GT.1.)) THEN
    CALL abort(__STAMP__,'ERROR Surface chemistry: EnergyAccommodation must be between 0 and 1 for reaction: '//TRIM(hilf))
  END IF

  ! Coverage and temperature windows must not be inverted
  IF (SurfChemReac(iReac)%MinCovTotal.GT.SurfChemReac(iReac)%MaxCovTotal) THEN
    CALL abort(__STAMP__,'ERROR Surface chemistry: MinimumTotalCoverage > MaximumTotalCoverage for reaction: '//TRIM(hilf))
  END IF
  IF (SurfChemReac(iReac)%MinTemp.GT.SurfChemReac(iReac)%MaxTemp) THEN
    CALL abort(__STAMP__,'ERROR Surface chemistry: MinimumTemperature > MaximumTemperature for reaction: '//TRIM(hilf))
  END IF
  IF (SurfChemReac(iReac)%CovDep) THEN
    DO iSpec = 1, nSpecies
      IF (SurfChemReac(iReac)%MinCov(iSpec).GT.SurfChemReac(iReac)%MaxCov(iSpec)) THEN
        CALL abort(__STAMP__,'ERROR Surface chemistry: MinimumCoverage > MaximumCoverage for reaction '//TRIM(hilf)//&
                   ' and species: ',IntInfoOpt=iSpec)
      END IF
    END DO
  END IF

  ! --- Type-specific derivation ------------------------------------------------------------------------------
  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

  CASE('A')
    !---- slot convention
    !     Products(1): adsorbed species (0 -> the impacting species itself, resolved at run time)
    !     Products(2): species released back into the gas phase -> dissociative adsorption
    !     Products(3): unused
    IF (SurfChemReac(iReac)%Products(3).NE.0) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: only two product slots are allowed for adsorption '//&
                 '(Product1 = adsorbate, Product2 = gas phase fragment). Reaction: '//TRIM(hilf))
    END IF

    SurfChemReac(iReac)%AdsorbedProduct = SurfChemReac(iReac)%Products(1)
    SurfChemReac(iReac)%GasProduct      = SurfChemReac(iReac)%Products(2)
    SurfChemReac(iReac)%DissociativeAds = (SurfChemReac(iReac)%GasProduct.NE.0)

    !---with Product1 = 0 the split of the molecule between surface and gas phase would be undefined.
    !           The replaced check tested GasProduct twice and therefore never caught this.
    IF (SurfChemReac(iReac)%DissociativeAds.AND.(SurfChemReac(iReac)%AdsorbedProduct.EQ.0)) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: dissociative adsorption requires Product1 (adsorbate) '//&
                 'to be set explicitly. Reaction: '//TRIM(hilf))
    END IF

    ! Adjacent free sites required: 2 for dissociative adsorption, 1 otherwise
    IF (SurfChemReac(iReac)%DissOrder.LT.0.) THEN
      SurfChemReac(iReac)%DissOrder = MERGE(2., 1., SurfChemReac(iReac)%DissociativeAds)
    END IF

    ! Enters the Kisliuk model as S_0*(1+K*(1/Theta-1))**(-1)
    IF ((SurfChemReac(iReac)%S_initial.LT.0.).OR.(SurfChemReac(iReac)%S_initial.GT.1.)) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: StickingCoefficient must be between 0 and 1 for reaction: '//TRIM(hilf))
    END IF

  CASE('D')
    ! Desorption order of the Polanyi-Wigner equation. Only 1 and 2 are supported, because PureSurfChemistry
    ! applies a separate prefactor conversion for the associative case (DissOrder = 2).
    IF (SurfChemReac(iReac)%DissOrder.LT.0.) SurfChemReac(iReac)%DissOrder = 1.
    IF ((SurfChemReac(iReac)%DissOrder.NE.1.).AND.(SurfChemReac(iReac)%DissOrder.NE.2.)) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: DissOrder must be 1 or 2 for desorption. Reaction: '//TRIM(hilf))
    END IF
    IF (ALL(SurfChemReac(iReac)%Products(:).EQ.0)) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: no product defined for desorption. Reaction: '//TRIM(hilf))
    END IF

  CASE('LH','LHD','ER')
    !--- not used by these types, but it must not stay undefined
    SurfChemReac(iReac)%DissOrder = 1.
    IF (ALL(SurfChemReac(iReac)%Products(:).EQ.0)) THEN
      CALL abort(__STAMP__,'ERROR Surface chemistry: no product defined for reaction type '//&
                 TRIM(SurfChemReac(iReac)%ReactType)//'. Reaction: '//TRIM(hilf))
    END IF

  CASE('P')
    ! Products are handled through SurfChem%EventProbInfo
    SurfChemReac(iReac)%DissOrder = 1.

  CASE DEFAULT
    CALL abort(__STAMP__,'ERROR Surface chemistry: unknown reaction type '//&
               TRIM(SurfChemReac(iReac)%ReactType)//' for reaction: '//TRIM(hilf))
  END SELECT

  ! Common lower bound, after all type-specific defaults have been applied
  IF (SurfChemReac(iReac)%DissOrder.LT.1.) THEN
    CALL abort(__STAMP__,'ERROR Surface chemistry: DissOrder must be >= 1 for reaction: '//TRIM(hilf))
  END IF
END DO

DEALLOCATE(ReadFromParameterFile)

! -----------------------------------------------------------------------------------------------------------------------------------
! 7.) Sanity check: the total reaction probability of a single species at a boundary must not exceed 1
! -----------------------------------------------------------------------------------------------------------------------------------
DO iSpec = 1, nSpecies
  IF (SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths.EQ.0) CYCLE
  ReacProbTest = 0.
  ! Loop over the reaction paths
  DO iReac = 1, SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths
    iReac2 = SurfChem%EventProbInfo(iSpec)%ReactionIndex(iReac)
    ! Loop over the boundaries where the reaction occurs
    DO iVal = 1, SurfChemReac(iReac2)%NumOfBounds
      iPartBound = SurfChemReac(iReac2)%Boundaries(iVal)
      ReacProbTest(iPartBound) = ReacProbTest(iPartBound) + SurfChem%EventProbInfo(iSpec)%ReactionProb(iReac)
      IF (ReacProbTest(iPartBound).GT.1.) THEN
        CALL abort(__STAMP__,'ERROR: Total probability above unity for species: ',IntInfoOpt=iSpec)
      END IF
    END DO
  END DO
END DO

CONTAINS

!===================================================================================================================================
!> Read an optional REAL attribute of the current species database dataset, or fall back to a default value.
!> Uses host association for file_id_specdb, dsetname and Attr_Exists.
!===================================================================================================================================
SUBROUTINE GetDBReal(AttribName,Value,DefaultValue)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=*),INTENT(IN) :: AttribName
REAL,INTENT(OUT)            :: Value
REAL,INTENT(IN)             :: DefaultValue
!===================================================================================================================================
CALL AttributeExists(file_id_specdb,TRIM(AttribName),TRIM(dsetname),AttrExists=Attr_Exists)
IF (Attr_Exists) THEN
  CALL ReadAttribute(file_id_specdb,TRIM(AttribName),1,DatasetName=dsetname,RealScalar=Value)
ELSE
  Value = DefaultValue
END IF
END SUBROUTINE GetDBReal


END SUBROUTINE InitializeVariablesSurfaceChemistry


SUBROUTINE InitSurfaceModelChemistry()
!===================================================================================================================================
! Allocation of side-specific arrays for chemistry modelling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode, PartBound, nSurfSample, nComputeNodeSurfTotalSides, SurfSide2GlobalSide
USE MOD_SurfaceModel_Vars       ,ONLY: ChemSampWall, ChemDesorpWall, ChemWallProp
! USE MOD_Particle_Surfaces_Vars
#if USE_MPI
USE MOD_MPI_Shared
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
USE MOD_SurfaceModel_Vars       ,ONLY: ChemSampWall_Shared, ChemSampWall_Shared_Win
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared, ChemWallProp_Shared_Win
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iSide, iSpec, iBC, SideID
!===================================================================================================================================

IF(.NOT.SurfTotalSideOnNode) RETURN

ALLOCATE(ChemSampWall(1:nSpecies+1,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemSampWall = 0.0
ALLOCATE(ChemDesorpWall(1:nSpecies,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemDesorpWall = 0.0

#if USE_MPI
CALL Allocate_Shared((/nSpecies+1,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),ChemSampWall_Shared_Win,ChemSampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,ChemSampWall_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  ChemSampWall_Shared = 0.
END IF
CALL BARRIER_AND_SYNC(ChemSampWall_Shared_Win,MPI_COMM_SHARED)

CALL Allocate_Shared((/nSpecies+1,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),ChemWallProp_Shared_Win,ChemWallProp_Shared)
CALL MPI_WIN_LOCK_ALL(0,ChemWallProp_Shared_Win,IERROR)
ChemWallProp => ChemWallProp_Shared
IF (myComputeNodeRank.EQ.0) THEN
  ChemWallProp = 0.
  DO iSide = 1, nComputeNodeSurfTotalSides
    ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
    SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
    iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
    DO iSpec = 1, nSpecies
      ! Initial surface coverage
      ChemWallProp(iSpec,:,:,iSide) = PartBound%CoverageIni(iBC, iSpec)
    END DO
  END DO
END IF
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(ChemWallProp(1:nSpecies+1,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemWallProp = 0.0
DO iSide = 1, nComputeNodeSurfTotalSides
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  DO iSpec = 1, nSpecies
  ! Initial surface coverage
    ChemWallProp(iSpec,:,:,iSide) = PartBound%CoverageIni(iBC, iSpec)
  END DO
END DO

#endif /*USE_MPI*/

END SUBROUTINE InitSurfaceModelChemistry


!===================================================================================================================================
!> Selection and execution of a catalytic gas-surface interaction
!> 0.) Determine the surface parameters: Coverage and number of surface molecules
!> 1.) Calculate the sticking coefficient by the Kisliuk model (adsorption)
!> 2.) Calculate the reaction probability by the Arrhenius equation (bias-free for multiple channels)
!> 3.) Choose the occurring pathway by comparison with a random number
!> 4.) Perform the chosen process
!>   a.) Adsorption: delete the incoming particle and update the surface values, for the special case of dissociative adsorption,
!>       the dissociated half is inserted in the gas phase
!>   b.) ER: delete the incoming particle, update the surface values and create the gas phase products
!===================================================================================================================================
SUBROUTINE SurfaceModelChemistry(PartID,SideID,GlobalElemID,n_Loc,PartPosImpact)
! MODULES
! ROUTINES / FUNCTIONS
USE MOD_Globals                   ,ONLY: abort,UNITVECTOR,OrthoNormVec,DOTPRODUCT,UNIT_StdOut,myRank
USE MOD_part_operations           ,ONLY: RemoveParticle, CreateParticle
USE MOD_part_tools                ,ONLY: GetParticleWeight
USE MOD_SurfaceModel_Tools        ,ONLY: MaxwellScattering, CalcPostWallCollVelo, SurfaceModelEnergyAccommodation
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
! VARIABLES
USE MOD_Globals_Vars              ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars             ,ONLY: PartSpecies,Species,usevMPF,WriteMacroSurfaceValues,nSpecies,PartState
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, GlobalSide2SurfSide, SurfSideArea
USE MOD_SurfaceModel_Vars         ,ONLY: SurfChem, SurfChemReac, ChemWallProp, ChemSampWall
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared, ElemMidPoint_Shared
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)
INTEGER,INTENT(IN) :: PartID, SideID
INTEGER,INTENT(IN) :: GlobalElemID        !< Global element ID of the particle impacting the surface
REAL,INTENT(IN)    :: PartPosImpact(1:3)  !< Position of impact of the bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: locBCID, SurfSideID, speciesID
INTEGER            :: iReac, i, iSpec, iCand, nCand
INTEGER            :: iValProd, iProd, iAdsProd, iGasProd, iReactant
INTEGER            :: iCoads, iCoadsSpec
INTEGER            :: NewPartID, CNElemID
INTEGER            :: SubP, SubQ, nGasProd
LOGICAL            :: DoReac
CHARACTER(LEN=5)   :: InteractionType
REAL               :: RanNum
REAL               :: Coverage, MaxCoverage, TotalCoverage, TestCoverage, Theta
REAL               :: CoAds_Coverage, CoAds_MaxCov
REAL               :: S_0, StickCoeff, EqConstant, DissOrder
REAL               :: WallTemp, VeloSquare
REAL               :: nu, E_act, Rate, Prob
REAL               :: SurfMol, InvSurfMol, AdsDens
REAL               :: AdsHeat, ReacHeat, BetaCoeff
REAL               :: partWeight
REAL               :: NewPos(1:3), NewVelo(1:3)
REAL               :: tang1(1:3), tang2(1:3), WallVelo(1:3)
INTEGER            :: CandReac(SurfChem%NumOfReact)  !< reaction index of the candidate
INTEGER            :: CandSpec(SurfChem%NumOfReact)  !< resolved adsorbate ('ER') or adsorbed product ('A')
REAL               :: CandProb(SurfChem%NumOfReact)  !< probability of the candidate
REAL               :: CandCov (SurfChem%NumOfReact)  !< coverage the probability was evaluated with
REAL               :: CandBeta(SurfChem%NumOfReact)  !< energy accommodation coefficient of the candidate
REAL               :: SumProb, ProbNone, ProbTot, ProbCum
LOGICAL,SAVE       :: SumProbWarnDone = .FALSE.      !< issue the SUM(p)>1 warning only once per rank
REAL,PARAMETER     :: eps  = 1e-6
REAL,PARAMETER     :: eps2 = 1.0-eps
REAL               :: ETrans, CovIncrement
!===================================================================================================================================
! -----------------------------------------------------------------------------------------------------------------------------------
! 0.) Determine the surface parameters: coverage and number of surface molecules
! -----------------------------------------------------------------------------------------------------------------------------------
locBCID     = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID  = GlobalSide2SurfSide(SURF_SIDEID,SideID)
CNElemID    = GetCNElemID(GlobalElemID)
speciesID   = PartSpecies(PartID)
WallTemp    = PartBound%WallTemp(locBCID)

InteractionType = 'None'
nCand           = 0
iReactant       = 0
iAdsProd        = 0
BetaCoeff       = 1.0

SubP = TrackInfo%p
SubQ = TrackInfo%q

! Particle weight (MacroParticleFactor)
partWeight = GetParticleWeight(PartID)
IF (.NOT.usevMPF) partWeight = partWeight * Species(speciesID)%MacroParticleFactor
ETrans = 0.5*Species(speciesID)%MassIC*DOTPRODUCT(PartState(4:6,PartID))

! Number of surface molecules of the sub-surface element
IF (PartBound%LatticeVec(locBCID).GT.0.) THEN
  ! Number of surface molecules depending on the occupancy of the unit cell
  SurfMol = PartBound%MolPerUnitCell(locBCID) * SurfSideArea(SubP,SubQ,SurfSideID) &
          / (PartBound%LatticeVec(locBCID)*PartBound%LatticeVec(locBCID))
ELSE
  ! Alternative: average number of surface molecules per area for a monolayer
  SurfMol = 10.**19 * SurfSideArea(SubP,SubQ,SurfSideID)
END IF

IF (SurfMol.LE.0.) CALL abort(__STAMP__,'ERROR SurfaceModelChemistry: SurfMol <= 0 for surf side ',IntInfoOpt=SurfSideID)
InvSurfMol = 1.0/SurfMol

TotalCoverage = SUM(ChemWallProp(1:nSpecies,SubP,SubQ,SurfSideID)) &
              + SUM(ChemSampWall(1:nSpecies,SubP,SubQ,SurfSideID))*InvSurfMol

! -----------------------------------------------------------------------------------------------------------------------------------
! 1.) First pass: evaluate the probability of every eligible reaction channel.
! -----------------------------------------------------------------------------------------------------------------------------------
DO iReac = 1, SurfChem%NumOfReact

  IF (.NOT.SurfChem%PSMap(locBCID)%PureSurfReac(iReac)) CYCLE
  ! Skip channels for which the impacting species is not a reactant
  IF (.NOT.ANY(SurfChemReac(iReac)%Reactants(:).EQ.speciesID)) CYCLE

  ! Check the coverage window of the reaction
  IF (SurfChemReac(iReac)%CovDep) THEN
    IF (TotalCoverage.LT.SurfChemReac(iReac)%MinCovTotal) CYCLE
    IF (TotalCoverage.GT.SurfChemReac(iReac)%MaxCovTotal) CYCLE
    DoReac = .TRUE.
    DO iSpec = 1,nSpecies
      TestCoverage = ChemWallProp(iSpec,SubP,SubQ,SurfSideID) + ChemSampWall(iSpec,SubP,SubQ,SurfSideID)*InvSurfMol
      IF (TestCoverage.LT.SurfChemReac(iReac)%MinCov(iSpec)) DoReac = .FALSE.
      IF (TestCoverage.GT.SurfChemReac(iReac)%MaxCov(iSpec)) DoReac = .FALSE.
    END DO
    IF (.NOT.DoReac) CYCLE
  END IF

  ! Check the temperature window of the reaction
  IF (SurfChemReac(iReac)%TempDep) THEN
    IF (WallTemp.LT.SurfChemReac(iReac)%MinTemp) CYCLE
    IF (WallTemp.GT.SurfChemReac(iReac)%MaxTemp) CYCLE
  END IF

  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

  ! ---------------------------------------------------------------------------------------------------------------------------------
  ! 1a.) Adsorption: sticking coefficient from the Kisliuk model
  ! ---------------------------------------------------------------------------------------------------------------------------------
  CASE('A')
    iAdsProd = SurfChemReac(iReac)%AdsorbedProduct
    IF (iAdsProd.EQ.0) iAdsProd = speciesID

    ! Absolute coverage in terms of the number of surface molecules (including the current time step)
    Coverage    = ChemWallProp(iAdsProd,SubP,SubQ,SurfSideID) + ChemSampWall(iAdsProd,SubP,SubQ,SurfSideID)*InvSurfMol
    MaxCoverage = PartBound%MaxCoverage(locBCID,iAdsProd)

    DissOrder  = SurfChemReac(iReac)%DissOrder
    S_0        = SurfChemReac(iReac)%S_initial
    EqConstant = SurfChemReac(iReac)%EqConstant

    ! Theta = fraction of free surface sites AVAILABLE for the adsorption. The number of adjacent free sites
    ! REQUIRED by the process enters below through the exponent Theta**DissOrder.
    IF (MaxCoverage.LE.0.) THEN
      Theta = 0.0
    ELSE
      Theta = 1.0 - Coverage/MaxCoverage
    END IF

    IF (SurfChemReac(iReac)%Promotion) THEN
      DO iCoads = 1,SIZE(SurfChemReac(iReac)%Promotors)
        iCoadsSpec = SurfChemReac(iReac)%Promotors(iCoads)
        IF (iCoadsSpec.EQ.0) CYCLE
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        IF (CoAds_MaxCov.LE.0.) CYCLE
        CoAds_Coverage = ChemWallProp(iCoadsSpec,SubP,SubQ,SurfSideID) &
                       + ChemSampWall(iCoadsSpec,SubP,SubQ,SurfSideID)*InvSurfMol
        Theta = Theta + CoAds_Coverage/CoAds_MaxCov
      END DO
    END IF

    IF (SurfChemReac(iReac)%Inhibition) THEN
      DO iCoads = 1,SIZE(SurfChemReac(iReac)%Inhibitors)
        iCoadsSpec = SurfChemReac(iReac)%Inhibitors(iCoads)
        IF (iCoadsSpec.EQ.0) CYCLE
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        IF (CoAds_MaxCov.LE.0.) CYCLE
        CoAds_Coverage = ChemWallProp(iCoadsSpec,SubP,SubQ,SurfSideID) &
                       + ChemSampWall(iCoadsSpec,SubP,SubQ,SurfSideID)*InvSurfMol
        Theta = Theta - CoAds_Coverage/CoAds_MaxCov
      END DO
    END IF

    Theta = MIN(MAX(Theta,0.0),1.0)

    ! The adsorbate units created per molecule enter the total coverage check as well
    IF (SurfChemReac(iReac)%DissociativeAds) THEN
      CovIncrement = partWeight*InvSurfMol
    ELSE
      CovIncrement = SurfChemReac(iReac)%DissOrder*partWeight*InvSurfMol
    END IF
    IF ((Theta.GT.0.0) .AND. ((TotalCoverage + CovIncrement).LE.PartBound%MaxTotalCoverage(locBCID))) THEN
      Theta = Theta**DissOrder
      ! Kisliuk model (for EqConstant=1 and MaxCoverage=1 this reduces to the Langmuir model, StickCoeff=Theta)
      StickCoeff = S_0 * (1.0 + EqConstant*(1.0/Theta - 1.0))**(-1.0)
      StickCoeff = MIN(MAX(StickCoeff,0.0),1.0)
    ELSE
      StickCoeff = 0.0
    END IF

    ! Register the channel as a candidate
    IF (StickCoeff.GT.0.0) THEN
      nCand           = nCand + 1
      CandReac(nCand) = iReac
      CandProb(nCand) = StickCoeff
      CandSpec(nCand) = iAdsProd
      CandCov (nCand) = Coverage
      CandBeta(nCand) = SurfChemReac(iReac)%HeatAccommodation
    END IF

  ! ---------------------------------------------------------------------------------------------------------------------------------
  ! 1b.) Eley-Rideal: reaction probability from the Arrhenius equation
  ! ---------------------------------------------------------------------------------------------------------------------------------
  CASE('ER')
    ! Translational energy of the impacting particle and the corresponding effective temperature.
    nu    = SurfChemReac(iReac)%Prefactor
    E_act = SurfChemReac(iReac)%ArrheniusEnergy

    iReactant = 0
    DO i = 1,SIZE(SurfChemReac(iReac)%Reactants(:))
      IF (SurfChemReac(iReac)%Reactants(i).EQ.0)         CYCLE
      IF (SurfChemReac(iReac)%Reactants(i).EQ.speciesID) CYCLE
      iReactant = SurfChemReac(iReac)%Reactants(i)
    END DO
    IF (iReactant.EQ.0) iReactant = speciesID

    IF (iReactant.EQ.SurfChem%SurfSpecies) THEN
      ! Involvement of the surface bulk species: always fully available
      Coverage = 1.0
    ELSE
      Coverage = ChemWallProp(iReactant,SubP,SubQ,SurfSideID) + ChemSampWall(iReactant,SubP,SubQ,SurfSideID)*InvSurfMol
      ! A coverage above a monolayer still enters the rate equation as 1
      Coverage = MIN(MAX(Coverage,0.0),1.0)
    END IF

    ! Absolute particle density of the sub-surface element
    AdsDens = Coverage * SurfMol / SurfSideArea(SubP,SubQ,SurfSideID)

    ! Reaction rate per area divided by the impingement rate -> dimensionless per-collision probability.
    ! nu is a bimolecular rate coefficient [m^3/s], SQRT(2*PI*m/(k*T)) the inverse Hertz-Knudsen flux per
    ! unit density, so the gas number density cancels.
    IF (ETrans.GT.E_act*BoltzmannConst) THEN
      Rate = nu * AdsDens * (1.0 - E_act*BoltzmannConst/ETrans)   ! activation energy given in [K]
      Prob = SQRT(2.*PI*Species(speciesID)%MassIC/(BoltzmannConst*WallTemp)) * Rate
      Prob = MIN(MAX(Prob,0.0),1.0)
    ELSE
      ! Below the threshold the reaction cannot occur
      Prob = 0.0
    END IF

    ! There must be enough adsorbate on the sub-surface element to consume one particle weight
    IF (partWeight.GT.(Coverage*SurfMol)) Prob = 0.0

    ! Register the channel as a candidate
    IF (Prob.GT.0.0) THEN
      nCand           = nCand + 1
      CandReac(nCand) = iReac
      CandProb(nCand) = Prob
      CandSpec(nCand) = iReactant
      CandCov (nCand) = Coverage
      CandBeta(nCand) = SurfChemReac(iReac)%HeatAccommodation
    END IF

  CASE DEFAULT
    ! Reaction types not handled by this model ('D','LH','LHD',...) are treated in PureSurfChemistry
  END SELECT

END DO ! iReac

! -----------------------------------------------------------------------------------------------------------------------------------
! 2.) Decide whether a reaction occurs and select the channel proportional to its probability
!     P_tot     = 1 - PRODUCT(1-p_i)  : probability that at least one of the competing processes fires.
!                                       Saturates smoothly towards 1 and reduces to SUM(p_i) for small p_i.
!     branching = p_i / SUM(p_i)
! -----------------------------------------------------------------------------------------------------------------------------------
IF (nCand.GT.0) THEN
  SumProb  = 0.0
  ProbNone = 1.0
  DO iCand = 1,nCand
    SumProb  = SumProb  + CandProb(iCand)
    ProbNone = ProbNone * (1.0 - CandProb(iCand))
  END DO
  ProbTot = 1.0 - ProbNone

  IF ((SumProb.GT.1.0).AND.(.NOT.SumProbWarnDone)) THEN
    SumProbWarnDone = .TRUE.
    IPWRITE(UNIT_StdOut,*) 'WARNING SurfaceModelChemistry: SUM of the surface reaction probabilities > 1: ',SumProb
    IPWRITE(UNIT_StdOut,*) '        Boundary index: ',locBCID,' impacting species: ',speciesID
    IPWRITE(UNIT_StdOut,*) '        Consider reducing the time step. This warning is shown only once per rank.'
  END IF

  CALL RANDOM_NUMBER(RanNum)
  IF (RanNum.LT.ProbTot) THEN
    ! Draw the channel with a second random number (roulette wheel over p_i/SUM(p_i))
    CALL RANDOM_NUMBER(RanNum)
    RanNum  = RanNum * SumProb
    ProbCum = 0.0
    iCand   = nCand
    DO i = 1,nCand
      ProbCum = ProbCum + CandProb(i)
      IF (RanNum.LT.ProbCum) THEN
        iCand = i
        EXIT
      END IF
    END DO

    ! Restore the state of the selected channel
    iReac           = CandReac(iCand)
    Coverage        = CandCov (iCand)
    BetaCoeff       = CandBeta(iCand)
    InteractionType = TRIM(SurfChemReac(iReac)%ReactType)
    IF (TRIM(InteractionType).EQ.'A') THEN
      iAdsProd = CandSpec(iCand)
      ! Heat of adsorption as a function of the coverage [J]
      AdsHeat  = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst
    ELSE
      iReactant = CandSpec(iCand)
      ! Reaction heat as a function of the coverage [J]
      ReacHeat  = (SurfChemReac(iReac)%EReact - Coverage*SurfChemReac(iReac)%EScale) * BoltzmannConst
    END IF
  END IF
END IF

! -----------------------------------------------------------------------------------------------------------------------------------
! 3.) Perform the chosen process
! -----------------------------------------------------------------------------------------------------------------------------------
SELECT CASE(TRIM(InteractionType))

! -----------------------------------------------------------------------------------------------------------------------------------
! 3a.) Adsorption: bind the adsorbate, release the gas fragment for dissociative adsorption and delete the
!      incoming particle. Both cases share everything except the optional gas phase product, so they are
!      handled in a single branch (v3, slot convention).
! -----------------------------------------------------------------------------------------------------------------------------------
CASE('A')
  ! Heat flux onto the surface created by the adsorption
  ChemSampWall(nSpecies+1,SubP,SubQ,SurfSideID) = ChemSampWall(nSpecies+1,SubP,SubQ,SurfSideID) + AdsHeat*partWeight

  IF (SurfChemReac(iReac)%DissociativeAds) THEN
    ! One fragment binds to the surface, the other one is released into the gas phase below -> one unit
    ChemSampWall(iAdsProd,SubP,SubQ,SurfSideID) = ChemSampWall(iAdsProd,SubP,SubQ,SurfSideID) + partWeight
  ELSE
    ! The whole molecule is bound and dissociates into DissOrder adsorbate units on the surface,
    ! e.g. O2(gas) -> 2 O(ads) for DissOrder = 2
    ChemSampWall(iAdsProd,SubP,SubQ,SurfSideID) = ChemSampWall(iAdsProd,SubP,SubQ,SurfSideID) &
                                                + SurfChemReac(iReac)%DissOrder*partWeight
  END IF

  ! Dissociative adsorption: the second product slot is released back into the gas phase.
  iGasProd = SurfChemReac(iReac)%GasProduct
  IF (iGasProd.NE.0) THEN
    WallVelo   = PartBound%WallVelo(1:3,locBCID)
    VeloSquare = 2.0*ETrans/Species(iGasProd)%MassIC
    CALL OrthoNormVec(n_loc,tang1,tang2)

    NewVelo(1:3) = CalcPostWallCollVelo(iGasProd,VeloSquare,WallTemp,BetaCoeff)
    NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
    NewPos(1:3)  = eps*ElemMidPoint_Shared(1:3,CNElemID) + eps2*PartPosImpact(1:3)

    CALL CreateParticle(iGasProd,NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0., &
                        NewPartID=NewPartID,NewMPF=partWeight)

    ! Energy transferred from the surface into the internal energies of the new particle
    CALL SurfaceModelEnergyAccommodation(NewPartID,locBCID,WallTemp)

    ! Sampling of the newly created particle
    IF (DSMC%CalcSurfaceVal.AND.(SamplingActive.OR.WriteMacroSurfaceValues)) &
      CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc, PartPosImpact_opt=NewPos)
  END IF

  ! Remove the impinging particle from the gas phase
  CALL RemoveParticle(PartID)

! -----------------------------------------------------------------------------------------------------------------------------------
! 3b.) Eley-Rideal: delete the incoming particle, consume the adsorbate and create the gas phase products
! -----------------------------------------------------------------------------------------------------------------------------------
CASE('ER')
  ! Heat flux onto the surface created by the reaction (accommodated fraction)
  ChemSampWall(nSpecies+1,SubP,SubQ,SurfSideID) = ChemSampWall(nSpecies+1,SubP,SubQ,SurfSideID) &
                                                + ReacHeat*partWeight*BetaCoeff

  WallVelo = PartBound%WallVelo(1:3,locBCID)
  CALL OrthoNormVec(n_loc,tang1,tang2)
  nGasProd = COUNT(SurfChemReac(iReac)%Products(:).GT.0)
  IF (nGasProd.LE.0) CALL abort(__STAMP__,'ERROR SurfaceModelChemistry: ER reaction without product, iReac = ',IntInfoOpt=iReac)

  DO iValProd = 1,SIZE(SurfChemReac(iReac)%Products(:))
    IF (SurfChemReac(iReac)%Products(iValProd).EQ.0) CYCLE
    iProd = SurfChemReac(iReac)%Products(iValProd)
    VeloSquare = 2.0*(ETrans + ReacHeat)/(Species(iProd)%MassIC*REAL(nGasProd))

    NewVelo(1:3) = CalcPostWallCollVelo(iProd,VeloSquare,WallTemp,BetaCoeff)
    NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
    NewPos(1:3)  = eps*ElemMidPoint_Shared(1:3,CNElemID) + eps2*PartPosImpact(1:3)

    CALL CreateParticle(iProd,NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0., &
                        NewPartID=NewPartID,NewMPF=partWeight)

    ! Energy transferred from the surface into the internal energies of the new particle
    CALL SurfaceModelEnergyAccommodation(NewPartID,locBCID,WallTemp)

    ! Sampling of the newly created particles
    IF (DSMC%CalcSurfaceVal.AND.(SamplingActive.OR.WriteMacroSurfaceValues)) &
      CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc, PartPosImpact_opt=NewPos)
  END DO

  IF (iReactant.NE.SurfChem%SurfSpecies) THEN
    ChemSampWall(iReactant,SubP,SubQ,SurfSideID) = ChemSampWall(iReactant,SubP,SubQ,SurfSideID) - partWeight
  END IF

  ! Remove the impinging particle from the gas phase
  CALL RemoveParticle(PartID)

! -----------------------------------------------------------------------------------------------------------------------------------
! 3c.) No reaction: regular wall interaction
! -----------------------------------------------------------------------------------------------------------------------------------
CASE DEFAULT
  CALL MaxwellScattering(PartID,SideID,n_Loc)

END SELECT

END SUBROUTINE SurfaceModelChemistry


!===================================================================================================================================
!> Perform a simple surface reaction based on a fixed probability
!> 1.) Check whether species has any reactions to perform at the boundary and select reaction path
!> 2.) Perform the selected reaction path
!===================================================================================================================================
SUBROUTINE SurfaceModelEventProbability(PartID,SideID,GlobalElemID,n_loc,PartPosImpact)
! MODULES
! ROUTINES / FUNCTIONS
USE MOD_Globals
USE MOD_SurfaceModel_Tools        ,ONLY: MaxwellScattering, CalcPostWallCollVelo, CalcRotWallVelo
USE MOD_SurfaceModel_Tools        ,ONLY: SurfaceModelEnergyAccommodation
USE MOD_part_operations           ,ONLY: CreateParticle, RemoveParticle
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_Mesh_Tools                ,ONLY: GetCNElemID
! VARIABLES
USE MOD_Globals_Vars              ,ONLY: BoltzmannConst
USE MOD_Particle_Vars             ,ONLY: PartSpecies, PartState, usevMPF, PartMPF, WriteMacroSurfaceValues
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, GlobalSide2SurfSide
USE MOD_SurfaceModel_Vars         ,ONLY: SurfChem, SurfChemReac
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared, ElemMidPoint_Shared
USE MOD_DSMC_Vars                 ,ONLY: DSMC, SamplingActive, BGGas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: PartID, SideID
INTEGER,INTENT(IN) :: GlobalElemID          !< Global element ID of the particle impacting the surface
REAL,INTENT(IN)    :: n_loc(1:3)
REAL,INTENT(IN)    :: PartPosImpact(1:3)    !< Charge and position of impact of bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: locBCID, SurfSideID, CNElemID
INTEGER            :: SpecID, ProdSpecID, NewPartID
INTEGER            :: iPath, iProd, ReacTodo, PathTodo
INTEGER            :: NumProd, NumReac
REAL               :: RanNum, WallTemp, TransACC, VeloSquare, TotalProb, OldMPF
REAL               :: tang1(1:3), tang2(1:3), WallVelo(3), NewVelo(3), NewPos(1:3)
REAL,PARAMETER     :: eps=1e-6, eps2=1.0-eps
!===================================================================================================================================
locBCID     = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID  = GlobalSide2SurfSide(SURF_SIDEID,SideID)
CNElemID    = GetCNElemID(GlobalElemID)
WallTemp    = PartBound%WallTemp(locBCID)
WallVelo    = PartBound%WallVelo(1:3,locBCID)
SpecID      = PartSpecies(PartID)
NumReac     = 0
ReacTodo    = 0
PathTodo    = 0
TotalProb   = 0.
IF(PartBound%RotVelo(locBCID)) THEN
  WallVelo(1:3) = CalcRotWallVelo(locBCID,PartPosImpact)
END IF
! ----------------------------------------------------------------------------------------------------------------------------------
! 1.) Check whether species has any reactions to perform at the boundary
IF(SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths.EQ.0) THEN
  PathTodo = 0
ELSE
! 1a.) Determine which reaction path to follow
  CALL RANDOM_NUMBER(RanNum)
  DO iPath = 1, SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths
    ! Check if the reaction is allowed at the boundary
    IF(.NOT.ANY(SurfChemReac(SurfChem%EventProbInfo(SpecID)%ReactionIndex(iPath))%Boundaries(:).EQ.locBCID)) CYCLE
    ! Sum up the probabilities of the reaction paths (sanity check during initialization to ensure that the sum is below 1.0)
    TotalProb = TotalProb + SurfChem%EventProbInfo(SpecID)%ReactionProb(iPath)
    ! Decide which reaction path to follow
    IF(TotalProb.GT.RanNum) THEN
      PathTodo = iPath
      EXIT
    END IF
  END DO
END IF

! 2.) Perform the selected reaction path
IF(PathTodo.GT.0) THEN
  ReacTodo = SurfChem%EventProbInfo(SpecID)%ReactionIndex(PathTodo)
  NumProd = COUNT(SurfChemReac(ReacTodo)%Products(:).GT.0)
  ! Create products if any have been defined
  IF(NumProd.GT.0) THEN
    CALL OrthoNormVec(n_loc,tang1,tang2)
    VeloSquare = DOTPRODUCT(PartState(4:6,PartID)) / NumProd
    IF(SurfChem%EventProbInfo(SpecID)%ProdTransACC(PathTodo).EQ.-1) THEN
      TransACC = PartBound%TransACC(locBCID)
    ELSE
      TransACC = SurfChem%EventProbInfo(SpecID)%ProdTransACC(PathTodo)
    END IF
    DO iProd = 1, NumProd
      ProdSpecID = SurfChemReac(ReacTodo)%Products(iProd)
      ! Do not emit background gas species (but consider them in the energy distribution in VeloSquare)
      IF(BGGas%BackgroundSpecies(ProdSpecID)) CYCLE
      ! Calculate the velocity based on the accommodation coefficient
      NewVelo(1:3) = CalcPostWallCollVelo(ProdSpecID,VeloSquare,WallTemp,TransACC)
      ! Perform vector transformation from the local to the global coordinate system and add wall velocity
      NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
      ! Create new position by using POI and moving the particle by eps in the direction of the element center
      NewPos(1:3) = eps*ElemMidPoint_Shared(1:3,CNElemID) + eps2*PartPosImpact(1:3)
      IF(usevMPF)THEN
        ! Get MPF of old particle
        OldMPF = PartMPF(PartID)
        ! New particle acquires the MPF of the impacting particle (not necessarily the MPF of the newly created particle species)
        CALL CreateParticle(ProdSpecID,NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=OldMPF)
      ELSE
        ! New particle acquires the MPF of the new particle species
        CALL CreateParticle(ProdSpecID,NewPos(1:3),GlobalElemID,GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
      END IF ! usevMPF
      ! Adding the energy that is transferred from the surface onto the internal energies of the particle
      CALL SurfaceModelEnergyAccommodation(NewPartID,locBCID,WallTemp)
      ! Sampling of newly created particles
      IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
        CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc,PartPosImpact_opt=NewPos)
    END DO
  END IF
  ! Remove original reactant
  CALL RemoveParticle(PartID,BCID=locBCID)
ELSE
  CALL MaxwellScattering(PartID,SideID,n_loc)
END IF

END SUBROUTINE SurfaceModelEventProbability


SUBROUTINE SurfChemCoverage()
!===================================================================================================================================
!> calculation of the surface coverage
!> 1) calculation of the coverage (needed without MPI as well)
!> 2) compute-node leaders ensure synchronization of shared arrays on their node
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound
USE MOD_Particle_Vars           ,ONLY: nSpecies,WriteMacroSurfaceValues
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp
USE MOD_Particle_Boundary_vars  ,ONLY: SurfSideArea, SurfSide2GlobalSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode, nComputeNodeSurfTotalSides
USE MOD_DSMC_Vars               ,ONLY: DSMC,SamplingActive
#if USE_MPI
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED, nComputeNodeProcessors
USE MOD_MPI_Shared_Vars         ,ONLY: myComputeNodeRank
USE MOD_SurfaceModel_Vars       ,ONLY: ChemSampWall_Shared, ChemSampWall_Shared_Win
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared_Win
#else
USE MOD_SurfaceModel_Vars       ,ONLY: ChemSampWall
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSide, firstSide, lastSide, GlobalSideID, locBCID, iSpec
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

#if USE_MPI
ASSOCIATE(ChemSampWall => ChemSampWall_Shared)
firstSide = INT(REAL( myComputeNodeRank   *nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))+1
lastSide  = INT(REAL((myComputeNodeRank+1)*nComputeNodeSurfTotalSides)/REAL(nComputeNodeProcessors))
#else
firstSide = 1
lastSide  = nComputeNodeSurfTotalSides
#endif /*USE_MPI*/

! calculate the coverage from the sampled values (also required in the MPI=OFF case) and nullify the ChemSampWall array
! in the case of MPI, ChemSampWall is an associate to the _Shared variant
DO iSide = firstSide, lastSide
  GlobalSideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  locBCID = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,GlobalSideID))
  DO iSpec = 1, nSpecies
    IF (PartBound%LatticeVec(locBCID).GT.0.) THEN
      ! update the surface coverage (direct calculation of the number of surface atoms)
      ChemWallProp(iSpec,:,:,iSide) = MAX(0., ChemWallProp(iSpec,:,:,iSide) + ChemSampWall(iSpec,:,:,iSide) * PartBound%LatticeVec(locBCID)* &
      PartBound%LatticeVec(locBCID)/(PartBound%MolPerUnitCell(locBCID)*SurfSideArea(:,:,iSide)))
    ELSE
      ! update the surface coverage (calculation with a surface monolayer)
      ChemWallProp(iSpec,:,:,iSide) = MAX(0., ChemWallProp(iSpec,:,:,iSide) + ChemSampWall(iSpec,:,:,iSide) / &
      (10.**(19)*SurfSideArea(:,:,iSide)))
    END IF
  END DO
  ! calculate the heat flux on the surface subside
  IF (DSMC%CalcSurfaceVal.AND.(SamplingActive.OR.WriteMacroSurfaceValues)) THEN
    ChemWallProp(nSpecies+1,:,:,iSide) = ChemWallProp(nSpecies+1,:,:,iSide) + ChemSampWall(nSpecies+1,:,:,iSide)
  END IF
  ChemSampWall(:,:,:,iSide) = 0.0
END DO

#if USE_MPI
! ensure synchronization on compute node
CALL BARRIER_AND_SYNC(ChemSampWall_Shared_Win         ,MPI_COMM_SHARED)
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win         ,MPI_COMM_SHARED)
END ASSOCIATE
#endif /*USE_MPI*/

END SUBROUTINE SurfChemCoverage


#if USE_MPI
SUBROUTINE ExchangeSurfChemCoverage()
!===================================================================================================================================
!> 1) compute-node leaders communicate the calculated coverage to the halo sides
!> 2) compute-node leaders ensure synchronization of shared arrays on their node
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode
USE MOD_MPI_Shared              ,ONLY: BARRIER_AND_SYNC
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED,MPI_COMM_LEADERS_SURF
USE MOD_MPI_Shared_Vars         ,ONLY: nSurfLeaders,myComputeNodeRank,mySurfRank
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp_Shared, ChemWallProp_Shared_Win
USE MOD_Particle_Boundary_Vars  ,ONLY: nSurfSample
USE MOD_Particle_Boundary_Vars  ,ONLY: GlobalSide2SurfSide, SurfMapping
USE MOD_Particle_MPI_Vars       ,ONLY: SurfSendBuf,SurfRecvBuf
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iSpec,iProc,SideID,iPos,p,q
INTEGER                         :: MessageSize,iSurfSide,SurfSideID
INTEGER                         :: nValues, SurfChemSampSize
TYPE(MPI_Request)               :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

! communication of the coverage values to the halo region of compute nodes
! swap the receive/send array, send coverage values to the sides, where you would receive from
IF (myComputeNodeRank.EQ.0) THEN
  SurfChemSampSize = nSpecies + 1
  nValues = SurfChemSampSize*nSurfSample**2
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we would have sent to this leader node (thus SendSurfSides)
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nSendSurfSides * nValues
    CALL MPI_IRECV( SurfSendBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , RecvRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! build message
  DO iProc = 0,nSurfLeaders-1
    ! Ignore myself
    IF (iProc .EQ. mySurfRank) CYCLE

    ! Only assemble message if we would have received from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Nullify everything
    iPos = 0
    SurfRecvBuf(iProc)%content = 0.

    DO iSurfSide = 1,SurfMapping(iProc)%nRecvSurfSides
      ! Get the right side id through the receive global id mapping
      SideID     = SurfMapping(iProc)%RecvSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      ! Assemble message
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          DO iSpec =1, nSpecies+1
            SurfRecvBuf(iProc)%content(iPos+1) = ChemWallProp_Shared(iSpec,p,q,SurfSideID)
            iPos = iPos + 1
          END DO
        END DO ! p=0,nSurfSample
      END DO ! q=0,nSurfSample
    END DO ! iSurfSide=1,SurfMapping(iProc)%nRecvSurfSides
  END DO

  ! send message
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we are expecting sides from this leader node
    IF (SurfMapping(iProc)%nRecvSurfSides.EQ.0) CYCLE

    ! Message is sent on MPI_COMM_LEADERS_SURF, so rank is indeed iProc
    MessageSize = SurfMapping(iProc)%nRecvSurfSides * nValues
    CALL MPI_ISEND( SurfRecvBuf(iProc)%content                   &
                  , MessageSize                                  &
                  , MPI_DOUBLE_PRECISION                         &
                  , iProc                                        &
                  , 1209                                         &
                  , MPI_COMM_LEADERS_SURF                        &
                  , SendRequest(iProc)                           &
                  , IERROR)
  END DO ! iProc

  ! Finish received number of sampling surfaces
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    IF (SurfMapping(iProc)%nRecvSurfSides.NE.0) THEN
      CALL MPI_WAIT(SendRequest(iProc),MPI_STATUS_IGNORE,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPI_STATUS_IGNORE,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF
  END DO ! iProc

  ! add data do my list
  DO iProc = 0,nSurfLeaders-1
    ! ignore myself
    IF (iProc.EQ.mySurfRank) CYCLE

    ! Only open recv buffer if we would have sent this leader node
    IF (SurfMapping(iProc)%nSendSurfSides.EQ.0) CYCLE

    iPos=0
    DO iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      SideID     = SurfMapping(iProc)%SendSurfGlobalID(iSurfSide)
      SurfSideID = GlobalSide2SurfSide(SURF_SIDEID,SideID)
      ! Store values in the halo regions
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          DO iSpec =1, nSpecies+1
            ChemWallProp_Shared(iSpec,p,q,SurfSideID) = SurfSendBuf(iProc)%content(iPos+1)
            iPos = iPos + 1
          END DO
        END DO ! p=0,nSurfSample
      END DO ! q=0,nSurfSample
    END DO ! iSurfSide = 1,SurfMapping(iProc)%nSendSurfSides
      ! Nullify buffer
    SurfSendBuf(iProc)%content = 0.
  END DO ! iProc
END IF

! ensure synchronization on compute node
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win         ,MPI_COMM_SHARED)

END SUBROUTINE ExchangeSurfChemCoverage
#endif /*USE_MPI*/

END MODULE MOD_SurfaceModel_Chemistry
