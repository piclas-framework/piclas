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
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Products'  &
                                           ,'Products of Reaction[$] (Product1, Product2, Product3)', '0 , 0, 0' &
                                           , numberedmulti=.TRUE.)
 CALL prms%CreateRealOption(     'Surface-Reaction[$]-ReactHeat', &
                                    'Heat flux to or from the surface due to the reaction [K]', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-HeatScaling', &
                                    'Linear dependence of the heat flux on the coverage', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EnergyAccommodation', &
                                    'Energy accommodation coefficient', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Inhibition','Inhibition/Coadsorption behaviour due to other reactions', &
                                '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Promotion','Promotion/Coadsorption behaviour due to other reactions', &
                                '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-StickingCoefficient','Ratio of adsorbed to impinging particles on a\n' //&
                                'reactive surface, Langmuir or Kisluik model', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-DissOrder',  &
                                    'Associative = 1, dissociative = 2', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EqConstant',  &
                                    'Equilibrium constant between the adsorption and desorption (K), Langmuir: K=1', '1.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Reaction[$]-DissociativeAdsorption', 'Special adsorption case, only one half of the molecule is' //&
                                'adsorbed and the other remains in the gas-phase', '.FALSE.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-AdsorptionProduct','Species that stays adsorbed on the surface', &
                                    '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-GasPhaseProduct','Species that is desorbed into the gas-phase', &
                                '0', numberedmulti=.TRUE.)
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
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies, SpeciesDatabase
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, nPartBound
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChem, SurfChemReac, DoChemSurface
! USE MOD_Particle_Surfaces_Vars
USE MOD_io_hdf5
USE MOD_HDF5_input              ,ONLY:ReadAttribute, DatasetExists, AttributeExists
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
INTEGER(HID_T)        :: file_id_specdb                       ! File identifier
LOGICAL               :: DataSetFound
LOGICAL               :: Attr_Exists
CHARACTER(LEN=32)     :: hilf
INTEGER               :: iReac, iReac2, iBound, iVal, err
INTEGER               :: ReadInNumOfReact
INTEGER               :: iSpec, SpecID, ReactionPathPerSpecies(nSpecies), BCID
!===================================================================================================================================

IF(SurfChem%NumOfReact.LE.0) RETURN

ReadInNumOfReact = SurfChem%NumOfReact
LBWRITE(*,*) '| Number of considered reaction paths on Surfaces: ', SurfChem%NumOfReact
!----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(SurfChemReac(ReadInNumOfReact))
! Surface map
ALLOCATE(SurfChem%BoundIsChemSurf(nPartBound))
SurfChem%BoundIsChemSurf = .FALSE.
ALLOCATE(SurfChem%PSMap(nPartBound))
SurfChem%CatBoundNum = 0
DO iBound=1, nPartBound
  ALLOCATE(SurfChem%PSMap(iBound)%PureSurfReac(ReadInNumOfReact))
  SurfChem%PSMap(iBound)%PureSurfReac = .FALSE.
END DO

! Probability based surface chemistry model
ALLOCATE(SurfChem%EventProbInfo(nSpecies))
SurfChem%EventProbInfo(:)%NumOfReactionPaths = 0
ReactionPathPerSpecies = 0

! Get the reaction names, reactive species and boundaries
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  SurfChemReac(iReac)%Reactants(:)           = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
  SurfChemReac(iReac)%Products(:)            = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',3,'0,0,0')
  SurfChemReac(iReac)%ReactType             = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))
  SurfChemReac(iReac)%NumOfBounds           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries')
  IF (SurfChemReac(iReac)%NumOfBounds.EQ.0) THEN
    CALL abort(__STAMP__,'ERROR: At least one boundary must be defined for each surface reaction!',iReac)
  END IF

  SurfChemReac(iReac)%Boundaries = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Boundaries',SurfChemReac(iReac)%NumOfBounds)
  ! Define the surface model
  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))
  CASE('P')
    ! Simple probability based surface model
    PartBound%SurfaceModel(SurfChemReac(iReac)%Boundaries) = 2
    IF(SurfChemReac(iReac)%Reactants(2).NE.0) CALL abort(__STAMP__,' ERROR: Probability based model only supports one reactant!')
    SpecID = SurfChemReac(iReac)%Reactants(1)
    SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths = SurfChem%EventProbInfo(SpecID)%NumOfReactionPaths + 1
  CASE('A','D','LH','LHD','ER')
    SurfChemReac(iReac)%CatName              = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-SurfName'))
    ! Check if a surface model is already defined, if not set the boundary to reactive
    IF (ANY(PartBound%SurfaceModel(SurfChemReac(iReac)%Boundaries).GT.0)) THEN
      IF (ANY(PartBound%SurfaceModel(SurfChemReac(iReac)%Boundaries).NE.20)) THEN
        SWRITE(*,*) 'WARNING: The surface model for boundary ', SurfChemReac(iReac)%Boundaries, ' is set to catalytic.'
      END IF
    ELSE
      PartBound%SurfaceModel(SurfChemReac(iReac)%Boundaries) = 20
      PartBound%Reactive(SurfChemReac(iReac)%Boundaries) = .TRUE.
    END IF
    DO iReac2 = 1, SurfChemReac(iReac)%NumOfBounds
      SurfChem%BoundIsChemSurf(SurfChemReac(iReac)%Boundaries(iReac2)) = .TRUE.
    END DO
    DoChemSurface = .TRUE.
    ! Select pure surface reactions
    DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
      iBound = SurfChemReac(iReac)%Boundaries(iVal)
      SurfChem%PSMap(iBound)%PureSurfReac(iReac) = .TRUE.
    END DO
  CASE DEFAULT
    SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac(iReac)%ReactType)
    CALL abort(__STAMP__,' ERROR: Surface Reaction Type does not exist!')
  END SELECT
END DO

! Determine the number of boundaries with a surface reaction with a surface flux on them
DO iBound = 1, nPartBound
  IF (SurfChem%BoundIsChemSurf(iBound)) THEN
    SurfChem%CatBoundNum = SurfChem%CatBoundNum + 1
  END IF
END DO

DO iSpec = 1, nSpecies
  IF(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths.GT.0) THEN
    ! Allocate the species specific type with the number of the possible reaction paths
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ReactionIndex(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ReactionIndex = 0
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ReactionProb(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ReactionProb = 0.
    ALLOCATE(SurfChem%EventProbInfo(iSpec)%ProdTransACC(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths))
    SurfChem%EventProbInfo(iSpec)%ProdTransACC = -1.
  END IF
END DO

! Bulk species involved in the reactions
SurfChem%SurfSpecies                   = GETINT('Surface-Species','0')

! Diffusion
SurfChem%Diffusion                     = GETLOGICAL('Surface-Diffusion', '.FALSE.')
SurfChem%TotDiffusion                  = GETLOGICAL('Surface-TotalDiffusion', '.FALSE.')

! SpeciesDatabase
SurfChem%OverwriteCatParameters        = GETLOGICAL('OverwriteCatParameters', '.FALSE.')

IF (SpeciesDatabase.EQ.'none') THEN
  SurfChem%OverwriteCatParameters = .TRUE.
END IF

IF(SpeciesDatabase.NE.'none') THEN
  CALL H5OPEN_F(err)

  CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)

  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    dsetname = TRIM('/Surface-Chemistry/'//TRIM(SurfChemReac(iReac)%CatName))

    CALL DatasetExists(file_id_specdb,TRIM(dsetname),DataSetFound)
    IF(.NOT.DataSetFound)THEN
      SurfChem%OverwriteCatParameters = .TRUE.
      SWRITE(*,*) 'WARNING: DataSet not found: ['//TRIM(dsetname)//'] ['//TRIM(SpeciesDatabase)//']'
    ELSE
      CALL ReadAttribute(file_id_specdb,'Type',1,DatasetName = dsetname,StrScalar=SurfChemReac(iReac)%ReactType)
      CALL AttributeExists(file_id_specdb,'Inhibition',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'Inhibition',1,DatasetName = dsetname,IntScalar=SurfChemReac(iReac)%Inhibition)
      ELSE
        SurfChemReac(iReac)%Inhibition= 0
      END IF
      CALL AttributeExists(file_id_specdb,'Promotion',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'Promotion',1,DatasetName = dsetname,IntScalar=SurfChemReac(iReac)%Promotion)
      ELSE
        SurfChemReac(iReac)%Promotion= 0
      END IF
      CALL AttributeExists(file_id_specdb,'ReactHeat',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'ReactHeat',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%EReact)
      ELSE
        SurfChemReac(iReac)%EReact= 0.
      END IF
      CALL AttributeExists(file_id_specdb,'HeatScaling',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'HeatScaling',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%EScale)
      ELSE
        SurfChemReac(iReac)%EScale= 0.
      END IF
      CALL AttributeExists(file_id_specdb,'EnergyAccommodation',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'EnergyAccommodation',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%HeatAccommodation)
      ELSE
        SurfChemReac(iReac)%HeatAccommodation= 0.
      END IF

      SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))
      CASE('A')
        SurfChemReac(iReac)%DissociativeAds = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-DissociativeAdsorption', '.FALSE.')
        CALL AttributeExists(file_id_specdb,'StickingCoefficient',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'StickingCoefficient',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%S_initial)
        ELSE
          SurfChemReac(iReac)%S_initial= 1.
        END IF
        CALL AttributeExists(file_id_specdb,'EqConstant',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'EqConstant',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%EqConstant)
        ELSE
          SurfChemReac(iReac)%EqConstant= 1.
        END IF
        CALL AttributeExists(file_id_specdb,'DissOrder',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'DissOrder',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%DissOrder)
        ELSE
          SurfChemReac(iReac)%DissOrder= 1.
        END IF
        ! Special case: dissociative adsorption
        IF (SurfChemReac(iReac)%DissociativeAds) THEN
          CALL AttributeExists(file_id_specdb,'AdsorptionProduct',TRIM(dsetname), AttrExists=Attr_Exists)
          IF (Attr_Exists) THEN
            CALL ReadAttribute(file_id_specdb,'AdsorptionProduct',1,DatasetName = dsetname,IntScalar=SurfChemReac(iReac)%AdsorbedProduct)
          ELSE
            CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
          END IF
          CALL AttributeExists(file_id_specdb,'GasPhaseProduct',TRIM(dsetname), AttrExists=Attr_Exists)
          IF (Attr_Exists) THEN
            CALL ReadAttribute(file_id_specdb,'GasPhaseProduct',1,DatasetName = dsetname,IntScalar=SurfChemReac(iReac)%GasProduct)
          ELSE
            CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
          END IF
        END IF

      CASE('D')
        CALL AttributeExists(file_id_specdb,'LateralInteraction',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'LateralInteraction',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%W_interact)
        ELSE
          SurfChemReac(iReac)%W_interact= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Ca',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Ca',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%C_a)
        ELSE
          SurfChemReac(iReac)%C_a= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Cb',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Cb',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%C_b)
        ELSE
          SurfChemReac(iReac)%C_b= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%Prefactor)
        ELSE
          SurfChemReac(iReac)%Prefactor= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%E_initial)
        ELSE
          SurfChemReac(iReac)%E_initial= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'DissOrder',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'DissOrder',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%DissOrder)
        ELSE
          SurfChemReac(iReac)%DissOrder= 1.
        END IF

      CASE('LH')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%ArrheniusEnergy)
        ELSE
          SurfChemReac(iReac)%ArrheniusEnergy= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%Prefactor)
        ELSE
          SurfChemReac(iReac)%Prefactor= 1.
        END IF

      CASE('LHD')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%ArrheniusEnergy)
        ELSE
          SurfChemReac(iReac)%ArrheniusEnergy= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%Prefactor)
        ELSE
          SurfChemReac(iReac)%Prefactor= 1.
        END IF

      CASE('ER')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%ArrheniusEnergy)
        ELSE
          SurfChemReac(iReac)%ArrheniusEnergy = 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac(iReac)%Prefactor)
        ELSE
          SurfChemReac(iReac)%Prefactor = 1.
        END IF

      CASE DEFAULT
        SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac(iReac)%ReactType)
        CALL abort(__STAMP__,'Surface Reaction Type does not exist')
      END SELECT

    END IF !DatasetFound
  END DO !iReac
  ! Close the file.
  CALL H5FCLOSE_F(file_id_specdb, err)
  ! Close FORTRAN interface.
  CALL H5CLOSE_F(err)
END IF !SpeciesDatabase


IF (SurfChem%OverwriteCatParameters) THEN
  ! Loop over the surface reactions
  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    SurfChemReac(iReac)%EReact                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ReactHeat','0.')
    SurfChemReac(iReac)%EScale                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-HeatScaling','0.')
    SurfChemReac(iReac)%HeatAccommodation     = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EnergyAccommodation','1.')

    SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))
    CASE('A')
      SurfChemReac(iReac)%Inhibition          = GETINT('Surface-Reaction'//TRIM(hilf)//'-Inhibition','0')
      SurfChemReac(iReac)%Promotion           = GETINT('Surface-Reaction'//TRIM(hilf)//'-Promotion','0')
      SurfChemReac(iReac)%S_initial           = GETREAL('Surface-Reaction'//TRIM(hilf)//'-StickingCoefficient','1.')
      SurfChemReac(iReac)%EqConstant          = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EqConstant','1.')
      SurfChemReac(iReac)%DissOrder           = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.')
      SurfChemReac(iReac)%DissociativeAds     = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-DissociativeAdsorption', '.FALSE.')
      ! Special case of the dissociative adsorption, half of the molecule is desorbed back into the gas-phase
      IF (SurfChemReac(iReac)%DissociativeAds) THEN
        SurfChemReac(iReac)%AdsorbedProduct   = GETINT('Surface-Reaction'//TRIM(hilf)//'-AdsorptionProduct','0')
        SurfChemReac(iReac)%GasProduct        = GETINT('Surface-Reaction'//TRIM(hilf)//'-GasPhaseProduct','0')
        IF ((SurfChemReac(iReac)%GasProduct.EQ.0).OR.(SurfChemReac(iReac)%GasProduct.EQ.0)) THEN
          CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
        END IF
      END IF

    CASE('D')
      SurfChemReac(iReac)%W_interact = GETREAL('Surface-Reaction'//TRIM(hilf)//'-LateralInteraction','0.')
      SurfChemReac(iReac)%C_a = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Ca','0.')
      SurfChemReac(iReac)%C_b = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Cb','0.')
      SurfChemReac(iReac)%Prefactor = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','0.')
      SurfChemReac(iReac)%E_initial = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac(iReac)%DissOrder = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.')

    CASE('LH')
      SurfChemReac(iReac)%ArrheniusEnergy = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac(iReac)%Prefactor = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE('LHD')
      SurfChemReac(iReac)%ArrheniusEnergy = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac(iReac)%Prefactor = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE('ER')
      SurfChemReac(iReac)%ArrheniusEnergy = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac(iReac)%Prefactor = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE('P')
      SpecID = SurfChemReac(iReac)%Reactants(1)
      ReactionPathPerSpecies(SpecID) = ReactionPathPerSpecies(SpecID) + 1
      SurfChem%EventProbInfo(SpecID)%ReactionIndex(ReactionPathPerSpecies(SpecID)) = iReac
      SurfChem%EventProbInfo(SpecID)%ReactionProb(ReactionPathPerSpecies(SpecID)) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EventProbability')
      SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ProductAccommodation')
      ! Sanity checks
      IF(SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).NE.-1.) THEN
        ! If a reaction-specific accommodation coefficient is used, check if it is between 0 and 1
        IF ((SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).LT.0.).OR. &
            (SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).GT.1.)) THEN
          CALL abort(__STAMP__,'Reaction-specific thermal accommodation must be between 0 and 1 for reaction: ', IntInfoOpt=iReac)
        END IF
        ! If the reaction-specific accommodation coefficient is greater than 0, check if a wall temperature has been defined
        IF(SurfChem%EventProbInfo(SpecID)%ProdTransACC(ReactionPathPerSpecies(SpecID)).GT.0.) THEN
          DO iVal = 1, SurfChemReac(iReac)%NumOfBounds
            BCID = SurfChemReac(iReac)%Boundaries(iVal)
            IF(PartBound%WallTemp(BCID).EQ.0.) THEN
              CALL abort(__STAMP__,'Reaction-specific thermal accommodation requires a wall temperature for boundary '//&
                        TRIM(PartBound%SourceBoundName(BCID))//' used for reaction: ', IntInfoOpt=iReac)
            END IF
          END DO
        END IF
      END IF
    END SELECT
  END DO
END IF

! Sanity check: Total reaction probability of a single species at the surface must not be above 1
DO iSpec = 1, nSpecies
  IF(SurfChem%EventProbInfo(iSpec)%NumOfReactionPaths.GT.0) THEN
    IF(SUM(SurfChem%EventProbInfo(iSpec)%ReactionProb(:)).GT.1.) THEN
      CALL abort(__STAMP__,'ERROR: Total probability above unity for species: ', IntInfoOpt=iSpec)
    END IF
  END IF
END DO

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

ALLOCATE(ChemSampWall(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemSampWall = 0.0
ALLOCATE(ChemDesorpWall(1:nSpecies,1,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemDesorpWall = 0.0

#if USE_MPI
CALL Allocate_Shared((/nSpecies,2,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),ChemSampWall_Shared_Win,ChemSampWall_Shared)
CALL MPI_WIN_LOCK_ALL(0,ChemSampWall_Shared_Win,IERROR)
IF (myComputeNodeRank.EQ.0) THEN
  ChemSampWall_Shared = 0.
END IF
CALL BARRIER_AND_SYNC(ChemSampWall_Shared_Win,MPI_COMM_SHARED)

CALL Allocate_Shared((/nSpecies,2,nSurfSample,nSurfSample,nComputeNodeSurfTotalSides/),ChemWallProp_Shared_Win,ChemWallProp_Shared)
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
      ChemWallProp(iSpec,1,:,:,iSide) = PartBound%CoverageIni(iBC, iSpec)
    END DO
  END DO
END IF
CALL BARRIER_AND_SYNC(ChemWallProp_Shared_Win,MPI_COMM_SHARED)
#else
ALLOCATE(ChemWallProp(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemWallProp = 0.0
DO iSide = 1, nComputeNodeSurfTotalSides
  ! get global SideID. This contains only nonUniqueSide, no special mortar treatment required
  SideID = SurfSide2GlobalSide(SURF_SIDEID,iSide)
  iBC = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
  DO iSpec = 1, nSpecies
  ! Initial surface coverage
    ChemWallProp(iSpec,1,:,:,iSide) = PartBound%CoverageIni(iBC, iSpec)
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
USE MOD_Globals                   ,ONLY: abort,UNITVECTOR,OrthoNormVec
USE MOD_DSMC_PolyAtomicModel      ,ONLY: DSMC_SetInternalEnr
USE MOD_part_operations           ,ONLY: RemoveParticle, CreateParticle
USE MOD_part_tools                ,ONLY: VeloFromDistribution, GetParticleWeight
USE MOD_SurfaceModel_Tools        ,ONLY: MaxwellScattering, CalcPostWallCollVelo
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
! VARIABLES
USE MOD_Globals_Vars              ,ONLY: PI, BoltzmannConst
USE MOD_Particle_Vars             ,ONLY: PartSpecies,Species,usevMPF, WriteMacroSurfaceValues
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, GlobalSide2SurfSide, SurfSideArea
USE MOD_SurfaceModel_Vars         ,ONLY: SurfChem, SurfChemReac , ChemWallProp, ChemSampWall
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared, BoundsOfElem_Shared
USE MOD_DSMC_Vars                 ,ONLY: DSMC, RadialWeighting, SamplingActive
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: n_loc(1:3)
INTEGER,INTENT(IN) :: PartID, SideID
INTEGER,INTENT(IN) :: GlobalElemID        !< Global element ID of the particle impacting the surface
REAL,INTENT(IN)    :: PartPosImpact(1:3)  !< Charge and position of impact of bombarding particle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ProductSpecNbr   !< number of emitted particles for ProductSpec(2)
INTEGER            :: locBCID, SurfSideID
CHARACTER(LEN=5)   :: InteractionType
REAL               :: RanNum, RanNum2
REAL               :: Coverage, MaxCoverage, TotalCoverage, Theta, MaxTotalCov
REAL               :: CoAds_Coverage, CoAds_MaxCov
REAL               :: S_0, StickCoeff
REAL               :: EqConstant, DissOrder
REAL               :: WallTemp
REAL               :: nu, E_act, Rate, Prob, Prob_new, Prob_Scaled
REAL               :: NewPos(1:3)
REAL               :: SurfMol, AdCountIter, AdsDens
REAL               :: tang1(1:3), tang2(1:3), WallVelo(1:3), BoundsOfElemCenter(1:3), NewVelo(3)
REAL               :: AdsHeat, ReacHeat, BetaCoeff
REAL               :: partWeight
REAL,PARAMETER     :: eps=1e-6
REAL,PARAMETER     :: eps2=1.0-eps
INTEGER            :: speciesID
INTEGER            :: iReac_Ads, iReac_ER, iReac_ER_new
INTEGER            :: iReac, iValProd, iProd, iReactant, iValReac
INTEGER            :: iCoadsReac, iCoadsSpec
INTEGER            :: NewPartID
INTEGER            :: SubP, SubQ
!===================================================================================================================================
! 0.) Determine the surface parameters: Coverage and number of surface molecules
locBCID     = PartBound%MapToPartBC(SideInfo_Shared(SIDE_BCID,SideID))
SurfSideID  = GlobalSide2SurfSide(SURF_SIDEID,SideID)
ProductSpecNbr = 0
SubP = TrackInfo%p
SubQ = TrackInfo%q
InteractionType = 'None'
StickCoeff = 0.0
Prob = 0.0
iReac_ER = 0
speciesID = PartSpecies(PartID)
! MacroParticleFactor
partWeight = GetParticleWeight(PartID)
IF(.NOT.(usevMPF.OR.RadialWeighting%DoRadialWeighting)) THEN
  partWeight = partWeight * Species(speciesID)%MacroParticleFactor
END IF

IF(PartBound%LatticeVec(locBCID).GT.0.) THEN
  ! Number of surface molecules in dependence of the occupancy of the unit cell
  SurfMol = PartBound%MolPerUnitCell(locBCID) * SurfSideArea(SubP, SubQ,SurfSideID) &
                              /(PartBound%LatticeVec(locBCID)*PartBound%LatticeVec(locBCID))
ELSE
  ! Alternative calculation by the average number of surface molecules per area for a monolayer
  SurfMol = 10.**19 * SurfSideArea(SubP, SubQ,SurfSideID)
END IF ! LatticeVec.GT.0

DO iReac = 1, SurfChem%NumOfReact
  SELECT CASE (TRIM(SurfChemReac(iReac)%ReactType))

  ! 1.) Calculate the sticking coefficient by the Kisliuk model (adsorption)
  CASE('A')
    IF(ANY(SurfChemReac(iReac)%Reactants(:).EQ.speciesID)) THEN
      iReac_Ads = iReac

      ! Absolute coverage in terms of the number of surface molecules
      ! Special case: dissociative adsorption
      IF (SurfChemReac(iReac)%DissociativeAds) THEN
        iProd = SurfChemReac(iReac)%AdsorbedProduct
        Coverage = ChemWallProp(iProd,1,SubP,SubQ,SurfSideID)
      ELSE
        IF(ANY(SurfChemReac(iReac)%Products(:).NE.0)) THEN
          DO iValProd=1, SIZE(SurfChemReac(iReac)%Products(:))
            IF(SurfChemReac(iReac)%Products(iValProd).NE.0) THEN
              iProd = SurfChemReac(iReac)%Products(iValProd)
              Coverage = ChemWallProp(iProd,1,SubP,SubQ,SurfSideID)
            END IF
          END DO
        ELSE
          Coverage = ChemWallProp(speciesID,1,SubP,SubQ,SurfSideID)
        END IF
      END IF ! Dissociative ADS

      ! Definition of the variables
      MaxCoverage = PartBound%MaxCoverage(locBCID,speciesID)
      TotalCoverage = SUM(ChemWallProp(:,1,SubP, SubQ, SurfSideID))
      MaxTotalCov = PartBound%MaxTotalCoverage(locBCID)
      DissOrder = SurfChemReac(iReac)%DissOrder
      S_0 = SurfChemReac(iReac)%S_initial
      EqConstant = SurfChemReac(iReac)%EqConstant
      StickCoeff = SurfChemReac(iReac)%StickCoeff

      ! Determine the heat of adsorption in dependence of the coverage [J]
      AdsHeat = (SurfChemReac(iReac)%EReact - Coverage * SurfChemReac(iReac)%EScale) * BoltzmannConst

      ! Theta = free surface sites required for the adsorption
      ! Determination of possible coadsorption processes
      IF(SurfChemReac(iReac)%Inhibition.NE.0) THEN
        iCoadsReac = SurfChemReac(iReac)%Inhibition
        iCoadsSpec = SurfChemReac(iCoadsReac)%Reactants(1)
        CoAds_Coverage = ChemWallProp(iCoadsSpec,1,SubP, SubQ, SurfSideID)
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        Theta = 1.0 - Coverage/MaxCoverage - CoAds_Coverage/CoAds_MaxCov
      ELSE IF(SurfChemReac(iReac)%Promotion.NE.0) THEN
        iCoadsReac = SurfChemReac(iReac)%Promotion
        iCoadsSpec = SurfChemReac(iCoadsReac)%Reactants(1)
        CoAds_Coverage = ChemWallProp(iCoadsSpec,1,SubP, SubQ, SurfSideID)
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        Theta = 1.0 - Coverage/MaxCoverage + CoAds_Coverage/CoAds_MaxCov
      ELSE
        Theta = 1.0 - Coverage/MaxCoverage
      END IF

      ! Check whether the maximum coverage value is reached:
      IF(Theta.GT.0.0 .AND. TotalCoverage.LT.MaxTotalCov) THEN
        Theta = Theta**DissOrder
      ! Kisliuk model (for EqConstant=1 and MaxCoverage=1: Langmuir model)
        StickCoeff = S_0 * (1.0 + EqConstant * (1.0/Theta - 1.0))**(-1.0)
      ELSE
        StickCoeff = 0.0
      END IF

    END IF

  ! 2.) Calculate the reaction probability by the Arrhenius equation (bias-free for multiple channels)
  CASE('ER')
    IF(ANY(SurfChemReac(iReac)%Reactants(:).EQ.speciesID)) THEN

      ! Definition of the variables
      WallTemp = PartBound%WallTemp(locBCID)
      nu = SurfChemReac(iReac)%Prefactor
      E_act = SurfChemReac(iReac)%ArrheniusEnergy
      Rate = SurfChemReac(iReac)%Rate
      BetaCoeff = SurfChemReac(iReac)%HeatAccommodation

      ! Check for the coverage values of the reactant adsorbed on the surface
      IF(ANY(SurfChemReac(iReac)%Reactants(:).NE.speciesID)) THEN
        DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
          IF(SurfChemReac(iReac)%Reactants(iValReac).NE.speciesID .AND. SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
            iReactant = SurfChemReac(iReac)%Reactants(iValReac)
            IF(iReactant.NE.SurfChem%SurfSpecies) THEN
              Coverage = ChemWallProp(iReactant,1,SubP, SubQ, SurfSideID)
              AdCountIter = ChemSampWall(iReactant, 1,SubP,SubQ, SurfSideID)
            ELSE  ! Involvement of the bulk species
              Coverage = 1.
              AdCountIter = 0.
            END IF
          END IF
        END DO
      ELSE
        Coverage = ChemWallProp(speciesID,1,SubP, SubQ, SurfSideID)
        AdCountIter = ChemSampWall(speciesID, 1,SubP,SubQ, SurfSideID)
      END IF
      ! Absolute particle density of the element
      AdsDens = Coverage * SurfMol / SurfSideArea(SubP, SubQ,SurfSideID)

      ! Determine the reaction heat in dependence of the coverage [J]
      ReacHeat = (SurfChemReac(iReac)%EReact - Coverage * SurfChemReac(iReac)%EScale) * BoltzmannConst

      ! Bias free calculation for multiple reaction channels
      IF(iReac_ER.EQ.0) THEN
        iReac_ER = iReac
        Rate = nu * AdsDens * exp(-E_act/WallTemp) ! Energy in K
        Prob = SQRT(2.*PI*Species(speciesID)%MassIC/(BoltzmannConst*WallTemp)) * Rate

        ! Comparison of adsorbate numbers to reactant particle weights
        IF(partWeight.GT.(Coverage*SurfMol)) THEN
          Prob = 0.0
        ! Test for the changes during the iteration
        ELSE IF ((Coverage*SurfMol + AdCountIter - partWeight).LT.0.0) THEN
          Prob = 0.0
        END IF

      ELSE ! iReac.NE.0
        iReac_ER_new = iReac
        Rate = nu * AdsDens * exp(-E_act/WallTemp) ! Energy in K
        Prob_new = SQRT(2.*PI*Species(speciesID)%MassIC/(BoltzmannConst*WallTemp)) * Rate

        ! Comparison of adsorbate numbers to reactant particle weights
        IF(partWeight.GT.(Coverage*SurfMol)) THEN
          Prob_new = 0.0
        ! Test for the changes during the iteration
        ELSE IF ((Coverage*SurfMol + AdCountIter - partWeight).LT.0.0) THEN
          Prob_new = 0.0
        END IF

        ! determine most likely reaction channel
        CALL RANDOM_NUMBER(RanNum)

        IF(Prob_new.GT.Prob) THEN
          iReac_ER = iReac_ER_new
          Prob = Prob_new
        ELSE IF(Prob_new.EQ.Prob) THEN
          IF(RanNum.GT.0.5) THEN
            iReac_ER = iReac_ER_new
            Prob = Prob_new
          END IF
        END IF
      END IF !iReac.EQ.0
    END IF !iReac.EQ.speciesID

  CASE DEFAULT
  END SELECT !Surface reaction
END DO !iReac

! ----------------------------------------------------------------------------------------------------------------------------------
! 3.) Choose the occurring pathway by comparison with a random number
! Rescale the probability (ER-Reaction) and the sticking coefficient (adsorption)
IF ((Prob+StickCoeff).GT.0.) THEN
  CALL RANDOM_NUMBER(RanNum)
  CALL RANDOM_NUMBER(RanNum2)
  Prob_Scaled = Prob/(Prob + StickCoeff)
  IF(Prob_Scaled.GT.RanNum) THEN
    IF(Prob.GT.RanNum2) THEN
      InteractionType = 'ER'
      iReac = iReac_ER
    END IF
  ELSE
    IF (StickCoeff.GT.RanNum2) THEN
      InteractionType = 'A'
      iReac = iReac_Ads
    END IF
  END IF
END IF

! ----------------------------------------------------------------------------------------------------------------------------------
! 4.) Perform the chosen process
SELECT CASE(TRIM(InteractionType))
! 4a.) Adsorption: delete the incoming particle and update the surface values
CASE('A')
  ! Dissociative adsorption
  IF (SurfChemReac(iReac)%DissociativeAds) THEN
    CALL RemoveParticle(PartID)

    ! Heat flux on the surface created by the adsorption
    ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) + AdsHeat * partWeight
    ! Update the number of adsorbed molecules by binding half of the molecule to the surface
    iProd = SurfChemReac(iReac)%AdsorbedProduct
    ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) = ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) + DissOrder * partWeight

    ! Re-insert the other half of the molecule into the gas-phase
    WallVelo = PartBound%WallVelo(1:3,locBCID)
    CALL OrthoNormVec(n_loc,tang1,tang2)

    ! Get Elem Center
    BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobalElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,2,GlobalElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,3,GlobalElemID)) /) / 2.

    iProd = SurfChemReac(iReac)%GasProduct

    NewVelo(1:3) = CalcPostWallCollVelo(iProd,0.,WallTemp,BetaCoeff)
    NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
    NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*PartPosImpact(1:3)

    CALL CreateParticle(iProd,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=partWeight)

    CALL DSMC_SetInternalEnr(iProd,locBCID,NewPartID,4)

    IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
    CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)

  ELSE
    CALL RemoveParticle(PartID)

    ! Heat flux on the surface created by the adsorption
    ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) + AdsHeat * partWeight
    ! Update the number of adsorbed molecules
    IF(ANY(SurfChemReac(iReac)%Products(:).NE.0)) THEN
      DO iValProd=1, SIZE(SurfChemReac(iReac)%Products(:))
        IF(SurfChemReac(iReac)%Products(iValProd).NE.0) THEN
          iProd = SurfChemReac(iReac)%Reactants(iValProd)
          ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) = ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) + DissOrder * partWeight
        END IF
      END DO
    ELSE
      ChemSampWall(speciesID, 1,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 1,SubP,SubQ, SurfSideID) + DissOrder * partWeight
    END IF
  END IF

! 4b.) ER: delete the incoming particle, update the surface values and create the gas phase products
CASE('ER')
  CALL RemoveParticle(PartID)

  ! Heat flux on the surface created by the reaction
  ChemSampWall(speciesID,2,SubP,SubQ,SurfSideID) = ChemSampWall(speciesID,2,SubP,SubQ,SurfSideID) + ReacHeat*partWeight*BetaCoeff

  ! Create the Eley-Rideal reaction product
  ! Incomplete energy accommodation: remaining energy is added to the product
  WallVelo = PartBound%WallVelo(1:3,locBCID)
  CALL OrthoNormVec(n_loc,tang1,tang2)

  ! Get Elem Center
  BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobalElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,2,GlobalElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,3,GlobalElemID)) /) / 2.

  DO iValProd=1, SIZE(SurfChemReac(iReac)%Products(:))
    IF(SurfChemReac(iReac)%Products(iValProd).NE.0) THEN
      iProd = SurfChemReac(iReac)%Products(iValProd)

      NewVelo(1:3) = CalcPostWallCollVelo(iProd,0.,WallTemp,BetaCoeff)
      NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
      NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*PartPosImpact(1:3)

      CALL CreateParticle(iProd,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=partWeight)

      CALL DSMC_SetInternalEnr(iProd,locBCID,NewPartID,4)

      ! Sampling of newly created particles
      IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
      CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
    END IF
  END DO

  ! Update the number of adsorbed molecules
  IF(ANY(SurfChemReac(iReac)%Reactants(:).NE.speciesID)) THEN
    DO iValReac=1, SIZE(SurfChemReac(iReac)%Reactants(:))
      IF(SurfChemReac(iReac)%Reactants(iValReac).NE.speciesID .AND. SurfChemReac(iReac)%Reactants(iValReac).NE.0) THEN
        iReactant = SurfChemReac(iReac)%Reactants(iValReac)
        IF(iReactant.NE.SurfChem%SurfSpecies) THEN
          ChemSampWall(iReactant, 1,SubP,SubQ, SurfSideID) = ChemSampWall(iReactant, 1,SubP,SubQ, SurfSideID) - partWeight
        END IF
      END IF
    END DO
  ELSE
    ChemSampWall(speciesID, 1,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 1,SubP,SubQ, SurfSideID) - partWeight
  END IF

CASE DEFAULT
  CALL MaxwellScattering(PartID,SideID,n_Loc)
END SELECT !Interaction Type

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
USE MOD_SurfaceModel_Tools        ,ONLY: SurfaceModelParticleEmission, MaxwellScattering, CalcPostWallCollVelo, CalcRotWallVelo
USE MOD_SurfaceModel_Tools        ,ONLY: SurfaceModelEnergyAccommodation, DiffuseReflection
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
    TotalProb = TotalProb + SurfChem%EventProbInfo(SpecID)%ReactionProb(iPath)
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
        CALL CreateParticle(ProdSpecID,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=OldMPF)
      ELSE
        ! New particle acquires the MPF of the new particle species
        CALL CreateParticle(ProdSpecID,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID)
      END IF ! usevMPF
      ! Adding the energy that is transferred from the surface onto the internal energies of the particle
      CALL SurfaceModelEnergyAccommodation(NewPartID,locBCID,WallTemp)
      ! Sampling of newly created particles
      IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
        CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
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
USE MOD_Particle_Vars           ,ONLY: nSpecies
USE MOD_SurfaceModel_Vars       ,ONLY: ChemWallProp
USE MOD_Particle_Boundary_vars  ,ONLY: SurfSideArea, SurfSide2GlobalSide
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: SurfTotalSideOnNode, nComputeNodeSurfTotalSides
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
      ChemWallProp(iSpec,1,:,:,iSide) = ChemWallProp(iSpec,1,:,:,iSide) + ChemSampWall(iSpec,1,:,:,iSide) * PartBound%LatticeVec(locBCID)* &
      PartBound%LatticeVec(locBCID)/(PartBound%MolPerUnitCell(locBCID)*SurfSideArea(:,:,iSide))
    ELSE
      ! update the surface coverage (calculation with a surface monolayer)
      ChemWallProp(iSpec,1,:,:,iSide) = ChemWallProp(iSpec,1,:,:,iSide) + ChemSampWall(iSpec,1,:,:,iSide) / &
      (10.**(19)*SurfSideArea(:,:,iSide))
    END IF
    ! calculate the heat flux on the surface subside
    ChemWallProp(iSpec,2,:,:,iSide) = ChemWallProp(iSpec,2,:,:,iSide) + ChemSampWall(iSpec,2,:,:,iSide)
  END DO
  ChemSampWall(:,:,:,:,iSide) = 0.0
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
INTEGER                         :: nValues, SurfChemVarNum, SurfChemSampSize
INTEGER                         :: RecvRequest(0:nSurfLeaders-1),SendRequest(0:nSurfLeaders-1)
!===================================================================================================================================
! nodes without sampling surfaces do not take part in this routine
IF (.NOT.SurfTotalSideOnNode) RETURN

! communication of the coverage values to the halo region of compute nodes
! swap the receive/send array, send coverage values to the sides, where you would receive from
IF (myComputeNodeRank.EQ.0) THEN
  SurfChemVarNum   = 2
  SurfChemSampSize = SurfChemVarNum * nSpecies
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
          DO iSpec =1, nSpecies
            SurfRecvBuf(iProc)%content(iPos+1:iPos+SurfChemVarNum) = ChemWallProp_Shared(iSpec,:,p,q,SurfSideID)
            iPos = iPos + SurfChemVarNum
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
      CALL MPI_WAIT(SendRequest(iProc),MPIStatus,IERROR)
      IF (IERROR.NE.MPI_SUCCESS) CALL ABORT(__STAMP__,' MPI Communication error',IERROR)
    END IF

    IF (SurfMapping(iProc)%nSendSurfSides.NE.0) THEN
      CALL MPI_WAIT(RecvRequest(iProc),MPIStatus,IERROR)
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
          DO iSpec =1, nSpecies
            ChemWallProp_Shared(iSpec,:,p,q,SurfSideID) = SurfSendBuf(iProc)%content(iPos+1:iPos+SurfChemVarNum)
            iPos = iPos + SurfChemVarNum
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
