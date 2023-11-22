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
PUBLIC :: DefineParametersSurfaceChemistry, InitializeVariablesSurfaceChemistry, SurfaceModel_Chemistry_Init
PUBLIC :: SurfaceModelChemistry
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
                                'No default, options are A (adsorption), D (desorption), ER (Eley-Rideal), LH (Langmuir-Hinshelwood), LHD', &
                                 numberedmulti=.TRUE.)
CALL prms%CreateStringOption(   'Surface-Reaction[$]-SurfName', 'none' ,numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'OverwriteCatParameters', 'Flag to set catalytic parameters manually', '.FALSE.')
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2)', '0 , 0' , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Products'  &
                                           ,'Products of Reaction[$] (Product1, Product2)', '0 , 0' &
                                           , numberedmulti=.TRUE.)
 CALL prms%CreateRealOption(     'Surface-Reaction[$]-ReactHeat', &
                                    'Heat flux to or from the surface due to the reaction', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-HeatScaling', &
                                    'Linear dependence of the heat flux on the coverage', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EnergyAccommodation', &
                                    'Energy accommodation coefficient', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Inhibition','Inhibition/Coadsorption behaviour due to other reactions', &
                                '0', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Promotion','Promotion/Coadsorption behaviour due to other reactions', &
                                '0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-StickingCoefficient', &
                                    'Ratio of adsorbed to impinging particles on a reactive surface', '0.' , numberedmulti=.TRUE.)
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
                                    'Arrhenius energy for the reaction/desorption', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'Surface-Diffusion', 'Diffusion along the surface', '.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-TotalDiffusion', 'Diffusion along all possible surface', '.FALSE.')
CALL prms%CreateIntOption(      'Surface-Reaction[$]-NumOfBoundaries', 'Num of boundaries for surface reaction.', &
                                    numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Boundaries', 'Array of boundary indices of surface reaction.', &
                                    numberedmulti=.TRUE., no=0)
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
USE MOD_SurfaceModel_Vars       ,ONLY: SurfChemReac, DoChemSurface
! USE MOD_Particle_Surfaces_Vars
USE MOD_io_hdf5
USE MOD_HDF5_input,         ONLY:ReadAttribute, DatasetExists, AttributeExists
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
REAL, ALLOCATABLE     :: StoichCoeff(:,:)
!===================================================================================================================================

IF(SurfChemReac%NumOfReact.LE.0) THEN
  RETURN
END IF
ReadInNumOfReact = SurfChemReac%NumOfReact
IF(ReadInNumOfReact.GT.0) THEN
  DoChemSurface = .TRUE.
END IF
SWRITE(*,*) '| Number of considered reaction paths on Surfaces: ', SurfChemReac%NumOfReact
!----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(StoichCoeff(nSpecies,2))
ALLOCATE(SurfChemReac%ReactType(SurfChemReac%NumOfReact))
SurfChemReac%ReactType = '0'
ALLOCATE(SurfChemReac%CatName(SurfChemReac%NumOfReact))
SurfChemReac%CatName = '0'
ALLOCATE(SurfChemReac%Reactants(SurfChemReac%NumOfReact,2))
SurfChemReac%Reactants = 0
ALLOCATE(SurfChemReac%Products(SurfChemReac%NumOfReact,2))
SurfChemReac%Products = 0
ALLOCATE(SurfChemReac%Inhibition(SurfChemReac%NumOfReact))
SurfChemReac%Inhibition = 0
ALLOCATE(SurfChemReac%Promotion(SurfChemReac%NumOfReact))
SurfChemReac%Promotion = 0
ALLOCATE(SurfChemReac%BoundisChemSurf(nPartBound))
SurfChemReac%BoundisChemSurf = .FALSE.
ALLOCATE(SurfChemReac%NumOfBounds(SurfChemReac%NumOfReact))
SurfChemReac%NumOfBounds = 0
ALLOCATE(SurfChemReac%BoundMap(SurfChemReac%NumOfReact))
ALLOCATE(SurfChemReac%PSMap(nPartBound))
SurfChemReac%CatBoundNum = 0

! Surface map
DO iBound=1, nPartBound
  ALLOCATE(SurfChemReac%PSMap(iBound)%PureSurfReac(ReadInNumOfReact))
  SurfChemReac%PSMap(iBound)%PureSurfReac = .FALSE.
END DO

! Adsorption parameter
ALLOCATE(SurfChemReac%S_initial(SurfChemReac%NumOfReact))
SurfChemReac%S_initial = 0.0
ALLOCATE(SurfChemReac%EqConstant(SurfChemReac%NumOfReact))
SurfChemReac%EqConstant = 1.0
ALLOCATE(SurfChemReac%DissOrder(SurfChemReac%NumOfReact))
SurfChemReac%DissOrder = 0.0
ALLOCATE(SurfChemReac%StickCoeff(SurfChemReac%NumOfReact))
SurfChemReac%StickCoeff = 0.0
ALLOCATE(SurfChemReac%DissociativeAds(SurfChemReac%NumOfReact))
SurfChemReac%DissociativeAds = .FALSE.
ALLOCATE(SurfChemReac%AdsorbedProduct(SurfChemReac%NumOfReact))
SurfChemReac%AdsorbedProduct = 0.0
ALLOCATE(SurfChemReac%GasProduct(SurfChemReac%NumOfReact))
SurfChemReac%GasProduct = 0.0

! Desorption parameter
ALLOCATE(SurfChemReac%E_initial(SurfChemReac%NumOfReact))
SurfChemReac%E_initial = 0.0
ALLOCATE(SurfChemReac%W_interact(SurfChemReac%NumOfReact))
SurfChemReac%W_interact = 0.0
ALLOCATE(SurfChemReac%C_a(SurfChemReac%NumOfReact))
SurfChemReac%C_a = 0.0
ALLOCATE(SurfChemReac%C_b(SurfChemReac%NumOfReact))
SurfChemReac%C_b = 0.0

! General reaction parameter
ALLOCATE(SurfChemReac%Rate(SurfChemReac%NumOfReact))
SurfChemReac%Rate = 0.0
ALLOCATE(SurfChemReac%Prob(SurfChemReac%NumOfReact))
SurfChemReac%Prob = 0.0
ALLOCATE(SurfChemReac%ArrheniusEnergy(SurfChemReac%NumOfReact))
SurfChemReac%ArrheniusEnergy = 0.0
ALLOCATE(SurfChemReac%Prefactor(SurfChemReac%NumOfReact))
SurfChemReac%Prefactor = 0.0
ALLOCATE(SurfChemReac%EReact(SurfChemReac%NumOfReact))
SurfChemReac%EReact = 0.0
ALLOCATE(SurfChemReac%EScale(SurfChemReac%NumOfReact))
SurfChemReac%EScale = 0.0
ALLOCATE(SurfChemReac%HeatAccommodation(SurfChemReac%NumOfReact))
SurfChemReac%HeatAccommodation = 0.0

! Get the reaction names, reactive species and boundaries
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  SurfChemReac%CatName(iReac)              = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-SurfName'))
  SurfChemReac%Reactants(iReac,:)           = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
  SurfChemReac%Products(iReac,:)            = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',2,'0,0')

  SurfChemReac%NumOfBounds(iReac)           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries')
  IF (SurfChemReac%NumOfBounds(iReac).EQ.0) THEN
      CALL abort(__STAMP__,'ERROR: At least one boundary must be defined for each surface reaction!',iReac)
  END IF

  SurfChemReac%BoundMap(iReac)%Boundaries = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Boundaries', &
                                            SurfChemReac%NumOfBounds(iReac))
  ! Define the surface model
  PartBound%SurfaceModel(SurfChemReac%BoundMap(iReac)%Boundaries) = 20

  DO iReac2 = 1, SurfChemReac%NumOfBounds(iReac)
    SurfChemReac%BoundisChemSurf(SurfChemReac%BoundMap(iReac)%Boundaries(iReac2)) = .TRUE.
  END DO

  ! Select pure surface reactions
  DO iVal = 1, SurfChemReac%NumOfBounds(iReac)
    iBound = SurfChemReac%BoundMap(iReac)%Boundaries(iVal)
    SurfChemReac%PSMap(iBound)%PureSurfReac(iReac) = .TRUE.
  END DO
END DO

! Determine the number of boundaries with a surface reaction on them
DO iBound = 1, nPartBound
  IF (SurfChemReac%BoundisChemSurf(iBound)) THEN
    SurfChemReac%CatBoundNum = SurfChemReac%CatBoundNum + 1
  END IF
END DO

! Bulk species involved in the reactions
SurfChemReac%SurfSpecies                   = GETINT('Surface-Species','0')

! Diffusion
SurfChemReac%Diffusion                     = GETLOGICAL('Surface-Diffusion', '.FALSE.')
SurfChemReac%TotDiffusion                  = GETLOGICAL('Surface-TotalDiffusion', '.FALSE.')

! SpeciesDatabase
SurfChemReac%OverwriteCatParameters        = GETLOGICAL('OverwriteCatParameters', '.FALSE.')

IF (SpeciesDatabase.EQ.'none') THEN
  SurfChemReac%OverwriteCatParameters = .TRUE.
END IF

IF(SpeciesDatabase.NE.'none') THEN
  CALL H5OPEN_F(err)

  CALL H5FOPEN_F (TRIM(SpeciesDatabase), H5F_ACC_RDONLY_F, file_id_specdb, err)

  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    dsetname = TRIM('/Surface-Chemistry/'//TRIM(SurfChemReac%CatName(iReac)))

    CALL DatasetExists(file_id_specdb,TRIM(dsetname),DataSetFound)
    IF(.NOT.DataSetFound)THEN
      SurfChemReac%OverwriteCatParameters = .TRUE.
      SWRITE(*,*) 'WARNING: DataSet not found: ['//TRIM(dsetname)//'] ['//TRIM(SpeciesDatabase)//']'
    ELSE
      CALL ReadAttribute(file_id_specdb,'Type',1,DatasetName = dsetname,StrScalar=SurfChemReac%ReactType(iReac))
      CALL AttributeExists(file_id_specdb,'Inhibition',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'Inhibition',1,DatasetName = dsetname,IntScalar=SurfChemReac%Inhibition(iReac))
      ELSE
        SurfChemReac%Inhibition(iReac)= 0
      END IF
      CALL AttributeExists(file_id_specdb,'Promotion',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'Promotion',1,DatasetName = dsetname,IntScalar=SurfChemReac%Promotion(iReac))
      ELSE
        SurfChemReac%Promotion(iReac)= 0
      END IF
      CALL AttributeExists(file_id_specdb,'ReactHeat',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'ReactHeat',1,DatasetName = dsetname,RealScalar=SurfChemReac%EReact(iReac))
      ELSE
        SurfChemReac%EReact(iReac)= 0.
      END IF
      CALL AttributeExists(file_id_specdb,'HeatScaling',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'HeatScaling',1,DatasetName = dsetname,RealScalar=SurfChemReac%EScale(iReac))
      ELSE
        SurfChemReac%EScale(iReac)= 0.
      END IF
      CALL AttributeExists(file_id_specdb,'EnergyAccommodation',TRIM(dsetname), AttrExists=Attr_Exists)
      IF (Attr_Exists) THEN
        CALL ReadAttribute(file_id_specdb,'EnergyAccommodation',1,DatasetName = dsetname,RealScalar=SurfChemReac%HeatAccommodation(iReac))
      ELSE
        SurfChemReac%HeatAccommodation(iReac)= 0.
      END IF

      SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))
      CASE('A')
        SurfChemReac%DissociativeAds(iReac) = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-DissociativeAdsorption', '.FALSE.')
        CALL AttributeExists(file_id_specdb,'StickingCoefficient',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'StickingCoefficient',1,DatasetName = dsetname,RealScalar=SurfChemReac%S_initial(iReac))
        ELSE
          SurfChemReac%S_initial(iReac)= 1.
        END IF
        CALL AttributeExists(file_id_specdb,'EqConstant',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'EqConstant',1,DatasetName = dsetname,RealScalar=SurfChemReac%EqConstant(iReac))
        ELSE
          SurfChemReac%EqConstant(iReac)= 1.
        END IF
        CALL AttributeExists(file_id_specdb,'DissOrder',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'DissOrder',1,DatasetName = dsetname,RealScalar=SurfChemReac%DissOrder(iReac))
        ELSE
          SurfChemReac%DissOrder(iReac)= 1.
        END IF
        ! Special case: dissociative adsorption
        IF (SurfChemReac%DissociativeAds(iReac)) THEN
          CALL AttributeExists(file_id_specdb,'AdsorptionProduct',TRIM(dsetname), AttrExists=Attr_Exists)
          IF (Attr_Exists) THEN
            CALL ReadAttribute(file_id_specdb,'AdsorptionProduct',1,DatasetName = dsetname,IntScalar=SurfChemReac%AdsorbedProduct(iReac))
          ELSE
            CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
          END IF
          CALL AttributeExists(file_id_specdb,'GasPhaseProduct',TRIM(dsetname), AttrExists=Attr_Exists)
          IF (Attr_Exists) THEN
            CALL ReadAttribute(file_id_specdb,'GasPhaseProduct',1,DatasetName = dsetname,IntScalar=SurfChemReac%GasProduct(iReac))
          ELSE
            CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
          END IF
        END IF

      CASE('D')
        CALL AttributeExists(file_id_specdb,'LateralInteraction',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'LateralInteraction',1,DatasetName = dsetname,RealScalar=SurfChemReac%W_interact(iReac))
        ELSE
          SurfChemReac%W_interact(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Ca',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Ca',1,DatasetName = dsetname,RealScalar=SurfChemReac%C_a(iReac))
        ELSE
          SurfChemReac%C_a(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Cb',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Cb',1,DatasetName = dsetname,RealScalar=SurfChemReac%C_b(iReac))
        ELSE
          SurfChemReac%C_b(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac%Prefactor(iReac))
        ELSE
          SurfChemReac%Prefactor(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac%E_initial(iReac))
        ELSE
          SurfChemReac%E_initial(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'DissOrder',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'DissOrder',1,DatasetName = dsetname,RealScalar=SurfChemReac%DissOrder(iReac))
        ELSE
          SurfChemReac%DissOrder(iReac)= 1.
        END IF

      CASE('LH')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac%ArrheniusEnergy(iReac))
        ELSE
          SurfChemReac%ArrheniusEnergy(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac%Prefactor(iReac))
        ELSE
          SurfChemReac%Prefactor(iReac)= 1.
        END IF

      CASE('LHD')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac%ArrheniusEnergy(iReac))
        ELSE
          SurfChemReac%ArrheniusEnergy(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac%Prefactor(iReac))
        ELSE
          SurfChemReac%Prefactor(iReac)= 1.
        END IF

      CASE('ER')
        CALL AttributeExists(file_id_specdb,'Energy',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Energy',1,DatasetName = dsetname,RealScalar=SurfChemReac%ArrheniusEnergy(iReac))
        ELSE
          SurfChemReac%ArrheniusEnergy(iReac)= 0.
        END IF
        CALL AttributeExists(file_id_specdb,'Prefactor',TRIM(dsetname), AttrExists=Attr_Exists)
        IF (Attr_Exists) THEN
          CALL ReadAttribute(file_id_specdb,'Prefactor',1,DatasetName = dsetname,RealScalar=SurfChemReac%Prefactor(iReac))
        ELSE
          SurfChemReac%Prefactor(iReac)= 1.
        END IF

      CASE DEFAULT
        SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac%ReactType(iReac))
        CALL abort(__STAMP__,'Surface Reaction Type does not exist')
      END SELECT

    END IF !DatasetFound
  END DO !iReac
  ! Close the file.
  CALL H5FCLOSE_F(file_id_specdb, err)
  ! Close FORTRAN interface.
  CALL H5CLOSE_F(err)
END IF !SpeciesDatabase


IF (SurfChemReac%OverwriteCatParameters) THEN
  ! Loop over the surface reactions
  DO iReac = 1, ReadInNumOfReact
    WRITE(UNIT=hilf,FMT='(I0)') iReac
    SurfChemReac%ReactType(iReac)             = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))

    SurfChemReac%Inhibition(iReac)            = GETINT('Surface-Reaction'//TRIM(hilf)//'-Inhibition','0')
    SurfChemReac%Promotion(iReac)             = GETINT('Surface-Reaction'//TRIM(hilf)//'-Promotion','0')
    SurfChemReac%EReact(iReac)                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ReactHeat','0.')
    SurfChemReac%EScale(iReac)                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-HeatScaling','0.')
    SurfChemReac%HeatAccommodation(iReac)     = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EnergyAccommodation','1.')

    SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))
    CASE('A')
      SurfChemReac%S_initial(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-StickingCoefficient','1.')
      SurfChemReac%EqConstant(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EqConstant','1.')
      SurfChemReac%DissOrder(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.')
      SurfChemReac%DissociativeAds(iReac) = GETLOGICAL('Surface-Reaction'//TRIM(hilf)//'-DissociativeAdsorption', '.FALSE.')
      ! Special case of the dissociative adsorption, half of the molecule is desorbed back into the gas-phase
      IF (SurfChemReac%DissociativeAds(iReac)) THEN
        SurfChemReac%AdsorbedProduct(iReac) = GETINT('Surface-Reaction'//TRIM(hilf)//'-AdsorptionProduct','0')
        SurfChemReac%GasProduct(iReac)      = GETINT('Surface-Reaction'//TRIM(hilf)//'-GasPhaseProduct','0')
        IF ((SurfChemReac%GasProduct(iReac).EQ.0).OR.(SurfChemReac%GasProduct(iReac).EQ.0)) THEN
          CALL abort(__STAMP__,'Product not defined for the dissociative-adsorption')
        END IF
      END IF

    CASE('D')
      SurfChemReac%W_interact(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-LateralInteraction','0.')
      SurfChemReac%C_a(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Ca','0.')
      SurfChemReac%C_b(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Cb','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','0.')
      SurfChemReac%E_initial(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%DissOrder(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.')

    CASE('LH')
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE('LHD')
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE('ER')
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE DEFAULT
      SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac%ReactType(iReac))
      CALL abort(__STAMP__,'Surface Reaction Type does not exist')
    END SELECT
  END DO

END IF

END SUBROUTINE InitializeVariablesSurfaceChemistry


SUBROUTINE SurfaceModel_Chemistry_Init()
!===================================================================================================================================
! Allocation of side-specific arrays for chemistry modelling
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  ,ONLY: PartBound, nSurfSample, nComputeNodeSurfTotalSides, SurfSide2GlobalSide
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

ALLOCATE( ChemSampWall(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
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

END SUBROUTINE SurfaceModel_Chemistry_Init


!===================================================================================================================================
!> Selection and execution of a catalytic gas-surface interaction
!> 0.) Determine the surface parameters: Coverage and number of surface molecules
!> 1.) Calculate the sticking coefficient by the Kisliuk model (adsorption)
!> 2.) Calculate the reaction probability by the Arrhenius equation (bias-free for multiple channels)
!> 3.) Choose the occuring pathway by comparison with a random number
!> 4.) Perform the chosen process
!>   a.) Adsorption: delete the incoming particle and update the surface values
!>   b.) ER: delete the incoming particle, update the surface values and create the gas phase products
!===================================================================================================================================
SUBROUTINE SurfaceModelChemistry(PartID,SideID,GlobalElemID,n_Loc,PartPosImpact)
! MODULES
! ROUTINES / FUNCTIONS
USE MOD_Globals                   ,ONLY: abort,UNITVECTOR,OrthoNormVec
USE MOD_DSMC_PolyAtomicModel      ,ONLY: DSMC_SetInternalEnr
USE MOD_part_operations           ,ONLY: RemoveParticle, CreateParticle
USE MOD_part_tools                ,ONLY: VeloFromDistribution, GetParticleWeight
USE MOD_SurfaceModel_Tools        ,ONLY: MaxwellScattering
! VARIABLES
USE MOD_Globals_Vars              ,ONLY: PI, BoltzmannConst
USE MOD_TimeDisc_Vars             ,ONLY: dt
USE MOD_Particle_Vars             ,ONLY: PartSpecies,Species,usevMPF,PartMPF, WriteMacroSurfaceValues
USE MOD_Particle_Tracking_Vars    ,ONLY: TrackInfo
USE MOD_Particle_Boundary_Vars    ,ONLY: PartBound, GlobalSide2SurfSide, dXiEQ_SurfSample,SurfSideArea_Shared
USE MOD_SurfaceModel_Vars         ,ONLY: nPorousBC, SurfChemReac , ChemWallProp, ChemSampWall
USE MOD_Particle_Mesh_Vars        ,ONLY: SideInfo_Shared, BoundsOfElem_Shared
USE MOD_Particle_Vars             ,ONLY: PDM, LastPartPos
USE MOD_DSMC_Vars                 ,ONLY: RadialWeighting, DSMC, SamplingActive
USE MOD_Particle_Boundary_Tools   ,ONLY: CalcWallSample
USE MOD_SurfaceModel_Tools        ,ONLY: CalcPostWallCollVelo
USE MOD_Particle_Mesh_Vars        ,ONLY: BoundsOfElem_Shared
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
REAL               :: TempErgy         !< temperature, energy or velocity used for VeloFromDistribution
INTEGER            :: locBCID
INTEGER            :: iBC, SurfSideID
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
INTEGER            :: NewPartID, iNewPart
INTEGER            :: SurfNumOfReac
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
SurfNumOfReac = SurfChemReac%NumOfReact
! MacroParticleFactor
partWeight = GetParticleWeight(PartID)
IF(.NOT.(usevMPF.OR.RadialWeighting%DoRadialWeighting)) THEN
  partWeight = partWeight * Species(speciesID)%MacroParticleFactor
END IF

IF(PartBound%LatticeVec(locBCID).GT.0.) THEN
  ! Number of surface molecules in dependence of the occupancy of the unit cell
  SurfMol = PartBound%MolPerUnitCell(locBCID) * SurfSideArea_Shared(SubP, SubQ,SurfSideID) &
                              /(PartBound%LatticeVec(locBCID)*PartBound%LatticeVec(locBCID))
ELSE
  ! Alternative calculation by the average number of surface molecules per area for a monolayer
  SurfMol = 10.**19 * SurfSideArea_Shared(SubP, SubQ,SurfSideID)
END IF ! LatticeVec.GT.0

DO iReac = 1, SurfNumOfReac
  SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))

  ! 1.) Calculate the sticking coefficient by the Kisliuk model (adsorption)
  CASE('A')
    IF(ANY(SurfChemReac%Reactants(iReac,:).EQ.speciesID)) THEN
      iReac_Ads = iReac

      ! Absolute coverage in terms of the number of surface molecules
      ! Special case: dissociative adsorption
      IF (SurfChemReac%DissociativeAds(iReac)) THEN
        iProd = SurfChemReac%AdsorbedProduct(iReac)
        Coverage = ChemWallProp(iProd,1,SubP,SubQ,SurfSideID)
      ELSE
        IF(ANY(SurfChemReac%Products(iReac,:).NE.0)) THEN
          DO iValProd=1, SIZE(SurfChemReac%Products(iReac,:))
            IF(SurfChemReac%Products(iReac,iValProd).NE.0) THEN
              iProd = SurfChemReac%Products(iReac,iValProd)
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
      DissOrder = SurfChemReac%DissOrder(iReac)
      S_0 = SurfChemReac%S_initial(iReac)
      EqConstant = SurfChemReac%EqConstant(iReac)
      StickCoeff = SurfChemReac%StickCoeff(iReac)

      ! Determine the heat of adsorption in dependence of the coverage [J]
      AdsHeat = (SurfChemReac%EReact(iReac) - Coverage * SurfChemReac%EScale(iReac)) * BoltzmannConst

      ! Theta = free surface sites required for the adsorption
      ! Determination of possible coadsorption processes
      IF(SurfChemReac%Inhibition(iReac).NE.0) THEN
        iCoadsReac = SurfChemReac%Inhibition(iReac)
        iCoadsSpec = SurfChemReac%Reactants(iCoadsReac,1)
        CoAds_Coverage = ChemWallProp(iCoadsSpec,1,SubP, SubQ, SurfSideID)
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        Theta = 1.0 - Coverage/MaxCoverage - CoAds_Coverage/CoAds_MaxCov
      ELSE IF(SurfChemReac%Promotion(iReac).NE.0) THEN
        iCoadsReac = SurfChemReac%Promotion(iReac)
        iCoadsSpec = SurfChemReac%Reactants(iCoadsReac,1)
        CoAds_Coverage = ChemWallProp(iCoadsSpec,1,SubP, SubQ, SurfSideID)
        CoAds_MaxCov = PartBound%MaxCoverage(locBCID,iCoadsSpec)
        Theta = 1.0 - Coverage/MaxCoverage + CoAds_Coverage/CoAds_MaxCov
      ELSE
        Theta = 1.0 - Coverage/MaxCoverage
      END IF

      ! Check whether the maximum coverage value is reached:
      IF(Theta.GE.0.0 .AND. TotalCoverage.LT.MaxTotalCov) THEN
        Theta = Theta**DissOrder
      ! Kisliuk model (for EqConstant=1 and MaxCoverage=1: Langmuir model)
        StickCoeff = S_0 * (1.0 + EqConstant * (1.0/Theta - 1.0))**(-1.0)
      ELSE
        StickCoeff = 0.0
      END IF

    END IF

  ! 2.) Calculate the reaction probability by the Arrhenius equation (bias-free for multiple channels)
  CASE('ER')
    IF(ANY(SurfChemReac%Reactants(iReac,:).EQ.speciesID)) THEN

      ! Definition of the variables
      WallTemp = PartBound%WallTemp(locBCID)
      nu = SurfChemReac%Prefactor(iReac)
      E_act = SurfChemReac%ArrheniusEnergy(iReac)
      Rate = SurfChemReac%Rate(iReac)
      BetaCoeff = SurfChemReac%HeatAccommodation(iReac)

      ! Check for the coverage values of the reactant adsorbed on the surface
      IF(ANY(SurfChemReac%Reactants(iReac,:).NE.speciesID)) THEN
        DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
          IF(SurfChemReac%Reactants(iReac,iValReac).NE.speciesID .AND. SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
            iReactant = SurfChemReac%Reactants(iReac,iValReac)
            IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
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
      AdsDens = Coverage * SurfMol / SurfSideArea_Shared(SubP, SubQ,SurfSideID)

      ! Determine the reaction heat in dependence of the coverage [J]
      ReacHeat = (SurfChemReac%EReact(iReac) - Coverage * SurfChemReac%EScale(iReac)) * BoltzmannConst

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

! 3.) Choose the occuring pathway by comparison with a random number
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

! 4.) Perform the chosen process
SELECT CASE(TRIM(InteractionType))

! 4a.) Adsorption: delete the incoming particle and update the surface values
CASE('A')
  ! Dissociative adsorption
  IF (SurfChemReac%DissociativeAds(iReac)) THEN
    CALL RemoveParticle(PartID)

    ! Heat flux on the surface created by the adsorption
    ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) + AdsHeat * partWeight
    ! Update the number of adsorbed molecules by binding half of the molecule to the surface
    iProd = SurfChemReac%AdsorbedProduct(iReac)
    ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) = ChemSampWall(iProd, 1,SubP,SubQ, SurfSideID) + DissOrder * partWeight

    ! Re-insert the other half of the molecule into the gas-phase
    WallVelo = PartBound%WallVelo(1:3,locBCID)
    CALL OrthoNormVec(n_loc,tang1,tang2)

    ! Get Elem Center
    BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobalElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,2,GlobalElemID)), &
                                SUM(BoundsOfElem_Shared(1:2,3,GlobalElemID)) /) / 2.

    iProd = SurfChemReac%GasProduct(iReac)

    NewVelo(1:3) = CalcPostWallCollVelo(iProd,0.,WallTemp,BetaCoeff)
    NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
    NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*PartPosImpact(1:3)

    CALL CreateParticle(iProd,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=partWeight)

    CALL DSMC_SetInternalEnr(iProd,locBCID,NewPartID,4,iReac)

    IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
    CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)

  ELSE
    CALL RemoveParticle(PartID)

    ! Heat flux on the surface created by the adsorption
    ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) + AdsHeat * partWeight
    ! Update the number of adsorbed molecules
    IF(ANY(SurfChemReac%Products(iReac,:).NE.0)) THEN
      DO iValProd=1, SIZE(SurfChemReac%Products(iReac,:))
        IF(SurfChemReac%Products(iReac,iValProd).NE.0) THEN
          iProd = SurfChemReac%Reactants(iReac,iValProd)
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
  ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) = ChemSampWall(speciesID, 2,SubP,SubQ, SurfSideID) + ReacHeat * partWeight * BetaCoeff

  ! Create the Eley-Rideal reaction product
  ! Incomplete energy accomodation: remaining energy is added to the product
  WallVelo = PartBound%WallVelo(1:3,locBCID)
  CALL OrthoNormVec(n_loc,tang1,tang2)

  ! Get Elem Center
  BoundsOfElemCenter(1:3) = (/SUM(BoundsOfElem_Shared(1:2,1,GlobalElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,2,GlobalElemID)), &
                              SUM(BoundsOfElem_Shared(1:2,3,GlobalElemID)) /) / 2.

  DO iValProd=1, SIZE(SurfChemReac%Products(iReac,:))
    IF(SurfChemReac%Products(iReac,iValProd).NE.0) THEN
      iProd = SurfChemReac%Products(iReac,iValProd)

      NewVelo(1:3) = CalcPostWallCollVelo(iProd,0.,WallTemp,BetaCoeff)
      NewVelo(1:3) = tang1(1:3)*NewVelo(1) + tang2(1:3)*NewVelo(2) - n_loc(1:3)*NewVelo(3) + WallVelo(1:3)
      NewPos(1:3) = eps*BoundsOfElemCenter(1:3) + eps2*PartPosImpact(1:3)

      CALL CreateParticle(iProd,NewPos(1:3),GlobalElemID,NewVelo(1:3),0.,0.,0.,NewPartID=NewPartID, NewMPF=partWeight)

      CALL DSMC_SetInternalEnr(iProd,locBCID,NewPartID,4,iReac)

      ! Sampling of newly created particles
      IF((DSMC%CalcSurfaceVal.AND.SamplingActive).OR.(DSMC%CalcSurfaceVal.AND.WriteMacroSurfaceValues)) &
      CALL CalcWallSample(NewPartID,SurfSideID,'new',SurfaceNormal_opt=n_loc)
    END IF
  END DO

  ! Update the number of adsorbed molecules
  IF(ANY(SurfChemReac%Reactants(iReac,:).NE.speciesID)) THEN
    DO iValReac=1, SIZE(SurfChemReac%Reactants(iReac,:))
      IF(SurfChemReac%Reactants(iReac,iValReac).NE.speciesID .AND. SurfChemReac%Reactants(iReac,iValReac).NE.0) THEN
        iReactant = SurfChemReac%Reactants(iReac,iValReac)
        IF(iReactant.NE.SurfChemReac%SurfSpecies) THEN
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

END MODULE MOD_SurfaceModel_Chemistry
