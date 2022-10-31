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
PUBLIC :: DefineParametersSurfaceChemistry, SurfaceModel_Chemistry_Init
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
CALL prms%CreateStringOption(   'Surface-Database', 'none')
CALL prms%CreateStringOption(   'Surface-Reaction[$]-SurfName', 'none' ,numberedmulti=.TRUE.)
CALL prms%CreateLogicalOption(  'OverwriteCatParameters', 'Flag to set catalytic parameters manually', '.FALSE.')
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2)', '0 , 0' , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Products'  &
                                           ,'Products of Reaction[j] (Product1, Product2)', '0 , 0' &
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
                                    'Equilibrium constant between the adsorption and desorption (K), Langmuir: K=1', '0.' , numberedmulti=.TRUE.)
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

SUBROUTINE SurfaceModel_Chemistry_Init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES 
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  
USE MOD_SurfaceModel_Vars     
USE MOD_Particle_Surfaces_Vars
USE MOD_io_hdf5
USE MOD_HDF5_input,         ONLY:ReadAttribute, DatasetExists, AttributeExists
#if USE_MPI
USE MOD_MPI_Shared_Vars         ,ONLY: MPI_COMM_SHARED, myComputeNodeRank
USE MOD_MPI_Shared
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=64)     :: dsetname
CHARACTER(LEN=64)     :: AttribName
INTEGER(HID_T)        :: file_id_specdb                       ! File identifier
INTEGER(HID_T)        :: dset_id_specdb                       ! Dataset identifier
LOGICAL               :: DataSetFound
LOGICAL               :: Attr_Exists
INTEGER(HID_T)        :: Loc_ID, Attr_ID
CHARACTER(LEN=32)      :: hilf, hilf2
INTEGER               :: iReac, iReac2, iSpec, iBound, iVal, err
INTEGER               :: ReadInNumOfReact
INTEGER               :: iSide, SideID, iBC
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
SurfChemReac%EqConstant = 0.0
ALLOCATE(SurfChemReac%DissOrder(SurfChemReac%NumOfReact))
SurfChemReac%DissOrder = 0.0
ALLOCATE(SurfChemReac%StickCoeff(SurfChemReac%NumOfReact))
SurfChemReac%StickCoeff = 0.0

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

  SurfChemReac%NumOfBounds(iReac)           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries','0')
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

!SpeciesDatabase
SpeciesDatabase                            = GETSTR('Surface-Database', 'none')
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

      ! Convert the prefactor from absolute to coverage values for associative desorption
        IF(SurfChemReac%DissOrder(iReac).EQ.2) THEN
          SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)
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
        SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)

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
        SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)

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
    SurfChemReac%HeatAccommodation(iReac)      = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EnergyAccommodation','1.')

    SELECT CASE (TRIM(SurfChemReac%ReactType(iReac)))
    CASE('A')
      SurfChemReac%S_initial(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-StickingCoefficient','1.')
      SurfChemReac%EqConstant(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EqConstant','1.')
      SurfChemReac%DissOrder(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.')

    CASE('D')
      SurfChemReac%W_interact(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-LateralInteraction','0.')
      SurfChemReac%C_a(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Ca','0.')
      SurfChemReac%C_b(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Cb','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','0.')
      SurfChemReac%E_initial(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%DissOrder(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.') 

    ! Convert the prefactor from absolute to coverage values for associative desorption
      IF(SurfChemReac%DissOrder(iReac).EQ.2) THEN
        SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)
      END IF

    CASE('LH') 
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')
      ! Convert the prefactor to coverage dependent values
      SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)

    CASE('LHD') 
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')
      ! Convert the prefactor to coverage dependent values
      SurfChemReac%Prefactor(iReac) = SurfChemReac%Prefactor(iReac) * 10.0**(15)

    CASE('ER')
      SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
      SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

    CASE DEFAULT
      SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac%ReactType(iReac))
      CALL abort(__STAMP__,'Surface Reaction Type does not exist')
    END SELECT
  END DO

END IF

ALLOCATE( ChemSampWall(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemSampWall = 0.0
ALLOCATE(ChemDesorpWall(1:nSpecies,1,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemDesorpWall = 0.0
! ALLOCATE(ChemCountReacWall(1:ReadInNumOfReact,1,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
! ChemCountReacWall = 0.0

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

END MODULE MOD_SurfaceModel_Chemistry
