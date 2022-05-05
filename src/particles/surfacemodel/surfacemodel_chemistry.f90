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
! Initialization of chemical module
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
!> Define parameters for DSMC (Direct Simulation Monte Carlo)
!==================================================================================================================================
SUBROUTINE DefineParametersSurfaceChemistry()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection("Surface Chemistry")
CALL prms%CreateIntOption(      'Surface-NumOfReactions','Number of chemical Surface reactions')
!CALL prms%CreateIntOption(      'Surface-NumOfRecombinations','Number of catalytic Surface recombinations')
CALL prms%CreateStringOption(   'Surface-Reaction[$]-Type',  &
                                'No default', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Reactants'  &
                                           ,'Reactants of Reaction[$]\n'//&
                                            '(SpecNumOfReactant1,\n'//&
                                            'SpecNumOfReactant2)', '0 , 0' , numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Products'  &
                                           ,'Products of Reaction[j] (Product1, Product2)', '0 , 0' &
                                           , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Inhibition','Inhibition/Coadsorption behaviour due to other reactions', &
                                '0', numberedmulti=.TRUE.)
!CALL prms%CreateRealOption(     'Surface-Reaction[$]-ReactProbability', &
!                                    'TODO', '0.' , numberedmulti=.TRUE.)
! CALL prms%CreateRealOption(     'Surface-Reaction[$]-RecombinationCoefficient', &
!                                 'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-StickingCoefficient', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-DissOrder',  &
                                    'TODO', '0.' , numberedmulti=.TRUE.) 
 CALL prms%CreateRealOption(     'Surface-Reaction[$]-EqConstant',  &
                                    'TODO', '0.' , numberedmulti=.TRUE.)                                                                         
CALL prms%CreateRealOption(     'Surface-Reaction[$]-LateralInteraction', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Ca', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Cb', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Prefactor', &
                                    'TODO', '0.' , numberedmulti=.TRUE.) 
CALL prms%CreateRealOption(     'Surface-Reaction[$]-Energy', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)                                
CALL prms%CreateRealOption(     'Surface-Reaction[$]-FormationEnergy', &
                                    'TODO', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-NumOfBoundaries', &
                                          'Num of boundaries for surface reaction.', &
                                            numberedmulti=.TRUE.) 
CALL prms%CreateIntArrayOption( 'Surface-Reaction[$]-Boundaries'  &
                                           ,'Array of boundary indices of surface reaction.' &
                                           ,numberedmulti=.TRUE.)
END SUBROUTINE DefineParametersSurfaceChemistry

SUBROUTINE SurfaceModel_Chemistry_Init()
!===================================================================================================================================
! Readin of variables and definition of reaction cases
!===================================================================================================================================
! MODULES 
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_DSMC_Vars               ,ONLY: SpecDSMC
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies
USE MOD_Particle_Boundary_Vars  ,ONLY: nPartBound, PartBound
USE MOD_SurfaceModel_Vars     
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iReac, iReac2, iSpec, iPart
INTEGER               :: ReadInNumOfReact, ReadInNumOfRecomb
REAL, ALLOCATABLE     :: StoichCoeff(:,:)
!===================================================================================================================================

SurfChemReac%NumOfReact = GETINT('Surface-NumOfReactions', '0')
IF(SurfChemReac%NumOfReact.LE.0) THEN
  RETURN
END IF
ReadInNumOfReact = SurfChemReac%NumOfReact
IF(ReadInNumOfReact.GT.0) THEN
  DoChemSurface = .TRUE.
END IF
SWRITE(*,*) '| Number of considered reaction paths on Surfaces: ', SurfChemReac%NumOfReact

! RecombModel%NumOfReact = GETINT('Surface-NumOfRecombinations', '0')
! IF(RecombModel%NumOfReact.LE.0) THEN
!   RETURN
! END IF
! ReadInNumOfRecomb = RecombModel%NumOfReact
! IF(ReadInNumOfRecomb.GT.0) THEN
!   DoCatSurface = .TRUE.
! END IF
! SWRITE(*,*) '| Number of considered recombiantion paths on Surfaces: ', RecombModel%NumOfReact
!----------------------------------------------------------------------------------------------------------------------------------
ALLOCATE(StoichCoeff(nSpecies,2))
ALLOCATE(SurfChemReac%ReactType(SurfChemReac%NumOfReact))
SurfChemReac%ReactType = '0'
ALLOCATE(SurfChemReac%Reactants(SurfChemReac%NumOfReact,2))
SurfChemReac%Reactants = 0
ALLOCATE(SurfChemReac%Products(SurfChemReac%NumOfReact,2))
SurfChemReac%Products = 0
ALLOCATE(SurfChemReac%Inhibition(SurfChemReac%NumOfReact))
SurfChemReac%Inhibition = 0
!ALLOCATE(SurfChemReac%ReactProb(SurfChemReac%NumOfReact))
!SurfChemReac%ReactProb = 0.0
ALLOCATE(SurfChemReac%EForm(SurfChemReac%NumOfReact))
SurfChemReac%EForm = 0.0
ALLOCATE(SurfChemReac%BoundisChemSurf(nPartBound))
ALLOCATE(SurfChemReac%NumOfBounds(SurfChemReac%NumOfReact))
ALLOCATE(SurfChemReac%BoundMap(SurfChemReac%NumOfReact))

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

! ! Simple recombination model
! ALLOCATE(RecombModel%ReactType(RecombModel%NumOfReact))
! RecombModel%ReactType = '0'
! ALLOCATE(RecombModel%Reactants(RecombModel%NumOfReact,2))
! RecombModel%Reactants = 0
! ALLOCATE(RecombModel%Products(RecombModel%NumOfReact,2))
! RecombModel%Products = 0
! ALLOCATE(RecombModel%RecombCoeff(RecombModel%NumOfReact))
! RecombModel%RecombCoeff = 0.0
! ALLOCATE(RecombModel%BoundisCatSurf(nPartBound))
! ALLOCATE(RecombModel%NumOfBounds(RecombModel%NumOfReact))
! ALLOCATE(RecombModel%BoundMap(RecombModel%NumOfReact))

 
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  SurfChemReac%ReactType(iReac)             = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))

  SurfChemReac%Reactants(iReac,:)           = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
  SurfChemReac%Products(iReac,:)            = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',2,'0,0') 
  SurfChemReac%Inhibition(iReac)            = GETINT('Surface-Reaction'//TRIM(hilf)//'-Inhibition','0') 
  SurfChemReac%NumOfBounds(iReac)           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries','0')
  IF (SurfChemReac%NumOfBounds(iReac).EQ.0) THEN
      CALL abort(&
    __STAMP__&
    ,'ERROR: At least one boundary must be defined for each surface reaction!',iReac)
  END IF
  ALLOCATE(SurfChemReac%BoundMap(iReac)%nSurfacefluxBCs(SurfChemReac%NumOfBounds(iReac)))
  SurfChemReac%BoundMap(iReac)%Boundaries = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Boundaries', &
                                            SurfChemReac%NumOfBounds(iReac))
  ! Define the surface model
  PartBound%SurfaceModel(SurfChemReac%BoundMap(iReac)%Boundaries) = 20

  DO iReac2 = 1, SurfChemReac%NumOfBounds(iReac)   
    SurfChemReac%BoundisChemSurf(SurfChemReac%BoundMap(iReac)%Boundaries(iReac2)) = .TRUE.                                   
  END DO

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
    SurfChemReac%E_initial(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','1.')
    SurfChemReac%DissOrder(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-DissOrder','1.') 

  CASE('LH') 
    SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
    SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

  CASE('ER')
    SurfChemReac%ArrheniusEnergy(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Energy','0.')
    SurfChemReac%Prefactor(iReac) = GETREAL('Surface-Reaction'//TRIM(hilf)//'-Prefactor','1.')

  CASE DEFAULT
    SWRITE(*,*) ' Reaction Type does not exists: ', TRIM(SurfChemReac%ReactType(iReac))
    CALL abort(__STAMP__,'Surface Reaction Type does not exist')
  END SELECT

  ! Total Collision Energy: Arrhenius-based chemistry model
  ! SurfChemReac%ReactProb(iReac)   = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ReactProbability','0')
  SurfChemReac%EForm(iReac)  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-FormationEnergy','0')
  ! Filling up ChemReac-Array for the given non-reactive dissociation/electron-impact ionization partners

  ! Calculation of stoichiometric coefficients and calculation of the heat of formation
  StoichCoeff = 0
  !SurfChemReac%EForm(iReac) = 0.0
  DO iSpec=1, nSpecies
    ! Reactants
    DO iPart = 1,2
      IF(SurfChemReac%Reactants(iReac,iPart).EQ.iSpec) THEN
          StoichCoeff(iSpec,1) = StoichCoeff(iSpec,1) + 1
      END IF
    END DO
    ! Products
    DO iPart = 1,2
      IF(SurfChemReac%Products(iReac,iPart).EQ.iSpec) THEN
        StoichCoeff(iSpec,2) = StoichCoeff(iSpec,2) + 1
      END IF
    END DO
    ! Calculation of the enthalpy of reaction by the summation of the enthalpies of formation of the respective species
    ! (ionization energy of ionized species was already added in dsmc_init.f90)
    SurfChemReac%EForm(iReac) = SurfChemReac%EForm(iReac) &
                            - StoichCoeff(iSpec,2)*SpecDSMC(iSpec)%HeatOfFormation  &
                            + StoichCoeff(iSpec,1)*SpecDSMC(iSpec)%HeatOfFormation
  END DO
  
END DO

! ! Simple recombination model
! DO iReac = 1, ReadInNumOfRecomb
!   WRITE(UNIT=hilf,FMT='(I0)') iReac
!   RecombModel%ReactType(iReac)             = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))

!   RecombModel%Reactants(iReac,:)           = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
!   RecombModel%Products(iReac,:)            = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',2,'0,0') 
!   RecombModel%NumOfBounds(iReac)           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries','0')
!   IF (RecombModel%NumOfBounds(iReac).EQ.0) THEN
!       CALL abort(&
!     __STAMP__&
!     ,'ERROR: At least one boundary must be defined for each surface reaction!',iReac)
!   END IF
!   ALLOCATE(RecombModel%BoundMap(iReac)%Boundaries(RecombModel%NumOfBounds(iReac)))
!   RecombModel%BoundMap(iReac)%Boundaries = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Boundaries', &
!                                             RecombModel%NumOfBounds(iReac))
!   ! Define the surface model
!   PartBound%SurfaceModel(RecombModel%BoundMap(iReac)%Boundaries) = 30

!   DO iReac2 = 1, RecombModel%NumOfBounds(iReac)   
!     RecombModel%BoundisCatSurf(RecombModel%BoundMap(iReac)%Boundaries(iReac2)) = .TRUE.                                   
!   END DO

!   RecombModel%RecombCoeff(iReac)  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-RecombinationCoefficient','0')
  
! END DO

END SUBROUTINE SurfaceModel_Chemistry_Init

END MODULE MOD_SurfaceModel_Chemistry
