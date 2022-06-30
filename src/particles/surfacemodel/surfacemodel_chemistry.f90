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
CALL prms%CreateStringOption(   'Surface-Reaction[$]-Type',  &
                                'No default, options are A (adsorption), D (desorption), ER (Eley-Rideal), LH (Langmuir-Hinshelwood), LHD', &
                                 numberedmulti=.TRUE.)
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
CALL prms%CreateRealOption(     'Surface-Reaction[$]-EnergyAccomodation', &
                                    'Energy accomodation coefficient', '0.' , numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Reaction[$]-Inhibition','Inhibition/Coadsorption behaviour due to other reactions', &
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
USE MOD_PARTICLE_Vars           ,ONLY: nSpecies
USE MOD_Mesh_Vars               ,ONLY: SideToElem, offsetElem
USE MOD_Particle_Mesh_Tools     ,ONLY: GetGlobalNonUniqueSideID
USE MOD_Particle_Mesh_Vars      ,ONLY: SideInfo_Shared
USE MOD_Particle_Boundary_Vars  
USE MOD_SurfaceModel_Vars     
USE MOD_Particle_Surfaces_Vars
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
CHARACTER(LEN=3)      :: hilf
INTEGER               :: iReac, iReac2, iSpec, iBound, iVal
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
ALLOCATE(SurfChemReac%Reactants(SurfChemReac%NumOfReact,2))
SurfChemReac%Reactants = 0
ALLOCATE(SurfChemReac%Products(SurfChemReac%NumOfReact,2))
SurfChemReac%Products = 0
ALLOCATE(SurfChemReac%Inhibition(SurfChemReac%NumOfReact))
SurfChemReac%Inhibition = 0
ALLOCATE(SurfChemReac%EForm(SurfChemReac%NumOfReact))
SurfChemReac%EForm = 0.0
ALLOCATE(SurfChemReac%EReact(SurfChemReac%NumOfReact))
SurfChemReac%EReact = 0.0
ALLOCATE(SurfChemReac%EScale(SurfChemReac%NumOfReact))
SurfChemReac%EScale = 0.0
ALLOCATE(SurfChemReac%HeatAccomodation(SurfChemReac%NumOfReact))
SurfChemReac%HeatAccomodation = 0.0
ALLOCATE(SurfChemReac%BoundisChemSurf(nPartBound))
ALLOCATE(SurfChemReac%NumOfBounds(SurfChemReac%NumOfReact))
ALLOCATE(SurfChemReac%BoundMap(SurfChemReac%NumOfReact))
ALLOCATE(SurfChemReac%PSMap(nPartBound))

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
 
DO iReac = 1, ReadInNumOfReact
  WRITE(UNIT=hilf,FMT='(I0)') iReac
  SurfChemReac%ReactType(iReac)             = TRIM(GETSTR('Surface-Reaction'//TRIM(hilf)//'-Type'))

  SurfChemReac%Reactants(iReac,:)           = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Reactants',2,'0,0')
  SurfChemReac%Products(iReac,:)            = GETINTARRAY('Surface-Reaction'//TRIM(hilf)//'-Products',2,'0,0') 
  SurfChemReac%Inhibition(iReac)            = GETINT('Surface-Reaction'//TRIM(hilf)//'-Inhibition','0') 
  SurfChemReac%EReact(iReac)                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-ReactHeat','0.') 
  SurfChemReac%EScale(iReac)                = GETREAL('Surface-Reaction'//TRIM(hilf)//'-HeatScaling','0.') 
  SurfChemReac%HeatAccomodation(iReac)      = GETREAL('Surface-Reaction'//TRIM(hilf)//'-EnergyAccomodation','0.')
  SurfChemReac%NumOfBounds(iReac)           = GETINT('Surface-Reaction'//TRIM(hilf)//'-NumOfBoundaries','0')
  IF (SurfChemReac%NumOfBounds(iReac).EQ.0) THEN
      CALL abort(&
    __STAMP__&
    ,'ERROR: At least one boundary must be defined for each surface reaction!',iReac)
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
  SurfChemReac%EForm(iReac)  = GETREAL('Surface-Reaction'//TRIM(hilf)//'-FormationEnergy','0')
END DO


ALLOCATE( ChemSampWall(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ALLOCATE(ChemDesorpWall(1:nSpecies,2,1:nSurfSample,1:nSurfSample,1:nComputeNodeSurfTotalSides))
ChemDesorpWall = 0.0
ChemSampWall = 0.0
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
