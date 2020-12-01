!==================================================================================================================================
! Copyright (c) 2015-2019 Wladimir Reschke
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

MODULE MOD_SurfaceModel_Init
!===================================================================================================================================
!> Module for initialization of surface models
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
PUBLIC :: DefineParametersSurfModel
PUBLIC :: InitSurfaceModel
PUBLIC :: FinalizeSurfaceModel
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for surface model
!==================================================================================================================================
SUBROUTINE DefineParametersSurfModel()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("SurfaceModel")

#if (PP_TimeDiscMethod==42)
CALL prms%CreateLogicalOption(  'Surface-Adsorption-LateralInactive'&
  , 'Flag for disabling lateral innteractions. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-Adsorption-CoverageReduction'&
  , 'Flag for reducing coverage by time dependant interval. Interval=coverage_start/num_iter. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateLogicalOption(  'Surface-Adsorption-doTPD'&
  , 'Flag for TPD spectrum calculation. Only for TD=42 (RESERVOIR)','.FALSE.')
CALL prms%CreateRealOption(     'Surface-Adsorption-TPD-Beta'&
  , 'Temperature slope used for TPD calculation [K/s]. Surface temperature used for desorption probability is increased'//&
    'each iteration','0.')
#endif
CALL prms%CreateStringOption(  'Particles-SurfCoverageFile'&
  , 'Give relative path to Surface-state-file ([Projectname]_DSMCSurfState***.h5) used for coverage init. '//&
    'File must be of same format as calculation (same mesh, same amount of species, cells and surfacemodel).\n'//&
    'If no file specified:\n'// &
    'Define Part-Species[$]-PartBound[$]-InitialCoverage for initial coverage different than 0.' &
  , 'none')

CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-InitialCoverage'&
  , 'Initial coverage used for species [$] and surfaces of boundary [$] in case of no surface-state-file init' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-ReactionAccomodation'&
  , 'Define energy accomodation coefficient.\n'//&
    'Describes the percentage of reaction enthalpy of surface reaction transmitted to surface.','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBond', 'Bond dissociation energy (K) needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly1'&
                                         , 'Dissociation energy (K) of first bond needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology for Poly atomic species' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-EDissBondPoly2'&
                                         , 'Dissociation energy (K) of second bond needed for calculation of'//&
                                           'adsorption enthalpy in bond order methodology for Poly atomic species' &
                                         , '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Surface-Adsorption-CalcTST'&
                                          , 'Define, which rate coefficient for adsorption are used:\n'//&
                                          ' 0: only those, specified in INI\n'//&
                                          ' 1: those not define in INI are calculated with TST\n'//&
                                          ' 2: calculate every coefficient with TST','2')
CALL prms%CreateIntOption(      'Surface-MaxDissNum', 'Define maximum number of Dissociations per species','0')
CALL prms%CreateIntArrayOption( 'Part-Species[$]-SurfDiss[$]-Products'&
                                          , 'Define product species for surface dissociation number of considered species' &
                                          , '0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-EDissBond'&
                                         , 'Bond dissociation energy (K) for dissociation into product species of'//&
                                          ' surface dissociation','0.', numberedmulti=.TRUE.)

CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Adsorption-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Powerfactor'&
                                         , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-SurfDiss[$]-Prefactor'&
                                         , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Powerfactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Surf-ER[$]-Prefactor'&
                                          , 'TODO-DEFINE-PARAMETER','-1.', numberedmulti=.TRUE.)

CALL prms%CreateIntOption(      'Surface-Nbr-DissocReactions'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                          'Resulting species for given dissoc (2,MaxDissNum,nSpecies)','0')
CALL prms%CreateIntOption(      'Surface-Nbr-ExchangeReactions'&
                                          , 'TODO-DEFINE-PARAMETER','0')
CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateIntArrayOption( 'Surface-ExchReact[$]-Products'&
                                          , 'TODO-DEFINE-PARAMETER','0 , 0', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Reactants'&
                                          , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)
CALL prms%CreateRealArrayOption('Surface-ExchReact[$]-DissBond_K-Products'&
                                         , 'TODO-DEFINE-PARAMETER','0. , 0.', numberedmulti=.TRUE.)

CALL prms%SetSection("Surfacemodel1 [currently not working]")

CALL prms%CreateRealOption(     'Part-Species[$]-MaximumCoverage'&
  , 'Maximum coverage of surfaces with species [$] (used for surfacemodel=1)','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-InitialStick'&
  , 'Initial sticking coefficient (S_0) of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PrefactorStick'&
  , 'Prefactor of sticking coefficient of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '0.', numberedmulti=.TRUE.)
CALL prms%CreateIntOption(      'Part-Species[$]-Adsorbexp'&
  , 'Adsorption exponent of species [$] for surfaces (used for Kisliuk model, surfacemodel=1)' &
  , '1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Nu-a'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                              'Nu exponent a for surface n','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Nu-b'&
                                         , 'TODO-DEFINE-PARAMETER\n'//&
                                              'Nu exponent b for surface n','0.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Desorption-Energy-K'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                            'Desorption energy (K) for surface n','1.', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-Intensification-K'&
                                          , 'TODO-DEFINE-PARAMETER\n'//&
                                             'Intensification energy (K) for surface n','0.', numberedmulti=.TRUE.)


CALL prms%SetSection("Surfacemodel2")

CALL prms%CreateIntOption(     'Part-Species[$]-PartBound[$]-ResultSpec'&
                               ,'Resulting recombination species (one of nSpecies)','-1', numberedmulti=.TRUE.)
CALL prms%CreateRealOption(     'Part-Species[$]-PartBound[$]-RecombinationCoeff'&
                                          , 'TODO-DEFINE-PARAMETER','0.', numberedmulti=.TRUE.)


END SUBROUTINE DefineParametersSurfModel


SUBROUTINE InitSurfaceModel()
!===================================================================================================================================
!> Initialize surface model variables
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars              ,ONLY: nSpecies
USE MOD_ReadInTools                ,ONLY: GETINT
USE MOD_Particle_Boundary_Vars     ,ONLY: nPartBound, PartBound
USE MOD_SurfaceModel_Vars          ,ONLY: Adsorption, SurfModel
!-----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(32)                    :: hilf, hilf2
INTEGER                          :: iSpec, iPartBound
!===================================================================================================================================
IF (.NOT.(ANY(PartBound%Reactive))) RETURN
! allocate info and constants
ALLOCATE( SurfModel%Info(1:nSpecies))
DO iSpec = 1,nSpecies
  SurfModel%Info(iSpec)%WallCollCount = 0
  SurfModel%Info(iSpec)%NumOfAds = 0
  SurfModel%Info(iSpec)%NumOfDes = 0
END DO

! initialize model specific variables
DO iSpec = 1,nSpecies
  WRITE(UNIT=hilf,FMT='(I0)') iSpec
  DO iPartBound=1,nPartBound
    IF(.NOT.PartBound%Reactive(iPartBound)) CYCLE
    WRITE(UNIT=hilf2,FMT='(I0)') iPartBound
    hilf2=TRIM(hilf)//'-PartBound'//TRIM(hilf2)
    SELECT CASE(PartBound%SurfaceModel(iPartBound))
!-----------------------------------------------------------------------------------------------------------------------------------
    CASE(5,6,7)
!-----------------------------------------------------------------------------------------------------------------------------------
      IF (.NOT.ALLOCATED(Adsorption%ResultSpec)) ALLOCATE( Adsorption%ResultSpec(1:nPartBound,1:nSpecies))
      Adsorption%ResultSpec(iPartBound,iSpec) = GETINT('Part-Species'//TRIM(hilf2)//'-ResultSpec')
!-----------------------------------------------------------------------------------------------------------------------------------
    END SELECT
!-----------------------------------------------------------------------------------------------------------------------------------
  END DO
END DO

END SUBROUTINE InitSurfaceModel


SUBROUTINE FinalizeSurfaceModel()
!===================================================================================================================================
!> Deallocate surface model vars
!===================================================================================================================================
! MODULES
USE MOD_SurfaceModel_Vars         ,ONLY: Adsorption, SurfModel
USE MOD_SurfaceModel_Analyze_Vars ,ONLY: SurfModelAnalyzeInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SurfModelAnalyzeInitIsDone=.FALSE.
! variables used if particles are kept after adsorption (currentyl not working)
SDEALLOCATE(SurfModel%Info)
SDEALLOCATE(Adsorption%ResultSpec)

END SUBROUTINE FinalizeSurfaceModel

END MODULE MOD_SurfaceModel_Init
