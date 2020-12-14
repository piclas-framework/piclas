!==================================================================================================================================
! Copyright (c) 2010 - 2019 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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

MODULE MOD_part_emission_tools
!===================================================================================================================================
! module for particle emission
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

! no interface allowed (do not remove this comment)
!INTERFACE IntegerDivide
  !MODULE PROCEDURE IntegerDivide
!END INTERFACE

INTERFACE SetParticleChargeAndMass
  MODULE PROCEDURE SetParticleChargeAndMass
END INTERFACE

INTERFACE SetParticleMPF
  MODULE PROCEDURE SetParticleMPF
END INTERFACE

INTERFACE CalcVelocity_maxwell_lpn
  MODULE PROCEDURE CalcVelocity_maxwell_lpn
END INTERFACE

INTERFACE CalcVelocity_taylorgreenvortex
  MODULE PROCEDURE CalcVelocity_taylorgreenvortex
END INTERFACE

INTERFACE SamplePoissonDistri
  MODULE PROCEDURE SamplePoissonDistri
END INTERFACE

INTERFACE BessK
  MODULE PROCEDURE BessK
END INTERFACE

INTERFACE DEVI
  MODULE PROCEDURE DEVI
END INTERFACE

INTERFACE SYNGE
  MODULE PROCEDURE SYNGE
END INTERFACE

INTERFACE QUASIREL
  MODULE PROCEDURE QUASIREL
END INTERFACE

INTERFACE SetCellLocalParticlePosition
  MODULE PROCEDURE SetCellLocalParticlePosition
END INTERFACE

INTERFACE InsideExcludeRegionCheck
  MODULE PROCEDURE InsideExcludeRegionCheck
END INTERFACE

INTERFACE CalcNbrOfPhotons
  MODULE PROCEDURE CalcNbrOfPhotons
END INTERFACE

INTERFACE CalcIntensity_Gaussian
  MODULE PROCEDURE CalcIntensity_Gaussian
END INTERFACE

INTERFACE CalcVelocity_FromWorkFuncSEE
  MODULE PROCEDURE CalcVelocity_FromWorkFuncSEE
END INTERFACE

INTERFACE DSMC_SetInternalEnr_LauxVFD
  MODULE PROCEDURE DSMC_SetInternalEnr_LauxVFD
END INTERFACE

#if CODE_ANALYZE
INTERFACE CalcVectorAdditionCoeffs
  MODULE PROCEDURE CalcVectorAdditionCoeffs
END INTERFACE
#endif /*CODE_ANALYZE*/

!===================================================================================================================================
PUBLIC :: CalcVelocity_taylorgreenvortex, CalcVelocity_gyrotroncircle
PUBLIC :: IntegerDivide,SetParticleChargeAndMass,SetParticleMPF,CalcVelocity_maxwell_lpn,SamplePoissonDistri
PUBLIC :: BessK,DEVI,SYNGE,QUASIREL
PUBLIC :: SetCellLocalParticlePosition,InsideExcludeRegionCheck
PUBLIC :: SetParticlePositionPoint, SetParticlePositionEquidistLine, SetParticlePositionLine, SetParticlePositionDisk
PUBLIC :: SetParticlePositionCircle, SetParticlePositionGyrotronCircle, SetParticlePositionCuboidCylinder
PUBLIC :: SetParticlePositionSphere, SetParticlePositionSinDeviation, SetParticleTimeStep
PUBLIC :: CalcNbrOfPhotons, CalcPhotonEnergy
PUBLIC :: CalcIntensity_Gaussian
PUBLIC :: CalcVelocity_FromWorkFuncSEE, DSMC_SetInternalEnr_LauxVFD
PUBLIC :: SetParticlePositionPhotonSEEDisc, SetParticlePositionPhotonCylinder
#if CODE_ANALYZE
PUBLIC :: CalcVectorAdditionCoeffs
#endif /*CODE_ANALYZE*/
!===================================================================================================================================
CONTAINS


SUBROUTINE IntegerDivide(Ntot,length,Ai,Ni)
!===================================================================================================================================
! Divide the Integer Ntot into separate Ni inside different "areas" Ai (attention: old Ni is counted up -> needs to be initialized!)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: Ntot, length
REAL,INTENT(IN)                  :: Ai(1:length)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)            :: Ni(1:length)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iN, iRan, Nitemp, Nrest, Ntot0
REAL            :: Atot, Bi(0:length), RandVal1, A2i(1:length), A2tot !,Error,Nrel(1:length),Arel(1:length)
!===================================================================================================================================

IF(Ntot.EQ.0) RETURN

Atot=0.
Ntot0=0
DO iN=1,length
  Atot=Atot+Ai(iN)
  Ntot0=Ntot0+Ni(iN)
END DO
!print*,Ai/Atot

!-- divide into INT-parts
Nrest=Ntot
A2tot=0.
Bi(:)=0.
DO iN=1,length
  Nitemp=INT(REAL(Ai(iN))/REAL(Atot)*Ntot) !INT-part
  Ni(iN)=Ni(iN)+Nitemp
  Nrest=Nrest-Nitemp !remaining number
  A2i(iN)=REAL(Ai(iN))/REAL(Atot)*Ntot - Nitemp !elem weight for remaining number
  A2tot=A2tot+A2i(iN)
  Bi(iN)=A2tot !elem upper limit for remaining number
END DO

!-- distribute remaining number
IF (Nrest.LT.0) THEN
  CALL abort(&
__STAMP__&
,'ERROR 1 in IntegerDivide!')
ELSE IF (Nrest.GT.0) THEN
  DO iN=1,length
    Bi(iN)=Bi(iN)/A2tot !normalized upper limit
  END DO
  DO iRan=1,Nrest
    CALL RANDOM_NUMBER(RandVal1)
    DO iN=1,length
      IF( Bi(iN-1).LT.RandVal1 .AND. RandVal1.LE.Bi(iN) ) THEN
        Ni(iN)=Ni(iN)+1
        EXIT
      END IF
    END DO
  END DO
END IF

!-- test if remaining number was distributed
Nrest=Ntot+Ntot0
DO iN=1,length
  Nrest=Nrest-Ni(iN)
END DO
IF (Nrest.NE.0) THEN
  IPWRITE(*,*) 'Ntot: ',Ntot
  IPWRITE(*,*) 'Ntot0: ',Ntot0
  IPWRITE(*,*) 'Nrest: ',Nrest
  CALL abort(&
__STAMP__&
,'ERROR 2 in IntegerDivide!')
END IF

!Error=0
!DO iN=1,length
!  Nrel(iN)=REAL(Ni(iN))/REAL(Ntot)
!  Arel(iN)=Ai(iN)      /Atot
!  Error=Error+(Nrel(iN)-Arel(iN))**2
!END DO
!IPWRITE(*,*)'Error=',Error

END SUBROUTINE IntegerDivide


SUBROUTINE SetParticleTimeStep(NbrOfParticle)
!===================================================================================================================================
!> Set the particle time step based on its position (2D/axisymmetric simulations). Can only be used for particle not yet added to
!> the particle vector, loops over the total number of particles and the indices in the nextFreePosition array.
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars           ,ONLY: PDM, VarTimeStep, PEM, PartState
USE MOD_Particle_VarTimeStep    ,ONLY: CalcVarTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iPart, PositionNbr
!===================================================================================================================================
DO iPart=1, NbrOfParticle
  PositionNbr = PDM%nextFreePosition(iPart+PDM%CurrentNextFreePosition)
  VarTimeStep%ParticleTimeStep(PositionNbr) = &
                CalcVarTimeStep(PartState(1,PositionNbr), PartState(2,PositionNbr),PEM%LocalElemID(PositionNbr))
END DO

END SUBROUTINE SetParticleTimeStep


SUBROUTINE SetParticleChargeAndMass(FractNbr,NbrOfParticle)
!===================================================================================================================================
!> Set the mass and charge of the particles based on the input species index.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,    ONLY : PDM, PartSpecies
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr
!===================================================================================================================================
DO i=1, NbrOfParticle
  PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
  IF (PositionNbr .NE. 0) THEN
    PartSpecies(PositionNbr) = FractNbr
  ELSE
    CALL abort(&
    __STAMP__&
    ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
  END IF
END DO

END SUBROUTINE SetParticleChargeAndMass


SUBROUTINE SetParticleMPF(FractNbr,NbrOfParticle)
!===================================================================================================================================
! finally, set the MPF
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars           ,ONLY : PDM, PartMPF, Species, PartState
USE MOD_DSMC_Vars               ,ONLY : RadialWeighting
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                       :: FractNbr
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                    :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                  :: i,PositionNbr
!===================================================================================================================================
i = 1
DO WHILE (i .le. NbrOfParticle)
  PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
  IF (PositionNbr .NE. 0) THEN
    IF(RadialWeighting%DoRadialWeighting) THEN
      PartMPF(PositionNbr) = CalcRadWeightMPF(PartState(2,PositionNbr),FractNbr,PositionNbr)
    ELSE
      PartMPF(PositionNbr) = Species(FractNbr)%MacroParticleFactor
    END IF
  ELSE
    CALL abort(&
    __STAMP__&
    ,'ERROR in SetParticlePosition:ParticleIndexNbr.EQ.0 - maximum nbr of particles reached?')
  END IF
  i = i + 1
END DO

END SUBROUTINE SetParticleMPF


SUBROUTINE CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit, Temperature)
!===================================================================================================================================
!> Calculate the velocity vector for a particle of a certain species based on the Maxwell distribution function.
!> Input: Index of the init or a temperature. Routine is suitable for low particle numbers (lpn).
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY: BoltzmannConst
USE MOD_Particle_Vars,          ONLY: Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: FractNbr
INTEGER,INTENT(IN), OPTIONAL    :: iInit
REAL,INTENT(IN), OPTIONAL       :: Temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: RandVal(3), Velo1, Velo2, Velosq, Tx, ty, Tz, v_drift(3)
!===================================================================================================================================
IF(PRESENT(iInit).AND.PRESENT(Temperature)) CALL abort(&
  __STAMP__&
  ,'CalcVelocity_maxwell_lpn: iInit and Temperature cannot both be input arguments!')
IF(PRESENT(iInit))THEN
  Tx=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Ty=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Tz=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  v_drift=Species(FractNbr)%Init(iInit)%VeloIC * Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
ELSE IF(PRESENT(Temperature))THEN
  Tx=Temperature
  Ty=Temperature
  Tz=Temperature
  v_drift=0.0
ELSE
  CALL abort(&
        __STAMP__&
        ,'CalcVelocity_maxwell_lpn: No iInit or temperature given for species: ', FractNbr)
END IF

Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1 = 2.*RandVal(1) - 1.
  Velo2 = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
! x component
Vec3D(1) = Velo1*SQRT(-2*BoltzmannConst*Tx / Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)
! y component
Vec3D(2) = Velo2*SQRT(-2*BoltzmannConst*Ty / Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)
Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1 = 2.*RandVal(1) - 1.
  Velo2 = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
! z component
Vec3D(3) = Velo1*SQRT(-2*BoltzmannConst*Tz / Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)
! Adding the bulk drift velocity
Vec3D(1:3) = Vec3D(1:3) + v_drift

END SUBROUTINE CalcVelocity_maxwell_lpn


SUBROUTINE DSMC_SetInternalEnr_LauxVFD(iSpecies, iInit, iPart, init_or_sf)
!===================================================================================================================================
!> Energy distribution according to dissertation of Laux (diatomic)
!===================================================================================================================================
! MODULES
USE MOD_Globals                 ,ONLY: abort
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_DSMC_Vars               ,ONLY: PartStateIntEn, SpecDSMC, DSMC
USE MOD_Particle_Vars           ,ONLY: Species, PEM, Adaptive_MacroVal
USE MOD_DSMC_ElectronicModel    ,ONLY: InitElectronShell
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)             :: iSpecies, iInit, iPart, init_or_sf
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: iRan
INTEGER                         :: iQuant
REAL                            :: TVib                       ! vibrational temperature
REAL                            :: TRot                       ! rotational temperature
INTEGER                         :: ElemID
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! Set internal energies (vibrational and rotational)
!-----------------------------------------------------------------------------------------------------------------------------------
ElemID = PEM%LocalElemID(iPart)
IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
  SELECT CASE (init_or_sf)
  CASE(1) !iInit
    TVib=SpecDSMC(iSpecies)%Init(iInit)%TVib
    TRot=SpecDSMC(iSpecies)%Init(iInit)%TRot
  CASE(2) !SurfaceFlux
    IF(Species(iSpecies)%Surfaceflux(iInit)%Adaptive) THEN
      SELECT CASE(Species(iSpecies)%Surfaceflux(iInit)%AdaptiveType)
        CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
          TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
          TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
        CASE(2) ! adaptive Outlet/freestream
          TVib = Species(iSpecies)%Surfaceflux(iInit)%AdaptivePressure &
                  / (BoltzmannConst * Adaptive_MacroVal(DSMC_NUMDENS,ElemID,iSpecies))
          TRot = TVib
        CASE DEFAULT
          CALL abort(&
          __STAMP__&
          ,'Wrong adaptive type for Surfaceflux in int_energy -> lauxVDF!')
      END SELECT
    ELSE
      TVib=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TVib
      TRot=SpecDSMC(iSpecies)%Surfaceflux(iInit)%TRot
    END IF
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'neither iInit nor Surfaceflux defined as reference!')
  END SELECT
  ! Set vibrational energy
  CALL RANDOM_NUMBER(iRan)
  iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
  DO WHILE (iQuant.GE.SpecDSMC(iSpecies)%MaxVibQuant)
    CALL RANDOM_NUMBER(iRan)
    iQuant = INT(-LOG(iRan)*TVib/SpecDSMC(iSpecies)%CharaTVib)
  END DO
  PartStateIntEn( 1,iPart) = (iQuant + DSMC%GammaQuant)*SpecDSMC(iSpecies)%CharaTVib*BoltzmannConst
  ! Set rotational energy
  CALL RANDOM_NUMBER(iRan)
  PartStateIntEn( 2,iPart) = -BoltzmannConst*TRot*LOG(iRan)
ELSE
  ! Nullify energy for atomic species
  PartStateIntEn( 1,iPart) = 0
  PartStateIntEn( 2,iPart) = 0
END IF
!-----------------------------------------------------------------------------------------------------------------------------------
! Set electronic energy
!-----------------------------------------------------------------------------------------------------------------------------------
IF (DSMC%ElectronicModel) THEN
  IF((SpecDSMC(iSpecies)%InterID.NE.4).AND.(.NOT.SpecDSMC(iSpecies)%FullyIonized)) THEN
    CALL InitElectronShell(iSpecies,iPart,iInit,init_or_sf)
  ELSE
    PartStateIntEn( 3,iPart) = 0.
  END IF
ENDIF

END SUBROUTINE DSMC_SetInternalEnr_LauxVFD

SUBROUTINE CalcVelocity_taylorgreenvortex(FractNbr, Vec3D, iInit, Element)
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: BoltzmannConst
USE MOD_Particle_Vars      ,ONLY: Species
USE MOD_Mesh_Vars          ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN)               :: iInit
INTEGER                          :: Element
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RandVal(3), Velo1, Velo2, Velosq, v_drift(3)
REAL                             :: T  ! temperature
REAL                             :: p  ! pressure
REAL                             :: p0 ! base pressure
!===================================================================================================================================

! V0 = Ma*c_s
!   Ma := 0.3
! c_s = sqrt(gamma*R*T/M)
!   gamma := 1.4
!   R     := 8.3144598
!   T     := 273.15
!   M     := 28.0134e-3
ASSOCIATE( V0   => 101.0694686816                       ,& !Species(FractNbr)%Init(iInit)%VeloIC ,&
           x    => ElemBaryNGeo(1,Element)              ,&
           y    => ElemBaryNGeo(2,Element)              ,&
           z    => ElemBaryNGeo(3,Element)              ,&
           L    => GEO%xmaxglob                         ,&
           rho0 => 1.25                                 ,&
           R_N2 => 296.8                                 & ! unit of R_N2 is [J*kg^-1K^-1]
           )

  v_drift(1) =  V0*SIN(x/L)*COS(y/L)*COS(z/L)
  v_drift(2) = -V0*COS(x/L)*SIN(y/L)*COS(z/L)
  v_drift(3) = 0.

  p0 = rho0 * R_N2 * Species(FractNbr)%Init(iInit)%MWTemperatureIC
  p  = p0 + (rho0*V0**2/16.)*( COS(2*x/L)+COS(2*y/L) )*( COS(2*z/L)+2 )
  T  = p / (rho0*R_N2)

END ASSOCIATE

Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1  = 2.*RandVal(1) - 1.
  Velo2  = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
Vec3D(1) = Velo1*SQRT(-2*BoltzmannConst*T/Species(FractNbr)%MassIC*LOG(Velosq)/Velosq) !x-Komponente
Vec3D(2) = Velo2*SQRT(-2*BoltzmannConst*T/Species(FractNbr)%MassIC*LOG(Velosq)/Velosq) !y-Komponente
Velosq = 2
DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
  CALL RANDOM_NUMBER(RandVal)
  Velo1  = 2.*RandVal(1) - 1.
  Velo2  = 2.*RandVal(2) - 1.
  Velosq = Velo1**2 + Velo2**2
END DO
Vec3D(3) = Velo1*SQRT(-2*BoltzmannConst*T/Species(FractNbr)%MassIC*LOG(Velosq)/Velosq) !z-Komponente


Vec3D(1:3) = Vec3D(1:3) + v_drift

END SUBROUTINE CalcVelocity_taylorgreenvortex

SUBROUTINE CalcVelocity_gyrotroncircle(FractNbr, Vec3D, iInit, iPart)
!===================================================================================================================================
! Subroutine to sample current cell values (partly copied from 'LD_DSMC_Mean_Bufferzone_A_Val' and 'dsmc_analyze')
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars      ,ONLY: Species, PartState
USE MOD_PICInterpolation_vars ,ONLY: externalField
INTEGER,INTENT(IN)               :: FractNbr, iPart
INTEGER,INTENT(IN)               :: iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                 :: r1, r2, x_1, y_1, x_2, y_2, a, b, e, f, g, x_01, x_02, y_01, y_02, RandVal1, Radius(3), n_vec(3), tan_vec(3)
REAL                 :: NormalIC(1:3), RadiusIC, RadiusICGyro, Alpha, GyroVecDirSIGN, VeloIC
!===================================================================================================================================
!! Position of particle on gyro circle changed in SetParticlePosition.F90: Problem
!! We don't have the radius-vector any more. Thus transport the radius vector from there to here.
! Or do Alternative way: Hack the radius by intersecting two circles (big IC and small gyro circle)
  IF (externalField(6).NE.0) THEN
    GyroVecDirSIGN = -externalField(6)/(ABS(externalField(6)))
  ELSE
    GyroVecDirSIGN = -1.
  END IF
  NormalIC=Species(FractNbr)%Init(iInit)%NormalIC(1:3)
  RadiusIC=Species(FractNbr)%Init(iInit)%RadiusIC
  RadiusICGyro=Species(FractNbr)%Init(iInit)%RadiusICGyro
  Alpha=Species(FractNbr)%Init(iInit)%alpha
  VeloIC=Species(FractNbr)%Init(iInit)%VeloIC
  r1 = RadiusIC
  r2 = RadiusICGyro
  x_1 = 0.
  y_1 = 0.
  x_2 = PartState(1,iPart)
  y_2 = PartState(2,iPart)
  IF (x_1 .eq. x_2) THEN
    a = (x_1 - x_2)/(y_2-y_1)
    b = ((r1**2-r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
    & /(2.*y_2-2.*y_1)
    e = (a**2+1.)
    f = (2.*a*(b-y_1))-2.*x_1
    g = (b-y_1)**2-r1**2+x_1**2
    ! intersection points
    x_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
    x_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
    y_01 = x_01 * a + b
    y_02 = x_02 * a + b
  ELSE
    a = (y_1 - y_2)/(x_2-x_1)
    b = ((r1**2 - r2**2)-(x_1**2-x_2**2)-(y_1**2-y_2**2))&
    & /(2.*x_2 - 2. * x_1)
    e = (a**2 + 1.)
    f = 2. * a * (b - x_1) -2 *y_1
    g = (b-x_1)**2 - r1**2 + y_1**2
    y_01 = (-f + SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
    y_02 = (-f - SQRT(ABS(f**2 - 4. * e * g)))/(2.*e) ! the term in SQRT can be -0.0 , therefore the ABS
    x_01 = y_01 * a + b
    x_02 = y_02 * a + b
  END IF
  CALL RANDOM_NUMBER(RandVal1)
  IF (RandVal1 .ge. 0.5) THEN
    Radius(1) = PartState(1,iPart) - x_01
    Radius(2) = PartState(2,iPart) - y_01
  ELSE
    Radius(1) = PartState(1,iPart) - x_02
    Radius(2) = PartState(2,iPart) - y_02
  END IF

  Radius(3) = 0.
  !Check if Radius has correct length
  IF ((SQRT(Radius(1)**2+Radius(2)**2)-r1).ge.1E-15) THEN
    IPWRITE(UNIT_stdOut,*)"Error in setparticle velocity, gyrotron circle. &
    & Radius too big after intersection."
  END IF
  !  Normal Vector of circle
  n_vec(1:3) = NormalIC(1:3)
  Radius(1:3) = Radius(1:3) / SQRT(Radius(1)**2+Radius(2)**2+Radius(3)**2)
  !  Vector Product rxn
  tan_vec(1) = Radius(2)*n_vec(3) * GyroVecDirSIGN - Radius(3)*n_vec(2)
  tan_vec(2) = Radius(3)*n_vec(1) - Radius(1)*n_vec(3) * GyroVecDirSIGN
  tan_vec(3) = Radius(1)*n_vec(2) - Radius(2)*n_vec(1)
  ! If Gyrotron resonator: Add velocity in normal direction!
  IF (Alpha .gt. 0.) THEN
    n_vec = n_vec * ( 1. / Alpha )
  ELSE
    n_vec = 0.
  END IF

  Vec3D(1:3) = (tan_vec(1:3) + n_vec(1:3)) * VeloIC

  IF (ABS(SQRT(Vec3D(1)*Vec3D(1) + Vec3D(2)*Vec3D(2))- VeloIC) .GT. 10.) THEN
    SWRITE(*,'(A,3(E21.14,X))') 'Velocity=', PartState(4:6,iPart)
    CALL abort(&
    __STAMP__&
    ,'ERROR in gyrotron_circle spaceIC!',iPart)
  END If
  IF (Vec3D(1).NE.Vec3D(1).OR.Vec3D(2).NE.Vec3D(2).OR.Vec3D(3).NE.Vec3D(3)) THEN
    SWRITE(*,'(A,3(E21.14,X))') 'WARNING:! NaN: Velocity=', Vec3D(1:3)
  END If

END SUBROUTINE CalcVelocity_gyrotroncircle

FUNCTION BessK(ord,arg)
!===================================================================================================================================
! Modified Bessel function of second kind and integer order (currently only 2nd...) and real argument,
! required for Maxwell-Juettner distribution
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,    ONLY: PI,EuMas
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,   INTENT(IN)  :: arg
INTEGER,INTENT(IN)  :: ord
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL                :: BessK
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: BessI0, BessI1, BessK0, BessK1, BessK0_old
REAL     :: rr, eps, ct, w0
REAL     :: set_a(12), set_b(12), set_c(8)
INTEGER  :: kk, k0
!===================================================================================================================================

  !em = 0.577215664901533_8        ! Euler–Mascheroni constant
  eps= 1E-15_8

  set_a = (/0.125E0_8, 7.03125E-2_8,                  &
          7.32421875E-2_8, 1.1215209960938E-1_8,      &
          2.2710800170898E-1_8, 5.7250142097473E-1_8, &
          1.7277275025845E0_8, 6.0740420012735E0_8,    &
          2.4380529699556E01_8, 1.1001714026925E02_8, &
          5.5133589612202E02_8, 3.0380905109224E03_8/)

  set_b = (/-0.375E0_8, -1.171875E-1_8,                 &
          -1.025390625E-1_8, -1.4419555664063E-1_8,     &
          -2.7757644653320E-1_8, -6.7659258842468E-1_8, &
          -1.9935317337513E0_8, -6.8839142681099E0_8,   &
          -2.7248827311269E01_8, -1.2159789187654E02_8, &
          -6.0384407670507E02_8, -3.3022722944809E03_8/)

  set_c = (/0.125E0_8, 0.2109375E0_8,                 &
          1.0986328125E0_8, 1.1775970458984E01_8,     &
          2.1461706161499E2_8, 5.9511522710323E03_8,  &
          2.3347645606175E05_8, 1.2312234987631E07_8/)


!==========================================================================================!
! Compute I_0(x) and I_1(x)
!==========================================================================================!
  IF (arg .EQ. 0.) THEN
    BessI1 = 0.
    BessI0 = 1.

  ELSE IF (arg .LE. 18.) THEN
    BessI0 = 1.
    rr     = 1.
    kk     = 0
    DO WHILE ((rr/BessI0) .GT. eps)
      kk = kk+1
      rr = .25*rr*arg*arg/(kk*kk)
      BessI0 = BessI0 + rr
    END DO
!     WRITE(*,*) 'BessI0:', BessI0
!     WRITE(*,*) kk
    BessI1 = 1.
    rr     = 1.
    kk     = 0
    DO WHILE ((rr/BessI1) .GT. eps)
      kk = kk+1
      rr = .25*rr*arg*arg/(kk*(kk+1))
      BessI1 = BessI1 + rr
    END DO
    BessI1 = 0.5*arg*BessI1
!     WRITE(*,*) 'BessI1:', BessI1

  ELSE
    IF      (arg .LT. 35.) THEN
      k0 = 12
    ELSE IF (arg .LT. 50.) THEN
      k0 =  9
    ELSE
      k0 =  7
    END IF
    BessI0 = 1._8
    DO kk = 1,k0
      BessI0 = BessI0 + set_a(kk)*arg**(-kk)
    END DO
    BessI0 = exp(arg)/sqrt(2._8*pi*arg)*BessI0
!     WRITE(*,*) 'BessI0: ', BessI0
    BessI1 = 1._8
    DO kk = 1,k0
      BessI1 = BessI1 + set_b(kk)*arg**(-kk)
    END DO
    BessI1 = exp(arg)/sqrt(2._8*pi*arg)*BessI1
!     WRITE(*,*) 'BessI1: ', BessI1
  END IF

!==========================================================================================!
! Compute K_0(x)
!==========================================================================================!
  IF (arg .LE. 0.) THEN
    CALL abort(&
__STAMP__&
,' mod. Bessel function of second kind requries pos arg:')
  ELSE IF (arg .LE. 9.) THEN
    kk = 1
    ct = -log(arg/2.)-EuMas
    w0 = 1._8
    rr = 0.25*arg*arg
    BessK0 = rr*(w0+ct)
    BessK0_old = 1.E20
    DO WHILE (abs((BessK0-BessK0_old)/BessK0) .GT. eps)
      kk = kk+1
      BessK0_old = BessK0
      w0 = w0+1._8/kk
      rr = 0.25*rr*arg*arg/(kk*kk)
      BessK0 = BessK0 + rr*(w0+ct)
    END DO
    BessK0 = BessK0 + ct
  ELSE
    BessK0 = 1._8
    DO kk = 1,8
      BessK0 = BessK0 + set_c(kk)*arg**(-2._8*kk)
    END DO
    BessK0 = BessK0/(2._8*arg*BessI0)
!     WRITE(*,*) 'BessK0: ', BessK0
  END IF

!==========================================================================================!
! Compute K_1(x) and K_n(x)
!==========================================================================================!
  BessK1 = (1._8/arg-BessI1*BessK0)/BessI0
  BessK = 2._8*(ord-1._8)*BessK1/arg + BessK0

END FUNCTION BessK


PURE FUNCTION DEVI(mass, temp, gamma)
!===================================================================================================================================
! derivative to find max of function
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,  ONLY: BoltzmannConst,c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: mass, temp, gamma
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL                :: DEVI
!===================================================================================================================================
  DEVI = mass*c2/(BoltzmannConst*temp)* &
           gamma*(gamma*gamma-1._8)-5._8*gamma*gamma+3._8
END FUNCTION DEVI


PURE FUNCTION SYNGE(velabs, temp, mass, BK2)
!===================================================================================================================================
! Maxwell-Juettner distribution according to Synge Book p.48
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,   ONLY: BoltzmannConst,c_inv,c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)   :: velabs, temp, mass, BK2
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL              :: SYNGE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
REAL              :: gamma
!===================================================================================================================================
gamma = 1./sqrt(1.-(velabs*c_inv)*(velabs*c_inv))
SYNGE = velabs*velabs*gamma**5/BK2*exp(-mass*c2*gamma/(BoltzmannConst*temp))
END FUNCTION SYNGE


PURE FUNCTION QUASIREL(velabs, temp, mass)
!===================================================================================================================================
! discard gamma in the prefactor, maintain it in the computation of the energy
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars,   ONLY: BoltzmannConst,c_inv,c2
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL ,INTENT(IN)    :: velabs, temp, mass
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLE
REAL     :: QUASIREL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE
REAL     :: gamma
!===================================================================================================================================
  gamma = 1/sqrt(1-(velabs*c_inv)*(velabs*c_inv))
  QUASIREL = velabs*velabs*gamma**5._8* &
               exp((1._8-gamma)*mass*c2/(BoltzmannConst*temp))
END FUNCTION QUASIREL


SUBROUTINE SamplePoissonDistri(RealTarget,IntSample,Flag_opt)
!===================================================================================================================================
! Sample IntSample from Poisson-Distri around RealTarget (if Flag present it will be turned off at sample limit, otherwise abort)
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                :: RealTarget
LOGICAL,INTENT(INOUT),OPTIONAL :: Flag_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: IntSample
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL         :: Flag
INTEGER         :: Npois
REAL            :: Tpois, RandVal1
!===================================================================================================================================

IF (PRESENT(Flag_opt)) THEN
  Flag=Flag_opt
ELSE
  Flag=.FALSE.
END IF

Npois=0
Tpois=1.0
CALL RANDOM_NUMBER(RandVal1)
DO
  Tpois=RandVal1*Tpois
  IF (Tpois.LT.TINY(Tpois)) THEN
    IF (Flag) THEN !Turn off Poisson Sampling and "sample" by random-rounding
      IPWRITE(*,*)'WARNING: target is too large for poisson sampling: switching now to Random rounding...'
      IntSample = INT(RealTarget + RandVal1)
      Flag = .FALSE.
      EXIT
    ELSE !Turning off not allowed: abort (RealTarget must be decreased ot PoissonSampling turned off manually)
      CALL abort(&
__STAMP__&
,'ERROR in SamplePoissonDistri: RealTarget (e.g. flux) is too large for poisson sampling!')
    END IF
  END IF
  IF (Tpois.GT.EXP(-RealTarget)) THEN
    Npois=Npois+1
    CALL RANDOM_NUMBER(RandVal1)
  ELSE
    IntSample = Npois
    EXIT
  END IF
END DO

END SUBROUTINE SamplePoissonDistri


SUBROUTINE SetCellLocalParticlePosition(chunkSize,iSpec,iInit,UseExactPartNum)
!===================================================================================================================================
!> routine for inserting particles positions locally in every cell
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars              ,ONLY: RadialWeighting
USE MOD_part_tools              ,ONLY: CalcRadWeightMPF
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Localization  ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemEpsOneCell
USE MOD_Particle_Mesh_Vars     ,ONLY: BoundsOfElem_Shared,ElemVolume_Shared,ElemMidPoint_Shared
USE MOD_Mesh_Tools             ,ONLY: GetCNElemID
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Tracking_Vars ,ONLY: TrackingMethod
USE MOD_Particle_Vars          ,ONLY: Species, PDM, PartState, PEM, Symmetry, VarTimeStep
USE MOD_Particle_VarTimeStep   ,ONLY: CalcVarTimeStep
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)              :: iSpec
INTEGER, INTENT(IN)              :: iInit
LOGICAL, INTENT(IN)              :: UseExactPartNum
!-----------------------------------------------------------------------------------------------------------------------------------
! INOUTPUT VARIABLES
INTEGER, INTENT(INOUT)           :: chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: iElem, ichunkSize
INTEGER                          :: iPart,  nPart
REAL                             :: iRan, RandomPos(3)
REAL                             :: PartDens
LOGICAL                          :: InsideFlag
REAL                             :: Det(6,2)
REAL                             :: RefPos(1:3)
INTEGER                          :: CellChunkSize(1+offsetElem:nElems+offsetElem)
INTEGER                          :: chunkSize_tmp, ParticleIndexNbr
REAL                             :: adaptTimestep
!-----------------------------------------------------------------------------------------------------------------------------------
  IF (UseExactPartNum) THEN
    IF(chunkSize.GE.PDM%maxParticleNumber) THEN
      CALL abort(&
__STAMP__,&
'ERROR in SetCellLocalParticlePosition: Maximum particle number reached! max. particles needed: ',chunksize)
    END IF
    CellChunkSize(:)=0
    CALL IntegerDivide(chunkSize,nElems,ElemVolume_Shared(1+offsetElem:nElems+offsetElem) &
        ,CellChunkSize(1+offsetElem:nElems+offsetElem))
  ELSE
    PartDens = Species(iSpec)%Init(iInit)%PartDensity / Species(iSpec)%MacroParticleFactor   ! numerical Partdensity is needed
    IF(RadialWeighting%DoRadialWeighting) PartDens = PartDens * 2. / (RadialWeighting%PartScaleFactor)
    chunkSize_tmp = INT(PartDens * LocalVolume)
    IF(chunkSize_tmp.GE.PDM%maxParticleNumber) THEN
      CALL abort(&
__STAMP__,&
'ERROR in SetCellLocalParticlePosition: Maximum particle number during sanity check! max. particles needed: ',chunkSize_tmp)
    END IF
  END IF

  ichunkSize = 1
  ParticleIndexNbr = 1
  DO iElem = 1+offsetElem, nElems+offsetElem
    ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
      IF (UseExactPartNum) THEN
        nPart = CellChunkSize(iElem)
      ELSE
        IF(RadialWeighting%DoRadialWeighting) THEN
          PartDens = Species(iSpec)%Init(iInit)%PartDensity / CalcRadWeightMPF(ElemMidPoint_Shared(2,GetCNElemID(iElem)), iSpec)
        END IF
        CALL RANDOM_NUMBER(iRan)
        IF(VarTimeStep%UseVariableTimeStep) THEN
          adaptTimestep = CalcVarTimeStep(ElemMidPoint_Shared(1,GetCNElemID(iElem)), ElemMidPoint_Shared(2,GetCNElemID(iElem)), iElem)
          nPart = INT(PartDens / adaptTimestep * ElemVolume_Shared(GetCNElemID(iElem)) + iRan)
        ELSE
          nPart = INT(PartDens * ElemVolume_Shared(GetCNElemID(iElem)) + iRan)
        END IF
      END IF
      DO iPart = 1, nPart
        ParticleIndexNbr = PDM%nextFreePosition(iChunksize + PDM%CurrentNextFreePosition)
        IF (ParticleIndexNbr .ne. 0) THEN
          InsideFlag=.FALSE.
          DO WHILE(.NOT.InsideFlag)
            CALL RANDOM_NUMBER(RandomPos)
            IF(Symmetry%Axisymmetric.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
              ! Treatment of axisymmetry without weighting
              RandomPos(1) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
              RandomPos(2) = SQRT(RandomPos(2)*(Bounds(2,2)**2-Bounds(1,2)**2)+Bounds(1,2)**2)
            ELSE
              RandomPos = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
            END IF
            IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
            IF(Symmetry%Order.LE.1) RandomPos(2) = 0.

            SELECT CASE(TrackingMethod)
              CASE(TRIATRACKING)
                CALL ParticleInsideQuad3D(RandomPos,iElem,InsideFlag,Det)

              CASE(TRACING)
                CALL PartInElemCheck(RandomPos,iPart,iElem,InsideFlag)

              CASE(REFMAPPING)
                CALL GetPositionInRefElem(RandomPos,RefPos,iElem)
                IF (MAXVAL(ABS(RefPos)).LE.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
            END SELECT
          END DO
          PartState(1:3,ParticleIndexNbr) = RandomPos(1:3)
          PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
          PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
          PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
          PEM%GlobalElemID(ParticleIndexNbr) = iElem
          ichunkSize = ichunkSize + 1
        ELSE
          WRITE(UNIT_stdOut,*) ""
          IPWRITE(UNIT_stdOut,*) "ERROR:"
          IPWRITE(UNIT_stdOut,*) "                iPart :", iPart
          IPWRITE(UNIT_stdOut,*) "PDM%maxParticleNumber :", PDM%maxParticleNumber
          CALL abort(&
              __STAMP__&
              ,'ERROR in SetCellLocalParticlePosition: Maximum particle number reached during inserting! --> ParticleIndexNbr.EQ.0')
        END IF
      END DO
    END ASSOCIATE
  END DO
  chunkSize = ichunkSize - 1

END SUBROUTINE SetCellLocalParticlePosition


SUBROUTINE SetParticlePositionPoint(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3)
INTEGER                 :: i
!===================================================================================================================================
Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC
DO i=1,chunkSize
   particle_positions(i*3-2) = Particle_pos(1)
   particle_positions(i*3-1) = Particle_pos(2)
   particle_positions(i*3  ) = Particle_pos(3)
END DO
END SUBROUTINE SetParticlePositionPoint


SUBROUTINE SetParticlePositionEquidistLine(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), VectorGap(3)
INTEGER                 :: i
!===================================================================================================================================
  IF(chunkSize.EQ.1)THEN
     Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + 0.5 * Species(FractNbr)%Init(iInit)%BaseVector1IC
  ELSE
    VectorGap = Species(FractNbr)%Init(iInit)%BaseVector1IC/(REAL(chunkSize)-1.)
    DO i=1,chunkSize
      Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (i-1)*VectorGap
      particle_positions(i*3-2) = Particle_pos(1)
      particle_positions(i*3-1) = Particle_pos(2)
      particle_positions(i*3  ) = Particle_pos(3)
    END DO
  END IF
END SUBROUTINE SetParticlePositionEquidistLine


SUBROUTINE SetParticlePositionLine(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), iRan
INTEGER                 :: i
!===================================================================================================================================
  DO i=1,chunkSize
    CALL RANDOM_NUMBER(iRan)
    Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC*iRan
    particle_positions(i*3-2) = Particle_pos(1)
    particle_positions(i*3-1) = Particle_pos(2)
    particle_positions(i*3  ) = Particle_pos(3)
  END DO
END SUBROUTINE SetParticlePositionLine


SUBROUTINE SetParticlePositionDisk(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_SuperB_Tools           ,ONLY: FindLinIndependentVectors, GramSchmidtAlgo
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), RandVec(2), lineVector(3), lineVector2(3), radius
INTEGER                 :: i
!===================================================================================================================================
  CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  DO i=1,chunkSize
   radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
   DO WHILE(radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC)
      CALL RANDOM_NUMBER(RandVec)
      RandVec = RandVec * 2. - 1.
      Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
               (RandVec(1) * lineVector + RandVec(2) *lineVector2)

      radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
                     (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
                     (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
                     (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
                     (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
                     (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
   END DO
   particle_positions(i*3-2) = Particle_pos(1)
   particle_positions(i*3-1) = Particle_pos(2)
   particle_positions(i*3  ) = Particle_pos(3)
  END DO
END SUBROUTINE SetParticlePositionDisk


SUBROUTINE SetParticlePositionCircle(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_SuperB_Tools           ,ONLY: FindLinIndependentVectors, GramSchmidtAlgo
USE MOD_Globals_Vars           ,ONLY: Pi
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), iRan, lineVector(3), lineVector2(3), radius, Phi
INTEGER                 :: i
!===================================================================================================================================
  CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  radius = Species(FractNbr)%Init(iInit)%RadiusIC
  DO i=1,chunkSize
    IF(TRIM(Species(FractNbr)%Init(iInit)%SpaceIC).EQ.'circle') THEN
      CALL RANDOM_NUMBER(iRan)
      Phi = 2.*Pi*iRan
    ELSE
      Phi = 2.*Pi*REAL(i)/ REAL(chunkSize)
    END IF
    Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC +        &
                  linevector * COS(Phi) * radius +  &
                  linevector2 * SIN(Phi) * radius
    particle_positions(i*3-2) = Particle_pos(1)
    particle_positions(i*3-1) = Particle_pos(2)
    particle_positions(i*3  ) = Particle_pos(3)
  END DO
END SUBROUTINE SetParticlePositionCircle


SUBROUTINE SetParticlePositionGyrotronCircle(FractNbr,iInit,chunkSize,particle_positions, NbrOfParticle)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_SuperB_Tools           ,ONLY: FindLinIndependentVectors, GramSchmidtAlgo
USE MOD_Globals_Vars           ,ONLY: Pi
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac, dt
USE MOD_PICInterpolation_vars  ,ONLY: useVariableExternalField, VariableExternalField
USE MOD_PICInterpolation       ,ONLY: InterpolateVariableExternalField
USE MOD_Globals_Vars           ,ONLY: c_inv
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize, NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), iRan, lineVector(3), lineVector2(3), radius, Phi, n(3), radius_vec(3)
REAL                    :: JJ(3,3), II(3,3), NN(3,3), rgyrate, Bintpol
INTEGER                 :: i, j
!===================================================================================================================================
  CALL FindLinIndependentVectors(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  CALL GramSchmidtAlgo(Species(FractNbr)%Init(iInit)%NormalIC(1:3), lineVector(1:3), lineVector2(1:3))
  radius = Species(FractNbr)%Init(iInit)%RadiusIC
  DO i=1,chunkSize
     CALL RANDOM_NUMBER(iRan)
     Phi = 2.*Pi*iRan
     Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + (linevector * COS(Phi) + linevector2 * SIN(Phi)) * radius
     ! Change position of particle on the small gyro circle
     ! take normal vecotr of the circle
     n(1:3) = Species(FractNbr)%Init(iInit)%NormalIC(1:3)
     ! generate radius vector (later it will be multiplied by the length of the
     ! gyro circles. For now we just need the vector)
     radius_vec(1) = Particle_pos(1) - Species(FractNbr)%Init(iInit)%BasePointIC(1)
     radius_vec(2) = Particle_pos(2) - Species(FractNbr)%Init(iInit)%BasePointIC(2)
     radius_vec(3) = Particle_pos(3) - Species(FractNbr)%Init(iInit)%BasePointIC(3)
     !rotate radius vector with random angle
     CALL RANDOM_NUMBER(iRan)
     Phi=2.*Pi*iRan
     JJ(1,1:3) = (/   0.,-n(3), n(2)/)
     JJ(2,1:3) = (/ n(3),   0.,-n(1)/)
     JJ(3,1:3) = (/-n(2), n(1),   0./)
     II(1,1:3) = (/1.,0.,0./)
     II(2,1:3) = (/0.,1.,0./)
     II(3,1:3) = (/0.,0.,1./)
     FORALL(j=1:3) NN(:,j) = n(:)*n(j)

     ! 1. determine the z-position in order to get the interpolated curved B-field
     CALL RANDOM_NUMBER(iRan)
     IF (NbrOfParticle.EQ.Species(FractNbr)%Init(iInit)%initialParticleNumber) THEN
       particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                       + iRan * Species(FractNbr)%Init(iInit)%CuboidHeightIC
     ELSE
       particle_positions(i*3  ) = Species(FractNbr)%Init(iInit)%BasePointIC(3) &
                                       + iRan * dt*RKdtFrac &
                                       * Species(FractNbr)%Init(iInit)%VeloIC/Species(FractNbr)%Init(iInit)%alpha
     END IF

     ! 2. calculate curved B-field at z-position in order to determine size of gyro radius
     IF (useVariableExternalField) THEN
        IF(particle_positions(i*3).LT.VariableExternalField(1,1))THEN ! assume particles travel in positive z-direction
          CALL abort(__STAMP__,'SetParticlePosition: particle_positions(i*3) cannot be smaller than VariableExternalField(1,1). Fix *.csv data or emission!')
        END IF
        Bintpol = InterpolateVariableExternalField(particle_positions(i*3))
        rgyrate = 1./ SQRT ( 1. - (Species(FractNbr)%Init(iInit)%VeloIC**2 * (1. + 1./Species(FractNbr)%Init(iInit)%alpha**2)) &
                            * c_inv * c_inv ) * Species(FractNbr)%MassIC * Species(FractNbr)%Init(iInit)%VeloIC / &
                  ( Bintpol * abs( Species(FractNbr)%ChargeIC) )
     ELSE
       rgyrate =  Species(FractNbr)%Init(iInit)%RadiusICGyro
     END IF

     radius_vec = MATMUL( NN+cos(Phi)*(II-NN)+sin(Phi)*JJ , radius_vec )
     radius_vec(1:3) = radius_vec(1:3) / SQRT(radius_vec(1)**2+radius_vec(2)**2+radius_vec(3)**2) &
                   * rgyrate !Species(1)%RadiusICGyro
     ! Set new particles position:
     particle_positions(i*3-2) = Particle_pos(1) + radius_vec(1)
     particle_positions(i*3-1) = Particle_pos(2) + radius_vec(2)
     !particle_positions(i*3  )=0.
  END DO
END SUBROUTINE SetParticlePositionGyrotronCircle


SUBROUTINE SetParticlePositionCuboidCylinder(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species, Symmetry
USE MOD_Timedisc_Vars          ,ONLY: RKdtFrac, dt
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)  :: chunkSize
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), RandVal(3), lineVector(3), radius
INTEGER                 :: i, chunkSize2
LOGICAL                 :: insideExcludeRegion
!===================================================================================================================================
  lineVector(1) = Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2)
  lineVector(2) = Species(FractNbr)%Init(iInit)%BaseVector1IC(3) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(3)
  lineVector(3) = Species(FractNbr)%Init(iInit)%BaseVector1IC(1) * Species(FractNbr)%Init(iInit)%BaseVector2IC(2) - &
    Species(FractNbr)%Init(iInit)%BaseVector1IC(2) * Species(FractNbr)%Init(iInit)%BaseVector2IC(1)
  IF ((lineVector(1).eq.0).AND.(lineVector(2).eq.0).AND.(lineVector(3).eq.0)) THEN
    CALL abort(__STAMP__,'BaseVectors are parallel!')
  ELSE
    lineVector = lineVector / SQRT(lineVector(1) * lineVector(1) + lineVector(2) * lineVector(2) + &
      lineVector(3) * lineVector(3))
  END IF
  i=1
  chunkSize2=0
  DO WHILE (i .LE. chunkSize)
    SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%SpaceIC))
    CASE ('cuboid')
      CALL RANDOM_NUMBER(RandVal)
      Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%BaseVector1IC * RandVal(1)
      Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BaseVector2IC * RandVal(2)
      IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
      ELSE
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CuboidHeightIC * RandVal(3)
      END IF
      IF(Symmetry%Order.EQ.1) Particle_pos(2:3) = 0.
    CASE ('cylinder')
      radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1.
      DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC) .OR.(radius.LT.Species(FractNbr)%Init(iInit)%Radius2IC))
         CALL RANDOM_NUMBER(RandVal)
         Particle_pos = Species(FractNbr)%Init(iInit)%BaseVector1IC * (RandVal(1)*2.-1.) &
                      + Species(FractNbr)%Init(iInit)%BaseVector2IC * (RandVal(2)*2.-1.)
         radius = SQRT( Particle_pos(1) * Particle_pos(1) + &
                        Particle_pos(2) * Particle_pos(2) + &
                        Particle_pos(3) * Particle_pos(3) )
      END DO
      Particle_pos = Particle_pos + Species(FractNbr)%Init(iInit)%BasePointIC
      IF (Species(FractNbr)%Init(iInit)%CalcHeightFromDt) THEN !directly calculated by timestep
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%VeloIC * dt*RKdtFrac * RandVal(3)
      ELSE
        Particle_pos = Particle_pos + lineVector * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal(3)
      END IF
    END SELECT
    IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
      IF (insideExcludeRegion) THEN
        i=i+1
        CYCLE !particle is in excluded region
      END IF
    END IF
    particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
    particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
    particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
    i=i+1
    chunkSize2=chunkSize2+1
  END DO
  chunkSize = chunkSize2
END SUBROUTINE SetParticlePositionCuboidCylinder


SUBROUTINE SetParticlePositionSphere(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Part_tools             ,ONLY: DICEUNITVECTOR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(INOUT)  :: chunkSize
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3), iRan, radius
INTEGER                 :: i, chunkSize2
LOGICAL                 :: insideExcludeRegion
!===================================================================================================================================
  i=1
  chunkSize2=0
  DO WHILE (i .LE. chunkSize)
    CALL RANDOM_NUMBER(iRan)
    radius = Species(FractNbr)%Init(iInit)%RadiusIC*iRan**(1./3.)
    Particle_pos = DICEUNITVECTOR()*radius + Species(FractNbr)%Init(iInit)%BasePointIC
    IF (Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions.GT.0) THEN
      CALL InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
      IF (insideExcludeRegion) THEN
        i=i+1
        CYCLE !particle is in excluded region
      END IF
    END IF
    particle_positions((chunkSize2+1)*3-2) = Particle_pos(1)
    particle_positions((chunkSize2+1)*3-1) = Particle_pos(2)
    particle_positions((chunkSize2+1)*3  ) = Particle_pos(3)
    i=i+1
    chunkSize2=chunkSize2+1
  END DO
  chunkSize = chunkSize2
END SUBROUTINE SetParticlePositionSphere


SUBROUTINE SetParticlePositionSinDeviation(FractNbr,iInit,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO
USE MOD_Globals_Vars           ,ONLY: Pi
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: xlen, ylen, zlen, pilen, x_step, y_step, z_step, x_pos, y_pos
INTEGER                 :: i, iPart, j, k
!===================================================================================================================================
  IF(Species(FractNbr)%Init(iInit)%initialParticleNumber.NE. &
      (Species(FractNbr)%Init(iInit)%maxParticleNumberX * Species(FractNbr)%Init(iInit)%maxParticleNumberY &
      * Species(FractNbr)%Init(iInit)%maxParticleNumberZ)) THEN
   SWRITE(*,*) 'for species ',FractNbr,' does not match number of particles in each direction!'
   CALL abort(__STAMP__,'ERROR: Number of particles in init / emission region',iInit)
  END IF
  xlen = ABS(GEO%xmaxglob  - GEO%xminglob)
  ylen = ABS(GEO%ymaxglob  - GEO%yminglob)
  zlen = ABS(GEO%zmaxglob  - GEO%zminglob)
  pilen=2.0*PI/xlen
  x_step = xlen/Species(FractNbr)%Init(iInit)%maxParticleNumberX
  y_step = ylen/Species(FractNbr)%Init(iInit)%maxParticleNumberY
  z_step = zlen/Species(FractNbr)%Init(iInit)%maxParticleNumberZ
  iPart = 1
  DO i=1,Species(FractNbr)%Init(iInit)%maxParticleNumberX
    x_pos = (i * x_step - x_step*0.5)
    x_pos = GEO%xminglob + x_pos + Species(FractNbr)%Init(iInit)%Amplitude &
            * SIN(Species(FractNbr)%Init(iInit)%WaveNumber * pilen * x_pos)
    DO j=1,Species(FractNbr)%Init(iInit)%maxParticleNumberY
      y_pos =  GEO%yminglob + j * y_step - y_step * 0.5
      DO k=1,Species(FractNbr)%Init(iInit)%maxParticleNumberZ
        particle_positions(iPart*3-2) = x_pos
        particle_positions(iPart*3-1) = y_pos
        particle_positions(iPart*3  ) = GEO%zminglob &
                                  + k * z_step - z_step * 0.5
        iPart = iPart + 1
      END DO
    END DO
  END DO
END SUBROUTINE SetParticlePositionSinDeviation


SUBROUTINE InsideExcludeRegionCheck(FractNbr, iInit, Particle_pos, insideExcludeRegion)
!===================================================================================================================================
! Subroutine for checking if calculated particle position would be inside user-defined ExcludeRegion (cuboid or cylinder)
!===================================================================================================================================
! MODULES
USE MOD_Globals,                ONLY : abort
USE MOD_Particle_Vars,          ONLY : Species
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
REAL,INTENT(IN)                  :: Particle_pos(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)              :: insideExcludeRegion
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: VecExclude(3), DistExclude
INTEGER                          :: iExclude
!===================================================================================================================================

insideExcludeRegion=.FALSE.
DO iExclude=1,Species(FractNbr)%Init(iInit)%NumberOfExcludeRegions
  VecExclude = Particle_pos - Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BasePointIC
  SELECT CASE (TRIM(Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%SpaceIC))
  CASE ('cuboid')
    !--check normal direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CuboidHeightIC) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check BV1 direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector1IC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(1)**2) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check BV2 direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%BaseVector2IC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%ExcludeBV_lenghts(2)**2) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
      RETURN !particle is inside current ExcludeRegion based an all dimensions
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
  CASE ('cylinder')
    !--check normal direction
    DistExclude = VecExclude(1)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(1) &
      + VecExclude(2)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(2) &
      + VecExclude(3)*Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%NormalIC(3)
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%CylinderHeightIC) &
      .AND. (DistExclude .GE. 0.) ) THEN
      insideExcludeRegion = .TRUE.
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
    !--check radial direction
    DistExclude = SQRT( VecExclude(1)**2 + VecExclude(2)**2 + VecExclude(3)**2 - DistExclude**2 )
    IF ( (DistExclude .LE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%RadiusIC) &
      .AND. (DistExclude .GE. Species(FractNbr)%Init(iInit)%ExcludeRegion(iExclude)%Radius2IC) ) THEN
      insideExcludeRegion = .TRUE.
      RETURN !particle is inside current ExcludeRegion based an all dimensions
    ELSE
      insideExcludeRegion = .FALSE.
      CYCLE
    END IF
  CASE DEFAULT
    CALL abort(&
__STAMP__&
,'wrong SpaceIC for ExcludeRegion!')
  END SELECT
END DO

END SUBROUTINE InsideExcludeRegionCheck


#if CODE_ANALYZE
PURE FUNCTION CalcVectorAdditionCoeffs(point,Vector1,Vector2)
!===================================================================================================================================
! robust calculation of Coeffs C(1) and C(2) from point = C(1)*Vector1 + C(2)*Vector2
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)         :: point(3), Vector1(3), Vector2(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL                     :: CalcVectorAdditionCoeffs(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL                     :: denom(3)
!===================================================================================================================================
denom(1)=Vector2(2)*Vector1(3)-Vector2(3)*Vector1(2)
denom(2)=Vector2(1)*Vector1(2)-Vector2(2)*Vector1(1)
denom(3)=Vector2(3)*Vector1(1)-Vector2(1)*Vector1(3)

IF (ABS(denom(1)).GT.ABS(denom(2)) .AND. ABS(denom(1)).GT.ABS(denom(3))) THEN
  CalcVectorAdditionCoeffs(2)=(point(2)*Vector1(3)-point(3)*Vector1(2))/denom(1)
ELSE IF (ABS(denom(2)).GT.ABS(denom(1)) .AND. ABS(denom(2)).GT.ABS(denom(3))) THEN
  CalcVectorAdditionCoeffs(2)=(point(1)*Vector1(2)-point(2)*Vector1(1))/denom(2)
ELSE
  CalcVectorAdditionCoeffs(2)=(point(3)*Vector1(1)-point(1)*Vector1(3))/denom(3)
END IF

IF (ABS(Vector1(1)).GT.ABS(Vector1(2)) .AND. ABS(Vector1(1)).GT.ABS(Vector1(3))) THEN
  CalcVectorAdditionCoeffs(1)=(point(1)-CalcVectorAdditionCoeffs(2)*Vector2(1))/Vector1(1)
ELSE IF (ABS(Vector1(2)).GT.ABS(Vector1(1)) .AND. ABS(Vector1(2)).GT.ABS(Vector1(3))) THEN
  CalcVectorAdditionCoeffs(1)=(point(2)-CalcVectorAdditionCoeffs(2)*Vector2(2))/Vector1(2)
ELSE
  CalcVectorAdditionCoeffs(1)=(point(3)-CalcVectorAdditionCoeffs(2)*Vector2(3))/Vector1(3)
END IF

END FUNCTION CalcVectorAdditionCoeffs
#endif /*CODE_ANALYZE*/

PURE FUNCTION CalcIntensity_Gaussian(x,x_norm)
!===================================================================================================================================
!> Calculates an exponential function of the Gaussian form for a given input variable (e.g. time) and a normalization variable
!> (e.g. pulse duration)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)         :: x, x_norm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: CalcIntensity_Gaussian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CalcIntensity_Gaussian = EXP(-(x/x_norm)**2)

END FUNCTION CalcIntensity_Gaussian


SUBROUTINE CalcNbrOfPhotons(i,iInit,NbrOfPhotons)
!===================================================================================================================================
!> Calculation of number of photons based on an intensity function integral of timestep and radius
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars  ,ONLY: PI
USE MOD_Particle_Vars ,ONLY: Species
USE MOD_Timedisc_Vars ,ONLY: dt,time
#ifdef LSERK
USE MOD_Timedisc_Vars ,ONLY: iStage, RK_c, nRKStages
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(IN)   :: iInit
! INOUTPUT VARIABLES
REAL, INTENT(OUT)     :: NbrOfPhotons
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                  :: t_1, t_2
REAL                  :: E_Intensity
INTEGER               :: NbrOfRepetitions
!===================================================================================================================================

ASSOCIATE( tau         => Species(i)%Init(iInit)%PulseDuration      ,&
           tShift      => Species(i)%Init(iInit)%tShift             ,&
           w_b         => Species(i)%Init(iInit)%WaistRadius        ,&
           Radius      => Species(i)%Init(iInit)%RadiusIC           ,&
           I_0         => Species(i)%Init(iInit)%IntensityAmplitude ,&
           lambda      => Species(i)%Init(iInit)%WaveLength         ,&
           NbrOfPulses => Species(i)%Init(iInit)%NbrOfPulses        ,&
           Period      => Species(i)%Init(iInit)%Period              &
          )

! Calculate the current pulse
NbrOfRepetitions = INT(Time/Period)

! Temporal bound of integration
#ifdef LSERK
IF (iStage.EQ.1) THEN
t_1 = Time
t_2 = Time + RK_c(2) * dt
ELSE
  IF (iStage.NE.nRKStages) THEN
    t_1 = Time + RK_c(iStage) * dt
    t_2 = Time + RK_c(iStage+1) * dt
  ELSE
    t_1 = Time + RK_c(iStage) * dt
    t_2 = Time + dt
  END IF
END IF
#else
t_1 = Time
t_2 = Time + dt
#endif
! Add arbitrary time shift (-4 sigma_t) so that I_max is not at t=0s
! Note that sigma_t = tau / sqrt(2)

t_1 = t_1 - tShift - NbrOfRepetitions * Period
t_2 = t_2 - tShift - NbrOfRepetitions * Period
! check if t_2 is outside of the pulse
IF(t_2.GT.2.0*tShift) t_2 = 2.0*tShift

! Integral of I(r,t) = I_0 exp(-(t/tau)**2)exp(-(r/w_b)**2)
! Integrate I(r,t)*r dr dt dphi and don't forget the Jacobian
!   dr : from 0 to R
!   dt : from t1 to t2
! dphi : from 0 to 2*PI
E_Intensity = 0.5 * I_0 * PI**(3.0/2.0) * w_b**2 * tau &
            * (1.0-EXP(-Radius**2/w_b**2)) &
            * (ERF(t_2/tau)-ERF(t_1/tau))
NbrOfPhotons = E_Intensity / CalcPhotonEnergy(lambda)

END ASSOCIATE

END SUBROUTINE CalcNbrOfPhotons


PURE FUNCTION CalcPhotonEnergy(lambda)
!===================================================================================================================================
!> Calculation of photon energy based on wavelength
!===================================================================================================================================
! MODULES
USE MOD_Globals_Vars          ,ONLY: PlanckConst,c
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)         :: lambda
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                     :: CalcPhotonEnergy
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

CalcPhotonEnergy = PlanckConst * c / lambda

END FUNCTION CalcPhotonEnergy


SUBROUTINE CalcVelocity_FromWorkFuncSEE(FractNbr, Vec3D, iInit)
!===================================================================================================================================
!> Subroutine to sample photon SEE electrons velocities from given energy distribution based on a work function.
!> Perform ARM for the energy distribution and a second ARM for the emission angle (between the impacting photon and the emitting
!> electron)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars            ,ONLY: PI, ElementaryCharge
USE MOD_Particle_Vars           ,ONLY: Species
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN), OPTIONAL     :: iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: RandVal
REAL               :: E_temp, E_max, VeloABS
REAL               :: Theta, Chi!, Psi_temp
REAL               :: PDF_temp, PDF_max
REAL, PARAMETER    :: PDF_max2=4./ACOS(-1.)
REAL               :: VeloVec_norm(3), RotationAxi(3)
LOGICAL            :: ARM_SEE_PDF
REAL               :: Theta_temp
!===================================================================================================================================

ASSOCIATE( W     => Species(FractNbr)%Init(iInit)%WorkFunctionSEE    ,&
           m     => Species(FractNbr)%MassIC                         ,&
           !beta  => Species(FractNbr)%Init(iInit)%AngularBetaSEE     ,&
           t_vec => Species(FractNbr)%Init(iInit)%BaseVector1IC      ,&
           n_vec => Species(FractNbr)%Init(iInit)%NormalIC            )

! ARM for energy distribution
E_max = 50.0 ! in eV (arbitrary)
PDF_max = 81.0 / (128.0 * W)  ! PDF_max at E = W/3 (derivation of 6W^2E/(E+W)^4 == 0)
ARM_SEE_PDF=.TRUE.
DO WHILE(ARM_SEE_PDF)
  CALL RANDOM_NUMBER(RandVal)
  E_temp = RandVal * E_max
  PDF_temp = 6 * W**2 * E_temp / (E_temp + W)**4
  CALL RANDOM_NUMBER(RandVal)
  IF ((PDF_temp/PDF_max).GT.RandVal) ARM_SEE_PDF = .FALSE.
END DO
VeloABS = SQRT(2.0 * E_temp * ElementaryCharge / m)

! ARM for angular distribution
CALL RANDOM_NUMBER(RandVal)
Chi = RandVal * 2.0 * PI
!PDF_max = -(2.0*(beta+4.0)) / (PI * (beta-8.0)) ! Henke 1977
!PDF_max2 = 4. / PI
ARM_SEE_PDF=.TRUE.
DO WHILE(ARM_SEE_PDF)
  CALL RANDOM_NUMBER(RandVal)
  !Psi_temp =0.5 * PI * (1.0 + RandVal)
  !Psi_temp = 0.5 * PI * RandVal
  !PDF_temp = ( (1.0 - beta / 2.0) + 3.0/4.0*beta*SIN(Psi_temp)**2 ) / (PI-PI*beta/8.0) ! Henke 1977
  Theta_temp = RandVal * 0.5 * PI
  PDF_temp = 4.0 / PI * COS(Theta_temp)**2
  CALL RANDOM_NUMBER(RandVal)
  IF ((PDF_temp/PDF_max2).GT.RandVal) ARM_SEE_PDF = .FALSE.
END DO
!Theta = PI - Psi_temp ! Henke 1977
Theta = Theta_temp

! Construct norm. VeloVec based on n_vec, t_vec, Theta and Chi
! first:  rotation of t_vec about n_vec with Chi (anzimuthal)
RotationAxi = n_vec
VeloVec_norm = t_vec * COS(Chi) + CROSS(RotationAxi,t_vec) * SIN(Chi)
! second: rotation of VeloVec_norm about RotationAxi with Theta (pi/2-theta)
RotationAxi = CROSS(VeloVec_norm,n_vec)
VeloVec_norm = VeloVec_norm * SIN(Theta) + CROSS(RotationAxi,VeloVec_norm) * COS(Theta)

! Calc VeloVec
Vec3D = VeloVec_norm * VeloABS

END ASSOCIATE

END SUBROUTINE CalcVelocity_FromWorkFuncSEE


SUBROUTINE SetParticlePositionPhotonSEEDisc(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),radius
INTEGER                 :: i
REAL                    :: lineVector(3),lineVector2(3)
LOGICAL                 :: ARM_Gauss
REAL                    :: RandVal(2),RandVal1
!===================================================================================================================================
! Use the base vectors BaseVector1IC and BaseVector2IC as coordinate system (they must be perpendicular)
lineVector = UNITVECTOR(Species(FractNbr)%Init(iInit)%BaseVector1IC(1:3))
lineVector2 = UNITVECTOR(Species(FractNbr)%Init(iInit)%BaseVector2IC(1:3))

DO i=1,chunkSize
  radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
  ARM_Gauss = .TRUE.
  DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC).OR.(ARM_Gauss))
    CALL RANDOM_NUMBER(RandVal)
    ! Check if particles are to be inserted in the first quadrant only, otherwise map R [0,1] -> R [-1,1]
    IF(.NOT.Species(FractNbr)%Init(iInit)%FirstQuadrantOnly) RandVal = RandVal * 2. - 1.
    Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
        (RandVal(1) * lineVector + RandVal(2) *lineVector2)

    radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
        (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
        (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
        (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
        (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
        (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
    ! Start ARM for Gauss distribution
    CALL RANDOM_NUMBER(RandVal1)
    IF(CalcIntensity_Gaussian(radius,Species(FractNbr)%Init(iInit)%WaistRadius).GT.RandVal1) ARM_Gauss = .FALSE.
    ! End ARM for Gauss distribution
  END DO
  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO
END SUBROUTINE SetParticlePositionPhotonSEEDisc


SUBROUTINE SetParticlePositionPhotonCylinder(FractNbr,iInit,chunkSize,particle_positions)
!===================================================================================================================================
! Set particle position
!===================================================================================================================================
! modules
USE MOD_Globals
USE MOD_Particle_Vars          ,ONLY: Species
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)     :: FractNbr, iInit, chunkSize
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)       :: particle_positions(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: Particle_pos(3),radius
INTEGER                 :: i
REAL                    :: lineVector(3),lineVector2(3)
LOGICAL                 :: ARM_Gauss
REAL                    :: RandVal(2),RandVal1
!===================================================================================================================================
! Use the base vectors BaseVector1IC and BaseVector2IC as coordinate system (they must be perpendicular)
lineVector = UNITVECTOR(Species(FractNbr)%Init(iInit)%BaseVector1IC(1:3))
lineVector2 = UNITVECTOR(Species(FractNbr)%Init(iInit)%BaseVector2IC(1:3))

DO i=1,chunkSize
  radius = Species(FractNbr)%Init(iInit)%RadiusIC + 1
  ARM_Gauss = .TRUE.
  DO WHILE((radius.GT.Species(FractNbr)%Init(iInit)%RadiusIC).OR.(ARM_Gauss))
    CALL RANDOM_NUMBER(RandVal)
    ! Check if particles are to be inserted in the first quadrant only, otherwise map R [0,1] -> R [-1,1]
    IF(.NOT.Species(FractNbr)%Init(iInit)%FirstQuadrantOnly) RandVal = RandVal * 2. - 1.
    Particle_pos = Species(FractNbr)%Init(iInit)%BasePointIC + Species(FractNbr)%Init(iInit)%RadiusIC * &
        (RandVal(1) * lineVector + RandVal(2) *lineVector2)

    radius = SQRT( (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) * &
        (Particle_pos(1)-Species(FractNbr)%Init(iInit)%BasePointIC(1)) + &
        (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) * &
        (Particle_pos(2)-Species(FractNbr)%Init(iInit)%BasePointIC(2)) + &
        (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) * &
        (Particle_pos(3)-Species(FractNbr)%Init(iInit)%BasePointIC(3)) )
    ! Start ARM for Gauss distribution
    CALL RANDOM_NUMBER(RandVal1)
    IF(CalcIntensity_Gaussian(radius,Species(FractNbr)%Init(iInit)%WaistRadius).GT.RandVal1) ARM_Gauss = .FALSE.
    ! End ARM for Gauss distribution
  END DO
  CALL RANDOM_NUMBER(RandVal1)
  Particle_pos = Particle_pos &
      + Species(FractNbr)%Init(iInit)%NormalIC * Species(FractNbr)%Init(iInit)%CylinderHeightIC * RandVal1
  particle_positions(i*3-2) = Particle_pos(1)
  particle_positions(i*3-1) = Particle_pos(2)
  particle_positions(i*3  ) = Particle_pos(3)
END DO
END SUBROUTINE SetParticlePositionPhotonCylinder


END MODULE MOD_part_emission_tools
