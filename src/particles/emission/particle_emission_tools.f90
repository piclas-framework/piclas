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

INTERFACE CalcVelocity_emmert
  MODULE PROCEDURE CalcVelocity_emmert
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

#if CODE_ANALYZE
INTERFACE CalcVectorAdditionCoeffs
  MODULE PROCEDURE CalcVectorAdditionCoeffs
END INTERFACE
#endif /*CODE_ANALYZE*/

!===================================================================================================================================
PUBLIC :: CalcVelocity_taylorgreenvortex, CalcVelocity_emmert
PUBLIC :: IntegerDivide,SetParticleChargeAndMass,SetParticleMPF,CalcVelocity_maxwell_lpn,SamplePoissonDistri
PUBLIC :: BessK,DEVI,SYNGE,QUASIREL
PUBLIC :: SetCellLocalParticlePosition,InsideExcludeRegionCheck
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


SUBROUTINE SetParticleChargeAndMass(FractNbr,NbrOfParticle)
!===================================================================================================================================
! And partilces mass and charge
!===================================================================================================================================
! MODULES
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

IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
  NbrOfParticle = PDM%maxParticleNumber
END IF
i = 1
DO WHILE (i .le. NbrOfParticle)
  PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
  IF (PositionNbr .ne. 0) THEN
    PartSpecies(PositionNbr) = FractNbr
  END IF
  i = i + 1
END DO

END SUBROUTINE SetParticleChargeAndMass


SUBROUTINE SetParticleMPF(FractNbr,NbrOfParticle)
!===================================================================================================================================
! finally, set the MPF
!===================================================================================================================================
! MODULES
USE MOD_Particle_Vars,    ONLY : PDM, PartMPF, Species
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

IF(NbrOfParticle.gt.PDM%maxParticleNumber)THEN
  NbrOfParticle = PDM%maxParticleNumber
END IF
i = 1
DO WHILE (i .le. NbrOfParticle)
  PositionNbr = PDM%nextFreePosition(i+PDM%CurrentNextFreePosition)
  IF (PositionNbr .ne. 0) THEN
    PartMPF(PositionNbr) = Species(FractNbr)%MacroParticleFactor
  END IF
  i = i + 1
END DO

END SUBROUTINE SetParticleMPF


SUBROUTINE CalcVelocity_maxwell_lpn(FractNbr, Vec3D, iInit, Element, Temperature)
!===================================================================================================================================
! Subroutine to sample current cell values (partly copied from 'LD_DSMC_Mean_Bufferzone_A_Val' and 'dsmc_analyze')
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species!, DoZigguratSampling
!USE Ziggurat,                   ONLY : rnor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN), OPTIONAL     :: iInit
INTEGER, OPTIONAL                :: Element !for BGG from VTK
REAL,INTENT(IN), OPTIONAL        :: Temperature
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RandVal(3), Velo1, Velo2, Velosq, Tx, ty, Tz, v_drift(3)
!===================================================================================================================================
IF(PRESENT(iInit).AND.PRESENT(Temperature))CALL abort(&
__STAMP__&
,'CalcVelocity_maxwell_lpn. iInit and Temperature cannot both be input arguments!')
IF(PRESENT(iInit).AND..NOT.PRESENT(Element))THEN
  Tx=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Ty=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  Tz=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  v_drift=Species(FractNbr)%Init(iInit)%VeloIC *Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
ELSE IF (PRESENT(Element)) THEN
  IF (Species(FractNbr)%Init(iInit)%ElemTemperatureFileID.GT.0) THEN
    Tx=Species(FractNbr)%Init(iInit)%ElemTemperatureIC(1,Element)
    Ty=Species(FractNbr)%Init(iInit)%ElemTemperatureIC(2,Element)
    Tz=Species(FractNbr)%Init(iInit)%ElemTemperatureIC(3,Element)
  ELSE
    Tx=Species(FractNbr)%Init(iInit)%MWTemperatureIC
    Ty=Species(FractNbr)%Init(iInit)%MWTemperatureIC
    Tz=Species(FractNbr)%Init(iInit)%MWTemperatureIC
  END IF
  IF (Species(FractNbr)%Init(iInit)%ElemVelocityICFileID.GT.0) THEN
    v_drift=Species(FractNbr)%Init(iInit)%ElemVelocityIC(1:3,Element)
  ELSE
    v_drift=Species(FractNbr)%Init(iInit)%VeloIC *Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
  END IF
ELSE IF(PRESENT(Temperature))THEN
  Tx=Temperature
  Ty=Temperature
  Tz=Temperature
  v_drift=0.0
ELSE
CALL abort(&
__STAMP__&
,'PO: force temperature!!')
END IF

!IF (.NOT.DoZigguratSampling) THEN !polar method
  Velosq = 2
  DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    CALL RANDOM_NUMBER(RandVal)
    Velo1 = 2.*RandVal(1) - 1.
    Velo2 = 2.*RandVal(2) - 1.
    Velosq = Velo1**2 + Velo2**2
  END DO
  Vec3D(1) = Velo1*SQRT(-2*BoltzmannConst*Tx/ &
    Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !x-Komponente
  Vec3D(2) = Velo2*SQRT(-2*BoltzmannConst*Ty/ &
  Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !y-Komponente
  Velosq = 2
  DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    CALL RANDOM_NUMBER(RandVal)
    Velo1 = 2.*RandVal(1) - 1.
    Velo2 = 2.*RandVal(2) - 1.
    Velosq = Velo1**2 + Velo2**2
  END DO
  Vec3D(3) = Velo1*SQRT(-2*BoltzmannConst*Tz/ &
    Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !z-Komponente
!ELSE !ziggurat method
!  Velo1 = rnor()
!  Vec3D(1) = Velo1*SQRT(BoltzmannConst*Tx/Species(FractNbr)%MassIC)             !x-Komponente
!  Velo1 = rnor()
!  Vec3D(2) = Velo1*SQRT(BoltzmannConst*Ty/Species(FractNbr)%MassIC)             !y-Komponente
!  Velo1 = rnor()
!  Vec3D(3) = Velo1*SQRT(BoltzmannConst*Tz/Species(FractNbr)%MassIC)             !z-Komponente
!END IF
Vec3D(1:3) = Vec3D(1:3) + v_drift

END SUBROUTINE CalcVelocity_maxwell_lpn


SUBROUTINE CalcVelocity_taylorgreenvortex(FractNbr, Vec3D, iInit, Element)
!===================================================================================================================================
! Subroutine to sample current cell values (partly copied from 'LD_DSMC_Mean_Bufferzone_A_Val' and 'dsmc_analyze')
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars       ,ONLY: BoltzmannConst
USE MOD_Particle_Vars      ,ONLY: Species
USE MOD_Mesh_Vars          ,ONLY: ElemBaryNGeo
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
INTEGER,INTENT(IN)               :: FractNbr
INTEGER,INTENT(IN), OPTIONAL     :: iInit
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


SUBROUTINE CalcVelocity_emmert(FractNbr, iInit, Vec3D)
!===================================================================================================================================
! Subroutine to sample particle velos in VecIC from distri by Emmert et al. [Phys. Fluids 23, 803 (1980)] and in normal dir. from MB
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : Species!, DoZigguratSampling
!USE Ziggurat,                   ONLY : rnor
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: FractNbr, iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                 :: Vec3D(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: RandVal(3), Velo1, Velo2, Velosq, T, v_dir(3), vec_t1(3), vec_t2(3), v_d
!===================================================================================================================================

T=Species(FractNbr)%Init(iInit)%MWTemperatureIC
v_dir=Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
v_d=Species(FractNbr)%Init(iInit)%VeloIC

!--build arbitrary vectors normal to v_dir
IF (.NOT.ALMOSTEQUAL(v_dir(3),0.)) THEN
  vec_t1(1) = 1.0
  vec_t1(2) = 1.0
  vec_t1(3) = -(v_dir(1)+v_dir(2))/v_dir(3)
  vec_t2(1) = v_dir(2) * vec_t1(3) - v_dir(3)
  vec_t2(2) = v_dir(3) - v_dir(1) * vec_t1(3)
  vec_t2(3) = v_dir(1) - v_dir(2)
  vec_t1 = vec_t1 / SQRT(2.0 + vec_t1(3)*vec_t1(3))
ELSE
  IF (.NOT.ALMOSTEQUAL(v_dir(2),0.)) THEN
    vec_t1(1) = 1.0
    vec_t1(3) = 1.0
    vec_t1(2) = -(v_dir(1)+v_dir(3))/v_dir(2)
    vec_t2(1) = v_dir(2) - v_dir(3) * vec_t1(2)
    vec_t2(2) = v_dir(3) - v_dir(1)
    vec_t2(3) = v_dir(1) * vec_t1(2) - v_dir(2)
    vec_t1 = vec_t1 / SQRT(2.0 + vec_t1(2)*vec_t1(2))
  ELSE
    IF (.NOT.ALMOSTEQUAL(v_dir(1),0.)) THEN
      vec_t1(2) = 1.0
      vec_t1(3) = 1.0
      vec_t1(1) = -(v_dir(2)+v_dir(3))/v_dir(1)
      vec_t2(1) = v_dir(2) - v_dir(3)
      vec_t2(2) = v_dir(3) * vec_t1(1) - v_dir(1)
      vec_t2(3) = v_dir(1) - v_dir(2) * vec_t1(1)
      vec_t1 = vec_t1 / SQRT(2.0 + vec_t1(1)*vec_t1(1))
    ELSE
      CALL abort(&
__STAMP__&
,'Error in CalcVelocity_emmert, VeloVecIC is zero!')
    END IF
  END IF
END IF
vec_t2 = vec_t2 / SQRT(vec_t2(1)*vec_t2(1) + vec_t2(2)*vec_t2(2) + vec_t2(3)*vec_t2(3))

!--sample velocities
!IF (.NOT.DoZigguratSampling) THEN !polar method
  Velosq = 2
  DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
    CALL RANDOM_NUMBER(RandVal)
    Velo1 = 2.*RandVal(1) - 1.
    Velo2 = 2.*RandVal(2) - 1.
    Velosq = Velo1**2 + Velo2**2
  END DO
  Vec3D(1:3) =              vec_t1(1:3)*Velo1*SQRT(-2*BoltzmannConst*T/ &
    Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !n1-Komponente (maxwell_lpn)
  Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3)*Velo2*SQRT(-2*BoltzmannConst*T/ &
    Species(FractNbr)%MassIC*LOG(Velosq)/Velosq)                                !n2-Komponente (maxwell_lpn)
!ELSE !ziggurat method
!  Velo1=rnor()
!  Vec3D(1:3) =              vec_t1(1:3)*Velo1*SQRT(BoltzmannConst*T/ &
!    Species(FractNbr)%MassIC)                                !n1-Komponente (maxwell_lpn)
!  Velo2=rnor()
!  Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3)*Velo2*SQRT(BoltzmannConst*T/ &
!    Species(FractNbr)%MassIC)                                !n2-Komponente (maxwell_lpn)
!END IF
Vec3D(1:3) = Vec3D(1:3) + v_dir(1:3)*SQRT(BoltzmannConst*T/Species(FractNbr)%MassIC)* &                ! (emmert)
  sign(1.d0,RandVal(3)-0.5d0)*SQRT(-2*log(1-sign(1.d0,RandVal(3)-0.5d0)*(2*RandVal(3)-1)))

Vec3D(1:3) = Vec3D(1:3) + v_dir(1:3)*v_d

END SUBROUTINE CalcVelocity_emmert


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
USE MOD_Globals_Vars,  ONLY: BoltzmannConst
USE MOD_Equation_Vars,  ONLY: c2
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
USE MOD_Globals_Vars,   ONLY: BoltzmannConst

USE MOD_Equation_Vars,  ONLY: c_inv,c2
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
USE MOD_Globals_Vars,  ONLY: BoltzmannConst
USE MOD_Equation_Vars,  ONLY: c_inv,c2
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
USE MOD_DSMC_Symmetry2D        ,ONLY: CalcRadWeightMPF
USE MOD_Eval_xyz               ,ONLY: GetPositionInRefElem
USE MOD_MacroBody_Vars         ,ONLY: UseMacroBody
USE MOD_MacroBody_Tools        ,ONLY: INSIDEMACROBODY
USE MOD_Mesh_Vars              ,ONLY: nElems,offsetElem
USE MOD_Particle_Localization  ,ONLY: PartInElemCheck
USE MOD_Particle_Mesh_Vars     ,ONLY: LocalVolume
USE MOD_Particle_Mesh_Vars     ,ONLY: GEO, ElemEpsOneCell
USE MOD_Particle_Mesh_Tools    ,ONLY: ParticleInsideQuad3D
USE MOD_Particle_Tracking_Vars ,ONLY: DoRefMapping, TriaTracking
USE MOD_Particle_Vars          ,ONLY: Species, PDM, PartState, PEM, Symmetry2D, Symmetry2DAxisymmetric, VarTimeStep, PartMPF
USE MOD_Particle_VarTimeStep   ,ONLY: CalcVarTimeStep
#if USE_MPI
USE MOD_MPI_Shared_Vars        ,ONLY: BoundsOfElem_Shared,ElemVolume_Shared
#else
USE MOD_Mesh_Vars              ,ONLY: BoundsOfElem_Shared,ElemVolume_Shared
#endif /*USE_MPI*/
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
INTEGER                          :: CellChunkSize(1:nElems)
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
    IF (Species(iSpec)%Init(iInit)%ElemPartDensityFileID.EQ.0) THEN
      CALL IntegerDivide(chunkSize,nElems,ElemVolume_Shared(:),CellChunkSize(:))
    ELSE
      CALL IntegerDivide(chunkSize,nElems,Species(iSpec)%Init(iInit)%ElemPartDensity(:)*ElemVolume_Shared(:),CellChunkSize(:))
    END IF
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
  DO iElem = 1, nElems
    !ASSOCIATE( Bounds => GEO%BoundsOfElem(1:2,1:3,iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
    ASSOCIATE( Bounds => BoundsOfElem_Shared(1:2,1:3,offsetElem+iElem) ) ! 1-2: Min, Max value; 1-3: x,y,z
      IF (UseExactPartNum) THEN
        nPart = CellChunkSize(iElem)
      ELSE
        IF(RadialWeighting%DoRadialWeighting) THEN
          PartDens = Species(iSpec)%Init(iInit)%PartDensity / CalcRadWeightMPF(GEO%ElemMidPoint(2,iElem), iSpec)
        END IF
        CALL RANDOM_NUMBER(iRan)
        IF(VarTimeStep%UseVariableTimeStep) THEN
          adaptTimestep = CalcVarTimeStep(GEO%ElemMidPoint(1,iElem), GEO%ElemMidPoint(2,iElem), iElem)
          nPart = INT(PartDens / adaptTimestep * ElemVolume_Shared(iElem) + iRan)
        ELSE
          nPart = INT(PartDens * ElemVolume_Shared(iElem) + iRan)
        END IF
      END IF
      DO iPart = 1, nPart
        ParticleIndexNbr = PDM%nextFreePosition(iChunksize + PDM%CurrentNextFreePosition)
        IF (ParticleIndexNbr .ne. 0) THEN
          InsideFlag=.FALSE.
          DO WHILE(.NOT.InsideFlag)
            CALL RANDOM_NUMBER(RandomPos)
            IF(Symmetry2DAxisymmetric.AND.(.NOT.RadialWeighting%DoRadialWeighting)) THEN
              ! Treatment of axisymmetry without weighting
              RandomPos(1) = Bounds(1,1) + RandomPos(1)*(Bounds(2,1)-Bounds(1,1))
              RandomPos(2) = SQRT(RandomPos(2)*(Bounds(2,2)**2-Bounds(1,2)**2)+Bounds(1,2)**2)
            ELSE
              RandomPos = Bounds(1,:) + RandomPos*(Bounds(2,:)-Bounds(1,:))
            END IF
            IF(Symmetry2D) RandomPos(3) = 0.
            IF (DoRefMapping) THEN
              CALL GetPositionInRefElem(RandomPos,RefPos,iElem)
              IF (MAXVAL(ABS(RefPos)).GT.ElemEpsOneCell(iElem)) InsideFlag=.TRUE.
            ELSE
              IF (TriaTracking) THEN
                CALL ParticleInsideQuad3D(RandomPos,iElem,InsideFlag,Det)
              ELSE
                CALL PartInElemCheck(RandomPos,iPart,iElem,InsideFlag)
              END IF
            END IF
          END DO
          IF (UseMacroBody) THEN
            IF (INSIDEMACROBODY(RandomPos)) THEN
              CYCLE !particle is inside MacroParticle
            END IF
          END IF
          PartState(1:3,ParticleIndexNbr) = RandomPos(1:3)
          PDM%ParticleInside(ParticleIndexNbr) = .TRUE.
          PDM%IsNewPart(ParticleIndexNbr)=.TRUE.
          PDM%dtFracPush(ParticleIndexNbr) = .FALSE.
          PEM%Element(ParticleIndexNbr) = iElem
          ichunkSize = ichunkSize + 1
          IF (VarTimeStep%UseVariableTimeStep) THEN
            VarTimeStep%ParticleTimeStep(ParticleIndexNbr) = &
              CalcVarTimeStep(PartState(1,ParticleIndexNbr), PartState(2,ParticleIndexNbr),iElem)
          END IF
          IF(RadialWeighting%DoRadialWeighting) THEN
            PartMPF(ParticleIndexNbr) = CalcRadWeightMPF(PartState(2,ParticleIndexNbr),1,ParticleIndexNbr)
          END IF
        ELSE
          CALL abort(&
              __STAMP__&
              ,'ERROR in SetCellLocalParticlePosition: Maximum particle number reached during inserting! --> ParticleIndexNbr.EQ.0')
        END IF
      END DO
    END ASSOCIATE
  END DO
  chunkSize = ichunkSize - 1

END SUBROUTINE SetCellLocalParticlePosition


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


END MODULE MOD_part_emission_tools
