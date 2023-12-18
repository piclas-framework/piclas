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

MODULE MOD_DSMC_AmbipolarDiffusion
!===================================================================================================================================
! Module for use of a background gas for the simulation of trace species (if number density of bg gas is multiple orders of
! magnitude larger than the trace species)
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
PUBLIC :: InitializeVariablesAmbipolarDiff, AD_SetInitElectronVelo, AD_InsertParticles, AD_DeleteParticles, AD_SetSFElectronVelo
!===================================================================================================================================

CONTAINS

SUBROUTINE InitializeVariablesAmbipolarDiff()
!===================================================================================================================================
!> Ambipolar Diffusion: Electrons are attached to and move with the ions, but still have their own velocity vector for collisions
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Globals_Vars        ,ONLY: ElementaryCharge
USE MOD_Particle_Vars       ,ONLY: nSpecies,Species, PDM
USE MOD_DSMC_Vars           ,ONLY: useDSMC, DSMC, AmbipolElecVelo
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: iSpec
!===================================================================================================================================
IF(useDSMC) THEN
  DSMC%DoAmbipolarDiff = GETLOGICAL('Particles-DSMC-AmbipolarDiffusion')
  IF (DSMC%DoAmbipolarDiff) THEN
    DSMC%AmbiDiffElecSpec = 0
    DO iSpec = 1, nSpecies
      IF (Species(iSpec)%ChargeIC.GT.0.0) CYCLE
      IF(NINT(Species(iSpec)%ChargeIC/(-ElementaryCharge)).EQ.1) DSMC%AmbiDiffElecSpec=iSpec
    END DO
    IF(DSMC%AmbiDiffElecSpec.EQ.0) THEN
      CALL abort(__STAMP__&
          ,'ERROR: No electron species found for ambipolar diffusion: ' &
          ,IntInfoOpt=DSMC%AmbiDiffElecSpec)
    END IF
    IF(.NOT.ALLOCATED(AmbipolElecVelo)) ALLOCATE(AmbipolElecVelo(PDM%maxParticleNumber))
  END IF
END IF

END SUBROUTINE InitializeVariablesAmbipolarDiff


SUBROUTINE AD_SetInitElectronVelo(FractNbr,iInit,NbrOfParticle)
!===================================================================================================================================
!> Initialize the electron velocity vector during the initial particle insertion
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars
USE MOD_Globals_Vars            ,ONLY: BoltzmannConst
USE MOD_part_emission_tools     ,ONLY: CalcVelocity_maxwell_lpn
USE MOD_part_tools              ,ONLY: BuildTransGaussNums, GetNextFreePosition
USE MOD_DSMC_Vars               ,ONLY: DSMC, AmbipolElecVelo
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: FractNbr,iInit
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)           :: NbrOfParticle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! LOCAL VARIABLES
INTEGER                         :: i, PositionNbr
CHARACTER(30)                   :: velocityDistribution
REAL                            :: VeloIC, VeloVecIC(3), maxwellfac, VeloVecNorm
REAL                            :: iRanPart(3, NbrOfParticle), Vec3D(3)
!===================================================================================================================================
IF(NbrOfParticle.LT.1) RETURN
IF(Species(FractNbr)%ChargeIC.LE.0.0) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PDM%maxParticleNumber!')
END IF

velocityDistribution=Species(FractNbr)%Init(iInit)%velocityDistribution
VeloIC=Species(FractNbr)%Init(iInit)%VeloIC
VeloVecIC=Species(FractNbr)%Init(iInit)%VeloVecIC(1:3)
VeloVecNorm = VECNORM(VeloVecIC(1:3))
IF (VeloVecNorm.GT.0.0) THEN
  VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
END IF

SELECT CASE(TRIM(velocityDistribution))
CASE('constant')
  DO i = 1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
    AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloVecIC(1:3) * VeloIC
  END DO
CASE('maxwell_lpn')
  DO i = 1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    CALL CalcVelocity_maxwell_lpn(DSMC%AmbiDiffElecSpec, Vec3D, Temperature=Species(FractNbr)%Init(iInit)%MWTemperatureIC)
    IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
    AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + Vec3D(1:3)
  END DO
CASE('maxwell')
  CALL BuildTransGaussNums(NbrOfParticle, iRanPart)
  maxwellfac = SQRT(BoltzmannConst*Species(FractNbr)%Init(iInit)%MWTemperatureIC/Species(DSMC%AmbiDiffElecSpec)%MassIC)
  DO i = 1,NbrOfParticle
    PositionNbr = GetNextFreePosition(i)
    IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
    AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = VeloIC *VeloVecIC(1:3) + iRanPart(1:3,i)*maxwellfac
  END DO
CASE DEFAULT
  CALL abort(&
__STAMP__&
,'Velo-Distri not implemented for ambipolar diffusion!')
END SELECT

END SUBROUTINE AD_SetInitElectronVelo


SUBROUTINE AD_SetSFElectronVelo(iSpec,iSFIon,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns,particle_xis)
!===================================================================================================================================
!> Initialize the electron velocity vector during the surface flux particle insertion
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Globals_Vars,           ONLY : PI, BoltzmannConst
USE MOD_Particle_Vars
USE MOD_Particle_Surfaces_Vars, ONLY : SurfMeshSubSideData, TriaSurfaceFlux
USE MOD_Particle_Surfaces,      ONLY : CalcNormAndTangBezier
USE MOD_Particle_Sampling_Vars  ,ONLY: AdaptBCMapElemToSample, AdaptBCMacroVal
USE MOD_DSMC_Vars               ,ONLY: AmbiPolarSFMapping, AmbipolElecVelo, DSMC
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)               :: iSpec,iSFIon,iSample,jSample,iSide,BCSideID,SideID,ElemID,NbrOfParticle,PartIns
REAL,INTENT(IN)                  :: particle_xis(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,PositionNbr,envelope,currentBC,SampleElemID, iSF,iPart
REAL                             :: Vec3D(3), vec_nIn(1:3), vec_t1(1:3), vec_t2(1:3)
REAL                             :: a,zstar,RandVal1,RandVal2(2),RandVal3(3),u,RandN,RandN_save,Velo1,Velo2,Velosq,T,beta,z
LOGICAL                          :: RandN_in_Mem
REAL                             :: projFak                          ! VeloVecIC projected to inwards normal of tria
REAL                             :: Velo_t1                          ! Velo comp. of first orth. vector in tria
REAL                             :: Velo_t2                          ! Velo comp. of second orth. vector in tria
REAL                             :: VeloIC
REAL                             :: VeloVec(1:3)
REAL                             :: VeloVecIC(1:3),v_thermal, pressure
!===================================================================================================================================

IF(PartIns.LT.1) RETURN
IF(Species(iSpec)%ChargeIC.LE.0.0) RETURN
IF(NbrOfParticle.GT.PDM%maxParticleNumber)THEN
     CALL abort(&
__STAMP__&
,'NbrOfParticle > PDM%maxParticleNumber!')
END IF

iSF = AmbiPolarSFMapping(iSpec,iSFIon)
RandN_in_Mem=.FALSE.
envelope=-1
currentBC = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%BC

IF (.NOT.Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
  vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
END IF !.NOT.VeloIsNormal

IF(.NOT.Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%Adaptive) THEN
  VeloIC = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIC
  T = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%MWTemperatureIC
  a = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%a_nIn
  projFak = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%projFak
  Velo_t1 = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t1
  Velo_t2 = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%SurfFluxSubSideData(iSample,jSample,iSide)%Velo_t2
ELSE !Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%Adaptive
  SampleElemID = AdaptBCMapElemToSample(ElemID)
  SELECT CASE(Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%AdaptiveType)
  CASE(1,3,4) ! Pressure and massflow inlet (pressure/massflow, temperature const)
    T =  Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%MWTemperatureIC
  CASE(2) ! adaptive Outlet/freestream
    pressure = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%AdaptivePressure
    T = pressure / (BoltzmannConst * AdaptBCMacroVal(4,SampleElemID,DSMC%AmbiDiffElecSpec))
  CASE DEFAULT
    CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong adaptive type for Surfaceflux velocities!')
  END SELECT
  VeloVec(1) = AdaptBCMacroVal(DSMC_VELOX,SampleElemID,DSMC%AmbiDiffElecSpec)
  VeloVec(2) = AdaptBCMacroVal(DSMC_VELOY,SampleElemID,DSMC%AmbiDiffElecSpec)
  VeloVec(3) = AdaptBCMacroVal(DSMC_VELOZ,SampleElemID,DSMC%AmbiDiffElecSpec)
  vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
  VeloVec(1:3) = DOT_PRODUCT(VeloVec,vec_nIn)*vec_nIn(1:3)
  VeloIC = SQRT(DOT_PRODUCT(VeloVec,VeloVec))
  IF (ABS(VeloIC).GT.0.) THEN
    VeloVecIC = VeloVec / VeloIC
  ELSE
    VeloVecIC = (/1.,0.,0./)
  END IF
  projFak = DOT_PRODUCT(vec_nIn,VeloVecIC) !VeloVecIC projected to inwards normal
  v_thermal = SQRT(2.*BoltzmannConst*T/Species(DSMC%AmbiDiffElecSpec)%MassIC) !thermal speed
  IF ( ALMOSTEQUAL(v_thermal,0.)) THEN
    v_thermal = 1.
  END IF
  a = VeloIC * projFak / v_thermal !speed ratio proj. to inwards n (can be negative!)
  Velo_t1 = VeloIC * DOT_PRODUCT(vec_t1,VeloVecIC) !v in t1-dir
  Velo_t2 = VeloIC * DOT_PRODUCT(vec_t2,VeloVecIC) !v in t2-dir
END IF !Adaptive SurfaceFlux

! Set velocities
SELECT CASE(TRIM(Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%velocityDistribution))
CASE('constant')
  IF (.NOT.Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
    VeloVecIC(1:3) = Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloVecIC(1:3)
    VeloVecIC(1:3) = VeloVecIC(1:3) / VECNORM(VeloVecIC(1:3))
  END IF
  iPart = 0
  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    iPart = iPart + 1
    PositionNbr = GetNextFreePosition(i)
    ! In case of side-normal velocities: calc n-vector at particle position, xi was saved in PartState(4:5)
    IF (Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      vec_t1(1:3) = 0. !dummy
      vec_t2(1:3) = 0. !dummy
    ELSE IF (Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      ! CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
      CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),xi=particle_xis(2*(iPart-1)+1),eta=particle_xis(2*(iPart-1)+2),SideID=SideID )
      vec_nIn(1:3) = -vec_nIn(1:3)
      vec_t1(1:3) = 0. !dummy
      vec_t2(1:3) = 0. !dummy
    ELSE
      vec_nIn(1:3) = VeloVecIC(1:3)
    END IF !VeloIsNormal
    ! Build complete velo-vector
    Vec3D(1:3) = vec_nIn(1:3) * Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIC
    ! PartState(4:6,PositionNbr) = Vec3D(1:3)
    IF (PositionNbr.GT.0) THEN
      IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
      ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
      AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = Vec3D(1:3)
    END IF
  END DO !i = ...NbrOfParticle
CASE('maxwell','maxwell_lpn')
  !-- determine envelope for most efficient ARM [Garcia and Wagner 2006, JCP217-2]
  IF (ALMOSTZERO(VeloIC*projFak)) THEN
    ! Rayleigh distri
    envelope = 0
  ELSE IF (-0.4.LT.a .AND. a.LT.1.3) THEN
    ! low speed flow
    IF (a.LE.0.) THEN
      envelope = 1
    ELSE
      envelope = 3
    END IF !choose envelope based on flow direction
  ELSE
    ! high speed / general flow
    IF (a.LT.0.) THEN
      envelope = 2
    ELSE
      envelope = 4
    END IF !choose envelope based on flow direction
  END IF !low speed / high speed / rayleigh flow

  iPart = 0
  DO i = NbrOfParticle-PartIns+1,NbrOfParticle
    iPart = iPart + 1
    PositionNbr = GetNextFreePosition(i)
    !-- 0a.: In case of side-normal velocities: calc n-/t-vectors at particle position, xi was saved in PartState(4:5)
    IF (Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal .AND. TriaSurfaceFlux) THEN
      vec_nIn(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_nIn(1:3)
      vec_t1(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t1(1:3)
      vec_t2(1:3) = SurfMeshSubSideData(iSample,jSample,BCSideID)%vec_t2(1:3)
    ELSE IF (Species(DSMC%AmbiDiffElecSpec)%Surfaceflux(iSF)%VeloIsNormal) THEN
      ! CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),tang1=vec_t1(1:3),tang2=vec_t2(1:3) &
      !   ,xi=PartState(4,PositionNbr),eta=PartState(5,PositionNbr),SideID=SideID )
      CALL CalcNormAndTangBezier( nVec=vec_nIn(1:3),tang1=vec_t1(1:3),tang2=vec_t2(1:3) &
        ,xi=particle_xis(2*(iPart-1)+1),eta=particle_xis(2*(iPart-1)+2),SideID=SideID )
      vec_nIn(1:3) = -vec_nIn(1:3)
    END IF !VeloIsNormal
    !-- 1.: determine zstar (initial generation of potentially too many RVu is for needed indentities of RVu used multiple times!
    SELECT CASE(envelope)
    CASE(0)
      CALL RANDOM_NUMBER(RandVal1)
      zstar = -SQRT(-LOG(RandVal1))
    CASE(1)
      DO
        CALL RANDOM_NUMBER(RandVal2)
        zstar = -SQRT(a*a-LOG(RandVal2(1)))
        IF ( -(a-zstar)/zstar .GT. RandVal2(2)) THEN
          EXIT
        END IF
      END DO
    CASE(2)
      z = 0.5*(a-SQRT(a*a+2.))
      beta  = a-(1.0-a)*(a-z)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        IF (EXP(-(beta*beta))/(EXP(-(beta*beta))+2.0*(a-z)*(a-beta)*EXP(-(z*z))).GT.RandVal3(1)) THEN
          zstar=-SQRT(beta*beta-LOG(RandVal3(2)))
          IF ( -(a-zstar)/zstar .GT. RandVal3(3)) THEN
            EXIT
          END IF
        ELSE
          zstar=beta+(a-beta)*RandVal3(2)
          IF ( (a-zstar)/(a-z)*EXP(z*z-(zstar*zstar)) .GT. RandVal3(3)) THEN
            EXIT
          END IF
        END IF
      END DO
    CASE(3)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        u = RandVal3(1)
        IF ( a*SQRT(PI)/(a*SQRT(PI)+1+a*a) .GT. u) THEN
!            IF (.NOT.DoZigguratSampling) THEN !polar method
            IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
              RandN = RandN_save
              RandN_in_Mem=.FALSE.
            ELSE
              Velosq = 2
              DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
                CALL RANDOM_NUMBER(RandVal2)
                Velo1 = 2.*RandVal2(1) - 1.
                Velo2 = 2.*RandVal2(2) - 1.
                Velosq = Velo1**2 + Velo2**2
              END DO
              RandN = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_in_Mem=.TRUE.
            END IF
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
          zstar = -1./SQRT(2.)*ABS(RandN)
          EXIT
        ELSE IF ( (a*SQRT(PI)+1.)/(a*SQRT(PI)+1+a*a) .GT. u) THEN
          zstar = -SQRT(-LOG(RandVal3(2)))
          EXIT
        ELSE
          zstar = (1.0-SQRT(RandVal3(2)))*a
          IF (EXP(-(zstar*zstar)).GT.RandVal3(3)) THEN
            EXIT
          END IF
        END IF
      END DO
    CASE(4)
      DO
        CALL RANDOM_NUMBER(RandVal3)
        IF (1.0/(2.0*a*SQRT(PI)+1.0).GT.RandVal3(1)) THEN
          zstar=-SQRT(-LOG(RandVal3(2)))
        ELSE
!            IF (.NOT.DoZigguratSampling) THEN !polar method
            IF (RandN_in_Mem) THEN !reusing second RandN form previous polar method
              RandN = RandN_save
              RandN_in_Mem=.FALSE.
            ELSE
              Velosq = 2
              DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
                CALL RANDOM_NUMBER(RandVal2)
                Velo1 = 2.*RandVal2(1) - 1.
                Velo2 = 2.*RandVal2(2) - 1.
                Velosq = Velo1**2 + Velo2**2
              END DO
              RandN = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_save = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
              RandN_in_Mem=.TRUE.
            END IF
!            ELSE !ziggurat method
!              RandN=rnor()
!            END IF
          zstar = 1./SQRT(2.)*RandN
        END IF
        IF ( (a-zstar)/a .GT. RandVal3(3)) THEN
          EXIT
        END IF
      END DO
    CASE DEFAULT
      CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong envelope in SetSurfacefluxVelocities!')
    END SELECT
    !-- 2.: sample normal directions and build complete velo-vector
    Vec3D(1:3) = vec_nIn(1:3) * SQRT(2.*BoltzmannConst*T/Species(DSMC%AmbiDiffElecSpec)%MassIC)*(a-zstar)
!      IF (.NOT.DoZigguratSampling) THEN !polar method
      Velosq = 2
      DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal2)
        Velo1 = 2.*RandVal2(1) - 1.
        Velo2 = 2.*RandVal2(2) - 1.
        Velosq = Velo1**2 + Velo2**2
      END DO
      Velo1 = Velo1*SQRT(-2*LOG(Velosq)/Velosq)
      Velo2 = Velo2*SQRT(-2*LOG(Velosq)/Velosq)
!      ELSE !ziggurat method
!        Velo1=rnor()
!        Velo2=rnor()
!      END IF
    Vec3D(1:3) = Vec3D(1:3) + vec_t1(1:3) * ( Velo_t1+Velo1*SQRT(BoltzmannConst*T/Species(DSMC%AmbiDiffElecSpec)%MassIC) )
    Vec3D(1:3) = Vec3D(1:3) + vec_t2(1:3) * ( Velo_t2+Velo2*SQRT(BoltzmannConst*T/Species(DSMC%AmbiDiffElecSpec)%MassIC) )
    IF (ALLOCATED(AmbipolElecVelo(PositionNbr)%ElecVelo)) DEALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo)
    ALLOCATE(AmbipolElecVelo(PositionNbr)%ElecVelo(3))
    AmbipolElecVelo(PositionNbr)%ElecVelo(1:3) = Vec3D(1:3)
  END DO !i = ...NbrOfParticle
CASE DEFAULT
  CALL abort(__STAMP__,'ERROR in SurfaceFlux: Wrong velocity distribution!')
END SELECT

END SUBROUTINE AD_SetSFElectronVelo


SUBROUTINE AD_InsertParticles(iPartIndx_Node, nPart, iPartIndx_NodeTotalAmbi, TotalPartNum)
!===================================================================================================================================
!> Creating electrons for each actual ion simulation particle, using the stored velocity vector for the electron
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_DSMC_Vars               ,ONLY: BGGas, CollisMode, DSMC, PartStateIntEn, AmbipolElecVelo, RadialWeighting
USE MOD_DSMC_Vars               ,ONLY: DSMCSumOfFormedParticles, newAmbiParts, iPartIndx_NodeNewAmbi
USE MOD_PARTICLE_Vars           ,ONLY: PDM, PartSpecies, PartState, PEM, Species, PartMPF, Symmetry, usevMPF
USE MOD_PARTICLE_Vars           ,ONLY: UseVarTimeStep, PartTimeStep
USE MOD_Particle_Tracking       ,ONLY: ParticleInsideCheck
USE MOD_Part_Tools              ,ONLY: GetNextFreePosition
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(INOUT)             :: nPart, TotalPartNum
INTEGER,INTENT(INOUT)             :: iPartIndx_Node(1:nPart)
INTEGER,INTENT(INOUT),ALLOCATABLE :: iPartIndx_NodeTotalAmbi(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNewPart, iPart, PositionNbr, iLoop, nNewElectrons, IonIndX(nPart), iElem, PartNum
REAL              :: MaxPos(3), MinPos(3), Vec3D(3), RandomPos(3)
LOGICAL           :: InsideFlag
!===================================================================================================================================

MaxPos = -HUGE(MaxPos)
MinPos = HUGE(MinPos)
iNewPart=0
PositionNbr = 0
nNewElectrons = 0

DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  MaxPos(1) = MAX(MaxPos(1),PartState(1,iPart))
  MaxPos(2) = MAX(MaxPos(2),PartState(2,iPart))
  MaxPos(3) = MAX(MaxPos(3),PartState(3,iPart))
  MinPos(1) = MIN(MinPos(1),PartState(1,iPart))
  MinPos(2) = MIN(MinPos(2),PartState(2,iPart))
  MinPos(3) = MIN(MinPos(3),PartState(3,iPart))
  IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
    nNewElectrons = nNewElectrons + 1
    IonIndX(nNewElectrons) = iPart
  END IF
END DO
ALLOCATE(iPartIndx_NodeTotalAmbi(nPart + nNewElectrons))
iPartIndx_NodeTotalAmbi(1:nPart) = iPartIndx_Node(1:nPart)
TotalPartNum = nPart

DO iLoop = 1, nNewElectrons
  DSMCSumOfFormedParticles = DSMCSumOfFormedParticles + 1
  PositionNbr = GetNextFreePosition()
  InsideFlag=.FALSE.
  iElem = PEM%GlobalElemID(iPartIndx_Node(1))
  DO WHILE(.NOT.InsideFlag)
    CALL RANDOM_NUMBER(Vec3D(1:3))
    RandomPos(1:3) = MinPos(1:3)+Vec3D(1:3)*(MaxPos(1:3)-MinPos(1:3))
    IF(Symmetry%Order.LE.2) RandomPos(3) = 0.
    IF(Symmetry%Order.LE.1) RandomPos(2) = 0.
    InsideFlag = ParticleInsideCheck(RandomPos,iPart,iElem)
  END DO
  PartState(1:3,PositionNbr) = RandomPos(1:3)
  PartSpecies(PositionNbr) = DSMC%AmbiDiffElecSpec
  PEM%GlobalElemID(PositionNbr) = iElem
  PDM%ParticleInside(PositionNbr) = .TRUE.
  PDM%isNewPart(PositionNbr) = .TRUE.
  PartState(4:6,PositionNbr) = AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(1:3)
  DEALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)
  IF ((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
    PartStateIntEn( 1,PositionNbr) = 0.
    PartStateIntEn( 2,PositionNbr) = 0.
    IF (DSMC%ElectronicModel.GT.0)   PartStateIntEn( 3,PositionNbr) = 0.
  END IF
  IF(usevMPF.OR.RadialWeighting%DoRadialWeighting) PartMPF(PositionNbr) = PartMPF(IonIndX(iLoop))
  IF(UseVarTimeStep) PartTimeStep(PositionNbr) = PartTimeStep(IonIndX(iLoop))
  iPartIndx_NodeTotalAmbi(nPart+iLoop) = PositionNbr
END DO
! Output variable
TotalPartNum = nPart + nNewElectrons
! Variable used for allocation
PartNum = TotalPartNum
! In case of the background gas, additional particles will be added
IF(BGGas%NumberOfSpecies.GT.0) PartNum = PartNum + TotalPartNum

newAmbiParts = 0
IF (ALLOCATED(iPartIndx_NodeNewAmbi)) DEALLOCATE(iPartIndx_NodeNewAmbi)
ALLOCATE(iPartIndx_NodeNewAmbi(PartNum))

END SUBROUTINE AD_InsertParticles


SUBROUTINE AD_DeleteParticles(iPartIndx_Node, nPart_opt)
!===================================================================================================================================
!> Deletes all electron created for the collision process and saves their velocity vector to the ion they are attached to
!> 1) Counting the ions/electrons within particles that already existed within the cell (but may have changed their species)
!> 2) Counting the ions/electrons within particles that were newly created during the time step
!> 3) Assinging each ion an electron, saving its velocity vector and deleting it
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PARTICLE_Vars       ,ONLY: PDM, PartSpecies, Species, PartState, PEM
USE MOD_DSMC_Vars           ,ONLY: AmbipolElecVelo, newAmbiParts, iPartIndx_NodeNewAmbi
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: iPartIndx_Node(:), nPart_opt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iPart, iLoop, nElectron, nIon, nPart
INTEGER, ALLOCATABLE        :: ElecIndx(:), IonIndX(:)
!===================================================================================================================================
nElectron =0; nIon = 0

IF(PRESENT(nPart_opt)) THEN
  ALLOCATE(ElecIndx(2*nPart_opt), IonIndX(2*nPart_opt))
  nPart = nPart_opt
ELSE
  ALLOCATE(ElecIndx(newAmbiParts), IonIndX(newAmbiParts))
  nPart = 0
END IF

! 1) Treating the particles that already existed within the cell (but may have changed their species)
DO iLoop = 1, nPart
  iPart = iPartIndx_Node(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    IF(PARTISELECTRON(iPart)) THEN
      nElectron = nElectron + 1
      ElecIndx(nElectron) = iPart
    END IF
    IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
      nIon = nIon + 1
      IonIndX(nIon) = iPart
    END IF
  END IF
END DO
! 2) Treating particles that were newly created during the time step
DO iLoop = 1, newAmbiParts
  iPart = iPartIndx_NodeNewAmbi(iLoop)
  IF (PDM%ParticleInside(iPart)) THEN
    IF(PARTISELECTRON(iPart)) THEN
      nElectron = nElectron + 1
      ElecIndx(nElectron) = iPart
    END IF
    IF(Species(PartSpecies(iPart))%ChargeIC.GT.0.0) THEN
      nIon = nIon + 1
      IonIndX(nIon) = iPart
    END IF
  END IF
END DO
! Sanity check: number of electrons and ions has to be identical
IF(nIon.NE.nElectron) THEN
  IPWRITE(*,*) 'Initial number of particles:', nPart, 'Number of new particles', newAmbiParts
  DO iLoop = 1, nPart
    iPart = iPartIndx_Node(iLoop)
    IPWRITE(*,*) 'Index, Element ID, Particle Inside?:', iPart, PDM%ParticleInside(iPart), PEM%GlobalElemID(iPart)
    IPWRITE(*,*) 'Species, Charge:', PartSpecies(iPart), Species(PartSpecies(iPart))%ChargeIC
  END DO
  IPWRITE(*,*) 'List of new particles:'
  DO iLoop = 1, newAmbiParts
    iPart = iPartIndx_NodeNewAmbi(iLoop)
    IPWRITE(*,*) 'Index, Element ID, Particle Inside?:', iPart, PDM%ParticleInside(iPart), PEM%GlobalElemID(iPart)
    IPWRITE(*,*) 'Species, Charge:', PartSpecies(iPart), Species(PartSpecies(iPart))%ChargeIC
  END DO
  CALL abort(__STAMP__,&
      'ERROR: Number of electrons and ions is not equal for ambipolar diffusion: ',IntInfoOpt=nIon-nElectron)
END IF

! 3) Assinging each ion an electron, saving its velocity vector and deleting it
DO iLoop = 1, nElectron
  IF (ALLOCATED(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)) DEALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo)
  ALLOCATE(AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(3))
  AmbipolElecVelo(IonIndX(iLoop))%ElecVelo(1:3) = PartState(4:6,ElecIndx(iLoop))
  PDM%ParticleInside(ElecIndx(iLoop)) = .FALSE.
END DO

DEALLOCATE(iPartIndx_NodeNewAmbi)

END SUBROUTINE AD_DeleteParticles

END MODULE MOD_DSMC_AmbipolarDiffusion
