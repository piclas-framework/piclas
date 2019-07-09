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
MODULE MOD_LD_internal_Temp

!===================================================================================================================================
! module for determination of cell quantities
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

PUBLIC :: CalcInternalTemp_LD_first, CalcInternalTemp_LD_second, CalcInternalTemp_LD_third
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcInternalTemp_LD_first(iElem)

USE MOD_LD_Vars
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : PEM, Species, nSpecies, PartSpecies
USE MOD_Particle_Mesh_Vars,     ONLY : GEO
USE MOD_DSMC_Vars,              ONLY : SpecDSMC, PartStateIntEn, DSMC
USE MOD_TimeDisc_Vars,          ONLY : dt

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                  !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER                :: iPart, nPart, iPartIndx, iSpec
  REAL                   :: NumDensOfPartSpec(nSpecies)
  REAL                   :: CollFreq, CollNum, RotRelaxProbLD, VibRelaxProbLD, iRan, AFactor, BFactor, TauMW, TauV
  REAL                   :: PreRotEnergy, PreVibEnergy, PostRotEnergy, PostVibEnergy, NumDensTot, DoFConstV
  INTEGER                       :: iQuant
  INTEGER, INTENT(IN)           :: iElem

!--------------------------------------------------------------------------------------------------!

  nPart = PEM%pNumber(iElem)
  NumDensOfPartSpec(1:nSpecies) = 0.0

  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    NumDensOfPartSpec(PartSpecies(iPartIndx)) = NumDensOfPartSpec(PartSpecies(iPartIndx)) + 1.0
    iPartIndx = PEM%pNext(iPartIndx)
  END DO
  NumDensOfPartSpec(1:nSpecies) = NumDensOfPartSpec(1:nSpecies) * Species(1:nSpecies)%MacroParticleFactor / GEO%Volume(iElem)
  NumDensTot = SUM(NumDensOfPartSpec)

  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    PreVibEnergy  = PartStateIntEn(iPartIndx,1)
    PreRotEnergy  = PartStateIntEn(iPartIndx,2)
    CollFreq = 0.0
    DO iSpec = 1, nSpecies
        CollFreq = CollFreq + 0.5 * (SpecDSMC(iSpec)%DrefVHS + SpecDSMC(PartSpecies(iPartIndx))%DrefVHS)**2.0 &
                 * NumDensOfPartSpec(iSpec) &
                 * ( 2.0 * 3.14159265359 * BoltzmannConst * SpecDSMC(PartSpecies(iPartIndx))%TrefVHS &
                 * (Species(iSpec)%MassIC + Species(PartSpecies(iPartIndx))%MassIC) &
                 / (Species(iSpec)%MassIC * Species(PartSpecies(iPartIndx))%MassIC) )**0.5 &
                 * (PartStateBulkValues(iPartIndx,4) / SpecDSMC(PartSpecies(iPartIndx))%TrefVHS)**(0.5 - SpecDSMC(iSpec)%omegaVHS)
    END DO
    CollNum = 18.1 / ( 1.0 + 0.5*3.14159265359**1.5*(91.5/ PartStateBulkValues(iPartIndx,4) )**0.5 &
            + (3.14159265359 + 0.25*3.14159265359**2.0)*91.5 / PartStateBulkValues(iPartIndx,4) )
    RotRelaxProbLD = 3.0 / 5.0 * (1.0 - EXP(-5.0/3.0 * CollFreq/CollNum * dt))
    CALL RANDOM_NUMBER(iRan)

    IF(RotRelaxProbLD.GT.iRan) THEN
      CALL RANDOM_NUMBER(iRan)
      PartStateIntEn(iPartIndx,2) = -BoltzmannConst * PartStateBulkValues(iPartIndx,4) * LOG(iRan)
    END IF

    AFactor = 0.00116 * (0.5 * 6.02214129E26 * Species(PartSpecies(iPartIndx))%MassIC)**0.5 &
            * SpecDSMC(PartSpecies(iPartIndx))%CharaTVib**(4.0/3.0)
    BFactor = -0.015 * AFactor * (0.5 * 6.02214129E26 * Species(PartSpecies(iPartIndx))%MassIC)**0.25 - 18.42
    TauMW = EXP(AFactor * PartStateBulkValues(iPartIndx,4)**(-1.0/3.0) + BFactor) &
          / (NumDensTot * BoltzmannConst * PartStateBulkValues(iPartIndx,4)  )!* 1E-5)
    TauV = TauMW + ( 3.14159265359 * Species(PartSpecies(iPartIndx))%MassIC &
         / (8.0*BoltzmannConst*PartStateBulkValues(iPartIndx,4)) )**0.5 &
         / (5.81E-21*NumDensTot)

    DoFConstV = 1.0 + 2./3.*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/PartStateBulkValues(iPartIndx,4) &
              / (EXP(SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/PartStateBulkValues(iPartIndx,4)) - 1.0)
    VibRelaxProbLD = 1/DoFConstV * (1.0 - EXP(-DoFConstV*dt/TauV))

    CALL RANDOM_NUMBER(iRan)
    IF(VibRelaxProbLD.GT.iRan) THEN
      CALL RANDOM_NUMBER(iRan)
      iQuant = INT(-LOG(iRan)*PartStateBulkValues(iPartIndx,4)/SpecDSMC(PartSpecies(iPartIndx))%CharaTVib)
      DO WHILE (iQuant.GE.SpecDSMC(PartSpecies(iPartIndx))%MaxVibQuant)
        CALL RANDOM_NUMBER(iRan)
        iQuant = INT(-LOG(iRan)*PartStateBulkValues(iPartIndx,4)/SpecDSMC(PartSpecies(iPartIndx))%CharaTVib)
      END DO
      PartStateIntEn(iPartIndx, 1) = (iQuant + DSMC%GammaQuant)*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib*BoltzmannConst
    END IF
    PostVibEnergy = PartStateIntEn(iPartIndx,1)
    PostRotEnergy = PartStateIntEn(iPartIndx,2)
    PartStateBulkValues(iPartIndx,4) = PartStateBulkValues(iPartIndx,4) &
                                      + (PreRotEnergy+PreVibEnergy-PostRotEnergy-PostVibEnergy) &
                                      / (1.5 * BoltzmannConst)
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

END SUBROUTINE CalcInternalTemp_LD_first

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcInternalTemp_LD_second(iElem)

USE MOD_LD_Vars
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : PEM, Species, PartSpecies, PDM
USE MOD_Particle_Mesh_Vars,     ONLY : GEO
USE MOD_DSMC_Vars,              ONLY : SpecDSMC, CollInf, PartStateIntEn, DSMC
USE MOD_TimeDisc_Vars,          ONLY : dt
USE MOD_part_tools,             ONLY : DiceUnitVector

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                  !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER               :: iPart, nPart, iPartIndx, iPType
  INTEGER               :: nPair, iPair, cPart1, cPart2, cSpec1, cSpec2, iCase
  REAL                   :: PreTransEnergy, PostTransEnergy,iRan
  REAL                   :: Velosq, Velo1, Velo2
  REAL                   :: RandVal(2)
  REAL, ALLOCATABLE    :: TempPartVelo(:,:) ! (Loc_iPart,3)  velox; veloy; veloz
  INTEGER, ALLOCATABLE :: iPartIndxArray(:)
  REAL                   :: ProbabilityFactor, SpecNum1, SpecNum2   ! see LD Paper for LD relaxation (P_max)
TYPE tLDPairData
  REAL              :: CRela2     ! squared relative velo of the particles in a pair
  REAL              :: Prob       ! colision probability
  INTEGER           :: iPart_p1   ! first particle of the pair
  INTEGER           :: iPart_p2   ! second particle of the pair
  INTEGER           :: PairType   ! type of pair (=iCase, CollInf%Coll_Case)
  REAL              :: Ec         ! Collision Energy
  REAL              :: TempPairVelo1(3)   ! Temporal velocity (x; y; z) of first particle
  REAL              :: TempPairVelo2(3)   ! Temporal velocity (x; y; z) of first particle
END TYPE tLDPairData

TYPE(tLDPairData), ALLOCATABLE    :: LD_Coll_pData(:)            ! LD Data of collision pairs into a cell (nPair)

  REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
  REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
  REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
  LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2   ! Check whether rot or vib relax is performed
  REAL                          :: Xi_rel, Xi, FakXi                ! Factors of DOF
  INTEGER                       :: iQuaMax, iQua                    ! Quantum Numbers
  REAL                          :: MaxColQua                        ! Max. Quantum Number

  INTEGER, INTENT(IN)           :: iElem
  REAL                          :: RanVec(3)

!--------------------------------------------------------------------------------------------------!

ALLOCATE(TempPartVelo(PDM%maxParticleNumber,3))

#if (PP_TimeDiscMethod==1001)
  IF((BulkValues(iElem)%CellType.EQ.3).OR.(BulkValues(iElem)%CellType.EQ.4)) THEN  ! --- LD Cell ?
#endif
  nPart = PEM%pNumber(iElem)

!--------------------------------------------------------------------------------------------------!
! Init some varibles for pairing, calcprob etc.
!--------------------------------------------------------------------------------------------------!

  ALLOCATE(iPartIndxArray(nPart))
  nPair                     = INT(nPart/2)
  ALLOCATE(LD_Coll_pData(nPair))
  CollInf%Coll_CaseNum      = 0
  CollInf%Coll_SpecPartNum  = 0
  PreTransEnergy  = 0.0
  PostTransEnergy = 0.0

!--------------------------------------------------------------------------------------------------!
! resample particl velocity from a Maxwellian distribution at CellBulkTemperature
!--------------------------------------------------------------------------------------------------!

  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart

      Velosq = 2
      DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal)
        Velo1 = 2.*RandVal(1) - 1.
        Velo2 = 2.*RandVal(2) - 1.
        Velosq = Velo1**2 + Velo2**2
      END DO
      TempPartVelo(iPartIndx,1) = Velo1*SQRT(-2*BoltzmannConst*BulkValues(iElem)%BulkTemperature/ &
        Species(PartSpecies(iPartIndx))%MassIC*LOG(Velosq)/Velosq)                                !x-Komponente
      TempPartVelo(iPartIndx,2) = Velo2*SQRT(-2*BoltzmannConst*BulkValues(iElem)%BulkTemperature/ &
        Species(PartSpecies(iPartIndx))%MassIC*LOG(Velosq)/Velosq)                                !y-Komponente
      Velosq = 2
      DO WHILE ((Velosq .GE. 1.) .OR. (Velosq .EQ. 0.))
        CALL RANDOM_NUMBER(RandVal)
        Velo1 = 2.*RandVal(1) - 1.
        Velo2 = 2.*RandVal(2) - 1.
        Velosq = Velo1**2 + Velo2**2
      END DO
      TempPartVelo(iPartIndx,3) = Velo1*SQRT(-2*BoltzmannConst*BulkValues(iElem)%BulkTemperature/ &
        Species(PartSpecies(iPartIndx))%MassIC*LOG(Velosq)/Velosq)                                !z-Komponente

      PreTransEnergy = PreTransEnergy + 0.5 * Species(PartSpecies(iPartIndx))%MassIC * &
                     ( TempPartVelo(iPartIndx,1)**2 + TempPartVelo(iPartIndx,2)**2 + TempPartVelo(iPartIndx,3)**2 )

    CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx)) = CollInf%Coll_SpecPartNum(PartSpecies(iPartIndx)) + 1
    iPartIndxArray(iPart) = iPartIndx
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

!--------------------------------------------------------------------------------------------------!
! END: resample particl velocity from a Maxwellian distribution at CellBulkTemperature
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Calculation of average pre-translation-energie
!--------------------------------------------------------------------------------------------------!

  PreTransEnergy = PreTransEnergy / nPart

!--------------------------------------------------------------------------------------------------!
! Pairing
!--------------------------------------------------------------------------------------------------!

  DO iPair = 1, nPair
    CALL RANDOM_NUMBER(iRan)
    cPart1 = INT( 1 + nPart * iRan )
    LD_Coll_pData(iPair)%iPart_p1 = iPartIndxArray(cPart1)
    iPartIndxArray(cPart1) = iPartIndxArray(nPart)
    nPart = nPart - 1
    cPart2 = INT( 1 + nPart * iRan )
    LD_Coll_pData(iPair)%iPart_p2 = iPartIndxArray(cPart2)
    iPartIndxArray(cPart2) = iPartIndxArray(nPart)
    nPart = nPart - 1

    cSpec1 = PartSpecies(LD_Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(LD_Coll_pData(iPair)%iPart_p2) !spec of particle 2

    iCase = CollInf%Coll_Case(cSpec1, cSpec2)
    CollInf%Coll_CaseNum(iCase) = CollInf%Coll_CaseNum(iCase) + 1 !sum of coll case (Sab)

    LD_Coll_pData(iPair)%TempPairVelo1(1) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,1)
    LD_Coll_pData(iPair)%TempPairVelo1(2) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,2)
    LD_Coll_pData(iPair)%TempPairVelo1(3) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,3)

    LD_Coll_pData(iPair)%TempPairVelo2(1) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,1)
    LD_Coll_pData(iPair)%TempPairVelo2(2) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,2)
    LD_Coll_pData(iPair)%TempPairVelo2(3) = TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,3)

    LD_Coll_pData(iPair)%CRela2 = (TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,1) - &
                                   TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,1))**2 &
                                + (TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,2) - &
                                   TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,2))**2 &
                                + (TempPartVelo(LD_Coll_pData(iPair)%iPart_p1,3) - &
                                   TempPartVelo(LD_Coll_pData(iPair)%iPart_p2,3))**2
    LD_Coll_pData(iPair)%PairType = iCase
  END DO

  IF (nPart .NE. 0) THEN   ! Don`t forget the last particle without a friend!!!
    PostTransEnergy = PostTransEnergy + 0.5 * Species(PartSpecies(iPartIndxArray(nPart)))%MassIC * &
                    ( TempPartVelo(iPartIndxArray(nPart),1)**2 &
                    + TempPartVelo(iPartIndxArray(nPart),2)**2 &
                    + TempPartVelo(iPartIndxArray(nPart),3)**2 )
  END IF

  nPart = PEM%pNumber(iElem)
!--------------------------------------------------------------------------------------------------!
! END: Pairing
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Calculaton of Collision probability
!--------------------------------------------------------------------------------------------------!
  DO iPair = 1, nPair
    iPType = SpecDSMC(PartSpecies(LD_Coll_pData(iPair)%iPart_p1))%InterID &
           + SpecDSMC(PartSpecies(LD_Coll_pData(iPair)%iPart_p2))%InterID !definition of col  type (atom-atom; atom-mol...)
    SpecNum1 = CollInf%Coll_SpecPartNum(PartSpecies(LD_Coll_pData(iPair)%iPart_p1)) !number of spec1 particels
    SpecNum2 = CollInf%Coll_SpecPartNum(PartSpecies(LD_Coll_pData(iPair)%iPart_p2)) !number of spec2 particels
    cSpec1 = PartSpecies(LD_Coll_pData(iPair)%iPart_p1) !spec of particle 1
    cSpec2 = PartSpecies(LD_Coll_pData(iPair)%iPart_p2) !spec of particle 2

    IF (iPType .EQ. 2) THEN
      PostTransEnergy = PostTransEnergy + 0.5 * Species(cSpec1)%MassIC * &
                      ( LD_Coll_pData(iPair)%TempPairVelo1(1)**2 &
                      + LD_Coll_pData(iPair)%TempPairVelo1(2)**2 &
                      + LD_Coll_pData(iPair)%TempPairVelo1(3)**2 ) &
                      + 0.5 * Species(cSpec2)%MassIC * &
                      ( LD_Coll_pData(iPair)%TempPairVelo2(1)**2 &
                      + LD_Coll_pData(iPair)%TempPairVelo2(2)**2 &
                      + LD_Coll_pData(iPair)%TempPairVelo2(3)**2 )
      CYCLE
    END IF

    IF (SpecDSMC(PartSpecies(LD_Coll_pData(iPair)%iPart_p1))%InterID .EQ. 1) THEN ! first collision partner is an atom
      ProbabilityFactor = 1.0 !SpecDSMC(cSpec2)%RotRelaxProb ! muss erweitert werden, wennRotRelxa nicht mehr konst.
    ELSE
      ProbabilityFactor = 1.0 !SpecDSMC(cSpec1)%RotRelaxProb ! muss erweitert werden, wennRotRelxa nicht mehr konst.
    END IF

    LD_Coll_pData(iPair)%Prob = SpecNum1*SpecNum2/(1 + CollInf%KronDelta(LD_Coll_pData(iPair)%PairType))  &
              * CollInf%Cab(LD_Coll_pData(iPair)%PairType)                                               & ! Cab species comb fac
              * Species(PartSpecies(LD_Coll_pData(iPair)%iPart_p1))%MacroParticleFactor                  &
                      ! weighting Fact, here only one MPF is used!!!
              / CollInf%Coll_CaseNum(LD_Coll_pData(iPair)%PairType)                                      & !sum of coll cases Sab
              * LD_Coll_pData(iPair)%CRela2 ** (0.5-SpecDSMC(cSpec1)%omegaVHS) &
                      ! relative velo to the power of (1 -2omega) !! only one omega is used!!
              * dt / GEO%Volume(iElem)

    LD_Coll_pData(iPair)%Prob = LD_Coll_pData(iPair)%Prob * ProbabilityFactor

!--------------------------------------------------------------------------------------------------!
! END: Calculaton of Collision probability
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Desicion if Collision and Rotation, Vibration Relaxation of particles is performed and do it
!--------------------------------------------------------------------------------------------------!
    IF(LD_Coll_pData(iPair)%Prob .GT. 1.0) print*,LD_Coll_pData(iPair)%Prob

!!    LD_Coll_pData(iPair)%Prob = 1.0

    CALL RANDOM_NUMBER(iRan)
    IF (LD_Coll_pData(iPair)%Prob .GT. iRan) THEN    ! Collision???

      DoRot1  = .FALSE.
      DoRot2  = .FALSE.
      DoVib1  = .FALSE.
      DoVib2  = .FALSE.

      Xi_rel = 2*(2 - SpecDSMC(cSpec1)%omegaVHS)
        ! DOF of relative motion in VHS model, only for one omega!!

      LD_Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(LD_Coll_pData(iPair)%PairType)*LD_Coll_pData(iPair)%CRela2

      Xi = Xi_rel !Xi are all DOF in the collision

      CALL RANDOM_NUMBER(iRan)
      IF(SpecDSMC(cSpec1)%InterID.EQ.2) THEN
        IF(SpecDSMC(cSpec1)%RotRelaxProb / ProbabilityFactor .GT.iRan) THEN     ! Rotation???
          DoRot1 = .TRUE.
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec + PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,2) ! adding rot energy
                                                                                                              ! to coll energy
          Xi = Xi + SpecDSMC(cSpec1)%Xi_Rot
          IF(SpecDSMC(cSpec1)%VibRelaxProb / ProbabilityFactor .GT.iRan) DoVib1 = .TRUE.     ! Vibration???
        END IF
      END IF

      CALL RANDOM_NUMBER(iRan)
      IF(SpecDSMC(cSpec2)%InterID.EQ.2) THEN
        IF(SpecDSMC(cSpec2)%RotRelaxProb / ProbabilityFactor .GT.iRan) THEN     ! Rotation???
          DoRot2 = .TRUE.
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec + PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,2) ! adding rot energy
                                                                                                              ! to coll energy
          Xi = Xi + SpecDSMC(cSpec2)%Xi_Rot
          IF(SpecDSMC(cSpec2)%VibRelaxProb / ProbabilityFactor .GT.iRan) DoVib2 = .TRUE.     ! Vibration???
        END IF
      END IF

      FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!

      IF(DoVib1) THEN
        LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec + PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,1) ! adding vib energy
                                                                                                            ! to coll energy
        MaxColQua = LD_Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(cSpec1)%CharaTVib)  &
                  - DSMC%GammaQuant
        iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(cSpec1)%MaxVibQuant)
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)
        CALL RANDOM_NUMBER(iRan)
        DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
          !laux diss page 31
          CALL RANDOM_NUMBER(iRan)
          iQua = INT(iRan * iQuaMax)
          CALL RANDOM_NUMBER(iRan)
        END DO
          PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(cSpec1)%CharaTVib
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec - PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,1)
      END IF

      IF(DoVib2) THEN
        LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec + PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,1) ! adding vib energy
                                                                                                            ! to coll energy
        MaxColQua = LD_Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(cSpec2)%CharaTVib)  &
                  - DSMC%GammaQuant
        iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(cSpec2)%MaxVibQuant)
        CALL RANDOM_NUMBER(iRan)
        iQua = INT(iRan * iQuaMax)
        CALL RANDOM_NUMBER(iRan)
        DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
         !laux diss page 31
         CALL RANDOM_NUMBER(iRan)
         iQua = INT(iRan * iQuaMax)
         CALL RANDOM_NUMBER(iRan)
        END DO
          PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(cSpec2)%CharaTVib
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec - PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,1)
      END IF

!--------------------------------------------------------------------------------------------------!
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------!

      IF(DoRot1) THEN
        CALL RANDOM_NUMBER(iRan)
          PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,2) = LD_Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec - PartStateIntEn(LD_Coll_pData(iPair)%iPart_p1,2)
          FakXi = FakXi - 0.5*SpecDSMC(cSpec2)%Xi_Rot
      END IF

      IF(DoRot2) THEN
        CALL RANDOM_NUMBER(iRan)
          PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,2) = LD_Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
          LD_Coll_pData(iPair)%Ec = LD_Coll_pData(iPair)%Ec - PartStateIntEn(LD_Coll_pData(iPair)%iPart_p2,2)
      END IF

!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

      FracMassCent1 = CollInf%FracMassCent(PartSpecies(LD_Coll_pData(iPair)%iPart_p1), LD_Coll_pData(iPair)%PairType)
      FracMassCent2 = CollInf%FracMassCent(PartSpecies(LD_Coll_pData(iPair)%iPart_p2), LD_Coll_pData(iPair)%PairType)

      !Calculation of velo from center of mass
      VeloMx = FracMassCent1 * LD_Coll_pData(iPair)%TempPairVelo1(1) &
             + FracMassCent2 * LD_Coll_pData(iPair)%TempPairVelo2(1)
      VeloMy = FracMassCent1 * LD_Coll_pData(iPair)%TempPairVelo1(2) &
             + FracMassCent2 * LD_Coll_pData(iPair)%TempPairVelo2(2)
      VeloMz = FracMassCent1 * LD_Coll_pData(iPair)%TempPairVelo1(3) &
             + FracMassCent2 * LD_Coll_pData(iPair)%TempPairVelo2(3)

      !calculate random vec and new squared velocities
      LD_Coll_pData(iPair)%CRela2 = 2 * LD_Coll_pData(iPair)%Ec/CollInf%MassRed(LD_Coll_pData(iPair)%PairType)
      RanVec(1:3) = DiceUnitVector()

      RanVelox = SQRT(LD_Coll_pData(iPair)%CRela2) * RanVec(1)
      RanVeloy = SQRT(LD_Coll_pData(iPair)%CRela2) * RanVec(2)
      RanVeloz = SQRT(LD_Coll_pData(iPair)%CRela2) * RanVec(3)

      ! deltaV particle 1
      LD_Coll_pData(iPair)%TempPairVelo1(1) = VeloMx + FracMassCent2*RanVelox
      LD_Coll_pData(iPair)%TempPairVelo1(2) = VeloMy + FracMassCent2*RanVeloy
      LD_Coll_pData(iPair)%TempPairVelo1(3) = VeloMz + FracMassCent2*RanVeloz
     ! deltaV particle 2
      LD_Coll_pData(iPair)%TempPairVelo2(1) = VeloMx - FracMassCent1*RanVelox
      LD_Coll_pData(iPair)%TempPairVelo2(2) = VeloMy - FracMassCent1*RanVeloy
      LD_Coll_pData(iPair)%TempPairVelo2(3) = VeloMz - FracMassCent1*RanVeloz

      PostTransEnergy = PostTransEnergy + 0.5 * Species(cSpec1)%MassIC * &
                     ( LD_Coll_pData(iPair)%TempPairVelo1(1)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo1(2)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo1(3)**2 ) &
                     + 0.5 * Species(cSpec2)%MassIC * &
                     ( LD_Coll_pData(iPair)%TempPairVelo2(1)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo2(2)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo2(3)**2 )
    ELSE    ! no Collision???
      PostTransEnergy = PostTransEnergy + 0.5 * Species(cSpec1)%MassIC * &
                     ( LD_Coll_pData(iPair)%TempPairVelo1(1)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo1(2)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo1(3)**2 ) &
                     + 0.5 * Species(cSpec2)%MassIC * &
                     ( LD_Coll_pData(iPair)%TempPairVelo2(1)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo2(2)**2 &
                     + LD_Coll_pData(iPair)%TempPairVelo2(3)**2 )
    END IF    ! end Collision???

!--------------------------------------------------------------------------------------------------!
! Desicion if Collision and Rotation, Vibration Relaxation of particles is performed and do it
!--------------------------------------------------------------------------------------------------!

  END DO  ! loop over all pairs

!--------------------------------------------------------------------------------------------------!
! Calculation of new CellBulkTemperature out of post- and pre-translation-energie
!--------------------------------------------------------------------------------------------------!

  PostTransEnergy = PostTransEnergy / nPart

  BulkValues(iElem)%BulkTemperature = BulkValues(iElem)%BulkTemperature + 2.0 / (3.0 * BoltzmannConst) &
                                    * (PostTransEnergy - PreTransEnergy)

#if (PP_TimeDiscMethod==1001)
  END IF  ! --- END LD Cell or Bufferzone_A?
#endif
  DEALLOCATE(iPartIndxArray)
  DEALLOCATE(LD_Coll_pData)

DEALLOCATE(TempPartVelo)

END SUBROUTINE CalcInternalTemp_LD_second

!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcInternalTemp_LD_third(iElem)

USE MOD_LD_Vars
USE MOD_Globals_Vars,           ONLY : BoltzmannConst
USE MOD_Particle_Vars,          ONLY : PEM, Species, nSpecies, PartSpecies
USE MOD_Particle_Mesh_Vars,     ONLY : GEO
USE MOD_DSMC_Vars,              ONLY : SpecDSMC, CollInf, PartStateIntEn, DSMC
USE MOD_TimeDisc_Vars,          ONLY : dt

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                  !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER               :: iPart, nPart, iPartIndx, iSpec, jSpec
  REAL                   :: VibTemp, Xi_vib
  REAL                   :: NewTransEnergy, NewRotEnergy, NewVibEnergy, NewVibTemp, EquXi_vib
  REAL                   :: SumEnergy, EquEnergy, MeanVibQua, NewMeanVibQua
  REAL                   :: VibRelaPart(nSpecies), RotRelaPart(nSpecies), NumOfPartSpec(nSpecies)
  REAL                   :: PreRotEnergyPerSpec(nSpecies), PreVibEnergyPerSpec(nSpecies), Xi_vibPerSpec(nSpecies)
  REAL                   :: PreTransEnergy, PreRotEnergy, PreVibEnergy
  REAL                   :: CollPartion, SumOfCombination, PartionFactor, MolNum, Xi_vibPerPart
  INTEGER                :: PairType
  REAL                   :: CRela2, GasMixFak

  INTEGER, INTENT(IN)           :: iElem

!--------------------------------------------------------------------------------------------------!

  nPart = PEM%pNumber(iElem)
  PreVibEnergyPerSpec(1:nSpecies) = 0.0
  PreRotEnergyPerSpec(1:nSpecies) = 0.0
  NumOfPartSpec(1:nSpecies) = 0.0
  PreVibEnergy = 0.0
  PreRotEnergy = 0.0
  Xi_vib = 0.0
  NewRotEnergy = 0.0
  NewVibEnergy = 0.0

  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    NumOfPartSpec(PartSpecies(iPartIndx)) = NumOfPartSpec(PartSpecies(iPartIndx)) + 1.0
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

  DO iSpec = 1, nSpecies
    CollPartion = 0.0
    SumOfCombination = 0.0
    DO jSpec = 1, nSpecies
      IF (jSpec .EQ. iSpec) THEN
        SumOfCombination = SumOfCombination + 0.5 * NumOfPartSpec(iSpec) * (NumOfPartSpec(iSpec) - 1.0)
      ELSE
        SumOfCombination = SumOfCombination + NumOfPartSpec(iSpec) * NumOfPartSpec(jSpec)
      END IF
    END DO

    DO jSpec = 1, nSpecies
      PairType = CollInf%Coll_Case(iSpec,jSpec)
      CRela2 = 8.0 * BoltzmannConst * BulkValues(iElem)%BulkTemperature &
             / (3.1415926536 * CollInf%MassRed(PairType))
      IF (jSpec .EQ. iSpec) THEN
        PartionFactor = 0.5 * NumOfPartSpec(iSpec) * (NumOfPartSpec(iSpec) - 1.0) / SumOfCombination
      ELSE
        PartionFactor = NumOfPartSpec(iSpec) * NumOfPartSpec(jSpec) / SumOfCombination
      END IF
      CollPartion = CollPartion + PartionFactor &
                  * (nPart - 1.0) &
                  * CollInf%Cab(PairType) & ! Cab species comb fac
                  * Species(iSpec)%MacroParticleFactor                  &
                          ! weighting Fact, here only one MPF is used!!!
                  * CRela2 ** (0.5-SpecDSMC(iSpec)%omegaVHS) &
                          ! relative velo to the power of (1 -2omega) !! only one omega is used!!
                  * dt / GEO%Volume(iElem)
    END DO
    RotRelaPart(iSpec) = SpecDSMC(iSpec)%RotRelaxProb * CollPartion
    VibRelaPart(iSpec) = SpecDSMC(iSpec)%VibRelaxProb * CollPartion
  END DO
  PreTransEnergy = 1.0/2.0 * BoltzmannConst * BulkValues(iElem)%BulkTemperature ! per DOF
  MolNum = 0.0
  iPartIndx = PEM%pStart(iElem)
  Xi_vibPerSpec(:) = 0.
  DO ipart = 1, nPart
    IF (SpecDSMC(PartSpecies(iPartIndx))%InterID .EQ. 2) THEN
      PreVibEnergyPerSpec(PartSpecies(iPartIndx)) = PreVibEnergyPerSpec(PartSpecies(iPartIndx)) + PartStateIntEn(iPartIndx,1)
      PreRotEnergyPerSpec(PartSpecies(iPartIndx)) = PreRotEnergyPerSpec(PartSpecies(iPartIndx)) + PartStateIntEn(iPartIndx,2)

      MeanVibQua = PartStateIntEn(iPartIndx,1)/(BoltzmannConst*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib)
      VibTemp = SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/LOG(1 + 1/(MeanVibQua-DSMC%GammaQuant))

      IF (MeanVibQua.LE.DSMC%GammaQuant) THEN
        Xi_vibPerPart = 0.0
      ELSE
        Xi_vibPerPart = 2.0 * SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/VibTemp &
               / (EXP(SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/VibTemp)-1) &
               + 2.0*DSMC%GammaQuant*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/VibTemp
      END IF

      PreVibEnergy = PreVibEnergy + PartStateIntEn(iPartIndx,1)
      PreRotEnergy = PreRotEnergy + PartStateIntEn(iPartIndx,2)
      Xi_vibPerSpec(PartSpecies(iPartIndx)) = Xi_vibPerSpec(PartSpecies(iPartIndx)) + Xi_vibPerPart
      Xi_vib = Xi_vib + Xi_vibPerPart
      MolNum = MolNum + 1.0
    END IF
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

  DO iSpec = 1, nSpecies
    Xi_vibPerSpec(iSpec) = Xi_vibPerSpec(iSpec) / NumOfPartSpec(iSpec)
    PreVibEnergyPerSpec(iSpec) = PreVibEnergyPerSpec(iSpec) / (Xi_vibPerSpec(iSpec) * NumOfPartSpec(iSpec))
    PreRotEnergyPerSpec(iSpec) = PreRotEnergyPerSpec(iSpec) / (2.0 * NumOfPartSpec(iSpec))
  END DO
  PreVibEnergy = PreVibEnergy / MolNum
  PreRotEnergy = PreRotEnergy / MolNum
  Xi_vib = Xi_vib / MolNum
  PreVibEnergy = PreVibEnergy / Xi_vib ! per DOF
  PreRotEnergy = PreRotEnergy / 2.0 ! per DOF

  SumEnergy  = 3.0*PreTransEnergy + 2.0*PreRotEnergy + Xi_vib*PreVibEnergy
  EquEnergy = ( 3.0*PreTransEnergy + 2.0*PreRotEnergy + Xi_vib*PreVibEnergy ) / (5.0 + Xi_vib)

  NewVibEnergy = 0.0
  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    IF (SpecDSMC(PartSpecies(iPartIndx))%InterID .EQ. 2) THEN

      NewMeanVibQua = EquEnergy*Xi_vib /(BoltzmannConst*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib)

      NewVibTemp = SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/LOG(1 + 1/(NewMeanVibQua-DSMC%GammaQuant))

      IF (NewMeanVibQua.LE.DSMC%GammaQuant) THEN
        EquXi_vib = 0.0
      ELSE
        EquXi_vib = 2.0 * SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/NewVibTemp &
                  / (EXP(SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/NewVibTemp)-1) &
                  + 2.0*DSMC%GammaQuant*SpecDSMC(PartSpecies(iPartIndx))%CharaTVib/NewVibTemp
      END IF
      PartStateIntEn(iPartIndx,1) = (1.0-VibRelaPart(PartSpecies(iPartIndx))) &
                                  * PreVibEnergyPerSpec(PartSpecies(iPartIndx)) * Xi_vibPerSpec(PartSpecies(iPartIndx)) &
                                  + VibRelaPart(PartSpecies(iPartIndx)) * EquEnergy * EquXi_vib
      NewVibEnergy = NewVibEnergy + PartStateIntEn(iPartIndx,1)
    END IF
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

  NewVibEnergy = NewVibEnergy / MolNum
  EquEnergy = (SumEnergy - NewVibEnergy) / 5.0

  iPartIndx = PEM%pStart(iElem)
  DO ipart = 1, nPart
    IF (SpecDSMC(PartSpecies(iPartIndx))%InterID .EQ. 2) THEN
      PartStateIntEn(iPartIndx,2) = (1.0-RotRelaPart(PartSpecies(iPartIndx))) &
                                  * PreRotEnergyPerSpec(PartSpecies(iPartIndx)) * 2.0 &
                                  + RotRelaPart(PartSpecies(iPartIndx)) * EquEnergy * 2.0
      NewRotEnergy = NewRotEnergy + PartStateIntEn(iPartIndx,2)
    END IF
    iPartIndx = PEM%pNext(iPartIndx)
  END DO

  NewRotEnergy = NewRotEnergy / MolNum
  NewTransEnergy = (SumEnergy - NewRotEnergy - NewVibEnergy) ! sum of DOF
  GasMixFak = MolNum / nPart
  BulkValues(iElem)%BulkTemperature = (1.0 - GasMixFak) * BulkValues(iElem)%BulkTemperature &
                                    + GasMixFak * 2.0/3.0 * NewTransEnergy / BoltzmannConst
END SUBROUTINE CalcInternalTemp_LD_third

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!   SUBROUTINE CalcInternalTemp_LD_third_back(iElem)
!
!   USE MOD_LD_Vars
!   USE MOD_Globals_Vars,           ONLY : BoltzmannConst
!   USE MOD_Particle_Vars,          ONLY : PEM, Species
!   USE MOD_Particle_Mesh_Vars,     ONLY : GEO
!   USE MOD_DSMC_Vars,              ONLY : SpecDSMC, CollInf, PartStateIntEn, DSMC
!   USE MOD_TimeDisc_Vars,          ONLY : dt
!
!   !--------------------------------------------------------------------------------------------------!
!      IMPLICIT NONE                                                                                  !
!   !--------------------------------------------------------------------------------------------------!
!   ! argument list declaration                                                                        !
!   ! Local variable declaration                                                                       !
!     INTEGER               :: iPart, nPart, iPartIndx, iSpec
!     REAL                   :: PreTransEnergy, PreRotEnergy, PreVibEnergy, VibTemp, Xi_vib
!     REAL                   :: NewTransEnergy, NewRotEnergy, NewVibEnergy, NewVibTemp, NewXi_vib, EquXi_vib
!     REAL                   :: SumEnergy, EquEnergy, MeanVibQua, NewMeanVibQua
!     REAL                   :: VibRelaPart, RotRelaPart, CollPartion
!     INTEGER                :: PairType
!     REAL                   :: CRela2
!
!     INTEGER, INTENT(IN)           :: iElem
!
!   !--------------------------------------------------------------------------------------------------!
!
!     nPart = PEM%pNumber(iElem)
!     PreVibEnergy = 0.0
!     PreRotEnergy = 0.0
!     iSpec = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ACHTUNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     PairType = CollInf%Coll_Case(1,1)
!     CRela2 = 8.0 * BoltzmannConst * BulkValues(iElem)%BulkTemperature &
!            / (3.1415926536 * CollInf%MassRed(PairType))
!
!     CollPartion = nPart*nPart/2  &
!                 * CollInf%Cab(PairType)   & ! Cab species comb fac
!                 * Species(iSpec)%MacroParticleFactor                  &
!                         ! weighting Fact, here only one MPF is used!!!
!                 / (nPart * 0.5) & !sum of coll cases Sab
!                 * CRela2 ** (0.5-SpecDSMC(iSpec)%omegaVHS) &
!                         ! relative velo to the power of (1 -2omega) !! only one omega is used!!
!                 * dt / GEO%Volume(iElem)
!
!     RotRelaPart = SpecDSMC(iSpec)%RotRelaxProb * CollPartion
!     VibRelaPart = SpecDSMC(iSpec)%VibRelaxProb * CollPartion
!
!     PreTransEnergy = 1.0/2.0 * BoltzmannConst * BulkValues(iElem)%BulkTemperature ! per DOF
!
!     iPartIndx = PEM%pStart(iElem)
!     DO ipart = 1, nPart
!       PreVibEnergy = PreVibEnergy + PartStateIntEn(iPartIndx,1)
!       PreRotEnergy = PreRotEnergy + PartStateIntEn(iPartIndx,2)
!       iPartIndx = PEM%pNext(iPartIndx)
!     END DO
!     PreVibEnergy = PreVibEnergy / REAL(nPart)
!     PreRotEnergy = PreRotEnergy / REAL(nPart)
!     MeanVibQua = PreVibEnergy /(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
!     VibTemp = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(MeanVibQua-DSMC%GammaQuant))
!
!     IF (MeanVibQua.LE.DSMC%GammaQuant) THEN
!       Xi_vib = 0.0
!     ELSE
!       Xi_vib = 2.0 * SpecDSMC(iSpec)%CharaTVib/VibTemp / (EXP(SpecDSMC(iSpec)%CharaTVib/VibTemp)-1) &
!              + 2.0*DSMC%GammaQuant*SpecDSMC(iSpec)%CharaTVib/VibTemp
!     END IF
!
!     PreVibEnergy = PreVibEnergy / Xi_vib ! per DOF
!     PreRotEnergy = PreRotEnergy / 2.0 ! per DOF
!
!
!     SumEnergy  = 3.0*PreTransEnergy + 2.0*PreRotEnergy + Xi_vib*PreVibEnergy
!     EquEnergy = ( 3.0*PreTransEnergy + 2.0*PreRotEnergy + Xi_vib*PreVibEnergy ) / (5.0 + Xi_vib)
!
!     NewMeanVibQua = EquEnergy*Xi_vib /(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
!
!     NewVibTemp = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(NewMeanVibQua-DSMC%GammaQuant))
!
!     IF (NewMeanVibQua.LE.DSMC%GammaQuant) THEN
!       EquXi_vib = 0.0
!     ELSE
!       EquXi_vib = 2.0 * SpecDSMC(iSpec)%CharaTVib/NewVibTemp / (EXP(SpecDSMC(iSpec)%CharaTVib/NewVibTemp)-1) &
!              + 2.0*DSMC%GammaQuant*SpecDSMC(iSpec)%CharaTVib/NewVibTemp
!     END IF
!     NewVibEnergy = (1.0-VibRelaPart) * PreVibEnergy * Xi_vib + VibRelaPart * EquEnergy * EquXi_vib
!
!     NewMeanVibQua = NewVibEnergy* Xi_vib /(BoltzmannConst*SpecDSMC(iSpec)%CharaTVib)
!     NewVibTemp = SpecDSMC(iSpec)%CharaTVib/LOG(1 + 1/(NewMeanVibQua-DSMC%GammaQuant))
!     IF (NewMeanVibQua.LE.DSMC%GammaQuant) THEN
!       NewXi_vib = 0.0
!     ELSE
!       NewXi_vib = 2.0 * SpecDSMC(iSpec)%CharaTVib/NewVibTemp / (EXP(SpecDSMC(iSpec)%CharaTVib/NewVibTemp)-1) &
!              + 2.0*DSMC%GammaQuant*SpecDSMC(iSpec)%CharaTVib/NewVibTemp
!     END IF
!
!     EquEnergy = (SumEnergy - NewVibEnergy) / 5.0
!
!     NewRotEnergy = (1.0-RotRelaPart) * PreRotEnergy * 2.0 + RotRelaPart * EquEnergy * 2.0
!
!     NewTransEnergy = (SumEnergy - NewRotEnergy - NewVibEnergy)  ! sum of DOF
!
!
!
!     BulkValues(iElem)%BulkTemperature = 2.0/3.0 * NewTransEnergy / BoltzmannConst
!     iPartIndx = PEM%pStart(iElem)
!     DO ipart = 1, nPart
!       PartStateIntEn(iPartIndx,1) = NewVibEnergy
!       PartStateIntEn(iPartIndx,2) = NewRotEnergy
!       iPartIndx = PEM%pNext(iPartIndx)
!     END DO
!
!   END SUBROUTINE CalcInternalTemp_LD_third_back
!   !-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_LD_internal_Temp
