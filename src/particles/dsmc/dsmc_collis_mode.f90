MODULE MOD_DSMC_Collis
!===================================================================================================================================
! module including collisions
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

PUBLIC :: DSMC_perform_collision
!===================================================================================================================================

CONTAINS

!--------------------------------------------------------------------------------------------------!
SUBROUTINE DSMC_Elastic_Col(iPair, iElem)
  
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, PairE_vMPF
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, RandomVec, NumRanVec, PartState, usevMPF, PartMPF, GEO
  USE MOD_vmpf_collision,         ONLY : vMPF_PostVelo

!--------------------------------------------------------------------------------------------------!
! perform collision
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
INTEGER                       :: iVec
REAL                          :: iRan
REAL                          :: Epre, Epost                      ! Energy check
! input variable declaration                                                                       !
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(INOUT)        :: iElem
!--------------------------------------------------------------------------------------------------! 

  FracMassCent1 = CollInf%FracMassCent(PartSpecies(Coll_pData(iPair)%iPart_p1), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(Coll_pData(iPair)%iPart_p2), Coll_pData(iPair)%PairType)
  
!Coll E check
!  Epre = 0.5*(Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC & 
!       * (PartState(Coll_pData(iPair)%iPart_p1, 4)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p1, 5)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p1, 6)**2) &
!       + Species(PartSpecies(Coll_pData(iPair)%iPart_p2))%MassIC & 
!       * (PartState(Coll_pData(iPair)%iPart_p2, 4)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p2, 5)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p2, 6)**2)) 


  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 4) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 4)
  VeloMy = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 5) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 5)
  VeloMz = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 6) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 6)

  IF(usevMPF) THEN
    IF (iPair.EQ.PairE_vMPF(1)) THEN  
      Coll_pData(iPair)%CRela2 = Coll_pData(iPair)%CRela2 + 2 * GEO%DeltaEvMPF(iElem) &
                               / CollInf%MassRed(Coll_pData(iPair)%PairType) &
                               / PartMPF(PairE_vMPF(2))
    END IF
  END IF

  !calculate random vec
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)
  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)
  
 ! deltaV particle 1 
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,1) = VeloMx + FracMassCent2*RanVelox &
          - PartState(Coll_pData(iPair)%iPart_p1, 4)
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,2) = VeloMy + FracMassCent2*RanVeloy &
          - PartState(Coll_pData(iPair)%iPart_p1, 5)
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,3) = VeloMz + FracMassCent2*RanVeloz &
          - PartState(Coll_pData(iPair)%iPart_p1, 6)
 ! deltaV particle 2
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,1) = VeloMx - FracMassCent1*RanVelox &
          - PartState(Coll_pData(iPair)%iPart_p2, 4)
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,2) = VeloMy - FracMassCent1*RanVeloy &
          - PartState(Coll_pData(iPair)%iPart_p2, 5)
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,3) = VeloMz - FracMassCent1*RanVeloz &
          - PartState(Coll_pData(iPair)%iPart_p2, 6)

  IF(usevMPF) CALL vMPF_PostVelo(iPair, iElem)

!Echeck
!  Epost = 0.5* (Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC &
!        * ((VeloMx + FracMassCent2*RanVelox)**2 & 
!        + (VeloMy + FracMassCent2*RanVeloy)**2 &
!        + (VeloMz + FracMassCent2*RanVeloz)**2) &
!        + Species(PartSpecies(Coll_pData(iPair)%iPart_p2))%MassIC & 
!        * ((VeloMx - FracMassCent1*RanVelox)**2 & 
!        + (VeloMy - FracMassCent1*RanVeloy)**2 &
!        + (VeloMz - FracMassCent1*RanVeloz)**2)) 

!print*, Epre, Epost

END SUBROUTINE DSMC_Elastic_Col

!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
SUBROUTINE DSMC_Relax_Col_LauxTSHO(iPair, iElem)
  
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, DSMC, &
                                         SpecDSMC, PartStateIntEn, PairE_vMPF
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, RandomVec, NumRanVec, &
                                         PartState, BoltzmannConst, usevMPF, GEO, PartMPF
  USE MOD_vmpf_collision,         ONLY : vMPF_PostVelo 
  USE MOD_DSMC_ElectronicModel,   ONLY : ElectronicEnergyExchange, TVEEnergyExchange
  USE MOD_DSMC_PolyAtomicModel,   ONLY : DSMC_RotRelaxPoly, DSMC_VibRelaxPoly

!--------------------------------------------------------------------------------------------------!
! perform collision
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
REAL                          :: FracMassCent1, FracMassCent2     ! mx/(mx+my)
REAL                          :: VeloMx, VeloMy, VeloMz           ! center of mass velo
REAL                          :: RanVelox, RanVeloy, RanVeloz     ! random relativ velo
INTEGER                       :: iVec
REAL (KIND=8)                 :: iRan, iRan2
REAL                          :: Epre, Epost                      ! Energy check
LOGICAL                       :: DoRot1, DoRot2, DoVib1, DoVib2   ! Check whether rot or vib relax is performed
REAL (KIND=8)                 :: Xi_rel, Xi, FakXi                ! Factors of DOF
INTEGER                       :: iQuaMax, iQua                    ! Quantum Numbers
REAL                          :: MaxColQua                        ! Max. Quantum Number
REAL                          :: PartStateIntEnTemp, Phi, DeltaPartStateIntEn ! temp. var for inertial energy (needed for vMPF)
! input variable declaration                                                                       !
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(INOUT)        :: iElem
!--------------------------------------------------------------------------------------------------!
! variables for electronic level relaxation and transition
LOGICAL                       :: DoElec1, DoElec2
INTEGER                       :: iQuaold, MaxElecQuant, iQuaMax2, ii, iQuaMax3
REAL                          :: gtemp, gmax
!--------------------------------------------------------------------------------------------------!

  DoRot1  = .FALSE.
  DoRot2  = .FALSE.
  DoVib1  = .FALSE.
  DoVib2  = .FALSE.
  DoElec1 = .FALSE.
  DoElec2 = .FALSE.
  
!    Epre = 0.5*(Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC  &
!       * (PartState(Coll_pData(iPair)%iPart_p1, 4)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p1, 5)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p1, 6)**2) &
!       + Species(PartSpecies(Coll_pData(iPair)%iPart_p2))%MassIC &
!       * (PartState(Coll_pData(iPair)%iPart_p2, 4)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p2, 5)**2 &
!       + PartState(Coll_pData(iPair)%iPart_p2, 6)**2)) & 


  Xi_rel = 2*(2 - SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%omegaVHS) 
    ! DOF of relative motion in VHS model, only for one omega!!
 
  Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2

!  IF((usevMPF).AND.(iPair.EQ.PairE_vMPF(1))) THEN         ! adding energy lost due to vMPF
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + GEO%DeltaEvMPF(iElem) / PartMPF(PairE_vMPF(2))
!    GEO%DeltaEvMPF(iElem) = 0.0
!  END IF

  Xi = Xi_rel !Xi are all DOF in the collision

!--------------------------------------------------------------------------------------------------!
! Desicion if Rotation, Vibration and Electronic Relaxation of particles is performed
!--------------------------------------------------------------------------------------------------!

  CALL RANDOM_NUMBER(iRan)
  IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.EQ.2) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%RotRelaxProb.GT.iRan) THEN
      DoRot1 = .TRUE.
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) ! adding rot energy to coll energy
      Xi = Xi + SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%Xi_Rot
      IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%VibRelaxProb.GT.iRan) DoVib1 = .TRUE.
    END IF
  END IF
  IF ( DSMC%ElectronicState ) THEN
    ! step as TRUE
    IF ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID .ne. 4 ) THEN
      IF ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%ElecRelaxProb .GT. iRan ) THEN
        DoElec1 = .TRUE.
      END IF
    END IF
  END IF

  CALL RANDOM_NUMBER(iRan)
  IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%InterID.EQ.2) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%RotRelaxProb.GT.iRan) THEN
      DoRot2 = .TRUE.
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2) ! adding rot energy to coll energy
      Xi = Xi + SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%Xi_Rot
      IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%VibRelaxProb.GT.iRan) DoVib2 = .TRUE.
    END IF
  END IF
  IF ( DSMC%ElectronicState ) THEN
    ! step as TRUE
    IF ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID .ne. 4 ) THEN
      IF ( SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%ElecRelaxProb .GT. iRan ) THEN
        DoElec2 = .TRUE.
      END IF
    END IF
  END IF

  FakXi = 0.5*Xi  - 1  ! exponent factor of DOF, substitute of Xi_c - Xi_vib, laux diss page 40

!--------------------------------------------------------------------------------------------------!
! Electronic Relaxation / Transition
!--------------------------------------------------------------------------------------------------!

  IF(DoElec1.AND.DoVib1) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,3)  &
                         +    PartStateIntEn(Coll_pData(iPair)%iPart_p1,1)
#if ( PP_TimeDiscMethod == 42)
    CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p1,FakXi,Coll_pData(iPair)%iPart_p2)
#else
    CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p1,FakXi)
#endif
    DoElec1=.false.
    DoVib1=.false.
  END IF

  IF(DoElec2.AND.DoVib2) THEN
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)  &
                         +    PartStateIntEn(Coll_pData(iPair)%iPart_p2,1)
#if ( PP_TimeDiscMethod == 42)
    CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi,Coll_pData(iPair)%iPart_p1)
#else
    CALL TVEEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi)
#endif
    DoElec2=.false.
    DoVib2=.false.
  END IF

!    IF ( DSMC%ElectronicState ) THEN
  !    IF (usevMPF) THEN
  !      STOP 'vMPF not implemented for electronic relaxation'
  !    END IF
      ! Relaxation of first particle
  IF ( DoElec1 ) THEN
      ! calculate energy for electronic relaxation of particle 1
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,3)
#if ( PP_TimeDiscMethod == 42)
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p1,FakXi,Coll_pData(iPair)%iPart_p2)
#else
    IF (usevMPF) THEN
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p1,FakXi,Coll_pData(iPair)%iPart_p2,iElem)
    ELSE
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p1,FakXi )
    END IF
#endif
  END IF
    ! Electronic relaxation of second particle
  IF ( DoElec2 ) THEN
    ! calculate energy for electronic relaxation of particle 1
    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,3)
#if ( PP_TimeDiscMethod == 42)
    CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi,Coll_pData(iPair)%iPart_p1)
#else
    IF (usevMPF) THEN
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi,Coll_pData(iPair)%iPart_p1,iElem)
    ELSE
      CALL ElectronicEnergyExchange(Coll_pData(iPair)%Ec,Coll_pData(iPair)%iPart_p2,FakXi )
    END IF
#endif
  END IF
!  END IF

#if (PP_TimeDiscMethod==42)
  ! for TimeDisc 42 & only transition counting: prohibit relaxation and energy exchange
  IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif


!--------------------------------------------------------------------------------------------------!
! Vibrational Relaxation
!--------------------------------------------------------------------------------------------------!
!  IF (DoRot1) THEN
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2)
!    Xi = Xi + SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%Xi_Rot
!  END IF
!  IF (DoRot2) THEN 
!    Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
!    Xi = Xi + SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%Xi_Rot
!  END IF
!  FakXi = 0.5*Xi  - 1

  IF(DoVib1) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%PolyatomicMol) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1)
      CALL DSMC_VibRelaxPoly(Coll_pData(iPair)%Ec,PartSpecies(Coll_pData(iPair)%iPart_p1),Coll_pData(iPair)%iPart_p1,FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,1)
    ELSE
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) ! adding vib energy to coll energy
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib)  &
                - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
      DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
       !laux diss page 31
       CALL RANDOM_NUMBER(iRan)
       iQua = INT(iRan * iQuaMax)
       CALL RANDOM_NUMBER(iRan)
      END DO

      IF (usevMPF) THEN
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
    !      DeltaPartStateIntEn = 0.0
          Phi = PartMPF(Coll_pData(iPair)%iPart_p2) / PartMPF(Coll_pData(iPair)%iPart_p1)
          PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
          PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + Phi * PartStateIntEnTemp
          ! searche for new vib quant
          iQua = INT(PartStateIntEnTemp/ &
                 (BoltzmannConst*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib) - DSMC%GammaQuant)
          CALL RANDOM_NUMBER(iRan)
          IF(iRan .LT. PartStateIntEnTemp/(BoltzmannConst &
                     * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib) - DSMC%GammaQuant &
                     - DBLE(iQua)) THEN
            iQua = iQua + 1
          END IF
          PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib
          DeltaPartStateIntEn = PartMPF(Coll_pData(iPair)%iPart_p1) &
                              * (PartStateIntEnTemp - PartStateIntEn(Coll_pData(iPair)%iPart_p1,1))
    !      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2)) 
          GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
        END IF
      ELSE
        PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%CharaTVib
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,1)
      END IF
    END IF
  END IF

  IF(DoVib2) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%PolyatomicMol) THEN
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1)
      CALL DSMC_VibRelaxPoly(Coll_pData(iPair)%Ec,PartSpecies(Coll_pData(iPair)%iPart_p2),Coll_pData(iPair)%iPart_p2,FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,1)
    ELSE
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) ! adding vib energy to coll energy
      MaxColQua = Coll_pData(iPair)%Ec/(BoltzmannConst*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib)  &
                - DSMC%GammaQuant
      iQuaMax = MIN(INT(MaxColQua) + 1, SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%MaxVibQuant)
      CALL RANDOM_NUMBER(iRan)
      iQua = INT(iRan * iQuaMax)
      CALL RANDOM_NUMBER(iRan)
      DO WHILE (iRan.GT.(1 - iQua/MaxColQua)**FakXi)
       !laux diss page 31
       CALL RANDOM_NUMBER(iRan)
       iQua = INT(iRan * iQuaMax)
       CALL RANDOM_NUMBER(iRan)
      END DO

      IF (usevMPF) THEN
        IF (PartMPF(Coll_pData(iPair)%iPart_p2).GT.PartMPF(Coll_pData(iPair)%iPart_p1)) THEN
    !      DeltaPartStateIntEn = 0.0
          Phi = PartMPF(Coll_pData(iPair)%iPart_p1) / PartMPF(Coll_pData(iPair)%iPart_p2)
          PartStateIntEnTemp = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
          PartStateIntEnTemp = (DBLE(1)-Phi) * PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) + Phi * PartStateIntEnTemp
          ! searche for new vib quant
          iQua = INT(PartStateIntEnTemp/ &
                 (BoltzmannConst*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib) - DSMC%GammaQuant)
          CALL RANDOM_NUMBER(iRan)
          IF(iRan .LT. PartStateIntEnTemp/(BoltzmannConst &
                     * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib) - DSMC%GammaQuant &
                     - DBLE(iQua)) THEN
            iQua = iQua + 1
          END IF
          PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                        * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib
          DeltaPartStateIntEn = PartMPF(Coll_pData(iPair)%iPart_p2) &
                              * (PartStateIntEnTemp - PartStateIntEn(Coll_pData(iPair)%iPart_p2,1))
    !      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + DeltaPartStateIntEn / PartMPF(PairE_vMPF(2)) 
          GEO%DeltaEvMPF(iElem) = GEO%DeltaEvMPF(iElem) + DeltaPartStateIntEn
        END IF
      ELSE
        PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) = (iQua + DSMC%GammaQuant) * BoltzmannConst &
                      * SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%CharaTVib
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) 
      END IF
    END IF
  END IF

!--------------------------------------------------------------------------------------------------! 
! Rotational Relaxation
!--------------------------------------------------------------------------------------------------! 
  IF(DoRot1) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%PolyatomicMol) THEN
      FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, Coll_pData(iPair)%iPart_p1, FakXi)      
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      IF (usevMPF) THEN
        IF (PartMPF(Coll_pData(iPair)%iPart_p1).GT.PartMPF(Coll_pData(iPair)%iPart_p2)) THEN
          Phi = PartMPF(Coll_pData(iPair)%iPart_p2) / PartMPF(Coll_pData(iPair)%iPart_p1)
          PartStateIntEnTemp = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
          FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%Xi_Rot
          PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) = (1-Phi) * PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) &
                                                       + Phi * PartStateIntEnTemp
        END IF
      ELSE
        PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p1,2)
        FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%Xi_Rot
      END IF
    END IF
  END IF

  IF(DoRot2) THEN
    IF(SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%PolyatomicMol) THEN
      FakXi = FakXi - 0.5*SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p2))%Xi_Rot
      CALL DSMC_RotRelaxPoly(Coll_pData(iPair)%Ec, Coll_pData(iPair)%iPart_p2, FakXi)
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
    ELSE
      CALL RANDOM_NUMBER(iRan)
      IF (usevMPF) THEN
        IF (PartMPF(Coll_pData(iPair)%iPart_p2).GT.PartMPF(Coll_pData(iPair)%iPart_p1)) THEN
          Phi = PartMPF(Coll_pData(iPair)%iPart_p1) / PartMPF(Coll_pData(iPair)%iPart_p2)
          PartStateIntEnTemp = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
          Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEnTemp
          PartStateIntEn(Coll_pData(iPair)%iPart_p2,2) = (1-Phi) * PartStateIntEn(Coll_pData(iPair)%iPart_p2,2) &
                                                       + Phi * PartStateIntEnTemp
        END IF
      ELSE
        PartStateIntEn(Coll_pData(iPair)%iPart_p2,2) = Coll_pData(iPair)%Ec * (1.0 - iRan**(1.0/FakXi))
        Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec - PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
      END IF
    END IF
  END IF


!--------------------------------------------------------------------------------------------------!
! Calculation of new particle velocities
!--------------------------------------------------------------------------------------------------!

  IF(usevMPF) THEN
    IF (iPair.EQ.PairE_vMPF(1)) THEN         ! adding energy lost due to vMPF
      Coll_pData(iPair)%Ec = Coll_pData(iPair)%Ec + GEO%DeltaEvMPF(iElem) / PartMPF(PairE_vMPF(2))
      GEO%DeltaEvMPF(iElem) = 0.0
    END IF
  END IF

  FracMassCent1 = CollInf%FracMassCent(PartSpecies(Coll_pData(iPair)%iPart_p1), Coll_pData(iPair)%PairType)
  FracMassCent2 = CollInf%FracMassCent(PartSpecies(Coll_pData(iPair)%iPart_p2), Coll_pData(iPair)%PairType)

  !Calculation of velo from center of mass
  VeloMx = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 4) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 4)
  VeloMy = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 5) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 5)
  VeloMz = FracMassCent1 * PartState(Coll_pData(iPair)%iPart_p1, 6) &
         + FracMassCent2 * PartState(Coll_pData(iPair)%iPart_p2, 6)

  !calculate random vec and new squared velocities
  Coll_pData(iPair)%CRela2 = 2 * Coll_pData(iPair)%Ec/CollInf%MassRed(Coll_pData(iPair)%PairType)
  CALL RANDOM_NUMBER(iRan)
  iVec = INT(NumRanVec * iRan + 1)

  RanVelox = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,1)
  RanVeloy = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,2)
  RanVeloz = SQRT(Coll_pData(iPair)%CRela2) * RandomVec(iVec,3)

  ! deltaV particle 1
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,1) = VeloMx + FracMassCent2*RanVelox &
          - PartState(Coll_pData(iPair)%iPart_p1, 4)
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,2) = VeloMy + FracMassCent2*RanVeloy &
          - PartState(Coll_pData(iPair)%iPart_p1, 5)
  DSMC_RHS(Coll_pData(iPair)%iPart_p1,3) = VeloMz + FracMassCent2*RanVeloz &
          - PartState(Coll_pData(iPair)%iPart_p1, 6)
 ! deltaV particle 2
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,1) = VeloMx - FracMassCent1*RanVelox &
          - PartState(Coll_pData(iPair)%iPart_p2, 4)
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,2) = VeloMy - FracMassCent1*RanVeloy &
          - PartState(Coll_pData(iPair)%iPart_p2, 5)
  DSMC_RHS(Coll_pData(iPair)%iPart_p2,3) = VeloMz - FracMassCent1*RanVeloz &
          - PartState(Coll_pData(iPair)%iPart_p2, 6)

  IF(usevMPF) CALL vMPF_PostVelo(iPair, iElem)

!  Epost = 0.5*(Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC &
!        * ((VeloMx + FracMassCent2*RanVelox)**2 & 
!        + (VeloMy + FracMassCent2*RanVeloy)**2 &
!        + (VeloMz + FracMassCent2*RanVeloz)**2) &
!        + Species(PartSpecies(Coll_pData(iPair)%iPart_p1))%MassIC &
!        * ((VeloMx - FracMassCent1*RanVelox)**2 & 
!        + (VeloMy - FracMassCent1*RanVeloy)**2 &
!        + (VeloMz - FracMassCent1*RanVeloz)**2)) &
#if (PP_TimeDiscMethod==42)
  ! for TimeDisc 42 & only transition counting: prohibit relaxation and energy exchange
  END IF
# endif

END SUBROUTINE DSMC_Relax_Col_LauxTSHO


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE DSMC_perform_collision(iPair, iElem, NodeVolume, NodePartNum)

  USE MOD_DSMC_Vars,          ONLY : CollisMode , DSMC

!--------------------------------------------------------------------------------------------------!
! perform collision
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
! input variable declaration                                                                       !
  INTEGER, INTENT(IN)           :: iPair
  INTEGER, INTENT(INOUT)        :: iElem
  REAL, INTENT(IN), OPTIONAL    :: NodeVolume
  INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
  LOGICAL                       :: RelaxToDo
!--------------------------------------------------------------------------------------------------!
  SELECT CASE(CollisMode)
    CASE(1) ! elastic collision
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
        CALL DSMC_Elastic_Col(iPair, iElem)
#if (PP_TimeDiscMethod==42)
      END IF
# endif
    CASE(2) ! collision with relaxation
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
        CALL DSMC_Relax_Col_LauxTSHO(iPair, iElem)
#if (PP_TimeDiscMethod==42)
      END IF
# endif
    CASE(3) ! chemical reactions
      RelaxToDo = .TRUE.
      IF (PRESENT(NodeVolume).AND.PRESENT(NodePartNum)) THEN
        CALL ReactionDecicson(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
      ELSE
        CALL ReactionDecicson(iPair, RelaxToDo, iElem)
      END IF
#if (PP_TimeDiscMethod==42)
      ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
      IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
        IF (RelaxToDo) CALL DSMC_Relax_Col_LauxTSHO(iPair, iElem)
#if (PP_TimeDiscMethod==42)
      END IF
# endif
    CASE DEFAULT
      PRINT*, 'ERROR in DSMC_collis: Wrong Reaction Mode ', CollisMode
      STOP 
  END SELECT


END SUBROUTINE DSMC_perform_collision


!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

SUBROUTINE ReactionDecicson(iPair, RelaxToDo, iElem, NodeVolume, NodePartNum)
  USE MOD_DSMC_Vars,              ONLY : Coll_pData, CollInf, DSMC_RHS, DSMC, &
                                         SpecDSMC, PartStateIntEn, ChemReac
  USE MOD_Particle_Vars,          ONLY : Species, PartSpecies, PartState , BoltzmannConst, PEM, GEO, usevMPF
  USE MOD_DSMC_ChemReact,         ONLY : ElecImpactIoni, MolecDissoc, MolecExch, AtomRecomb 
  USE MOD_Globals,                ONLY : Unit_stdOut
  USE MOD_Equation_Vars,          ONLY : Pi
  USE MOD_vmpf_collision,         ONLY : AtomRecomb_vMPF
  USE MOD_DSMC_QK_PROCEDURES,     ONLY : QK_dissociation, QK_recombination, QK_exchange, &
                                         QK_ImpactIonization
!--------------------------------------------------------------------------------------------------!
! choose reactions
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
INTEGER                       :: CaseOfReaction, iReac, SpecToExec, PartToExec, PartReac2, iPart_p3!
INTEGER                       :: PartToExecSec, PartReac2Sec, iReac2
INTEGER                       :: nPartNode, nPair, PairForRec, iQua1, iQua2, ii
REAL                          :: ReactionProb, JToEv, EZeroPoint, ReactionProb2, Volume
REAL (KIND=8)                 :: iRan
!-------------------
! input variable declaration                                                                       !
INTEGER, INTENT(IN)           :: iPair
INTEGER, INTENT(INOUT)        :: iElem
LOGICAL, INTENT(INOUT)        :: RelaxToDo
REAL, INTENT(IN), OPTIONAL    :: NodeVolume
INTEGER, INTENT(IN), OPTIONAL :: NodePartNum
!--------------------------------------------------------------------------------------------------!
  JToEv = 1.602176565E-19

  IF (ChemReac%NumOfReact.EQ.0) THEN
  ELSE
    CaseOfReaction = ChemReac%ReactCase(PartSpecies(Coll_pData(iPair)%iPart_p1),PartSpecies(Coll_pData(iPair)%iPart_p2))
  END IF
  SELECT CASE(CaseOfReaction)
!--------------------------------------------------------------------------------------------------!
    CASE(1)! recombination of two atoms
!--------------------------------------------------------------------------------------------------!
      ! searching third collison partner
      IF(ChemReac%RecombParticle.EQ. 0) THEN
        IF(ChemReac%nPairForRec.GT. 1) THEN
          PairForRec = iPair
          DO WHILE(PairForRec.EQ.iPair)
            CALL RANDOM_NUMBER(iRan)
            PairForRec = 1 + INT(iRan * ChemReac%nPairForRec)
          END DO
          Coll_pData(PairForRec)%NeedForRec = .TRUE.
          iPart_p3 = Coll_pData(PairForRec)%iPart_p1
          ChemReac%RecombParticle = Coll_pData(PairForRec)%iPart_p2
          ChemReac%nPairForRec = ChemReac%nPairForRec - 1
        ELSE
          iPart_p3 = 0
        END IF
      ELSE
        iPart_p3 = ChemReac%RecombParticle
        ChemReac%RecombParticle = 0
      END IF
!--------------------------------------------------------------------------------------------------!
      IF ( iPart_p3 .GT. 0 ) THEN
        iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), &
                                  PartSpecies(Coll_pData(iPair)%iPart_p2), &
                                  PartSpecies(iPart_p3))
        IF ( ChemReac%QKProcedure(iReac)  ) THEN
          ! QK - model
          CALL QK_recombination(iPair,iReac,iPart_p3,RelaxToDo,iElem,NodeVolume,NodePartNum)
        ELSE
!--------------------------------------------------------------------------------------------------!
        ! traditional Recombination
          IF (PRESENT(NodeVolume)) THEN
            Volume = NodeVolume
          ELSE
            Volume = GEO%Volume(PEM%Element(iPart_p3))
          END IF
          IF (PRESENT(NodePartNum)) THEN
            nPartNode = NodePartNum
          ELSE
            nPartNode = PEM%pNumber(PEM%Element(iPart_p3))
          END IF
          Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2         &
                               + 0.5 * Species(PartSpecies(iPart_p3))%MassIC                                      * &
                               ( PartState(iPart_p3,4)**2 + PartState(iPart_p3,5)**2 + PartState(iPart_p3,6)**2 )   &
                              + PartStateIntEn(iPart_p3,1) + PartStateIntEn(iPart_p3,2)
!
!           iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2) &
!                 , PartSpecies(iPart_p3))
!
          EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
          IF((Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)).GE.EZeroPoint) THEN
            ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Rec_Arrhenius(ChemReac%MeanEVibQua_PerIter(PartSpecies(iPart_p3))) &
                         * nPartNode * Species(PartSpecies(iPart_p3))%MacroParticleFactor                                    &
                         / Volume * Coll_pData(iPair)%Ec**(ChemReac%Arrhenius_Powerfactor(iReac)                             &
                          - 0.5 + SpecDSMC(PartSpecies(iPart_p3))%omegaVHS)
          ELSE
            ReactionProb = 0.0
          END IF
#if (PP_TimeDiscMethod==42)
          IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
          END IF
# endif
          CALL RANDOM_NUMBER(iRan)
          IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
            IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
              IF (usevMPF) THEN
                CALL AtomRecomb_vMPF(iReac, iPair, iPart_p3, iElem)
              ELSE
                CALL AtomRecomb(iReac, iPair, iPart_p3)
              END IF
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
            END IF
# endif
            RelaxToDo = .FALSE.
          END IF
        END IF
      END IF
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
    CASE(2)! only dissociation
      iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
      IF ( ChemReac%QKProcedure(iReac)  ) THEN
        CALL QK_dissociation(iPair,iReac,RelaxToDo)
!--------------------------------------------------------------------------------------------------!
      ELSE
      ! reaction propability based on equilibrium method
      ! Reaction based on probability and rate coefficient
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                 &
                            + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                            + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
!         iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
          ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius(                                                         &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))                                               &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))                                               &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac))                                                     &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))            &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)                                            &
                          * Coll_pData(iPair)%Ec                                                                                &
                          ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor                                  &
                          - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))            &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)                                            &
                          * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
        ELSE
          ReactionProb = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
        END IF
# endif
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
          IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
            CALL MolecDissoc(iReac, iPair)
#if (PP_TimeDiscMethod==42)
          END IF
          IF ( DSMC%ReservoirRateStatistic  ) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
          END IF
# endif
          RelaxToDo = .FALSE.
        END IF
      END IF
!--------------------------------------------------------------------------------------------------!
    CASE(3)! only exchange reactions
      iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
      IF ( ChemReac%QKProcedure(iReac) ) THEN
        CALL QK_exchange(iPair,iReac,RelaxToDo)
!--------------------------------------------------------------------------------------------------!
      ELSE
    ! equilibrium and continuum based exchange reaction
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                  &
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
!         iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
! would be the equation if phi3 would not be zero for exchange reactions
!        ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
!                          ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
!                        * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
!                        ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
!                        + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * Coll_pData(iPair)%Ec &
!                        ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
!                        - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
          ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(                                                         &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)))                                              &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac))                                                     &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))            &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)                                            &
                          * Coll_pData(iPair)%Ec                                                                                &
                          ** (1.0 - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))    &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)
        ELSE
          ReactionProb = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
        END IF
# endif
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
          IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
            CALL MolecExch(iReac, iPair)
#if (PP_TimeDiscMethod==42)
          END IF
          IF ( DSMC%ReservoirRateStatistic  ) THEN
            ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
          END IF
# endif
          RelaxToDo = .FALSE.
        END IF
      END IF
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
    CASE(4)! exchange reaction and diss reaction possible
      iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
      IF ( ChemReac%QKProcedure(iReac)  ) THEN
        ! first check, if the the molecule dissociate, afterwards, check if an exchange reaction is possible
        CALL QK_dissociation(iPair,iReac,RelaxToDo)

        IF ( RelaxToDo  ) THEN
        ! exchange reactions
          iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
          IF ( ChemReac%QKProcedure(iReac)  ) THEN
            CALL QK_exchange(iPair,iReac,RelaxToDo)
          ELSE
            ! Arrhenius based Exchange Reaction
            Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                  &
                                 + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                                 + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
!         iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
            IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
              PartToExec = Coll_pData(iPair)%iPart_p1
              PartReac2 = Coll_pData(iPair)%iPart_p2
            ELSE
              PartToExec = Coll_pData(iPair)%iPart_p2
              PartReac2 = Coll_pData(iPair)%iPart_p1
            END IF
            EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
            IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
! would be the equation if phi3 would not be zero for exchange reactions
!        ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
!                          ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
!                        * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
!                        ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
!                        + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * Coll_pData(iPair)%Ec &
!                        ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
!                        - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
              ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius(                                                     &
                             ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)))                                             &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac))                                                     &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))            &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)                                            &
                          * Coll_pData(iPair)%Ec                                                                                &
                          ** (1.0 - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))    &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)
            ELSE
              ReactionProb = 0.0
            END IF
#if (PP_TimeDiscMethod==42)
            IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
            END IF
# endif
            CALL RANDOM_NUMBER(iRan)
            IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
              IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
                CALL MolecExch(iReac, iPair)
#if (PP_TimeDiscMethod==42)
              END IF
              IF ( DSMC%ReservoirRateStatistic  ) THEN
                ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
              END IF
# endif
              RelaxToDo = .FALSE.
            END IF
          END IF
        END IF
      ELSE
!--------------------------------------------------------------------------------------------------!
        ! purley Arrehnius rate based reactions
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2                  &
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
        ! calculation of dissociation probability
        iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF   
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
          ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * Coll_pData(iPair)%Ec &
                          ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                          - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
        ELSE
          ReactionProb = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
        END IF
# endif
        ! calculation of exchange reaction probability
        iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(ChemReac%DefinedReact(iReac,2,1))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
! would be the equation if phi3 would not be zero for exchange reactions
!        ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
!                          ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
!                        * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
!                        ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
!                        + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * Coll_pData(iPair)%Ec &
!                        ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
!                        - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
!                        , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
!                        * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
          ReactionProb2 = ChemReac%ReactInfo(iReac)%Beta_Exch_Arrhenius( &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec))) &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * Coll_pData(iPair)%Ec &
                          ** (1.0 - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2)
        ELSE
          ReactionProb2 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb2  ! for calculation of reactionrate coeficient
        END IF
# endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
            IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
              PartToExec = Coll_pData(iPair)%iPart_p1
              PartReac2 = Coll_pData(iPair)%iPart_p2
            ELSE
              PartToExec = Coll_pData(iPair)%iPart_p2
              PartReac2 = Coll_pData(iPair)%iPart_p1
            END IF
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
#if (PP_TimeDiscMethod==42)
            IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
              CALL MolecDissoc(iReac, iPair)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
            END IF
# endif
          ELSE
            iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
            IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
              PartToExec = Coll_pData(iPair)%iPart_p1
              PartReac2 = Coll_pData(iPair)%iPart_p2
            ELSE
              PartToExec = Coll_pData(iPair)%iPart_p2
              PartReac2 = Coll_pData(iPair)%iPart_p1
            END IF
#if (PP_TimeDiscMethod==42)
          ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
            IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
              CALL MolecExch(iReac, iPair)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
!--------------------------------------------------------------------------------------------------!
    CASE(5)! diss reaction and diss reaction possible
      iReac  = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
      iReac2 = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
      IF ( ChemReac%QKProcedure(iReac) .AND. ChemReac%QKProcedure(iReac2) ) THEN ! both Reaction QK
        ! collision energy without internal energy
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2
        ! first pseudo reaction probability
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        ReactionProb = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1)    - &
                         SpecDSMC(PartSpecies(PartToExec))%Ediss_eV * JToEv  )  / &
                       ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) )
        IF ( 0 .gt. ReactionProb ) THEN
          ReactionProb = 0
        END IF
        IF (ChemReac%DefinedReact(iReac2,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExecSec = Coll_pData(iPair)%iPart_p1
          PartReac2Sec = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExecSec = Coll_pData(iPair)%iPart_p2
          PartReac2Sec = Coll_pData(iPair)%iPart_p1
        END IF
        ! pseudo probability for second reaction
        ReactionProb2 = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1 )  - &
                          SpecDSMC(PartSpecies(PartToExecSec))%Ediss_eV * JToEv )  / &
                        ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1) )
        IF ( 0 .gt. ReactionProb2 ) THEN
          ReactionProb2 = 0
        END IF
        IF ( ReactionProb .gt. 0 ) THEN
          ! determine if first molecule dissociate
          CALL RANDOM_NUMBER(iRan)
          IF ( ReactionProb / ( ReactionProb + ReactionProb2) .gt. iRan) THEN
            CALL QK_dissociation(iPair,iReac,RelaxToDo)
            ! first molecule does not dissociate. what is the second particle doing?
          ELSE IF ( ReactionProb2 .gt. 0 ) THEN
          ! dissociation second molecule
            CALL QK_dissociation(iPair,iReac2,RelaxToDo)
          END IF
        ! ReactionProb = 0, check if second dissociation is possible
        ELSE IF ( ReactionProb2 .gt. 0 ) THEN
!           ! dissociationof second molecule
            CALL QK_dissociation(iPair,iReac2,RelaxToDo)
        END IF
      ELSE IF ( ChemReac%QKProcedure(iReac) .OR. ChemReac%QKProcedure(iReac2) ) THEN ! only one Reaction QK
        ! collision energy without internal energy
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2
        ! set QK reaction as first tested reaction
        IF (ChemReac%QKProcedure(iReac2)) THEN
          iReac  = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
          iReac2 = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
        END IF
        ! first pseude reaction probability
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        ReactionProb = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1)    - &
                         SpecDSMC(PartSpecies(PartToExec))%Ediss_eV * JToEv  )  / &
                       ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExec,1) )
        IF ( 0 .gt. ReactionProb ) THEN
          ReactionProb = 0
        END IF
        ! second reaction probability
        IF (ChemReac%DefinedReact(iReac2,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExecSec = Coll_pData(iPair)%iPart_p1
          PartReac2Sec = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExecSec = Coll_pData(iPair)%iPart_p2
          PartReac2Sec = Coll_pData(iPair)%iPart_p1
        END IF
        ! pseudo probability for second reaction
        ReactionProb2 = ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1 )  - &
                          SpecDSMC(PartSpecies(PartToExecSec))%Ediss_eV * JToEv )  / &
                        ( Coll_pData(iPair)%Ec + PartStateIntEn(PartToExecSec,1) )
        IF ( 0 .gt. ReactionProb2 ) THEN
          ReactionProb2 = 0
        END IF
        IF ( ReactionProb .gt. 0 ) THEN
          ! determine if first molecule dissociate
          CALL RANDOM_NUMBER(iRan)
          IF ( ReactionProb / ( ReactionProb + ReactionProb2) .gt. iRan) THEN
              CALL QK_dissociation(iPair,iReac,RelaxToDo)
          ELSE
          ! performe Arrehnius dissociation
            Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 &
                                  + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                                  + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)

            EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2Sec))%CharaTVib
            IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac2)) THEN
              ReactionProb = ChemReac%ReactInfo(iReac2)%Beta_Diss_Arrhenius(                                                     &
                             ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))                                            &
                             , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))                                          &
                             * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac2))                                                  &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac2) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac2,1,1))%omegaVHS&
                             + ChemReac%ReactInfo(iReac2)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))      &
                            , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))/2)                                        &
                             * Coll_pData(iPair)%Ec                                                                              &
                             ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac2,1,1))%VFD_Phi3_Factor                               &
                                 - ChemReac%ReactInfo(iReac2)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExecSec))  &
                                  , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2Sec)))/2)                                  &
                                 * PartStateIntEn(PartToExecSec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac2,1,1))%VFD_Phi3_Factor
            ELSE
              ReactionProb = 0.0
            END IF
#if (PP_TimeDiscMethod==42)
            IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
              ChemReac%NumReac(iReac2) = ChemReac%NumReac(iReac2) + ReactionProb  ! for calculation of reactionrate coeficient
            END IF
# endif
            CALL RANDOM_NUMBER(iRan)
            IF (ReactionProb.GT.iRan) THEN
#if (PP_TimeDiscMethod==42)
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
              IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
                CALL MolecDissoc(iReac2, iPair)
#if (PP_TimeDiscMethod==42)
              END IF
              IF ( DSMC%ReservoirRateStatistic  ) THEN
                ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
              END IF
# endif
              RelaxToDo = .FALSE.
            END IF
          END IF
        END IF
!--------------------------------------------------------------------------------------------------!
      ELSE ! both reactions Arrhenius
!--------------------------------------------------------------------------------------------------!
        ! Arrehnius rate based reations
        Coll_pData(iPair)%Ec = 0.5 * CollInf%MassRed(Coll_pData(iPair)%PairType)*Coll_pData(iPair)%CRela2 & 
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,1) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,1) &
                             + PartStateIntEn(Coll_pData(iPair)%iPart_p1,2) + PartStateIntEn(Coll_pData(iPair)%iPart_p2,2)
        ! calculation of dissociation probability
        iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
          ReactionProb = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * Coll_pData(iPair)%Ec &
                          ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                          - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
        ELSE
          ReactionProb = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb  ! for calculation of reactionrate coeficient
        END IF
#endif
        ! calculation of dissociation probability
        iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
        IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
          PartToExec = Coll_pData(iPair)%iPart_p1
          PartReac2 = Coll_pData(iPair)%iPart_p2
        ELSE
          PartToExec = Coll_pData(iPair)%iPart_p2
          PartReac2 = Coll_pData(iPair)%iPart_p1
        END IF
        EZeroPoint = DSMC%GammaQuant*BoltzmannConst*SpecDSMC(PartSpecies(PartReac2))%CharaTVib
        IF((Coll_pData(iPair)%Ec-EZeroPoint).GE.ChemReac%EActiv(iReac)) THEN
          ReactionProb2 = ChemReac%ReactInfo(iReac)%Beta_Diss_Arrhenius( &
                            ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2))) &
                          * (Coll_pData(iPair)%Ec - ChemReac%EActiv(iReac)) &
                          ** (ChemReac%Arrhenius_Powerfactor(iReac) - 1.5 + SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%omegaVHS &
                          + ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * Coll_pData(iPair)%Ec &
                          ** (1.0 - SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor &
                          - ChemReac%ReactInfo(iReac)%Xi_Total(ChemReac%MeanEVibQua_PerIter(PartSpecies(PartToExec)) &
                          , ChemReac%MeanEVibQua_PerIter(PartSpecies(PartReac2)))/2) &
                          * PartStateIntEn(PartToExec,1) ** SpecDSMC(ChemReac%DefinedReact(iReac,1,1))%VFD_Phi3_Factor
        ELSE
          ReactionProb2 = 0.0
        END IF
#if (PP_TimeDiscMethod==42)
        IF ( DSMC%ReservoirRateStatistic .EQV. .FALSE. ) THEN
          ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + ReactionProb2  ! for calculation of reactionrate coeficient
        END IF
#endif
        CALL RANDOM_NUMBER(iRan)
        IF ((ReactionProb + ReactionProb2).GT.iRan) THEN
          CALL RANDOM_NUMBER(iRan)
          IF((ReactionProb/(ReactionProb + ReactionProb2)).GT.iRan) THEN
            iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
            IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
              PartToExec = Coll_pData(iPair)%iPart_p1
              PartReac2 = Coll_pData(iPair)%iPart_p2
            ELSE
              PartToExec = Coll_pData(iPair)%iPart_p2
              PartReac2 = Coll_pData(iPair)%iPart_p1
            END IF
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
#if (PP_TimeDiscMethod==42)
            IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
              CALL MolecDissoc(iReac, iPair)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
            END IF
# endif
          ELSE
            iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 2)
            IF (ChemReac%DefinedReact(iReac,1,1).EQ.PartSpecies(Coll_pData(iPair)%iPart_p1)) THEN
              PartToExec = Coll_pData(iPair)%iPart_p1
              PartReac2 = Coll_pData(iPair)%iPart_p2
            ELSE
              PartToExec = Coll_pData(iPair)%iPart_p2
              PartReac2 = Coll_pData(iPair)%iPart_p1
            END IF
            ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
#if (PP_TimeDiscMethod==42)
            IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
              CALL MolecDissoc(iReac, iPair)
#if (PP_TimeDiscMethod==42)
            END IF
            IF ( DSMC%ReservoirRateStatistic  ) THEN
              ChemReac%NumReac(iReac) = ChemReac%NumReac(iReac) + 1  ! for calculation of reactionrate coeficient
            END IF
# endif
          END IF
          RelaxToDo = .FALSE.
        END IF
      END IF
    CASE(6) !only ionization
      iReac = ChemReac%ReactNum(PartSpecies(Coll_pData(iPair)%iPart_p1), PartSpecies(Coll_pData(iPair)%iPart_p2), 1)
      IF ( ChemReac%QKProcedure(iReac)  ) THEN
        IF ( .NOT. DSMC%ElectronicState ) THEN
            WRITE(*,*) 'ERROR! Atomic electron shell has to be initalized.'
            STOP
        END IF
        CALL QK_ImpactIonization(iPair,iReac,RelaxToDo)
      END IF
!---------------------------------------------------------------------------------------------------------------------------
      IF (.NOT. ChemReac%QKProcedure(iReac) ) THEN
      !.... Check if the total collision energy is greater than the ionization
      !     energy as a requirement for a possible reaction
        IF (SpecDSMC(PartSpecies(Coll_pData(iPair)%iPart_p1))%InterID.EQ.4) THEN
          SpecToExec = PartSpecies(Coll_pData(iPair)%iPart_p2)
        ELSE
          SpecToExec = PartSpecies(Coll_pData(iPair)%iPart_p1)
        END IF
        IF ((Coll_pData(iPair)%Ec).GE.SpecDSMC(SpecToExec)%Eion_eV*JToEv) THEN 
          ReactionProb = Coll_pData(iPair)%Sigma(2)/Coll_pData(iPair)%Sigma(0)
        ELSE
          ReactionProb = 0.0
        END IF
        CALL RANDOM_NUMBER(iRan)
        IF (ReactionProb.GT.iRan) THEN
        ! Reservoir simulation for obtaining the reaction rate at one given point does not require to performe the reaction
#if (PP_TimeDiscMethod==42)
          IF ( DSMC%ReservoirSimuRate .EQV. .FALSE. ) THEN
# endif
            CALL ElecImpactIoni(iReac, iPair)
#if (PP_TimeDiscMethod==42)
          END IF
# endif
        RelaxToDo = .FALSE.
      END IF
    END IF
  END SELECT


END SUBROUTINE ReactionDecicson

!--------------------------------------------------------------------------------------------------!
END MODULE MOD_DSMC_Collis
