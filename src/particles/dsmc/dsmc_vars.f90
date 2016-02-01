MODULE MOD_DSMC_Vars
!===================================================================================================================================
! Contains the DSMC variables
!===================================================================================================================================
! MODULES
#ifdef MPI
USE MOD_Particle_MPI_Vars, ONLY: tPartMPIConnect
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL                          :: Debug_Energy(2)=0.0        ! debug variable energy conservation
INTEGER                       :: DSMCSumOfFormedParticles   !number of formed particles per iteration in chemical reactions
                                                            ! for counting the nextfreeparticleposition

REAL  , ALLOCATABLE           :: DSMC_RHS(:,:)              ! RHS of the DSMC Method/ deltaV (npartmax, direction)

INTEGER                       :: CollisMode                 ! Mode of Collision:, ini_1
                                                            !    0: No Collisions (=free molecular flow with DSMC-Sampling-Routines)
                                                            !    1: Elastic Collision
                                                            !    2: Relaxation + Elastic Collision
                                                            !    3: Chemical Reactions

INTEGER                       :: PairE_vMPF(2)              ! 1: Pair chosen for energy redistribution
                                                            ! 2: partical with minimal MPF of this Pair
LOGICAL                       :: useDSMC
REAL    , ALLOCATABLE         :: PartStateIntEn(:,:)        ! (npartmax,1:3) with 2nd index: Evib, Erot, Eel
INTEGER , ALLOCATABLE	      :: PartElecQua(:)

INTEGER                         :: LD_MultiTemperaturMod   ! Modell choice for MultiTemperature
                                                              ! 0 = no MultiTemperature Modeling
                                                              ! 1 = LD1 see Paper
                                                              ! 2 = LD2
                                                              ! 3 = LD3

TYPE tSpecInit
  REAL                        :: TVib                       ! vibrational temperature, ini_1
  REAL                        :: TRot                       ! rotational temperature, ini_1
  REAL                        :: Telec                      ! electronic temperature, ini_1
END TYPE tSpecInit


TYPE tSpeciesDSMC                                           ! DSMC Species Param
  TYPE(tSpecInit),ALLOCATABLE :: Init(:) !   =>NULL()
  LOGICAL                     :: PolyatomicMol              ! Species is a polyatomic molecule
  INTEGER                     :: SpecToPolyArray            ! 
  CHARACTER(LEN=64)           :: Name                       ! Species Name, required for DSMCSpeciesElectronicDatabase
  INTEGER                     :: InterID                    ! Identification number (e.g. for DSMC_prob_calc), ini_2
                                                            !     1   : Atom
                                                            !     2   : Molecule
                                                            !     4   : Electron
                                                            !     10  : Atomic ion
                                                            !     15  : Atomic CEX/MEX ion
                                                            !     20  : Molecular ion
                                                            !     40  : Excited atom
                                                            !     100 : Excited atomic ion
                                                            !     200 : Excited molecule
                                                            !     400 : Excited molecular ion
  REAL                        :: TrefVHS                    ! VHS reference temp, ini_2
  REAL                        :: DrefVHS                    ! VHS reference diameter, ini_2
  REAL                        :: omegaVHS                   ! VHS exponent omega, ini_2
  INTEGER                     :: NumOfPro                   ! Number of Protons, ini_2
  REAL                        :: Eion_eV                    ! Energy of Ionisation in eV, ini_2
  REAL                        :: RelPolarizability          ! relative polarizability, ini_2
  INTEGER                     :: NumEquivElecOutShell       ! number of equivalent electrons in outer shell, ini2
  INTEGER                     :: Xi_Rot                     ! Rotational DOF 
  REAL                        :: CharaTVib                  ! Charac vibrational Temp, ini_2
  REAL                        :: Ediss_eV                   ! Energy of Dissosiation in eV, ini_2
  INTEGER                     :: MaxVibQuant                ! Max vib quantum number + 1
  INTEGER                     :: MaxElecQuant               ! Max elec quantum number + 1
  REAL                        :: RotRelaxProb               ! rotational relaxation probability
                                                            !this should be a value for every pair, and not fix!
  REAL                        :: VibRelaxProb               ! vibrational relaxation probability,
                                                            !this should be a value for every pair, and not fix!
  REAL                        :: ElecRelaxProb              ! electronic relaxation probability
                                                            !this should be a value for every transition, and not fix!
  REAL                        :: VFD_Phi3_Factor            ! Factor of Phi3 in VFD Method: Phi3 = 0 => VFD -> TCE, ini_2
#if (PP_TimeDiscMethod==42)
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: levelcounter         ! counter for electronic levels; only debug
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: dtlevelcounter       ! counter for produced electronic levels per timestep; only debug
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ElectronicTransition ! counter for electronic transition from state i to j
#endif
  REAL,ALLOCATABLE,DIMENSION(:,:) :: ElectronicState        ! Array with electronic State for each species
                                                            ! first  index: 1 - degeneracy & 2 - char. Temp,el
                                                            ! second index: energy level
END TYPE tSpeciesDSMC

TYPE(tSpeciesDSMC), ALLOCATABLE     :: SpecDSMC(:)          ! Species DSMC params (nSpec)

TYPE tDSMC 
  REAL                          :: GammaQuant               ! GammaQuant for zero point energy in Evib (perhaps also Erot), 
                                                            ! should be 0.5 or 0
  INTEGER(KIND=8), ALLOCATABLE  :: NumColl(:)               ! Number of Collision for each case + entire Collision number
  REAL                          :: TimeFracSamp             ! %-of simulation time for sampling
  INTEGER                       :: SampNum                  ! number of Samplingsteps
  INTEGER                       :: NumOutput                ! number of Outputs
  REAL                          :: DeltaTimeOutput          ! Time intervall for Output
  LOGICAL                       :: ReservoirSimu            ! Flag for reservoir simulation
  LOGICAL                       :: ReservoirSimuRate        ! Does not performe the collision.
                                                            ! Switch to enable to create reaction rates curves
  LOGICAL                       :: ReservoirRateStatistic   ! if false, calculate the reaction coefficient rate by the probability
                                                            ! Default Value is false
  INTEGER                       :: VibEnergyModel           ! Model for vibration Energy: 
                                                            !       0: SHO (default value!)
                                                            !       1: TSHO 
  INTEGER                       :: PartNumOctreeNode        ! Max Number of Particles per Octree Node
  LOGICAL                       :: UseOctree                ! Flag for Octree
  LOGICAL                       :: CalcSurfaceVal           ! Flag for calculation of surfacevalues like heatflux or force at walls
  LOGICAL                       :: CalcSurfaceTime          ! Flag for sampling in time-domain or iterations
  REAL                          :: CalcSurfaceSumTime       ! Flag for sampling in time-domain or iterations
  INTEGER                       :: CalcSurfCollis_NbrOfSpecies     ! Nbr. of Species to be counted for wall collisions (def. 0: all)
  LOGICAL,ALLOCATABLE           :: CalcSurfCollis_SpeciesFlags(:)  ! Species counted for wall collisions (def.: all species=T)
  LOGICAL                       :: CalcSurfCollis_OnlySwaps        ! count only wall collisions being SpeciesSwaps (def. F)
  LOGICAL                       :: CalcSurfCollis_Only0Swaps       ! count only wall collisions being delete-SpeciesSwaps (def. F)
  LOGICAL                       :: CalcSurfCollis_Output           ! Print sums of all counted wall collisions (def. F)
  REAL                          :: CollMean                 ! output of mean Collision Probability
  INTEGER                       :: CollMeanCount            ! counter for mean Collision Probability max
  REAL, ALLOCATABLE            :: CollProbOut(:,:)           ! Collision probability per cell for output
                                                            ! (1: Maximal collision prob, 2: Time-averaged mean collision prob)
  REAL, ALLOCATABLE            :: CollProbSamp(:)       ! Sampling of mean collision probability per cell
  LOGICAL                       :: ElectronicState          ! Flag for Electronic State of atoms and molecules
  CHARACTER(LEN=64)             :: ElectronicStateDatabase  ! Name of Electronic State Database | h5 file
  INTEGER                       :: NumPolyatomMolecs        ! Number of polyatomic molecules
  LOGICAL                       :: OutputMeshInit           ! Write Outputmesh (for const. pressure BC) at Init.
  LOGICAL                       :: OutputMeshSamp           ! Write Outputmesh (for const. pressure BC) 
                                                            ! with sampling values at t_analyze
END TYPE tDSMC

TYPE(tDSMC)                        :: DSMC
  
TYPE tBGGas
  INTEGER                       :: BGGasSpecies             ! Number which Species is Background Gas
  REAL                          :: BGGasDensity             ! Density of Background Gas
  REAL                          :: BGColl_SpecPartNum       ! PartNum of BGGas per cell   
  INTEGER                       :: BGMeanEVibQua            ! Mean EVib qua number for dissociation probability    
END TYPE tBGGas

TYPE(tBGGas)                        :: BGGas

TYPE tPairData
  REAL              :: CRela2                               ! squared relative velo of the particles in a pair
  REAL              :: Prob                                 ! colision probability
  INTEGER           :: iPart_p1                             ! first particle of the pair
  INTEGER           :: iPart_p2                             ! second particle of the pair
  INTEGER           :: PairType                             ! type of pair (=iCase, CollInf%Coll_Case)
  REAL, ALLOCATABLE :: Sigma(:)                             ! cross sections sigma of the pair
                                                            !       0: sigma total
                                                            !       1: sigma elast
                                                            !       2: sigma ionization
                                                            !       3: sigma exciation
  REAL              :: Ec                                   ! Collision Energy
  LOGICAL           :: NeedForRec                           ! Flag if pair is need for Recombination
END TYPE tPairData

TYPE(tPairData), ALLOCATABLE    :: Coll_pData(:)            ! Data of collision pairs into a cell (nPair)

TYPE tCollInf             ! informations of collision                                              
  INTEGER       , ALLOCATABLE    :: Coll_Case(:,:)          ! Case of species combination (Spec1, Spec2)
  INTEGER                        :: NumCase                 ! Number of possible collision combination
  INTEGER       , ALLOCATABLE    :: Coll_CaseNum(:)         ! number of species combination per cell Sab (number of cases)
  INTEGER       , ALLOCATABLE    :: Coll_SpecPartNum(:)     ! number of particle of species n per cell (nSpec)
  REAL          , ALLOCATABLE    :: Cab(:)                  ! species factor for cross section (number of case)
  INTEGER       , ALLOCATABLE    :: KronDelta(:)            ! (number of case)
  REAL          , ALLOCATABLE    :: FracMassCent(:,:)       ! mx/(my+mx) (nSpec, number of cases)
  REAL          , ALLOCATABLE    :: MassRed(:)              ! reduced mass (number of cases)
END TYPE

TYPE(tCollInf)               :: CollInf

TYPE tSampDSMC             ! DSMC sample                                              
  REAL                           :: PartV(3), PartV2(3)     ! Velocity, Velocity^2 (vx,vy,vz)
  REAL                           :: PartNum                 ! Particle Number
  INTEGER                        :: SimPartNum
  REAL                           :: ERot                    ! Rotational  Energy
  REAL                           :: EVib                    ! Vibrational Energy
  REAL                           :: EElec                   ! electronic  Energy
END TYPE

TYPE(tSampDSMC), ALLOCATABLE     :: SampDSMC(:,:)           ! DSMC sample array (number of Elements, nSpec)

TYPE tMacroDSMC           ! DSMC output 
  REAL                            :: PartV(4), PartV2(3)    ! Velocity, Velocity^2 (vx,vy,vz,|v|)
  REAL                            :: PartNum                ! Particle Number
  REAL                            :: Temp(4)                ! Temperature (Tx, Ty, Tz, Tt)
  REAL                            :: NumDens                ! Particle density        
  REAL                            :: TVib                   ! Vibrational Temp
  REAL                            :: TRot                   ! Rotational Temp
  REAL                            :: TElec                  ! Electronic Temp
END TYPE

TYPE(tMacroDSMC), ALLOCATABLE     :: MacroDSMC(:,:)         ! DSMC sample array (number of Elements, nSpec)

TYPE tReactInfo  
   REAL,  ALLOCATABLE             :: Xi_Total(:,:)          ! Total DOF of Reaction (quant num part1, quant num part2)  
   REAL,  ALLOCATABLE             :: Beta_Diss_Arrhenius(:,:) ! Beta_d for calculation of the Dissociation probability 
                                                            ! (quant num part1, quant num part2)   
   REAL,  ALLOCATABLE             :: Beta_Exch_Arrhenius(:) ! Beta_d for calculation of the Excchange reaction probability 
                                                            ! (quant num part1) 
   REAL,  ALLOCATABLE             :: Beta_Rec_Arrhenius(:)  ! Beta_d for calculation of the Recombination reaction probability 
                                                            ! (quant num part3)
END TYPE   

TYPE tChemReactions
  INTEGER                         :: NumOfReact             ! Number of possible Reactions
  LOGICAL, ALLOCATABLE            :: QKProcedure(:)         ! Defines if QK Procedure is selected
  INTEGER, ALLOCATABLE            :: QKMethod(:)            ! Recombination method for Q-K model (1 by Bird / 2 by Gallis)
  REAL,ALLOCATABLE,DIMENSION(:,:) :: QKCoeff                ! QKRecombiCoeff for Birds method
  REAL, ALLOCATABLE               :: NumReac(:)             ! Number of occured reactions for each reaction number
  INTEGER, ALLOCATABLE            :: ReacCount(:)           ! Counter of chemical reactions for the determination of rate
                                                            ! coefficient based on the reaction probabilities
!  INTEGER(KIND=8), ALLOCATABLE    :: NumReac(:)            ! Number of occured reactions for each reaction number
  LOGICAL                         :: MeanEVib_Necc          ! Flag if the MeanEVibQua_PerIter is necessary to calculate
  CHARACTER(LEN=1),ALLOCATABLE    :: ReactType(:)           ! Type of Reaction (reaction num)
                                                            !    i (electron impact ionization)
                                                            !    R (molecular recombination
                                                            !    D (molecular dissociation)
                                                            !    E (molecular exchange reaction)
                                                            !    x (simple charge exchange reaction)
  INTEGER, ALLOCATABLE            :: DefinedReact(:,:,:)    ! Defined Reaction 
                                                            !(reaction num; 1:reactant, 2:product;
                                                            !  1-3 spezieses of reactants and products,
                                                            ! 0: no spezies -> only 2 reactants or products)
  INTEGER, ALLOCATABLE            :: ReactCase(:,:)         ! Case of reaction in combination of (spec1, spec2)
  INTEGER, ALLOCATABLE            :: ReactNum(:,:,:)        ! Number of Reaction of (spec1, spec2, 
                                                            ! Case 1: Recomb: func. of species 3
                                                            ! Case 2: dissociation, only 1
                                                            ! Case 3: exchange reaction, only 1
                                                            ! Case 4: RN of 1. dissociation
                                                            !               2. exchange
                                                            ! Case 5: RN of 1. dissociation 1
                                                            !               2. dissociation 2
                                                            ! Case 6: ionization, only 1
                                                            ! Case 7: simple CEX, only 1
   REAL,  ALLOCATABLE             :: Arrhenius_Prefactor(:)     ! pre-exponential factor af of Arrhenius ansatz (nReactions)
   REAL,  ALLOCATABLE             :: Arrhenius_Powerfactor(:)   ! powerfactor bf of temperature in Arrhenius ansatz (nReactions)
   REAL,  ALLOCATABLE             :: EActiv(:)              ! activation energy (relative to k) (nReactions)
   REAL,  ALLOCATABLE             :: EForm(:)               ! heat of formation  (relative to k) (nReactions)
   REAL,  ALLOCATABLE             :: MeanEVib_PerIter(:)    ! MeanEVib per iteration for calculation of 
   INTEGER,  ALLOCATABLE          :: MeanEVibQua_PerIter(:) ! MeanEVib per iteration for calculation of 
                                                            ! xi_vib per cell (nSpecies)
   REAL,  ALLOCATABLE             :: CEXa(:)                ! CEX log-factor (g-dep. cross section in Angstrom (nReactions)
   REAL,  ALLOCATABLE             :: CEXb(:)                ! CEX const. factor (g-dep. cross section in Angstrom (nReactions)
   REAL,  ALLOCATABLE             :: MEXa(:)                ! MEX log-factor (g-dep. cross section in Angstrom (nReactions)
   REAL,  ALLOCATABLE             :: MEXb(:)                ! MEX const. factor (g-dep. cross section in Angstrom (nReactions)
   INTEGER                       :: RecombParticle = 0      ! P. Index for Recombination, if zero -> no recomb particle avaible
   INTEGER                       :: nPairForRec
   INTEGER                       :: nPartForRec
   TYPE(tReactInfo), ALLOCATABLE  :: ReactInfo(:)           ! Informations of Reactions (nReactions)   
END TYPE

TYPE tTreeNode
!  TYPE (tTreeNode), POINTER       :: One, Two, Three, Four, Five, Six, Seven, Eight !8 Childnodes of Octree Treenode
  TYPE (tTreeNode), POINTER       :: ChildNode              !8 Childnodes of Octree Treenode
  REAL                            :: MidPoint(1:3)          ! approx Middle Point of Treenode
  INTEGER                         :: PNum_Node              ! Particle Number of Treenode
  INTEGER, ALLOCATABLE            :: iPartIndx_Node(:)      ! Particle Index List of Treenode
  INTEGER                         :: PairNum_Node           ! Number of Particle Pairs            
END TYPE

TYPE(tChemReactions)              :: ChemReac

TYPE tSampWall             ! DSMC sample for Wall                                             
  REAL                           :: Energy(9)               ! 1-3 E_tra (pre, wall, re),
                                                            ! 4-6 E_rot (pre, wall, re),
                                                            ! 7-9 E_vib (pre, wall, re)
  REAL                           :: Force(3)                ! x, y, z direction
  REAL, ALLOCATABLE              :: Counter(:)              ! Wall-Collision counter
END TYPE

TYPE(tSampWall), ALLOCATABLE     :: SampWall(:)             ! Wall sample array (number of BC-Sides)
#ifdef MPI
TYPE(tSampWall), ALLOCATABLE     :: SampWallHaloCell(:)     ! Wall sample array (number of BC-HALO-Sides)
#endif

TYPE tMacroSurfaceVal                                       ! DSMC sample for Wall    
  REAL                           :: Heatflux                ! 
  REAL                           :: Force(3)                ! x, y, z direction
  REAL, ALLOCATABLE              :: Counter(:)              ! Wall-Collision counter of all Species
  REAL                           :: CounterOut              ! Wall-Collision counter for Output
END TYPE

TYPE(tMacroSurfaceVal), ALLOCATABLE     :: MacroSurfaceVal(:) ! Wall sample array (number of BC-Sides)

 
REAL                              :: realtime               ! realtime of simulation

TYPE tPolyatomMolDSMC !DSMC Species Param
  LOGICAL                         :: LinearMolec            ! Is a linear Molec?
  INTEGER                         :: NumOfAtoms             ! Number of Atoms in Molec
  INTEGER                         :: VibDOF                 ! DOF in Vibration, equals number of independent SHO's
  REAL, ALLOCATABLE               :: CharaTVibDOF(:)        ! Chara TVib for each DOF
  INTEGER,ALLOCATABLE             :: LastVibQuantNums(:)    ! Last Quantum Numbers for vibrational insering
  INTEGER, ALLOCATABLE            :: MaxVibQuantDOF(:)      ! Max Vib Quant for each DOF
END TYPE

TYPE (tPolyatomMolDSMC), ALLOCATABLE    :: PolyatomMolDSMC(:)        ! Infos for Polyatomic Molecule

TYPE tPolyatomMolVibQuant !DSMC Species Param
  INTEGER, ALLOCATABLE               :: Quants(:)            ! Vib quants of each DOF for each particle
END TYPE

TYPE (tPolyatomMolVibQuant), ALLOCATABLE    :: VibQuantsPar(:)

INTEGER                           :: nOutput                 ! output counter for DSMC
INTEGER(KIND=8)                   :: iter_loc, iter_macvalout, istep

INTEGER                           :: nSurfSample             ! polynomial degree of surface supersampling
REAL,ALLOCATABLE                  :: XiEq_Surf(:)            ! position of equidistant interpolation points on surface
REAL                              :: deltaXiEQ_Surf              ! delta of equidistant surface sampling

TYPE tSampleCartmesh_VolWe
  REAL                                  :: BGMdeltas(3)       ! Backgroundmesh size in x,y,z
  REAL                                  :: FactorBGM(3)       ! Divider for BGM (to allow real numbers)
  REAL                                  :: BGMVolume          ! Volume of a BGM Cell
  INTEGER,ALLOCATABLE               :: GaussBGMIndex(:,:,:,:,:) ! Background mesh index of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
  REAL,ALLOCATABLE                  :: GaussBGMFactor(:,:,:,:,:) ! BGM factor of gausspoints (1:3,PP_N,PP_N,PP_N,nElems)
  INTEGER                               :: BGMminX            ! Local minimum BGM Index in x
  INTEGER                               :: BGMminY            ! Local minimum BGM Index in y
  INTEGER                               :: BGMminZ            ! Local minimum BGM Index in z
  INTEGER                               :: BGMmaxX            ! Local maximum BGM Index in x
  INTEGER                               :: BGMmaxY            ! Local maximum BGM Index in y
  INTEGER                               :: BGMmaxZ            ! Local maximum BGM Index in z
  INTEGER, ALLOCATABLE                  :: PeriodicBGMVectors(:,:)           ! = periodic vectors in backgroundmesh coords
  LOGICAL                               :: SelfPeriodic
  REAL, ALLOCATABLE                    :: BGMVolumes(:,:,:)
  REAL, ALLOCATABLE                    :: BGMVolumes2(:,:,:)
  LOGICAL, ALLOCATABLE                 :: isBoundBGCell(:,:,:)
  INTEGER                               :: OrderVolInt
  REAL, ALLOCATABLE                    :: x_VolInt(:)
  REAL, ALLOCATABLE                    :: w_VolInt(:)
#ifdef MPI
  TYPE(tPartMPIConnect)        , ALLOCATABLE :: MPIConnect(:)             ! MPI connect for each process
#endif
END TYPE

TYPE (tSampleCartmesh_VolWe) DSMCSampVolWe

TYPE tDSMCSampNearInt
  REAL,ALLOCATABLE                      :: GaussBorder(:)     ! 1D coords of gauss points in -1|1 space
END TYPE

TYPE (tDSMCSampNearInt) DSMCSampNearInt

TYPE tDSMCSampCellVolW
  REAL,ALLOCATABLE                      :: xGP(:)     
END TYPE

TYPE (tDSMCSampCellVolW) DSMCSampCellVolW

TYPE tHODSMC
  LOGICAL                 :: HODSMCOutput         !High Order DSMC Output
  INTEGER                 :: nOutputDSMC          !HO DSMC output order
  REAL,ALLOCATABLE        :: DSMC_xGP(:,:,:,:,:)  ! XYZ positions (first index 1:3) of the volume Gauss Point
  REAL,ALLOCATABLE        :: DSMC_wGP(:)
  CHARACTER(LEN=256)      :: SampleType
  CHARACTER(LEN=256)      :: NodeType
  REAL,ALLOCATABLE        :: sJ(:,:,:,:)
END TYPE tHODSMC

TYPE(tHODSMC)             :: HODSMC
REAL,ALLOCATABLE          :: DSMC_HOSolution(:,:,:,:,:,:) !1:3 v, 4:6 v^2, 7 dens, 8 Evib, 9 erot, 10 eelec
!===================================================================================================================================
END MODULE MOD_DSMC_Vars
