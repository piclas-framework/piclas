# PICLas Regression Testing

PICLas utilizes the Reggie2.0 toolbox for regression testing. A detailed documentation on its usage is available [at this repository](https://gitlab.com/reggie2.0/reggie2.0/blob/master/README.md). A list detailing the test cases and which features are tested is given below.

# List of Cases

## Check-in

Overview of the test cases performed after a commit.

| **No.** |           **Case**           |                 **CMAKE-CONFIG**                 |                              **Feature**                               | **Execution**  |     **Comparing**      |                        **Readme**                         |
| :-----: | :--------------------------: | :----------------------------------------------: | :--------------------------------------------------------------------: | :------------: | :--------------------: | :-------------------------------------------------------: |
|   01    |          run_basic           |                   maxwell,RK4                    |                              DG-Operator                               | nProcs=1,2,5,8 |        L2,Linf         |                                                           |
|   02    |         CHE_maxwell          |               pure maxwell DG,RK4                |                              DG-Operator                               | nProcs=1,2,5,8 |        L2,Linf         |                                                           |
|   03    |         CHE_poisson          |                   Poisson,RK3                    |                              DG-Operator                               | nProcs=1,2,5,8 |        L2,Linf         |                                                           |
|   05    | CHE_PIC_gyrotron_variable_Bz |                 maxwell,PIC,RK4                  |                            PIC, variable Bz                            |   nProcs=1,2   | Database.csv, relative |                                                           |
|   06    | CHE_PIC_single_particle_PML  |                 maxwell,PIC,RK4                  |                                PIC, PML                                |    particle    |   nProcs=1,2,5,8,10    |                                                           |
|   07    |     CHE_PIC_IMD_coupling     |                 maxwell,PIC,RK4                  |                       mapping from IMP to PICLas                       |    nProcs=1    |    PartPata in Box     |                                                           |
|   08    |    [CHE_DSMC](#che_dsmc)     |                       DSMC                       |                                                                        |                |                        |                                                           |
|   09    |  CHE_PIC_maxwell_implicitBC  |              maxwell,PIC,ImplicitO4              |                    Implicit reflective particle BC                     |    nProcs=1    |   Particle Position    |                                                           |
|   10    |           CHE_BGK            | [BGK](regressioncheck/checks/CHE_BGK/builds.ini) | Relax to thermal equi. with ESBGK/SBGK, continuous/quantized vibration |    nProcs=1    |  T_rot,T_vib,T_trans   | [Link](regressioncheck/checks/CHE_BGK/RELAX_N2/readme.md) |

#### CHE_DSMC

Small test cases to check DSMC features: [Link to build](regressioncheck/checks/CHE_DSMC/builds.ini).

| **No.** |               **Case**                | **CMAKE-CONFIG** |            **Feature**            | **Execution** |           **Comparing**            | **Readme** |
| :-----: | :-----------------------------------: | :--------------: | :-------------------------------: | :-----------: | :--------------------------------: | :--------: |
|  08-1   | BC_surfaceflux_adaptive_constPressure |                  | SurfaceFlux with AdaptiveType=1/2 |   nProcs=4    |        Integrated mass flux        |            |
|  08-2   | BC_surfaceflux_adaptive_constMassflow |                  | SurfaceFlux with AdaptiveType=3,4 |   nProcs=1    |        Integrated mass flux        |            |
|  08-3   |              BC_porousBC              |                  | PorousBC as a pump with 2 species |   nProcs=3    | Total # of removed part through BC |            |
|  08-4   |                 cube                  |                  |          Collismode=2,3           |   nProcs=2    |                                    |            |

## Nightly

Overview of the test cases performed during the nightly regression testing.

| **No.** |                **Case**                 |    **CMAKE-CONFIG**     |           **Feature**            |                **Execution**                |         **Comparing**          | **Readme** |
| :-----: | :-------------------------------------: | :---------------------: | :------------------------------: | :-----------------------------------------: | :----------------------------: | :--------: |
|   10    |     [NIG_Reservoir](#nig_reservoir)     |      maxwell, DSMC      | Relaxation, (Surface-) Chemistry |                                             |                                |            |
|   11    | [NIG_tracking_DSMC](#nig_tracking_dsmc) |      maxwell, DSMC      |             Tracking             |                                             |                                |            |
|   12    |         NIG_PIC_maxwell_bgfield         |     maxwell,PIC,RK4     |   External Background-field,h5   |                  nProcs=2                   |          DG_Solution           |            |
|   13    |      feature_poisson_powerdensity       | Poisson, Crank-Nicolson |      Implicit, CalcTimeAvg       |         DoRefMapping=T/F, nProcs=2          |     Final TimeAvg, h5diff      |            |
|   14    |      feature_poisson_powerdensity       |       poisson,CN        |           CalcTimeAvg            |   DoRefMapping=1,2, nProcs=2, CN implicit   |            TimeAvg             |            |
|   15    |        feature_emission_gyrotron        |       maxwell,RK4       |       Part-Inflow,TimeDep        | N=1,3,6,9,10, nProcs=1,2,10,25, gyro-circle |   LineIntegration of nPartIn   |            |
|   16    |        feature_TWT_recordpoints         |       maxwell,RK4       |          RPs, ExactFlux          |     nProcs=1,4, RPs, interior TE-Inflow     |       RP_State, RP_Daata       |            |
|   17    |       feature_PIC_HDG_plasma_wave       |     poisson,RK4,CN      |  Poisson-PIC,Shape-Function-1D   |            nProcs=2, Imex for CN            | W_el LineIntegration over 2Per |            |

#### NIG_Reservoir

Testing more complex DSMC routines with reservoir (heat bath) simulations: [Link to build](regressioncheck/checks/NIG_Reservoir/builds.ini).

| **No.** |            **Case**            | **CMAKE-CONFIG** |                            **Feature**                            | **Execution** | **Comparing** |                                      **Readme**                                       |
| :-----: | :----------------------------: | :--------------: | :---------------------------------------------------------------: | :-----------: | :-----------: | :-----------------------------------------------------------------------------------: |
|  10-x   |  CHEM_RATES_dissocication_CH4  |                  |       TCE rates for a dissociation: CH4 + M -> CH3 + H + M        |   nProcs=1    |               |  [Link](regressioncheck/checks/NIG_Reservoir/CHEM_RATES_dissocication_CH4/readme.md)  |
|  10-x   |   CHEM_RATES_exchange_CH4_H    |                  |          TCE rates for an exchange: CH4 + H <-> CH3 + H2          |   nProcs=1    |               |    [Link](regressioncheck/checks/NIG_Reservoir/CHEM_RATES_exchange_CH3/readme.md)     |
|  10-x   |     CHEM_RATES_recomb_CH4      |                  |           TCE rates for a recombination: CH3 + H -> CH4           |   nProcs=1    |               |     [Link](regressioncheck/checks/NIG_Reservoir/CHEM_RATES_recomb_CH4/readme.md)      |
|  10-x   | CHEM_RATES_ionization-recomb_H |                  | QK rates for ionization and recombination: H + e <-> HIon + e + e |   nProcs=1    |               | [Link](regressioncheck/checks/NIG_Reservoir/CHEM_RATES_ionization-recomb_H/readme.md) |
|  10-x   | CHEM_multi-ionization_C_to_C6+ |                  |         Impact ionization, from neutral to fully ionized          |   nProcs=1    |               | [Link](regressioncheck/checks/NIG_Reservoir/CHEM_multi-ionization_C_to_C6+/readme.md) |
|  10-x   |            RELAX_N2            |                  |          Rotational, vibrational, electronic relaxation           |   nProcs=1    |               |            [Link](regressioncheck/checks/NIG_Reservoir/RELAX_N2/readme.md)            |
|  10-x   |           RELAX_CO2            |                  |                Rotational, vibrational relaxation                 |   nProcs=1    |               |           [Link](regressioncheck/checks/NIG_Reservoir/RELAX_CO2/readme.md)            |
|  10-x   |          RELAX_N2Ion           |                  |          Rotational, vibrational, electronic relaxation           |   nProcs=1    |               |          [Link](regressioncheck/checks/NIG_Reservoir/RELAX_N2Ion/readme.md)           |

#### NIG_tracking_DSMC

Testing of different tracking routines with DSMC: [Link to build](regressioncheck/checks/NIG_tracking_DSMC/builds.ini).

| **No.** |     **Case**      | **CMAKE-CONFIG** |          **Feature**          |                 **Execution**                 |          **Comparing**           | **Readme** |
| :-----: | :---------------: | :--------------: | :---------------------------: | :-------------------------------------------: | :------------------------------: | :--------: |
|  11-1   |     ANSA box      |                  |                               |         DoRefMapping=T,F, nProcs=1,2          | PartInt, PartPos in bounding box |            |
|  11-2   |      curved       |                  |                               |         DoRefMapping=T  , nProcs=1,2          | PartInt with relative tolerance  |            |
|  11-3   |      mortar       |                  |                               |         DoRefMapping=T,F, nProcs=1,2          | PartInt, PartPos in bounding box |            |
|  11-4   |     periodic      |                  |                               |       DoRefMapping=T,F, nProcs=1,2,5,10       | PartInt, PartPos in bounding box |            |
|  11-5   |  periodic_2cells  |                  |                               |  DoRefMapping=T,F;TriaTracking=T,F, nProcs=1  |     PartPos in bounding box      |            |
|  11-6   |    semicircle     |                  |                               |         DoRefMapping=T,F, nProcs=1,2          |     PartPos in bounding box      |            |
|  11-7   |    sphere_soft    |                  |                               | DoRefMapping=T;RefMappingGuess=1,3,nProcs=1,2 |     PartPos in bounding box      |            |
|  11-8   | tracing_cylinder1 |                  | mortar,curved,single particle |           DoRefMapping=F, nProcs=1            |    PartPos-X in bounding box     |            |
|  11-9   | tracing_cylinder2 |                  | mortar,curved,single particle |           DoRefMapping=F, nProcs=1            |    PartPos-X in bounding box     |            |

## Weekly

Overview of the testcases performed every week.

| **No.** |            **Case**             |    **CMAKE-CONFIG**    |                **Feature**                 |         **Execution**         |         **Comparing**          |                                  **Readme**                                   |
| :-----: | :-----------------------------: | :--------------------: | :----------------------------------------: | :---------------------------: | :----------------------------: | :---------------------------------------------------------------------------: |
|   18    | feature_PIC_maxwell_plasma_wave | maxwell,RK4,ImplicitO4 |       Maxwell-PIC,SF1D, FastPeriodic       | nProcs=2, IMEX for ImplicitO4 | W_el LineIntegration over 2Per |                                                                               |
|    x    |     CHEM_EQUI_ionization_H      |     DSMC Reservoir     | Relaxation into equilibrium with chemistry |           nProcs=1            |      PartAnalyze_ref.csv       | [Link](regressioncheck/checks/WEK_Reservoir/CHEM_EQUI_ionization_H/readme.md) |
