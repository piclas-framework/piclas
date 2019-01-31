## PICLas Regressioncheck 

## How to execute

## Which examples?

## Reggie2.0 (new reggie)

| **No.** |     **Check**     | **When** |  **CMAKE-CONFIG**   |    **Examples**    | **Feature**  | **Execution** | **Comparing** |
| :-----: | :---------------: | :------: | :-----------------: | :----------------: | :----------: | :-----------: | :-----------: |
|    1    | run_basic (flexi) | checkin  |       default       |   freestream_2D    | DG-Operator  |    MPI=1,2    |      L2       |
|         |                   |          |                     |   freestream_3D    | DG-Operator  |    MPI=1,2    |      L2       |
|    2    | convtest (flexi)  | nighlty  |     FLEXI_2D=ON     |        h_2D        | h-convergece |    single     |      L2       |
|         |                   |          |                     |        h_3D        | h-convergece |    single     |      L2       |
|         |                   |          |     FLEXI_FV=ON     |      h_3D_FV       | h-convergece |    single     |      L2       |
|         |                   |          |                     |    h_3D_mortar     | h-convergece |    single     |      L2       |
|         |                   |          | FLEXI_PARABOLIC=OFF | h_3D_parabolic_off | h-convergece |    single     |      L2       |
|         |                   |          |     FLEXI_2D=ON     |        p_2D        | p-convergece |    single     |      L2       |
|         |                   |          |                     |        p_3D        | p-convergece |    single     |      L2       |
|         |                   |          |                     |    p_3D_mortar     | p-convergece |    single     |      L2       |
|         |                   |          | FLEXI_PARABOLIC=OFF | p_3D_parabolic_off | p-convergece |    single     |      L2       |
|         |                   |          |                     |                    |              |               |               |
|         |                   |          |                     |                    |              |               |               |
|         |                   |          |                     |                    |              |               |               |


## Analyze routines

see [the reggie repository](https://gitlab.iag.uni-stuttgart.de/reggie/reggie/blob/master/README.md)

# List of Examples

## Check-in

| **No.** |         **Example**          |    **CMAKE-CONFIG**    |           **Feature**           | **Execution**  |     **Comparing**      |
| :-----: | :--------------------------: | :--------------------: | :-----------------------------: | :------------: | :--------------------: |
|   01    |          run_basic           |      maxwell,RK4       |           DG-Operator           | nProcs=1,2,5,8 |        L2,Linf         |
|   02    |         CHE_maxwell          |  pure maxwell DG,RK4   |           DG-Operator           | nProcs=1,2,5,8 |        L2,Linf         |
|   03    |         CHE_poisson          |      Poisson,RK3       |           DG-Operator           | nProcs=1,2,5,8 |        L2,Linf         |
|   05    | CHE_PIC_gyrotron_variable_Bz |    maxwell,PIC,RK4     |        PIC, variable Bz         |   nProcs=1,2   | Database.csv, relative |
|   06    | CHE_PIC_single_particle_PML  |    maxwell,PIC,RK4     |            PIC, PML             |    particle    |   nProcs=1,2,5,8,10    |
|   07    |     CHE_PIC_IMD_coupling     |    maxwell,PIC,RK4     |   mapping from IMP to PICLas    |    nProcs=1    |    PartPata in Box     |
|   08    |        CHE_DSMC_check        |          DSMC          |                                 |                |                        |
|   09    |  CHE_PIC_maxwell_implicitBC  | maxwell,PIC,ImplicitO4 | Implicit reflective particle BC |    nProcs=1    |   Particle Position    |

### CHE_DSMC_check

| **No.** |              **Example**              | **CMAKE-CONFIG** |            **Feature**            | **Execution** |           **Comparing**            | **Readme** |
| :-----: | :-----------------------------------: | :--------------: | :-------------------------------: | :-----------: | :--------------------------------: | :--------: |
|  08-1   | BC_surfaceflux_adaptive_constPressure |                  | SurfaceFlux with AdaptiveType=1/2 |   nProcs=4    |        Integrated mass flux        |            |
|  08-2   | BC_surfaceflux_adaptive_constMassflow |                  | SurfaceFlux with AdaptiveType=3,4 |   nProcs=1    |        Integrated mass flux        |            |
|  08-3   |              BC_porousBC              |                  | PorousBC as a pump with 2 species |   nProcs=3    | Total # of removed part through BC |            |
|  08-4   |                 cube                  |                  |          Collismode=2,3           |   nProcs=2    |                                    |            |

## Nightly

Overview of the test cases performed during the nightly regression testing

| **No.** |         **Example**          |    **CMAKE-CONFIG**     |           **Feature**            |                **Execution**                |         **Comparing**          | **Readme** |
| :-----: | :--------------------------: | :---------------------: | :------------------------------: | :-----------------------------------------: | :----------------------------: | :--------: |
|   10    |        [NIG_Reservoir](https://gitlab.com/piclas/piclas/blob/update.ionization.routines/REGGIE.md#nig_reservoir)         |      maxwell, DSMC      | Relaxation, (Surface-) Chemistry |                                             |                                |            |
|   11    |      [NIG_tracking_DSMC](https://gitlab.com/piclas/piclas/blob/update.ionization.routines/REGGIE.md#nig_tracking_dsmc)       |      maxwell, DSMC      |             Tracking             |                                             |                                |            |
|   12    |   NIG_PIC_maxwell_bgfield    |     maxwell,PIC,RK4     |   External Background-field,h5   |                  nProcs=2                   |          DG_Solution           |            |
|   13    | feature_poisson_powerdensity | Poisson, Crank-Nicolson |      Implicit, CalcTimeAvg       |         DoRefMapping=T/F, nProcs=2          |     Final TimeAvg, h5diff      |            |
|   14    | feature_poisson_powerdensity |       poisson,CN        |           CalcTimeAvg            |   DoRefMapping=1,2, nProcs=2, CN implicit   |            TimeAvg             |            |
|   15    |  feature_emission_gyrotron   |       maxwell,RK4       |       Part-Inflow,TimeDep        | N=1,3,6,9,10, nProcs=1,2,10,25, gyro-circle |   LineIntegration of nPartIn   |            |
|   16    |   feature_TWT_recordpoints   |       maxwell,RK4       |          RPs, ExactFlux          |     nProcs=1,4, RPs, interior TE-Inflow     |       RP_State, RP_Daata       |            |
|   17    | feature_PIC_HDG_plasma_wave  |     poisson,RK4,CN      |  Poisson-PIC,Shape-Function-1D   |            nProcs=2, Imex for CN            | W_el LineIntegration over 2Per |            |

### NIG_Reservoir

Reservoir (heat bath) simulations *Link to build*

| **No.** |            **Example**            | **CMAKE-CONFIG** |        **Feature**         | **Execution** | **Comparing** | **Readme** |
| :-----: | :-------------------------------: | :--------------: | :------------------------: | :-----------: | :-----------: | :--------: |
|  10-1   |    CHEM_dissocication_rate_CH4    |                  |   Dissociation reactions   |   nProcs=1    |               |            |
|  10-2   |      CHEM_exchange_rate_CH3       |                  |     Exchange reactions     |   nProcs=1    |               |            |
|  10-3   |      CHEM_multi_ionization_N      |                  | Electron-impact ionization |   nProcs=1    |               |            |
|  10-4   | CHEM_recombination_rate_CH3_and_H |                  |       Recombination        |   nProcs=1    |               |            |

### NIG_tracking_DSMC

Testing of different tracking routines *Link to build*

| **No.** |    **Example**    | **CMAKE-CONFIG** |          **Feature**          |                 **Execution**                 |          **Comparing**           | **Readme** |
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

| **No.** |           **Example**           |    **CMAKE-CONFIG**    |          **Feature**           |         **Execution**         |         **Comparing**          |
| :-----: | :-----------------------------: | :--------------------: | :----------------------------: | :---------------------------: | :----------------------------: |
|   18    | feature_PIC_maxwell_plasma_wave | maxwell,RK4,ImplicitO4 | Maxwell-PIC,SF1D, FastPeriodic | nProcs=2, IMEX for ImplicitO4 | W_el LineIntegration over 2Per |