## PICLas Regressioncheck 

## How to execute

## Which examples?

## Reggie2.0 (new reggie)

| **No.** |            **Check**              | **When** |     **CMAKE-CONFIG**    |      **Examples**              |      **Feature**               |         **Execution**                       |           **Comparing**          |
|:-------:|:----------------------------------|:--------:|:-----------------------:|:------------------------------:|:------------------------------:|:-------------------------------------------:|:--------------------------------:|
|    1    | run_basic (flexi)                 | checkin  | default                 | freestream_2D                  |  DG-Operator                   |  MPI=1,2                                    | L2                               |
|         |                                   |          |                         | freestream_3D                  |  DG-Operator                   |  MPI=1,2                                    | L2                               |
|    2    | convtest (flexi)                  | nighlty  | FLEXI_2D=ON             | h_2D                           |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | h_3D                           |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_FV=ON             | h_3D_FV                        |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | h_3D_mortar                    |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_PARABOLIC=OFF     | h_3D_parabolic_off             |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_2D=ON             | p_2D                           |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | p_3D                           |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | p_3D_mortar                    |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_PARABOLIC=OFF     | p_3D_parabolic_off             |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         |                                |                                |                                             |                                  |
|         |                                   |          |                         |                                |                                |                                             |                                  |
|         |                                   |          |                         |                                |                                |                                             |                                  |


## Analyze routines

see [the reggie repository](https://gitlab.iag.uni-stuttgart.de/reggie/reggie/blob/master/README.md)

## List of Examples

| **No.** |            **Example**            | **When** |     **CMAKE-CONFIG**    |      **Feature**                |         **Execution**                          |           **Comparing**          |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|   01    |             run_basic             |  checkin |       maxwell,RK4       |      DG-Operator                |        nProcs=1,2,5,8                          |              L2,Linf             |
|   02    |             CHE_maxwell           |  checkin | pure maxwell DG,RK4     |      DG-Operator                |        nProcs=1,2,5,8                          |              L2,Linf             |
|   03    |             CHE_poisson           |  checkin |       Poisson,RK3       |      DG-Operator                |        nProcs=1,2,5,8                          |              L2,Linf             |
|   05    |    CHE_PIC_gyrotron_variable_Bz   |  checkin |    maxwell,PIC,RK4      |  PIC, variable Bz               |        nProcs=1,2                              |    Database.csv, relative        |
|   06    |    CHE_PIC_single_particle_PML    |  checkin |    maxwell,PIC,RK4      |  PIC, PML + particle            |        nProcs=1,2,5,8,10                       |    DG_Solution in State          |
|   07    |    CHE_PIC_IMD_coupling           |  checkin |    maxwell,PIC,RK4      |   mapping from IMP to PICLas    |        nProcs=1                                |    PartPata in Box               |
|   08    |    CHE_DSMC_check                 |  checkin |      DSMC               |  Collismode=2,3                 |        nProcs=2                                |                                  |
|   09    |    CHE_PIC_maxwell_implicitBC     |  checkin |  maxwell,PIC,ImplicitO4 | Implicit reflective particle BC |        nProcs=1                                |    Particle Position             |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|   10    |    NIG_Reservoir                  |  nightly |     DSMC                | Dissociation, recombination,    |        nProcs=1                                |    Particle Position             |
|         |                                   |          |                         | Exchange reactions              |        nProcs=1                                |    Database.csv                  |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|   11    |   NIG_tracking_DSMC               |  nightly |      maxwell, DSMC      |        Tracking                 |                                                |                                  |
|   11-1  |      ANSA box                     |          |                         |                                 | DoRefMapping=T,F, nProcs=1,2                   | PartInt, PartPos in bounding box |
|   11-2  |      curved                       |          |                         |                                 | DoRefMapping=T  , nProcs=1,2                   | PartInt with relative tolerance  |
|   11-3  |      mortar                       |          |                         |                                 | DoRefMapping=T,F, nProcs=1,2                   | PartInt, PartPos in bounding box |
|   11-4  |      periodic                     |          |                         |                                 | DoRefMapping=T,F, nProcs=1,2,5,10              | PartInt, PartPos in bounding box |
|   11-5  |      periodic_2cells              |          |                         |                                 | DoRefMapping=T,F;TriaTracking=T,F, nProcs=1    |  PartPos in bounding box         |
|   11-6  |      semicircle                   |          |                         |                                 | DoRefMapping=T,F, nProcs=1,2                   |  PartPos in bounding box         |
|   11-7  |      sphere_soft                  |          |                         |                                 | DoRefMapping=T;RefMappingGuess=1,3,nProcs=1,2  |  PartPos in bounding box         |
|   11-8  |      tracing_cylinder1            |          |                         |  mortar,curved,single particle  | DoRefMapping=F, nProcs=1                       |  PartPos-X in bounding box       |
|   11-9  |      tracing_cylinder2            |          |                         |  mortar,curved,single particle  | DoRefMapping=F, nProcs=1                       |  PartPos-X in bounding box       |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|   12    |    NIG_PIC_maxwell_bgfield        |  nightly |    maxwell,PIC,RK4      | External Background-field,h5    |        nProcs=2                                |    DG_Solution                   |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:-------------------------------:|:----------------------------------------------:|:--------------------------------:|
|    2    |    feature_poisson_powerdensity   |  nightly | Poisson, Crank-Nicolson | Implicit, CalcTimeAvg          |  DoRefMapping=T/F, nProcs=2                 |       Final TimeAvg, h5diff      |
|    8    | feature_poisson_powerdensity      |  nightly |       poisson,CN        |      CalcTimeAvg               | DoRefMapping=1,2, nProcs=2, CN implicit     |  TimeAvg                         |
|    9    |   feature_emission_gyrotron       |  nightly |       maxwell,RK4       | Part-Inflow,TimeDep            | N=1,3,6,9,10, nProcs=1,2,10,25, gyro-circle |  LineIntegration of nPartIn      |
|   10    |      feature_TWT_recordpoints     |  nightly |       maxwell,RK4       | RPs, ExactFlux                 | nProcs=1,4, RPs, interior TE-Inflow         |  RP_State, RP_Daata              |
|   11    |    feature_PIC_HDG_plasma_wave    |  nightly |       poisson,RK4,CN    | Poisson-PIC,Shape-Function-1D  | nProcs=2, Imex for CN                       |  W_el LineIntegration over 2Per  |
|   12    |  feature_PIC_maxwell_plasma_wave  |  weekly  | maxwell,RK4,ImplicitO4  | Maxwell-PIC,SF1D, FastPeriodic | nProcs=2, IMEX for ImplicitO4               |  W_el LineIntegration over 2Per  |

  stage: reggie_reservoir_nightly
