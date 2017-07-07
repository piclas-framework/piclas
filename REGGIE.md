## Boltzplatz Regressioncheck 

## How to execute

## Which examples?

| **No.** |            **Example**            | **When** |     **CMAKE-CONFIG**    |      **Feature**      |         **Execution**                   |           **Comparing**          |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:---------------------:|:---------------------------------------:|:--------------------------------:|
|    1    |             run_basic             |  checkin |       maxwell,RK4       |      DG-Operator      |        nProcs=1,2,5,8                   |              L2,Linf             |
|    2    |    feature_poisson_powerdensity   |  nightly | poisson, Crank-Nicolson | Implicit, CalcTimeAvg |  DoRefMapping=T/F, nProcs=2             |       Final TimeAvg, h5diff      |
|    3    |   feature_tracking_DSMC_ANSA_box  |  nightly |      maxwell, DSMC      |        Tracking       | DoRefMapping=T,F, nProcs=1,2            | PartInt, PartPos in Bounding Box |
|    4    |    feature_tracking_DSMC_curved   |  nightly |       maxwell,DSMC      |        Tracking       |  DoRefMapping=T, nProcs=1,2             |              PartInt             |
|    5    |    feature_tracking_DSMC_mortar   |  nightly |       maxwell,DSMC      |  Tracking with Mortar | DoRefMapping=T,F, nProcs=1,2            | PartInt, PartPos in Bounding Box |
|    6    |   feature_tracking_DSMC_periodic  |  nightly |       maxwell,DSMC      |        Tracking       | DoRefMapping=T,F, nProcs=1,2            | PartInt, PartPos in Bounding Box |
|    7    | feature_tracking_DSMC_sphere_soft |  nightly |       maxwell,DSMC      |        Tracking       |                                         |                                  |
|    8    | feature_poisson_powerdensity      |  nightly |       poisson,RK4,CN    |        Poisson-PIC    | DoRefMapping=1,2, nProcs=2, CN implicit |  TimeAvg                         |