# Code Extension

This section shall describe how to extend the code in a general way e.g. to implement new input and output parameters.

## Surface Sampling & Output

The surface sampling values are stored in the `SampWallState` array and the derived output variables in `MacroSurfaceVal`, which can be extended to include new **optional** variables and to exploit the already implemented MPI communication. The variables `SurfSampSize` and `SurfOutputSize` define the size of the arrays is set in `InitParticleBoundarySampling` in `piclas/src/particles/boundary/particle_boundary_sampling.f90` with default values as given in `piclas.h`

    SurfSampSize = SAMPWALL_NVARS+nSpecies
    SurfOutputSize = MACROSURF_NVARS

To add optional variables, you need to increase the sampling and/or output indices as shown in the example

    IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
      SurfSampSize = SurfSampSize + 1
      SWIStickingCoefficient = SurfSampSize
      SurfOutputSize = SurfOutputSize + 1
    END IF

To be able to store the new sampling variable at the correct position make sure to define the index (SWI: SampWallIndex) as well. The index variable `SWIStickingCoefficient` is defined in the `MOD_Particle_Boundary_Vars` and can be later utilized to write and access the `SampWallState` array at the correct position, e.g.

    SampWallState(SWIStickingCoefficient,SubP,SubQ,SurfSideID) = SampWallState(SWIStickingCoefficient,SubP,SubQ,SurfSideID) + Prob

The calculation & output of the sampled values is performed in `CalcSurfaceValues` in `piclas/src/particles/dsmc/dsmc_analyze.f90` through the array `MacroSurfaceVal`. In a loop over all `nComputeNodeSurfSides` the sampled variables can be averaged (or manipulated in any other way). The variable `nVarCount` guarantees that you do not overwrite other variables

    IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
      nVarCount = nVarCount + 1
      IF(CounterSum.GT.0) MacroSurfaceVal(nVarCount,p,q,OutputCounter) = SampWallState(SWIStickingCoefficient,p,q,iSurfSide) / CounterSum
    END IF

Finally, the `WriteSurfSampleToHDF5` routine in `piclas/src/particles/boundary/particle_boundary_sampling.f90` writes the prepared `MacroSurfaceVal` array to the `ProjectName_DSMCSurfState_Timestamp.h5` file. Here, you have define a variable name, which will be shown in the output (e.g. in ParaView)

    IF (ANY(PartBound%SurfaceModel.EQ.1)) CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Sticking_Coefficient')

The order of the variable names and their position in the `MacroSurfaceVal` array has to be the same. Thus, make sure to place the `AddVarName` call at the same position, where you placed the calculation and writing into the `MacroSurfaceVal` array, otherwise the names and values will be mixed up.