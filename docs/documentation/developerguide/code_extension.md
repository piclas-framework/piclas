# Code Extension

This section shall describe how to extend the code in a general way e.g. to implement new input and output parameters.

## Surface Sampling & Output

*Location: `piclas/src/particles/boundary/particle_boundary_sampling.f90`*

The surface sampling values are stored in the `SampWallState` array and the derived output variables in `MacroSurfaceVal`, which can be extended to include new **optional** variables and to exploit the already implemented MPI communication. The variables `SurfSampSize` and `SurfOutputSize` define the size of the arrays is set in `InitParticleBoundarySampling` with default values as given in `piclas.h`

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

The calculation & output of the sampled values is performed in `CalcSurfaceValues` through the array `MacroSurfaceVal`. In a loop over all `nComputeNodeSurfSides` the sampled variables can be averaged (or manipulated in any other way). The variable `nVarCount` guarantees that you do not overwrite other variables

    IF(ANY(PartBound%SurfaceModel.EQ.1)) THEN
      nVarCount = nVarCount + 1
      IF(CounterSum.GT.0) MacroSurfaceVal(nVarCount,p,q,OutputCounter) = SampWallState(SWIStickingCoefficient,p,q,iSurfSide) / CounterSum
    END IF

Finally, the `WriteSurfSampleToHDF5` routine writes the prepared `MacroSurfaceVal` array to the `ProjectName_DSMCSurfState_Timestamp.h5` file. Here, you have define a variable name, which will be shown in the output (e.g. in ParaView)

    IF (ANY(PartBound%SurfaceModel.EQ.1)) CALL AddVarName(Str2DVarNames,nVar2D_Total,nVarCount,'Sticking_Coefficient')

The order of the variable names and their position in the `MacroSurfaceVal` array has to be the same. Thus, make sure to place the `AddVarName` call at the same position, where you placed the calculation and writing into the `MacroSurfaceVal` array, otherwise the names and values will be mixed up.

## Arrays with size of PDM%maxParticleNumber

If an array is to store particle information, it is usually allocated with

    ALLOCATE(ParticleInformation(1:PDM%maxParticleNumber))

But since PDM%maxParticleNumber is dynamic, this behavior must also be captured by this array. Therefore its size has to be changed in the routines `IncreaseMaxParticleNumber` and `ReduceMaxParticleNumber` in `src/particles/particle_tools.f90`.

    IF(ALLOCATED(ParticleInformation)) CALL ChangeSizeArray(ParticleInformation,PDM%maxParticleNumber,NewSize, Default)

Default is an optional parameter if the new array memory is to be initialized with a specific value. The same must be done for TYPES of size PDM%maxParticleNumber. Please check both routines to see how to do it.

## Insert new particles

To add new particles, first create a new particle ID using the GetNextFreePosition function contained in `src/particles/particle_tools.f90`

    NewParticleID = GetNextFreePosition()

This directly increments the variable PDM%CurrentNextFreePosition by 1 and if necessary adjusts PDM%ParticleVecLength by 1. If this is not desired, it is possible to pass an offset. Then the two variables will not be incremented, which must be done later by the developer. This can happen if the particle generation process is divided into several functions, where each function contains a loop over all new particles (e.g. `src/particles/emission/particle_emission.f90`).

    DO iPart=1,nPart
        NewParticleID = GetNextFreePosition(iPart)
    END DO
    PDM%CurrentNextFreePosition = PDM%CurrentNextFreePosition + nPart
    PDM%ParticleVecLength = MAX(PDM%ParticleVecLength,GetNextFreePosition(0))

For the new particle to become a valid particle, the inside flag must be set to true and various other arrays must be filled with meaningful data. See SUBROUTINE CreateParticle in `src/particles/particle_operations.f90`. A basic example of the most important variables is given below:

    newParticleID = GetNextFreePosition()
    PDM%ParticleInside(newParticleID) = .TRUE.
    PDM%FracPush(newParticleID) = .FALSE.
    PDM%IsNewPart(newParticleID) = .TRUE.
    PEM%GlobalElemID(newParticleID) = GlobElemID
    PEM%LastGlobalElemID(newParticleID) = GlobElemID
    PartSpecies(newParticleID) = SpecID
    LastPartPos(1:3,newParticleID) = Pos(1:3)
    PartState(1:3,newParticleID) = Pos(1:3)
    PartState(4:6,newParticleID) = Velocity(1:3)