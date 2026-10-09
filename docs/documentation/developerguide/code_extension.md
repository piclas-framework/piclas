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

## Allocating particle data

If an array is to store particle information, it can be allocated with

    ALLOCATE(ParticleInformation(1:PDM%maxParticleNumber))

Since PDM%maxParticleNumber is dynamic, the new array has to be added to the routines `IncreaseMaxParticleNumber` and `ReduceMaxParticleNumber` in `src/particles/particle_tools.f90`.

    IF(ALLOCATED(ParticleInformation)) CALL ChangeSizeArray(ParticleInformation,PDM%maxParticleNumber,NewSize, Default)

Default is an optional parameter if the new array memory is to be initialized with a specific value. The same must be done for TYPES of size PDM%maxParticleNumber. Please check both routines to see how to do it. Alternatively, the internal particle data structure can be extended as well (refer to {ref}`sec:partinten`), especially, if the particle information depends on the particle type (e.g. molecular species).

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

## Remove particles

If a particle leaves the simulation domain (e.g. through open boundaries), disappears
due to chemical reactions, or is removed for any other reason, it is **not sufficient**
to only set `PDM%ParticleInside = .FALSE.`. Instead, the subroutine

    CALL RemoveParticle(iPart)

from the module `MOD_part_operations` **MUST** be used. This routine ensures that all the flags and associated arrays are properly deallocated and that internal particle data structures remain consistent.

(sec:partinten)=
## Particle internal data container

Each particle owns an internal data container defined in the module `MOD_DSMC_Vars`.
This container is implemented as the derived type `tPartIntEn`, whose size is defined
by `PDM%maxParticleNumber`. All particle-related internal properties are attached to this type and are **only
allocated if required by the respective simulation model**.

    TYPE tPartIntEn
      REAL, ALLOCATABLE    :: EVib(:)        ! Vibrational energy
      REAL, ALLOCATABLE    :: ERot(:)        ! Rotational energy
      REAL, ALLOCATABLE    :: EElec(:)       ! Electronic energy
      REAL, ALLOCATABLE    :: TSolid(:)      ! Temperature of solid particles
      INTEGER, ALLOCATABLE :: QVib(:)        ! Vibrational quantum numbers
      INTEGER, ALLOCATABLE :: QRot(:)        ! Rotational quantum numbers
      INTEGER, ALLOCATABLE :: QElec(:)       ! Electronic quantum numbers
      REAL, ALLOCATABLE    :: DistriFunc(:)  ! Electronic distribution function
      REAL, ALLOCATABLE    :: ElecVelo(:)    ! Electron velocity for ambipolar diffusion
    END TYPE tPartIntEn

Not all particles carry all internal properties. For example, vibrational energy (`EVib`) is only allocated for molecular particles.
For such particles, allocation is performed as:

    ALLOCATE(PartIntEn(iPart)%EVib(1))

Vibrational quantum numbers (`QVib`) are not required in all physical models.
The property `TSolid` is only allocated for solid particles.
All other internal properties follow the same conditional allocation principle.
It is essential that all required properties are **allocated correctly during particle
creation**, depending on particle type and the active physical models.

When extending `tPartIntEn` with additional particle properties, the following routines must be extended as part of the particle management in `MOD_part_tools`:

- `ChangePartID`
- `ReduceMaxParticleNumber`
- `IncreaseMaxParticleNumber`
- `RemoveParticle`

Existing `PartIntEn` operators can be copied and used as templates.

### MPI communication

All new particle properties that must be communicated across MPI ranks require
explicit handling. First, the size of the communication array must be extended. The size of the 2D array `PartMPIExchange%nPartsSend(:,:)` is managed by a global variable named `nPartMPIData`, located in the module `MOD_Particle_MPI_Vars`. The array is allocated as follows:

```fortran
ALLOCATE(PartMPIExchange%nPartsSend(nPartMPIData, 0:nExchangeProcessors-1))
```

The variable `nPartMPIData` defines the total number of particle properties exchanged via MPI. To include an additional property, increment the value of `nPartMPIData` within the `MOD_Particle_MPI_Vars` module (e.g., change it from 7 to 8). The following components are updated automatically in `MOD_Particle_MPI`:

* `InitParticleCommSize`: Defines the size of communication data and allocates the required arrays
* `IRecvNbOfParticles`: Opens of the receive buffer for the number of particles to be received

However, the actual additional particle data has be added to the message manually in the following routines in `MOD_Particle_MPI`:

- `SendNbOfParticles`: Determine the size of MPI message
- `MPIParticleSend`: Build and send the MPI message
- `MPIParticleRecv`: Receive and unroll the MPI message

Each particle property is treated independently. Existing communication structures
for `PartIntEn` can be copied and adapted to the new property.

## Useful Functions

The following functions and subroutines are available for re-use.

### General Functions and Subroutines

|   Function   |    Module     |     Input      |  Output   | Description                                                                                                |
| :----------: | :-----------: | :------------: | :-------: | :--------------------------------------------------------------------------------------------------------- |
| `UNITVECTOR` | `MOD_Globals` |   3D vector    | 3D vector | Normalizes a given vector by dividing all vectors entries by the vector's magnitude                        |
| `CROSSNORM`  | `MOD_Globals` | two 3D vectors | 3D vector | Computes the cross product of two 3-dimensional vectors: cross=v1 x v2 and normalizes the resulting vector |
|   `CROSS`    | `MOD_Globals` | two 3D vectors | 3D vector | Computes the cross product of two 3-dimensional vectors: cross=v1 x v2                                     |
| `VECNORM3D`  | `MOD_Globals` |   3D vector    |  `REAL`   | Computes the Euclidean norm (length) of a vector                                                           |
| `DOTPRODUCT` | `MOD_Globals` |   3D vector    |  `REAL`   | Computes the dot product of a vector with itself                                                           |

### Particle Functions and Subroutines

|                Function (Module)                 | Input                                                            |         Output         | Description                                                                                                                               |
| :----------------------------------------------: | :--------------------------------------------------------------- | :--------------------: | :---------------------------------------------------------------------------------------------------------------------------------------- |
|      `isChargedParticle` (`MOD_part_tools`)      | particle ID                                                      |       `LOGICAL`        | Check if particle has charge unequal to zero                                                                                              |
|         `PARTISELECTRON` (`MOD_globals`)         | particle ID                                                      |       `LOGICAL`        | Check if particle is an electron by checking if the charge is equal to 1.602176634e-19 (division and nearest integer)                     |
|      `isDepositParticle` (`MOD_part_tools`)      | particle ID                                                      |       `LOGICAL`        | Check if particle is to be deposited on the grid                                                                                          |
|       `isPushParticle` (`MOD_part_tools`)        | particle ID                                                      |       `LOGICAL`        | Check if particle is to be pushed (integrated in time)                                                                                    |
|    `isInterpolateParticle` (`MOD_part_tools`)    | particle ID                                                      |       `LOGICAL`        | Check if the field at a particle's is to be interpolated (accelerated)                                                                    |
|    `VeloFromDistribution` (`MOD_part_tools`)     | distribution type, Tempergy                                      |       3D vector        | **WIP**, Calculates a velocity vector from a defined velocity distribution and Tempergy (temperature [K] or energy [J] or velocity [m/s]) |
|       `DiceUnitVector` (`MOD_part_tools`)        | `None`                                                           |       3D vector        | Calculates a normalized vector in 3D (unit space) in random direction                                                                     |
| `DiceDeflectedVelocityVector` (`MOD_part_tools`) | cRela2(post-collison), alphaVSS(iSpecA,iSpecB)                   |       3D vector        | Calculates scaled post-collision relative velocity vector in center-of-mass frame                                                         |
|                                                  | if alphaVSS>1 (VSS) also: , cRelaX,cRelaY,cRelaZ (pre-collision) |                        | VSS case includes coordinate transformation due to anisotropic scattering                                                                 |
|       `CreateParticle` (`MOD_part_tools`)        | species ID, position, element ID, velocity and internal energies | particle ID (optional) | Creates a new particle at a given position and energetic state and return the new particle ID (optional)                                  |
|       `RemoveParticle` (`MOD_part_tools`)        | particle ID                                                      |                        | Removes the specified particle by setting all the required flags and deallocating internal data arrays                                    |
|      `GetParticleWeight` (`MOD_part_tools`)      | particle ID                                                      |         `REAL`         | Determines the weighting factor of a particle                                                                                             |

## CollectiveStop
When using the `CALL abort(__STAMP__,'ERROR ...')` subroutine, each MPI process that encounters this
call emits an output to std.out, which can result in a huge amount of output when running on hundreds
or thousands of processes, especially if an error in the .ini files produces the error in the initialization step.
To prevent excessive output to std.out, the `CollectiveStop` subroutine has been created.
```
!==================================================================================================================================
!> \brief Safely terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!>
!> Safely terminate program using a soft MPI_FINALIZE in the MPI case and writes the error message only on the root.
!> Terminate program using a soft MPI_FINALIZE in the MPI case and write the error message only on the root.
!> This routine can only be used if ALL processes are guaranteed to generate the same error at the same time!
!> Prime use is to exit PICLAS without MPI errors and with a single error message if some parameters are not set in the init
!> routines or a file is not found.
!>
!> Criteria where CollectiveStop may be used:
!> 0. In case of doubt stick with Abort, which is always safe!
!> 1. A routine is BY DESIGN (!) called by all processes, i.e. does not permit to be called by single processes or subgroups.
!> 2. The criteria for the CollectiveStop must be identical among all processors.
!> 3. The routine is only used during the init phase.
!> 4. The error must not originate from MPI errors (e.g. during MPI init)
!> 5. The error must not originate from checking roundof errors (e.g. accuracy of interpolation matrices)
!>
!==================================================================================================================================
```
An example where `CollectiveStop` should be used instead of `abort` is in `gradients.f90`

```
GradLimiterType=GETINT('Grad-LimiterType')
GradLimVktK=GETREAL('Grad-VktK')
SELECT CASE(GradLimiterType)
CASE(0)
  LBWRITE(UNIT_stdOut,*)'Limiter = 0 -> first order FV'
CASE(1) !minmax
  LBWRITE(UNIT_stdOut,*)'Using Barth-Jespersen Limiter'
CASE(4) !venkatakrishnan
  LBWRITE(UNIT_stdOut,*)'Using Venkatakrishnan limiter with K =', GradLimVktK
CASE(9) ! no limiter (central)
  LBWRITE(UNIT_stdOut,*)'Not using any limiter'
CASE DEFAULT
  CALL abort(__STAMP__,'Limiter type not implemented.')
END SELECT
```
which should read
```
GradLimiterType=GETINT('Grad-LimiterType')
GradLimVktK=GETREAL('Grad-VktK')
SELECT CASE(GradLimiterType)
CASE(0)
  LBWRITE(UNIT_stdOut,*)'Limiter = 0 -> first order FV'
CASE(1) !minmax
  LBWRITE(UNIT_stdOut,*)'Using Barth-Jespersen Limiter'
CASE(4) !venkatakrishnan
  LBWRITE(UNIT_stdOut,*)'Using Venkatakrishnan limiter with K =', GradLimVktK
CASE(9) ! no limiter (central)
  LBWRITE(UNIT_stdOut,*)'Not using any limiter'
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Limiter type not implemented.')
END SELECT
```
