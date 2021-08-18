\hypertarget{functions_and_subroutines}{}

# Useful Functions and Subroutines \label{chap:functions_and_subroutines}

This chapter contains a summary of useful functions and subroutines that might be re-used in the
future.

## General Functions and Subroutines

| Function     | Module        | Input          | Output     | Description                                                                                                |
| :----------: | :-----------: | :-----------:  | :--------: | :------------------------------------------------------------                                              |
| `UNITVECTOR` | `MOD_Globals` | 3D vector      | 3D vector  | Normalizes a given vector by dividing all vectors entries by the vector's magnitude                        |
| `CROSSNORM`  | `MOD_Globals` | two 3D vectors | 3D vector  | Computes the cross product of two 3-dimensional vectors: cross=v1 x v2 and normalizes the resulting vector |
| `CROSS`      | `MOD_Globals` | two 3D vectors | 3D vector  | Computes the cross product of two 3-dimensional vectors: cross=v1 x v2                                     |
| `VECNORM`    | `MOD_Globals` | 3D vector      | `REAL`     | Computes the Euclidean norm (length) of a vector                                                           |
| `DOTPRODUCT` | `MOD_Globals` | 3D vector      | `REAL`     | Computes the dot product of a vector with itself                                                           |

## Particle Functions and Subroutines

|                 Function (Module)                |                               Input                              |         Output         |                                                                Description                                                                |
|            :-------------------------:           |                      :---------------------:                     |      :----------:      |                                                     :--------------------------------                                                     |
|      `isChargedParticle` (`MOD_part_tools`)      |                            particle ID                           |        `LOGICAL`       |                                                Check if particle has charge unequal to zero                                               |
|         `PARTISELECTRON` (`MOD_globals`)         |                            particle ID                           |        `LOGICAL`       |           Check if particle is an electron by checking if the charge is equal to 1.602176634e-19 (division and nearest integer)           |
|      `isDepositParticle` (`MOD_part_tools`)      |                            particle ID                           |        `LOGICAL`       |                                              Check if particle is to be deposited on the grid                                             |
|        `isPushParticle` (`MOD_part_tools`)       |                            particle ID                           |        `LOGICAL`       |                                           Check if particle is to be pushed (integrated in time)                                          |
|    `isInterpolateParticle` (`MOD_part_tools`)    |                            particle ID                           |        `LOGICAL`       |                                   Check if the field at a particle's is to be interpolated (accelerated)                                  |
|     `VeloFromDistribution` (`MOD_part_tools`)    |                    distribution type, Tempergy                   |        3D vector       | **WIP**, Calculates a velocity vector from a defined velocity distribution and Tempergy (temperature [K] or energy [J] or velocity [m/s]) |
|        `DiceUnitVector` (`MOD_part_tools`)       |                              `None`                              |        3D vector       |                                   Calculates a normalized vector in 3D (unit space) in random direction                                   |
| `DiceDeflectedVelocityVector` (`MOD_part_tools`) |          cRela2(post-collison), alphaVSS(iSpecA,iSpecB)          |        3D vector       |                             Calculates scaled post-collision relative velocity vector in center-of-mass frame                             |
|                                                  | if alphaVSS>1 (VSS) also: , cRelaX,cRelaY,cRelaZ (pre-collision) |                        |                                 VSS case includes coordinate transformation due to anisotropic scattering                                 |
|        `CreateParticle` (`MOD_part_tools`)       | species ID, position, element ID, velocity and internal energies | particle ID (optional) |                  Creates a new particle at a given position and energetic state and return the new particle ID (optional)                 |
|      `GetParticleWeight` (`MOD_part_tools`)      |                            particle ID                           |         `REAL`         |                                               Determines the weighting factor of a particle                                               |
