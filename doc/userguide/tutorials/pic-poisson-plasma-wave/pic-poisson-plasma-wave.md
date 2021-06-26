## PIC Plasma Wave (Poisson's Equation)
\label{sec:tut_freestream}

The setup considers a 1D plasma oscillation, which is a common and simple electrostatic PIC benchmark [@Birdsall1991], [@Hockney1988],[@Jacobs2006b].
In PICLas it can be simulated either with the full Maxwell solver (DGSEM) or with the Poisson solver (HDGSEM), where the latter is
chosen for this this tutorial.
In this setup, electrons oscillate around the almost immobile ions, which creates a fluctuating electric field


Copy the `pic-poisson-plasma-wave` directory from the tutorial folder in the top level directory to a separate location

        cp -r $PICLAS_PATH/tutorials/pic-poisson-plasma-wave .
        cd pic-poisson-plasma-wave




### Mesh Generation with HOPR (pre-processing)

Before the actual simulation is conducted, a mesh file in the correct HDF5 format has to be supplied.
The mesh files used by **PICLas** are created by supplying an input file *hopr.ini* with the required information for a mesh that
has either been created by an external mesh generator or directly from block-structured information in the hopr.ini file itself.
Here, a block-structured grid is created directly from the information in the hopr.ini file.
To create the .h5 mesh file, simply run

    hopr hopr.ini

This creates the mesh file *plasma_wave_mesh.h5* in HDF5 format.
The size of the simulation domain is set to [$2\pi\times0.2\times0.2$] m$^{3}$ and is defined by the single block information in the line

    Corner         =   (/0.,0.,0.,,6.2831,0.,0.,,6.2831, ...........

The number of mesh elements in each direction can be adjusted by changing the line

    nElems         = (/60,1,1/)                ! number of elements in each directio

Alternatively, if you do not want to run **hopr** yourself, you can also use the provided mesh.


Contrary to the particle boundary conditions other than periodic, the field boundaries can directly be defined in the hopr.ini file

    !=============================================================================== !
    ! BOUNDARY CONDITIONS
    !=============================================================================== !
    nUserDefinedBoundaries=6
    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)
    BoundaryName = BC_periodicy+ ! Periodic (+vv2)
    BoundaryType = (/1,0,0,2/)   ! Periodic (+vv2)
    BoundaryName = BC_periodicy- ! Periodic (-vv2)
    BoundaryType = (/1,0,0,-2/)  ! Periodic (-vv2)
    BoundaryName = BC_periodicz+ ! Periodic (+vv3)
    BoundaryType = (/1,0,0,3/)   ! Periodic (+vv3)
    BoundaryName = BC_periodicz- ! Periodic (-vv3)
    BoundaryType = (/1,0,0,-3/)  ! Periodic (-vv3)

    nVV=3                     ! Number of periodic displacement vectors for each periodic BC
    VV=(/6.2831 , 0.  , 0./)  ! Displacement vector 1 (x-direction)
    VV=(/0.     , 0.2 , 0./)  ! Displacement vector 2 (y-direction)
    VV=(/0.     , 0.  , 0.2/) ! Displacement vector 3 (z-direction)

and in this case a fully periodic setup is chosen by defining periodic boundaries on all six sides of the block, reflecting each
positive and negative Cartesian coordinate. In x-direction,

    BoundaryName = BC_periodicx+ ! Periodic (+vv1)
    BoundaryType = (/1,0,0,1/)   ! Periodic (+vv1)
    BoundaryName = BC_periodicx- ! Periodic (-vv1)
    BoundaryType = (/1,0,0,-1/)  ! Periodic (-vv1)

where for each of the six boundaries, a name `BoundaryNamethe` and a type `BoundaryType` must be defined (in this order)
The first "1" in `BoundaryType` corresponds to the type "periodic" and the last entry, here, either "1" or "-1" corresponds to the
first periodic vector that is defined via `VV=(/6.2831 , 0.  , 0./)` that handles periodicity in the x-direction and gives the
orientation on the boundary for the vector.
Note that each periodic boundary must have one positive and one negative corresponding boundary for the same periodic vector.




![](tutorials/pic-poisson-plasma-wave/mesh/pic.pdf)

Figure: Mesh with $60\times1\times1$ elements and a size of [$2\pi\times0.2\times0.2$] m$^{3}$.




### PIC Simulation with PICLas

The simulation setup is defined in *parameter_piclas.ini*. The initial condition is selected via the variable vector
**RefState=(/1.225,1.,1.,1.,101325./)** which represents the solution vector $(\rho, u, v, w, p)^T$. **PICLas** allows for multiple
**RefState** vectors and numerates them respectively for them to be used by different functions. In this example a single **RefState** is
supplied and therefore is given the number **1**.


**IniRefState = 1** : the initial condition uses **RefState 1** for the initial flow field solution.

**IniExactFunc = 1** : the used exact function routine uses **RefState 1**, e.g., for the calculation of the $L_2$ error norm.

Constant flow properties like the gas constant are given in table \ref{tab:freestream_flow_prop}
and define the gas behavior in combination with the ideal gas law $p=\rho R T$.

Table: Numerical settings \label{tab:freestream_flow_prop}

| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           | 0.000018547 |
| ideal gas constant $R$          | R             |  276        |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |



#### Numerical Setup

\label{sec:tut_freestream_num_settings}

The DG solution on the mesh is represented by piecewise polynomials and the polynomial degree in this tutorial is chosen as $N=3$.

The default settings for these properties are displayed in table \ref{tab:freestream_num_set}.

Table: Numerical settings \label{tab:freestream_num_set}

| Variable        | Description                            | Value         |
| --------------- |:--------------------------------------:|:-------------:|
| N               | polynomial degree                      | 3             |
| MeshFile        |                                        |cartbox_mesh.h5|
| tend            |                                        | 1e-6          |
| Analyze_dt      |                                        | 1e-6          |
| nWriteData      |                                        | 1             |
| CFLscale        |                                        | 0.99          |
| DFLscale        |                                        | 0.4           |


The command

~~~~~~~
piclas parameter_piclas.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*.
If the run has completed successfully, which should take only a brief moment, the contents of the working folder should look like in figure \ref{fig:freestream_folder}

<!--![The folder contents after a successful run\label{fig:freestream_folder}](tutorials/00_freestream/freestream_folder.png)-->

Two additional files have been created, which are are named  **Projectname_State_Timestamp.h5**. They contain the solution vector of the conserved variables at each interpolation node at the given time, which corresponds to multiplies of **Analyze_dt**. If these files are not present, something went wrong during the execution of **PICLas**. In that case, check the _std.out_ file for an error message.

After a successful completion, the last lines in this files should look like in figure \ref{fig:freestream_stdout}

<!--![The _std.out_ file after a successful run\label{fig:freestream_stdout}](tutorials/00_freestream/freestream_stdout.png)-->



### Visualization (post-processing)

To visualize the solution, the *State*-files must be converted into a format suitable for **ParaView**. Issue the command

~~~~~~~
posti_visu parameter_postiVisu.ini parameter_piclas.ini cartbox_State_0000000.00000*
~~~~~~~
to generate the corresponding *vtu*-files, which can then be loaded into **ParaView**.
