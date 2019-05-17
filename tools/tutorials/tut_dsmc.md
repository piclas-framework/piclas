This section describes selected examples as tutorials to simplify the first usage of the parameter.ini files.
1. 'Cube-Visu' covers the difference between the micro and the macro level sampling.
2. An analyze and sampling example will be added soon.
3. The 'cuboid with diffusion and two initial values' deals with the use of different initial values throughout the mesh. 
4. The 'quasi-2D tangential flow' demonstrates the different inflow mechanisms using surfaceflux or volume emission.
5. The scattering surface creates an insight into the interaction with the walls with TransACC,VibACC,RotACC, and MomentACC. 
6. The cone 
7. 70° arrow


1. #cubeVisu
##The hopr.ini
The mesh for the simulation is determined by the hopr.ini. 
This mesh consists of two cells with perfectly reflecting boundaries.
To build the mesh files needed for the calculation, execute

*[path_to_hopr]/hopr* *./hopr.ini*

As a result, a file named *cube_mesh.h5* appears in the folder. 
##The parameter.ini 
The parameter.ini contains all the parameters needed for the calculation. To get neat ini files, one can split it up and create a seperate DSMC.ini file containing the detailed species information.
For many parameters a default value will be set, in case no other input is provided. 

Execute 

*[path]/piclas --help >> help.txt* 

to create a file containing information about the input parameters, their modes and default values.

In this tutorial the focus lies on the comprehension the difference between the macro values of the simulation and the micro output.
 
Part-WriteMacroValues     = T
Part-IterationForMacroVal = 1 
As a result, a *HOState* *SurfState* per iteration is created.

The *[ProjectName]_DSMCHOState_[Timestep].h5* contains the macroscopic high order information of of the simulation. (high order...)
The *[ProjectName]_DSMCSurfState_[Timestep].h5* contains the macroscopic surface information of of the simulation.

The *[ProjectName]_State_[Timestep].h5* contains the microscopic values of each simulated particle and is therefore larger. 


##Run PICLas 
To run the simulation you need to execute PIClas from whereever you have localized it on your computer.
The console input looks similar to 
*[path_to_piclas]/piclas [path_to_cubeVisu]/parameter.ini [path_to_cubeVisu]/DSMC.ini*

or in case you put all the information into the parameter.ini the DSMC.ini drops.

*[path_to_piclas]/piclas [path_to_cubeVisu]/parameter.ini*

When the simulation is successfully finished you want to visualize the results.

##Visualization
To do so, you need to convert the output files with their *h5* format to *vtk* or *vtu* files, respectively.

*[path_to_paraviewposti]/[build_folder]/bin/visuPart [path_to_cubeVisu]/CubeVisu_State_\* 

creates *CubeVisu_visPart_[TimeStamp].vtk* files.

In paraview you can open the group of state files to use the time tool later on. 
Now you are able to view the data which includes the velocity of each particle per axis, rotational and vibrational information.
The particles are visible, however the mesh is not.
To view the mesh in each time step you might execute the visu3D tool respectively.

*[path_to_paraviewposti]/[build_folder]/bin/visu3D [path_to_cubeVisu]/CubeVisu_State_\* 

and receive *CubeVisu_visu3D_[TimeStamp].pvtu* files which include information about the electric field, magnetic field and other data of the mesh.
Opening the file group in paraview the mesh is displayed.

##Change parameters
By altering 'Part-Species1-SpaceIC' the method of inserting particles is changed. For the initally shown option 'cuboid', some extra parameters, as shown in the file, 
have to be defined.

In this example, try the different inserting methods to get an understanding what the differences are by viewing the output. Change 'cuboid' to 'cell_local', where no further 
variables must be set. Therefore, comment the cuboid specific parameters out.

When the calculation runs through and you have visualized the output successfully, try 'cylinder' as 'SpaceIC'. This option needs the following additional parameters


Part-Species1-RadiusIC              = 5.0E-4             ! Required for cylinder 
Part-Species1-Radius2IC             = 4.0E-4             ! Required for cylinder 
Part-Species1-NormalIC              = (/0.,0.,1./)       ! Required for cylinder 
Part-Species1-CylinderHeightIC      = 1e-4               ! Required for cylinder

Try out some other inserting methods by checking the list of 'SpaceIC' options in the earlier created *help.txt* 


For example: 
SpaceIC=cuboid_vpi will lead to an abort, if the ParticleEmissionType is neither emission rate in part/s (=1) nor part/iteration (=2)
##Why is the Code aborting?
In case your calculation aborts immediately,
you might
+ Check if all the variables needed for the SpaceIC are chosen. Should the help.txt file not help to fix the problem, you may look in [piclas]/src/particles/particle_init.f90 for "SpaceIC".
  This way it is possible to determine the exact reason for the program to abort and adapt the respective parameters.

+ Make sure the chosen name of the mesh file, defined within the hopr.ini, is the name written in the parameter.ini  as 'MeshFile'

##Good to know!

+ As you may know, Fortran is a case-insensitive programming language, which is why the parameter.ini file is case-insensitive, too.

+ The balance between calculating and writing output files with the processing power needs to be adapted for bigger simulations 

+ The Part-FIBGMdeltas for the background mesh should be set between the smallest and meanest cell dimensions.
  To know which value to choose, check the hopr.ini, where 'corner' defines the dimension of one cell and 'nElems' creates two identical cells.
#2. multipleInit
Use more than one initial value for the same species. Choose between surfaceflux or adaptive 
## Change the Insertion Method
 
#3. Quasi 2D case of application
In this example the nInit is greater 1.

PICTURE SURFACEFLUX
VOLUMEINSERTING

cylinder Surface
Choose between Surfaceflux or Volume inserting

#4. 
#5. 
#6. 
#7. 
