# Magnetic Background Field

Certain application cases allow the utilization of a constant magnetic background field.
This background field can either be supplied via .csv (1D) or .h5 (2D) file, however, in this case the field must be based on an equidistant
Cartesian mesh.
Another method is to use the built-in tool **superB**, which is also available as stand-alone executable to generate magnetic fields
based on magnets or coils, which results in the creation of a .h5 file containing the field data based on a PICLas (HOPR) mesh file.
The following two sections give an overview of using the different methods.

(sec:variableExternalField)=
## Variable External Field

One-, two- and three-dimensional magnetic fields can be used as fixed background fields for certain time discretization methods
(full Maxwell time discs and the Poisson Boris-Leapfrog scheme)
The read-in variable for either .csv (only for 1D along $z$-direction) or .h5 (only for 2D axis symmetric $z$-direction or fully 3D)
files is set via

    PIC-variableExternalField = X.csv, X.h5

Three examples are located within the regression test directory

    regressioncheck/CHE_PIC_maxwell_RK4/gyrotron_variable_Bz
    regressioncheck/CHE_PIC_maxwell_RK4/2D_variable_B
    regressioncheck/CHE_PIC_maxwell_RK4/3D_variable_B

for 1D, 2D and 3D fields, respectively. Note that 1D currently only allows magnetic fields of type $B_{z}(z)$ and 2D only allows the 
components $B_{r}(r,z)$ and $B_{z}(r,z)$ that comprise a rotationally symmetric vector field $\textbf{B}$.

The first example (1D and .csv file) uses data via

    PIC-variableExternalField = variable_Bz.csv

which is csv-based data in the form (the delimiter is actually not a comma)

    -0.00132 	2.7246060625
    -0.000217551020408	2.700481016
    0.0008848979591837	2.6762685135
    0.0019873469387755	2.6519260266
    0.0030897959183674	2.6274128336
    ....

and the second (2D and .h5 file)

    PIC-variableExternalField = reggie-linear-rot-symmetry.h5 

that is read from a .h5 file.
The data structure in the .h5 file must be of the form (dataset is labelled "data") and contains

    r1 z1 Br1 Bz1
    r2 z2 Br2 Bz2
    r3 z3 Br3 Bz3
    r4 z4 Br4 Bz4
    r5 z5 Br5 Bz5
    ....

where for each data point one row is required.
The ordering of the data is also important.
It is only allowed that the first $N$ rows have the same $r$ value and varying $z$-components ($r$ is the outer loop variable and
$z$ is the inner loop variable when unrolling the data into an array).
This is automatically checked by comparing the distances in $r$ and $z$ direction, which must be equidistant.
In addition, the attributes r, z, Br and Bz, which contain the indices of the corresponding column number in "data".

Three-dimensional fields must be supplied in the following format

    x1 y1 z1 Bx1 By1 Bz1
    x2 y2 z2 Bx2 By2 Bz2
    x3 y3 z3 Bx3 By3 Bz3
    x4 y4 z4 Bx4 By4 Bz4
    x5 y5 z5 Bx5 By5 Bz5
    ....

where the data (dataset is labelled "data" in the .h5 file) is sorted in lines in ascending coordinates.
For everything to work, the order must be like this

    x1 y1 z1 Bx1 By1 Bz1
    x2 y1 z1 Bx2 By2 Bz2
    x1 y2 z1 Bx3 By3 Bz3
    x2 y2 z1 Bx4 By4 Bz4
    x1 y1 z2 Bx5 By5 Bz5
    x2 y1 z2 Bx6 By6 Bz6
    ....

where first the $x$-coordinate changes, then $y$ and finally the $z$.

(sec:superB)=
## superB
The magnetic field resulting from certain types of coils and permanent magnets can be calculated during the initialization within 
PICLas or with the standalone tool **superB** (see Section {ref}`sec:compiler-options` for compilation), which can be used to solely
create a .h5 file that contains the B-field data via

    superB parameter_superB.ini

For usage in PICLas, the background field can be enabled by

    PIC-BG-Field = T

The first option is to use a previously calculated background field. It can be read-in with

    PIC-BGFileName     = BField.h5 ! Path to a .h5 file that contains the B-field data
    PIC-NBG            = 1         ! Polynomial degree of the B-field
    PIC-BGFieldScaling = 1.        ! Scaling factor for the B-field

Additionally, the polynomial degree for the background field can be set by ``PIC-NBG`` and might differ from the actually read-in
polynomial degree that is used to represent the numerical solution for the field solver. Optionally, the read-in field can be
scaled by the last of the three parameters above.

The second option is to calculate the magnetic field during the initialization, which will produce an output .h5 file of the field.
The field will automatically be calculated from the supplied parameters, if the corresponding .h5 file does not exist.
For this purpose, different coil and permanent magnet geometries can be defined. For visualization purposes, the geometry of the
respective coils and permanent magnets can be directly written out as a VTK with

    PIC-CalcBField-OutputVTK = T

In the following the parameters for different coils and permanent magnets based on the implementation by Hinsberger {cite}`Hinsberger2017` are presented.

### Magnetic Field by Permanent Magnets

First, the total number of permanent magnets has to be defined and the type selected. Options are `cuboid`, `sphere`, `cylinder` and `conic`.

    NumOfPermanentMagnets = 1
    PermanentMagnet1-Type = cuboid
                            sphere
                            cylinder ! also used for hollow cylinders
                            conic

All options require the input of a base/origin vector, a number of discretization nodes (results in a different number of total
points depending on the chosen geometry) and a magnetisation in [A/m]

    PermanentMagnet1-BasePoint = (/0.,0.,0./)
    PermanentMagnet1-NumNodes = 10
    PermanentMagnet1-Magnetisation = (/0.,0.,1./)

The geometries require different input parameters given below

    ! Three vectors spanning the cuboid
    PermanentMagnet1-BaseVector1 = (/1.,0.,0./)
    PermanentMagnet1-BaseVector2 = (/0.,1.,0./)
    PermanentMagnet1-BaseVector3 = (/0.,0.,1./)
    ! Radius required for a spherical, cylindrical and conical magnet
    PermanentMagnet1-Radius = 1.
    ! Height vector required for a cylindrical and conical magnet
    PermanentMagnet1-HeightVector = (/0.,0.,1./)
    ! Second radius only required for a conical magnet or a hollow cylinder with inner radius
    ! 'Radius2' and outer radius "Radius1'
    PermanentMagnet1-Radius2 = 1.

### Magnetic Field by Coils

The total number of coils and the respective type of the cross-section (`linear`,`circle`,`rectangle`,`custom`) is defined by

    NumOfCoils = 1
    Coil1-Type = linear
                 circle
                 rectangle
                 custom

All options require the input of a base/origin vector, a length vector (vector normal to the cross-section of the coil) and the
current in [A]

    Coil1-BasePoint = (/0.0,0.0,0.0/)
    Coil1-LengthVector = (/0.0,1.0,0.0/)
    Coil1-Current = 1.

The first option `linear` represents a simple linear conductor (e.g. a straight wire) and requires only the input of a number of
discretization points

    Coil1-NumNodes = 5

The other three types, which are actually coils, are described by the number of loops and the number of discretization points per loop

    Coil1-LoopNum = 10
    Coil1-PointsPerLoop = 10

The cross-section of the coil is defined in a plane normal to the `-LengthVector`. A circular coil cross-section requires simply
the input of a radius whereas a rectangular coil cross-section is spanned by two vectors (`-RectVec1` and `-RectVec2`) and an
additional vector, which must be orthogonal to the `-LengthVector` to define the orientation of the cross-section (`-AxisVec1`).
In these two cases, the base/origin vector defines the middle point of the cross-section.

    ! Circular coil cross-section
    Coil1-Radius = 1.
    ! Rectangular coil cross-section
    Coil1-RectVec1 = (/1.0,0.0/)
    Coil1-RectVec2 = (/0.0,1.0/)
    Coil1-AxisVec1 = (/0.0,0.0,1.0/)

The last cross-section type `custom` allows the definition of a cross-section as a combination of multiple linear (`line`) and
circular (`circle`) segments and also requires an additional vector to define the orientation of the cross-section (`-AxisVec1`)

    Coil1-AxisVec1 = (/0.0,0.0,1.0/)
    Coil1-NumOfSegments = 3
    ! Linear segment defined by
    Coil1-Segment1-SegmentType = line
    Coil1-Segment1-NumOfPoints = 5
    Coil1-Segment1-LineVector = (/1.0,1.0/)
    ! Circular segment connected to the previous segment
    Coil1-Segment2-SegmentType = circle
    Coil1-Segment2-NumOfPoints = 5
    Coil1-Segment2-Radius = 1.
    Coil1-Segment2-Phi1 = 90.
    Coil1-Segment2-Phi2 = 0.
    ! Linear segment connected to the previous segment, closing the cross-section
    Coil1-Segment3-SegmentType = line
    Coil1-Segment3-NumOfPoints = 5
    Coil1-Segment3-LineVector = (/-2.0,0.0/)

The `-NumOfPoints` controls the number of discretization points per segment. A linear segment is simply described by a vector in
the cross-section plane. The circular segment is defined by a radius and the initial as well as final angle of the segment.
It should be noted that the base point defines the start of the first segment as opposed to the circular and rectangular
cross-sections, where it is the middle point of the cross-section.

### Time-dependent Magnetic Coils
A time-dependent magnetic field can be created by a time-varying electric current running through a coil.
Time-dependent coils can be combined with an arbitrary number of permanent magnets and coils (with a constant current).
Currently, all time-dependent coils must use the same frequency but can have different phases.
The time-dependent settings are required in addition to the ones used for a standard coil

    Coil1-TimeDepCoil      = T
    Coil1-CurrentFrequency = 1e6
    Coil1-CurrentPhase     = 0.0
    nTimePoints            = 11

where the frequency and phase of the sin function that is used for the electric current as well as the number of points in time for
the interpolation of the current is required. In piclas, times between two interpolation points are determined by linear
interpolation from the stored solution.
