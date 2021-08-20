(sec:superB)=
# Magnetic Background Field (superB)

Certain application cases allow the utilization of a constant magnetic background field. The magnetic field resulting from certain
types of coils and permanent magnets can be calculated during the initialization within PICLas or with the standalone tool
**superB** (see Section {ref}`sec:compiler-options` for compilation), which can be used to solely create a .h5 file that contains
the B-field data via

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

## Magnetic Field by Permanent Magnets

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

## Magnetic Field by Coils

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

