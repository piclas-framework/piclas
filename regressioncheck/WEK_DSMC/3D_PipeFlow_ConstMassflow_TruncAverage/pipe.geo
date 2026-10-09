// ============================================================================
//  pipe.geo  -  3D hexahedral mesh of a cylindrical pipe
//
//  Axis      : x  (x_min = 0 ... x_max = 5)
//  Cross-sec : circle in the y-z plane, radius R = 1/sqrt(pi)
//  Topology  : O-grid / butterfly  ->  pure hexahedra
//  BC groups : inlet (x_min), outlet (x_max), wall (curved surface)
//
//  Generate the mesh with:   gmsh -3 pipe.geo -o pipe.msh
// ============================================================================

// ---------------------------------------------------------------------------
//  Geometry parameters
// ---------------------------------------------------------------------------
R      = 0.564189583547756;   // pipe radius (= 1/sqrt(pi), cross-section area = 1)
L      = 5.0;                 // pipe length along x
fInner = 0.5;                 // inner (core) square corner radius / R

// ---------------------------------------------------------------------------
//  Mesh resolution  (Medium)
// ---------------------------------------------------------------------------
nCirc  = 5;    // nodes along each core-square edge / arc  (= 8 cells)
nRad   = 3;    // nodes in the radial O-grid ring          (= 4 cells)
nLen   = 10;   // extrusion layers along x                 (= 20 cells)

// ---------------------------------------------------------------------------
//  Derived point coordinates (cross-section lives at x = 0, y-z plane)
// ---------------------------------------------------------------------------
s  = R / Sqrt(2.0);            // outer circle points on the 45-deg diagonals
si = fInner * R / Sqrt(2.0);   // inner core-square corners

// Center
Point(1) = {0,  0.0,  0.0};

// Outer ring points on the circle (45,135,225,315 deg)
Point(2) = {0,  s,  s};
Point(3) = {0, -s,  s};
Point(4) = {0, -s, -s};
Point(5) = {0,  s, -s};

// Inner core-square corners
Point(6) = {0,  si,  si};
Point(7) = {0, -si,  si};
Point(8) = {0, -si, -si};
Point(9) = {0,  si, -si};

// ---------------------------------------------------------------------------
//  Lines
// ---------------------------------------------------------------------------
// Inner core-square edges
Line(1) = {6, 7};
Line(2) = {7, 8};
Line(3) = {8, 9};
Line(4) = {9, 6};

// Outer circular arcs (each spans 90 deg, center = point 1)
Circle(5) = {2, 1, 3};
Circle(6) = {3, 1, 4};
Circle(7) = {4, 1, 5};
Circle(8) = {5, 1, 2};

// Radial connectors (inner corner -> outer ring point)
Line( 9) = {6, 2};
Line(10) = {7, 3};
Line(11) = {8, 4};
Line(12) = {9, 5};

// ---------------------------------------------------------------------------
//  Surfaces  (1 central quad + 4 outer ring quads)
// ---------------------------------------------------------------------------
Curve Loop(1) = {1, 2, 3, 4};            Plane Surface(1) = {1};   // core
Curve Loop(2) = {1, 10, -5,  -9};        Plane Surface(2) = {2};   // ring (arc 5)
Curve Loop(3) = {2, 11, -6, -10};        Plane Surface(3) = {3};   // ring (arc 6)
Curve Loop(4) = {3, 12, -7, -11};        Plane Surface(4) = {4};   // ring (arc 7)
Curve Loop(5) = {4,  9, -8, -12};        Plane Surface(5) = {5};   // ring (arc 8)

// ---------------------------------------------------------------------------
//  Structured (transfinite) meshing of the cross-section
// ---------------------------------------------------------------------------
Transfinite Curve {1, 2, 3, 4}     = nCirc;   // core-square edges
Transfinite Curve {5, 6, 7, 8}     = nCirc;   // arcs
Transfinite Curve {9, 10, 11, 12}  = nRad;    // radial lines

Transfinite Surface {1, 2, 3, 4, 5};
Recombine Surface   {1, 2, 3, 4, 5};          // quads -> ready for hex extrusion

// ---------------------------------------------------------------------------
//  Extrude the cross-section along x to obtain hexahedra
//  Each Extrude returns: out[0] = top surface, out[1] = volume,
//                        out[2+k] = side surface of the k-th loop curve.
// ---------------------------------------------------------------------------
core[] = Extrude {L, 0, 0} { Surface{1}; Layers{nLen}; Recombine; };
r1[]   = Extrude {L, 0, 0} { Surface{2}; Layers{nLen}; Recombine; };
r2[]   = Extrude {L, 0, 0} { Surface{3}; Layers{nLen}; Recombine; };
r3[]   = Extrude {L, 0, 0} { Surface{4}; Layers{nLen}; Recombine; };
r4[]   = Extrude {L, 0, 0} { Surface{5}; Layers{nLen}; Recombine; };

// For the ring surfaces the arc is the 3rd curve in the loop (k = 2),
// so the curved wall side surface is out[2 + 2] = out[4].

// ---------------------------------------------------------------------------
//  Physical groups (named) for HOPR / PICLas
// ---------------------------------------------------------------------------
Physical Surface("BC_xminus")  = {1, 2, 3, 4, 5};                       // x_min
Physical Surface("BC_xplus") = {core[0], r1[0], r2[0], r3[0], r4[0]}; // x_max
Physical Surface("BC_shell")   = {r1[4], r2[4], r3[4], r4[4]};          // curved surface
