// Tutorial: 3D 70degCone test case with gmsh
SetFactory("OpenCASCADE");
// Import provided 3D model as step file
v() = ShapeFromFile("70degCone_3D_model.step");
// Create a cylinder to represent the fluid domain
Cylinder(2) = {-50, 0, 0, 100, 0, 0, 50, Pi/6};
// Substract the 70degCone (tool) from the cylinder (object) = Object minus tool
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };
// Define boundary conditions
Physical Surface("IN", 29) = {4, 1};
//+
Physical Surface("SYM", 30) = {3, 5};
//+
Physical Surface("OUT", 31) = {2};
//+
Physical Surface("WALL", 32) = {7, 8, 9, 10, 11, 6};
//
Mesh.MeshSizeMin = 1;
//
Mesh.MeshSizeMax = 10;
//
Field[1] = MathEval;
//
Field[1].F = "0.2";
//
Field[2] = Restrict;
//
Field[2].SurfacesList = {7, 8, 9};
//
Background Field = 2;
//
Mesh.Algorithm = 1;
//
Mesh.Algorithm3D = 7;
//
Mesh.SubdivisionAlgorithm = 2;
//
Mesh.OptimizeNetgen = 1;
//
Mesh 3;
//
// Save all elements even if they are not part of a physical group, required to output volume elements
Mesh.SaveAll = 1;
// Save as ASCII format, Version 4
Mesh.Binary = 0;
Mesh.MshFileVersion = 4.1;
//
Save "70degCone_3D.msh";
