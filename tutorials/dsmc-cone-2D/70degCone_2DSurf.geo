// Gmsh project created on Tue Oct 18 18:46:36 2022
SetFactory("OpenCASCADE");
v() = ShapeFromFile("70degCone_2DSurf.step");
//
Translate {0, 0, -1.0} {
  Surface{1}; 
}
//
Physical Curve("SYMAXIS") = {3};
//
Physical Curve("WALL") = {1, 2};
//
Physical Curve("IN") = {4};
//
Physical Curve("OUT") = {5, 6};
//
Mesh.MeshSizeMin = 0;
//
Mesh.MeshSizeMax = 10;
//
Field[1] = MathEval;
//
Field[1].F = "0.2";
//
Field[2] = Restrict;
//
Field[2].EdgesList = {1, 2};
//
Background Field = 2;
//
Mesh.Algorithm = 1;
//
Mesh.RecombinationAlgorithm = 2;
//
Mesh.RecombineAll = 1;
//
Mesh 2;
//
Mesh.SaveAll = 1;
//
Mesh.Binary = 0;
Mesh.MshFileVersion = 4.1;
//
Save "70degCone_2D.msh";
