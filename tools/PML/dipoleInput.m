function [input]=default
input.PreprocDipole(:,1)={'Mode';'nZones';'Corner';'nElems';'BCIndex';'elemtype';'useCurveds';'SpaceQuandt';'ConformConnect';'nUserDefinedBoundaries';'BoundaryName';'BoundaryType';'NVisu'};
input.PreprocDipole(:,2)={1;1;5;5;[1 1 1 1 1];108;'F';1;'T';1;'BC_outflow';[3 0 0 0];7};
input.Dipole(:,1)={'TIMEDISCMETHOD';'IniExactFunc';'N';'GeometricNGeo';'NAnalyze';'tend';'Analyze_dt';'CFLscale';'AlphaShape';'r_cutoff';'c_corr';'zeta0';'xyzPhysicalMinMax';'zetaShape';'PMLspread'};
input.Dipole(:,2)={2;4;4;4;10;100E-9;1E-9;1;2;0.5;1.0;0.e8;[-1.0 1.0 -1.0 1.0 -1.0 1.0];0;0};
input.RecordPoints(:,1)={'Point1';'Point2';'Point3'};
input.RecordPoints(:,2)={[5.9 0.0 0.0];[5.9 5.9 0.0];[5.9 5.9 5.9]};