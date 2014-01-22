function input=calcQuantities(input)

input.omega                                                 = input.frequency*2*acos(-1); % dipole angular frequency = 628 Mhz: w=2*pi*f
input.PreprocDipole{strcmp(input.PreprocDipole,'NVisu'),2}  = input.Dipole{strcmp(input.Dipole,'N'),2}+1;         % Number of visualization points

input.CellRatioCore=input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2};
input.CellLength = 2*input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}/input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2};
input.Dipole{strcmp(input.Dipole,'xyzPhysicalMinMax'),2}    = input.Dipole{strcmp(input.Dipole,'xyzPhysicalMinMax'),2}*input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2};
switch input.ProjectName
  case 'SingleParticle'                                                                          
    input.RecordPoints{strcmp(input.RecordPoints,'Point1'),2}=[0.9E-3 0 0];
    input.RecordPoints{strcmp(input.RecordPoints,'Point2'),2}=[0.9E-3 0.9E-3 0];
    input.RecordPoints{strcmp(input.RecordPoints,'Point3'),2}=[0.9E-3 0.9E-3 0.9E-3];
end
input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2} = input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}+2*input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}/input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}*input.PML_Layer;
input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2} = input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}+2*input.PML_Layer;
input.DOF         = input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}^3*(input.Dipole{strcmp(input.Dipole,'N'),2}+1)^3;
input.delta = input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}-input.Dipole{strcmp(input.Dipole,'xyzPhysicalMinMax'),2}(2);

disp(['xyzPhysicalMinMax = ' num2str(input.Dipole{strcmp(input.Dipole,'xyzPhysicalMinMax'),2})]);
disp(['Corner            = ' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2})]);
disp(['nElems            = ' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2})]);
disp(['polynomial degree = ' num2str(input.Dipole{strcmp(input.Dipole,'N'),2})]);
disp(['DOF               = ' num2str(input.DOF)]);
disp(['delta             = ' num2str(input.delta) ' (PML thickness)']);
disp(['zeta              = ' num2str(input.zeta)]);
disp(['omega             = ' num2str(input.frequency/1E6) ' MHz']);
if input.frequency ~0
  input.lambda = 299792458/input.frequency;
  input.PPW    = (input.Dipole{strcmp(input.Dipole,'N'),2}+1)*input.lambda/input.CellLength;
else
  input.lambda = 0;
  input.PPW    = 0;
end
disp(['lambda            = ' num2str(input.lambda) ]);
disp(['PPW               = ' num2str(input.PPW) ' (DOF pro Wellenlänge)']);
end