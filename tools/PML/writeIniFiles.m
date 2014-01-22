function writeIniFiles(input,casefolder)

names=fieldnames(input);

%% write parameter_preproc.ini
fid = fopen([casefolder filesep 'parameter_preproctool.ini'],'w');

fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! MAKEFILE PARAMETER (put a "#" in front, NO blanks!)']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! This is only a dummy parameter needed for the regression check']);
fprintf(fid,'%s\n\n',['#MPI=']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! OUTPUT']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['  ProjectName   = ' input.ProjectName '                      ! name of the project (used for filenames)']);
fprintf(fid,'%s\n',['  Debugvisu     =F                           ! Write debug mesh to tecplot file']);
fprintf(fid,'%s\n\n',['  Logging       =F                           ! Write log files']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! MESH']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['  Mode          =' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Mode'),2}) '                           ! 1 Cartesian 2 gambit file 3 CGNS ']);
fprintf(fid,'%s\n',['  nZones        =' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'nZones'),2}) '                           ! number of zones']);
fprintf(fid,'%s\n',['  Corner        =(/-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ' ,,-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',,-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) '/) ! [-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ']x[-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ']x[-' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) ']']);
fprintf(fid,'%s\n',['  nElems        =(/' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}) '/)                ! Anzahl der Elemente in jede Richtung']);
fprintf(fid,'%s\n',['  BCIndex       =(/' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(1)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(2)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(3)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(4)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(4)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BCIndex'),2}(5)) '/)             ! Indices of UserDefinedBoundaries']);
fprintf(fid,'%s\n',['  elemtype      =' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'elemtype'),2}) '                         ! Elementform (108: Hexaeder)']);
fprintf(fid,'%s\n',['  useCurveds    =' input.PreprocDipole{strcmp(input.PreprocDipole,'useCurveds'),2} '                           ! T if curved boundaries defined']);
fprintf(fid,'%s\n',['  SpaceQuandt   =' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'SpaceQuandt'),2}) '                           ! characteristic length of the mesh']);
fprintf(fid,'%s\n\n',['  ConformConnect=' input.PreprocDipole{strcmp(input.PreprocDipole,'ConformConnect'),2} ]);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! BOUNDARY CONDITIONS']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['  nUserDefinedBoundaries=1']);
fprintf(fid,'%s\n',['    BoundaryName=BC_outflow                  ! Outflow: open (absorbing)   [for MAXWELL]']);
fprintf(fid,'%s\n\n',['    BoundaryType=(/' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BoundaryType'),2}(1)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BoundaryType'),2}(2)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BoundaryType'),2}(3)) ',' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'BoundaryType'),2}(4)) '/)                 ! Outflow: open (absorbing)   [for MAXWELL]']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! BASIS']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n\n',['  NVisu         = ' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'NVisu'),2}) ]);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['! SEARCH']);
fprintf(fid,'%s\n',['!=============================================================================== !']);
fprintf(fid,'%s\n',['!  nElemsNodeSearch=50']);
fprintf(fid,'%s\n',['!  RefineSideSearch=50']);

fclose(fid);






%% write parameter_preproc.ini
fid = fopen([casefolder filesep 'parameter_recordpoints.ini'],'w');


fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! OUTPUT ']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['ProjectName   = ' input.ProjectName ' ']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! INTERPOLATION']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['GeometricNGeo = 4  ! Degree of mesh representation']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! MESH']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['MeshFile      = ' input.ProjectName '_mesh.h5']);
fprintf(fid,'%s\n',['useCurveds    =']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! RECORDPOINTS DEFINITION']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['NSuper=15     ! resolution of equidistan. search grid in element']);
fprintf(fid,'%s\n',['! last entry is group number']);
fprintf(fid,'%s\n',['GroupName =Gruppe1  ! Group ID 1']);
fprintf(fid,'%s\n\n',['!GroupName =Line1']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! Order of RPs is changed in the RPSet.h5']);
fprintf(fid,'%s\n\n',['! =============================================================================== !']);

fprintf(fid,'%s\n',['Point_GroupID = 1']);
fprintf(fid,'%s\n\n',['Point_x =(/' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point3'),2}(1)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point3'),2}(2)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point3'),2}(3)) '/)']);

fprintf(fid,'%s\n',['Point_GroupID = 1']);
fprintf(fid,'%s\n\n',['Point_x=(/' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point2'),2}(1)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point2'),2}(2)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point2'),2}(3)) '/)']);

fprintf(fid,'%s\n',['Point_GroupID = 1']);
fprintf(fid,'%s\n\n',['Point_x=(/' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point1'),2}(1)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point1'),2}(2)) ',' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point1'),2}(3)) '/)']);

fclose(fid);






%% write parameter_flexi.ini
fid = fopen([casefolder filesep 'parameter_flexi.ini'],'w');

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! MAKEFILE PARAMETER (put a "#" in front, NO blanks!)']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! eqn: navierstokes, linearscalaradvection, maxwell']);
fprintf(fid,'%s\n',['#EQNSYS=maxwell']);
fprintf(fid,'%s\n',['! explicit time discretization : 1...RK3, 2...RK4, 3...Taylor']);
fprintf(fid,'%s\n',['#TIMEDISCMETHOD=' num2str(input.Dipole{strcmp(input.Dipole,'TIMEDISCMETHOD'),2}) '']);
fprintf(fid,'%s\n',['! node type: 1...Gauss, 2...Gauss-Lobatto']);
fprintf(fid,'%s\n',['#NODETYPE=1']);
fprintf(fid,'%s\n',['! Riemann solver: 1...LF, 2...HLLC, 3...Roe']);
fprintf(fid,'%s\n',['#RIEMANN=1']);
fprintf(fid,'%s\n',['! Parallel execution: EMPTY...Off, T...On (MPI)']);
fprintf(fid,'%s\n',['#MPI=']);
fprintf(fid,'%s\n',['! optional: fixed number of elements']);
fprintf(fid,'%s\n',['#NELEMS=']);
fprintf(fid,'%s\n',['! optional: fixed polynomial degree']);
fprintf(fid,'%s\n',['#N=']);
fprintf(fid,'%s\n',['! optimizations ignoring inner array bounds (EMPTY...Off, T...On)']);
fprintf(fid,'%s\n',['! (cause errors when using array bound checks, always switched of in debug mode)']);
fprintf(fid,'%s\n\n',['#OPTIMIZED=T']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! EQUATION (linearscalaradvection)']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',  ['IniExactFunc  = ' num2str(input.Dipole{strcmp(input.Dipole,'IniExactFunc'),2}) '       ! 4=dipole, 41=pulsed dipole']);
fprintf(fid,'%s\n\n',['omega         = ' sprintf('%1.2E',input.omega) ' ! dipole angular frequency (f = ' num2str(input.frequency/1E6) ' MHz)']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! DISCRETIZATION']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['N             = ' num2str(input.Dipole{strcmp(input.Dipole,'N'),2}) '  ! Polynomial degree']);
fprintf(fid,'%s\n',['GeometricNGeo = ' num2str(input.Dipole{strcmp(input.Dipole,'GeometricNGeo'),2}) '  ! Degree of mesh representation']);
fprintf(fid,'%s\n\n',['NAnalyze      = ' num2str(input.Dipole{strcmp(input.Dipole,'NAnalyze'),2}) ' ! Number of analyze points']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! MESH']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['MeshFile      = ' input.ProjectName '_mesh.h5 ']);
fprintf(fid,'%s\n\n',['useCurveds    = F']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! OUTPUT / VISUALIZATION']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['ProjectName     = ' input.ProjectName]);
fprintf(fid,'%s\n',['OutputFormat    = 1    ! 0...Tecplot (only PostProcTool)']);
fprintf(fid,'%s\n',['ContinuousVisu  = 0    ! 0 - False | 1 - True | 2 - Both']);
fprintf(fid,'%s\n',['NVisu           = ' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'NVisu'),2}) '    ! Number of visualization points']);
fprintf(fid,'%s\n',['WriteErrorFiles = F']);
fprintf(fid,'%s\n\n',['Logging         = F']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! CALCULATION']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['tend       = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'tend'),2}) ' ! End time']);
fprintf(fid,'%s\n',['Analyze_dt = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'Analyze_dt'),2}) ' ! Timestep of analyze outputs ! encounter BC at 6.67e-9']);
fprintf(fid,'%s\n',['CFLscale   = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'CFLscale'),2}) '  ! Scaling of theoretical CFL number']);
fprintf(fid,'%s\n',['AlphaShape = ' num2str(input.Dipole{strcmp(input.Dipole,'AlphaShape'),2}) ]);
fprintf(fid,'%s\n',['r_cutoff   = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'r_cutoff'),2}) ]);
fprintf(fid,'%s\n',['c_corr     = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'c_corr'),2}) ]);
fprintf(fid,'%s\n',['c0         = 299792458']);
fprintf(fid,'%s\n',['eps        = 8.85418782e-12']);
fprintf(fid,'%s\n\n',['mu         = 12.566370614e-7']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! PML (perfectly matched layer)']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['zeta0             = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'zeta0'),2}) '! damping constant in PML region, physical regions zeta=0']);
fprintf(fid,'%s\n',['zetaShape         = ' num2str(input.Dipole{strcmp(input.Dipole,'zetaShape'),2}) '        ! shape function for damping constant (0=const, 1=linear, 2=sinus)']);
fprintf(fid,'%s\n',['PMLspread         = ' num2str(input.Dipole{strcmp(input.Dipole,'PMLspread'),2}) '        ! spread=1 dont spread=0']);
fprintf(fid,'%s\n',['PMLwriteZeta      = ' num2str(input.PMLwriteZeta) '        ! write zeta field=1,  dont write zeta field=0']);
fprintf(fid,'%s\n\n',['xyzPhysicalMinMax = (/' sprintf('%1.2E,',input.Dipole{strcmp(input.Dipole,'xyzPhysicalMinMax'),2}) '/) ! lower/upper boarder ']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! PARTICLES']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['Part-nBounds              = 1']);
fprintf(fid,'%s\n',['Part-Boundary1-SourceName = BC_outflow']);
fprintf(fid,'%s\n\n',['Part-Boundary1-Condition  = open']);

if strcmp(input.ProjectName,'SingleParticle')
    fprintf(fid,'%s\n',['Part-maxParticleNumber=250']);
    fprintf(fid,'%s\n',['Part-nSpecies=1']);
    fprintf(fid,'%s\n',['!PIC-DoInterpolation=TRUE !FALSE']);
    fprintf(fid,'%s\n',['PIC-externalField=(/0.,0.,0.,0.,0.,0.,0.,0./)']);
    fprintf(fid,'%s\n',['!Part-Species1-SpaceIC=line_with_equidistant_distribution']);
    fprintf(fid,'%s\n',['!Part-Species1-BasePointIC=(/0.0,-9.62962963E-4,0.0/)']);
    fprintf(fid,'%s\n',['!Part-Species1-BaseVector1IC=(/0.0,9.62962963E-4,0.0/)']);
    fprintf(fid,'%s\n',['Part-Species1-SpaceIC=point']);
    fprintf(fid,'%s\n',['Part-Species1-BasePointIC=(/0.0,0.0,0.0/)']);
    fprintf(fid,'%s\n',['Part-Species1-velocityDistribution=constant']);
    fprintf(fid,'%s\n',['Part-Species1-initialParticleNumber=1']);
    fprintf(fid,'%s\n',['Part-Species1-ParticleEmissionType=1']);
    fprintf(fid,'%s\n',['Part-Species1-ParticleEmission=0']);
    fprintf(fid,'%s\n',['Part-Species1-VeloIC=0 ! geändert von 1E7']);
    fprintf(fid,'%s\n',['Part-Species1-VeloVecIC=(/0.,1.,0./)']);
    fprintf(fid,'%s\n',['Part-Species1-ChargeIC=1.60217653E-19']);
    fprintf(fid,'%s\n',['Part-Species1-MassIC=1']);
    fprintf(fid,'%s\n\n',['Part-Species1-MacroParticleFactor=1 ! geändert von 1E9']);
    
    fprintf(fid,'%s\n',['Part-Species2-SpaceIC=line_with_equidistant_distribution']);
    fprintf(fid,'%s\n',['Part-Species2-BasePointIC=(/-0.001,0.0,0.0/)']);
    fprintf(fid,'%s\n',['Part-Species2-BaseVector1IC=(/0.002,0.0,0.0/)']);
    fprintf(fid,'%s\n',['Part-Species2-velocityDistribution=constant']);
    fprintf(fid,'%s\n',['Part-Species2-initialParticleNumber=201']);
    fprintf(fid,'%s\n',['Part-Species2-ParticleEmissionType=1']);
    fprintf(fid,'%s\n',['Part-Species2-ParticleEmission=0']);
    fprintf(fid,'%s\n',['Part-Species2-VeloIC=0']);
    fprintf(fid,'%s\n',['Part-Species2-VeloVecIC=(/1.,0.,0./)']);
    fprintf(fid,'%s\n',['Part-Species2-ChargeIC=0']);
    fprintf(fid,'%s\n',['Part-Species2-MassIC=1']);
    fprintf(fid,'%s\n\n',['Part-Species2-MacroParticleFactor=1']);
    %                                         Part-nBounds=1
    %                                         Part-Boundary1-SourceName=BC_open
    %                                         Part-Boundary1-Condition=open
    
    fprintf(fid,'%s\n',['PIC-Interpolation-Type=particle_position']);
    fprintf(fid,'%s\n',['PIC-Deposition-Type=nearest_blurycenter']);
    fprintf(fid,'%s\n',['PIC-shapefunction-radius=0.02']);
    fprintf(fid,'%s\n',['PICshapefunction-alpha=2']);
    fprintf(fid,'%s\n',['Part-FIBGMdeltas=(/0.02,0.02,0.02/)']);
    fprintf(fid,'%s\n',['PIC-Depo-Periodic=FALSE']);
    fprintf(fid,'%s\n',['CalcCharge = FALSE']);
    fprintf(fid,'%s\n',['CalcPotentialEnergy=TRUE']);
    fprintf(fid,'%s\n',['CalcKineticEnergy=TRUE']);
    
    fprintf(fid,'%s\n',['Part-NumberOfRandomSeeds =2']);
    fprintf(fid,'%s\n',['Particles-RandomSeed1= 1180520427']);
    fprintf(fid,'%s\n',['Particles-RandomSeed2= 1708457652']);
end

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! Analysis']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['CalcPotentialEnergy      = TRUE']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! Recordpoints']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['RP_OutputInterval        = 100']);
fprintf(fid,'%s\n',['RP_SamplingOffset        = 1 ! only each x time step']);
fprintf(fid,'%s\n',['RP_DefFile               = ' input.ProjectName '_RPSet.h5']);


fclose(fid);


 










%% write parameter_visu3D.ini
fid = fopen([casefolder filesep 'parameter_visu3D.ini'],'w');

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! MAKEFILE PARAMETER (put a "#" in front, NO blanks!)']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! eqn system: navierstokes, linearscalaradvection, maxwell']);
fprintf(fid,'%s\n',['#EQNSYS=maxwell']);
fprintf(fid,'%s\n',['! EOS: equation of state']);
fprintf(fid,'%s\n',['#EOS=idealgas']);
fprintf(fid,'%s\n',['! Parallel execution: EMPTY...Off, T...On (MPI)']);
fprintf(fid,'%s\n',['#MPI=T']);
fprintf(fid,'%s\n\n',['!R=4']);

fprintf(fid,'%s\n',['NVisu         = ' num2str(input.PreprocDipole{strcmp(input.PreprocDipole,'NVisu'),2})]);
fprintf(fid,'%s\n',['NodeType_Visu = VISU']);
fprintf(fid,'%s\n',['useCurveds    = F']);
fprintf(fid,'%s\n',['!VarName       = ElectricFieldX']);
fprintf(fid,'%s\n',['!VarName       = Density']);
fprintf(fid,'%s\n',['!Varname       = MomentumX']);
fprintf(fid,'%s\n',['!Varname       = MomentumY']);
fprintf(fid,'%s\n',['!Varname       = MomentumZ']);
fprintf(fid,'%s\n',['!Varname       = EnergyStagnationDensity']);
fprintf(fid,'%s\n',['!Varname       = VelocityX']);
fprintf(fid,'%s\n',['!Varname       = VelocityY']);
fprintf(fid,'%s\n',['!Varname       = VelocityZ']);
fprintf(fid,'%s\n',['!Varname       = VelocityMagnitude']);
fprintf(fid,'%s\n',['!Varname       = Pressure']);
fprintf(fid,'%s\n',['!VarName       = VelocitySound']);
fprintf(fid,'%s\n',['!VarName       = Temperature']);
fprintf(fid,'%s\n',['!Varname       = EnergyStagnation']);
fprintf(fid,'%s\n',['!Varname       = EnthalpyStagnation']);
fprintf(fid,'%s\n',['!Varname       = Entropy']);
fprintf(fid,'%s\n',['!VarName       = Jacobian']);
fprintf(fid,'%s\n',['!VarName       = PotTemp']);
fprintf(fid,'%s\n',['!VarName       = VorticityX']);
fprintf(fid,'%s\n',['!VarName       = VorticityY']);
fprintf(fid,'%s\n',['!VarName       = VorticityZ']);
fprintf(fid,'%s\n',['!VarName       = VorticityMagnitude']);
fprintf(fid,'%s\n',['!VarName       = Helicity']);
fprintf(fid,'%s\n',['!VarName       = Lambda2']);
fprintf(fid,'%s\n',['!VarName       = Dilatation']);
fprintf(fid,'%s\n',['outputformat   = 4 ! 3: *.pvtu and *.vtu  ']);
fprintf(fid,'%s\n',['                   ! 4: *.vtk']);

fclose(fid);













%% write parameter_visit.py
fid = fopen([casefolder filesep 'parameter_visit.py'],'w');
switch input.ProjectName
  case 'Dipole'
    fprintf(fid,'%s\n',['print "=============================================================="']);
    fprintf(fid,'%s\n',['print "Starting Script .............................................."']);
    fprintf(fid,'%s\n',['SuppressMessages(3) # 1 suppresses all types of messages. ']);
    fprintf(fid,'%s\n',['                    # 2 suppresses Warnings and Messages but does NOT suppress Errors']);
    fprintf(fid,'%s\n',['                    # 3 suppresses Messages but does not suppress Warnings or Errors. ']);
    fprintf(fid,'%s\n',['                    # 4 does not suppress any messages. The default setting is 4.']);
    fprintf(fid,'%s\n',['s = SaveWindowAttributes()']);
    fprintf(fid,'%s\n',['s.format = s.BMP']);
    fprintf(fid,'%s\n',['s.fileName = "mybmpfile"']);
    fprintf(fid,'%s\n',['s.width, s.height = 1024,768']);
    fprintf(fid,'%s\n',['s.screenCapture = 0']);
    fprintf(fid,'%s\n\n',['SetSaveWindowAttributes(s)']);
    
    fprintf(fid,'%s\n',['# Fuer Isosurface']);
    fprintf(fid,'%s\n',['#v = GetView3D()']);
    fprintf(fid,'%s\n',['#v.viewNormal = (-0.670597, 0.352165, 0.652901)']);
    fprintf(fid,'%s\n',['#v.viewUp = (0.256093, 0.935924, -0.24179)']);
    fprintf(fid,'%s\n\n',['#SetView3D(v)']);
    
    fprintf(fid,'%s\n',['# Fuer ThreeSlice']);
    fprintf(fid,'%s\n',['v = GetView3D()']);
    fprintf(fid,'%s\n',['v.viewNormal = (0.589854, 0.692355, 0.415592)']);
    fprintf(fid,'%s\n',['v.viewUp = (-0.269786, -0.316121, 0.909551)']);
    fprintf(fid,'%s\n\n',['SetView3D(v)']);
    
    fprintf(fid,'%s\n',['pwd=os.getcwd() # get current directory path']);
    fprintf(fid,'%s\n',['import glob']);
    fprintf(fid,'%s\n',['List = glob.glob(pwd+"/visit_vtk_files/*.vtk")']);
    fprintf(fid,'%s\n',['List = sorted(List)']);
    
    fprintf(fid,'%s\n',['plotType = ' num2str(input.plotType) ' # 1 = ThreeSlice']);
    fprintf(fid,'%s\n',['             # 2 = Isosurface']);
    fprintf(fid,'%s\n\n',['             # 3 = Contour']);
    
    fprintf(fid,'%s\n',['for file in List:']);
    fprintf(fid,'%s\n',['    print file']);
    fprintf(fid,'%s\n',['    OpenDatabase(file)']);
    fprintf(fid,'%s\n',['    DefineScalarExpression("E_r", "ElectricFieldX*cos(cylindrical_theta(mesh)) + ElectricFieldY*sin(cylindrical_theta(mesh))")']);
    fprintf(fid,'%s\n\n',['    DefineScalarExpression("E_abs", "sqrt(ElectricFieldX^2+ElectricFieldY^2+ElectricFieldZ^2)")']);
    
    fprintf(fid,'%s\n',['    #ThreeSlice']);
    fprintf(fid,'%s\n',['    if plotType == 1:']);
    fprintf(fid,'%s\n',['      AddPlot("Pseudocolor","E_abs")']);
    fprintf(fid,'%s\n',['      p = PseudocolorAttributes()']);
    fprintf(fid,'%s\n',['      p.colorTableName = "hot_and_cold" # Fuer ThreeSlice']);
    fprintf(fid,'%s\n',['      #p.colorTableName = "rainbow"']);
    fprintf(fid,'%s\n',['      #p.opacity = 0.5']);
    fprintf(fid,'%s\n',['      p.min, p.minFlag = -1, 1']);
    fprintf(fid,'%s\n',['      p.max, p.maxFlag = 1, 1']);
    fprintf(fid,'%s\n',['      SetPlotOptions(p)']);
    fprintf(fid,'%s\n\n',['      AddOperator("ThreeSlice")']);
    
    fprintf(fid,'%s\n',['    #Isosurface']);
    fprintf(fid,'%s\n',['    elif plotType == 2:']);
    fprintf(fid,'%s\n',['      AddPlot("Pseudocolor","E_abs")']);
    fprintf(fid,'%s\n',['      p = PseudocolorAttributes()']);
    fprintf(fid,'%s\n',['      p.min, p.minFlag = -1, 1']);
    fprintf(fid,'%s\n',['      p.max, p.maxFlag = 1, 1']);
    fprintf(fid,'%s\n',['      p.opacity = 0.5']);
    fprintf(fid,'%s\n',['      SetPlotOptions(p)']);
    fprintf(fid,'%s\n',['      AddOperator("Isosurface")']);
    fprintf(fid,'%s\n',['      p2 = IsosurfaceAttributes()']);
    fprintf(fid,'%s\n',['      p2.min, p2.minFlag = -1, 1']);
    fprintf(fid,'%s\n',['      p2.max, p2.maxFlag = 1, 1']);
    fprintf(fid,'%s\n\n',['      SetOperatorOptions(p2)']);
    
    fprintf(fid,'%s\n',['    #Contour']);
    fprintf(fid,'%s\n',['    elif plotType == 3:']);
    fprintf(fid,'%s\n',['      AddPlot("Contour","E_abs")']);
    fprintf(fid,'%s\n',['      c = ContourAttributes()']);
    fprintf(fid,'%s\n',['      c.contourNLevels = 2']);
    fprintf(fid,'%s\n',['      c.min, c.minFlag = -1, 1']);
    fprintf(fid,'%s\n',['      c.max, c.maxFlag = 1, 1']);
    fprintf(fid,'%s\n',['      SetPlotOptions(c)']);
    fprintf(fid,'%s\n',['    else:']);
    fprintf(fid,'%s\n',['      print("Only single-digit numbers are allowed between 1 and 3");']);
    fprintf(fid,'%s\n\n',['      exit()']);
    
    fprintf(fid,'%s\n',['    #p2.contourMethod = p2.Value']);
    fprintf(fid,'%s\n',['    #p2.variable = "E_r"']);
    fprintf(fid,'%s\n\n',['    #p2.contourValue = 1']);

    fprintf(fid,'%s\n',['    anno=AnnotationAttributes()']);
    fprintf(fid,'%s\n',['    anno.axes3D.bboxFlag=0']);
    fprintf(fid,'%s\n',['    anno.axes3D.visible=0']);
    fprintf(fid,'%s\n',['    SetAnnotationAttributes(anno)']);
    fprintf(fid,'%s\n',['    DrawPlots()']);
    fprintf(fid,'%s\n',['    SaveWindow()']);
    fprintf(fid,'%s\n',['    DeleteAllPlots()']);
    fprintf(fid,'%s\n',['    CloseDatabase(file)']);
    fprintf(fid,'%s\n',['print "Done.........................................................."']);
    fprintf(fid,'%s\n',['print "=============================================================="']);
    fprintf(fid,'%s\n',['exit()']);
  case 'SingleParticle'
    fprintf(fid,'%s\n',['print "=============================================================="']);
    fprintf(fid,'%s\n',['print "Starting Script .............................................."']);
    fprintf(fid,'%s\n',['SuppressMessages(3) # 1 suppresses all types of messages. ']);
    fprintf(fid,'%s\n',['                    # 2 suppresses Warnings and Messages but does NOT suppress Errors']);
    fprintf(fid,'%s\n',['                    # 3 suppresses Messages but does not suppress Warnings or Errors. ']);
    fprintf(fid,'%s\n\n',['                    # 4 does not suppress any messages. The default setting is 4.']);
    
    fprintf(fid,'%s\n',['s = SaveWindowAttributes()']);
    fprintf(fid,'%s\n',['s.format = s.BMP']);
    fprintf(fid,'%s\n',['s.fileName = "mybmpfile"']);
    fprintf(fid,'%s\n',['s.width, s.height = 1024,768']);
    fprintf(fid,'%s\n',['s.screenCapture = 0']);
    fprintf(fid,'%s\n\n',['SetSaveWindowAttributes(s)']);
    
    fprintf(fid,'%s\n',['v = GetView3D()']);
    fprintf(fid,'%s\n',['v.viewNormal = (0.589854, 0.692355, 0.415592)']);
    fprintf(fid,'%s\n',['v.viewUp = (-0.269786, -0.316121, 0.909551)']);
    fprintf(fid,'%s\n\n',['SetView3D(v)']);
    
    fprintf(fid,'%s\n',['pwd=os.getcwd() # get current directory path']);
    fprintf(fid,'%s\n',['amplitude=0.1']);
    fprintf(fid,'%s\n',['print amplitude']);
    fprintf(fid,'%s\n',['import glob']);
    fprintf(fid,'%s\n',['List = glob.glob(pwd+"/visit_vtk_files/*.vtk")']);
    fprintf(fid,'%s\n',['List = sorted(List)']);
    fprintf(fid,'%s\n',['for file in List:']);
    fprintf(fid,'%s\n',['    print file']);
    fprintf(fid,'%s\n',['    OpenDatabase(file)']);
    fprintf(fid,'%s\n\n',['    DefineScalarExpression("E_abs", "sqrt(ElectricFieldX^2+ElectricFieldY^2+ElectricFieldZ^2)")']);
    
    fprintf(fid,'%s\n',['    AddPlot("Pseudocolor","E_abs")']);
    fprintf(fid,'%s\n',['    p = PseudocolorAttributes()']);
    fprintf(fid,'%s\n',['    #p.min, p.minFlag = -amplitude, 1']);
    fprintf(fid,'%s\n',['    p.max, p.maxFlag = amplitude*2, 1']);
    fprintf(fid,'%s\n',['    #p.opacity = 0.5']);
    fprintf(fid,'%s\n\n',['    SetPlotOptions(p)']);

    fprintf(fid,'%s\n',['    AddOperator("Isosurface")']);
    fprintf(fid,'%s\n',['    p2 = IsosurfaceAttributes()']);
    fprintf(fid,'%s\n',['    p2.contourNLevels=250']);
    fprintf(fid,'%s\n',['    p2.scaling=1']);
    fprintf(fid,'%s\n',['    #p2.min, p2.minFlag = -amplitude, 1']);
    fprintf(fid,'%s\n',['    p2.max, p2.maxFlag = amplitude*4, 1']);
    fprintf(fid,'%s\n\n',['    SetOperatorOptions(p2)']);

    fprintf(fid,'%s\n\n',['    AddOperator("ThreeSlice")']);

    fprintf(fid,'%s\n',['    anno=AnnotationAttributes()']);
    fprintf(fid,'%s\n',['    anno.axes3D.bboxFlag=0']);
    fprintf(fid,'%s\n',['    anno.axes3D.visible=0']);
    fprintf(fid,'%s\n',['    SetAnnotationAttributes(anno)']);
    fprintf(fid,'%s\n',['    DrawPlots()']);
    fprintf(fid,'%s\n',['    SaveWindow()']);
    fprintf(fid,'%s\n',['    DeleteAllPlots()']);
    fprintf(fid,'%s\n\n',['    CloseDatabase(file)']);
    
    fprintf(fid,'%s\n',['print "Done.........................................................."']);
    fprintf(fid,'%s\n',['print "=============================================================="']);
    fprintf(fid,'%s\n',['exit()']);
end
fclose(fid);





%% write parameter_postrec.ini
fid = fopen([casefolder filesep 'parameter_postrec.ini'],'w');

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! RECORDPOINTS POSTPROC']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n\n',['ProjectName   = ' input.ProjectName]);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! RP INFO          ']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['RP_DefFile=' input.ProjectName '_RPSet.h5']);
fprintf(fid,'%s\n\n',['GroupName=Gruppe1']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! TIME INTERVAL']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['!equiTimeSpacing  =F']);
fprintf(fid,'%s\n',['!CalcTimeAverage  =T']);
fprintf(fid,'%s\n',['!Line_LocalCoords =F']);
fprintf(fid,'%s\n',['!OutputTimeData    = T']);
fprintf(fid,'%s\n',['ZeroCrossing     =T']);
fprintf(fid,'%s\n',['OutputPureSignal =T ! further analysis']);
fprintf(fid,'%s\n',['!doFFT             =T']);
fprintf(fid,'%s\n',['FFTt0             = ' sprintf('%1.2E',input.Dipole{strcmp(input.Dipole,'Analyze_dt'),2}) '']);
fprintf(fid,'%s\n\n',['FFTtend           = ' sprintf('%1.2E',0.99*input.Dipole{strcmp(input.Dipole,'tend'),2}) '']);   % weniger als Tend, da es scheinbar nicht klappt, deshalb 99%

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! EQUATION     ']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['c_corr     = 1  ! divergence cleaning -> 1 badbad bad']);
fprintf(fid,'%s\n',['c0         = 299792458. ']);
fprintf(fid,'%s\n',['eps        = 8.8541878176E-12']);
fprintf(fid,'%s\n\n',['mu         = 12.566370614e-7 ']);

fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! OUTPUT VARIABLES']);
fprintf(fid,'%s\n',['! =============================================================================== !']);
fprintf(fid,'%s\n',['! conservative variables']);
fprintf(fid,'%s\n',['VarName=ElectricFieldX']);
fprintf(fid,'%s\n',['VarName=ElectricFieldY']);
fprintf(fid,'%s\n',['VarName=ElectricFieldZ']);
fprintf(fid,'%s\n',['VarName=MagneticFieldX']);
fprintf(fid,'%s\n',['VarName=MagneticFieldY']);
fprintf(fid,'%s\n',['VarName=MagneticFieldZ']);
fprintf(fid,'%s\n',['VarName=Phi']);
fprintf(fid,'%s\n',['VarName=Psi']);
fprintf(fid,'%s\n',['! derived variables ! only with ZeroCrossing']);
fprintf(fid,'%s\n',['VarName=ElectricFieldR']);
fprintf(fid,'%s\n',['VarName=ElectricFieldTheta']);
fprintf(fid,'%s\n',['VarName=MagneticFieldR']);
fprintf(fid,'%s\n',['VarName=MagneticFieldTheta']);
fprintf(fid,'%s\n',['VarName=MagElectric']);
fprintf(fid,'%s\n',['VarName=MagMagnetic']);

fclose(fid);
%field=getfield(input,names{1})
end


