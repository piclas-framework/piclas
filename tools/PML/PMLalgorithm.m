clear ; clc; close all;
cd(fileparts(which(mfilename))); %avoid problems if splashingalgorithm is already an added path; path directories may fail
matpath = pwd;
home = [pwd,'/..'];
addpath(matpath);
addpath([matpath filesep 'export_fig']);
disp('==============================================================');
disp('                    Running PML algorithm                     ');
disp('==============================================================');
[~ , machine ] = unix('whoami');
machine=strtrim(machine);
[~, nodes]=unix('cat machines | wc -l');
nodes=strtrim(nodes);
if length(nodes) > 30
  nodes='2';
end
disp(['nodes             = ' nodes]);
%% Load input file for preproc/mesh/dipole/parametervisu 
input=dipoleInput;

%% Setup
% Misc
debug                                                       = 'true';    % true = run flexi in single mode without MPI
create_plots                                                = 'false';    % Visu3D + Visit Plots + Video
input.ProjectName                                           = 'Dipole';   % Dipole or SingleParticle
copy                                                        = 1;          % Copy exe files (1) or not (0)
                                                                           % 1 = local storage (here)
                                                                           % 0 = remote stograge (non Linux partition)
flexi_source                                                = 'Flexi/Boltzplatz4'; % Source Folder
% Flexi 
input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2} = 12;        % Number of Cells in the physical domain +2*PML      particle: 20
input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2} = 6;         % Physical Domain Size   and Corner_new = Corner_old + 2*Corner_old/nElems*PML_Layer   particle:1E-3
input.Dipole{strcmp(input.Dipole,'N'),2}                    = 3;         % DG polynomial order
input.Dipole{strcmp(input.Dipole,'tend'),2}                 = 100E-9;% End Time Simulation in [ns]  particle: 200E-12
input.Dipole{strcmp(input.Dipole,'Analyze_dt'),2}           = 10E-9; % output time in [ns]         particle: 20E-12

input.frequency                                             = 100E6;     % frequency:[ default = 100MHz
input.Dipole{strcmp(input.Dipole,'CFLscale'),2}             = 0.75;      % CFL Number
input.Dipole{strcmp(input.Dipole,'AlphaShape'),2}           = 2;         % shape function type
input.Dipole{strcmp(input.Dipole,'r_cutoff'),2}             = 0.5;       % shape function cut off radius
input.Dipole{strcmp(input.Dipole,'GeometricNGeo'),2}        = 1;         % Degree of mesh representation
input.Dipole{strcmp(input.Dipole,'NAnalyze'),2}             = 10;        % Number of analyze points
input.Dipole{strcmp(input.Dipole,'TIMEDISCMETHOD'),2}       = 2;         % Time Discretisation: 1=RK3, 2=RK4 DON'T CHANGE !!! (PML in RK3 not implemented)
input.Dipole{strcmp(input.Dipole,'IniExactFunc'),2}         = 41;         % Exact Funktion (4  = contiunous dipole, 41 = pulsed dipole (t0=30ns/5 und w=30ns/20))
input.Dipole{strcmp(input.Dipole,'zetaShape'),2}            = 2;         % Damping Coefficient distribution in PML region (0 = constant, 1 = linear, 2 = sinusoidal)
input.Dipole{strcmp(input.Dipole,'PMLspread'),2}            = 0;         % copy zeta_i value to all zetas x,y,z (0 = don't spread, 1 = spread)
input.Dipole{strcmp(input.Dipole,'c_corr'),2}               = 1;         % chi: Divergence Cleaning Velocity
% Visit
input.plotType                                              = 1; %Visit Plot Type (1 = ThreeSlice, 2 = Isosurface, 3 = Contour)
% PML
input.PMLwriteZeta                                          = 1; % write zeta field=1,  dont write zeta field=0 (for debugging)
input.PML_Layer                                             = 2; % PML Layers
input.zeta=[0, 4E8];
%% Initialize
input=calcQuantities(input);

%%   
if input.frequency/1E6 ~= 100
  input.namefix=['_' num2str(input.frequency/1E6) 'MHz'];
else
  input.namefix =[''];
end
if input.Dipole{strcmp(input.Dipole,'zetaShape'),2} ~= 2
  switch input.Dipole{strcmp(input.Dipole,'zetaShape'),2}
    case 0
      input.namefix=[input.namefix '_constant'];
    case 1
      input.namefix=[input.namefix '_linear'];
  end
end

if input.Dipole{strcmp(input.Dipole,'c_corr'),2} ~= 1
  input.namefix=[input.namefix '_chi' num2str(input.Dipole{strcmp(input.Dipole,'c_corr'),2})];
end

disp(['RP01: ' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point1'),2})])
disp(['RP02: ' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point2'),2})])
disp(['RP03: ' num2str(input.RecordPoints{strcmp(input.RecordPoints,'Point3'),2})])
%
order=linspace(input.Dipole{strcmp(input.Dipole,'N'),2},input.Dipole{strcmp(input.Dipole,'N'),2},1);
%elements=linspace(6,18,7);
%% Create case folder 
for K=1:length(input.zeta)
  input.Dipole{strcmp(input.Dipole,'zeta0'),2} = input.zeta(K);
  disp(['zeta(' num2str(K) ')         = ' num2str(input.zeta(K))]);
  %input.Dipole{strcmp(input.Dipole,'N'),2}      = order(K);
  %input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}=elements(K);
  if input.Dipole{strcmp(input.Dipole,'zeta0'),2} == 0
    input.CellRatio=0;
  else
    input.CellRatio=(input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}^3-input.CellRatioCore^3)/input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}^3;
  end
  disp(['PML Cell %        = ' num2str(input.CellRatio*100)]);
  disp(['DOF in PML        = ' num2str(input.CellRatio*input.DOF)]);
  %input.DOF=(input.Dipole{strcmp(input.Dipole,'N'),2}+1)^3*input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}^3;
  if copy == 1
    casefolder=home;
  else
    casefolder='/media/Stephen 1TB_';
  end
  
  if input.Dipole{strcmp(input.Dipole,'zeta0'),2} == 0
    casefolder=[casefolder '/PMLalgorithm_cases' filesep 'domain' sprintf('%03.0f',input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) '_nElems' sprintf('%03.0f',input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}) '_order'  sprintf('%02.0f',input.Dipole{strcmp(input.Dipole,'N'),2}) '_PML'  sprintf('%02.0f',input.PML_Layer) '_zeta0E+00' input.namefix];
  else
    casefolder=[casefolder '/PMLalgorithm_cases' filesep 'domain' sprintf('%03.0f',input.PreprocDipole{strcmp(input.PreprocDipole,'Corner'),2}) '_nElems' sprintf('%03.0f',input.PreprocDipole{strcmp(input.PreprocDipole,'nElems'),2}) '_order'  sprintf('%02.0f',input.Dipole{strcmp(input.Dipole,'N'),2}) '_PML'  sprintf('%02.0f',input.PML_Layer) '_zeta' num2str(input.Dipole{strcmp(input.Dipole,'zeta0'),2}) input.namefix];
  end
  
  if  exist(casefolder,'dir')==7 % Wenn Ordner existiert
    rmdir(casefolder,'s') %Entfernen mit Subdirectories
  end
  mkdir(casefolder)
  
  
  %% write file for 1. preproc (mesh) 2. flexi (dipole) 3. visit (python script) & copy executables flexi/preproctool/visu3D
  writeIniFiles(input,casefolder);
  
  if copy ==1
    if strcmp(machine,'stephen')
      copyfile([pwd filesep '../' flexi_source '/bin/flexi'],[casefolder filesep '.']);
      %copyfile([pwd filesep 'flexi'],[casefolder filesep '.']);
    else
      copyfile([pwd filesep 'flexi'],[casefolder filesep '.']);
    end
    copyfile([pwd filesep 'preproctool'],[casefolder filesep '.']);
    copyfile([pwd filesep 'recordpoints'],[casefolder filesep '.']);
    copyfile([pwd filesep 'visu3D'],[casefolder filesep '.']);
    copyfile([pwd filesep 'postrec'],[casefolder filesep '.']);
    sourcepath='';
  else
    copyfile([pwd filesep '../' flexi_source '/bin/flexi'],[matpath filesep '.']);
    sourcepath=[matpath filesep];
  end  
  cd(casefolder);
  
  %% preproctool
  input.time_preproctool=run('preproctool',sourcepath,'./preproctool parameter_preproctool.ini');
  
  %% recordpoints
  input.time_recordpoints=run('recordpoints',sourcepath,['./recordpoints parameter_recordpoints.ini ' input.ProjectName '_mesh.h5']);
  
  %% flexi
  input.time_start=clock;
  if strcmp(machine,'stephen')
    if strcmp(debug,'true')
      input.time_flexi2=run('flexi',[matpath filesep],'./flexi_single parameter_flexi.ini 1>std.out 2>err.out');
    else
      input.time_flexi2=run('flexi',['mpirun -np 2 ' sourcepath],'./flexi parameter_flexi.ini 1>std.out 2>err.out');
    end
  else
    input.time_flexi2=run('flexi',sourcepath,['mpirun -np ' nodes ' ./flexi parameter_flexi.ini 1>std.out 2>err.out']);
  end
  input.time_end=clock;
  input.time_flexi=dot((input.time_end-input.time_start)',[0,0,86400,3600,60,1]'); %[year month day hour minute seconds]
  unix(['touch ' num2str(input.time_flexi)]);
  [~, ~] = unix(['mkdir ' input.ProjectName '_State && mv ' input.ProjectName '_State_* ' input.ProjectName '_State']); % move state files
  [~, ~] = unix(['mkdir ' input.ProjectName '_RP && mv ' input.ProjectName '_RP_* ' input.ProjectName '_RP']);          % move record points
  
  %% visu3D
  if strcmp(create_plots,'true') % Mannheim Cluster Visit File Creation does not work properly?!
    input.time_visu3D=run('visu3D',sourcepath,['./visu3D parameter_visu3D.ini ' input.ProjectName  '_State/' input.ProjectName '_State_00.0000*']);
    [~, ~] = unix('mkdir visit_vtk_files && mv *vtk visit_vtk_files');
    
    %% Visit
    if strcmp(machine,'stephen')
      input.time_visit=run('Visit','','/home/stephen/visit2_6_0.linux-x86_64/bin/./visit -nowin -cli -s parameter_visit.py');
    else
      input.time_visit=run('Visit','','visit -nowin -cli -s parameter_visit.py');
    end
    disp('-------------------------------------------------------------------');
    disp(['Running avconv']);
    disp(['Running Command: avconv -r 10 -i mybmpfile%4d.bmp -b:v 1000k ' input.ProjectName '.mp4'])
    tic
    [~, ~] = unix(['avconv -r 10 -i mybmpfile%4d.bmp -b:v 1000k ' input.ProjectName '.mp4']);
    input.time_avconv=toc;
    disp(['End avconv: ' num2str(input.time_avconv)]);
    disp('-------------------------------------------------------------------');
    unix('mkdir visit_pictures && mv mybmpfile* visit_pictures');
    
  end
  %% postrec
  disp('-------------------------------------------------------------------');
  disp(['Running postrec']);
  disp(['Running Command: ' sourcepath './postrec parameter_postrec.ini ' input.ProjectName '_RP/' input.ProjectName '_RP_0*']);
  tic
  [status, ~] = unix([sourcepath './postrec parameter_postrec.ini ' input.ProjectName '_RP/' input.ProjectName '_RP_0*']);
  input.time_postrec=toc;
  disp(['End postrec: ' num2str(input.time_postrec)]);
  disp('-------------------------------------------------------------------');
  if status ~=0
    fileList = getAllFiles([pwd filesep input.ProjectName '_RP/']);    % alle Dateien im Ordner in eine Liste setzen
    if ~isempty(fileList)
      A=cell2mat(fileList(end)); lastname=A(max(strfind(A,'/'))+1:end);
      [~, ~] = unix(['mv ' input.ProjectName '_RP/' lastname  ' moved_' lastname]);
      input.time_postrec=run('postrec',sourcepath,['./postrec parameter_postrec.ini ' input.ProjectName '_RP/' input.ProjectName '_RP_*']);
    end
  end
  unix('mkdir Probes &&  mv *RP0* Probes');
  
  %--------------------------------------------------------------------------
  % Record points -
  %--------------------------------------------------------------------------
  record_points(input,casefolder,matpath)
  
  %--------------------------------------------------------------------------
  % Potential Energy
  %--------------------------------------------------------------------------
  potential_energy(input,casefolder)
  
  %--------------------------------------------------------------------------
  % Clean up
  %--------------------------------------------------------------------------
  Clean_up(matpath,copy,input)
  
  %
  save([casefolder '/CaseInput.mat'],'input');
  disp('==============================================================');
  disp(['                        Done ... ' num2str(K) ' of ' num2str(length(input.zeta))]);
  disp('==============================================================');
end
disp('==============================================================');
disp('             PML algorithm completed successfully             ');
disp('==============================================================');

