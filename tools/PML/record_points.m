function record_points(input,casefolder,matpath)
order      = input.Dipole{strcmp(input.Dipole,'N'),2};
zeta       = input.Dipole{strcmp(input.Dipole,'zeta0'),2};
PML        = input.PML_Layer;
time       = input.time_flexi;
shape      = input.Dipole{strcmp(input.Dipole,'zetaShape'),2};
chi        = input.Dipole{strcmp(input.Dipole,'c_corr'),2};
PMLcells   = input.CellRatio;
PMLspread  = input.Dipole{strcmp(input.Dipole,'PMLspread'),2};
DOF        = input.DOF;
frequency  = input.frequency;
load([matpath '/Reference_Case.mat'])
switch input.Dipole{strcmp(input.Dipole,'IniExactFunc'),2}
  case 41 % pulsed dipole: in reference case data only order 2,3,4,5,6,7 are stored, therefore shift = -1 !!  p=2 -> RP(1), p=3 -> RP(2), p=4 -> RP(3) etc ...
    switch frequency/1E6
      case 100
        shift = -1; % kleinster polynomgrad 2
        orderShifted = order+shift;
      case 83
        switch order
          case 3
            orderShifted = case_start(3);
          case 7
            orderShifted = case_start(3)+3;
          otherwise
            orderShifted = case_end(3);
        end
      case 125
        switch order
          case 3
            orderShifted = case_start(4);
          case 7
            orderShifted = case_start(4)+3;
          otherwise
            orderShifted = case_end(4);
        end
      case 250
        switch order
          case 3
            orderShifted = case_start(5);
          case 7
            orderShifted = case_start(5)+3;
          otherwise
            orderShifted = case_end(5);
        end
      case 20
        switch order
          case 3
            orderShifted = case_start(6);
          case 7
            orderShifted = case_start(6)+3;
          otherwise
            orderShifted = case_end(6);
        end
    end
    disp(['orderShifted = ' num2str(orderShifted)])
    disp(['Pulsed dipole comparison ... ' num2str(frequency/1E6) ' MHz'])
  case 4 % continuous dipole needs different shift
    shift = -3; %eigentlich 4 % kleinster polynomgrad 5 für 5,6,7 was genau 3 ordner ergibt maximum = 10
    orderShifted = order+shift+case_start(2)-1;
    disp(['orderShifted = ' num2str(orderShifted)])
    disp(['Continuous dipole comparison ... ' num2str(frequency/1E6) ' MHz'])
  case 0
    switch order
      case 3
        orderShifted = case_start(7);
      case 4
        orderShifted = case_start(7)+1;
      case 6
        orderShifted = case_start(7)+2;
      otherwise % for order = 0 (ref case skript) or order = 7
        orderShifted = case_end(7);
    end
    disp(['orderShifted = ' num2str(orderShifted)])
    disp(['Single Particle comparison ... '])
end

%% ========================================================================
% E-Field absolute value
%==========================================================================
% Wenn Datei existiert: 'Time','ElectricFieldX','ElectricFieldY','ElectricFieldZ','MagneticFieldX', 'MagneticFieldY','MagneticFieldZ','Psi','Phi'
if  exist([casefolder '/Probes/RP01_ElectricMagneticField.dat'],'file')==2 
  %RP01 = textscan(fopen([casefolder '/Probes/RP01_ElectricMagneticField.dat'],'r'),'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
  %RP02 = textscan(fopen([casefolder '/Probes/RP02_ElectricMagneticField.dat'],'r'),'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
  %RP03 = textscan(fopen([casefolder '/Probes/RP03_ElectricMagneticField.dat'],'r'),'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
  %disp(['RecordPoint 1 = ' num2str(input.RecordPoints{1,2}) ]);
  %disp(['RecordPoint 2 = ' num2str(input.RecordPoints{2,2}) ]);
  %disp(['RecordPoint 3 = ' num2str(input.RecordPoints{3,2}) ]);
  for K=1:3
    fid=fopen([casefolder '/Probes/RP0' num2str(K) '_ElectricMagneticField.dat'],'r');
    a=textscan(fid,'%s %f %f %f',1,'HeaderLines',1)
    frewind(fid);
    for J = 1:3
      if abs(input.RecordPoints{J,2}-[a{1,2} a{1,3} a{1,4}])<=1E-6
        switch J
          case 1
            RP01 = textscan(fid,'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
          case 2
            RP02 = textscan(fid,'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
          case 3
            RP03 = textscan(fid,'%f %f %f %f %f %f %f %f %f', 'HeaderLines',3);
        end
      else
        
      end
      disp(['abs(input.RecordPoints{J,2}-[a{1,2} a{1,3} a{1,4}]) =' num2str(  abs(input.RecordPoints{J,2}-[a{1,2} a{1,3} a{1,4}])  )])
    end
    fclose(fid);
  end  
  
  t=RP01{1,1};
  % Plot E_abs
  f=figure; hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('E [V/m]');
  set(gca, 'ColorOrder', [0 0 1]);       RP01_E_abs=(RP01{1,2}.^2+RP01{1,3}.^2+RP01{1,4}.^2).^0.5;  E1_abs=plot(t*1e9,RP01_E_abs,'.-');
  set(gca, 'ColorOrder', [1 0 0]);       RP02_E_abs=(RP02{1,2}.^2+RP02{1,3}.^2+RP02{1,4}.^2).^0.5;  E2_abs=plot(t*1e9,RP02_E_abs,'.-');
  set(gca, 'ColorOrder', [0 0.4706 0]);  RP03_E_abs=(RP03{1,2}.^2+RP03{1,3}.^2+RP03{1,4}.^2).^0.5;  E3_abs=plot(t*1e9,RP03_E_abs,'.-');
  legend([E1_abs(1),E2_abs(1),E3_abs(1)],'Probe 1: E_{abs}','Probe 2: E_{abs}','Probe 3: E_{abs}','Location','EastOutside');
  export_fig([casefolder '/Probes/E_abs_RP.pdf']);
  close(gcf);
  %save([casefolder '/Probes/RP_order_'  sprintf('%02.0f',input.Dipole{find(strcmp(input.Dipole,'N')),2}) '.mat'],'t','RP01_E_abs','RP02_E_abs','RP03_E_abs');
  
  % Plot B_abs
  f=figure; hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('B [T]');
  set(gca, 'ColorOrder', [0 0 1]);       RP01_B_abs=(RP01{1,5}.^2+RP01{1,6}.^2+RP01{1,7}.^2).^0.5;  B1_abs=plot(t*1e9,RP01_B_abs,'.-');
  set(gca, 'ColorOrder', [1 0 0]);       RP02_B_abs=(RP02{1,5}.^2+RP02{1,6}.^2+RP02{1,7}.^2).^0.5;  B2_abs=plot(t*1e9,RP02_B_abs,'.-');
  set(gca, 'ColorOrder', [0 0.4706 0]);  RP03_B_abs=(RP03{1,5}.^2+RP03{1,6}.^2+RP03{1,7}.^2).^0.5;  B3_abs=plot(t*1e9,RP03_B_abs,'.-');
  legend([B1_abs(1),B2_abs(1),B3_abs(1)],'Probe 1: B_{abs}','Probe 2: B_{abs}','Probe 3: B_{abs}','Location','EastOutside');
  export_fig([casefolder '/Probes/B_abs_RP.pdf']);
  close(gcf);
  
  % Plot E_diff_percent
  RP01_E_abs_ref_temp=[RP01_ref{orderShifted}{1,1} (RP01_ref{orderShifted}{1,2}.^2+RP01_ref{orderShifted}{1,3}.^2+RP01_ref{orderShifted}{1,4}.^2).^0.5]; % RP01 Reference clean [t E_abs]
  RP02_E_abs_ref_temp=[RP02_ref{orderShifted}{1,1} (RP02_ref{orderShifted}{1,2}.^2+RP02_ref{orderShifted}{1,3}.^2+RP02_ref{orderShifted}{1,4}.^2).^0.5]; % RP02 Reference clean [t E_abs]
  RP03_E_abs_ref_temp=[RP03_ref{orderShifted}{1,1} (RP03_ref{orderShifted}{1,2}.^2+RP03_ref{orderShifted}{1,3}.^2+RP03_ref{orderShifted}{1,4}.^2).^0.5]; % RP03 Reference clean [t E_abs]
  RP01_E_abs_ref_interpol=zeros(length(t),1); % RP01 Reference interpolated [t E_abs]
  RP02_E_abs_ref_interpol=zeros(length(t),1); % RP01 Reference interpolated [t E_abs]
  RP03_E_abs_ref_interpol=zeros(length(t),1); % RP01 Reference interpolated [t E_abs]
  for J=1:length(RP01_E_abs)
    RP01_E_abs_ref_interpol(J)=interpolate(RP01_E_abs_ref_temp,t(J));
    RP02_E_abs_ref_interpol(J)=interpolate(RP02_E_abs_ref_temp,t(J));
    RP03_E_abs_ref_interpol(J)=interpolate(RP03_E_abs_ref_temp,t(J));
  end
  clear RP01_E_abs_ref_temp RP02_E_abs_ref_temp RP03_E_abs_ref_temp;
  E_diff_percent_uncorrected=zeros(length(RP01_E_abs),3);
  %max(RP01_E_abs_ref_interpol)
  %max(RP02_E_abs_ref_interpol)
  %max(RP03_E_abs_ref_interpol)
  E_diff_percent_uncorrected(:,1)=abs(  RP01_E_abs - RP01_E_abs_ref_interpol) ./ max(RP01_E_abs_ref_interpol);
  E_diff_percent_uncorrected(:,2)=abs(  RP02_E_abs - RP02_E_abs_ref_interpol) ./ max(RP02_E_abs_ref_interpol);
  E_diff_percent_uncorrected(:,3)=abs(  RP03_E_abs - RP03_E_abs_ref_interpol) ./ max(RP03_E_abs_ref_interpol);
  
  E_diff_percent_uncorrected(:,1)=abs(  RP01_E_abs - RP01_E_abs_ref_interpol) ./ (RP01_E_abs_ref_interpol);
  E_diff_percent_uncorrected(:,2)=abs(  RP02_E_abs - RP02_E_abs_ref_interpol) ./ (RP02_E_abs_ref_interpol);
  E_diff_percent_uncorrected(:,3)=abs(  RP03_E_abs - RP03_E_abs_ref_interpol) ./ (RP03_E_abs_ref_interpol);
  
  
  E_diff_percent_eps_uncorrected= 1/3*( max(E_diff_percent_uncorrected(:,1)) + max(E_diff_percent_uncorrected(:,2)) + max(E_diff_percent_uncorrected(:,3)) ) ;
  disp(['E_diff_percent_eps_uncorrected = ' num2str(E_diff_percent_eps_uncorrected)]);
  
  % Plot E_diff_percent_uncorrected: Abweichung zwischen dem aktuellen Fall und der Referenz
  f=figure;hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('(E_{abs}-E_{abs,ref})/max(E_{abs,ref}) [-]');
  set(gca, 'ColorOrder', [0 0 1]);       E_diff1=plot(t*1e9,E_diff_percent_uncorrected(:,1),'-');
  set(gca, 'ColorOrder', [1 0 0]);       E_diff2=plot(t*1e9,E_diff_percent_uncorrected(:,2),'-');
  set(gca, 'ColorOrder', [0 0.4706 0]);  E_diff3=plot(t*1e9,E_diff_percent_uncorrected(:,3),'-');
  E_diff_max=plot([0 t(end)]*1e9,[1 1]*E_diff_percent_eps_uncorrected,'--b');
  legend([E_diff1(1),E_diff2(1),E_diff3(1),E_diff_max],'Record Point 1: Deviation [%]','Record Point 2: Deviation [%]','Record Point 3: Deviation [%]','Averaged Deviation [%]','Location','NorthOutside');
  export_fig([casefolder '/Probes/E_abs_RP_diff_percent_uncorrected.pdf']);
  close(gcf);
  
  % remove possible spikes from E_diff_percent_uncorrected
  E_diff_percent= slopeCheckSpike(t, E_diff_percent_uncorrected, casefolder);
  E_diff_percent_eps= 1/3*( max(E_diff_percent(:,1)) + max(E_diff_percent(:,2)) + max(E_diff_percent(:,3)) ) ;
  disp(['E_diff_percent_eps = ' num2str(E_diff_percent_eps)]);
  save([pwd '/Error.mat'],'zeta','order','E_diff_percent_eps','PML','time','shape','chi','PMLcells','DOF','PMLspread','frequency');
  
   % Plot E_diff_percent: Abweichung zwischen dem aktuellen Fall und der Referenz
  f=figure;hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('(E_{abs}-E_{abs,ref})/max(E_{abs,ref}) [-]');
  set(gca, 'ColorOrder', [0 0 1]);       E_diff1=plot(t*1e9,E_diff_percent(:,1),'-');
  set(gca, 'ColorOrder', [1 0 0]);       E_diff2=plot(t*1e9,E_diff_percent(:,2),'-');
  set(gca, 'ColorOrder', [0 0.4706 0]);  E_diff3=plot(t*1e9,E_diff_percent(:,3),'-');
  E_diff_max=plot([0 t(end)]*1e9,[1 1]*E_diff_percent_eps,'--b');
  legend([E_diff1(1),E_diff2(1),E_diff3(1),E_diff_max],'Record Point 1: Deviation [%]','Record Point 2: Deviation [%]','Record Point 3: Deviation [%]','Averaged Deviation [%]','Location','NorthOutside');
  export_fig([casefolder '/Probes/E_abs_RP_diff_percent.pdf']);
  close(gcf);
  
  % Plot E and E_ref
  f=figure;hold on;  set(f, 'color', 'white');  xlabel('t [ns]'); ylabel('E [V/m]');
  set(gca, 'ColorOrder', [0 0 1]);       E1_abs=plot(t*1e9,RP01_E_abs,'-');  E1_abs_ref=plot(t*1e9,RP01_E_abs_ref_interpol,'-k');
  set(gca, 'ColorOrder', [1 0 0]);       E2_abs=plot(t*1e9,RP02_E_abs,'-');  E2_abs_ref=plot(t*1e9,RP02_E_abs_ref_interpol,'-k');
  set(gca, 'ColorOrder', [0 0.4706 0]);  E3_abs=plot(t*1e9,RP03_E_abs,'-');  E3_abs_ref=plot(t*1e9,RP03_E_abs_ref_interpol,'-k');
  legend([E1_abs(1),E2_abs(1),E3_abs(1),E1_abs_ref(1),E2_abs_ref(1),E3_abs_ref(1)],'Probe 1: E_{abs}','Probe 2: E_{abs}','Probe 3: E_{abs}','Probe 1: E_{abs,ref}','Probe 2: E_{abs,ref}','Probe 3: E_{abs,ref}','Location','EastOutside');
  export_fig([casefolder '/Probes/E_abs_RP_and_Reference.pdf']);
  close(gcf);
end
%% ========================================================================
% Frequencies
%==========================================================================
% Wenn Datei existiert: 'TimeEx','Frequency_Ex','TimeEy','Frequency_Ey','TimeEtheta','Frequency_Etheta'
if  exist([casefolder '/Probes/DominantFrequency_RP01.dat'],'file')==2 
  filename=([casefolder '/Probes/DominantFrequency_RP01.dat']);
  fid = fopen(filename,'r');
  Frequency_RP01 = textscan(fid,'%f %f %f %f %f %f', 'HeaderLines',3);
  filename=([casefolder '/Probes/DominantFrequency_RP02.dat']);
  fid = fopen(filename,'r');
  Frequency_RP02 = textscan(fid,'%f %f %f %f %f %f', 'HeaderLines',3);
  filename=([casefolder '/Probes/DominantFrequency_RP03.dat']);
  fid = fopen(filename,'r');
  Frequency_RP03 = textscan(fid,'%f %f %f %f %f %f', 'HeaderLines',3);
  
  f=figure;  set(f, 'Units', 'normalized', 'Position', [0.2, 0.2, 0.7, 0.7]); 
  set(f, 'color', 'white');
  xlabel('t [ns]'); ylabel('F [1/s]');
  set(gca, 'ColorOrder', [0 0 1]);
  F1_Ex=semilogy(Frequency_RP01{1,1}*1e9,Frequency_RP01{1,2},'o-');hold on;
  F1_Ey=semilogy(Frequency_RP01{1,3}*1e9,Frequency_RP01{1,4},'s-');
  
  set(gca, 'ColorOrder', [1 0 0]);
  F2_Ex=semilogy(Frequency_RP02{1,1}*1e9,Frequency_RP02{1,2},'o-');
  F2_Ey=semilogy(Frequency_RP02{1,3}*1e9,Frequency_RP02{1,4},'s-');
  
  set(gca, 'ColorOrder', [0 0.4706 0]);
  F3_Ex=semilogy(Frequency_RP03{1,1}*1e9,Frequency_RP03{1,2},'o-');
  F3_Ey=semilogy(Frequency_RP03{1,3}*1e9,Frequency_RP03{1,4},'s-');
  
  F_analytical=semilogy([0 input.Dipole{strcmp(input.Dipole,'tend'),2}]*1e9,[1 1]*input.frequency,'k:','Linewidth',2);
  if ~isempty(F1_Ex)
    legend([F1_Ex(1),F1_Ey(1),F2_Ex(1),F2_Ey(1),F3_Ex(1),F3_Ey(1),F_analytical(1)],...
      'Probe 1: Frequency E_{x}','Probe 1: Frequency E_{y}',...%'Probe 1: Frequency E_{\theta}',...%
      'Probe 2: Frequency E_{x}','Probe 2: Frequency E_{y}',...%'Probe 2: Frequency E_{\theta}',...
      'Probe 3: Frequency E_{x}','Probe 3: Frequency E_{y}',...%'Probe 3: Frequency E_{\theta}',...
      ['Dipole frequency f = ' num2str(input.frequency/1E6) ' MHz'],...
      'Location','EastOutside'); grid on;
    export_fig([casefolder '/Probes/Frequencies.pdf']);
    disp(['Frequency at end for E_x = [ ' num2str(Frequency_RP01{1,2}(end)/(1E6)) ' ' num2str(Frequency_RP02{1,2}(end)/(1E6)) ' ' num2str(Frequency_RP03{1,2}(end)/(1E6)) ' ] MHz'])
    disp(['Frequency at end for E_y = [ ' num2str(Frequency_RP01{1,4}(end)/(1E6)) ' ' num2str(Frequency_RP02{1,4}(end)/(1E6)) ' ' num2str(Frequency_RP03{1,4}(end)/(1E6)) ' ] MHz'])
  end
  close(gcf);
end
  
  end