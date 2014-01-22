%function compare_energy %(input,casefolder)
clear ; clc; close all;
addpath('/home/stephen/Diplomarbeit/PMLalgorithm/export_fig');
particle = 'false';


for K = 1:3
  case_number = K;
  
  switch case_number
    % 20MHz - p7
    case 1
      % 1PML - 20MHz
      reference_case = '/home/stephen/Diplomarbeit/PMLalgorithm/Reference_Case/Mannheim/20MHz/domain018_nElems030_order07_PML10_zeta0E+00_20MHz          ';
      clean_case     = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/20MHz_1PML/domain007_nElems012_order07_PML01_zeta0E+00_20MHz      ';
      corrected_case = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/20MHz_1PML/domain007_nElems012_order07_PML01_zeta300000000_20MHz  ';
      name='pulsed_20MHz_1PML';
    case 2
      % 5PML - 20MHz
      reference_case = '/home/stephen/Diplomarbeit/PMLalgorithm/Reference_Case/Mannheim/20MHz/domain018_nElems030_order07_PML10_zeta0E+00_20MHz          ';
      clean_case     = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/20MHz_5PML/domain012_nElems020_order07_PML05_zeta0E+00_20MHz      ';
      corrected_case = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/20MHz_5PML/domain012_nElems020_order07_PML05_zeta125000000_20MHz  ';
      name='pulsed_20MHz_5PML';
    case 3
      % continuous
      reference_case = '/home/stephen/Diplomarbeit/PMLalgorithm/Reference_Case/Mannheim/continuous_dipole/domain018_nElems030_order07_PML10_zeta0E+00                       ';
      clean_case     = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/continuous_100MHz_1PML/domain007_nElems012_order07_PML01_zeta0E+00_continuous        ';
      corrected_case = '/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/pulsed/continuous_100MHz_1PML/domain007_nElems012_order07_PML01_zeta287500000_continuous    ';
      name='continuous_100MHz_1PML';
  end
  
  %%
  load([strtrim(corrected_case) '/Error.mat'])
  A=[reference_case;clean_case;corrected_case];
  for I=1:3
    disp(['Opening: ' strtrim(A(I,:)) '/Database.csv']);
    filename=([strtrim(A(I,:)) '/Database.csv']); fid = fopen(filename,'r');
    switch particle
      case 'false'
        data = textscan(fid,'%f %s %f %s %f', 'HeaderLines',1);
      case 'true'
        data = textscan(fid,'%f %s %f %s %f %s %f', 'HeaderLines',1);
        switch I
          case 1
            E_kin_1 = data{1,7};
          case 2
            E_kin_2 = data{1,7};
          case 3
            E_kin_3 = data{1,7};
        end
    end
    switch I
      case 1
        t_ref     = data{1,1};
        W_el_ref  = data{1,3};
        W_mag_ref = data{1,5};
        W_abs_ref = W_el_ref+W_mag_ref;
      case 2
        t_cle     = data{1,1};
        W_el_cle  = data{1,3};
        W_mag_cle = data{1,5};
        W_abs_cle = W_el_cle+W_mag_cle;
      case 3
        t_cor     = data{1,1};
        W_el_cor  = data{1,3};
        W_mag_cor = data{1,5};
        W_abs_cor = W_el_cor+W_mag_cor;
    end
  end
  
  f=figure; set(f, 'color', 'white'); xlabel(['t [ns]']); ylabel('W_{abs} / MAX(W_{abs,ref}) [J]'); hold on;
  normalize = max(W_abs_ref);
  set(gca, 'ColorOrder', [0 0 1]);      ref_case=plot(t_ref*1E9,W_abs_ref/normalize,'-','MarkerSize',3,'LineWidth',1);
  set(gca, 'ColorOrder', [1 0 0]);      cle_case=plot(t_cle*1E9,W_abs_cle/normalize,'--','MarkerSize',4,'LineWidth',1);
  set(gca, 'ColorOrder', [0 0.4706 0]); cor_case=plot(t_cor*1E9,W_abs_cor/normalize,'-.','MarkerSize',2,'LineWidth',1);
  
  legend([ref_case cle_case cor_case],...
    ['Reference'],...
    ['Clean'],...
    ['PML, \zeta = ' num2str(zeta/1E8) 'E8'],...
    'Location','North');
  axis square; legend boxoff;
  
  
  export_fig(['/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/' name '_dipole.pdf']);
  %ylim([1E-6 4E-3])
  
  % ylim([0.85 1])
  % axis normal
  % set(f, 'Units', 'normalized', 'Position', [0, 0, 0.3, 0.5]); MS=4; LW=0.75;
  % export_fig(['/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/' name '_dipole_zoom.pdf']);
  
  E_diff_percent_eps
  close(gcf);
  
  
  difference_clean     = calcDifferenceInterpolation([ t_cle W_abs_cle ], [ t_ref W_abs_ref ]);
  difference_corrected = calcDifferenceInterpolation([ t_cor W_abs_cor ], [ t_ref W_abs_ref ]);
  f=figure; set(f, 'color', 'white');
  set(gca, 'ColorOrder', [1 0 0]);      cle_case=plot(difference_clean(:,1)    *1E9,difference_clean(:,2)    /normalize,'--','MarkerSize',3,'LineWidth',1);hold on;
  set(gca, 'ColorOrder', [0 0.4706 0]); cor_case=plot(difference_corrected(:,1)*1E9,difference_corrected(:,2)/normalize,'.-','MarkerSize',4,'LineWidth',1);
  
  legend([cle_case cor_case],...
    ['Clean'],...
    ['PML, \zeta = ' num2str(zeta/1E8) 'E8'],...
    'Location','NorthWest');
  axis square; legend boxoff;  xlabel(['t [ns]']); ylabel('(W_{abs}-W_{abs,ref}) / MAX(W_{abs,ref}) [J]');
  %ylim([10^-25 10^0])
  export_fig(['/home/stephen/Diplomarbeit/Auswertung/Energy_comparison/' name '_dipole_difference.pdf']);
  close(gcf);
  
  
end












return
filename=([casefolder '/Database.csv']);


t     = data{1,1};
W_el  = data{1,3};
W_mag = data{1,5};
W_abs = W_el+W_mag;

f=figure;hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('W [J]');

set(gca, 'ColorOrder', [0 0 1]);      W1=plot(t*1e9,W_el,'.-');
set(gca, 'ColorOrder', [1 0 0]);      W2=plot(t*1e9,W_mag,'.-');
set(gca, 'ColorOrder', [0 0.4706 0]); W3=plot(t*1e9,W_abs,'.-');

legend([W1(1),W2(1),W3(1)],'W_{el}','W_{mag}','W_{ges}','Location','EastOutside');
export_fig([casefolder '/Probes/W_el_mag.pdf']);
close(gcf);

if exist('E_kin','var')
  f=figure;hold on;
  set(f, 'color', 'white');
  xlabel('t [ns]'); ylabel('W [J]');
  set(gca, 'ColorOrder', [0 0 1]); E_kin=plot(t*1e9,E_kin,'.-');
  legend([E_kin(1)],'E_{kin}','Location','EastOutside');
  export_fig([casefolder '/Probes/E_kin.pdf']);
  close(gcf);
end