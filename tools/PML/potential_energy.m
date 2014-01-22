function potential_energy(input,casefolder)

filename=([casefolder '/Database.csv']);
fid = fopen(filename,'r');
switch input.ProjectName
  case 'Dipole'
    data = textscan(fid,'%f %s %f %s %f', 'HeaderLines',1);
  case 'SingleParticle'
    data = textscan(fid,'%f %s %f %s %f %s %f', 'HeaderLines',1);
    E_kin = data{1,7};
end
t     = data{1,1};
W_el  = data{1,3};
W_mag = data{1,5};
W_abs = W_el+W_mag;

f=figure;hold on; set(f, 'color', 'white'); xlabel('t [ns]'); ylabel('W [J]');


set(gca, 'ColorOrder', [0 0 1]);      W1=plot(t*1e9,W_el,'.-');
  strValuesel = strtrim(cellstr([num2str(W_el(end))])); %,'(%d,%d)'
  %text(t(end)*1e9*0.95,W_el(end),strValues, 'horizontal','center', 'vertical','top');
  
set(gca, 'ColorOrder', [1 0 0]);      W2=plot(t*1e9,W_mag,'.-');
  strValuesmag = strtrim(cellstr([num2str(W_mag(end))])); %,'(%d,%d)'
  %text(t(end)*1e9*0.95,W_mag(end),strValues, 'horizontal','center', 'vertical','top');
  
set(gca, 'ColorOrder', [0 0.4706 0]); W3=plot(t*1e9,W_abs,'.-');
  strValues = strtrim(cellstr(['W_{abs} = ' num2str(W_abs(end))])); %,'(%d,%d)'
  text(t(end)*1e9*0.95,W_abs(end)*1.07,strValues, 'horizontal','center', 'vertical','top');
  
  
legend([W1 W2 W3],['W_{el} = ' num2str(W_mag(end))],['W_{mag} = ' num2str(W_mag(end))],['W_{abs} = ' num2str(W_abs(end))],'Location','EastOutside');
export_fig([casefolder '/Probes/W_el_mag.pdf']);
close(gcf);

if exist('E_kin','var')
  f=figure;
  set(f, 'color', 'white');
  xlabel('t [ns]'); ylabel('W [J]');
  set(gca, 'ColorOrder', [0 0 1]); E_kin_plot=semilogy(t*1e9,E_kin,'.-'); hold on;
  legend([E_kin_plot],'E_{kin}','Location','EastOutside');
  ylim([0 max(E_kin(end),1E100)]);
  
  strValues = strtrim(cellstr(num2str(E_kin(end)))); %,'(%d,%d)'
  text(t(end)*1e9*0.95,E_kin(end),strValues, 'horizontal','center', 'vertical','top');
  
  export_fig([casefolder '/Probes/E_kin.pdf']);
  close(gcf);
end

end