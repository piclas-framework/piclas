function E_diff_percent_corrected = slopeCheckSpike(t, f, casefolder) 
slope_i1=zeros(length(t),3); slope_i2=zeros(length(t),3);
barrier = 1E7;
for K=1:3
  for I=2:length(t)-1
    
    x_i_minus_1 = t(I-1); 
    x_i = t(I); 
    x_i_plus_1 = t(I+1);
    
    y_i_minus_1 = f(I-1,K); 
    y_i = f(I,K); 
    y_i_plus_1 = f(I+1,K);
    
    slope_i1(I,K) = abs((y_i-y_i_minus_1)/(x_i-x_i_minus_1));
    slope_i2(I,K) = abs((y_i_plus_1-y_i)/(x_i_plus_1-x_i));
     
    if (slope_i1(I,K) > barrier) && (slope_i2(I,K) > barrier)
      disp('Spike detected !!!!! ')
      f(I,K)=(y_i_minus_1+y_i_plus_1)/2;
    end
  end
end

fig=figure; set(fig, 'color', 'white'); xlabel('t [ns]'); ylabel('slope [-]'); set(fig, 'Units', 'normalized', 'Position', [0, 0, 0.4, 0.8]); MS=8; LW=2; MA='-';
set(gca, 'ColorOrder', [0 0 1]);       slope_i1_RP1=semilogy(t*1e9,slope_i1(:,1),'^-'); hold on; slope_i2_RP1=semilogy(t*1e9,slope_i2(:,1),'o-');
set(gca, 'ColorOrder', [1 0 0]);       slope_i1_RP2=semilogy(t*1e9,slope_i1(:,2),'^-'); slope_i2_RP2=semilogy(t*1e9,slope_i2(:,2),'o-');
set(gca, 'ColorOrder', [0 0.4706 0]);  slope_i1_RP3=semilogy(t*1e9,slope_i1(:,3),'^-'); slope_i2_RP3=semilogy(t*1e9,slope_i2(:,3),'o-');
barrier_plot = semilogy([0 max(t*1e9)],[1 1]*barrier,'k--');
legend([slope_i1_RP1 slope_i2_RP1 slope_i1_RP2 slope_i2_RP2 slope_i1_RP3 slope_i2_RP3 barrier_plot],...
  'Record Point 1: slope left','Record Point 1: slope right',...
  'Record Point 2: slope left','Record Point 2: slope right',...
  'Record Point 2: slope left','Record Point 3: slope right',...
  'Barrier slope spike detection',...
  'Location','NorthOutside');
ylim([1E-15 1E15])
export_fig([casefolder '/Probes/E_abs_RP_diff_percent_slope.pdf']);

close(gcf);
E_diff_percent_corrected=f;



end