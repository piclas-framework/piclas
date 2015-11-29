function PlotConvergence(ClipConvergenceXiEta,MaxClipIterations,Tolerance)
% plot convergence rate - tolerance
plotClipConvergenceXiEta=figure;  hold on;set(gcf, 'color', 'white');
plotxi=plot(ClipConvergenceXiEta(:,1),'rs--');
ploteta=plot(ClipConvergenceXiEta(:,2),'bo:');
plottol=plot([0 MaxClipIterations],Tolerance*[1 1],'k--');
set(gca, 'yscale', 'log')
title('convergence rate of the xi/eta clipping algorithm')
legend([plotxi ploteta plottol],...
  'xi-convergence',...
  'eta-convergence',...
  'break tolerance')
ylim([Tolerance/1000 1])


% plot convergence rate: smax-smin=? (X%)
plotClipConvergenceXiEta2=figure;  hold on;set(gcf, 'color', 'white');
plotxi=plot(ClipConvergenceXiEta(ClipConvergenceXiEta(:,3)>0,3),'rs--');
ploteta=plot(ClipConvergenceXiEta(ClipConvergenceXiEta(:,4)>0,4),'bo:');
title('convergence rate of the xi/eta clipping algorithm')
legend([plotxi ploteta ],...
  'xi-convergence',...
  'eta-convergence')
%ylim([Tolerance/10 1])


end