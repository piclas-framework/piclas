function P_elevated=ElevateBezierSurface(N,index,P,elevation,elevation_direction,ElevationMatrix)
% ==================================================================
% Description
% ==================================================================
%

% ==================================================================
% Input
% ==================================================================
%                   ::

% ==================================================================
% Output
% ==================================================================
%                   ::

% ==================================================================
% Function Start
% ==================================================================
global BezierSurfaceFigure
p_elevated=N+elevation;
temp=zeros(p_elevated+1,N+1,3);
P_elevated   =zeros(p_elevated+1,p_elevated+1,3);
switch elevation_direction
  case 1
    %%
    for dim=1:3
      for j=1:N+1 % for each N+1 line of the original patch
        [~,temp(:,j,dim)]=ElevateBezierPolynomial(N,P(:,j,dim,index),elevation,ElevationMatrix);
        %plotCP_elevated=plot(Xi_NGeo_elevated,P_elevated,'^m--');
      end
    end
    
    for dim=1:3
      for i=1:p_elevated+1
        [~,P_elevated(i,:,dim)]=ElevateBezierPolynomial(N,temp(i,:,dim),elevation,ElevationMatrix);
        %plotCP_elevated=plot(Xi_NGeo_elevated,P_elevated,'^m--');
      end
    end
  case 2
    %%
    for dim=1:3
      for j=1:N+1 % for each N+1 line of the original patch
        [~,temp(:,j,dim)]=ElevateBezierPolynomial(N,P(j,:,dim,index),elevation,ElevationMatrix);
        %plotCP_elevated=plot(Xi_NGeo_elevated,P_elevated,'^m--');
      end
    end
    for dim=1:3
      for i=1:p_elevated+1
        [~,P_elevated(i,:,dim)]=ElevateBezierPolynomial(N,temp(i,:,dim),elevation,ElevationMatrix);
        %plotCP_elevated=plot(Xi_NGeo_elevated,P_elevated,'^m--');
      end
    end
end

%% plotting
for j=1:p_elevated
  for i=1:p_elevated
    PlotWirframe(P_elevated(i:i+1,j:j+1,:))
  end
end
for j=1:N+1
  for i=1:p_elevated+1
    plot3(temp(i,j,1),temp(i,j,2),temp(i,j,3),'r^','MarkerSize',10,'LineWidth',2)
  end
end
for j=1:p_elevated+1
  for i=1:p_elevated+1
    plot3(P_elevated(i,j,1),P_elevated(i,j,2),P_elevated(i,j,3),'k+','MarkerSize',10,'LineWidth',2)
  end
end
view(66,42)
%export_fig(BezierSurfaceFigure,['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/3D_plot_elevated_direction' num2str(elevation_direction) '_export_fig.pdf']);
% ==================================================================
% Function End
% ==================================================================
end