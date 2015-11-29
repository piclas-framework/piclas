function plot_2D_clipped_Bezier(N,P,boundaries) 
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
% plot clipped bezier curve (without new control points)
global Bezier2d;global n_Ls;global n_Ls_inv;global projectedFace
N_plot=5;
xi=linspace(boundaries(1,1),boundaries(1,2),N_plot);
eta=linspace(boundaries(2,1),boundaries(2,2),N_plot);
Bezier2d=zeros(length(xi),length(eta),2);
for q=1:N
  for p=1:N
    h1=plot([P(p,q,1) P(p+1,q,1) P(p+1,q+1,1) P(p,q+1,1) P(p,q,1)],[P(p,q,2) P(p+1,q,2) P(p+1,q+1,2) P(p,q+1,2) P(p,q,2)],'k-','Parent',get(projectedFace, 'children'));
  end
end
for I=1:length(xi)
  for J=1:length(eta)
    for q=1:N+1
      for p=1:N+1
        %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
        Bezier2d(I,J,:)=Bezier2d(I,J,:)+B(N,p-1,xi(I))*B(N,q-1,eta(J))*P(p,q,:);
      end
    end
  end
end
plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'r-',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r-','Parent',get(projectedFace, 'children'));
end