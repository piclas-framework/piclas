function notused(N,P) 
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
% 2D projected bezier surface
global Bezier2d
global n_Ls
global n_Ls_inv
N_plot=10;
xi=linspace(-1,1,N_plot);eta=linspace(-1,1,N_plot);
Bezier2d=zeros(length(xi),length(eta),2);
for q=1:N
  for p=1:N
    h1=plot([P(p,q,1) P(p+1,q,1) P(p+1,q+1,1) P(p,q+1,1) P(p,q,1)],[P(p,q,2) P(p+1,q,2) P(p+1,q+1,2) P(p,q+1,2) P(p,q,2)],'k-');
  end
end
for I=1:length(xi)
  for J=1:length(eta)
    for q=1:N+1
      for p=1:N+1
        % plot the control points of the 2D bezier polynomial
        plot(P(p,q,1),P(p,q,2),'bo','MarkerSize',10,'LineWidth',3)
        Bezier2d(I,J,:)=Bezier2d(I,J,:)+B(N,p-1,xi(I))*B(N,q-1,eta(J))*P(p,q,:);
      end
    end
  end
end
% plot origin as circle (the ray pierces the berzier surface at this position)
plot(0,0,'ko','MarkerSize',10,'LineWidth',5)
plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'b-',Bezier2d(:,:,1)',Bezier2d(:,:,2)','b-');
%% get Line L_s (vector n_Ls and n_Ls_inv)
axis equal
calcLsLines(P)
end