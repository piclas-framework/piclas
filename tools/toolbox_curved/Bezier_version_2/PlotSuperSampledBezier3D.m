function PlotBezier3D=PlotSuperSampledBezier3D(N,P,SideID,sampPts,dimension)
% ==================================================================
% Description
% ==================================================================
% 3D plot of the bezier surface: evaulates the bernstein polynomial numerous times ! costly? yes, but yolo !
% 1.) plot3 plot: control points and connectivity in 3D
% 2.) surf plot : the compelte bezier surface polynomial

% ==================================================================
% Input
% ==================================================================
% sampPts           ::

% ==================================================================
% Output
% ==================================================================
%                   ::

% ==================================================================
% Function Start
% ==================================================================
switch dimension
  case 2
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
    % plot wireframe
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'b-',Bezier2d(:,:,1)',Bezier2d(:,:,2)','b-');

    axis equal
  case 3
    global ParticleBezierSurface Face
    xi=linspace(-1,1,sampPts);eta=linspace(-1,1,sampPts);
    Bezier=zeros(length(xi),length(eta),3);

    
    %% 1.) super sample the the Bezier polynomial and
    %  2.) plot the control points
    for I=1:length(xi)
      for J=1:length(eta)
        for q=1:N+1
          for p=1:N+1
            %disp(['p=' num2str(p) ' q=' num2str(q) ])
            %disp(['n=' num2str(N) ' k=' num2str(p-1) '(p)'])
            % disp(['n=' num2str(N) ' k=' num2str(q-1) '(q)'])
            plot3(P(p,q,1),P(p,q,2),P(p,q,3),'bo','MarkerSize',10,'LineWidth',2)
            Bezier(I,J,:)=Bezier(I,J,:)+B(N,p-1,xi(I))*B(N,q-1,eta(J))*P(p,q,:);
          end
        end
      end
    end
    
    %% plot the connectivity of the control points (plots a mesh)
    %h1=surf(P(:,:,1),P(:,:,2),P(:,:,3));set(h1,'FaceAlpha',0);
    for j=1:N
      for i=1:N
        PlotWirframe(P(i:i+1,j:j+1,:))
      end
    end
    
    %% surf plot of the Bezier polynomial
    if SideID==Face
      % P_i,j = [Px Py Pz]^T
      % x-coordinates: Px = Bezier(:,:,1)
      % y-coordinates: Py = Bezier(:,:,2)
      % z-coordinates: Pz = Bezier(:,:,3)
      ParticleBezierSurface=surf(Bezier(:,:,1),Bezier(:,:,2),Bezier(:,:,3));
      set(ParticleBezierSurface,'FaceAlpha',0.9)
    else
      h=surf(Bezier(:,:,1),Bezier(:,:,2),Bezier(:,:,3));
      set(h,'FaceAlpha',0.25)
    end
    
  otherwise
    error('fail')
end
% ==================================================================
% Function End
% ==================================================================
end
