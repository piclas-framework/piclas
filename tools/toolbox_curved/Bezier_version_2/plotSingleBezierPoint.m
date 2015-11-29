function plotSingleBezierPoint(N,P,Xi_Eta,dimension)
% ==================================================================
% Description
% ==================================================================
% 

% ==================================================================
% Input
% ==================================================================
% N                 :: polynomial degree
% P                 :: Control Points of the Bézier polynomial
% Xi_Eta(1)         :: xi-variable
% Xi_Eta(2)         :: eta-variable
% dimension         :: switch for 2D or 3D evaluation and plot 

% ==================================================================
% Output
% ==================================================================
%                   :: 

% ==================================================================
% Function Start
% ==================================================================
switch dimension
  case 2
    global projectedFace
    Bezier2d=0;
    %Xi_Eta(1)
    %Xi_Eta(2)
    for q=1:N+1
      for p=1:N+1
        Bezier2d=Bezier2d+B(N,p-1,Xi_Eta(1))*B(N,q-1,Xi_Eta(2))*P(p,q,:);
      end
    end
    %plot(Bezier2d(1),Bezier2d(2),'ro-','MarkerSize',20);
    plot(Bezier2d(1),Bezier2d(2),'bo-','MarkerSize',20,'Parent',get(projectedFace, 'children'));
  case 3
    global BezierSurfaceFigure
    Bezier3d=0;
    % Xi_Eta(1)=0
    % Xi_Eta(2)=0
    for q=1:N+1
      for p=1:N+1
        Bezier3d=Bezier3d+B(N,p-1,Xi_Eta(1))*B(N,q-1,Xi_Eta(2))*P(p,q,:);
      end
    end
    plot3(Bezier3d(1),Bezier3d(2),Bezier3d(3),'k+','MarkerSize',20,'Parent',get(BezierSurfaceFigure, 'children'));
  otherwise
    error('plotSingleBezierPoint(N,P,Xi_Eta,dimension): wrong dimension')
end

end
