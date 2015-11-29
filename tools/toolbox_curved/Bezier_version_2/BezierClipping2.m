function BezierClipping=BezierClipping(direction,smin,smax,N,BezierControlPoints2D,method,createPlots)
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
%%
global projectedFace;

switch direction
  case 1 % u direction
    color='k+';
    Bezier2d=zeros(N+1,N+1,2);
    Bezier2d_temp=BezierControlPoints2D;
    %% 1. cut off the top part [-1,+1] -> [-1,smax]
    s=smax;
    for K=1:N+1
      for J=1:N+1
        for L=1:K
          Bezier2d(K,J,:)=Bezier2d(K,J,:)+Bezier2d_temp(L,J,:)*B(K-1,L-1,s);
        end
      end
    end
    % plot first cufoff
    %plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'rs',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r+','Parent',get(projectedFace, 'children'));
    
    
    % Bedingungen: tut!
    % Bezier2d(1,1,:) == [-0.2838, -0.3946]T <-- bleibt bei yi=-1
    % Bezier2d(3,3,:) == [ 0.0918,  0.1031]T <-- smax verschiebung
    %% 1. cut off the bottom part [-1,smax] -> [smin,smax]
    Bezier2d_temp=Bezier2d;
    Bezier2d=zeros(N+1,N+1,2); % Bezier2d=0 geht nicht !!!!!
    
    s=(smin+1)/(smax+1); %[-1, +1]
    s=2*(1-s)-1;
    for J=1:N+1
      for K=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          %disp(['N+2-L = ' num2str(N+2-L)])
          Bezier2d(K,J,:)=Bezier2d(K,J,:)+Bezier2d_temp(N+2-L,J,:)*B(K-1,L-1,s);
          %plot(Bezier2d(K,J,1),Bezier2d(K,J,2),'go','MarkerSize',10,'LineWidth',5)
        end
      end
    end
    
    % plot second cufoff
    % plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'go',Bezier2d(:,:,1)',Bezier2d(:,:,2)','g+','Parent',get(projectedFace, 'children'));
    
  %% flip the Xi-direction (the first index. it is somehow wronly directed)
    Bezier2d_flipped=Bezier2d; % create temp matrix for flipping
    for i=1:N+1
      Bezier2d_flipped(i,:,1)=Bezier2d(N+2-i,:,1);
      Bezier2d_flipped(i,:,2)=Bezier2d(N+2-i,:,2);
    end
    Bezier2d=Bezier2d_flipped;
    
    % plot second cufoff
    %plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'go',Bezier2d(:,:,1)',Bezier2d(:,:,2)','g+','Parent',get(projectedFace, 'children'));

    % Bedingungen:
    % Bezier2d(1,1,:) == [-0.3524, -0.09689]T  <-- smin verschiebung
    % Bezier2d(3,3,:) == [ 0.0918,  0.1031]T   <-- smax verschiebung
    
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),color,Bezier2d(:,:,1)',Bezier2d(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    BezierClipping=Bezier2d;
    %BezierClipping=BezierControlPoints2D

    
  case 2 % v direction
    color='rs';
    %         BezierControlPoints2D(:,:,1)=transpose(BezierControlPoints2D(:,:,1));
    %         BezierControlPoints2D(:,:,2)=transpose(BezierControlPoints2D(:,:,2));
    Bezier2d=zeros(N+1,N+1,2);
    Bezier2d_temp=BezierControlPoints2D;
    %plot(BezierControlPoints2D(:,:,1),BezierControlPoints2D(:,:,2),'ro',BezierControlPoints2D(:,:,1)',BezierControlPoints2D(:,:,2)','ro','Parent',get(projectedFace, 'children'));
    %% 1. top
    s=smax;
    for K=1:N+1
      for I=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          Bezier2d(I,K,:)=Bezier2d(I,K,:)+Bezier2d_temp(I,L,:)*B(K-1,L-1,s);
        end
      end
    end
    %plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'ks',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r+','Parent',get(projectedFace, 'children'));
    %% 2. bottom
    %  Bezier2d_temp(:,:,1)=transpose(Bezier2d(:,:,1));% altes array mit neuen punkten überschreiben für 2. clip
    %   Bezier2d_temp(:,:,2)=transpose(Bezier2d(:,:,2));
    Bezier2d_temp=Bezier2d;
    Bezier2d=zeros(N+1,N+1,2);
    %s=smin/smax; [0, 1]
    s=(smin+1)/(smax+1); %[-1, +1]
    s=2*(1-s)-1;
    for I=1:N+1
      for K=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          Bezier2d(I,K,:)=Bezier2d(I,K,:)+Bezier2d_temp(I,N+2-L,:)*B(K-1,L-1,s);
        end
      end
    end
    
    %% flip the Eta-direction (the first index. it is somehow wronly directed)
    Bezier2d_flipped=Bezier2d; % create temp matrix for flipping
    for i=1:N+1
      Bezier2d_flipped(:,i,1)=Bezier2d(:,N+2-i,1);
      Bezier2d_flipped(:,i,2)=Bezier2d(:,N+2-i,2);
    end
    Bezier2d=Bezier2d_flipped;
    
    
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),color,Bezier2d(:,:,1)',Bezier2d(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    BezierClipping=Bezier2d;
    
end


if and(K==1,true(createPlots))
  export_fig(projectedFace,['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/2D_projected_bezier_' num2str(J) '_clip_export_fig.pdf']);
end
end