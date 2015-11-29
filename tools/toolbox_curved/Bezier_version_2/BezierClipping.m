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

% switch method
%   case 0
%     result=zeros(4,2);
%     for L=1:4
%       %L=6
%       sol=zeros(N,N,2);
%       N2=N;
%       BezierControlPoints2D_temp=BezierControlPoints2D;
%       switch L
%         case 1
%           u=1/2*(umin+1);v=1/2*(vmin+1);
%           %u=0;v=0;
%         case 2
%           %u=umax;v=vmin;
%           u=1/2*(umax+1);v=1/2*(vmin+1);
%           %u=1;v=0;
%         case 3
%           %u=umax;v=vmax;
%           u=1/2*(umax+1);v=1/2*(vmax+1);
%           %u=0;v=1;
%         case 4
%           %u=umin;v=vmax;
%           u=1/2*(umin+1);v=1/2*(vmax+1);
%           %u=1;v=1;
%         case 5
%           u=1/2*((umin+umax)/2+1);v=1/2*(vmax+1);
%         case 6
%           u=1/2*(umin+1);v=1;
%       end
%       for I=1:N %de casteljau schritt
%         for q=1:N2
%           for p=1:N2
%             A=[BezierControlPoints2D_temp(p  ,q,:) BezierControlPoints2D_temp(p  ,q+1,:);...
%               BezierControlPoints2D_temp(p+1,q,:) BezierControlPoints2D_temp(p+1,q+1,:)];
%             %         A=[BezierControlPoints2D_temp(p  ,q,:) BezierControlPoints2D_temp(p+1  ,q,:);...
%             %           BezierControlPoints2D_temp(p,q+1,:) BezierControlPoints2D_temp(p+1,q+1,:)];
%             for K=1:2
%               sol(p,q,K)=[1-u u]*A(:,:,K)*[1-v;v];
%             end
%           end
%         end
%         N2=N2-1;
%         BezierControlPoints2D_temp=sol;
%         sol=zeros(N2,N2,2);
%         %
%         %axes1 = axes('Parent',projectedFace)
%         plot(BezierControlPoints2D_temp(:,:,1),BezierControlPoints2D_temp(:,:,2),'bo','Parent',get(projectedFace, 'children'))
%       end
%       result(L,:)=BezierControlPoints2D_temp;
%       plot(result(L,1),result(L,2),'sr-','Parent',get(projectedFace, 'children'))
%     end
%     patch=[result;result(1,:)];
%     plot(patch(:,1),patch(:,2),'r-','Parent',get(projectedFace, 'children'))
%
%   case 1
%%

switch direction
  case 1 % u direction
    color='k+';
    Bezier2d=zeros(N+1,N+1,2);
    Bezier2d_temp=BezierControlPoints2D;
    %plot(BezierControlPoints2D(:,:,1),BezierControlPoints2D(:,:,2),'ro',BezierControlPoints2D(:,:,1)',BezierControlPoints2D(:,:,2)','ro','Parent',get(projectedFace, 'children'));
    %% 1. top
    s=smax;
    for J=1:N+1
      for K=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          Bezier2d(K,J,:)=Bezier2d(K,J,:)+Bezier2d_temp(L,J,:)*B(K-1,L-1,s);
        end
      end
    end
    %plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'ks',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r+','Parent',get(projectedFace, 'children'));
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'rs',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r+','Parent',get(projectedFace, 'children'));
    %% 2. bottom
    %  Bezier2d_temp(:,:,1)=transpose(Bezier2d(:,:,1));% altes array mit neuen punkten überschreiben für 2. clip
    %   Bezier2d_temp(:,:,2)=transpose(Bezier2d(:,:,2));
    Bezier2d_temp=Bezier2d;
    Bezier2d=zeros(N+1,N+1,2);
    %s=smin/smax; [0, 1]
    s=(smin+1)/(smax+1); %[-1, +1]
    s=2*(1-s)-1;
    for J=1:N+1
      for K=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          Bezier2d(K,J,:)=Bezier2d(K,J,:)+Bezier2d_temp(N+2-L,J,:)*B(K-1,L-1,s);
        end
      end
    end
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'go',Bezier2d(:,:,1)',Bezier2d(:,:,2)','g+','Parent',get(projectedFace, 'children'));
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
    for I=1:N+1
      for K=1:N+1
        for L=1:K
          %plot(P(p,q,1),P(p,q,2),'go','MarkerSize',10,'LineWidth',5)
          Bezier2d(I,K,:)=Bezier2d(I,K,:)+Bezier2d_temp(I,L,:)*B(K-1,L-1,s);
        end
      end
    end
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'rs',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r+','Parent',get(projectedFace, 'children'));
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
    %plot(Bezier2d(:,:,1),Bezier2d(:,:,2),color,Bezier2d(:,:,1)',Bezier2d(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'go',Bezier2d(:,:,1)',Bezier2d(:,:,2)','g+','Parent',get(projectedFace, 'children'));
    BezierClipping=Bezier2d;
    
end


if and(K==1,true(createPlots))
  export_fig(projectedFace,['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/2D_projected_bezier_' num2str(J) '_clip_export_fig.pdf']);
end
end