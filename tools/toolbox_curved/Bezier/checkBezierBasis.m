function plot_cases
clc; clear all; close all;
global Bezier2d;global n_Ls;global n_Ls_inv;global projectedFace;global Face;global BezierSurface
NGeo=2;

load faces.dat
%load /home/stephen/PMLalgorithm_cases/Testing_Bezier_domain010_nElems002_order04_PML00_zeta0E+00_polynom_CFL0.5_N10_Parts500_DoPML_False/BezierControlPoints.dat
load BezierControlPoints.dat
%factorial(5)
faces=BezierControlPoints;
dimension = size(faces);

N=dimension(2)-1;
N=NGeo;
nSides=dimension(1)/(N+1)^2;
%%
P=zeros(N+1,N+1,3,nSides);
for i=1:nSides
  index1=(i-1)*(N+1)^2;
  for q=1:N+1
    for p=1:N+1
      %disp([(q-1)*(N+1)+p+index1])
      %faces((q-1)*(N+1)+p,:);
      P(p,q,:,i)=faces((q-1)*(N+1)+p+index1,:);
      %disp([P(:,p,q)])
    end
  end
end
%%

Tolerance=1E-2
%pvec=[0.4 0.2 0.001];v=[0 1 0];Face=9;
%pvec=[0.6 0.2 0.001];v=[1 -1 1];Face=9; %miss
%pvec=[0.6 0.05 0.2];v=[-1 1 -1];Face=9;pvec=pvec+v*0.1 %double
%pvec=[0.6 0.2 0.001];v=[1 -1 0];Face=9; %miss
pvec=[0.4 0.2 0.001];v=[0.7 0.9 1];Face=9;
%pvec=[0.4 0.2 0.001];v=[0 0 1];Face=10;
%pvec=[0.3 0.3 0.001];v=[1 1 0];Face=9;
%pvec=[0.3 0.3 0.001];v=[0 0 1];Face=10;
%pvec=[0.3 0.3 0.001];v=[0.2 -1 0.2];Face=8;
%%

BezierSurface=figure;  hold on;set(gcf, 'color', 'white');view(29,60);grid on;xlabel('x');ylabel('y');zlabel('z');
%plot3(faces(:,1),faces(:,2),faces(:,3),'ro','MarkerSize',10,'LineWidth',5)
Xi_NGeo=zeros(1,N+1);
for I=1:N+1
  Xi_NGeo(I) = 2/N*(I-1)-1;
end
count=0;
for i=1:nSides
  disp([num2str(i)])
  %if i==10
  if i==7||i==9||i==3||i==2
    index1=(i-1)*(N+1)^2+1;
    index2=i*(N+1)^2;
    %plot3(faces(index1:index2,1),faces(index1:index2,2),faces(index1:index2,3),'ro','MarkerSize',10,'LineWidth',5)
    g_h(N,P(:,:,:,i),i);
    count=count+1;
    index=i;
    %return
  end
end
if count==111
  x=P(:,:,1,index);y=P(:,:,2,index);z=P(:,:,3,index);
  %   b=max([max(x),max(y),max(z)])
  %   a=min([min(x),min(y),min(z)])
  axis equal
  percent=0.3;
  xlim([min(min(x))-(1-percent)*abs(min(min(x))) max(max(x))+(1+percent)*abs(max(max(x)))]);
  ylim([min(min(y))-(1-percent)*abs(min(min(y))) max(max(y))+(1+percent)*abs(max(max(y)))]);
  zlim([min(min(z))-(1-percent)*abs(min(min(z))) max(max(z))+(1+percent)*abs(max(max(z)))]);
  
end

%% particle

SideID=9;dimensions=size(P);
BezierControlPoints=P;
BezierControlPoints(:,:,:,SideID);

v=1/(v*v')*v;
p2=pvec+0.3*v;
plot3(pvec(1),pvec(2),pvec(3),'or','LineWidth',4,'MarkerSize',10)
plot3(p2(1),p2(2),p2(3),'or','LineWidth',4,'MarkerSize',10)
quiver3(pvec(1),pvec(2),pvec(3),v(1),v(2),v(3),0,'LineWidth',2,'MarkerSize',10)
view(35,26)
view(v) % view in direction of the vector
%return
%% projection 3D->2D
%1. create perpendicular vectors
if abs(v(3))<1E-6
  n1=[-v(2)-v(3),v(1),v(1)];
else
  n1=[v(3),v(3),-v(1)-v(2)];
end
n1=1/sqrt(n1*n1')*n1;
quiver3(pvec(1),pvec(2),pvec(3),n1(1),n1(2),n1(3),0,'LineWidth',2,'MarkerSize',10)
n2=cross(v,n1);
quiver3(pvec(1),pvec(2),pvec(3),n2(1),n2(2),n2(3),0,'LineWidth',2,'MarkerSize',10)
%dimensions=size(BezierControlPoints);
d1ij=zeros(dimensions(1),dimensions(2),dimensions(4));
d2ij=zeros(dimensions(1),dimensions(2),dimensions(4));
for SideID=1:nSides
  for q=1:NGeo+1
    for p=1:NGeo+1
      d1ij(p,q,SideID)=([BezierControlPoints(p,q,1,SideID) BezierControlPoints(p,q,2,SideID) BezierControlPoints(p,q,3,SideID)]-pvec)*n1';
      d2ij(p,q,SideID)=([BezierControlPoints(p,q,1,SideID) BezierControlPoints(p,q,2,SideID) BezierControlPoints(p,q,3,SideID)]-pvec)*n2';
    end
  end
end
%%
monitor='small';
set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');
position = get(gcf,'Position');
outerpos = get(gcf,'OuterPosition');
borders = outerpos - position;
edge = -borders(1)/2;
if strcmp(monitor,'large'),pos1=[outerpos(1)*5.5,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos2=[outerpos(1)*6.7,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos3=[outerpos(1)*5.5,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];
  pos4=[outerpos(1)*6.7,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];end;

if strcmp(monitor,'small'),pos1=[1400,550,scnsize(3)/6,scnsize(4)/2.2];
  pos2=[1950,550,scnsize(3)/6,scnsize(4)/2.2];
  pos3=[1400,50,scnsize(3)/7,scnsize(4)/2.2];
  pos4=[1800,50,scnsize(3)/7,scnsize(4)/2.2];end;
set(gcf,'OuterPosition',pos1)
axis equal
addpath(['/home/stephen/MATLAB_Programme/export_fig/']);
%export_fig([pwd '/3D_plot.pdf']);
%%
%pause(3)
BezierControlPoints2D=zeros(dimensions(1),dimensions(2),2,dimensions(4));
for SideID=1:nSides
  for q=1:NGeo+1
    for p=1:NGeo+1
      BezierControlPoints2D(p,q,:,SideID)=[d1ij(p,q,SideID) d2ij(p,q,SideID)];
    end
  end
end
projectedFace = figure;  hold on;set(gcf, 'color', 'white');set(gcf,'OuterPosition',pos2)
%g_h(N,P(:,:,:,i));
BezierControlPoints2Dorg=BezierControlPoints2D(:,:,:,Face);
g_h2D(N,BezierControlPoints2D(:,:,:,Face));
%%
BezierControlPoints1D=zeros(dimensions(1),dimensions(2),dimensions(4));
xi=linspace(-1,1,NGeo+1);eta=linspace(-1,1,NGeo+1);
minmax=zeros(2,NGeo+1)-2;
%disp(['n_Ls_inv*n_Ls=' num2str(n_Ls_inv*n_Ls') ])

boundaries=[-1 1; -1 1];
DoClipping=ones(2,1);
sarray=zeros(25,2,2);
sarray(:,:,1)=-1.0;
sarray(:,:,2)=1.0;
iiter=0;
jiter=0;
for K=1:5
  for J=1:2
    if true(DoClipping(J))
      switch J
        case 1
          iiter=iiter+1;
          vector = [n_Ls_inv(1) n_Ls_inv(2)];posX=pos3;
        case 2
          jiter=jiter+1;
          vector = [n_Ls(1) n_Ls(2)];posX=pos4;
      end
      f=figure;hold on;set(gcf, 'color', 'white');grid on;xlabel('x');ylabel('y');plot([-1 1],[0 0],'k-');set(gcf,'OuterPosition',posX)
      for SideID=1:nSides
        for q=1:NGeo+1
          for p=1:NGeo+1
            BezierControlPoints1D(p,q,SideID)=([BezierControlPoints2D(p,q,1,SideID) BezierControlPoints2D(p,q,2,SideID)])*vector';
          end
        end
      end
      if J==1
        for I=1:NGeo+1
          plot(xi,BezierControlPoints1D(:,I,Face),'ko-')
          minmax(2,I)=max(BezierControlPoints1D(I,:,Face)); % Upper
          minmax(1,I)=min(BezierControlPoints1D(I,:,Face)); % Lower
        end
      else
        for I=1:NGeo+1
          plot(xi,BezierControlPoints1D(I,:,Face),'ko-')
          minmax(2,I)=max(BezierControlPoints1D(:,I,Face)); % Upper
          minmax(1,I)=min(BezierControlPoints1D(:,I,Face)); % Lower
        end
      end
        
      smin= 1.5;
      smax=-1.5;
      %check streckenzug upper/lower
      for I=1:NGeo
        if  minmax(2,I)*minmax(2,I+1)<=0 % Upper
          m=(minmax(2,I+1)-minmax(2,I))/(xi(I+1)-xi(I));
          smin_proposed=xi(I)-minmax(2,I)/m;
          smin=min(smin_proposed,smin);
        end
        if  minmax(1,I)*minmax(1,I+1)<=0 % Upper
          m=(minmax(1,I+1)-minmax(1,I))/(xi(I+1)-xi(I));
          smax_proposed=xi(I)-minmax(1,I)/m;
          smax=max(smax_proposed,smax);
        end
      end
      %% streckenzug:
      if  minmax(2,1)*minmax(2,NGeo+1)<=0 %1. zeile anfang/ende -> upper
        m=(minmax(2,NGeo+1)-minmax(2,1))/(xi(NGeo+1)-xi(1));
        smin_proposed=xi(1)-minmax(2,1)/m;
        smin=min(smin_proposed,smin);
      end
      if  minmax(1,1)*minmax(1,NGeo+1)<=0 %2. zeile anfang/ende -> lower
        m=(minmax(1,NGeo+1)-minmax(1,1))/(xi(NGeo+1)-xi(1));
        smax_proposed=xi(1)-minmax(1,1)/m;
        smax=max(smax_proposed,smax);
      end
      %% streckenzug:
      if  minmax(1,1)*minmax(2,1)<=0 % 1 spalte anfang/ende
        smin_proposed=xi(1);
        smin=min(smin_proposed,smin);
      end
      if  minmax(1,2)*minmax(2,2)<=0 % 2 spalte anfang/ende
        smax_proposed=xi(NGeo+1);
        smax=max(smax_proposed,smax);
      end
      %%
      if smin== 1.5, smin=-1;end;
      if smax==-1.5, smax= 1;end;
      %       sarray(K,J,1)=smin;
      %       sarray(K,J,2)=smax;
      sarray(K,J,1)=smin;
      sarray(K,J,2)=smax;
      
      
      % disp([' k= ',num2str(K),' J= ',num2str(J),' smin= ',num2str(smin),' smax= ',num2str(smax)])
      plot([smin smax],[0 0],'r+')
      
      %%
      boundaries(J,:)=[smin smax];
      %% CLIP
      %       if J==2
      %         smax=0.5
      %         smin=-0.5
      %       end
      x_sol=mean(mean(BezierControlPoints2D(:,:,1,Face)));
      y_sol=mean(mean(BezierControlPoints2D(:,:,2,Face)));
      switch J
        case 1
          direction='u';
        case 2
          direction='v';
      end
      disp(['K=' num2str(K) ' J=' num2str(J) ' Clip ' direction '  smin=' sprintf('%1.2e',smin) ', smax=' sprintf('%1.2e',smax) ' smax-smin=' sprintf('%1.2e',smax-smin) ' (' sprintf('%3.2f',(smax-smin)/2*100) '%) delta=' sprintf('%1.2e',sqrt(x_sol^2+y_sol^2))])
      %if (smax-smin)/2*100>40,disp('gebiet muss geteilt werden!');end
      BezierControlPoints2D(:,:,:,Face)=BezierClipping(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2); % J={u,v}
      calcLsLines(BezierControlPoints2D(:,:,:,Face))
      %% check epsilon
      
      if abs(smax-smin)<1e-4
        DoClipping(J)=0;
      end
      %str=input('Press ENTER to continue or type "y" to exit: ','s');
      %if strcmp(str,'y'),error('You have terminated the program!');end;
      
      if sqrt(x_sol^2+y_sol^2)<Tolerance
        disp(['sqrt(x_sol^2+y_sol^2)<' num2str(Tolerance) ': origin found!'])
        plot(x_sol,y_sol, 'g+','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        plot(x_sol,y_sol, 'go','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        % compute the correct smin and smax value
        
        umean=0.5*sarray(K,1,1)+0.5*sarray(K,1,2);
        vmean=0.5*sarray(K,2,1)+0.5*sarray(K,2,2);
        umean=0.;
        vmean=0.;
%         disp([' umean: ',num2str(umean),' vmean: ',num2str(vmean)])
        for i=K:-1:1
          [umean,vmean]=LinIntPol2D(sarray(i,1,1),sarray(i,1,2),...
                                    sarray(i,2,1),sarray(i,2,2),...
                                    umean,vmean);
        end
        [umean,vmean]=LinIntPol2D(-1.,1.,...
                                  -1.,1.,...
                                  umean,vmean);
        disp([' umean: ',num2str(umean),' vmean: ',num2str(vmean)])
         smean=zeros(2,1);
         for L=1:2
           switch L
             case 1
               niter=iiter;
             case 2
               niter=jiter;
           end
           smin=sarray(niter,L,1);
           smax=sarray(niter,L,2);
           smean(L)=0.5*smin+0.5*smax;
           %smean(L)=0.; % center of deepest patch
           for i=niter-1:-1:1
              %disp(['i ',num2str(i)]);
              smin=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smin);
              smax=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smax);
              smean(L)=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smean(L));
           end
           %smean(L)=LinIntPol1D(-1,1,smean(L));
           %smean(L)=0.5*smin+0.5*smax;
           %          smean(J)=0.5*smin+0.5*smax;
           %smean(J)=0.5*(smin+smax);
           disp(['L: ',num2str(L),' smin: ',num2str(smin),' smax: ',num2str(smax),' mean: ',num2str(smean(L))])
         end
            
        %           smean(J)=0.5*smin+0.5*smax;
        %         for J=1:2
        %           switch J
        %             case 1
        %               niter=iiter;
        %             case 2
        %               niter=jiter;
        %           end
        %           smin=sarray(niter,J,1);
        %           smax=sarray(niter,J,2);
        %           smean(J)=0.5*smin+0.5*smax;
        %           for i=niter-1:-1:1
        % %             smin=LinIntPol1D(sarray(i,J,1),sarray(i,J,2),smin);
        % %             smax=LinIntPol1D(sarray(i,J,1),sarray(i,J,2),smax);
        %              smean(J)=LinIntPol1D(sarray(i,J,1),sarray(i,J,2),smean(J));
        %           end
        % %          smean(J)=0.5*smin+0.5*smax;
        %           %smean(J)=0.5*(smin+smax);
        %           disp(['J: ',num2str(J),' smin: ',num2str(smin),' smax: ',num2str(smax),' mean: ',num2str(smean(J))])
        %         end
        %plotBezierPoint2D(N,BezierControlPoints2Dorg,[0;0])
        plotBezierPoint2D(N,BezierControlPoints2Dorg,[smean(1);smean(2)])
        plotBezierPoint2D(N,BezierControlPoints2Dorg,[umean;vmean])
        %plotBezierPoint2D(N,BezierControlPoints2Dorg,[vmean;umean])
        %plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[0.;0.3])
        plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[smean(1);smean(2)])
        plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[umean;vmean])
        return
      end
    end
  end
  disp(' ')
end
plot(x_sol,y_sol, 'r+','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
plot(x_sol,y_sol, 'ro','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
plot_2D_clipped_Bezier(N,BezierControlPoints2D(:,:,:,Face),[-1 1; -1 1])



return
boundaries(:,1)=boundaries(:,1)-1E-4;
boundaries(:,2)=boundaries(:,2)+1E-4;

plot_2D_clipped_Bezier(N,BezierControlPoints2D(:,:,:,Face),boundaries);


%%
addpath(['/home/stephen/MATLAB_Programme/export_fig/']);
%export_fig([pwd '/2D_plot.pdf']);
end










function g_h=g_h(N,P,SideID) % 3D
global ParticleBezierSurface;global Face
N_plot=5;
xi=linspace(-1,1,N_plot);eta=linspace(-1,1,N_plot);
Bezier=zeros(length(xi),length(eta),3);
h1=surf(P(:,:,1),P(:,:,2),P(:,:,3));set(h1,'FaceAlpha',0)
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
%Bezier
if SideID==Face
  ParticleBezierSurface=surf(Bezier(:,:,1),Bezier(:,:,2),Bezier(:,:,3));
  set(ParticleBezierSurface,'FaceAlpha',0.9)
else
  h=surf(Bezier(:,:,1),Bezier(:,:,2),Bezier(:,:,3));
  set(h,'FaceAlpha',0.25)
end

end

function g_h2D=g_h2D(N,P) % 2D
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
        plot(P(p,q,1),P(p,q,2),'bo','MarkerSize',10,'LineWidth',3)
        Bezier2d(I,J,:)=Bezier2d(I,J,:)+B(N,p-1,xi(I))*B(N,q-1,eta(J))*P(p,q,:);
      end
    end
  end
end
plot(0,0,'ko','MarkerSize',10,'LineWidth',5)
plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'b-',Bezier2d(:,:,1)',Bezier2d(:,:,2)','b-');
%% get Line L_s (vector n_Ls and n_Ls_inv)
axis equal
calcLsLines(P)
end



function plotBezierPoint3D(N,P,boundaries)
global BezierSurface
Bezier3d=0;
% boundaries(1)=0
% boundaries(2)=0
for q=1:N+1
  for p=1:N+1
    Bezier3d=Bezier3d+B(N,p-1,boundaries(1))*B(N,q-1,boundaries(2))*P(p,q,:);
  end
end
plot3(Bezier3d(1),Bezier3d(2),Bezier3d(3),'ko-','MarkerSize',10,'Parent',get(BezierSurface, 'children'));
end

function plotBezierPoint2D(N,P,boundaries)
global projectedFace
Bezier2d=0;
%boundaries(1)
%boundaries(2)
for q=1:N+1
  for p=1:N+1
    Bezier2d=Bezier2d+B(N,p-1,boundaries(1))*B(N,q-1,boundaries(2))*P(p,q,:);
  end
end
plot(Bezier2d(1),Bezier2d(2),'bo-','MarkerSize',20,'Parent',get(projectedFace, 'children'));
end





function plot_2D_clipped_Bezier(N,P,boundaries) % plot clipped bezier curve (without new control points)
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

function B=B(n,k,x)
%disp(['n=' num2str(n) ' k=' num2str(k)])
B=(1/(2^n))*choose(n,k)*(x+1).^k.*(1-x).^(n-k);
end

function choose=choose(n,k)
if or(k==0,n==k)
  choose = 1;
else
  choose = factorial(n) / (factorial(k) * factorial(n-k));
end
end

function BezierClipping=BezierClipping(direction,smin,smax,N,BezierControlPoints2D,method)
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
    for I=1:N+1
      for K=1:N+1
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
    plot(Bezier2d(:,:,1),Bezier2d(:,:,2),color,Bezier2d(:,:,1)',Bezier2d(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    BezierClipping=Bezier2d;
    
end



end

function result=LinIntPol1D(y1,y2,x)
m=(y2-y1)/2;
result=m*(x+1)+y1;
end

function [unew,vnew]=LinIntPol2D(xmin,xmax,ymin,ymax,u,v)
P(1,1)=xmin;
P(1,2)=ymin;
P(2,1)=xmax;
P(2,2)=ymin;
P(3,1)=xmax;
P(3,2)=ymax;
P(4,1)=xmin;
P(4,2)=ymax;

unew=0.25*(P(1,1)*(1-v)*(1-u)+P(2,1)*(1-v)*(1+u)+...
  P(3,1)*(1+v)*(1+u)+P(4,1)*(1+v)*(1-u));

vnew=0.25*(P(1,2)*(1-v)*(1-u)+P(2,2)*(1-v)*(1+u)+...
  P(3,2)*(1+v)*(1+u)+P(4,2)*(1+v)*(1-u));

end


function  calcLsLines(P)
global n_Ls
global n_Ls_inv
global projectedFace
%% 1st direction
v1=P(1,end,:)-P(1,1,:);
v2=P(end,end,:)-P(end,1,:);
n_Ls=v1+v2;
n_Ls=(1/sqrt([n_Ls(1) n_Ls(2)]*[n_Ls(1);n_Ls(2)]))*n_Ls;
%% 2nd direction
u1=P(end,1,:)-P(1,1,:);
u2=P(end,end,:)-P(1,end,:);
n_Ls_inv=u1+u2;
n_Ls_inv=(1/sqrt([n_Ls_inv(1) n_Ls_inv(2)]*[n_Ls_inv(1);n_Ls_inv(2)]))*n_Ls_inv;
%% plot result
plot([0 n_Ls(1)*0.25], [0 n_Ls(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
plot([0 n_Ls_inv(1)*0.25], [0 n_Ls_inv(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
end