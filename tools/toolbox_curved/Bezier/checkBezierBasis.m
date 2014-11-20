function plot_cases
clc; clear all; close all;
global Bezier2d;global n_Ls;global n_Ls_inv;global projectedFace;global Face;global BezierSurface;global GeometricRepresentation;
NGeo=2;
onlyPlotData=1;
setWindows=1;
load faces.dat
%load /home/stephen/PMLalgorithm_cases/Testing_Bezier_domain010_nElems002_order04_PML00_zeta0E+00_polynom_CFL0.5_N10_Parts500_DoPML_False/BezierControlPoints.dat
load /scratch/iagortwe/Testcases/semicircle/BezierControlPoints.dat
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
      P(p,q,:,i)=faces((q-1)*(N+1)+p+index1,:);
    end
  end
end
%%
GeometricRepresentation=5;
MaxIter=10;
Tolerance=1E-5;
alpha0=0.5;
%pvec=[0.4 0.2 0.001];v=[0 1 0];Face=9;
pvec=[0.6 0.2 0.001];v=[1 -1 1];Face=9; %miss
%pvec=[0.6 0.05 0.2];v=[-1 1 -1];Face=9;pvec=pvec+v*0.1 %multiple
%intersection
%pvec=[0.6 0.2 0.001];v=[1 -1 0];Face=9; %miss
%pvec=[0.2 0.1 0.001];v=[0.7 0.9 0.3];Face=9;
%pvec=[0.4 0.2 0.001];v=[0 0 1];Face=10;
%pvec=[0.3 0.3 0.001];v=[1 1 0];Face=9; % center of face
% pvec=[sqrt(9/2)/10 sqrt(9/2)/10 0.001];v=[0 0 1];Face=10;  % center of upper face
pvec=[-0.3519 0.1426 0]; v=[0.87 -0.5 0];Face=5;
pvec=[0.3 0.3 0.001];v=[0.2 -1 0.2];Face=8;
% pvec=[0.25 0.25 0.001];v=[0 1 0];Face=6;
pvec=[0.25 0.25 0.000];v=[0 0 1];Face=9;
% pvec=[ 0.25 0.25 0.000 ];v=[-1 0 0];Face=9; %9 hier ist 4 im code
%%

BezierSurface=figure;  hold on;set(gcf, 'color', 'white');view(29,60);grid on;xlabel('x');ylabel('y');zlabel('z');
%plot3(faces(:,1),faces(:,2),faces(:,3),'ro','MarkerSize',10,'LineWidth',5)
Xi_NGeo=linspace(-1,1,N+1);
delta=2/N
%
% Xi_NGeo=zeros(1,N+1);
% for I=1:N+1
%   Xi_NGeo(I) = 2/N*(I-1)-1;
% end
count=0;
for i=1:nSides
  disp([num2str(i)])
  if i==Face
    %if i==7||i==9||i==3||i==2
    index1=(i-1)*(N+1)^2+1;
    index2=i*(N+1)^2;
    %plot3(faces(index1:index2,1),faces(index1:index2,2),faces(index1:index2,3),'ro','MarkerSize',10,'LineWidth',5)
    plotBezierSurface3D(N,P(:,:,:,i),i);
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

v=1/(v*v')*v; %direction

p2=pvec+alpha0*v;

plot3(pvec(1),pvec(2),pvec(3),'or','LineWidth',4,'MarkerSize',10)
plot3(p2(1),p2(2),p2(3),'or','LineWidth',4,'MarkerSize',10)
quiver3(pvec(1),pvec(2),pvec(3),v(1),v(2),v(3),0,'LineWidth',2,'MarkerSize',10)
view(35,26)
%view(v) % view in direction of the vector
if true(onlyPlotData),return;end;
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
monitor='0';set(0,'Units','pixels');scnsize = get(0,'ScreenSize');position = get(gcf,'Position');outerpos = get(gcf,'OuterPosition');borders = outerpos - position;edge = -borders(1)/2;
if strcmp(monitor,'0'),pos1=[0,1,scnsize(3)/10,scnsize(4)/2];
                       pos2=[375,1,scnsize(3)/6,scnsize(4)/2];
                       pos3=[1000,200,scnsize(3)/10,scnsize(4)/3];
                       pos4=[1000,1000,scnsize(3)/10,scnsize(4)/3];end;
if strcmp(monitor,'1'),pos1=[outerpos(1)*5.5,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos2=[outerpos(1)*6.7,scnsize(4),scnsize(3)/3-edge,scnsize(4)/1.25];
  pos3=[outerpos(1)*5.5,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];
  pos4=[outerpos(1)*6.7,-outerpos(1)/1.5,scnsize(3)/3-edge,scnsize(4)/2];end;
if strcmp(monitor,'2'),pos1=[1400,550,scnsize(3)/9,scnsize(4)/2.2];
  pos2=[1770,550,scnsize(3)/6,scnsize(4)/2.2];
  pos3=[1400,50,scnsize(3)/9,scnsize(4)/2.2];
  pos4=[1770,50,scnsize(3)/9,scnsize(4)/2.2];end;
if true(setWindows),set(gcf,'OuterPosition',pos1);end;
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
projectedFace = figure;  hold on;set(gcf, 'color', 'white');if true(setWindows),set(gcf,'OuterPosition',pos2);end;
%plotBezierSurface3D(N,P(:,:,:,i));
BezierControlPoints2Dorg=BezierControlPoints2D(:,:,:,Face);
plotProjectedBezierSurface2D(N,BezierControlPoints2D(:,:,:,Face));
%%
BezierControlPoints1D=zeros(dimensions(1),dimensions(2),dimensions(4));
xi=linspace(-1,1,NGeo+1);eta=linspace(-1,1,NGeo+1);
minmax=zeros(2,NGeo+1)-2;
%disp(['n_Ls_inv*n_Ls=' num2str(n_Ls_inv*n_Ls') ])

boundaries=[-1 1; -1 1];
DoClipping=ones(2,1);
sarray=zeros(MaxIter,2,2);
sarray(:,:,1)=-1.0;
sarray(:,:,2)=1.0;
iiter=0;
jiter=0;

for K=1:MaxIter
  for J=1:2
    if true(DoClipping(J))
      %% get Line L_s (vector n_Ls and n_Ls_inv)
      calcLsLines(BezierControlPoints2D(:,:,:,Face)) % get n_Ls and n_Ls_inv vectors, attention: n_Ls_inv is used for Ls and vice versa
      switch J
        case 1
          iiter=iiter+1;
          vector = [n_Ls_inv(1) n_Ls_inv(2)];posX=pos3;
        case 2
          jiter=jiter+1;
          vector = [n_Ls(1) n_Ls(2)];posX=pos4;
      end
      f=figure;hold on;set(gcf, 'color', 'white');grid on;xlabel('x');ylabel('y');plot([-1 1],[0 0],'k-');if true(setWindows),set(gcf,'OuterPosition',posX);end;
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
      % ! check upper and lower intersection with convex hull
      %%  1.) check traverse line UPPER/LOWER
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
      %% 2.) check BEGINNING/END upper convex hull
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
      
      %% 3.) check vertical line LEFT/RIGHT of convex hull
      if  minmax(1,1)*minmax(2,1)<=0 % 1 spalte anfang/ende
        smin_proposed=xi(1);
        smin=min(smin_proposed,smin);
      end
      if  minmax(1,end)*minmax(2,end)<=0 % 2 spalte anfang/ende
        smax_proposed=xi(NGeo+1);
        smax=max(smax_proposed,smax);
      end
      %%
      if or(smin==1.5,smax==-1.5)
        disp('Miss!')
        return
      end
      sarray(K,J,1)=smin;
      sarray(K,J,2)=smax;
      % disp([' k= ',num2str(K),' J= ',num2str(J),' smin= ',num2str(smin),' smax= ',num2str(smax)])
      plot([smin smax],[0 0],'r+')
      boundaries(J,:)=[smin smax];
      %% CLIP
      switch J
        case 1
          direction='u';
        case 2
          direction='v';
      end
      %if (smax-smin)/2*100>40,disp('gebiet muss geteilt werden!');end
      BezierControlPoints2D(:,:,:,Face)=BezierClipping(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face)); % J={u,v}
      %calcLsLines(BezierControlPoints2D(:,:,:,Face)) %moved to top
      x_sol=mean(mean(BezierControlPoints2D(:,:,1,Face)));
      y_sol=mean(mean(BezierControlPoints2D(:,:,2,Face)));
      disp(['K=' num2str(K) ' J=' num2str(J) ' Clip ' direction '  smin=' sprintf('%1.2e',smin) ', smax=' sprintf('%1.2e',smax) ' smax-smin=' sprintf('%1.2e',smax-smin) ' (' sprintf('%3.2f',(smax-smin)/2*100) '%) delta=' sprintf('%1.2e',sqrt(x_sol^2+y_sol^2))])
      %       disp(['sarray=' num2str(sarray(K,J,1)),'   ',num2str(sarray(K,J,2))])
      %% check epsilon
      if abs(smax-smin)<1e-4
        DoClipping(J)=0;
      end
      str=input('Press ENTER to continue or type "y" to exit: ','s');
      if strcmp(str,'y'),error('You have terminated the program!');end;
      if sqrt(x_sol^2+y_sol^2)<Tolerance
        disp(['sqrt(x_sol^2+y_sol^2)<' num2str(Tolerance) ': origin found!'])
        plot(x_sol,y_sol, 'g+','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        plot(x_sol,y_sol, 'go','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        % compute the correct smin and smax value
        umean=0.5*sarray(K,1,1)+0.5*sarray(K,1,2);
        vmean=0.5*sarray(K,2,1)+0.5*sarray(K,2,2);
        %umean=0.;
        %vmean=0.;
        disp([' umean: ',num2str(umean),' vmean: ',num2str(vmean)])
        for i=K-1:-1:1
          [umean,vmean]=LinIntPol2D(sarray(i,1,1),sarray(i,1,2),...
            sarray(i,2,1),sarray(i,2,2),...
            umean,vmean);
        end
        %         [umean,vmean]=LinIntPol2D(-1.,1.,...
        %           -1.,1.,...
        %           umean,vmean);
        disp([' umean: ',num2str(umean),' vmean: ',num2str(vmean)])
        %         % assuming that the center of patch is intersection with ray
        %         %smean=zeros(2,1);
        %         %smean(1)=mean(sarray(iiter,1,:));
        %         s=1;
        %         for K=iiter:-1:1
        %           s=1+sarray(K,1,1)+(sarray(K,1,2)-sarray(K,1,1))*0.5*(s);
        %         end
        %         smean(1)=-1+s;
        %         s=1;
        %         for K=jiter:-1:1
        %           s=1+sarray(K,2,1)+(sarray(K,2,2)-sarray(K,2,1))*0.5*(s);
        %         end
        %         smean(2)=-1+s;
        %
        % %  smean(L)=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smean(L));
        s=mean(sarray(iiter,1,:));
        disp(['                               U=' num2str(s)  ])
        for K=iiter-1:-1:1
          s=sarray(K,1,1)+(sarray(K,1,2)-sarray(K,1,1))*0.5*(s+1);
          disp(['Umin=' num2str(sarray(K,1,1)) '   Umax=' num2str(sarray(K,1,2)) '   U=' num2str(s)  ])
        end
        smean(1)=s;
        s=mean(sarray(jiter,2,:));
        disp(['                               V=' num2str(s)  ])
        for K=jiter-1:-1:1
          s=sarray(K,2,1)+(sarray(K,2,2)-sarray(K,2,1))*0.5*(s+1);
          disp(['Vmin=' num2str(sarray(K,1,1)) '   Vmax=' num2str(sarray(K,1,2)) '   V=' num2str(s)  ])
        end
        smean(2)=s;
        disp(['U/Vmean ', num2str(smean(1)),' ',num2str(smean(2))])
        %smean(1:2)=0
        %%% old and wrong point
        %         smean=zeros(2,1);
        %         for L=1:2
        %           switch L
        %             case 1
        %               niter=iiter;
        %               iRow1=1;
        %               iRow2=2;
        %             case 2
        %               niter=jiter;
        %               iRow1=2;
        %               iRow2=1;
        %           end
        %           smin=sarray(niter,L,1);
        %           smax=sarray(niter,L,2);
        %           smean(L)=0.5*smin+0.5*smax;
        %           %smean(L)=0.; % center of deepest patch
        %           for i=niter-1:-1:1
        %             %disp(['i ',num2str(i)]);
        %             smin=LinIntPol1D(sarray(i,L,iRow1),sarray(i,L,iRow2),smin);
        %             smax=LinIntPol1D(sarray(i,L,iRow2),sarray(i,L,iRow1),smax);
        %             smean(L)=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smean(L));
        %           end
        %           %smean(L)=LinIntPol1D(-1,1,smean(L));
        %           %smean(L)=0.5*smin+0.5*smax;
        %           %          smean(J)=0.5*smin+0.5*smax;
        %           %smean(J)=0.5*(smin+smax);
        %           disp(['L: ',num2str(L),' smin: ',num2str(smin),' smax: ',num2str(smax),' mean: ',num2str(smean(L))])
        %         end
        
        %plotBezierPoint2D(N,BezierControlPoints2Dorg,[0;0])
        plotBezierPoint2D(N,BezierControlPoints2Dorg,[smean(1);smean(2)],'.')
        plotBezierPoint2D(N,BezierControlPoints2Dorg,[umean;vmean],'x')
        %plotBezierPoint2D(N,BezierControlPoints2Dorg,[vmean;umean])
        %plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[0.;0.3])
        
        a=plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[umean;vmean]);
        
        a=plotBezierPoint3D(N,BezierControlPoints(:,:,:,Face),[smean(1);smean(2)]);
        a=[a(1) a(2) a(3)]
        v0=a-pvec;
        %alpha=sqrt(v0*v0');
        alpha=v0*v';
        disp(['alpha0 = ' num2str(alpha0) '    alpha = ' num2str(alpha)])
        if alpha<0
          disp('Particle oposite direction !')
        else
          if alpha <= alpha0
            disp('Particle Hit !')
          else
            disp('Particle MISS !')
          end
        end
        %% calc distance
        return
      end % sqrt(x_sol^2+y_sol^2)<Tolerance
    end % true(DoClipping(J))
  end % J=1:2
  disp(' ')
end % K=1:MaxIter
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










function plotBezierSurface3D(N,P,SideID) % 3D plot using Bezier control points and bernstein polynomials
global ParticleBezierSurface;global Face;global GeometricRepresentation;
xi=linspace(-1,1,GeometricRepresentation);eta=linspace(-1,1,GeometricRepresentation);
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

function plotProjectedBezierSurface2D(N,P) % 2D plot using Bezier control points and bernstein polynomials
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
axis equal
end



function result=plotBezierPoint3D(N,P,boundaries)
global BezierSurface

% boundaries(1)=0
% boundaries(2)=0
method='Bernstein';
method='DeCasteljau';
switch method
  case 'Bernstein'
    Bezier3d=0;
    for q=1:N+1
      for p=1:N+1
        Bezier3d=Bezier3d+B(N,p-1,boundaries(1))*B(N,q-1,boundaries(2))*P(p,q,:);
      end
    end
  case 'DeCasteljau'
    u=0.5*(boundaries(1)+1);
    v=0.5*(boundaries(2)+1);
    N2=N;
    BezierControlPoints3D_temp=P;
    for I=1:N %de casteljau schritti
      for q=1:N2
        for p=1:N2
          A=[BezierControlPoints3D_temp(p  ,q,:) BezierControlPoints3D_temp(p  ,q+1,:);...
             BezierControlPoints3D_temp(p+1,q,:) BezierControlPoints3D_temp(p+1,q+1,:)];
          for K=1:3
            BezierControlPoints3D_temp(p,q,K)=[1-u u]*A(:,:,K)*[1-v;v];
          end
        end
      end
      N2=N2-1;
    end
end
plot3(BezierControlPoints3D_temp(1,1,1),BezierControlPoints3D_temp(1,1,2),BezierControlPoints3D_temp(1,1,3),'ko','LineWidth',5,'MarkerSize',10,'Parent',get(BezierSurface, 'children'));
plot3(BezierControlPoints3D_temp(1,1,1),BezierControlPoints3D_temp(1,1,2),BezierControlPoints3D_temp(1,1,3),'y.','LineWidth',5,'MarkerSize',10,'Parent',get(BezierSurface, 'children'));
result=BezierControlPoints3D_temp(1,1,:);
end

function plotBezierPoint2D(N,P,boundaries,symbol)
global projectedFace
Bezier2d=0;
%boundaries(1)
%boundaries(2)
for q=1:N+1
  for p=1:N+1
    Bezier2d=Bezier2d+B(N,p-1,boundaries(1))*B(N,q-1,boundaries(2))*P(p,q,:);
  end
end
plot(Bezier2d(1),Bezier2d(2),'bo','MarkerSize',20,'Parent',get(projectedFace, 'children'));
plot(Bezier2d(1),Bezier2d(2),['b' symbol],'MarkerSize',20,'Parent',get(projectedFace, 'children'));
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
        Bezier2d(I,J,:)=Bezier2d(I,J,:)+B(N,p-1,xi(I))*B(N,q-1,eta(J))*P(p,q,:);
      end
    end
  end
end
plot(Bezier2d(:,:,1),Bezier2d(:,:,2),'r-',Bezier2d(:,:,1)',Bezier2d(:,:,2)','r-','Parent',get(projectedFace, 'children'));
end

function B=B(n,k,x)
%disp(['n=' num2str(n) '   k=' num2str(k)])
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

function BezierClipping=BezierClipping(direction,smin,smax,N,BezierControlPoints2D)
global projectedFace;
switch direction
  % ===================================================================================================
  case 1 % u direction
    color='k.';
    Bezier2d_temp=zeros(N+1,N+1,2);;
    % 1. top xi
    s=smax;
    for q=1:N+1
      for p=1:N+1
        for l=1:p
          Bezier2d_temp        (p,q,:)=...
          Bezier2d_temp        (p,q,:)+...
          BezierControlPoints2D(l,q,:)*B(p-1,l-1,s);
        end
      end
    end
    % 1. bottom xi
    BezierControlPoints2D=Bezier2d_temp;
    Bezier2d_temp=zeros(N+1,N+1,2);
    s=(smin+1)/(smax+1); %[-1, +1]
    s=2*(1-s)-1;
    for q=1:N+1
      for p=1:N+1
        for l=1:p
          Bezier2d_temp        (N+2-p,q,:)=...
          Bezier2d_temp        (N+2-p,q,:)+...
          BezierControlPoints2D(N+2-l,q,:)*B(p-1,l-1,s);
        end
      end
    end
    plot(Bezier2d_temp(:,:,1),Bezier2d_temp(:,:,2),color,Bezier2d_temp(:,:,1)',Bezier2d_temp(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    for J=1:N+1
      text(Bezier2d_temp(:,J,1),Bezier2d_temp(:,J,2),'u','FontSize',22,'HorizontalAlignment','center','VerticalAlignment','bottom','Parent',get(projectedFace, 'children'))
      text(Bezier2d_temp(:,J,1),Bezier2d_temp(:,J,2),[num2str([Bezier2d_temp(:,J,1) Bezier2d_temp(:,J,2)],'(%1.2e,%1.2e)')],'HorizontalAlignment','center','VerticalAlignment','top','Parent',get(projectedFace, 'children'))
    end
    BezierClipping=Bezier2d_temp;
    % ===================================================================================================
  case 2 % v direction
    color='r.';
    Bezier2d_temp=zeros(N+1,N+1,2);
    % 2. top eta
    s=smax;
    for q=1:N+1
      for p=1:N+1
        for l=1:p
          Bezier2d_temp        (q,p,:)=...
          Bezier2d_temp        (q,p,:)+...
          BezierControlPoints2D(q,l,:)*B(p-1,l-1,s);
        end
      end
    end
    BezierControlPoints2D=Bezier2d_temp;
    % 2. bottom eta
    Bezier2d_temp=zeros(N+1,N+1,2);
    s=(smin+1)/(smax+1); %[-1, +1]
    s=2*(1-s)-1;
    for q=1:N+1
      for p=1:N+1 %N+1:-1:1
        for l=1:N+2-p
          Bezier2d_temp        (q,p,:)=...
          Bezier2d_temp        (q,p,:)+...
          BezierControlPoints2D(q,N+2-l,:)*B(N+1-p,l-1,s);
        end
      end
    end
    plot(Bezier2d_temp(:,:,1),Bezier2d_temp(:,:,2),color,Bezier2d_temp(:,:,1)',Bezier2d_temp(:,:,2)',color,'Parent',get(projectedFace, 'children'));
    
    
    for J=1:N+1
      text(Bezier2d_temp(:,J,1),Bezier2d_temp(:,J,2),'v','HorizontalAlignment','center','VerticalAlignment','bottom','Parent',get(projectedFace, 'children'))
      text(Bezier2d_temp(:,J,1),Bezier2d_temp(:,J,2),[num2str([Bezier2d_temp(:,J,1) Bezier2d_temp(:,J,2)],'(%1.2e,%1.2e)')],'HorizontalAlignment','center','VerticalAlignment','top','Parent',get(projectedFace, 'children'))
    end
%     sprintf('%1.2e'
%     strValues = strtrim(cellstr(num2str([X(:) Y(:)],'(%d,%d)')));
%     text(X,Y,strValues,'VerticalAlignment','bottom');
    BezierClipping=Bezier2d_temp;
    % ===================================================================================================
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
v1=P(1  ,end,:)-P(1  ,1,:);
v2=P(end,end,:)-P(end,1,:);
n_Ls=v1+v2;
n_Ls=(1/sqrt([n_Ls(1) n_Ls(2)]*[n_Ls(1);n_Ls(2)]))*n_Ls;
%% 2nd direction
u1=P(end,  1,:)-P(1,  1,:);
u2=P(end,end,:)-P(1,end,:);
n_Ls_inv=u1+u2;
n_Ls_inv=(1/sqrt([n_Ls_inv(1) n_Ls_inv(2)]*[n_Ls_inv(1);n_Ls_inv(2)]))*n_Ls_inv;
%% plot result
plot([0 n_Ls(1)*0.25], [0 n_Ls(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
plot([0 n_Ls_inv(1)*0.25], [0 n_Ls_inv(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
%n_Ls
%n_Ls_inv
end