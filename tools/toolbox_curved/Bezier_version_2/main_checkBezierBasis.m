% ==================================================================
% Description
% ==================================================================
% CAUTION: this algorithm only converges for certain problems, where
% a surface splitting is not needed (this would require the algortihm to be recursive)
addpath(['/home/stephen/Documents/MATLAB/export_fig/']);
clc; clear all; close all;
% global Bezier2d;global n_Ls;global n_Ls_inv;
global projectedFace;global Face;global BezierSurfaceFigure
global Volume_box
NGeo=2;
load faces.dat
%load /home/stephen/PMLalgorithm_cases/Testing_Bezier_domain010_nElems002_order04_PML00_zeta0E+00_polynom_CFL0.5_N10_Parts500_DoPML_False/BezierControlPoints.dat
load BezierControlPoints.dat
%factorial(5)
faces                    = BezierControlPoints; % externally specified faces order in vector form (z-curve type)
dimension                = size(faces);         % number of external faces and coordinates (x,y and z)
ZeroRadiusTolerance      = 1E-9;                % clipping tolerance, i.e., the minimum distance from intersection point to middle of clipped bézier polynomial
DiffSminSmaxTolerance    = 1E-9;                % set DoClip=0 for xi/eta direction if this occurs
MaxClipIterations        = 6;                   % number of iterations for clipping
N_plot                   = 5;                  % plotting number of points for 3D surface plot of the bezier surface

createPlots              = 0;                   % export plots to .pdf files
create1Dplots            = 0;                   % make plots for 1D projected cliped bezier

VerifyControlPoints      = 0;                   % interpolate bezier basis and use lagrange onde (CL) 
                                                % to change basis back to the bezier basis
Points2DSequence         = 0;                   % plot Pi,j: press "ENTER" for single output or "y" for complete output

ElevateControlPoints     = 1;                   % Perform degree elevation in 3D
elevation                = 1;                   % elevate from p->p+elevation
elevation_direction      = 1;                   % 1:xiber
                                                % 2:eta
                                                
OrientedSlabBox          = 0;                   % use bezier surface and create a bounding box (oriented-slab box)

OrientedSlabBoxCell      = 0;                   % use bezier surfaces from the cell faces and create a bounding box (oriented-slab box)
cell                     = [6,7,8,9,10,11];     % face IDs for the cell -> for cell bounding box
cell=9

adjustScreen             = 0;                   % adjust plot window sizes



% points
N=dimension(2)-1;
N=NGeo;
nSides=dimension(1)/(N+1)^2;
Volume_box=zeros(1,nSides);
if(ZeroRadiusTolerance>1e-5),disp('your ZeroRadiusTolerance might be set too high! Continuing'); end;
%% assign control points from loaded values in "faces" variable
% ("load BezierControlPoints.dat" and "faces=BezierControlPoints")
% re-shape values from vector to a 4 dimensional tensor "P(p,q,:,:)"
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
%% specify the particle origin "BasePointIC" and direction "v" (the direction will be normed v/|v|->v)


%BasePointIC=[0.6 0.05 0.2];v=[-1 1 -1];Face=9;BasePointIC=BasePointIC+v*0.1 %double
%BasePointIC=[0.4 0.2 0.001];v=[0 0 1];Face=10;
%BasePointIC=[0.3 0.3 0.001];v=[0 0 1];Face=10;
%BasePointIC=[0.3 0.3 0.001];v=[0.2 -1 0.2];Face=8;
 


% hit --------------------------------------------------------------------
  BasePointIC=[0.4 0.2 0.001];v=[0 1 0];Face=9;
  BasePointIC=[0.4 0.2 0.001];v=[2 1 0.4];Face=9;
 
% miss --------------------------------------------------------------------
%BasePointIC=[0.6 0.2 0.001];v=[1 -1 0];Face=9; % misses correctly
%BasePointIC=[0.6 0.2 0.001];v=[1 -1 1];Face=9; % misses correctly: long convergence 

% should be hit, but does not converge (no splitting implemented) ---------
%BasePointIC=[0.4 0.2 0.001];v=[0.7 0.9 1];Face=9; 
%BasePointIC=[0.3 0.3 0.001];v=[1 1 0];Face=9;
%BasePointIC=[0.6 0.05 0.2];v=[-1 1 -1];Face=9;BasePointIC=BasePointIC+v*0.1 %double intersection
% BasePointIC=[0.8 0.5 0.4];v=[-0.9 -0.4 -0.5];Face=9; % does not converge: (0.87%)

% new, wrong point found?! ------------------------------------------------
% -> not resolved: resulting intersection point does not fit
%

% for picture creation
BasePointIC=[0.8 0.55 0.4];v=[-0.9 -0.4 -0.5];Face=9; % does not converge: (0.87%)
velocity=0.3;
SetView=[-145.0375,38.9157];

% for PICASSO checking
%BasePointIC=[-0.1, 0.15, -0.1];v=[0.9, 0.4, 0.5];Face=9; % does not converge: (0.87%)
velocity=1;

%%
BezierSurfaceFigure=figure;  hold on;set(gcf, 'color', 'white');view(29,60);grid on;xlabel('x');ylabel('y');zlabel('z');
%plot3(faces(:,1),faces(:,2),faces(:,3),'ro','MarkerSize',10,'LineWidth',5)
Xi_NGeo=zeros(1,N+1);
for I=1:N+1
  Xi_NGeo(I) = 2/N*(I-1)-1;
end
count=0;
for i=1:nSides
  %if i>0%==9 % only plot these faces
    %if i==7||i==9||i==3||i==2 % only plot these faces
  %if i==
  switch i
    case num2cell(cell) % complete cell: all faces in "cell"
      index1=(i-1)*(N+1)^2+1; % shifting in data file where all sides are located above each other
      index2=i*(N+1)^2;
      %plot3(faces(index1:index2,1),faces(index1:index2,2),faces(index1:index2,3),'ro','MarkerSize',10,'LineWidth',5)
      PlotSuperSampledBezier3D(N,P(:,:,:,i),i,N_plot,3);
      count=count+1;
      index=i;
      fprintf(2,'Face: %3.0f \n',i);
      %fprintf(2,'Face: %3.0f \n',i);
    otherwise
      %else
      fprintf(1,'Face: %3.0f \n',i);
  end
end

% what is this ?
if count==111
  x=P(:,:,1,index);y=P(:,:,2,index);z=P(:,:,3,index);
  %   b=max([max(x),max(y),max(z)])
  %   a=min([min(x),min(y),min(z)])
  percent=0.3;
  xlim([min(min(x))-(1-percent)*abs(min(min(x))) max(max(x))+(1+percent)*abs(max(max(x)))]);
  ylim([min(min(y))-(1-percent)*abs(min(min(y))) max(max(y))+(1+percent)*abs(max(max(y)))]);
  zlim([min(min(z))-(1-percent)*abs(min(min(z))) max(max(z))+(1+percent)*abs(max(max(z)))]);
  
end
axis equal

%%
% use control points to interpolate to a lagrange basis and use it the
% inverse of the Vandermonde to change the basis back to the bezier control
% points
if true(VerifyControlPoints)
  checkControlPoints(N,P(:,:,:,Face),N_plot)
end




%% particle plot

%Face=9;
dimensions=size(P);
BezierControlPoints=P;
BezierControlPoints(:,:,:,Face);

v=1/(v*v')*v;
VeloVecIC=velocity*v;
p2=BasePointIC+VeloVecIC;
plot3(BasePointIC(1),BasePointIC(2),BasePointIC(3),'or','LineWidth',4,'MarkerSize',10);
plot3(p2(1),p2(2),p2(3),'or','LineWidth',4,'MarkerSize',10);
quiver3(BasePointIC(1),BasePointIC(2),BasePointIC(3),v(1),v(2),v(3),0,'LineWidth',2,'MarkerSize',10);
view(35,26)
view(-v) % view in direction of the vector
axis tight
pause(0.1)


%% Bounding Box Face: Oriented Slab Box
if true(OrientedSlabBox)
  CalcOrientedSlabBox(BasePointIC,VeloVecIC,BezierControlPoints(:,:,:,Face),NGeo,0)
end

%% Bounding Box Cell: Oriented Slab Box
if true(OrientedSlabBoxCell)
  CalcOrientedSlabBoxCell(BasePointIC,VeloVecIC,BezierControlPoints(:,:,:,cell),NGeo)
end


%% Elevation
if true(ElevateControlPoints)
  ElevationMatrix=CalcElevationMatrix(N,elevation);
  BezierSurfaceFigure=figure;  hold on;set(gcf, 'color', 'white');view(29,60);grid on;xlabel('x');ylabel('y');zlabel('z');
  PlotSuperSampledBezier3D(N,P(:,:,:,Face),Face,N_plot,3);
  
  plot3(BasePointIC(1),BasePointIC(2),BasePointIC(3),'or','LineWidth',4,'MarkerSize',10);
  plot3(p2(1),p2(2),p2(3),'or','LineWidth',4,'MarkerSize',10);
  quiver3(BasePointIC(1),BasePointIC(2),BasePointIC(3),v(1),v(2),v(3),0,'LineWidth',2,'MarkerSize',10);
  
  
  P_elevated=ElevateBezierSurface(N,Face,P,elevation,elevation_direction,ElevationMatrix);
  if true(OrientedSlabBox)
    CalcOrientedSlabBox(BasePointIC,VeloVecIC,P_elevated,NGeo+elevation,Face)
  end
   
end

 
%% not working properly, layers are not arranged correctly
if true(createPlots)
  view(SetView);
  %saveas(gcf, '/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/3D_plot_in_vector_direction_with_Control_points_saveas', 'pdf');
  %print(BezierSurfaceFigure,'/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/3D_plot_in_vector_direction_with_Control_points_print','-depsc')
  export_fig(BezierSurfaceFigure,['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/3D_plot_with_Control_points_export_fig.pdf']);
end
%% projection 3D->2D
%1. create perpendicular vectors
if abs(v(3))<1E-6
  n1=[-v(2)-v(3),v(1),v(1)];
else
  n1=[v(3),v(3),-v(1)-v(2)];
end
n1=1/sqrt(n1*n1')*n1;
quiver3(BasePointIC(1),BasePointIC(2),BasePointIC(3),n1(1),n1(2),n1(3),0,'LineWidth',2,'MarkerSize',10)
n2=cross(v,n1);
quiver3(BasePointIC(1),BasePointIC(2),BasePointIC(3),n2(1),n2(2),n2(3),0,'LineWidth',2,'MarkerSize',10)
%dimensions=size(BezierControlPoints);
d1ij=zeros(dimensions(1),dimensions(2),dimensions(4));
d2ij=zeros(dimensions(1),dimensions(2),dimensions(4));
for SideID=1:nSides
  for q=1:NGeo+1
    for p=1:NGeo+1
      d1ij(p,q,SideID)=([BezierControlPoints(p,q,1,SideID) BezierControlPoints(p,q,2,SideID) BezierControlPoints(p,q,3,SideID)]-BasePointIC)*n1';
      d2ij(p,q,SideID)=([BezierControlPoints(p,q,1,SideID) BezierControlPoints(p,q,2,SideID) BezierControlPoints(p,q,3,SideID)]-BasePointIC)*n2';
    end
  end
end



% plotSingleBezierPoint(N,P(:,:,:,9),[1,-1],3)
% return
%% not working properly, layers are not arranged correctly
if true(createPlots)
  view(-v) % view in direction of the vector
  export_fig(BezierSurfaceFigure,['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/3D_plot_in_vector_direction_with_Control_points_export_fig.pdf']);
end
%%
if true(adjustScreen)
  AdjustScreenSize()
end
axis equal
%export_fig([pwd '/3D_plot.pdf']);
%% projection to 2D: PLOT 2D PROJECTED PLANE -> PlotSuperSampledBezier3D(...)
BezierControlPoints2D    = CalcBezierControlPoints2D(dimensions,nSides,NGeo,d1ij,d2ij,N,N_plot,createPlots);
BezierControlPoints2Dorg = BezierControlPoints2D(:,:,:,Face); % save the original points (BezierControlPoints2D will be updated within each clipping step)
if true(Points2DSequence)
  PlotBezierControlPoints2DSequence(BezierControlPoints2D(:,:,:,Face),NGeo)
end
% get Line L_s and plot (vector n_Ls and n_Ls_inv): PLOT THESE VECTORS IN 2D PROJECTED PLANE
[v1,v2]=calcLsLines(BezierControlPoints2D(:,:,:,Face));

xMin=min(min(BezierControlPoints2D(:,:,1,Face)));
xMin=xMin-abs(xMin*0.1);
xMax=max(max(BezierControlPoints2D(:,:,1,Face)));
xMax=xMax+abs(xMax*0.1);

yMin=min(min(BezierControlPoints2D(:,:,2,Face)));
yMin=yMin-abs(yMin*0.1);
yMax=max(max(BezierControlPoints2D(:,:,2,Face)));
yMax=yMax+abs(yMax*0.1);

ylim([min(xMin,yMin) max(xMax,yMax)])
xlim([min(xMin,yMin) max(xMax,yMax)])

plot(xlim,[0 0],'k-');plot([0 0],ylim,'k-');
plot(xlim,max(ylim)*[1 1],'k-');plot(max(xlim)*[1 1],ylim,'k-');


%%
BezierControlPoints1D=zeros(dimensions(1),dimensions(2),dimensions(4));
xi=linspace(-1,1,NGeo+1);eta=linspace(-1,1,NGeo+1);
minmax=zeros(2,NGeo+1)-2;
%disp(['n_Ls_inv*n_Ls=' num2str(n_Ls_inv*n_Ls') ])

boundaries           = [-1 1; -1 1]; % 
DoClipping           = ones(2,1);      % switch for both clipping directions, will be turned false is further clipping is not needed
sarray               = zeros(25,2,2);              % old
sarray(:,:,1)        = -1.0;                       % old
sarray(:,:,2)        =  1.0;                       % old
ClipConvergenceXiEta = zeros(MaxClipIterations,4);
XiArray              = [-1*ones(25,1) ones(25,1)]; % new
EtaArray             = [-1*ones(25,1) ones(25,1)]; % new

iiter=0;
jiter=0;
for K=1:MaxClipIterations
  for J=1:2
    if true(DoClipping(J))
      switch J
        case 1
          iiter=iiter+1;
          
          vector=v1;
          plot([0 v1(1)*0.25], [0 v1(2)*0.25],'r-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
          text(v1(1)*0.25,v1(2)*0.25,'v1','Parent',get(projectedFace, 'children'))
        case 2
          jiter=jiter+1;
          %vector = [n_Ls(1) n_Ls(2)]          %posX=pos4;
          vector=v2;
          plot([0 v2(1)*0.25], [0 v2(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
          text(v2(1)*0.25,v2(2)*0.25,'v2','Parent',get(projectedFace, 'children'))
      end
      
      %set(gcf,'OuterPosition',posX)
      for SideID=1:nSides
        for q=1:NGeo+1
          for p=1:NGeo+1
            BezierControlPoints1D(p,q,SideID)=([BezierControlPoints2D(p,q,1,SideID) BezierControlPoints2D(p,q,2,SideID)])*vector;
            if (SideID==Face)
              plot(BezierControlPoints2D(p,q,1,SideID)*[1 1], BezierControlPoints2D(p,q,2,SideID)*[1 1],'g+')
              BezierControlPoints1D(p,q,SideID)
            end
          end
        end
      end
      
      
      
      [smin,smax]   = CalcSminSmax(J,NGeo,xi,BezierControlPoints1D,K,createPlots,create1Dplots);   % [smin,smax]=CalcSminSmax()
      sarray(K,J,1) = smin;
      sarray(K,J,2) = smax;
      switch J
        case 1
          XiArray(K,1) = smin;
          XiArray(K,2) = smax;
        case 2
          EtaArray(K,1) = smin;
          EtaArray(K,2) = smax;
      end
      % disp([' k= ',num2str(K),' J= ',num2str(J),' smin= ',num2str(smin),' smax= ',num2str(smax)])


      boundaries(J,:)=[smin smax];
      %% CLIP
      %       if J==2
      %         smax=0.5
      %         smin=-0.5
      %       end
      
      % moved from here ************************

      
      if (smax-smin)/2*100>40,disp('gebiet muss geteilt werden!');end
      


      % clip and plot the new ploynomial (old version)
      %BezierControlPoints2D(:,:,:,Face) = BezierClipping(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2,createPlots); % J={u,v}: lalt
      % clip and plot the new ploynomial (new version)
      %BezierControlPoints2D(:,:,:,Face) = BezierClipping2(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2,createPlots);
      
      % plot clipped surface (old)
      %BezierClipping(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2,createPlots);
      
      % replot with new version (new)
      %CalcBezierControlPoints2D(dimensions,nSides,NGeo,d1ij,d2ij,N,N_plot,createPlots);
      BezierControlPoints2D(:,:,:,Face)=BezierClipping2(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2,createPlots);
      
      if true(Points2DSequence)
        PlotBezierControlPoints2DSequence(BezierControlPoints2D(:,:,:,Face),N)
      end
      
%        BezierControlPoints2D_temp = BezierClipping(J,smin,smax,N,BezierControlPoints2D(:,:,:,Face),2,createPlots); % J={u,v}
%       BezierControlPoints2D(:,:,1,Face)=BezierControlPoints2D_temp(:,:,1)';
%       BezierControlPoints2D(:,:,2,Face)=BezierControlPoints2D_temp(:,:,2)';
      % moved to here ************************
      x_sol=mean(mean(BezierControlPoints2D(:,:,1,Face)));
      %x_sol=sum(sum(BezierControlPoints2D(:,:,1,Face)))/((N+1)*(N+1));
      y_sol=mean(mean(BezierControlPoints2D(:,:,2,Face)));
      
      
      % get Line L_s and plot (vector n_Ls and n_Ls_inv)
      %calcLsLines(BezierControlPoints2D(:,:,:,Face))
      [v1,v2]=calcLsLines(BezierControlPoints2D(:,:,:,Face));
%       if mod(K,2)
%         vector=-vector
%       end
      

      
      delta=sqrt(x_sol^2+y_sol^2);
      switch J
        case 1
          direction='u';
          ClipConvergenceXiEta(K,1)=delta;
          ClipConvergenceXiEta(K,3)=(smax-smin)/2*100;
        case 2
          direction='v';
          ClipConvergenceXiEta(K,2)=delta;
          ClipConvergenceXiEta(K,4)=(smax-smin)/2*100;
      end
      disp(['K=' num2str(K) ' J=' num2str(J) ' Clip ' direction '  smin=' sprintf('%1.2e',smin) ', smax=' sprintf('%1.2e',smax) ' smax-smin=' sprintf('%1.2e',smax-smin) ' (' sprintf('%3.2f',(smax-smin)/2*100) '%) ZeroRadius=' sprintf('%1.2e',delta)])
      
      
      
      %% check epsilon - finish search if true !
      if abs(smax-smin)<DiffSminSmaxTolerance
        DoClipping(J)=0;
      end
      % % inquire for user input to continue each step
      %str=input('Press ENTER to continue or type "y" to exit: ','s');
      %if strcmp(str,'y'),error('You have terminated the program!');end;
      
      
      if sqrt(x_sol^2+y_sol^2)<ZeroRadiusTolerance
        disp(['sqrt(x_sol^2+y_sol^2)<' num2str(ZeroRadiusTolerance) ': origin found!'])
        plot(x_sol,y_sol, 'g+','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        plot(x_sol,y_sol, 'go','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
        disp('__________Xi_________     _______Eta________')
        disp([num2str([XiArray(1:max(iiter,jiter),:) EtaArray(1:max(iiter,jiter),:)])])
        
        
        % switch sign of even directions: i=2..4..6..8
       % XiArray(~mod(1:length(XiArray),2),:) =XiArray(~mod(1:length(XiArray),2),:)*(-1);
       % EtaArray(~mod(1:length(XiArray),2),:)=EtaArray(~mod(1:length(EtaArray),2),:)*(-1);
       % disp('__________Xi_________     _______Eta________')
       % disp([num2str([XiArray(1:max(iiter,jiter),:) EtaArray(1:max(iiter,jiter),:)])])
        
        
        %% compute the correct smin and smax value
        %% method 1
        umean=0.5*(sarray(K,1,1)+sarray(K,1,2)); % innerstes level
        vmean=0.5*(sarray(K,2,1)+sarray(K,2,2)); % innerstes level
        
        umean=mean(XiArray(K,:));
        vmean=mean(EtaArray(K,:));
        
        disp(' ')
        disp('=====================================================')
        disp('method 1')
        % alt
        % vorletzter eintrag wird verwendet
        for i=K-1:-1:1
          [umean,vmean]=LinIntPol2D(sarray(i,1,1),sarray(i,1,2),...
                                    sarray(i,2,1),sarray(i,2,2),...
                                    umean,vmean);
        end
        % neu
        % vorletzter eintrag wird verwendet
%         for i=K-1:-1:1
%           [umean,vmean]=LinIntPol2D(XiArray(i,1) , XiArray(i,2),...
%                                     EtaArray(i,1),EtaArray(i,2),...
%                                     umean,vmean);
%         end

        
        disp([' umean: ',sprintf('%3.5f',umean),' vmean: ',sprintf('%3.5f',vmean)])
        disp('=====================================================')
        %% compute the correct smin and smax value
        %% method 2
        name=[' xi';'eta'];
        smean=zeros(2,1);
        for L=1:2 % xi/eta direction
          switch L
            case 1
              disp(' ')
              disp('=====================================================')
              disp('method 2')
              niter=iiter;
            case 2
              niter=jiter;
          end
          smin=sarray(niter,L,1);
          smax=sarray(niter,L,2);
          switch L
            case 1
              smin=XiArray(niter,1);
              smax=XiArray(niter,2);
            case 2
              smin=EtaArray(niter,1);
              smax=EtaArray(niter,2);
          end
          
          
          smean(L)=0.5*(smin+smax);
          %if ~mod(niter,2) % -> iseven?
          %  smean(L)=-smean(L);
          %end
          for i=niter-1:-1:1
            % smin=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smin);
            % smax=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smax);
            % alt
           % smean(L)=LinIntPol1D(sarray(i,L,1),sarray(i,L,2),smean(L));
            % neu
            switch L
              case 1
                smean(L)=LinIntPol1D(XiArray(i,1),XiArray(i,2),smean(L));
                
              case 2
                smean(L)=LinIntPol1D(EtaArray(i,1),EtaArray(i,2),smean(L));
            end
            if ~mod(i,2) % -> iseven?
             % smean(L)=-smean(L);
            end
          end
          %smean(L)=0.5*smin+0.5*smax;
          %          smean(J)=0.5*smin+0.5*smax;
          %smean(J)=0.5*(smin+smax); 
          disp(['L=' num2str(L) ': (', name(L,:),') smin: ', sprintf('%3.4f',smin),' smax: ', sprintf('%3.4f',smax)])
        end
        disp([' umean: ', sprintf('%3.5f',smean(1)), ' vmean: ', sprintf('%3.5f',smean(2))])
        disp('=====================================================')
        
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
        
        % plot points on bezier curve 2D
        
        %plotSingleBezierPoint(N,BezierControlPoints2Dorg,[umean;vmean],2);       % method 1
        plotSingleBezierPoint(N,BezierControlPoints2Dorg,[smean(1);smean(2)],2); % method 2
        
        % plot points on bezier curve 3D
        %plotSingleBezierPoint(N,BezierControlPoints(:,:,:,Face),[smin;smax],3)
        plotSingleBezierPoint(N,BezierControlPoints(:,:,:,Face),[umean;vmean],3)
        
        
        return
        PlotConvergence(ClipConvergenceXiEta,MaxClipIterations,ZeroRadiusTolerance)

        return
      end
    end
  end
  disp(' ')
end

if MaxClipIterations>0
  plot(x_sol,y_sol, 'r+','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
  plot(x_sol,y_sol, 'ro','Parent',get(projectedFace, 'children'),'MarkerSize',20,'LineWidth',2)
  plot_2D_clipped_Bezier(N,BezierControlPoints2D(:,:,:,Face),[-1 1; -1 1])
  
  PlotConvergence(ClipConvergenceXiEta,MaxClipIterations,ZeroRadiusTolerance)
end





return
boundaries(:,1)=boundaries(:,1)-1E-4;
boundaries(:,2)=boundaries(:,2)+1E-4;

plot_2D_clipped_Bezier(N,BezierControlPoints2D(:,:,:,Face),boundaries);


%%
addpath(['/home/stephen/MATLAB_Programme/export_fig/']);
%export_fig([pwd '/2D_plot.pdf']);








