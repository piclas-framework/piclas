function CalcOrientedSlabBox(BasePointIC,VeloVecIC,P,N,boxVol)


global Volume_box_base
global Volume_box
global Face
monitor='large'; % 'small' or 'large'
global input
global x1;global x2;global x3;global x4;global x5;global x6;global x7;global x8
% global points
% global export
% global v1; global v2;  global v3; global nRoot;
global x_p1;global y_p1;global z_p1;global x_p2;global y_p2;global z_p2
input.removeAx=0;
epsilon=1E-5;
% x1=[0.3,0.3,0.3];x2=[2.3,0.3,0.3];x3=[2.3,2.3,0.3];x4=[0.3,2.3,0.3];
% x5=[0.3,0.3,2.3];x6=[2.3,0.3,2.3];x7=[2.3,2.3,2.3];x8=[0.3,2.3,2.3];
% x1=[0.,0.,0.];x2=[1.,0.,0.];x3=[1.,1.,0.];x4=[0.,1.,0.];
% x5=[0.,0.,1.];x6=[1.,0.,1.];x7=[1.,1.,1.];x8=[0.,1.,1.];
x_p1=0.25;y_p1=0.25;z_p1=0;
x_p2=0.25;y_p2=0.25;z_p2=1;





%% Side 2

% given bounding box
O     = reshape(P(1,1,:),3,1);  % origin of local coordinate system: 1st control point
n_xi  = reshape( (P(end,1,:)-P(1,1,:)) + (P(end,end,:)-P(1,end,:)) ,3,1);
n_eta = reshape( (P(1,end,:)-P(1,1,:)) + (P(end,end,:)-P(end,1,:)) ,3,1);

n1=n_xi/sqrt(n_xi'*n_xi);              % vector 1
n2=cross(n1,n_eta/sqrt(n_eta'*n_eta)); % vector 2
n3=cross(n2,n1);                       % vector 3


if(abs(n1'*n1-1)>1e-10), error('n1 is wrong');end;
if(abs(n2'*n2-1)>1e-10), error('n2 is wrong');end;
if(abs(n3'*n3-1)>1e-10), error('n3 is wrong');end;

if(n1'*n2), error('n1 not perpendicular to n2');end;
if(n1'*n3), error('n1 not perpendicular to n3');end;
if(n2'*n3), error('n2 not perpendicular to n3');end;

scale=0.1;
quiver3(O(1),O(2),O(3),n1(1)*scale,n1(2)*scale,n1(3)*scale,0,'LineWidth',2,'MarkerSize',10);
quiver3(O(1),O(2),O(3),n2(1)*scale,n2(2)*scale,n2(3)*scale,0,'LineWidth',2,'MarkerSize',10);
quiver3(O(1),O(2),O(3),n3(1)*scale,n3(2)*scale,n3(3)*scale,0,'LineWidth',2,'MarkerSize',10);

alpha_temp=zeros(N+1,N+1,3);
for k=1:3
  switch k
    case 1
      n_k=n1;
    case 2
      n_k=n2;
    case 3
      n_k=n3;
  end
  for q=1:N+1
    for p=1:N+1
      alpha_temp(p,q,k) = reshape(P(p,q,:)-P(1,1,:),3,1)'*n_k;
    end
  end
end


alpha_k=zeros(3,2);
for k=1:3
  alpha_k(k,:)=[min(min(alpha_temp(:,:,k))) max(max(alpha_temp(:,:,k)))];
end


beta1=[0 0.5757];   beta2=[0 0];        beta3=[-0.2828 0.4242];   % ranges for the 3 vectors
beta1=[alpha_k(1,:)];L1=beta1(2)-beta1(1);
beta2=[alpha_k(2,:)];L2=beta2(2)-beta2(1);
beta3=[alpha_k(3,:)];L3=beta3(2)-beta3(1);
switch boxVol
  case 0 %non-elevated box
    Volume_box_base=L1*L2*L3;
    disp(['Box volume = ' num2str(Volume_box_base)])
  otherwise
    Volume_box(Face)=L1*L2*L3;
    disp(['Box volume = ' num2str(Volume_box(Face)) ': Reduction elevated(' num2str(N) ')/NGeo(p=2) = ' num2str(Volume_box(Face)/Volume_box_base)])
end


%beta1=[0.1 1]; beta2=[0.1 1]; beta3=[0.1 1];

% given bounding box


%% Particle
%BasePointIC=[ 0.102668968989024       0.113049102130101      -0.247914810958384  ];
%VeloVecIC  =[-0.536702757124147       0.813970319052967      -0.222266664608417  ];

%Face=6;
E=BasePointIC'; % LastPartPos

x_p1=BasePointIC(1);    y_p1=BasePointIC(2);    z_p1=BasePointIC(3);
x_p2=x_p1+VeloVecIC(1); y_p2=y_p1+VeloVecIC(2); z_p2=z_p1+VeloVecIC(3);
%%
x1=O+beta1(1)*n1+beta2(1)*n2+beta3(1)*n3;
x2=O+beta1(2)*n1+beta2(1)*n2+beta3(1)*n3;
x3=O+beta1(2)*n1+beta2(2)*n2+beta3(1)*n3;
x4=O+beta1(1)*n1+beta2(2)*n2+beta3(1)*n3;
x5=O+beta1(1)*n1+beta2(1)*n2+beta3(2)*n3;
x6=O+beta1(2)*n1+beta2(1)*n2+beta3(2)*n3;
x7=O+beta1(2)*n1+beta2(2)*n2+beta3(2)*n3;
x8=O+beta1(1)*n1+beta2(2)*n2+beta3(2)*n3;
%%
%f=figure;hold on;set(gcf, 'color', 'white');
x=(x1(1)+x2(1)+x3(1)+x4(1)+x5(1)+x6(1)+x7(1)+x8(1))/8;y=(x1(2)+x2(2)+x3(2)+x4(2)+x5(2)+x6(2)+x7(2)+x8(2))/8;z=(x1(3)+x2(3)+x3(3)+x4(3)+x5(3)+x6(3)+x7(3)+x8(3))/8;
getPlotFramework(x1,x2,x3,x4,x5,x6,x7,x8)
axis tight





%view(0,0)
%view(-107,12)
view(36,28)
%axis equal
%delta=2;
%xlim([x-delta x+delta]);ylim([y-delta y+delta]);zlim([z-delta z+delta]);
%%
N=50;
for I=1231233 %N
  if 1==2
    switch I
      case 1
        x_p1=-0.5;y_p1=-0.5;z_p1=0.6;
        x_p2=1.25;y_p2=-0.5;z_p2=0.1;
      case 2
        x_p1=0.5;y_p1=-0.2;z_p1=0.7;
        x_p2=0.1;y_p2=0.2;z_p2=1.6;
      case 3
        x_p1=-0.5;y_p1=0;z_p1=0.5;
        x_p2=1.25;y_p2=2;z_p2=2;
      case 4
        x_p1=3;y_p1=0.5;z_p1=0.1;
        x_p2=0.1;y_p2=0.5;z_p2=2;
      case 5
        x_p2=3;y_p2=0.5;z_p2=0.1;
        x_p1=0.1;y_p1=0.5;z_p1=2;
      case 6
        x_p1=3.7;y_p1=2;z_p1=1;
        x_p2=0.2;y_p2=0;z_p2=-0.2;
      case 7
        x_p1=2.4;y_p1=2.3;z_p1=2.3;
        x_p2=3;y_p2=3;z_p2=3;
      case 8
        x_p1=2.3;y_p1=2.3;z_p1=2.3;
        x_p2=0.3;y_p2=2.3;z_p2=2.3;
      case 9
        x_p1=1;y_p1=0;z_p1=1;
        x_p2=1;y_p2=3;z_p2=1;
      otherwise
        return
    end
  end
  if I==1111
    x_p1=1;y_p1=0;z_p1=1;
    x_p2=2;y_p2=3;z_p2=1;
  end
  if I==2111
    x_p2=1;y_p2=0;z_p2=1;
    x_p1=2;y_p1=3;z_p1=1;
  end
  %      x_p1=1.3;y_p1=1.3;z_p1=4;
  %      x_p2=0.5+x_p1+cos(2*pi/(N)*I);y_p2=0.5+y_p1+sin(2*pi/(N)*I);z_p2=2.29;
  
  
  
  
  
  %%
  d=[x_p2-x_p1,y_p2-y_p1,z_p2-z_p1]';
  t_particle=sqrt(d'*d);
  d=d/sqrt(d'*d);
  q=quiver3(x_p1,y_p1,z_p1,d(1),d(2),d(3),0,'b','LineWidth',1,'MaxHeadSize',23);
  %adjust_quiver_arrowhead_size(q, 1.5)
  %  E=O
  %O=[x_p1,y_p1,z_p1];
  
  plot3([x_p1 x_p2],[y_p1 y_p2],[z_p1 z_p2],'ko-','LineWidth',2)
  quiver3(x_p1,y_p1,z_p1,x_p2-x_p1,y_p2-y_p1,z_p2-z_p1,0,'b','LineWidth',1)
  %   n1=x2-x1;n1=n1/sqrt(n1*n1');
  %   n2=x5-x1;n2=n2/sqrt(n2*n2');
  %   n3=x4-x1;n3=n3/sqrt(n3*n3');
  if I==1
    quiver3(O(1),O(2),O(3),n1(1),n1(2),n1(3),0,'k','LineWidth',2);text(O(1)+n1(1),O(2)+n1(2),(O(2)+n1(2))-0.2,'n_1');
    quiver3(O(1),O(2),O(3),n2(1),n2(2),n2(3),0,'k','LineWidth',2);text(O(1)+n2(1)-0.2,O(2)+n2(2),(O(2)+n2(3)),'n_2');
    quiver3(O(1),O(2),O(3),n3(1),n3(2),n3(3),0,'k','LineWidth',2);text(O(1)+n3(1),O(2)+n3(2),(O(2)+n3(3))-0.2,'n_3');
  end
  beta1min=beta1(1); %          0;
  beta1max=beta1(2); %            x2(1)-x1(1);
  beta2min=beta2(1); %            0;
  beta2max=beta2(2); %         x5(3)-x1(3);
  beta3min=beta3(1); %           0;
  beta3max=beta3(2); %         x4(2)-x1(2);
  %%
  dnk1=d'*n1;
  if dnk1 < 0
    disp('t1 is inverted')
    beta1=[beta1max beta1min];
  else
    beta1=[beta1min beta1max];
  end
  if abs(dnk1)<epsilon
    disp('dnk1=0: parallel to n1')
    dnk1=epsilon;
  end
  
  t1min=((O-E)'*n1+beta1(1))/dnk1;
  t1max=((O-E)'*n1+beta1(2))/dnk1;
  t1=((O-E)'*n1+[beta1min beta1max])/dnk1;
  Pt1min=E+t1min*d;
  plot3(Pt1min(1),Pt1min(2),Pt1min(3),'or','LineWidth',4)
  Pt1max=E+t1max*d;
  plot3(Pt1max(1),Pt1max(2),Pt1max(3),'or','LineWidth',4)
  %%
  dnk2=d'*n2;
  if dnk2 < 0
    disp('t2 is inverted')
    beta2=[beta2max beta2min];
  else
    beta2=[beta2min beta2max];
  end
  if abs(dnk2)<epsilon
    disp('dnk2=0: parallel to n2')
    dnk2=epsilon;
  end
  t2min=((O-E)'*n2+beta2(1))/dnk2;
  t2max=((O-E)'*n2+beta2(2))/dnk2;
  Pt2min=E+t2min*d;
  plot3(Pt2min(1),Pt2min(2),Pt2min(3),'og','LineWidth',4)
  Pt2max=E+t2max*d;
  plot3(Pt2max(1),Pt2max(2),Pt2max(3),'og','LineWidth',4)
  %%
  dnk3=d'*n3;
  if dnk3 < 0
    disp('t3 is inverted')
    beta3=[beta3max beta3min];
  else
    beta3=[beta3min beta3max];
  end
  if abs(dnk3)<epsilon
    disp('dnk3=0: parallel to n3')
    dnk3=1E-3;
  end
  t3min=((O-E)'*n3+beta3(1))/dnk3;
  t3max=((O-E)'*n3+beta3(2))/dnk3;
  Pt3min=E+t3min*d;
  plot3(Pt3min(1),Pt3min(2),Pt3min(3),'ob','LineWidth',4)
  Pt3max=E+t3max*d;
  plot3(Pt3max(1),Pt3max(2),Pt3max(3),'ob','LineWidth',4)
  %%
  if I==1
    set(0,'Units','pixels')
    scnsize = get(0,'ScreenSize');
    position = get(gcf,'Position');
    outerpos = get(gcf,'OuterPosition');
    borders = outerpos - position;
    edge = -borders(1)/2;
    if strcmp(monitor,'large'),pos1=[outerpos(1),scnsize(4),scnsize(3)/4-edge,scnsize(4)/1.25]; end;
    if strcmp(monitor,'small'),pos1=[outerpos(1)*0.8,scnsize(4)*(2/2),scnsize(3)/5-edge,scnsize(4)/1]; end;
    set(gcf,'OuterPosition',pos1)
  end
  %%
   axis tight
   
   axis equal
   
   
   
   
   %%
  f=figure;hold on;set(gcf, 'color', 'white');
  %g=subplot(10,1,8:9);hold on
  plot([t1min t1max],[0 0],'-or','LineWidth',2,'MarkerSize',10)
  plot([t2min t2max],[0.1 0.1],'-og','LineWidth',2,'MarkerSize',10)
  plot([t3min t3max],[0.2 0.2],'-ob','LineWidth',2,'MarkerSize',10)
  plot([t_particle t_particle],[0 0.2],'-ok','LineWidth',2,'MarkerSize',10)
  ylim([-0.5 0.5])
  
  
  %h=subplot(2,1,2);
  %set(h,'Position',[500 500 300 150]);
  %f = figure('Position',[200 200 400 150]);
  
  dat = [t1min t1max;...
    t2min t2max;...
    t3min t3max]';
  if I==I
    cnames = {'       t1 (red) n3/n2','     t2 (green) n1/n3','      t3 (blue) n1/n2'};
    rnames = {'Min','Max'};
    t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,...
      'RowName',rnames,'Position',[120 80 418 60]);
  end
  
  rechts=min(dat(2,:));
  links=max(dat(1,:));
  
  if links<=rechts
    plot([links links],[-0.1 0.3],'-k','LineWidth',2,'MarkerSize',10)
    plot([rechts rechts],[-0.1 0.3],'-k','LineWidth',2,'MarkerSize',10)
    %   if rechts==0
    %     if links ==0
    %     else
    %       xlim([0.1*links 1.5*rechts])
    %     end
    %   end
    if links<t_particle+epsilon && links+epsilon>0
      title('Particle HIT !!!')
      P=E+links*d;
      plot3(P(1),P(2),P(3),'+k','LineWidth',2,'MarkerSize',20)
      plot3(P(1),P(2),P(3),'ok','LineWidth',2,'MarkerSize',20)
      disp(['HIT!!!! Links=' num2str(links) ' rechts=' num2str(rechts) ' t_p= ' num2str(t_particle)])
    else
      title('particle trajectory HIT !!!')
    end
  else
    title('NO particle trajectory Intersection !!!')
  end
  
end


set( f,'toolbar','figure')
return
%%
if t1min>t2max & t2min<t1max
  if t3min>t2max
    disp('red in green')
  elseif t1max<t2min
    
  end
elseif t2min>min(t1) & t2min<max(t1)
  disp('green in red')
end

% intersection=1;
%
% if true(intersection)
%   title('HIT !!!')
% else
%   title('NO HIT !!!')
% end
% pos2 = [pos1(1) + edge,...
%         pos1(2)/7.5,...
%         pos1(3),...
%         pos1(4)];
% set(f,'OuterPosition',pos2)





end