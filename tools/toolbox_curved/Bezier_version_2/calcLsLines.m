function  [v1,v2]=calcLsLines(P)
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
global n_Ls
global n_Ls_inv
global projectedFace

%% 1st direction  (eigentlich eta)
v1=P(1,end,:)-P(1,1,:);
v2=P(end,end,:)-P(end,1,:);
n_Ls=v1+v2;
n_Ls=(1/sqrt([n_Ls(1) n_Ls(2)]*[n_Ls(1);n_Ls(2)]))*n_Ls;
%plot([0 n_Ls(1)*0.25], [0 n_Ls(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))

%% 2nd direction (eigentlich xi)
u1=P(end,1,:)-P(1,1,:);
u2=P(end,end,:)-P(1,end,:);
n_Ls_inv=u1+u2;
n_Ls_inv=(1/sqrt([n_Ls_inv(1) n_Ls_inv(2)]*[n_Ls_inv(1);n_Ls_inv(2)]))*n_Ls_inv;
%plot([0 n_Ls_inv(1)*0.25], [0 n_Ls_inv(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))


%% xi direction
v1=(P(end,1,:)-P(1,1,:)) + (P(end,end,:)-P(1,end,:));
v1=reshape(v1,2,1);
v1=v1/sqrt(v1'*v1);
% plot([0 v1(1)*0.25], [0 v1(2)*0.25],'r-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
% text(v1(1)*0.25,v1(2)*0.25,'v1','Parent',get(projectedFace, 'children'))

%% eta direction
v2=(P(1,end,:)-P(1,1,:)) + (P(end,end,:)-P(end,1,:));
v2=reshape(v2,2,1);
v2=v2/sqrt(v2'*v2);
% plot([0 v2(1)*0.25], [0 v2(2)*0.25],'k-','LineWidth',2,'MarkerSize',10,'Parent',get(projectedFace, 'children'))
% text(v2(1)*0.25,v2(2)*0.25,'v2','Parent',get(projectedFace, 'children'))
end






