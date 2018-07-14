clear all
close all


N=3
first=1
minmax=zeros(N,2);
xeq=linspace(-1,1,N);
% 
%minmax(:,1) = [0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];
%minmax(:,2) = [0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];

% % upper
% minmax(:,1)=[-0.05,-0.05,-0.05,-0.05,-0.05];
% minmax(:,2) = [0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];
% 
% % lower
% minmax(:,2)=[0.05,0.05,0.05,0.05,0.05];
% minmax(:,1)=[0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];

% 
% minmax(:,2) = [0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];
% minmax(:,1) = [0.003009825069583,-0.004454623041696,-0.004462530863543,0.003724950818699,0.018604733779872];
% minmax = -minmax;

%minmax(:,1)=[-2.8041277543718753E-003,-3.2909417463801694E-003,-3.7777557383888421E-003,-4.2645697303970755E-003,-4.7513837224056991E-003];
%minmax(:,2)=[ 1.8557442802926811E-003, 1.3689302882844891E-003, 8.8211629627565013E-004, 3.9530230426754014E-004,-9.1511687741141013E-005];

% check for upper
% minmax(:,1)=[-2,-2,-2,-2,-2];
% minmax(:,2)=[2,1.3,0.4,0.1,-2];

% check for lower
% minmax(:,2)=[3,3,3,3,3];
% minmax(:,1)=[2,1.3,0.4,0.1,-2];

% nishita fig.12
minmax(:,2)=[-8,2,7]
minmax(:,1)=[-10,-0.5,5];



plot(xeq,minmax(:,1),xeq,minmax(:,2));
hold on
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis

% now check the bounding box
smin=1.5;
smax=-1.5;

% convex-hull check by
% 1) check upper line for an intersection with x-axis
% 2) check lower line for an intersection with x-axis
% 3) check if both ends  have an intersection

% check upper/lower line for an intersection with the x-axis
for i=first:N
    for j=first:N
        if i==j
            continue
        end
        % check for upper line
        if(minmax(i,2)*minmax(j,2)<0)
            m=(minmax(j,2)-minmax(i,2))/(xeq(j)-xeq(i));
            tmp=xeq(j)-minmax(j,2)/m;
            if(m<0) % move right boundary
                smax=max(smax,tmp);
            else % move left boundary
                smin=min(smin,tmp);
            end
        end
        % check for lower line
        if(minmax(i,1)*minmax(j,1)<0)
            m=(minmax(j,1)-minmax(i,1))/(xeq(j)-xeq(i));
            tmp=xeq(j)-minmax(j,1)/m;
            if(m>0) % move right boundary
                smax=max(smax,tmp);
            else % move left boundary
                smin=min(smin,tmp);
            end
        end       
    end
end

% 2) get min element and index
[MinVal,Index] = min(minmax(:,1));


% 3) check left and right boundary
if(minmax(first,1)*minmax(first,2)<0.)
    smin=-1.;
end
if(minmax(N,1)*minmax(N,2)<0.)
    smax=1.;
end

Xmin=smin
Xmax=smax
