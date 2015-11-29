function checkControlPoints(p,ControlPointsOriginal,N_plot)
% ==================================================================
% Description
% ==================================================================
% use control points to interpolate to a lagrange basis and use it the
% inverse of the Vandermonde to change the basis back to the bezier control
% points

% ==================================================================
% Input
% ==================================================================
% p                     :: polynomial order
% ControlPointsOriginal :: original control points which are to be checked
% N_plot                :: number of points for plotting

global Face
% ==================================================================
% Output
% ==================================================================
%                   ::

% ==================================================================
% Function Start
% ==================================================================
%dim=size(ControlPointsOriginal);p=dim(1)-1;
LagrangeValues=zeros(p+1,p+1,3);

%% the surface values are created for equidistant points !!! 
%  (however, for p=2 they coincide with CL nodes)
%  for higher degrees, copy the function and change to CL nodes
for dim=1:3
  for j=1:p+1
    [~,LagrangeValues(:,j,dim)]=SuperSampleBezierPolynomial(p,ControlPointsOriginal(:,j,dim),3);
  end
end
plot3(LagrangeValues(:,:,1),LagrangeValues(:,:,2),LagrangeValues(:,:,3),'+r','LineWidth',4,'MarkerSize',10);

%% get vandermonde (CL nodes !!!)
i=0:p;
Xi_NGeo=2/p.*i-1;
XiCL_NGeo=-cos(i/p*acos(-1.));
[Vdm_Bezier,sVdm_Bezier] = BuildBezierVdm(p,XiCL_NGeo);

%% change basis - get the control points
% Lagrange -> Bezier
P=ChangeBasis(p,LagrangeValues,sVdm_Bezier,3);

%% plot the bezier curve
test=figure;  hold on;set(gcf, 'color', 'white');view(29,60);grid on;xlabel('x');ylabel('y');zlabel('z');
PlotSuperSampledBezier3D(p,P,Face,N_plot,3);
plot3(LagrangeValues(:,:,1),LagrangeValues(:,:,2),LagrangeValues(:,:,3),'+m','LineWidth',4,'MarkerSize',10);
axis equal

% ==================================================================
% Function End
% ==================================================================
end