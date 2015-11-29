% ==================================================================
% Description
% ==================================================================
% uses a pre-defined polynomial to create a bezier basis from a lagrange
% basis (CL nodes to equidistant nodes): create the Vandermode and use it
% do derive the control points

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
%541 DO iGP=0,N_in
%542   xGP(iGP)=-COS(iGP/REAL(N_in)*ACOS(-1.))
%543 END DO
close all; clear; clc

pMax=5;

YLim=[-0.3 3.5];

%Polynomial = @(x) -x^2+1.2;
Polynomial = @(x) -1.5*x^4+2*x^2+1;

for p=pMax:pMax
  
  %% 1. plot (bezier and lagrange basis + original polynomial)
  figure;set(gcf, 'color', 'white');
  subplot(1,2,1);hold on;
  %% calc nodes (equidistant and CL distributed)
  i=0:p;
  Xi_NGeo=2/p.*i-1;
  XiCL_NGeo=-cos(i/p*acos(-1.));
  
  for i=1:length(XiCL_NGeo)
    XCL_NGeo(i)=Polynomial(XiCL_NGeo(i));
  end
  for i=1:length(XiCL_NGeo)
    plot(XiCL_NGeo(i)*[1 1],XCL_NGeo(i)*[0 1],':k');
  end
  plotXi_NGeo=plot(Xi_NGeo,ones(1,length(Xi_NGeo))*(p-pMax),'bs');
  plotXiCL_NGeo=plot(XiCL_NGeo,ones(1,length(XiCL_NGeo))*(p-pMax),'r.');
  
  
  plotXCL_NGeo=plot(XiCL_NGeo,XCL_NGeo,'ro');
  
  
  
  
  fplot(Polynomial,[-1,1],'-k');  plotPolynomial=plot(0,1,'k-');
  ylim(YLim);
  %[Xi_NGeo' XiCL_NGeo']
  
  %% get vandermode
  [Vdm_Bezier,sVdm_Bezier] = BuildBezierVdm(p,XiCL_NGeo);
  P=zeros(1,length(XCL_NGeo));
  
  %% calculate control points: matrix-vector multiplication
%   for i=1:length(XCL_NGeo)
%     P(i)=sVdm_Bezier(i,:)*XCL_NGeo';
%   end
  P=matmul(sVdm_Bezier,XCL_NGeo); % P = V^-1 * X
  
  
  plot([-1 1],[0 0],'k-');        % plot 0-line
  plotCP=plot(Xi_NGeo,P,'^b--');  % plot control points
  
  title(['polynomial degree: p=' num2str(p)])
  legend([plotXi_NGeo plotXiCL_NGeo plotPolynomial plotXCL_NGeo plotCP],...
    'Equidistant points',...
    'Lobatto-Chebyshev points',...
    'Polynomial function',...
    'Lagrange nodes',...
    'Bezier control points')
  
  
  
  %% check control points (do they reproduce the polynomial?)
  [X Y]=SuperSampleBezierPolynomial(p,P,40);
  plot(X,Y,'r--');
  
  %% 2. plot (Elevation)
  %figure;hold on;set(gcf, 'color', 'white');
  subplot(1,2,2);hold on;
  plotCP=plot(Xi_NGeo,P,'^b--');
  % elevation
  elevation=20;
  [Xi_NGeo_elevated,P_elevated]=ElevateBezierPolynomial(p,P,elevation);
  plotCP_elevated=plot(Xi_NGeo_elevated,P_elevated,'^m--');
  
  
  % check elevated control points (do they reproduce the polynomial?)
  [X Y]=SuperSampleBezierPolynomial(p,P,40);
  fplot(Polynomial,[-1,1],'-k');
  plotPolynomial2=plot(0,1,'k-');
  %plot(X,Y,'r--'); %test is the new control points reproduce the plonymoial
  ylim(YLim);
  title(['polynomial degree elevated from p=' num2str(p) ' to p_2=' num2str(p+elevation)])
  legend([  plotPolynomial2  plotCP plotCP_elevated],...
    'Polynomial function',...
    'Bezier control points',...
    'Bezier control points - elevated')
  
end
disp('done.')
% ==================================================================
% Function End
% ==================================================================
