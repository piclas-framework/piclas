function [Vdm_Bezier,sVdm_Bezier] = BuildBezierVdm(N_In,xi_In)
% ==================================================================
% Description
% ==================================================================
% create the vandermode matrix for Lagrage(CL-nodes) -> Bézier(equidistant)

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
Vdm_Bezier =zeros(N_In+1);
sVdm_Bezier=zeros(N_In+1);


% DO i=0,N_In
%   DO j=0,N_In
%     CALL BernsteinPolynomial(N_In,j,xi_in(i),Vdm_Bezier(i,j))
%     ! array with binomial coeffs for bezier clipping
%     IF(i.GE.j)THEN!only calculate LU (for n >= k, else 0)
%       arrayNchooseK(i,j)=REAL(CHOOSE(i,j))
%       FacNchooseK(i,j) = (1.0/(2.0**REAL(i)))*ArrayNchooseK(i,j)
%     ELSE
%       arrayNchooseK(i,j)=0.
%     END IF
% !    print*,'i,j',i,j,arrayNChooseK(i,j),FacNChooseK(i,j),'komisch',REAL(1.0/(2.0**Real(i)))
% !    read*
%   END DO !i
% END DO !j

for i=1:N_In+1
  for j=1:N_In+1
    Vdm_Bezier(i,j)=B(N_In,j-1,xi_In(i));
  end
end

% Invert the Vandermode matrix
sVdm_Bezier=Vdm_Bezier^(-1);

% check inversion:  dummy=SUM(ABS(MATMUL(sVdm_Bezier,Vdm_Bezier)))-REAL(N_In+1)
dummy=sum(sum(sVdm_Bezier*Vdm_Bezier))-(N_In+1);
if(abs(dummy)>1E-12), disp(['0 = ' num2str(dummy)]); error('vandermode was not inverted properly !?'); end;
% ==================================================================
% Function End
% ==================================================================
end