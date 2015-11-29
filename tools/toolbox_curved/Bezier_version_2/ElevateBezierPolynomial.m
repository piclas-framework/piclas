function [Xi_NGeo_elevated P_elevated]=ElevateBezierPolynomial(p,P,t,ElevationMatrix)
% ==================================================================
% Description
% ==================================================================
% this function creates a new equidistantly distributed set of control
% points (bézier polynomial basis coefficients) based on the control points
% "P" with order "p" and elevates them to "P_elevated" on
% "Xi_NGeo_elevated" with order "p+t"
% ==================================================================
% Input
% ==================================================================
% p                 :: Polynomial order of the Bézier basis
% P                 :: Control points of the Bézier basis
% t                 :: Increase the polynomial order from "p" to "p+t"
% ElevationMatrix   :: precomputed elevation matrix "p" to "p+t"
% ==================================================================
% Output
% ==================================================================
% Xi_NGeo_elevated  :: Equidistantly distributed set of control points
% P_elevated        :: Elevated control points
% ==================================================================
% Function Start
% ==================================================================
% elevate the degree from p -> p+t
p_elevated=p+t;
P_elevated       = zeros(1,p_elevated+1);
i                = 1:p_elevated+1;
Xi_NGeo_elevated = (2/(p_elevated).*(i-1)-1)';

% algorithm originally from "The NURBS Book" by Les Piegl, Wayne Tiller (p.205)
% the first and last points remain the same!
P_elevated(1)    = P(1);
P_elevated(p+t+1)= P(end);

% P_2,...,P_p
for i=1:p+t-1
  jStart=max(0,i-t);
  jEnd=min(p,i);
  for j=jStart:jEnd
    % old
    %P_elevated(i+1)=P_elevated(i+1)+...
    %                                   choose(p,j)*choose(t,i-j)*P(j+1)/...
    %                                           choose(p_elevated,i);
    % new
    P_elevated(i+1)=P_elevated(i+1)+ElevationMatrix(i+1,j+1)*P(j+1);
  end
end
% ==================================================================
% Function End
% ==================================================================
end