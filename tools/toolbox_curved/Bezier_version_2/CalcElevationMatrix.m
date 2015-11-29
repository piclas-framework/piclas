function [M]=CalcElevationMatrix(p,t)
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

% ==================================================================
% Output
% ==================================================================
% Xi_NGeo_elevated  :: Equidistantly distributed set of control points
% P_elevated        :: Elevated control points


% ==================================================================
% Function Start
% ==================================================================
M=zeros(p+t+1,p+1);

M(1,1)    =1;
M(end,end)=1;

for i=1:p+t-1% from 0+1 to p_new-1 -> remove the edge points
  jStart=max(0,i-t);
  jEnd=min(p,i);
  
  for j=jStart:jEnd
    %P_elevated(i)=P_elevated(i)+choose(p,j)*choose(t,i-j)*P(j)/choose(p+t,i);  % original gleichung, d.h., von 0 bis p
    M(i+1,j+1) = choose(p,j)*choose(t,i-j) / choose(p+t,i);
%     choose(p,j)
%     choose(t,i-j)
%     choose(p+t,i)
%     clc
  end
  if(abs(sum(M(i+1,:))-1)>1e-12),error('The line of the elevation matrix does not sum to unity! exit.');exit;end;
end



% ==================================================================
% Function End
% ==================================================================
end