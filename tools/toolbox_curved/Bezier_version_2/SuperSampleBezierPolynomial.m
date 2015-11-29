function [Xi_sampPts Bezier]=SuperSampleBezierPolynomial(p,P,sampPts)
% ==================================================================
% Description
% ==================================================================
% Equidistantly super sample a Bézier polynomial of degree "p" from the
% control points "P" with "sampPts" number of sampling points.

% ==================================================================
% Input
% ==================================================================
% p                 :: polynomial degree
% P                 :: control points
% sampPts           :: number of sampling points

% ==================================================================
% Output
% ==================================================================
% Xi_sampPts        :: 
% Bezier            :: 


% ==================================================================
% Function Start
% ==================================================================
i=1:sampPts;
Xi_sampPts=(2/(sampPts-1).*(i-1)-1)';

Bezier=zeros(1,sampPts);

%h1=surf(P(:,:,1),P(:,:,2),P(:,:,3));set(h1,'FaceAlpha',0)
for I=1:sampPts
      for ip=1:p+1
        %disp(['p=' num2str(p) ' q=' num2str(q) ])
        %disp(['n=' num2str(N) ' k=' num2str(p-1) '(p)'])
        % disp(['n=' num2str(N) ' k=' num2str(q-1) '(q)'])
        %plot3(P(p,q,1),P(p,q,2),P(p,q,3),'bo','MarkerSize',10,'LineWidth',2)
        Bezier(I)=Bezier(I)+P(ip)*B(p,ip-1,Xi_sampPts(I));
      end
end

Bezier=Bezier';
% ==================================================================
% Function End
% ==================================================================
end
