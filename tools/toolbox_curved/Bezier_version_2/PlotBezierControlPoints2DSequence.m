function PlotBezierControlPoints2DSequence(P,N)
% ==================================================================
% Description
% ==================================================================
% check the sequence of the Bezier control points by plotting them and
% labelling them with their respective number

% ==================================================================
% Input
% ==================================================================
% sampPts           ::

% ==================================================================
% Output
% ==================================================================
%                   ::

% ==================================================================
% Function Start
% ==================================================================

K=1;
for q=1:N+1
  for p=1:N+1
    % plot the control points of the 2D bezier polynomial
    plot(P(p,q,1),P(p,q,2),'bs','MarkerSize',10,'LineWidth',1)
    text(P(p,q,1),P(p,q,2)*0.9,['P_{' num2str(p) ','  num2str(q) '}'])
    if K==1
      str=input('Enter "y" to break: ','s');
      if strcmp(str,'y'),K=0;end;
    end
    if and(q==N+1,p==N)
      K=0;
    end
  end
end

% ==================================================================
% Function End
% ==================================================================
end