function [smin,smax]=CalcSminSmax(J,NGeo,xi,BezierControlPoints1D,K,createPlots,create1Dplots)
global Face
if true(create1Dplots)
  f=figure;hold on;set(gcf, 'color', 'white');grid on;xlabel('x');ylabel('y');plot([-1 1],[0 0],'k-');
end
%%
if J==1
  for I=1:NGeo+1
    if true(create1Dplots)
      plot(xi,BezierControlPoints1D(:,I,Face),'ko-')
    end
    minmax(2,I)=max(BezierControlPoints1D(I,:,Face)); % Upper
    minmax(1,I)=min(BezierControlPoints1D(I,:,Face)); % Lower
  end
else
  for I=1:NGeo+1
    if true(create1Dplots)
      plot(xi,BezierControlPoints1D(I,:,Face),'ko-')
    end
    minmax(2,I)=max(BezierControlPoints1D(:,I,Face)); % Upper
    minmax(1,I)=min(BezierControlPoints1D(:,I,Face)); % Lower
  end
end

%%
smin= 1.5;
smax=-1.5;
%check streckenzug upper/lower
for I=1:NGeo
  if  minmax(2,I)*minmax(2,I+1)<=0 % Upper
    m=(minmax(2,I+1)-minmax(2,I))/(xi(I+1)-xi(I));
    smin_proposed=xi(I)-minmax(2,I)/m;
    smin=min(smin_proposed,smin);
  end
  if  minmax(1,I)*minmax(1,I+1)<=0 % Upper
    m=(minmax(1,I+1)-minmax(1,I))/(xi(I+1)-xi(I));
    smax_proposed=xi(I)-minmax(1,I)/m;
    smax=max(smax_proposed,smax);
  end
end
%% streckenzug:
if  minmax(2,1)*minmax(2,NGeo+1)<=0 %1. zeile anfang/ende -> upper
  m=(minmax(2,NGeo+1)-minmax(2,1))/(xi(NGeo+1)-xi(1));
  smin_proposed=xi(1)-minmax(2,1)/m;
  smin=min(smin_proposed,smin);
end
if  minmax(1,1)*minmax(1,NGeo+1)<=0 %2. zeile anfang/ende -> lower
  m=(minmax(1,NGeo+1)-minmax(1,1))/(xi(NGeo+1)-xi(1));
  smax_proposed=xi(1)-minmax(1,1)/m;
  smax=max(smax_proposed,smax);
end
%% streckenzug:
if  minmax(1,1)*minmax(2,1)<=0 % 1 spalte anfang/ende
  %if true: both values have different sign
  smin_proposed=xi(1);
  smin=min(smin_proposed,smin);
end
%
%if  minmax(1,2)*minmax(2,2)<=0 % 2 spalte anfang/ende -> ALT!!!!!!!!!!!!!
if  minmax(1,NGeo+1)*minmax(2,NGeo+1)<=0 % 2 spalte anfang/ende
  %if true: both values have different sign
  smax_proposed=xi(NGeo+1);
  smax=max(smax_proposed,smax);
end
%%
if smin== 1.5, smin=-1;end;
if smax==-1.5, smax= 1;end;
%       sarray(K,J,1)=smin;
%       sarray(K,J,2)=smax;
if true(create1Dplots)
  plot([smin smax],[0 0],'r+')
end
if and(K==1,true(createPlots))
  export_fig(['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/1D_projected_bezier_direction' num2str(J) '_projection_export_fig.pdf']);
end

end