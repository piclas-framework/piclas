function BezierControlPoints2D=CalcBezierControlPoints2D(dimensions,nSides,NGeo,d1ij,d2ij,N,N_plot,createPlots)
global Face
global projectedFace

%%
BezierControlPoints2D=zeros(dimensions(1),dimensions(2),2,dimensions(4));
for SideID=1:nSides
  for q=1:NGeo+1
    for p=1:NGeo+1
      BezierControlPoints2D(p,q,:,SideID)=[d1ij(p,q,SideID) d2ij(p,q,SideID)];
    end
  end
end
projectedFace = figure;  hold on;set(gcf, 'color', 'white');


%% plot 2D projection with control points
PlotSuperSampledBezier3D(N,BezierControlPoints2D(:,:,:,Face),0,N_plot,2);


if true(createPlots)
  export_fig(['/home/stephen/Documents/Projekte/CurvedPIC/Paper/graphics/2D_projected_bezier_export_fig.pdf']);
end


end