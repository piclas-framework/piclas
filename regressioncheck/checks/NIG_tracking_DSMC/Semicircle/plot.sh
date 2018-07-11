r=0.5
set parametric
set trange[0:pi]
fx(t)=r*cos(t)
fy(t)=r*sin(t)
r2=0.1
fx2(t)=r2*cos(t)
fy2(t)=r2*sin(t)
#set size equal
plot './ParticlePosition.csv' u 3:4 w l, './geo/geo.csv' u 1:2 w l lc 4, fx(t),fy(t) lc 3,fx2(t),fy2(t) lc 3
pause -1
reload
