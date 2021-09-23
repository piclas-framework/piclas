r=0.5
set parametric
set trange[0:pi]
fx(t)=r*cos(t)
fy(t)=r*sin(t)
r2=0.1
fx2(t)=r2*cos(t)
fy2(t)=r2*sin(t)
#set size equal
plot './MyRank000000_ParticlePosition.csv' u 3:4 w l, './geo/geo.csv' u 1:2 w l lc 3, fx(t),fy(t) lc 3,fx2(t),fy2(t) lc 3
pause -1
