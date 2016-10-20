set terminal gif animate
set output 'output.gif'
stats './ParticlePosition.csv' nooutput
r=0.5
set parametric
set trange[0:pi]
set size ratio -1
set xrange [-0.5:0.5]
set yrange [0.:0.5]
set nokey
fx(t)=r*cos(t)
fy(t)=r*sin(t)
r2=0.1
fx2(t)=r2*cos(t)
fy2(t)=r2*sin(t)

do for [i=1:int(STATS_blocks)] {
   plot for [IDX=0:i-1] './ParticlePosition.csv' index IDX u 3:4 w p lc rgb "red" pt 1, fx(t),fy(t) lc 3,fx2(t),fy2(t) lc 3
}
