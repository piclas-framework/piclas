set xlabel 't / ps'
set ylabel 'v / m/s'
#set xrange[-0.6:0.6]
#set yrange[-0.6:0.6]
#set xrange[-0.005:0.005]
#set yrange[-0.005:0.005]
#set xtics 0.1
#set ytics 0.1
set size square
set grid
set key horizontal
set key top
#plot './implicit/Database.csv' u (1e9*$1):4 w l t 'implicit', './explicit/Database.csv' u (1e9*$1):4 w l t 'explicit'
#plot './ParticlePosition.csv' u ($1*1e12):($3**2+$4**2+$5**2) w l
plot './ParticlePosition.csv' u ($1*1e12):($6) w l, \
     './ParticlePosition.csv' u ($1*1e12):($7) w l, \
     './ParticlePosition.csv' u ($1*1e12):($8) w l  
pause  1
replot
reread
