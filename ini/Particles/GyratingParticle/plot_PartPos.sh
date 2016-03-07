set xlabel 'x / m'
set ylabel 'y / m'
#set xrange[-0.6:0.6]
#set yrange[-0.6:0.6]
set xrange[-0.005:0.005]
set yrange[-0.005:0.005]
set xtics 0.001
set ytics 0.001
set size square
set grid
set key horizontal
set key top
#plot './implicit/Database.csv' u (1e9*$1):4 w l t 'implicit', './explicit/Database.csv' u (1e9*$1):4 w l t 'explicit'
plot './ParticlePosition.csv' u 3:4 w l
#plot './MyRank000001_ParticlePosition.csv' u 3:4 w l
pause 1
replot
reread
