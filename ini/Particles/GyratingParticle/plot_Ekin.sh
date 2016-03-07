set xlabel 't / ns'
set ylabel 'E_kin / J'
#set yrange[0:3e-14]
set key horizontal
set key top
#plot './implicit/Database.csv' u (1e9*$1):4 w l t 'implicit', './explicit/Database.csv' u (1e9*$1):4 w l t 'explicit'
plot './Database.csv' u (1e0*$1):4 w l t 'E_kin'
pause 1
replot
reread
