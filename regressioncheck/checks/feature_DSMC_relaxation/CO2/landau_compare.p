## gnuplot-Skript für den Vergleich von Relaxationszeiten mit der Landau-Teller Gleichung (Anpassen der Spaltennummer $6,$7, etc.)
## Benutzung mit "gnuplot -> load './landau_compare.p'"

## Einlesen von Variablen am Ende und Anfang der Database:
Tvinf=1.*system("tail -4 Database.csv | head -1 | awk '{print $6}' | sed -e s/,$//")      # Vibrationstemperatur im Gleichgewicht
Xivinf=1.*system("tail -4 Database.csv | head -1 | awk '{print $7}' | sed -e s/,$//")     # Vib. Freiheitsgrad im Gleichgewicht
Tv0=1.*system("head -2 Database.csv | tail -1 | awk '{print $6}' | sed -e s/,$//")        # Vibrationstemperatur beim Start
Xiv0=1.*system("head -2 Database.csv | tail -1 | awk '{print $7}' | sed -e s/,$//")       # Vib. Freiheitsgrad beim Start
deltaT=1.*system("head -3 Database.csv | tail -1 | awk '{print $1}' | sed -e s/,$//")     # Zeitschritt
Pmean=1.*system("head -3 Database.csv | tail -1 | awk '{print $11}' | sed -e s/,$//")     # Mittlere Kollisionswahrscheinlichkeit

## Notwendig für den Vergleich mit einer zweitens Database
#Tv02=1.*system("head -2 Database_PDR.csv | tail -1 | awk '{print $6}' | sed -e s/,$//")
#Xiv02=1.*system("head -2 Database_PDR.csv | tail -1 | awk '{print $7}' | sed -e s/,$//")

## Benutzte Relaxationswahrscheinlichkeit (CO2, N2: 0.05, CH4: 0.02)
Pvib = 0.02

## Berechnung der normalisierten Temperaturdifferenz
f(x) = (Xivinf*Tvinf - x) / (Xivinf*Tvinf - Xiv0*Tv0)

## Nur notwendig für den Vergleich mit einer zweiten Database
#h(x) = (Xivinf*Tvinf - x) / (Xivinf*Tvinf - Xiv02*Tv02)

## Landau-Teller Gleichung
g(x) = exp(-x*Pmean*Pvib/deltaT)

set xlabel "Time [s]"
set ylabel "Normalized Energy Difference [-]"
set yr [0:1]
set format x "%g"
set xr [0:6e-6]
set xtics (0,1e-6,2e-6,3e-6,4e-6,5e-6,6e-6)

## Vergleich einer Database mit Landau-Teller
plot './Database.csv' u ($1):(f($6*$7)) w l title "PICLas", './Database.csv' u ($1):(g($1)) w l title "Landau-Teller"

## Vergleich mit anderen Database's
#plot './Database.csv' u ($1):(f($6*$7)) w l title "PICLas: Multi-mode relaxation", './Database_PDR.csv' u ($1):(h($6*$7)) w l title "PICLas: Prohibiting double relaxation", './Database.csv' u ($1):(g($1)) w l title "Landau-Teller"

## Ausgabe in eine Datei

#set terminal postscript eps size 3.5,2.62 enhanced color \
#    font 'Helvetica,16' linewidth 2
#set output "LT.eps"
#replot
#set terminal x11
#set size 1,1
