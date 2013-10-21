#!/bin/bash
#
zeilen=$(cat $1 | wc -l)
for ((i=1; i<=$zeilen; i++ ))
do
  zeile=`sed -n $i'p' $1`
  laenge=`echo "$zeile" | wc -L`
  #echo "Zeile $zeile, Länge $laenge"
  if [ $laenge -gt 132 ]
  then
    echo "Achtung: Zeile $i ist $laenge Zeichen lang und damit über 132 Zeichen!"
    echo $zeile
    echo "---"
  fi
done
