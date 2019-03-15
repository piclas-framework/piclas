#!/bin/bash

# Get MPF from parameter.ini 
filename=parameter.ini
while read line; do
  line="${line// /}"
  if [[ $line == *"="* ]]; then
    if [[ ${line:0:1} == *"!"* ]]; then
      continue
    else
      if [[ $line == *"Part-Species1-MacroParticleFactor"* ]]; then
        #echo $line
        MPF=$(echo $line | cut -d '=' -f 2 )
      fi
    fi
  fi
done < $filename

filenameA=reggie_MPF"$MPF"_merged_csv_files_full.csv
filenameB=reggie_MPF"$MPF"_merged_csv_files.csv

append_line () { temperature=$(echo $1 | sed -r 's/.*_Tvib_([0-9]*)\_.*/\1/g') ;
                 mytail=$(tail -n 1 $1) ;
                 echo $temperature", "$mytail >> $filenameA ;
               }


myheader=$(head -n 1 Database_Tvib_10000_ref.csv)
echo 'temperature'", "$myheader > $filenameA

append_line Database_Tvib_10000_ref.csv 
append_line Database_Tvib_15000_ref.csv 
append_line Database_Tvib_20000_ref.csv 
append_line Database_Tvib_25000_ref.csv 
append_line Database_Tvib_30000_ref.csv 


# shift 92 to 93 beacause of added temperature
cut -d ',' -f 1,93- $filenameA > $filenameB
cat $filenameB




