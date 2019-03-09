#!/bin/bash

append_line () { temperature=$(echo $1 | sed -r 's/.*_Tvib_([0-9]*)\_.*/\1/g') ;
                 #cut -d ',' -f 92- $1 > temp.csv ;
                 #mytail=$(tail -n 1 temp.csv) ;
                 mytail=$(tail -n 1 $1) ;
                 echo $temperature", "$mytail >> reggie_merged_csv_files_full.csv ;
               }


#cut -d ',' -f 92- Database_Tvib_10000_ref.csv > temp.csv
#myheader=$(head -n 1 temp.csv)
myheader=$(head -n 1 Database_Tvib_10000_ref.csv)
echo 'temperature'", "$myheader > reggie_merged_csv_files_full.csv

append_line Database_Tvib_10000_ref.csv 
append_line Database_Tvib_15000_ref.csv 
append_line Database_Tvib_20000_ref.csv 
append_line Database_Tvib_25000_ref.csv 
append_line Database_Tvib_30000_ref.csv 


# shift 92 to 93 beacause of added temperature
cut -d ',' -f 1,93- reggie_merged_csv_files_full.csv > reggie_merged_csv_files.csv
cat reggie_merged_csv_files.csv

# clean-up
rm temp.csv
