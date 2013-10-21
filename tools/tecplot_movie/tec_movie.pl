#!/usr/bin/perl -w 

foreach my $file (`ls ./*_Solution_*.plt`) {
  chomp($file);
  `cp ./movie.mcr ./tmp.mcr`;
  `cp ./tecmovie.lay ./tmp.lay`;
  $current_folder = `pwd`;
  chomp($current_folder);  
  `perl -pi -e 's|FILE_PLACEHOLDER|$file|g' ./tmp.lay`;
  `perl -pi -e 's|FOLDER_PLACEHOLDER|$current_folder|g' ./tmp.mcr`;
  $filetemp = substr($file,0,-4);
  $filetemp = substr($filetemp,2);
  `perl -pi -e 's|OUTPUT_PLACEHOLDER|$filetemp|g' ./tmp.mcr`;
  system("tecplot -p ./tmp.mcr");
}

`mencoder "mf://*Solution*.png" -mf fps=6 -o video.avi -ovc lavc -lavcopts vcodec=wmv2`;
