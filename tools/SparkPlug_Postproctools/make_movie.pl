#!/usr/bin/perl -w
# Skript zur erstellung von Bildern in Visit. Wenn man eine neue session für das Skript abspeichert, musst man darin noch den Dateinamen in FILENAME_PLACEHOLDER ändern.

foreach my $file (`ls ../StateBackup/*_State_*.h5`) {
  chomp($file);
  $file =~ s/\.h5$/.vtk/;
  (my $cgnsFile = $file) =~ s/\.\.\/StateBackup\///g;
  $cgnsFile =~ s/\.vtk/\.cgns/;
  $cgnsFile =~ s/State/visu/;
  (my $time = $cgnsFile) =~ s/^(.+)0(0\.\d+)\.cgns$/$2/;
  $time *= 1.0E+12;
  my $timeString = sprintf("%3.0f",$time);
  print "\$file=$file -> $cgnsFile at time $time | $timeString\n";
  `cp ./movie_sample.session ./tmp.session`;
  `cp ./movie_sample.session.gui ./tmp.session.gui`;
  `perl -pi -e 's|TIME_PLACEHOLDER|$timeString|g' ./tmp.session`;
  `perl -pi -e 's|VTK_FILENAME_PLACEHOLDER|$file|g' ./tmp.session`;
  `perl -pi -e 's/CGNS_FILENAME_PLACEHOLDER/$cgnsFile/g' ./tmp.session`;
  `perl -pi -e 's|VTK_FILENAME_PLACEHOLDER|$file|g' ./tmp.session.gui`;
  `perl -pi -e 's/CGNS_FILENAME_PLACEHOLDER/$cgnsFile/g' ./tmp.session.gui`;
  print "*  Visit starten (Skript-Modus)...\n";
  system("~/scratch/visit/bin/visit -nowin -cli -s ./make_movie.py");
  print "*  ...fertig!\n";
}

#foreach my $tifFile (`ls *.tif | grep visit`) {
#  chomp($tifFile);
#  (my $pngFile = $tifFile) =~ s/\.tif/.png/;
#  print "convert $tifFile $pngFile\n";
#  `convert $tifFile $pngFile`;
#}

#`mencoder "mf://visit*.png" -mf fps=6 -o video.avi -ovc lavc -lavcopts vcodec=wmv2`;
