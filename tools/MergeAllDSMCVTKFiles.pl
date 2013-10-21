#!/usr/bin/perl -w
# Skript zum Mergen von Density (+Temp +...) vtk Ausgaben

# hier eines der verzeichnisse (sinnvollerweise das mit den meisten ausgaben) angeben und direkt drunter die buchstaben bis inkl zum \ zaehlen
$counter = 0;
foreach my $file (`ls DSMCOut_*_0001.vtk`) {
   $counter++;
}
foreach my $file1 (`ls DSMCOut_00000_*.vtk`) {
  chomp($file1);
  $number = substr($file1,14,4);
  $num = int $number;
  `./MergeDSMCVTKFiles $counter $num`;
  `cp DSMCOut_Merged.vtk DSMCOut_Merged$number.vtk`;
   print "$number done!\n";
  }
