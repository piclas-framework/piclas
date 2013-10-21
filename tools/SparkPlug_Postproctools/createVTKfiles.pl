#!/usr/bin/perl -w
use strict;

foreach my $file (`ls ../StateBackup/zuendkerze_bgg_fine_State_00.0000*60.h5`) {
  chomp($file);
  print "processing file $file\n";
  `../../Zuendkerze_Bosch/extractParticleData $file`;
  `../../Zuendkerze_Bosch/dat2vtk.pl $file.dat`;
}
