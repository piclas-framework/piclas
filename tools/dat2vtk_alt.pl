#!/usr/bin/perl -w
# Aufruf ./dat2vtl.pl myFile.dat
use strict;

my $fileNameOrig = $ARGV[0];
chomp($fileNameOrig);
(my $fileName = $fileNameOrig) =~ s/\.h5//g;        # strip ".h5"
(my $outFileName = $fileName)  =~ s/\.dat$/.vtk/g;  # replace ".dat" by ".vtk"

if ($fileNameOrig =~ /\.h5/) { `mv $fileNameOrig $fileName`; }

my $nParticles = `awk 'NR==4 {print;}' $fileName`;
chomp($nParticles);
$nParticles =~ s/I=\s*(\d+)[,\s].+/$1/g;
print "nParticles=$nParticles\n";

open(OUTFILE,">$outFileName");
print OUTFILE "# vtk DataFile Version 3.1\n";
print OUTFILE "This file is generated from an HDF5 data file\n";
print OUTFILE "ASCII\n";
print OUTFILE "DATASET UNSTRUCTURED_GRID\n\n";
print OUTFILE "POINTS ".sprintf("%d",$nParticles)." FLOAT\n";
close(OUTFILE);
my $minNR = 4;
my $maxNR = $minNR+$nParticles+1;
`awk '(NR>$minNR)&&(NR<$maxNR) {print \$1,\$2,\$3;}' $fileName >> $outFileName`;
open(OUTFILE,">>$outFileName");
print OUTFILE "\n";
print OUTFILE "CELLS ".sprintf("%d",$nParticles)." ".sprintf("%d",2*$nParticles)." \n";
for (my $i=0; $i<$nParticles; $i++){ print OUTFILE "1 $i\n"; }
print OUTFILE "\n";
print OUTFILE "CELL_TYPES ".sprintf("%d",$nParticles)."\n";
for (my $i=0; $i<$nParticles; $i++){ print OUTFILE "1 "; }
print OUTFILE "\n\n";
print OUTFILE "POINT_DATA ".sprintf("%d",$nParticles)."\n";
print OUTFILE "SCALARS species FLOAT\n";
print OUTFILE "LOOKUP_TABLE default\n";
close(OUTFILE);
`awk '(NR>$minNR)&&(NR<$maxNR) {print \$7;}' $fileName >> $outFileName`;
open(OUTFILE,">>$outFileName");
print OUTFILE "\n";
print OUTFILE "SCALARS vx FLOAT\n";
print OUTFILE "LOOKUP_TABLE default\n";
close(OUTFILE);
`awk '(NR>$minNR)&&(NR<$maxNR) {print \$4;}' $fileName >> $outFileName`;
open(OUTFILE,">>$outFileName");
print OUTFILE "\n";
print OUTFILE "SCALARS vy FLOAT\n";
print OUTFILE "LOOKUP_TABLE default\n";
close(OUTFILE);
`awk '(NR>$minNR)&&(NR<$maxNR) {print \$5;}' $fileName >> $outFileName`;
open(OUTFILE,">>$outFileName");
print OUTFILE "\n";
print OUTFILE "SCALARS vz FLOAT\n";
print OUTFILE "LOOKUP_TABLE default\n";
close(OUTFILE);
`awk '(NR>$minNR)&&(NR<$maxNR) {print \$6;}' $fileName >> $outFileName`;
open(OUTFILE,">>$outFileName");
print OUTFILE "\n";
close(OUTFILE);
