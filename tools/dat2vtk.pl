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
$nParticles =~ s/I=(\d+),.+/$1/g;

my $nnewPart = 0;
open(OUTFILE,">$outFileName");
print OUTFILE "# vtk DataFile Version 3.1\n";
print OUTFILE "This file is generated from an HDF5 data file\n";
print OUTFILE "ASCII\n";
print OUTFILE "DATASET UNSTRUCTURED_GRID\n\n";
my $minNR = 4;
my $maxNR = $minNR+$nParticles+1;
if ($#ARGV == 0){
  print "nParticles=$nParticles\n";
  print OUTFILE "POINTS ".sprintf("%d",$nParticles)." FLOAT\n";
  close(OUTFILE);
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
} else {
  my $specnum = $ARGV[1];
  my @specarray = `awk '(NR>$minNR)&&(NR<$maxNR) {print \$7;}' $fileName`;
  my @sortspecarray = 0;
  my $i=0;
  for ($i=0;$i<@specarray;$i++){
    if ($specarray[$i] == $specnum) {
      $sortspecarray[$i] = 1;
      $nnewPart++;
    } else {
      $sortspecarray[$i] = 0;
    }
  }
  print "nParticles=$nnewPart\n";
  my @xarray =`awk '(NR>$minNR)&&(NR<$maxNR) {print \$1;}' $fileName`;
  my @yarray =`awk '(NR>$minNR)&&(NR<$maxNR) {print \$2;}' $fileName`;
  my @zarray =`awk '(NR>$minNR)&&(NR<$maxNR) {print \$3;}' $fileName`;
  print OUTFILE "POINTS ".sprintf("%d",$nnewPart)." FLOAT\n";
  for (my $i=0; $i<$nParticles; $i++){ 
    if ($sortspecarray[$i]==1){
      #print OUTFILE "$xarray[$i] $yarray[$i] $zarray[$i]";
      print OUTFILE sprintf("%.8f",$xarray[$i])." ".sprintf("%.8f",$yarray[$i])." ".sprintf("%.8f",$zarray[$i])."\n";  
    }
  }
  print OUTFILE "\n";
  print OUTFILE "CELLS ".sprintf("%d",$nnewPart)." ".sprintf("%d",2*$nnewPart)." \n";
  for (my $i=0; $i<$nnewPart; $i++){ print OUTFILE "1 $i\n"; }
  print OUTFILE "\n";
  print OUTFILE "CELL_TYPES ".sprintf("%d",$nnewPart)."\n";
  for (my $i=0; $i<$nnewPart; $i++){ print OUTFILE "1 "; }
  print OUTFILE "\n\n";
  print OUTFILE "POINT_DATA ".sprintf("%d",$nnewPart)."\n";
  print OUTFILE "SCALARS species FLOAT\n";
  print OUTFILE "LOOKUP_TABLE default\n";
  for ($i=0; $i<$nParticles; $i++){ 
    if ($sortspecarray[$i]==1){
      print OUTFILE "$specarray[$i]"; 
    }
  }
}
print OUTFILE "\n";
close(OUTFILE);
