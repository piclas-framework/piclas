#!/usr/bin/perl -w 

if ($#ARGV == -1) {
  foreach my $file (`ls ./*_Solution_*.dat`) {
    chomp($file);
    $file = substr($file,0,-4);
    `cp ./writeDataBase.py ./tmp.py`;
    `perl -pi -e 's|FILENAME_PLACEHOLDER|$file|g' ./tmp.py`;
    system("~/Programme/visit2_4_1.linux-x86_64/bin/visit -cli -nowin -s ./tmp.py -assume_format Tecplot");
  }
} else {
    $file = $ARGV[0];
    $file = substr($file,0,-4);
    `cp ./writeDataBase.py ./tmp.py`;
    `perl -pi -e 's|FILENAME_PLACEHOLDER|$file|g' ./tmp.py`;
    system("~/Programme/visit2_4_1.linux-x86_64/bin/visit -cli -nowin -s ./tmp.py -assume_format Tecplot");
}
