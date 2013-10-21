#!/usr/bin/perl -w 

if ($#ARGV == -1) {
  foreach my $file (`ls ./*_Solution_*_Proc000000.dat`) {
    chomp($file);
   $placeholder_tecmerge=".tec\" DATA_PLACEHOLDER";
  `cp ./convTec09ToTec06.mcr ./tmp.mcr`;
   $current_folder = `pwd`;
   chomp($current_folder);  
  `perl -pi -e 's|FOLDER_PLACEHOLDER|$current_folder|g' ./tmp.mcr`;
    $filetmp = substr($file,0,-10);
    $file = substr($file,0,-10);
    foreach my $fileperproc (`ls ./$filetmp*.dat`) {
      chomp($fileperproc);
      $fileperproc = substr($fileperproc,0,-4);
      `cp ./convTec09ToTec06.py ./tmp.py`;
      `perl -pi -e 's|FILENAME_PLACEHOLDER|$fileperproc|g' ./tmp.py`;
      system("visit -cli -nowin -s ./tmp.py -assume_format Tecplot");
      $fileperproc = substr($fileperproc,1);
      $tempdata = "\"".$current_folder.$fileperproc.$placeholder_tecmerge;
      `perl -pi -e 's|DATA_PLACEHOLDER|$tempdata|g' ./tmp.mcr`;
    } 
  $tempdata = "";
  `perl -pi -e 's|DATA_PLACEHOLDER|$tempdata|g' ./tmp.mcr`;
  $file = substr($file,0,-5);
  $file = "/".$file;
  `perl -pi -e 's|OUTPUT_PLACEHOLDER|$file|g' ./tmp.mcr`;
  system("tecplot -b -p ./tmp.mcr");
  $file = substr($file,1);
  $file =$file."*.tec";
  system("rm $file");
  }
} else {
  $file = $ARGV[0];
  $file = substr($file,0,-10);
  $MFBD="|MFBD|";
  $placeholder_tecmerge=".tec\" DATA_PLACEHOLDER";
  `cp ./convTec09ToTec06.mcr ./tmp.mcr`;
   $current_folder = `pwd`;
   chomp($current_folder);  
  `perl -pi -e 's|FOLDER_PLACEHOLDER|$current_folder|g' ./tmp.mcr`;
  foreach my $fileperproc (`ls ./$file*.dat`) {
    chomp($fileperproc);
    $fileperproc = substr($fileperproc,0,-4);
    `cp ./convTec09ToTec06.py ./tmp.py`;
    `perl -pi -e 's|FILENAME_PLACEHOLDER|$fileperproc|g' ./tmp.py`;
    system("visit -cli -nowin -s ./tmp.py -assume_format Tecplot");
    $fileperproc = substr($fileperproc,1);
    $tempdata = "\"".$current_folder.$fileperproc.$placeholder_tecmerge;
    `perl -pi -e 's|DATA_PLACEHOLDER|$tempdata|g' ./tmp.mcr`;
    #`perl -pi -e  "DATA_PLACEHOLDER" $tempdata ./tmp.mcr`;
  }  
  $tempdata = "";
  `perl -pi -e 's|DATA_PLACEHOLDER|$tempdata|g' ./tmp.mcr`;
  $file = substr($file,0,-5);
  $file = "/".$file;
  `perl -pi -e 's|OUTPUT_PLACEHOLDER|$file|g' ./tmp.mcr`;
  system("tecplot -b -p ./tmp.mcr");
  $file = substr($file,1);
  $file =$file."*.tec";
  system("rm $file");
}
