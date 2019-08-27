These script help setting up an environment for development with piclas

------------------------------------------------------------
Order of how to setup environment

1. open console and run:
      sudo -i
2. change dir to script directory ( cd /path/to/piclas/tools )
3. ./InstallPackages.sh
4. ./InstallModules.sh
  5. reboot and maybe second time ./InstallModules.sh is needed
6. ./InstallCMake.sh
7. ./InstallGCC.sh
8. ./InstallMPIallCOMPILERS.sh
9. ./InstallHDF5.sh

! next one might currently not work
10. ./InstallParaview.sh

------------------------------------------------------------
all scripts (6-9) can be rerun with "-r" or "-rerun" argument
  this cleans the created module file and build directory of the to-be-build version and rebuilds it
