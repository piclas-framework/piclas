These scripts help setting up an environment for development with piclas

------------------------------------------------------------
Order of how to setup environment

0. change dir to script directory ( cd /path/to/piclas/tools )

1. Install the basic packages depending on the OS
   - ./InstallPackagesUbuntu16.sh
   - ./InstallPackagesUbuntu20.sh
   - ./InstallPackagesUbuntu21.sh
   and if you have a server that has been setup with only basic packages, the following might be required
   sudo ./InstallPackagesServer.sh

2. open a new terminal and run:
      sudo -s
   or run the following script with sudo ./Install...

3. ./InstallModules.sh
  4. reboot and maybe second time ./InstallModules.sh is needed

5. ./InstallCMake.sh
6. ./InstallGCC.sh
7. ./InstallMPIallCOMPILERS.sh
8. ./InstallHDF5.sh

! next one might currently not work
9. ./InstallParaview.sh

------------------------------------------------------------
all scripts (5-8) can be rerun with "-r" or "-rerun" argument
  this cleans the created module file and build directory of the to-be-build version and rebuilds it
