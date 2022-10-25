#!/bin/bash

#==============================================================================
# title       : InstallPackages.sh
# description : This script installs the software packages required for
#               the module env scripts for creating a software environment for
#               PICLas/FLEXI code frameworks
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallPackages.sh
# notes       :
#==============================================================================

# Get OS Distributor ID
LSBNAME=$(lsb_release -si)

# Check OS
if [[ ${LSBNAME} != "Ubuntu" ]]; then
  echo "This script currently only supports Ubuntu. Exit."
  exit
fi

# Check for updates
sudo apt-get update

# compiler
sudo apt-get install make cmake cmake-curses-gui gfortran g++ gcc -y

# Boost C++ libraries development files
sudo apt-get install  libboost-dev -y

# cmake: developer's libraries for ncurses
sudo apt-get install  libncurses-dev  libncurses5-dev -y

# ctags
sudo apt-get install  exuberant-ctags -y

# blas and lapack
sudo apt-get install  libblas-dev  liblapack3  liblapack-dev -y

# tcl (required for module install scripts)
sudo apt-get install  tcl  tcl8.6-dev -y

# zlib is a library implementing the deflate compression method found in gzip and PKZIP (development)
sudo apt-get install   zlib1g-dev -y

# Secure Sockets Layer toolkit - development files
sudo apt-get install  libssl-dev -y
