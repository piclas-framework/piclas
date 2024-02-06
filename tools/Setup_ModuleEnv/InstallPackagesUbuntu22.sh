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
sudo apt update

# compiler
sudo apt install cmake -y
sudo apt install cmake-curses-gui -y
sudo apt install gfortran -y
sudo apt install g++ -y

# Boost C++ libraries development files
sudo apt install  libboost-dev -y

# cmake: developer's libraries for ncurses
sudo apt install libncurses-dev -y
sudo apt install libncurses5-dev -y

# ctags
sudo apt install exuberant-ctags -y

# blas and lapack
sudo apt install libblas-dev -y
sudo apt install liblapack-dev -y

# tcl (required for module install scripts)
sudo apt install tcl8.6-dev -y

# Secure Sockets Layer toolkit - development files
sudo apt install libssl-dev -y

# git
sudo apt install git
