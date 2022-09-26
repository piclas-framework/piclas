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

# Boost C++ libraries development files
sudo apt-get install  libboost-dev -y

# cmake
sudo apt-get install  libncurses-dev -y

# ctags
sudo apt-get install  exuberant-ctags -y

# libphonon-dev (multimedia framework from KDE using Qt 4 - development files)
sudo apt-get install  libphonon-dev -y

# make
sudo apt-get install  make -y

# python
sudo apt-get install  libpython-stdlib  libpython-dev -y

# blas and lapack
sudo apt-get install  libblas-dev  liblapack-dev -y

# tcl (required for module install scripts)
sudo apt-get install  tcl  tcl8.6-dev -y

# mesa (Mesa-specific OpenGL extensions)
sudo apt-get install  mesa-common-dev  libgl1-mesa-dri -y

# libz (compression library - development)
sudo apt-get install  libz-dev -y

# libxt (toolkit intrinsics library - development headers)
sudo apt-get install  libxt-dev -y

# Binutils are a collection of binary tools
sudo apt-get install  binutils -y

# GNU C Library: Development Libraries and Header Files
sudo apt-get install  libc6-dev -y

# Reference for all the packages needed to compile a Debian package.
# It generally includes the GCC/g++ compilers and libraries and some other utilities.
sudo apt-get install  build-essential -y

# multiple precision complex floating-point library development package
sudo apt-get install  libmpc-dev -y

# Multiprecision arithmetic library developers tools
sudo apt-get install  libgmp3-dev -y

# multiple precision floating-point computation developers tools
sudo apt-get install  libmpfr-dev  libmpfr-doc  libmpfr6 -y
