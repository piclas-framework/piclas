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

# Check for updates
sudo apt-get update

# filesystem operations
sudo apt-get install  libboost-dev

# cmake
sudo apt-get install  libncurses-dev

# ctags 
sudo apt-get install  ctags

# paraview
sudo apt-get install  qt4-qmake  qt4-dev-tools  libqt4-opengl-dev  libqt4-dev 

# libphonon-dev (multimedia framework from KDE using Qt 4 - development files)
sudo apt-get install  libphonon-dev

# make
sudo apt-get install  make

# python
sudo apt-get install  libpython-stdlib  libpython-dev

# blas and lapack
sudo apt-get install  libblas-dev  liblapack-dev

# tcl (required for module install scripts)
sudo apt-get install  tcl  tcl8.6-dev

# mesa (Mesa-specific OpenGL extensions)
sudo apt-get install  mesa-common-dev  libgl1-mesa-dri

# libz (compression library - development)
sudo apt-get install  libz-dev

# libxt (toolkit intrinsics library - development headers)
sudo apt-get install  libxt-dev

# Binutils are a collection of binary tools
sudo apt-get install  binutils

# GNU C Library: Development Libraries and Header Files
sudo apt-get install  libc6-dev

# Reference for all the packages needed to compile a Debian package.
# It generally includes the GCC/g++ compilers and libraries and some other utilities.
sudo apt-get install  build-essential

# multiple precision complex floating-point library development package
sudo apt-get install  libmpc-dev

# Multiprecision arithmetic library developers tools
sudo apt-get install  libgmp3-dev

# multiple precision floating-point computation developers tools
sudo apt-get install  libmpfr-dev  libmpfr-doc  libmpfr6

# Create Installsources directory and copy the module templates
sudo mkdir -p /opt/Installsources
sudo cp -r moduletemplates /opt/Installsources/moduletemplates
