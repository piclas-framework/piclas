#!/bin/bash

#==============================================================================
# title       : InstallPackagesReggie.sh
# description : This script installs the software packages required for the
#               regression check tool "reggie2.0"
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallPackagesReggie.sh
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

# python
sudo apt-get install python python-numpy python3 python3-numpy python3-matplotlib python-matplotlib -y

# editors
sudo apt-get install gedit vim vim-runtime vim-common vim-tiny gdb vim-gui-common -y

# git && versions
sudo apt-get install git gitg qgit subversion -y

# Tecplot
sudo apt-get install libstdc++5 -y

# tools
sudo apt-get install gzip gimp htop meld gnuplot gnuplot-x11 vlc okular ddd gmsh unzip -y

# Secure Sockets Layer toolkit - development files
sudo apt install libssl-dev
sudo apt-get install openvpn openssl openssh-client -y

# for FLEXI/PICLas
sudo apt-get install liblapack3 liblapack-dev zlib1g-dev exuberant-ctags -y

# for documentation
sudo apt-get install texlive-base -y
sudo apt-get install texlive-latex-extra -y

# hdf5-file viewer
sudo apt-get install hdfview -y

# Install libs for reggie
sudo apt-get install python-h5py -y

# Further libs
# this is maybe not required (do not install them if it works without these packages)
# sudo apt-get install hdf5-tools libhdf5-dev -y

# Linux Standard Base (LSB): required on, e.g., Ubuntu Server that is equipped only thinly with pre-installed software
sudo apt-get install lsb -y
