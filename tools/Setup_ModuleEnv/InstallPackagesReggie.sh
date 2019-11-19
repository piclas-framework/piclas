#!/bin/bash
# --------------------------------------------------------------------------
# This script installs the software packages required for the regression 
# check tool "reggie2.0"
# --------------------------------------------------------------------------

# Check for updates
sudo apt-get update

# compiler
sudo apt-get install make cmake cmake-curses-gui gfortran g++ gcc 

# python
sudo apt-get install python python-numpy python3 python3-numpy python3-matplotlib python-matplotlib

# editors
sudo apt-get install gedit vim vim-runtime vim-common vim-tiny gdb vim-gui-common

# git && versions
sudo apt-get install git gitg qgit subversion

# paraview
sudo apt-get install libpython-dev libboost-dev libphonon-dev libphonon4 libxt-dev mesa-common-dev
#apt-get install qt4-default qt4-dev-tools libqt4-dev qt4-qmake libqt4-opengl-dev
sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libgl1-mesa-dev

# Tecplot
sudo apt-get install libstdc++5

# tools
sudo apt-get install gzip gimp htop meld gnuplot gnuplot-x11 vlc okular ddd gmsh unzip
sudo apt-get install openvpn openssl openssh-client

# for FLEXI/PICLas
sudo apt-get install liblapack3 liblapack-dev zlib1g-dev exuberant-ctags

# for documentation
sudo apt-get install pandoc pandoc-citeproc
sudo apt-get install texlive-full

# hdf5-file viewer
sudo apt-get install hdfview 

# Install libs for reggie
sudo apt-get install python-h5py

# Further libs
sudo apt-get install hdf5-tools libhdf5-dev # this is maybe not required (do not install them if it works without these packages)
