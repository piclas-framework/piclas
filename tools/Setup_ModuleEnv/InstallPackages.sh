#!/bin/bash
apt-get update
apt-get install \
  libgmp3-dev \
  libmpfr-dev \
  libmpfr-doc \
  libmpfr6 \
  libmpc-dev \
  build-essential \
  libc6-dev \
  binutils \
  libblas-dev \
  liblapack-dev \
  libpython-stdlib \
  libpython-dev \
  make \
  libphonon-dev \
  qt4-qmake \
  qt4-dev-tools \
  libqt4-opengl-dev \
  libqt4-dev \
  libxt-dev \
  mesa-common-dev \
  libz-dev \
  libncurses-dev \
  tcl \
  tcl8.6-dev \
  ctags \
  libboost-dev \
  libgl1-mesa-dri

mkdir -p /opt/Installsources
cp -r moduletemplates /opt/Installsources/moduletemplates
