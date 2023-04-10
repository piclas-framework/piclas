#!/bin/bash

#==============================================================================
# title       : InstallPackages.sh
# description : This script installs the software packages required for
#               the module env script ParaView
# date        : Oct 21, 2021
# version     : 1.0
# usage       : bash InstallPackagesParaView.sh.sh
# notes       :
#==============================================================================

# Get OS Revision
LSBREV=$(lsb_release -sr)
LSBREVMAJOR=$(echo ${LSBREV%.*})
LSBREVMINOR=$(echo ${LSBREV##*.})
# Get OS Distributor ID
LSBNAME=$(lsb_release -si)

# Check OS
if [[ ${LSBNAME} != "Ubuntu" ]]; then
  echo "This script currently only supports Ubuntu. Exit."
  exit
fi


# For older versions
if [[ ${LSBREVMAJOR} -lt 21 ]]; then
  # from reggie (ubuntu 18)
  sudo apt-get install libpython-dev libphonon-dev libphonon4 -y
fi

sudo apt-get install libboost-dev -y
sudo apt-get install libxt-dev -y
sudo apt-get install mesa-common-dev -y

# Select revision-dependent packages (always assuming the contemporary ParaView version is used)
if [[ ${LSBREVMAJOR} == "16" ]]; then # Ubuntu 16
  #sudo apt-get install  qt4-qmake  qt4-dev-tools  libqt4-opengl-dev  libqt4-dev -y
  sudo apt-get install qt4-default qt4-dev-tools libqt4-dev qt4-qmake libqt4-opengl-dev -y
elif [[ ${LSBREVMAJOR} == "18" ]]; then # Ubuntu 18
  sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libgl1-mesa-dev -y
elif [[ ${LSBREVMAJOR} == "21" ]]; then # Ubuntu 21
  # ------------------------------------------------
  # This has been tested for ParaView version 5.9.1
  # ------------------------------------------------
  # qt5-default
  sudo apt install qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools -y
  # extras and mesa
  sudo apt install qttools5-dev -y
  sudo apt install libqt5x11extras5-dev -y
  sudo apt install libgl1-mesa-dev -y

  sudo apt install qtdeclarative5-dev -y

  # Could NOT find Qt5Svg (missing: Qt5Svg_DIR)
  sudo apt install libqt5svg5-dev -y

  # '/usr/lib/qt5/bin/xmlpatterns': No such file or directory
  # Tried the follwing, which did not help
  sudo apt install libqt5xmlpatterns5-dev -y
  # /usr/lib/qt5/bin/xmlpatterns is contained in package qtxmlpatterns5-dev-tools:
  sudo apt install qtxmlpatterns5-dev-tools -y

elif [[ ${LSBREVMAJOR} == "22" ]]; then # Ubuntu 22
  # ------------------------------------------------
  # This has been tested for ParaView version 5.9.1
  # ------------------------------------------------
  # python3 development
  sudo apt install python3-dev -y

  # extras and mesa
  sudo apt install libgl1-mesa-dev -y

  # Could NOT find Qt5Svg (missing: Qt5Svg_DIR)
  sudo apt install libqt5svg5-dev -y

  # '/usr/lib/qt5/bin/xmlpatterns': No such file or directory
  # Tried the follwing, which did not help
  sudo apt install libqt5xmlpatterns5-dev -y
  # /usr/lib/qt5/bin/xmlpatterns is contained in package qtxmlpatterns5-dev-tools:
  sudo apt install qtxmlpatterns5-dev-tools -y

fi
