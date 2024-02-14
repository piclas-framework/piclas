#!/bin/bash

#==============================================================================
# title       : InstallPackages.sh
# description : This script installs the software packages required for
#               the module env scripts for creating a software environment for
#               PICLas/FLEXI code frameworks
# date        : Oct 19, 2021
# version     : 1.0
# usage       : bash InstallPackagesServer.sh
# notes       :
#==============================================================================

# Get OS Distributor ID
LSBNAME=$(lsb_release -si)

# Check OS
if [[ ${LSBNAME} != "Ubuntu" ]]; then
  echo "This script currently only supports Ubuntu. Exit."
  exit
fi

# Linux Standard Base support package (LSB): required on, e.g., Ubuntu Server that is equipped only thinly with pre-installed software
sudo apt-get install  lsb -y
