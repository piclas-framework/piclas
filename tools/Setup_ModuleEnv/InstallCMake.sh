#!/bin/bash

#==============================================================================
# title       : InstallCMake.sh
# description : This script installs cmake with a specified version as given
#               below via CMAKEVERSION='X.XX.X'
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallCMake.sh
# notes       :
#==============================================================================

# Check privilege
if [[ "$EUID" -ne 0 ]]; then
  echo "Please run as root"
  exit 1
fi

# --------------------------------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------------------------------

if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    TERMCOLS=$(tput cols)
    BOLD="$(tput bold)"
    UNDERLINE="$(tput smul)"
    STANDOUT="$(tput smso)"
    NORMAL="$(tput sgr0)"
    NC="$(tput sgr0)"
    BLACK="$(tput setaf 0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
    BLUE="$(tput setaf 4)"
    MAGENTA="$(tput setaf 5)"
    CYAN="$(tput setaf 6)"
    WHITE="$(tput setaf 7)"
  fi
fi

# --------------------------------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------------------------------

INSTALLDIR=/opt
SOURCESDIR=/opt/sources
TEMPLATEDIR=/opt/sources/moduletemplates

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL CMAKE (example cmake-3.4.3)
# For current releases, see: https://github.com/Kitware/CMake/releases/
#CMAKEVERSION='3.4.3'
#CMAKEVERSION='3.13.3'
#CMAKEVERSION='3.15.3'
CMAKEVERSION='3.17.0'
CMAKEDIR=${INSTALLDIR}/cmake/${CMAKEVERSION}/standard
MODULEFILE=${INSTALLDIR}/modules/modulefiles/utilities/cmake/${CMAKEVERSION}-d

if [[ -n ${1} ]]; then
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
    rm ${MODULEFILE}
  fi
fi

if [ ! -e "${MODULEFILE}" ]; then
  echo "${GREEN}creating CMake-${CMAKEVERSION}${NC}"
  cd ${SOURCESDIR}
  if [ ! -e "${SOURCESDIR}/cmake-${CMAKEVERSION}.tar.gz" ]; then
    wget "https://github.com/Kitware/CMake/releases/download/v${CMAKEVERSION}/cmake-${CMAKEVERSION}.tar.gz"
  fi
  tar -xzf cmake-${CMAKEVERSION}.tar.gz
  if [ ! -d "${SOURCESDIR}/cmake-${CMAKEVERSION}/build" ]; then
    mkdir -p ${SOURCESDIR}/cmake-${CMAKEVERSION}/build
  fi
  if [[ ${1} =~ ^-r(erun)?$ ]] ; then
    rm ${SOURCESDIR}/cmake-${CMAKEVERSION}/build/*
  fi
  cd ${SOURCESDIR}/cmake-${CMAKEVERSION}/build
  ../bootstrap --prefix=${CMAKEDIR}
  make -j 2>&1 | tee make.out
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo "${RED}Failed: [make -j 2>&1 | tee make.out]${NC}"
    exit 1
  else
    make install 2>&1 | tee install.out
  fi

  if [ -e "${CMAKEDIR}/bin/cmake" ] && [ -e "${CMAKEDIR}/bin/ccmake" ]; then
    if [ ! -e "${INSTALLDIR}/modules/modulefiles/utilities/cmake" ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities/cmake
    fi
    cp ${TEMPLATEDIR}/utilities/cmake/cmake_temp ${MODULEFILE}
    sed -i 's/cmakeversion/'${CMAKEVERSION}'/g' ${MODULEFILE}
    sed -i 's/compilerversion/standard/g' ${MODULEFILE}
    rm -rf cmake-${CMAKEVERSION}.tar.gz
  else
    echo "${RED}ERROR in cmake installation, no modulefile created${NC}"
    if [ ! -e ${CMAKEDIR}/bin/cmake ]; then
      echo "${RED}ERROR: cmake not installed${NC}"
    fi
    if [ ! -e ${CMAKEDIR}/bin/ccmake ]; then
      echo "${RED}ERROR: cmake-curses-gui not installed${NC}"
    fi
  fi
else
  echo "${YELLOW}CMake-${CMAKEVERSION} already created (module file exists)${NC}"
fi
