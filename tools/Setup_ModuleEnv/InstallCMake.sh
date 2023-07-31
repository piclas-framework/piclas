#!/bin/bash -i

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

NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
INSTALLDIR=/opt
SOURCESDIR=/opt/sources
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/cmake/cmake_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL CMAKE (example cmake-3.4.3)
# For current releases, see: https://github.com/Kitware/CMake/releases/
#CMAKEVERSION='3.4.3'
#CMAKEVERSION='3.13.3'
#CMAKEVERSION='3.15.3'
#CMAKEVERSION='3.17.0'
#CMAKEVERSION='3.20.3'
#CMAKEVERSION='3.21.3'
#CMAKEVERSION='3.24.0'
#CMAKEVERSION='3.24.2'
CMAKEVERSION='3.26.4'

CMAKEDIR=${INSTALLDIR}/cmake/${CMAKEVERSION}/standard
MODULEFILE=${INSTALLDIR}/modules/modulefiles/utilities/cmake/${CMAKEVERSION}
BUILDDIR=${SOURCESDIR}/cmake-${CMAKEVERSION}/build
TARFILE=${SOURCESDIR}/cmake-${CMAKEVERSION}.tar.gz

# Remove INSTALL module directory during re-run
if [[ -n ${1} ]]; then
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
    #read -p "Delete ${MODULEFILE}?"
    rm ${MODULEFILE}
  fi
fi

if [ ! -e "${MODULEFILE}" ]; then

  # Check if module is available (not required here, but for following libs)
  if [[ -n $(module purge 2>&1) ]]; then
    echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read #! /bin/bash -i)${NC}"
    exit
  else
    #echo "MODULEPATH ="$MODULEPATH
    #echo "MODULESHOME="${MODULESHOME}
    #module av
    echo -e "${YELLOW}Warning: If the default gcc version is used for compiling CMake, the\nrequired GLIBCXX_3.4.XX (the default gcc11.2 compiler might have been wrongly labelled in ubuntu22 and actually 11.3 was installed) version might be\ntoo new and cmake will not work when another gcc version is loaded that does not support the required version.\n\nCheck with\n\n   strings /lib/i386-linux-gnu/libc.so.6 | grep GLIBC\n\nor with\n\n    ldd --version\n\nThe gcc compiler will now be loaded via 'module load gcc'.\n${NC}"
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
    module load gcc
    module li
  fi

  echo -e "This will install Cmake version ${GREEN}${CMAKEVERSION}${NC} with either the system gcc compiler or the module loaded version (note the warning above).\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
  read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  cd ${SOURCESDIR}

  # Download tar.gz file from FTP server
  if [ ! -f ${TARFILE} ]; then
    DOWNLOADPATH="https://github.com/Kitware/CMake/releases/download/v${CMAKEVERSION}/cmake-${CMAKEVERSION}.tar.gz"
    wget ${DOWNLOADPATH}
  fi

  # Check if tar.gz file was correctly downloaded, abort script if non-existent
  if [ ! -f ${TARFILE} ]; then
    echo "no cmake install-file downloaded for cmake-${CMAKEVERSION}"
    echo "check if ${DOWNLOADPATH} exists"
    exit
  fi

  # Extract tar.gz file
  tar -xzf cmake-${CMAKEVERSION}.tar.gz

  # Check if extraction failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " " && ls -l ${TARFILE}
    echo "${RED} Failed to extract: [tar -xzf ${TARFILE}]. Broken or failed download. Try removing ${TARFILE} before processing. Exit.${NC}"
    exit
  fi

  # Create build directory
  if [ ! -d ${BUILDDIR} ]; then
    mkdir -p ${BUILDDIR}
  fi
  # Remove SOURCE cmake-X.Y.Z/build/* directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] ; then
    #DELETE=$(echo ${SOURCESDIR}/cmake-${CMAKEVERSION}/build/*)
    #read -p "Delete ${DELETE} ?"
    rm ${BUILDDIR}/*
  fi

  # Change to build directory
  cd ${BUILDDIR}

  # Configure setup
  ../bootstrap --prefix=${CMAKEDIR}

  # Compile source files with NBROFCORES threads
  make -j${NBROFCORES} 2>&1 | tee make.out

  # Check if compilation failed
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
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/cmakeversion/'${CMAKEVERSION}'/g' ${MODULEFILE}
    sed -i 's/compilerversion/standard/g' ${MODULEFILE}

    # Remove SOURCE tar.gz file after successful installation
    if [[ -f ${TARFILE} ]]; then
      rm ${TARFILE}
    fi
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
  echo "${YELLOW}CMake-${CMAKEVERSION} already created (module file exists). Run with -r to remove and re-install.${NC}"
fi
