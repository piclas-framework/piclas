#!/bin/bash -i

#==============================================================================
# title       : InstallOpenBLAS.sh
# description : This script installs OpenBLAS with specific setting in a
#               pre-installed module env
# date        : Feb 5, 2026
# version     : 1.0
# usage       : bash InstallOpenBLAS.sh
# notes       : Bash in run interactively via "-i" to use "module load/purge"
#               commands
#==============================================================================

# Check privilege
if [[ "$EUID" -ne 0 ]]; then
  echo "Please run as root"
  exit 1
fi

# Check git
if [[ ! -x "$(command -v git)" ]]; then
  echo "Please install git"
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
# Functions
# --------------------------------------------------------------------------------------------------
# Function for checking modules
check_module () {
                 if [ -z "${2}" ]; then
                   echo "module for ${1} not found. Exit"
                   exit
                 else
                   echo "${1}: [${2}]"
                 fi
                }

# Function for loading modules and checking if they exist in a specific combination
load_module () {
                 ERROR=$(module load "$1" 2>&1 | grep -in "ERROR")
                 if [ ! -z "$ERROR" ]; then
                   echo " "
                   echo -e "$RED$ERROR$NC"
                   echo -e "$RED""Failed: [module load $1]. Exit$NC"
                   exit
                 else
                   module load "$1"
                 fi
                }

# --------------------------------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------------------------------
# Check command line arguments
RERUNMODE=0
UPDATEMODE=0
LOADMODULES=1
# default to openmpi
for ARG in "$@"
do

  if [ ${ARG} == "--update" ] || [ ${ARG} == "-u" ]; then
    UPDATEMODE=1
  fi

  if [ ${ARG} == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "Input arguments:"
    echo "--help/-h            print help information"
    echo "--modules/-m         use modules defined in this script by the user."
    echo "                     Otherwise, find modules automatically."
    exit
  fi

  if [ ${ARG} == "--modules" ] || [ ${ARG} == "-m" ]; then
    LOADMODULES=0
    # Set desired versions
    GCCVERSION=14.2.0
    GCCVERSION=15.2.0

  fi

  # Check if re-run mode is selected by the user
  if [[ ${ARG} == "--rerun" ]] || [[ ${ARG} =~ ^-r(erun)?$ ]]; then
    RERUNMODE=1
  fi

done

OPENBLASVERSION=v0.3.31
OPENBLASDOWNLOAD='https://github.com/OpenMathLib/OpenBLAS.git'

# --------------------------------------------------------------------------------------------------
# Check pre-requisites
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
#  Settings
# --------------------------------------------------------------------------------------------------

NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
#NBROFCORES=1 # set to fixed value when errors are encountered at random (not related to missing packages)
INSTALLDIR=/opt
SOURCESDIR=/opt/sources
MODULESDIR=/opt/modules/modulefiles
TEMPLATEPATH=$(echo `pwd`/moduletemplates/openblas/openblas_temp)
[ ! -f ${TEMPLATEPATH} ] && echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit." && exit
[ ! -d "${SOURCESDIR}" ] && mkdir -p "${SOURCESDIR}"

# Check if module is available (not required here, but for following libs)
if [[ -n $(module purge 2>&1) ]]; then
  echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
  exit
fi

# take the first gcc compiler
echo " "
if [[ $LOADMODULES -eq 1 ]]; then
  GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  echo -e "Modules found automatically.\n\nGCCVERSION=${GCCVERSION}\n\nWARNING: The combination might not be possible!"
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  fi
else
  echo "Modules defined by user. Check if the combination is possible!"
fi

check_module "gcc" "${GCCVERSION}"

OPENBLASMODULEFILEDIR=${MODULESDIR}/openblas/openblas/${OPENBLASVERSION}/gcc
MODULEFILE=${OPENBLASMODULEFILEDIR}/${GCCVERSION}

# if no OPENBLAS module for this compiler found, install OPENBLAS and create module
if [[ ! -e "${MODULEFILE}" || ${UPDATEMODE} -eq 1 ]]; then
  # Check if module already exists and update mode is enabled (this will remove the existing module)
  if [[ -e "${MODULEFILE}" && ${UPDATEMODE} -eq 1 ]]; then
    echo " "
    echo -e "${YELLOW}Update mode detected (--update or -u): This will overwrite the installed module version with the latest one.${NC}"
    echo -e "${YELLOW}The 'OPENBLAS-${OPENBLASVERSION}' module file already exists under [${MODULEFILE}] and will be deleted now.${NC}"
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
    rm -rf ${MODULEFILE}
  fi
  echo -e "$GREEN""creating OPENBLAS-${OPENBLASVERSION} for GCC-${GCCVERSION} under$NC"
  echo -e "$GREEN""$MODULEFILE$NC"
  echo " "
  module purge
  load_module "gcc/${GCCVERSION}"
  module load cmake
  module list
  echo " "
  echo -e "$GREEN""Important: If the compilation step fails, run the script again and if it still fails \n1) try compiling single, .i.e., remove -j from make -j or \n2) try make -j 2 (not all available threads)$NC"
  echo " "
  echo -e "This will install OPENBLAS version ${GREEN}${OPENBLASVERSION}${NC} from ${OPENBLASDOWNLOAD} \nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Have the correct modules been loaded? If yes, press [Enter] to continue or [Crtl+c] to abort!"
  fi

  # Install destination
  OPENBLASINSTALLDIR=/opt/openblas/${OPENBLASVERSION}/gcc-${GCCVERSION}

  # Create and change to install directory
  # mkdir -p ${OPENBLASINSTALLDIR}
  cd ${SOURCESDIR}

  # Check if repo is already downloaded
  CLONEDIR=${SOURCESDIR}/openblas

  # If repo already exists and update is performed, remove the directory first.
  if [[ -d ${CLONEDIR} && ${UPDATEMODE} -eq 1 ]]; then
    echo " "
    echo -e "${YELLOW}Update mode detected (--update or -u): This will remove the existing git repository and clone a new one.${NC}"
    echo -e "${YELLOW}The cloned directory already exists under [${CLONEDIR}] and will be deleted now.${NC}"
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
    rm -rf ${CLONEDIR}
  fi

  # Check if directory exists
  if [[ -d ${CLONEDIR} ]]; then
    # Inquiry: Continue the installation with the existing files OR remove them all and start fresh
    while true; do
      echo " "
      echo "${YELLOW}${CLONEDIR} already exists.${NC}"
      echo "${YELLOW}Do you want to continue the installation? Otherwise the whole directory [${CLONEDIR}] will be removed and a fresh installation will be performed.${NC}"
      # Inquiry
      if [[ ${RERUNMODE} -eq 0 ]]; then
        read -p "${YELLOW}Continue the installation without removing? [Y/n]${NC}" yn
      else
        yn=y
      fi
      # Select case
      case $yn in
          [Yy]* ) echo "Continuing... "; break;;
          [Nn]* ) rm -rf ${CLONEDIR} ; break;;
          * ) echo "Please answer yes or no.";;
      esac
    done
  fi

  # Clone git repository
  if [[ -d ${CLONEDIR} ]]; then
    cd ${CLONEDIR}
    git fetch --all
    git checkout ${OPENBLASVERSION}
  else
    git clone ${OPENBLASDOWNLOAD} ${CLONEDIR}
    cd ${CLONEDIR}
  fi

  # Check if extraction failed
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed to download [git clone ${OPENBLASDOWNLOAD}] or update [git checkout ${OPENBLASVERSION}]$NC"
    exit
  else
    OPENBLASBUILDDIR=${CLONEDIR}/build-gcc-${GCCVERSION}
    rm -rf ${OPENBLASBUILDDIR}
    mkdir -p ${OPENBLASBUILDDIR}
    cd ${OPENBLASBUILDDIR}
  fi

  # CMAKE COMPILE FLAGS
  cmake -DCMAKE_INSTALL_PREFIX=${OPENBLASINSTALLDIR} ${CLONEDIR}

  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed command: [cmake -DCMAKE_INSTALL_PREFIX=${OPENBLASINSTALLDIR} ${CLONEDIR}]$NC"
    exit
  else
    # Compile source files with NBROFCORES threads
    make -j${NBROFCORES} 2>&1 | tee make.out
  fi

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed command: [make -j${NBROFCORES} 2>&1 | tee make.out]$NC"
    exit
  else
    make install 2>&1 | tee install.out
  fi

  # create modulefile if the installation seems successful (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -e "${OPENBLASINSTALLDIR}/lib/libopenblas.a" ]; then
    if [ ! -d "${OPENBLASMODULEFILEDIR}" ]; then
      mkdir -p ${OPENBLASMODULEFILEDIR}
    fi
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/openblasversion/'${OPENBLASVERSION}'/gI' ${MODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${MODULEFILE}
    sed -i 's\OPENBLASTOPDIR\'${OPENBLASINSTALLDIR}'\gI' ${MODULEFILE}
  else
    echo -e "${RED}ERROR: No module file created for OPENBLAS-${OPENBLASVERSION} for GCC-${GCCVERSION}${NC}"
    echo -e "${RED}ERROR: no installation found in ${OPENBLASBUILDDIR}/bin${NC}"
  fi
else
  echo -e "${YELLOW}WARNING: OPENBLAS-${OPENBLASVERSION} already created: module file exists under ${MODULEFILE}${NC}"
fi