#!/bin/bash -i

#==============================================================================
# title       : InstallHOPR.sh
# description : This script installs HOPR with specific setting in a
#               pre-installed module env
# date        : Apr 18, 2022
# version     : 1.0
# usage       : bash InstallHOPR.sh
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
LOADMODULES=1
# default to openmpi
WHICHMPI='openmpi'
for ARG in "$@"
do

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
    #CMAKEVERSION=3.15.3-d
    #CMAKEVERSION=3.17.0-d
    #CMAKEVERSION=3.20.3
    #CMAKEVERSION=3.21.3
    #CMAKEVERSION=3.24.2
    CMAKEVERSION=3.26.4

    #GCCVERSION=9.2.0
    #GCCVERSION=9.3.0
    #GCCVERSION=10.1.0
    #GCCVERSION=10.2.0
    #GCCVERSION=11.2.0
    #GCCVERSION=13.1.0
    GCCVERSION=13.2.0

    #OPENMPIVERSION=3.1.4
    #OPENMPIVERSION=4.0.1
    #OPENMPIVERSION=4.0.2
    #OPENMPIVERSION=3.1.6
    #OPENMPIVERSION=4.1.1
    #OPENMPIVERSION=4.1.5

    MPICHVERSION=4.1.2

    # chose which mpi you want to have installed (openmpi or mpich), default is openmpi
    if [[ -n ${MPICHVERSION} ]]; then
      # Set mpich or mpich-debug
      # MPICH "debug", which uses MPICH installation with --with-device=ch3:sock.
      # This will use the older ch3:sock channel that does not busy poll.
      # This channel will be slower for intra-node communication, but it will perform much better in the oversubscription scenario.
      WHICHMPI='mpich'
      #WHICHMPI='mpich-debug'
      MPIVERSION=${MPICHVERSION}
    else
      if [[ -z ${OPENMPIVERSION} ]]; then
        echo "${RED}ERROR: Set either OPENMPIVERSION or MPICHVERSION in InstallPETSc.sh when running with '-m'${NC}. Exit."
        exit
      else
        MPIVERSION=${OPENMPIVERSION}
      fi
    fi

    #HDF5VERSION=1.10.5
    #HDF5VERSION=1.10.6
    #HDF5VERSION=1.12.1
    HDF5VERSION=1.14.0
  fi

  # Check if re-run mode is selected by the user
  if [[ ${ARG} == "--rerun" ]] || [[ ${ARG} =~ ^-r(erun)?$ ]]; then
    RERUNMODE=1
  fi

done

HOPRVERSION='master'
HOPRDOWNLOAD='https://github.com/hopr-framework/hopr.git'

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
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/hopr/hopr_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p "${SOURCESDIR}"
fi

# Check if module is available (not required here, but for following libs)
if [[ -n $(module purge 2>&1) ]]; then
  echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
  exit
fi

# take the first gcc compiler installed with first compatible openmpi/mpich and hdf5
echo " "
if [[ $LOADMODULES -eq 1 ]]; then
  CMAKEVERSION=$(ls ${MODULESDIR}/utilities/cmake/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  MPIVERSION=$(ls ${MODULESDIR}/MPI/${WHICHMPI}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  HDF5VERSION=$(ls ${MODULESDIR}/libraries/hdf5/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  echo -e "Modules found automatically.\n\nCMAKEVERSION=${CMAKEVERSION}\nGCCVERSION=${GCCVERSION}\n${WHICHMPI}-MPIVERSION=${MPIVERSION}\nHDF5VERSION=${HDF5VERSION}\n\nWARNING: The combination might not be possible!"
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  fi
else
  echo "Modules defined by user. Check if the combination is possible!"
fi

check_module "cmake" "${CMAKEVERSION}"
check_module "gcc" "${GCCVERSION}"
check_module "${WHICHMPI}" "${MPIVERSION}"
check_module "hdf5" "${HDF5VERSION}"

HOPRMODULEFILEDIR=${MODULESDIR}/utilities/hopr/${HOPRVERSION}/gcc/${GCCVERSION}/${WHICHMPI}/${MPIVERSION}/hdf5
MODULEFILE=${HOPRMODULEFILEDIR}/${HDF5VERSION}

# if no HOPR module for this compiler found, install HOPR and create module
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "$GREEN""creating HOPR-${HOPRVERSION} for GCC-${GCCVERSION} under$NC"
  echo -e "$GREEN""$MODULEFILE$NC"
  echo " "
  module purge
  load_module "cmake/${CMAKEVERSION}"
  load_module "gcc/${GCCVERSION}"
  load_module "${WHICHMPI}/${MPIVERSION}/gcc/${GCCVERSION}"
  load_module "hdf5/${HDF5VERSION}/gcc/${GCCVERSION}/${WHICHMPI}/${MPIVERSION}"
  module list
  echo " "
  echo -e "$GREEN""Important: If the compilation step fails, run the script again and if it still fails \n1) try compiling single, .i.e., remove -j from make -j or \n2) try make -j 2 (not all available threads)$NC"
  echo " "
  echo -e "This will install HOPR version ${GREEN}${HOPRVERSION}${NC} from ${HOPRDOWNLOAD} \nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Have the correct modules been loaded? If yes, press [Enter] to continue or [Crtl+c] to abort! Note that hopr will be installed via cmake with -DLIBS_USE_MPI=OFF"
  fi

  # Install destination
  HOPRINSTALLDIR=/opt/hopr/${HOPRVERSION}/gcc-${GCCVERSION}/${WHICHMPI}-${MPIVERSION}/hdf5-${HDF5VERSION}

  # Create and change to install directory
  mkdir -p ${HOPRINSTALLDIR}
  cd ${HOPRINSTALLDIR}

  # Check if repo is already downloaded
  CLONEDIR=${HOPRINSTALLDIR}/hopr
  if [[ -d ${CLONEDIR} ]]; then
    # Inquiry: Continue the installation with the existing files OR remove them all and start fresh
    while true; do
      echo " "
      echo "${YELLOW}${CLONEDIR} already exists.${NC}"
      echo "${YELLOW}Do you want to continue the installation (y/n)?${NC}"
      # Inquiry
      if [[ ${RERUNMODE} -eq 0 ]]; then
        read -p "${YELLOW}Otherwise the directory will be removed and a fresh installation will be performed. [Y/n]${NC}" yn
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

  # Check if git repo exists/is available
  # echo "Checking git repository ${HOPRDOWNLOAD}"
  # wget ${HOPRDOWNLOAD} --no-check-certificate -o /dev/null
  # ERRORCODE=$?
  # if [ ${ERRORCODE} -ne 0 ]; then
  #   echo " "
  #   echo -e "$RED""Failed to find git repository: [wget ${HOPRDOWNLOAD} --no-check-certificate -o /dev/null]$NC"
  #   exit
  # fi

  # Clone git repository
  if [[ -d ${CLONEDIR} ]]; then
    cd ${CLONEDIR}
    git fetch
    git pull origin ${HOPRVERSION}
  else
    git clone ${HOPRDOWNLOAD}
    cd ${CLONEDIR}
  fi

  # Check if extraction failed
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed to download [git clone ${HOPRDOWNLOAD}] or update [git pull origin ${HOPRVERSION}]$NC"
    exit
  else
    HOPRBUILDDIR=${CLONEDIR}/build
    rm -rf ${HOPRBUILDDIR}
    mkdir -p ${HOPRBUILDDIR}
    cd ${HOPRBUILDDIR}
  fi

  # CMAKE COMPILE FLAGS
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DLIBS_USE_MPI=OFF \
        -DLIBS_BUILD_HDF5=OFF \
        ${CLONEDIR}

  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed command: [cmake -DCMAKE_BUILD_TYPE=Release \
                       -DLIBS_USE_MPI=ON \
                       -DLIBS_BUILD_HDF5=OFF \
                       ${CLONEDIR}]$NC"
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


  #ERRORCODE=$?
  #echo ${ERRORCODE}
  #echo `pwd`
  #ls -l
    #read -p "Download okay?"


  # create modulefile if the installation seems successful (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -e "${HOPRBUILDDIR}/bin/hopr" ]; then
    if [ ! -d "${HOPRMODULEFILEDIR}" ]; then
      mkdir -p ${HOPRMODULEFILEDIR}
    fi
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/hoprversion/'${HOPRVERSION}'/gI' ${MODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${MODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${MODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${MPIVERSION}'/gI' ${MODULEFILE}
    sed -i 's/HDF5VERSIONFLAG/'${HDF5VERSION}'/gI' ${MODULEFILE}
    sed -i 's\HOPRTOPDIR\'${HOPRBUILDDIR}'\gI' ${MODULEFILE}
    sed -i 's\HOPRWHICHMPI\'${WHICHMPI}'\gI' ${MODULEFILE}
  else
    echo -e "$RED""No module file created for HOPR-${HOPRVERSION} for GCC-${GCCVERSION}$NC"
    echo -e "$RED""no installation found in ${HOPRBUILDDIR}/bin$NC"
  fi
else
  echo -e "$YELLOW""HOPR-${HOPRVERSION} already created: module file exists under ${MODULEFILE}$NC"
fi
