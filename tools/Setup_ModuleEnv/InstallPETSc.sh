#!/bin/bash -i

#==============================================================================
# title       : InstallPETSc.sh
# description : This script installs PETSc with specific setting in a
#               pre-installed module env
# date        : Apr 19, 2022
# version     : 1.0
# usage       : bash InstallPETSc.sh
# notes       : Bash in run interactively via "-i" to use "module load/purge"
#               commands
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

  BLAS_SUPPORT=''
  if [ ${ARG} == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "Input arguments:"
    echo "--help/-h            print help information"
    echo "--modules/-m         use modules defined in this script by the user."
    echo "                     Otherwise, find modules automatically."
    echo "--blas/-b            Install with BLAS (FBLASLAPACK) by running configure command with --download-fblaslapack=1"
    exit
  elif [ ${ARG} == "--blas" ] || [ ${ARG} == "-b" ]; then
    echo "BLAS (FBLASLAPACK) will be installed via --download-fblaslapack=1"
    BLAS_SUPPORT=' --download-fblaslapack=1 '
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
    #GCCVERSION=12.2.0
    #GCCVERSION=13.1.0
    GCCVERSION=13.2.0

    #OPENMPIVERSION=3.1.4
    #OPENMPIVERSION=4.0.1
    #OPENMPIVERSION=4.0.2
    #OPENMPIVERSION=3.1.6
    #OPENMPIVERSION=4.1.1
    #OPENMPIVERSION=4.1.4
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

  fi

  # Check if re-run mode is selected by the user
  if [[ ${ARG} == "--rerun" ]] || [[ ${ARG} =~ ^-r(erun)?$ ]]; then
    RERUNMODE=1
  fi

done

# DOWNLOAD and INSTALL PETSc (example PETSc-3.17.0)
#PETSCVERSION=3.17.0
#PETSCVERSION=3.18.4
PETSCVERSION=3.19.3

# Activate DEBUGGING MODE with ON/OFF
DEBUG=OFF

if [[ ${DEBUG} == 'ON' ]]; then
  DEBUGDIR='_debug'
  WITHDEBUG='yes'
  TESTCOL=${YELLOW}
else
  DEBUGDIR=''
  WITHDEBUG='0'
  TESTCOL=${GREEN}
fi
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
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/petsc/petsc_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi
TARFILE=${SOURCESDIR}/petsc-${PETSCVERSION}.tar.gz

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p "${SOURCESDIR}"
fi

# Check if module is available (not required here, but for following libs)
if [[ -n $(module purge 2>&1) ]]; then
  echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
  exit
fi

# take the first gcc compiler installed with first compatible openmpi or mpich
echo " "
if [[ $LOADMODULES -eq 1 ]]; then
  CMAKEVERSION=$(ls ${MODULESDIR}/utilities/cmake/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  MPIVERSION=$(ls ${MODULESDIR}/MPI/${WHICHMPI}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  echo -e "Modules found automatically.\n\nCMAKEVERSION=${CMAKEVERSION}\nGCCVERSION=${GCCVERSION}\n${WHICHMPI}-MPIVERSION=${MPIVERSION}\n\nWARNING: The combination might not be possible!"
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  fi
else
  echo "Modules defined by user. Check if the combination is possible!"
fi

check_module "cmake" "${CMAKEVERSION}"
check_module "gcc" "${GCCVERSION}"
check_module "${WHICHMPI}" "${MPIVERSION}"

PETSCMODULEFILEDIR=${MODULESDIR}/utilities/petsc/${PETSCVERSION}${DEBUGDIR}/gcc/${GCCVERSION}/${WHICHMPI}
MODULEFILE=${PETSCMODULEFILEDIR}/${MPIVERSION}

# if no PETSc module for this compiler found, install PETSc and create module
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "$GREEN""creating PETSc-${PETSCVERSION} for GCC-${GCCVERSION} under$NC"
  echo -e "$GREEN""$MODULEFILE$NC"
  echo " "
  module purge
  load_module "cmake/${CMAKEVERSION}"
  load_module "gcc/${GCCVERSION}"
  load_module "${WHICHMPI}/${MPIVERSION}/gcc/${GCCVERSION}"
  module list
  echo " "
  echo -e "$GREEN""Important: If the compilation step fails, run the script again and if it still fails \n1) try compiling single, .i.e., remove -j from make -j or \n2) try make -j 2 (not all available threads)$NC"
  echo " "
  echo -e "This will install PETSc version ${GREEN}${PETSCVERSION}${NC} with ${TESTCOL}--with-debugging=${WITHDEBUG}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Have the correct modules been loaded? If yes, press [Enter] to continue or [Crtl+c] to abort!"
  fi

  # Install destination
  PETSCINSTALLDIR=/opt/petsc/${PETSCVERSION}${DEBUGDIR}/gcc-${GCCVERSION}/${WHICHMPI}-${MPIVERSION}

  # Change to sources directors
  cd ${SOURCESDIR}

  # Check if tar file was already extracted
  CLONEDIR=${SOURCESDIR}/petsc-${PETSCVERSION}
  if [[ -d ${CLONEDIR} ]]; then
    # Inquiry: Continue the installation with the existing files OR remove them all and start fresh
    while true; do
      echo " "
      echo "${YELLOW}${CLONEDIR} already exists.${NC}"
      echo "${YELLOW}Do you want to continue the installation with the existing files under ${CLONEDIR} (y/n)?${NC}"
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

  # Download tar.gz file
  if [ ! -f ${TARFILE} ]; then
    wget --output-document=${TARFILE} "https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSCVERSION}.tar.gz"
  fi

  # Check if tar.gz file was correctly downloaded
  if [ ! -f ${TARFILE} ]; then
    echo -e "$RED""no source-file downloaded for petsc-${PETSCVERSION}$NC"
    echo -e "$RED""check https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSCVERSION}.tar.gz$NC"
    exit
  fi

  # Extract tar.gz file
  if [[ ! -d ${CLONEDIR} ]]; then
    echo "Extracting ${TARFILE} to ${SOURCESDIR} ..."
    tar -xzf ${TARFILE}

    # Check if extraction failed
    ERRORCODE=$?
    if [ ${ERRORCODE} -ne 0 ]; then
      echo " "
      echo -e "$RED""Failed: [tar -xzf ${TARFILE} petsc-${PETSCVERSION}]$NC"
      exit
    fi
  fi

  # Check if decompressed directory exists
  if [[ -d ${CLONEDIR} ]]; then
    cd ${CLONEDIR}
  else
    echo -e "$RED""Failed: Cannot find extracted directory ${CLONEDIR}$NC"
    exit
  fi

  # Configure
  MPIINSTALLDIR=${INSTALLDIR}/${WHICHMPI}/${MPIVERSION}/gcc/${GCCVERSION}
  if [[ ! -d ${MPIINSTALLDIR} ]]; then
    echo -e "$RED""Failed: Cannot find MPI directory ${MPIINSTALLDIR}$NC"
    exit
  fi
  ./configure PETSC_ARCH=arch-linux \
	      --prefix=${PETSCINSTALLDIR} \
	      --with-mpi-dir=${MPIINSTALLDIR} \
	      --with-debugging=${WITHDEBUG} \
	      COPTFLAGS='-O3 -march=native -mtune=native' \
	      CXXOPTFLAGS='-O3 -march=native -mtune=native' \
	      FOPTFLAGS='-O3 -march=native -mtune=native' \
	      --download-hypre \
	      --download-mumps \
	      --download-scalapack \
        ${BLAS_SUPPORT}

  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed command: [./configure PETSC_ARCH=arch-linux \
	      --prefix=${PETSCINSTALLDIR} \
	      --with-mpi-dir=${MPIINSTALLDIR} \
	      --with-debugging=${WITHDEBUG} \
	      COPTFLAGS='-O3 -march=native -mtune=native' \
	      CXXOPTFLAGS='-O3 -march=native -mtune=native' \
	      FOPTFLAGS='-O3 -march=native -mtune=native' \
	      --download-hypre \
	      --download-mumps \
	      --download-scalapack \
        ${BLAS_SUPPORT}]$NC"
    exit
  else
    # Compile source files with NBROFCORES threads
    make -j${NBROFCORES} PETSC_DIR=${CLONEDIR} PETSC_ARCH=arch-linux all 2>&1 | tee make.out
  fi

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed: [make -j${NBROFCORES} 2>&1 | tee make.out]$NC"
    exit
  else
    # Install include/lib/share directories to target directory
    make PETSC_DIR=${CLONEDIR} PETSC_ARCH=arch-linux install 2>&1 | tee install.out
  fi

  # create modulefile if the installation seems successful (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -d "${PETSCINSTALLDIR}/include" ]; then
    if [ ! -d "${PETSCMODULEFILEDIR}" ]; then
      mkdir -p ${PETSCMODULEFILEDIR}
    fi
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/petscversion/'${PETSCVERSION}'/gI' ${MODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${MODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${MODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${MPIVERSION}'/gI' ${MODULEFILE}
    sed -i 's\PETSCTOPDIR\'${PETSCINSTALLDIR}'\gI' ${MODULEFILE}
    sed -i 's\PETSCWHICHMPI\'${WHICHMPI}'\gI' ${MODULEFILE}
  else
    echo -e "$RED""No module file created for PETSc-${PETSCVERSION} for GCC-${GCCVERSION}$NC"
    echo -e "$RED""no installation found in ${PETSCINSTALLDIR}/include$NC"
  fi
else
  echo -e "$YELLOW""PETSc-${PETSCVERSION} already created: module file exists under ${MODULEFILE}$NC"
fi
