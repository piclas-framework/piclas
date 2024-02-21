#!/bin/bash -i

#==============================================================================
# title       : InstallMPIallCOMPILERS.sh
# description : This script installs openmpi or mpich in a pre-installed module
#               env for all compiler that are found and able to be loaded
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallMPIallCOMPILERS.sh
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
# Settings
# --------------------------------------------------------------------------------------------------
NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
# chose which mpi you want to have installed (openmpi or mpich)
#WHICHMPI=openmpi
WHICHMPI=mpich
#WHICHMPI=mpich-debug
# choose for which compilers mpi is build (gcc or intel)
WHICHCOMPILER=gcc

INSTALLDIR=/opt
SOURCESDIR=/opt/sources
MODULESDIR=/opt/modules/modulefiles
TEMPLATEPATH=$(echo `pwd`/moduletemplates/MPI/template)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

CONFIGSUFFIX=''
if [ "${WHICHMPI}" == "openmpi" ]; then
  # DOWNLOAD and INSTALL OPENMPI (example OpenMPI-2.1.6)
  #MPIVERSION=2.1.6
  #MPIVERSION=3.1.3
  #MPIVERSION=3.1.4
  #MPIVERSION=3.1.6
  #MPIVERSION=4.0.1
  #MPIVERSION=4.0.5
  #MPIVERSION=4.1.1
  #MPIVERSION=4.1.4
  MPIVERSION=4.1.5
elif [ "${WHICHMPI}" == "mpich" ]; then
  # DOWNLOAD and INSTALL MPICH (example mpich-3.2.0)
  #MPIVERSION=3.2
  MPIVERSION=4.1.2
elif [ "${WHICHMPI}" == "mpich-debug" ]; then
  # DOWNLOAD and INSTALL MPICH (example mpich-3.2.0)
  #MPIVERSION=3.2
  MPIVERSION=4.1.2
  CONFIGSUFFIX='--with-device=ch3:sock'
else
  echo -e "${RED}ERROR: Setting is neither 'openmpi' nor 'mpich'${NC}"
  echo -e "${RED}ERROR: no mpi installed will be installed. Exit.${NC}"
  exit
fi

WHICHMPISHORT=$(echo $WHICHMPI | cut -d "-" -f1)
WHICHMPIDEBUG=$(echo $WHICHMPI | cut -d "-" -f2)
MPIINSTALLDIR=${INSTALLDIR}/${WHICHMPI}/${MPIVERSION}
COMPILERPREFIX=compilers/ # required for modules 5.0.0
COMPILERPREFIX=

# --------------------------------------------------------------------------------------------------
# Install Module for MPI
# --------------------------------------------------------------------------------------------------

if [ "${WHICHCOMPILER}" == "gcc" ] || [ "${WHICHCOMPILER}" == "intel" ]; then

  if [ ! -d "${SOURCESDIR}" ]; then
    mkdir -p ${SOURCESDIR}
  fi

  # check how many ${WHICHCOMPILER} compilers are installed
  NCOMPILERS=$(ls ${MODULESDIR}/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | wc -l)
  # loop over all found compilers
  for i in $(seq 1 ${NCOMPILERS}); do
    # 'i'th compiler installed
    COMPILERVERSION=$(ls ${MODULESDIR}/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n ${i} | tail -n 1)
    MPIMODULEFILEDIR=${MODULESDIR}/MPI/${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}
    MPIMODULEFILE=${MPIMODULEFILEDIR}/${COMPILERVERSION}
    BUILDDIR=${SOURCESDIR}/${WHICHMPI}-${MPIVERSION}/build_${WHICHCOMPILER}-${COMPILERVERSION}
    TARFILE=${SOURCESDIR}/${WHICHMPI}-${MPIVERSION}.tar.gz
    TARFILESHORT=${SOURCESDIR}/${WHICHMPISHORT}-${MPIVERSION}.tar.gz

    # Remove INSTALL module directory during re-run
    if [[ -n ${1} ]]; then
      if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MPIMODULEFILE} ]]; then
        rm ${MPIMODULEFILE}
      fi
    fi
    # if no mpi module for this compiler found, install ${WHICHMPI} and create module
    if [ ! -e "${MPIMODULEFILE}" ]; then

      # Check if module is available (not required here, but for following libs)
      if [[ -n $(module purge 2>&1) ]]; then
        echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
        exit
      fi

      echo -e "$GREEN""creating ${WHICHMPI}-${MPIVERSION} for ${WHICHCOMPILER}-${COMPILERVERSION}${NC}"

      # Unload all possibly loaded modules and load specific modules for compilation of MPI
      module purge
      module av
      if [[ -n $(module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
        echo -e "${RED}module ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
        break
      fi
      module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION}

      echo ""
      module list
      echo ""
      echo "Have the correct modules been loaded?"
      echo -e "This will install ${WHICHMPI} version ${GREEN}${MPIVERSION}${NC} for ${WHICHCOMPILER} version ${GREEN}${COMPILERVERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
      read -p "If yes, press [Enter] to continue or [Crtl+c] to abort or [n] to skip this combination and go to the next!" ANSWER
      if [[ $ANSWER == "n" ]]; then
        continue
      fi

      # Change to sources directories
      cd ${SOURCESDIR}

      # Download tar.gz file
      if [[ ! -f ${TARFILESHORT} ]]; then
        if [ "${WHICHMPISHORT}" == "openmpi" ]; then
          wget "https://www.open-mpi.org/software/ompi/v${MPIVERSION%.*}/downloads/openmpi-${MPIVERSION}.tar.gz"
        elif [ "${WHICHMPISHORT}" == "mpich" ]; then
          wget "http://www.mpich.org/static/downloads/${MPIVERSION}/mpich-${MPIVERSION}.tar.gz"
        fi
      fi

      # Check if tar.gz file was correctly downloaded
      if [[ ! -f ${TARFILESHORT} ]]; then
        if [ "${WHICHMPISHORT}" == "openmpi" ]; then
          echo -e "${RED}no mpi install-file downloaded for OpenMPI-${MPIVERSION}${NC}"
          echo -e "${RED}check if https://www.open-mpi.org/software/ompi/v${MPIVERSION%.*}/downloads/openmpi-${MPIVERSION}.tar.gz exists${NC}"
          break
        elif [ "${WHICHMPISHORT}" == "mpich" ]; then
          echo -e "${RED}no mpi install-file downloaded for MPICH-${MPIVERSION}${NC}"
          echo -e "${RED}check if http://www.mpich.org/static/downloads/${MPIVERSION}/mpich-${MPIVERSION}.tar.gz exists${NC}"
          break
        fi
      fi

      # Extract tar.gz file
      #tar -xzf ${WHICHMPISHORT}-${MPIVERSION}.tar.gz
      if [ "${WHICHMPIDEBUG}" == "debug" ]; then
        mv ${TARFILESHORT} ${TARFILE}
      fi
      tar -xzf ${TARFILE}

      # Check if extraction failed
      if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo " " && ls -l ${TARFILE}
        echo "${RED} Failed to extract: [tar -xzf ${TARFILE}]. Broken or failed download. Try removing ${TARFILE} before processing. Exit.${NC}"
        exit
      fi

      if [ ! -d ${BUILDDIR} ]; then
        mkdir -p ${BUILDDIR}
      fi

      # Remove SOURCE ${BUILDDIR}/* directory during re-run
      if [[ ${1} =~ ^-r(erun)?$ ]] ; then
        #DELETE=$(echo ${BUILDDIR}/*)
        #read -p "Delete ${DELETE} ?"
        rm -rf ${BUILDDIR}/*
      fi

      # Change to build directory
      cd ${BUILDDIR}

      # Configure setup for openmpi or mpich using the same command
      if [ "${WHICHCOMPILER}" == "gcc" ]; then
        ../configure ${CONFIGSUFFIX} --prefix=${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION} CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
      elif [ "${WHICHCOMPILER}" == "intel" ]; then
        ../configure ${CONFIGSUFFIX} --prefix=${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION} CC=$(which icc) CXX=$(which icpc) FC=$(which ifort)
      fi

      # Compile source files with NBROFCORES threads
      make -j${NBROFCORES} 2>&1 | tee make.out

      # Check if compilation failed
      if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo " "
        echo -e "${RED}Failed: [make -j 2>&1 | tee make.out]${NC}"
        break
      else
        make install 2>&1 | tee install.out
      fi

      # Create module file if installation seems to have been successful (check if mpicc, mpicxx, mpifort exists in installdir)
      if [ -e "${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION}/bin/mpicc" ] && [ -e "${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION}/bin/mpicxx" ] && [ -e "${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION}/bin/mpifort" ]; then
        if [ ! -d "${MPIMODULEFILEDIR}" ]; then
          mkdir -p ${MPIMODULEFILEDIR}
        fi
        cp ${TEMPLATEPATH} ${MPIMODULEFILE}
        sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MPIMODULEFILE}
        sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MPIMODULEFILE}
        sed -i 's/whichmpi/'${WHICHMPI}'/gI' ${MPIMODULEFILE}
        sed -i 's/mpiversion/'${MPIVERSION}'/gI' ${MPIMODULEFILE}

        # Remove SOURCE tar.gz file after successful installation
        if [[ -f ${TARFILE} ]]; then
          rm ${TARFILE}
        fi
        if [[ -f ${TARFILESHORT} ]]; then
          rm ${TARFILESHORT}
        fi
      else
        echo -e "${RED}No module file created for ${WHICHMPI}-${MPIVERSION} for ${WHICHCOMPILER}-${COMPILERVERSION}${NC}"
        echo -e "${RED}No mpi found in ${MPIINSTALLDIR}/${WHICHCOMPILER}/${COMPILERVERSION}/bin${NC}"
      fi

    else
      echo -e "${YELLOW}${WHICHMPI}-${MPIVERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} already created (module file exists). Run with -r to remove and re-install.${NC}"
      continue
    fi
  done
  #cd ${SOURCESDIR}
  #rm -rf ${WHICHMPI}-${MPIVERSION}.tar.gz
else
  echo "WHICHCOMPILER-flag neither 'gcc' nor 'intel'"
  echo "no mpi installed"
  exit
fi
