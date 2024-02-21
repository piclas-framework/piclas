#!/bin/bash -i

#==============================================================================
# title       : InstallHDF5.sh
# description : This script installs HDF5 either for all compilers/mpi versions
#               found in a module env or for a specific compiler-mpi-combination
#               (also via module env)
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallHDF5.sh [--help] [--rerun] [--modules]
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
# Check command line arguments
# --------------------------------------------------------------------------------------------------
LOADMODULES=1
for arg in "$@"
do
  if [ "$arg" == "--help" ] || [ "$arg" == "-h" ]
  then
    echo ""
    echo "Input arguments:"
    echo "--help/-h            print this help information"
    echo "--rerun/-r           remove existing module files and re-install HDF5"
    echo "--modules/-m         use modules defined in script by the user."
    echo "                     Otherwise, find modules automatically and install"
    echo "                     all combinations possible."
    exit
  fi
  if [ "$arg" == "--modules" ] || [ "$arg" == "-m" ]
  then
    LOADMODULES=0
    # Set desired versions

    # GCC
    #USECOMPILERVERSION=13.1.0
    USECOMPILERVERSION=13.2.0

    # OpenMPI
    #MPINAMES='openmpi'
    #USEMPIVERSION=4.1.5

    # MPICH
    MPINAMES='mpich'
    USEMPIVERSION=4.1.2

    # MPICH "debug", which uses MPICH installation with --with-device=ch3:sock.
    # This will use the older ch3:sock channel that does not busy poll.
    # This channel will be slower for intra-node communication, but it will perform much better in the oversubscription scenario.
    #MPINAMES='mpich-debug'
    #USEMPIVERSION=4.1.2

    # Force --rerun via 'set'
    echo ""
    echo "Running '-m' with GCC $USECOMPILERVERSION and $MPINAMES $USEMPIVERSION"
    set -- -rerun
    break
  fi
done

if [[ $LOADMODULES -eq 1 ]]; then
  MPINAMES='openmpi mpich mpich-debug'
fi

NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
INSTALLDIR=/opt
SOURCESDIR=/opt/sources
TEMPLATEPATH=$(echo `pwd`/moduletemplates/libraries/hdf5)
if [[ ! -d ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL HDF5 for every compiler openmpi / mpich combination (example HDF5-1.8.18)
#HDF5VERSION=1.8.18
#HDF5VERSION=1.10.4
#HDF5VERSION=1.10.6
#HDF5VERSION=1.12.0 # CAUTION NIG_PIC_maxwell_RK4/TWT_recordpoints fails for: 1) gcc/11.2.0   2) cmake/3.21.3   3) openmpi/4.1.1/gcc/11.2.0   4) hdf5/1.12.0/gcc/11.2.0/openmpi/4.1.1
#HDF5VERSION=1.12.1
#HDF5VERSION=1.12.2
HDF5VERSION=1.14.0

COMPILERPREFIX=compilers/ # required for modules 5.0.0
MPIPREFIX=MPI/ # required for modules 5.0.0
COMPILERPREFIX=
MPIPREFIX=

HDF5DIR=${INSTALLDIR}'/hdf5/'${HDF5VERSION}
TARFILE=${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz

# Change to sources directors
cd ${SOURCESDIR}

echo -e "Download HF5 version ${GREEN}${HDF5VERSION}${NC}."
read -p "Press [Enter] to continue or [Crtl+c] to abort!"

# Download tar.gz file
if [ ! -f ${TARFILE} ]; then
  wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz"
fi

# Check if tar.gz file was correctly downloaded
if [ ! -f ${TARFILE} ]; then
  echo "${RED} no hdf5 install-file downloaded for HDF5-${HDF5VERSION}${NC}"
  echo "${RED} check if https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz exists${NC}"
  break
fi

# Extract tar.gz file
tar -xzf ${TARFILE}

# Check if extraction failed
if [ ${PIPESTATUS[0]} -ne 0 ]; then
  echo " " && ls -l ${TARFILE}
  echo "${RED} Failed to extract: [tar -xzf ${TARFILE}]. Broken or failed download. Try removing ${TARFILE} before processing. Exit.${NC}"
  exit
fi

# Loop gcc and intel compilers
COMPILERNAMES='gcc intel'
for WHICHCOMPILER in ${COMPILERNAMES}; do
  echo "${GREEN}$WHICHCOMPILER ------------------------------------------------------------------------------${NC}"
  if [ ! -d "${INSTALLDIR}/modules/modulefiles/compilers/${WHICHCOMPILER}" ]; then
    break
  fi
  NCOMPILERS=$(ls ${INSTALLDIR}/modules/modulefiles/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | wc -l)
  if [[ $LOADMODULES -eq 0 ]]; then
    NCOMPILERS=1 # limit to one compiler version
  fi
  for i in $(seq 1 ${NCOMPILERS}); do
    COMPILERVERSION=$(ls ${INSTALLDIR}/modules/modulefiles/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n ${i} | tail -n 1)
    if [[ $LOADMODULES -eq 0 ]]; then
      COMPILERVERSION=$USECOMPILERVERSION
    fi
    echo "${GREEN}  $COMPILERVERSION ------------------------------------------------------------------------------${NC}"
    if [ ! -e "${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}" ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}
    fi

    BUILDDIR=${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}

    # Check if module is available (not required here, but for following libs)
    if [[ -n $(module purge 2>&1) ]]; then
      echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
      exit
    fi

    # ============================================================================================================================================================================
    #--- build hdf5 in single
    # ============================================================================================================================================================================
    MODULEFILE=${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/single
    echo "${GREEN}      Installing under: ${MODULEFILE}${NC}"
    if [[ -n ${1} ]]; then
      # Remove INSTALL module directory during re-run
      if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
        rm ${MODULEFILE}
      fi
    fi
    if [ ! -e "${MODULEFILE}" ]; then
      echo "${GREEN}    creating HDF5-${HDF5VERSION} library for ${WHICHCOMPILER}-${COMPILERVERSION} single${NC}"

      # Unload all possibly loaded modules and load specific modules for compilation of HDF5
      module purge
      if [[ -n $(module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
        echo "${RED}      module ${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
        break
      fi
      module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION}

      echo ""
      module list
      echo ""
      echo -e "Compiling HDF5 SINGLE mode.\nHave the correct modules been loaded?"
      echo -e "This will install HF5 version ${GREEN}${HDF5VERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
      read -p "Press [Enter] to continue or [Crtl+c] to abort!"

      if [ ! -d "${BUILDDIR}/single" ]; then
        mkdir -p ${BUILDDIR}/single
      fi

      # Remove SOURCE ${BUILDDIR}/single/* directory during re-run
      if [[ ${1} =~ ^-r(erun)?$ ]] ; then
        #DELETE=$(echo ${BUILDDIR}/single/*)
        #read -p "Delete ${DELETE} ?"
        rm -rf ${BUILDDIR}/single/*
      fi

      # Change to build directory
      cd ${BUILDDIR}/single/

      # Configure setup
      if [ "${WHICHCOMPILER}" == "gcc" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
      elif [ "${WHICHCOMPILER}" == "intel" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which icc) CXX=$(which icpc) FC=$(which ifort)
      fi

      # Compile source files with NBROFCORES threads
      make -j${NBROFCORES} 2>&1 | tee make.out

      # Check if compilation failed
      if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo " "
        echo "${RED} Failed: [make -j 2>&1 | tee make.out]${NC}"
        exit
      else
        make install 2>&1 | tee install.out
      fi

      # Create modulefile if installation seems successful
      cp ${TEMPLATEPATH}/single_template ${MODULEFILE}
      sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
      sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
      sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE}
    else
      echo "${YELLOW}      HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} already created (module file exists). Run with -r to remove and re-install.${NC}"
    fi
    # ============================================================================================================================================================================
    #--- build hdf5 in single
    # ============================================================================================================================================================================


    # ============================================================================================================================================================================
    #--- build hdf5 with mpi
    # ============================================================================================================================================================================
    for WHICHMPI in ${MPINAMES}; do
      echo "${GREEN}    $WHICHMPI ------------------------------------------------------------------------------${NC}"
      if [ ! -d "${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}" ]; then
        mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}
      fi
      if [ ! -d "${INSTALLDIR}/modules/modulefiles/MPI/${WHICHMPI}" ]; then
        break
      fi
      NMPI=$(ls ${INSTALLDIR}/modules/modulefiles/MPI/${WHICHMPI} | sed 's/ /\n/g' | grep -i "[0-9]\." | wc -l)
      if [[ $LOADMODULES -eq 0 ]]; then
        NMPI=1 # limit to one openmpi version
      fi
      for j in $(seq 1 ${NMPI}); do
        MPIVERSION=$(ls ${INSTALLDIR}/modules/modulefiles/MPI/${WHICHMPI}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n ${j} | tail -n 1)
        if [[ $LOADMODULES -eq 0 ]]; then
          MPIVERSION=$USEMPIVERSION
        fi
        echo "${GREEN}    $MPIVERSION ------------------------------------------------------------------------------${NC}"
        MODULEFILE=${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
        echo "${GREEN}      Installing under: ${MODULEFILE}${NC}"

        if [[ -n ${1} ]]; then
          # Remove INSTALL module directory during re-run
          if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
            rm ${MODULEFILE}
          fi
        fi
        if [ ! -e "${MODULEFILE}" ]; then

          # Unload all possibly loaded modules and load specific modules for compilation of HDF5
          module purge
          if [[ -n $(module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
            echo "${RED}      module ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
            break
          fi
          module load ${COMPILERPREFIX}${WHICHCOMPILER}/${COMPILERVERSION}
          if [[ -n $(module load ${MPIPREFIX}${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
            echo "${RED}      module ${MPIPREFIX}${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
            break
          fi
          module load ${MPIPREFIX}${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION}

          echo ""
          module list
          echo ""
          echo -e "Compiling HDF5 with MPI.\nHave the correct modules been loaded?"
          read -p "If yes, press [Enter] to continue or [Crtl+c] to abort!"

          if [ ! -d "${BUILDDIR}/${WHICHMPI}/${MPIVERSION}" ]; then
            mkdir -p ${BUILDDIR}/${WHICHMPI}/${MPIVERSION}
          fi

          # Remove SOURCE ${BUILDDIR}//${WHICHMPI}/${MPIVERSION}/* directory during re-run
          if [[ ${1} =~ ^-r(erun)?$ ]] ; then
            #DELETE=$(echo ${BUILDDIR}/${WHICHMPI}/${MPIVERSION}/*)
            #read -p "Delete ${DELETE} ?"
            rm ${BUILDDIR}/${WHICHMPI}/${MPIVERSION}/*
          fi

          # Change to build directory
          cd ${BUILDDIR}/${WHICHMPI}/${MPIVERSION}

          # Configure setup
          ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION} --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel CC=$(which mpicc) CXX=$(which mpicxx) FC=$(which mpifort)

          # Compile source files with NBROFCORES threads
          make -j 2>&1 | tee make.out
          if [ ${PIPESTATUS[0]} -ne 0 ]; then
            echo " "
            echo "${RED}Failed: [make -j 2>&1 | tee make.out]${NC}"
            exit
          else
            make install 2>&1 | tee install.out
          fi

          # Create modulefile if installation seems successful
          cp ${TEMPLATEPATH}/mpi_template ${MODULEFILE}
          sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
          sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
          sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE}
          sed -i 's/whichmpi/'${WHICHMPI}'/gI' ${MODULEFILE}
          sed -i 's/mpiversion/'${MPIVERSION}'/gI' ${MODULEFILE}
        else
          echo "${YELLOW}      HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} and ${WHICHMPI}-${MPIVERSION} already created (module file exists). Run with -r to remove and re-install.${NC}"
          continue
        fi
      done # j in $(seq 1 ${NMPI}); do
    done # WHICHMPI in ${MPINAMES}; do
    # ============================================================================================================================================================================
    #--- build hdf5 with mpi
    # ============================================================================================================================================================================


  done # i in $(seq 1 ${NCOMPILERS}); do
done # WHICHCOMPILER in ${COMPILERNAMES}; do

# Remove SOURCE tar.gz file after successful installation
if [[ -f ${TARFILE} ]]; then
  rm ${TARFILE}
fi
