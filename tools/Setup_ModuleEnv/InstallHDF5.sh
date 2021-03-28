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
    USECOMPILERVERSION=9.3.0
    USEMPIVERSION=3.1.6
    # Force --rerun via 'set'
    echo ""
    echo "Running '-m' with GCC $USECOMPILERVERSION and  OpenMPI $USEMPIVERSION"
    set -- -rerun
    break
  fi
done

INSTALLDIR=/opt
SOURCESDIR=/opt/sources
TEMPLATEDIR=/opt/sources/moduletemplates

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL HDF5 for every compiler openmpi / mpich combination (example HDF5-1.8.18)
#HDF5VERSION=1.8.18
#HDF5VERSION=1.10.4
HDF5VERSION=1.10.6

HDF5DIR=${INSTALLDIR}'/hdf5/'${HDF5VERSION}

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

    #--- build hdf5 in single
    MODULEFILE=${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/single
    echo "${GREEN}      Installing under: ${MODULEFILE}${NC}"
    if [[ -n ${1} ]]; then
      if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
        rm ${MODULEFILE}
      fi
    fi
    if [ ! -e "${MODULEFILE}" ]; then
      echo "${GREEN}    creating HDF5-${HDF5VERSION} library for ${WHICHCOMPILER}-${COMPILERVERSION} single${NC}"
      module purge
      if [[ -n $(module load ${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
        echo "${RED}      module ${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
        break
      fi
      module load ${WHICHCOMPILER}/${COMPILERVERSION}
      module list

      echo ""
      read -p "Have the correct modules been loaded? If yes, press enter to continue!"

      cd ${SOURCESDIR}
      if [ ! -e "${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz" ]; then
        wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz"
      fi
      if [ ! -e "${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz" ]; then
        echo "${RED} no hdf5 install-file downloaded for HDF5-${HDF5VERSION}${NC}"
        echo "${RED} check if https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz exists${NC}"
        break
      fi
      tar -xzf hdf5-${HDF5VERSION}.tar.gz hdf5-${HDF5VERSION}
      if [ ! -d "${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single" ]; then
        mkdir -p ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single
      fi
      if [[ ${1} =~ ^-r(erun)?$ ]] ; then
        rm ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single/*
      fi
      cd ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single/
      if [ "${WHICHCOMPILER}" == "gcc" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
      elif [ "${WHICHCOMPILER}" == "intel" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which icc) CXX=$(which icpc) FC=$(which ifort)
      fi
      make -j 2>&1 | tee make.out
      if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo " "
        echo "${RED} Failed: [make -j 2>&1 | tee make.out]${NC}"
        exit
      else
        make install 2>&1 | tee install.out
      fi

      cp ${TEMPLATEDIR}/libraries/hdf5/single_template ${MODULEFILE}
      sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
      sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
      sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE}
    else
      echo "${YELLOW}      HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} already created (module file exists)${NC}"
    fi

    #--- build hdf5 with mpi
    MPINAMES='openmpi mpich'
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
          if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
            rm ${MODULEFILE}
          fi
        fi
        if [ ! -e "${MODULEFILE}" ]; then
          module purge
          if [[ -n $(module load ${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
            echo "${RED}      module ${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
            break
          fi
          module load ${WHICHCOMPILER}/${COMPILERVERSION}
          if [[ -n $(module load ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION} 2>&1) ]]; then
            echo "${RED}      module ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION} not found ${NC}"
            break
          fi
          module load ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION}

          if [ ! -d "${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}" ]; then
            mkdir -p ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
          fi
          if [[ ${1} =~ ^-r(erun)?$ ]] ; then
            rm ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}/*
          fi
          cd ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
          ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION} --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel CC=$(which mpicc) CXX=$(which mpicxx) FC=$(which mpifort)
          make -j 2>&1 | tee make.out
          if [ ${PIPESTATUS[0]} -ne 0 ]; then
            echo " "
            echo "${RED}Failed: [make -j 2>&1 | tee make.out]${NC}"
            exit
          else
            make install 2>&1 | tee install.out
          fi
          cp ${TEMPLATEDIR}/libraries/hdf5/mpi_template ${MODULEFILE}
          sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
          sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
          sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE}
          sed -i 's/whichmpi/'${WHICHMPI}'/gI' ${MODULEFILE}
          sed -i 's/mpiversion/'${MPIVERSION}'/gI' ${MODULEFILE}
        else
          echo "${YELLOW}      HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} and ${WHICHMPI}-${MPIVERSION} already created (module file exists)${NC}"
          continue
        fi
      done
    done

    #rm -rf ${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz
  done

done
