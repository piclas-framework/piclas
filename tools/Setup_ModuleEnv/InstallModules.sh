#!/bin/bash

#==============================================================================
# title       : InstallModules.sh
# description : This script installs  either for all compilers/mpi versions 
# description : This script installs the module env software package with a 
#               specified version as given below via MODULEVERSION='X.X.XX' 
#               from https://downloads.sourceforge.net/
# date        : Nov 27, 2019
# version     : 1.0   
# usage       : bash InstallModules.sh
# notes       : 
#==============================================================================

INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
MODULEVERSION='3.2.10'
MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-3.2.10/modules-3.2.10.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-3.2.10%2Fmodules-3.2.10.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959'
calcTrue() { awk 'BEGIN{printf "%d\n" , ('"$*"'?1:0)}';}

if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi

# download and install modules framework if no modules are present
if [ ! -d "${MODULESHOME}" ]; then
  if [ ! -d "${INSTALLDIR}/modules/${MODULEVERSION}" ]; then
    echo "creating module evironment with modules-${MODULEVERSION}"
    cd ${SOURCESDIR}
    if [ ! -e "${SOURCESDIR}/modules-${MODULEVERSION}.tar.gz" ]; then
      #echo "--output-document=modules-${MODULEVERSION}.tar.gz ${MODULEDLINK}"
      #exit
      wget --output-document=modules-${MODULEVERSION}.tar.gz "${MODULEDLINK}"
    fi
    tar -xzf modules-${MODULEVERSION}.tar.gz && rm -rf modules-${MODULEVERSION}.tar.gz
    if [ ! -e "${SOURCESDIR}/modules-${MODULEVERSION}.tar.gz" ]; then
      mkdir -p ${SOURCESDIR}/modules-${MODULEVERSION}/build
    fi
    cd ${SOURCESDIR}/modules-${MODULEVERSION}/build

    if [ `calcTrue "$(echo 'puts $tcl_version;exit 0' | tclsh) < 8.5"` -eq 1 ]; then
      ../configure --prefix=${INSTALLDIR}/modules/${MODULEVERSION} --with-module-path=${INSTALLDIR}/modules/modulefiles
    else
      CPPFLAGS="-DUSE_INTERP_ERRORLINE" ../configure --prefix=${INSTALLDIR}/modules/${MODULEVERSION} --with-module-path=${INSTALLDIR}/modules/modulefiles
    fi
    make 2>&1 | tee make.out
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " "
      echo "Failed: [make 2>&1 | tee make.out]"
      exit
    else
      make install 2>&1 | tee install.out
    fi

    # Copy initialization to /etc/profile
    #   /etc/profile: system-wide .profile file for the Bourne shell (sh(1))
    #   and Bourne compatible shells (bash(1), ksh(1), ash(1), ...).
    if [ -z "$(grep "if.*Modules.*${MODULEVERSION}.*init.*bash.*then" /etc/profile)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/Modules/'${MODULEVERSION}'/init/bash ]; then' /etc/profile
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/Modules/'${MODULEVERSION}'/init/bash' /etc/profile
      sed -i '3 i\fi' /etc/profile
    else
      echo "modules init already exists in /etc/profile"
      exit
    fi

    # Copy initialization to /etc/bash.bashrc
    #   System-wide .bashrc file for interactive bash(1) shells.
    if [ -z "$(grep "if.*Modules.*${MODULEVERSION}.*init.*bash.*then" /etc/bash.bashrc)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/Modules/'${MODULEVERSION}'/init/bash ]; then' /etc/bash.bashrc
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/Modules/'${MODULEVERSION}'/init/bash' /etc/bash.bashrc
      sed -i '3 i\fi' /etc/bash.bashrc
    else
      echo "modules init already exists in /etc/bash.bashrc"
      exit
    fi
    source /etc/profile

    # Change modulefiles path in init -> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    # comment everything in .modulespath
    sed -i 's/^/\# /' ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/utilities" >> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/MPI" >> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/libraries" >> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    # echo "/home/\$\{USER\}/modulefiles" >> ${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath
    mkdir -p ${INSTALLDIR}/modules/modulefiles/
    mkdir -p ${INSTALLDIR}/modules/modulefiles/compilers
    mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities
    mkdir -p ${INSTALLDIR}/modules/modulefiles/MPI
    mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries
    # cd /opt/modules/modulefiles/compilers/gcc
    # Check if .modulespath and bash exist
    if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/.modulespath" ]; then
      if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/Modules/${MODULEVERSION}/init/bash" ]; then
        echo "Modules correctly installed. System restart required."
      else
        echo "bash was not created correctly."
      fi
    else
      echo ".modulespath was not created correctly."
    fi
  else
    echo "module evironment (modules-${MODULEVERSION}) already created"
  fi
else
  echo "module evironment ($(module --version)) already existant"
  # Change modulefiles path in init -> ${MODULESHOME}/init/.modulespath
  # comment everything
  if [ ! -z "$(grep "${INSTALLDIR}/modules/modulefiles/" ${MODULESHOME}/init/.modulespath)" ]; then
    sed -i 's/^/\# /' ${MODULESHOME}/init/.modulespath
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/utilities" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/MPI" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/libraries" >> ${MODULESHOME}/init/.modulespath
    # echo "/home/\$\{USER\}/modulefiles" >> ${MODULESHOME}/init/.modulespath
    mkdir -p ${INSTALLDIR}/modules/modulefiles/
    mkdir -p ${INSTALLDIR}/modules/modulefiles/compilers
    mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities
    mkdir -p ${INSTALLDIR}/modules/modulefiles/MPI
    mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries
    # cd /opt/modules/modulefiles/compilers/gcc
  fi
fi

if [ -z "${MODULESHOME}" ] && [ -d "${INSTALLDIR}/modules/${MODULEVERSION}" ]; then
  echo "Modules installed. System restart required."
fi
