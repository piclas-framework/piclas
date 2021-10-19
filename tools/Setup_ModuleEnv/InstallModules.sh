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
MODULEVERSION='4.6.1'
#MODULEVERSION='5.0.0'

echo ""
echo -e "This will install Environment Modules version ${GREEN}${MODULEVERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
read -p "Press enter to continue!"

if [ "$MODULEVERSION" == "4.6.1" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-4.6.1/modules-4.6.1.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-4.6.1%2Fmodules-4.6.1.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959'
elif [ "$MODULEVERSION" == "5.0.0" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-5.0.0/modules-5.0.0.tar.gz?ts=gAAAAABhXF_Jxr8Tp_QSaLtNFIwqXte_JnuzMdO606UAbI0okB5uzbQFQ0B5NmlIQ-bgM2uJr_EAVFEtCm9GR0NuyWLthCNrrQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-5.0.0%2Fmodules-5.0.0.tar.gz%2Fdownload'
fi
calcTrue() { awk 'BEGIN{printf "%d\n" , ('"$*"'?1:0)}';}

# Create sources directory and copy the module templates
if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi
sudo cp -r moduletemplates /opt/sources/moduletemplates

# download and install modules framework if no modules are present
# Check ${MODULESHOME} variable which is set by module env if already installed
if [ ! -d "${MODULESHOME}" ]; then
  if [ ! -d "${INSTALLDIR}/modules/${MODULEVERSION}" ]; then
    echo "creating Module environment with modules-${MODULEVERSION}"
    cd ${SOURCESDIR}
    if [ ! -e "${SOURCESDIR}/modules-${MODULEVERSION}.tar.gz" ]; then
      #echo "--output-document=modules-${MODULEVERSION}.tar.gz ${MODULEDLINK}"
      #exit
      wget --output-document=modules-${MODULEVERSION}.tar.gz "${MODULEDLINK}"
    fi
    tar -xzf modules-${MODULEVERSION}.tar.gz && rm -rf modules-${MODULEVERSION}.tar.gz
    cd ${SOURCESDIR}/modules-${MODULEVERSION}

    if [ `calcTrue "$(echo 'puts $tcl_version;exit 0' | tclsh) < 8.5"` -eq 1 ]; then
      ./configure --prefix=${INSTALLDIR}/modules/${MODULEVERSION} --modulefilesdir=${INSTALLDIR}/modules/modulefiles --enable-dotmodulespath
    else
      CPPFLAGS="-DUSE_INTERP_ERRORLINE" ./configure --prefix=${INSTALLDIR}/modules/${MODULEVERSION} --modulefilesdir=${INSTALLDIR}/modules/modulefiles --enable-dotmodulespath
    fi
    make -j${NBROFCORES} 2>&1 | tee make.out
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " "
      echo "${RED}Failed: [make 2>&1 | tee make.out]${NC}"
      exit
    else
      make install 2>&1 | tee install.out
    fi

    # Copy initialization to /etc/profile
    #   /etc/profile: system-wide .profile file for the Bourne shell (sh(1))
    #   and Bourne compatible shells (bash(1), ksh(1), ash(1), ...).
    if [ -z "$(grep "if.*Modules.*${MODULEVERSION}.*init.*bash.*then" /etc/profile)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/init/bash ]; then' /etc/profile
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/init/bash' /etc/profile
      sed -i '3 i\fi' /etc/profile
    else
      echo "modules init already exists in /etc/profile"
      exit
    fi

    # Copy initialization to /etc/bash.bashrc
    #   System-wide .bashrc file for interactive bash(1) shells.
    if [ -z "$(grep "if.*Modules.*${MODULEVERSION}.*init.*bash.*then" /etc/bash.bashrc)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/init/bash ]; then' /etc/bash.bashrc
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/init/bash' /etc/bash.bashrc
      sed -i '3 i\fi' /etc/bash.bashrc
    else
      echo "modules init already exists in /etc/bash.bashrc"
      exit
    fi
    source /etc/profile

    # Change modulefiles path in init -> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    # comment everything in .modulespath
    if [ -f "${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath" ]; then
      sed -i 's/^/\# /' ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    fi
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/utilities" >> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/MPI" >> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    echo "/opt/modules/modulefiles/libraries" >> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    # echo "/home/\$\{USER\}/modulefiles" >> ${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
    mkdir -p ${INSTALLDIR}/modules/modulefiles/
    mkdir -p ${INSTALLDIR}/modules/modulefiles/compilers
    mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities
    mkdir -p ${INSTALLDIR}/modules/modulefiles/MPI
    mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries
    # cd /opt/modules/modulefiles/compilers/gcc
    # Check if .modulespath and bash exist
    if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath" ]; then
      if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/init/bash" ]; then
        echo "${GREEN}Modules correctly installed. System restart required.${NC}"
      else
        echo "${RED}bash was not created correctly.${NC}"
      fi
    else
      echo "${RED}.modulespath was not created correctly.${NC}"
    fi
  else
    echo "${YELLOW}module environment (modules-${MODULEVERSION}) already created${NC}"
  fi
else
  echo "${YELLOW}Module environment ($(module --version)) already existent${NC}"
  # Change modulefiles path in init -> ${MODULESHOME}/init/.modulespath
  # comment everything
  # Check if any line contains "/opt/modules/modulefiles/" in ${MODULESHOME}/init/.modulespath
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
  echo "${GREEN}Modules installed. System restart required.${NC}"
fi
