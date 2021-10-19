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
#MODULEVERSION='3.2.10'
MODULEVERSION='4.6.1'
#MODULEVERSION='5.0.0'
INSTALLPREFIX=${INSTALLDIR}/modules/${MODULEVERSION}
MODULESPATH=${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
INSTALLDIRMODULESFILES=${INSTALLDIR}/modules/modulefiles

echo ""
echo -e "This will install Environment Modules version ${GREEN}${MODULEVERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
read -p "Press enter to continue!"

if [ "$MODULEVERSION" == "4.6.1" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-4.6.1/modules-4.6.1.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-4.6.1%2Fmodules-4.6.1.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959'
elif [ "$MODULEVERSION" == "5.0.0" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-5.0.0/modules-5.0.0.tar.gz?ts=gAAAAABhXF_Jxr8Tp_QSaLtNFIwqXte_JnuzMdO606UAbI0okB5uzbQFQ0B5NmlIQ-bgM2uJr_EAVFEtCm9GR0NuyWLthCNrrQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-5.0.0%2Fmodules-5.0.0.tar.gz%2Fdownload'
fi
calcTrue() { awk 'BEGIN{printf "%d\n" , ('"$*"'?1:0)}';}

# Remove source directory during re-run
if [[ -n ${1} ]]; then
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d "${SOURCESDIR}/moduletemplates" ]]; then
    rm -rf "${SOURCESDIR}/moduletemplates"
    read -p "Delete ${SOURCESDIR}/moduletemplates?"
  fi
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLDIRMODULESFILES} ]]; then
    rm -rf ${INSTALLDIRMODULESFILES}
    read -p "Delete ${INSTALLDIRMODULESFILES}?"
  fi
fi

# Create sources directory and copy the module templates
if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi
sudo cp -r moduletemplates ${SOURCESDIR}/moduletemplates

# download and install modules framework if no modules are present
# Check if ${MODULESHOME} variable which is set by module env if already installed (possibly not loaded by this script?)
if [ ! -d "${MODULESHOME}" ]; then
  # Remove install directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLPREFIX} ]]; then
    rm -rf ${INSTALLPREFIX}
    read -p "Delete ${INSTALLPREFIX}?"
  fi

  # Check if module environment (modules-${MODULEVERSION}) already created
  if [ ! -d ${INSTALLPREFIX} ]; then
    echo "creating Module environment with modules-${MODULEVERSION}"

    # Change to source dir
    cd ${SOURCESDIR}

    # Remove modules-X.Y.Z directory during re-run
    if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d "${SOURCESDIR}/modules-${MODULEVERSION}" ]]; then
      rm -rf "${SOURCESDIR}/modules-${MODULEVERSION}"
      read -p "Delete ${SOURCESDIR}/modules-${MODULEVERSION}?"
    fi

    # Download tar.gz file
    if [ ! -e "${SOURCESDIR}/modules-${MODULEVERSION}.tar.gz" ]; then
      wget --output-document=modules-${MODULEVERSION}.tar.gz "${MODULEDLINK}"
    fi
    # Extract tar.gz file
    tar -xzf modules-${MODULEVERSION}.tar.gz && rm -rf modules-${MODULEVERSION}.tar.gz

    # Change directory
    cd ${SOURCESDIR}/modules-${MODULEVERSION}

    # Configure setup
    #
    #   --enable-dotmodulespath
    #
    #       Set the module paths defined by --with-modulepath option in a .modulespath file (following C version fashion)
    #       within the initialization directory defined by the --initdir option rather than within the modulerc file. (default=no)
    #
    if [ `calcTrue "$(echo 'puts $tcl_version;exit 0' | tclsh) < 8.5"` -eq 1 ]; then
      echo "Case A"
      ./configure --prefix=${INSTALLPREFIX} --modulefilesdir=${INSTALLDIR}/modules/modulefiles --enable-dotmodulespath
    else
      echo "Case B"
      CPPFLAGS="-DUSE_INTERP_ERRORLINE" ./configure --prefix=${INSTALLPREFIX} --modulefilesdir=${INSTALLDIR}/modules/modulefiles --enable-dotmodulespath
    fi

    make -j${NBROFCORES} 2>&1 | tee make.out

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " "
      echo "${RED}Failed: [make 2>&1 | tee make.out]${NC}"
      exit
    else
      make install 2>&1 | tee install.out
    fi

    # ----------------- THIS CHANGES /etc/profile -----------------
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
    # ----------------- THIS CHANGES /etc/profile -----------------


    # ----------------- THIS CHANGES /etc/bash.bashrc -----------------
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
    # ----------------- THIS CHANGES /etc/bash.bashrc -----------------


    source /etc/profile

    # Remove pre-installed modulefiles directory (also automatically installed stuff)
    # if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLDIRMODULESFILES} ]]; then
    #   rm -rf ${INSTALLDIRMODULESFILES}
    #   read -p "Delete ${INSTALLDIRMODULESFILES}?"
    # fi
    # Change modulefiles path in init -> ${MODULESPATH}
    if [ -f ${MODULESPATH} ]; then
      # Comment every line in .modulespath
      sed -i 's/^/\# /' ${MODULESPATH}
    else
      echo "${MODULESPATH} does not exit. Creating empty file for ${MODULESPATH}"
      touch ${MODULESPATH}
    fi
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/utilities" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/MPI" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/libraries" >> ${MODULESPATH}

    mkdir -p ${INSTALLDIRMODULESFILES}
    mkdir -p ${INSTALLDIRMODULESFILES}/compilers
    mkdir -p ${INSTALLDIRMODULESFILES}/utilities
    mkdir -p ${INSTALLDIRMODULESFILES}/MPI
    mkdir -p ${INSTALLDIRMODULESFILES}/libraries

    # Check if .modulespath and bash exist
    if [ -e ${MODULESPATH} ]; then
      if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/init/bash" ]; then
        echo "${GREEN}Modules correctly installed. System restart might be required.${NC}"
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
  MODULESPATH=${MODULESHOME}/init/.modulespath
  # Check if any line contains "/opt/modules/modulefiles/" in ${MODULESPATH}
  if [ ! -z "$(grep "${INSTALLDIR}/modules/modulefiles/" ${MODULESPATH})" ]; then
    # Comment every line in .modulespath
    sed -i 's/^/\# /' ${MODULESPATH}
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/utilities" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/MPI" >> ${MODULESPATH}
    echo "/opt/modules/modulefiles/libraries" >> ${MODULESPATH}

    mkdir -p ${INSTALLDIRMODULESFILES}
    mkdir -p ${INSTALLDIRMODULESFILES}/compilers
    mkdir -p ${INSTALLDIRMODULESFILES}/utilities
    mkdir -p ${INSTALLDIRMODULESFILES}/MPI
    mkdir -p ${INSTALLDIRMODULESFILES}/libraries
  fi
fi

if [ -z "${MODULESHOME}" ] && [ -d ${INSTALLPREFIX} ]; then
  echo "${GREEN}Modules installed. System restart might be required.${NC}"
fi
