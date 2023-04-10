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

calcTrue() { awk 'BEGIN{printf "%d\n" , ('"$*"'?1:0)}';}

# For current releases, see: https://sourceforge.net/projects/modules/files/Modules/
#MODULEVERSION='4.6.1'
#MODULEVERSION='5.0.0'
#MODULEVERSION='5.0.1'
MODULEVERSION='5.1.1'

NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
INSTALLDIR=/opt
SOURCESDIR=/opt/sources
INSTALLPREFIX=${INSTALLDIR}/modules/${MODULEVERSION}
INSTALLDIRMODULESFILES=${INSTALLDIR}/modules/modulefiles
MODULESPATH=${INSTALLDIR}/modules/${MODULEVERSION}/init/.modulespath
BUILDDIR=${SOURCESDIR}/modules-${MODULEVERSION}
TARFILE=${SOURCESDIR}/modules-${MODULEVERSION}.tar.gz

echo ""
echo -e "This will install Environment Modules version ${GREEN}${MODULEVERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
read -p "Press [Enter] to continue or [Crtl+c] to abort!"

#if [ "$MODULEVERSION" == "3.2.10" ]; then
  #MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-3.2.10/modules-3.2.10.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-3.2.10%2Fmodules-3.2.10.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959'
if [ "$MODULEVERSION" == "4.6.1" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-4.6.1/modules-4.6.1.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-4.6.1%2Fmodules-4.6.1.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959'
elif [ "$MODULEVERSION" == "5.0.0" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-5.0.0/modules-5.0.0.tar.gz?ts=gAAAAABhXF_Jxr8Tp_QSaLtNFIwqXte_JnuzMdO606UAbI0okB5uzbQFQ0B5NmlIQ-bgM2uJr_EAVFEtCm9GR0NuyWLthCNrrQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-5.0.0%2Fmodules-5.0.0.tar.gz%2Fdownload'
elif [ "$MODULEVERSION" == "5.0.1" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-5.0.1/modules-5.0.1.tar.gz?ts=gAAAAABhcU8vH1TG6SQext57ioXX7Ja-U8fAP4QzR_dzW9vk8_M4sH1kCEgrNjgUH1OYjtNI2bWUJLtIw3O0V3ClRMUhmqAoqw%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-5.0.1%2Fmodules-5.0.1.tar.gz%2Fdownload'
elif [ "$MODULEVERSION" == "5.1.1" ]; then
  MODULEDLINK='https://downloads.sourceforge.net/project/modules/Modules/modules-5.1.1/modules-5.1.1.tar.gz?ts=gAAAAABjLEKdzzvZbBrq1P1xVq0ZkSBd-sPTY61Zu6mneKWZNP111L-_WOIE-yEnfw_IVxvae25tU7tvcG4xoGEZIjp6Va7N6g%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2Flatest%2Fdownload'
fi

if [[ -z ${MODULEDLINK} ]]; then
  echo "ERROR: MODULEDLINK not set correctly. Exit."
  exit
fi


if [[ -n ${1} ]]; then
  # Remove INSTALL directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLDIRMODULESFILES} ]]; then
    #read -p "Delete ${INSTALLDIRMODULESFILES}?"
    rm -rf ${INSTALLDIRMODULESFILES}
  fi
fi

# Create sources directory and copy the module templates
if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi

# download and install modules framework if no modules are present
# Check if ${MODULESHOME} variable which is set by module env if already installed (possibly not loaded by this script?)
if [ ! -d "${MODULESHOME}" ]; then
  # Remove INSTALL directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLPREFIX} ]]; then
    #read -p "Delete ${INSTALLPREFIX}?"
    rm -rf ${INSTALLPREFIX}
  fi

  # Check if module environment (modules-${MODULEVERSION}) already created
  if [ ! -d ${INSTALLPREFIX} ]; then
    echo "creating Module environment with modules-${MODULEVERSION}"

    # Change to source dir
    cd ${SOURCESDIR}

    # Remove SOURCE modules-X.Y.Z directory during re-run
    if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d "${SOURCESDIR}/modules-${MODULEVERSION}" ]]; then
      #read -p "Delete ${SOURCESDIR}/modules-${MODULEVERSION}?"
      rm -rf "${SOURCESDIR}/modules-${MODULEVERSION}"
    fi

    # Download tar.gz file
    if [ ! -f ${TARFILE} ]; then
      wget --output-document=modules-${MODULEVERSION}.tar.gz "${MODULEDLINK}"
    fi

    # Check if tar.gz file was correctly downloaded, abort script if non-existent
    if [ ! -f ${TARFILE} ]; then
      echo "no environment modules install-file downloaded for modules-${MODULEVERSION}"
      echo "check if ${MODULEDLINK} exists"
      exit
    fi

    # Extract tar.gz file
    tar -xzf ${TARFILE}

    # Check if extraction failed
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " " && ls -l ${TARFILE}
      echo "${RED} Failed to extract: [tar -xzf ${TARFILE}]. Broken or failed download. Try removing ${TARFILE} before processing. Exit.${NC}"
      exit
    fi

#    # Create build directory
#    if [ ! -d ${BUILDDIR} ]; then
#      mkdir -p ${BUILDDIR}
#    fi
#
#    # Remove SOURCE cmake-X.Y.Z/build/* directory during re-run
#    if [[ ${1} =~ ^-r(erun)?$ ]] ; then
#      #DELETE=$(echo ${BUILDDIR}/*)
#      #read -p "Delete ${DELETE} ?"
#      rm -rf ${BUILDDIR}/*
#    fi

    # Change to build directory
    cd ${BUILDDIR}

    # Configure setup
    #
    #   --prefix=PREFIX
    #
    #       Installation root directory [/usr/local/Modules]
    #
    #   --with-initconf-in=VALUE
    #       Location where to install Modules initialization configuration files. Either initdir or etcdir (default=etcdir)
    #
    #       New in version 4.1.
    #       Changed in version 5.0: Configuration option default set to etcdir
    #
    #   --enable-modulepath
    #
    #       Set the module paths defined by --with-modulepath option in a .modulespath file (following C version fashion)
    #       within the initialization directory defined by the --initdir option rather than within the modulerc file. (default=no)
    #
    #   --initdir=DIR
    #
    #       Directory for the per-shell environment initialization scripts [PREFIX/init]
    #
    #   --with-modulepath=PATHLIST
    #       Default path list to setup as the default modulepaths. Each path in this list should be separated by ":".
    #       Defined value is registered in the initrc or modulespath configuration file, depending on the --enable-modulespath option.
    #       These files are respectively called modulerc and .modulespath if --with-initconf-in is set to initdir.
    #       The path list value is read at initialization time to populate the MODULEPATH environment variable.
    #       By default, this modulepath is composed of the directory set for the system modulefiles (default=PREFIX/modulefiles or
    #       BASEPREFIX/$MODULE_VERSION/modulefiles if versioning installation mode enabled)
    #
    #       New in version 4.0.
    #
    # registered in the initrc          ----> /opt/modules/5.0.0/etc/initrc
    # or modulespath configuration file ----> /opt/modules/5.0.0/etc/modulespath
    #
    #   --modulefilesdir=DIR
    #
    #       Directory of main modulefiles also called system modulefiles [PREFIX/modulefiles]
    #       S.C.: /op/modules/modulefiles
    #

    PATHLIST=/opt/modules/modulefiles/compilers:/opt/modules/modulefiles/utilities:/opt/modules/modulefiles/MPI:/opt/modules/modulefiles/libraries
    #PATHLIST=/opt/modules/modulefiles/compilersX

    # Check if TCL version is greater/equal 8.5
    TCLVERSION=$(echo 'puts $tcl_version;exit 0' | tclsh)
    RESULT=$(echo `calcTrue "$(echo 'puts $tcl_version;exit 0' | tclsh) < 8.5"`)
    if [[ -z ${RESULT} ]]; then
      echo "ERROR: could not retrieve [echo 'puts $tcl_version;exit 0' | tclsh]. Instead got [${RESULT}]. Exit."
      exit
    fi
    if [[ ${RESULT} -eq 0 ]]; then
      echo "TCL version $TCLVERSION is greater/equal 8.5: Setting CPPFLAGS=-DUSE_INTERP_ERRORLINE"
      export CPPFLAGS="-DUSE_INTERP_ERRORLINE"
    fi

    #./configure --prefix=${INSTALLPREFIX} --modulefilesdir=${INSTALLDIRMODULESFILES} --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --modulefilesdir=${INSTALLDIRMODULESFILES} --with-modulepath=${PATHLIST} --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --modulefilesdir=${INSTALLDIRMODULESFILES} --with-modulepath=${PATHLIST}
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST}
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --enable-modulespath --with-initconf-in=initdir
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --with-initconf-in=etcdir
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --with-initconf-in=initdir --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --with-initconf-in=initdir --enable-modulespath
    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST} --with-initconf-in=initdir
    #./configure --prefix=${INSTALLPREFIX} --enable-modulespath --with-modulepath=${PATHLIST} --with-initconf-in=initdir

    #./configure --prefix=${INSTALLPREFIX} --with-modulepath=${PATHLIST}
    # 5.0.0 creates initrc in etc/ without paths
    # 4.6.1 creates modulerc in init/ with paths

    #./configure --prefix=${INSTALLPREFIX} --enable-modulespath --with-modulepath=${PATHLIST}
    # 5.0.0 creates modulespath   initrc (paths are in initrc) in etc/
    # 4.6.1 creates .modulespath   modulerc (paths are in .modulespath and modulerc) in init/

    ./configure --prefix=${INSTALLPREFIX} --enable-modulespath --with-modulepath=${PATHLIST} --with-initconf-in=initdir
    # 5.0.1 creates .modulespath   modulerc (paths are in .modulespath) in init/
    # 5.0.0 creates .modulespath   modulerc (paths are in .modulespath) in init/
    # 4.6.1 creates .modulespath   modulerc (paths are in .modulespath) in init/


    # Check if configuration failed
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " "
      echo "${RED} Failed to configure: [${CPPFLAGS}]. Exit.${NC}"
      exit
    fi

    # Compile source files with NBROFCORES threads
    make -j${NBROFCORES} 2>&1 | tee make.out

    # Check if compilation failed
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
      echo " "
      echo "${RED}Failed: [make 2>&1 | tee make.out]${NC}"
      exit
    else
      make install 2>&1 | tee install.out
    fi

    # -----------------          other shells         -----------------
    echo ""
    echo "${GREEN}==============================================================================================================${NC}"
    echo "${GREEN}Initialization has been created in /etc/profile and /etc/profile but not for other shells.${NC}"
    echo "${GREEN}If you are using a different shell, e.g., zsh, then copy the following to your .zshrc.local or .zshrc to initialize the modules correctly${NC}"
    echo ""
    echo "if [ -f /opt/modules/'${MODULEVERSION}'/init/zsh ]; then"
    echo "  . /opt/modules/'${MODULEVERSION}'/init/zsh"
    echo "fi"
    echo ""
    echo "${GREEN}For other shells, note that 'zsh' must be exchanged for the specific shell type${NC}"
    echo "${GREEN}==============================================================================================================${NC}"
    echo ""
    # -----------------          other shells         -----------------


    # ----------------- THIS CHANGES /etc/profile -----------------
    # Copy initialization to /etc/profile
    #   /etc/profile: system-wide .profile file for the Bourne shell (sh(1))
    #   and Bourne compatible shells (bash(1), ksh(1), ash(1), ...).
    if [ -z "$(grep "if.*modules.*${MODULEVERSION}.*init.*bash.*then" /etc/profile)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/init/bash ]; then' /etc/profile
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/init/bash' /etc/profile
      sed -i '3 i\fi' /etc/profile
    else
      echo "modules init already exists in /etc/profile"
    fi
    # ----------------- THIS CHANGES /etc/profile -----------------


    # ----------------- THIS CHANGES /etc/bash.bashrc -----------------
    # Copy initialization to /etc/bash.bashrc
    #   System-wide .bashrc file for interactive bash(1) shells.
    if [ -z "$(grep "if.*modules.*${MODULEVERSION}.*init.*bash.*then" /etc/bash.bashrc)" ]; then
      sed -i '1 i\if [ -f /opt/modules/'${MODULEVERSION}'/init/bash ]; then' /etc/bash.bashrc
      sed -i '2 i\  . /opt/modules/'${MODULEVERSION}'/init/bash' /etc/bash.bashrc
      sed -i '3 i\fi' /etc/bash.bashrc
    else
      echo "modules init already exists in /etc/bash.bashrc"
    fi
    # ----------------- THIS CHANGES /etc/bash.bashrc -----------------

#    # Remove pre-installed modulefiles directory (also automatically installed stuff)
#    # if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -d ${INSTALLDIRMODULESFILES} ]]; then
#    #   rm -rf ${INSTALLDIRMODULESFILES}
#    #   read -p "Delete ${INSTALLDIRMODULESFILES}?"
#    # fi
#    # Change modulefiles path in init -> ${MODULESPATH}
#    if [ -f ${MODULESPATH} ]; then
#      # Comment every line in .modulespath
#      echo "${MODULESPATH} exists. Commenting out all lines in this file."
#      sed -i 's/^/\# /' ${MODULESPATH}
#    else
#      echo "${MODULESPATH} does not exist. Creating empty file for ${MODULESPATH}"
#      touch ${MODULESPATH}
#    fi
#    # add:
#    echo "/opt/modules/modulefiles/compilers" >> ${MODULESPATH}
#    echo "/opt/modules/modulefiles/utilities" >> ${MODULESPATH}
#    echo "/opt/modules/modulefiles/MPI" >> ${MODULESPATH}
#    echo "/opt/modules/modulefiles/libraries" >> ${MODULESPATH}
#
    mkdir -p ${INSTALLDIRMODULESFILES}
    mkdir -p ${INSTALLDIRMODULESFILES}/compilers
    mkdir -p ${INSTALLDIRMODULESFILES}/utilities
    mkdir -p ${INSTALLDIRMODULESFILES}/MPI
    mkdir -p ${INSTALLDIRMODULESFILES}/libraries

    source /etc/profile

    # Check if .modulespath and bash exist
    if [ -e ${MODULESPATH} ]; then
      if [ -e "${INSTALLDIR}/modules/${MODULEVERSION}/init/bash" ]; then
        echo "${GREEN}Modules correctly installed. System restart might be required.${NC}"

        # Remove SOURCE tar.gz file after successful installation
        if [[ -f ${TARFILE} ]]; then
          rm ${TARFILE}
        fi
      else
        echo "${RED}bash was not created correctly.${NC}"
      fi
    else
      echo "${RED}.modulespath was not created correctly in [${MODULESPATH}].${NC}"
    fi
  else
    # Module env not loaded, but folder already exists
    echo "${YELLOW}Module environment (modules-${MODULEVERSION}) already created (${INSTALLPREFIX} already exists). Reload your bash or re-install with -r!${NC}"
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

# Check if MODULESHOME is empty and INSTALLPREFIX is a directory
echo "MODULESHOME  :"${MODULESHOME}
echo "INSTALLPREFIX:"${INSTALLPREFIX}
if [ -z "${MODULESHOME}" ] && [ -d ${INSTALLPREFIX} ]; then
  echo "${GREEN}Modules installed. System restart might be required.${NC}"
fi
