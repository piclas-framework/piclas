#!/bin/bash -i

#==============================================================================
# title       : InstallAOTUS.sh
# description : This script installs AOTUS with specific setting in a
#               pre-installed module env
#
#               "The Aotus library provides a Fortran wrapper around the C-API
#               of the Lua scripting language, allowing a convenient usage of
#               Lua scripts as configuration files in Fortran applications."
#
# date        : Feb 9, 2023
# version     : 1.0
# usage       : bash InstallAOTUS.sh
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

# --------------------------------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------------------------------
# Check command line arguments
RERUNMODE=0
for ARG in "$@"
do

  if [ ${ARG} == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "no help available. Exit."
    echo "--help/-h            print help information"
    echo "--rerun/-r           Re-run mode: Try to continue an already failed installation attempt."
    exit
  fi

  # Check if re-run mode is selected by the user
  if [[ ${ARG} == "--rerun" ]] || [[ ${ARG} =~ ^-r(erun)?$ ]]; then
    RERUNMODE=1
  fi

done

# DOWNLOAD and INSTALL AOTUS
# f1d4a8743773   2022-11-18 21:24:50   Added support for the Fujitsu FRPTX compiler for ARM64FX
AOTUSVERSION=f1d4a8
FULLAOTUSVERSION=f1d4a87437734eaa40cc4df3f8fcb72df7771cec
#WEBSITE="https://osdn.net/projects/apes/scm/hg/aotus/tree/${FULLAOTUSVERSION}/"
WEBSITE="https://osdn.net/projects/apes/scm/hg/aotus/archive/${FULLAOTUSVERSION}/?format=tar.gz"

# --------------------------------------------------------------------------------------------------
# Check pre-requisites
# --------------------------------------------------------------------------------------------------
if [[ ! -x "$(command -v lua)" ]]; then
  read -p  "LUA not found. Install LUA via 'sudo apt install lua5.4'? Press [Enter] to continue installation or [Crtl+c] to abort!"
  sudo apt install lua5.4
else
  LUAVERSION=$(lua -v | cut -d " " -f2)
  echo "Detected lua version ${LUAVERSION}"
  LUAMAJOR=$(echo ${LUAVERSION} | cut -d "." -f1)
  LUAMINOR=$(echo ${LUAVERSION} | cut -d "." -f2)
  LUAPATCH=$(echo ${LUAVERSION} | cut -d "." -f3)
  if [[ ${LUAMAJOR} -lt 5 || ${LUAMINOR} -lt 4 ]]; then
    echo "Please install a newer LUA version >= 5.4"
    exit 1
  fi
fi

if [[ ! -x "$(command -v brave-browser)" ]]; then
  if [[ ! -x "$(command -v chromium)" ]]; then
    if [[ ! -x "$(command -v firefox)" ]]; then
      echo "Found neither brave-browser, nor chromium nor firefox. Exit."
      exit 1
    else
      BROWSER='firefox'
    fi
  else
    BROWSER='chromium'
  fi
else
  BROWSER='brave-browser'
fi

# --------------------------------------------------------------------------------------------------
#  Settings
# --------------------------------------------------------------------------------------------------

INSTALLDIR=/opt
SOURCESDIR=/opt/sources
MODULESDIR=/opt/modules/modulefiles
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/aotus/aotus_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi
# e.g. aotus-f1d4a8.tar.gz
TARFILE=${SOURCESDIR}/aotus-${AOTUSVERSION}.tar.gz

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p "${SOURCESDIR}"
fi

AOTUSMODULEFILEDIR=${MODULESDIR}/utilities/aotus
MODULEFILE=${AOTUSMODULEFILEDIR}/${AOTUSVERSION}

# if no AOTUS module for this compiler found, install AOTUS and create module
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "$GREEN""creating AOTUS-${AOTUSVERSION} under ${MODULEFILE}$NC"
  echo " "
  #echo -e "$GREEN""Important: If the compilation step fails, run the script again and if it still fails \n1) try compiling single, .i.e., remove -j from make -j or \n2) try make -j 2 (not all available threads)$NC"
  #echo " "
  echo -e "This will install AOTUS version ${GREEN}${AOTUSVERSION}${NC}.\n"
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  fi

  # Install destination
  AOTUSINSTALLDIR=/opt/aotus/${AOTUSVERSION}

  # Change to sources directors
  cd ${SOURCESDIR}

  # Check if tar file was already extracted
  CLONEDIR=${SOURCESDIR}/aotus-${AOTUSVERSION}
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

  # Download tar.gz file
  if [ ! -f ${TARFILE} ]; then
    echo ""
    echo "${YELLOW}${TARFILE} not found under ${SOURCESDIR}${NC}"
    echo "${YELLOW}Download the ${TARFILE} from ${WEBSITE} with ${BROWSER}.${NC}"
    echo "${YELLOW}After downloading, the file must be moved (with sudo rights) to ${SOURCESDIR} by hand, e.g.,  sudo mv ~/Desktop/${TARFILE} /opt/sources/${NC}"
    echo "${YELLOW}Afterwards, run this script again.${NC}"
    echo ""
    exit
  fi

  # Check if tar.gz file was correctly downloaded
  if [ ! -f ${TARFILE} ]; then
    echo -e "${RED}no source-file downloaded for aotus-${AOTUSVERSION}$NC"
    echo -e "${RED}check ${WEBSITE} again and try to download the .tar.gz file$NC"
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
      echo -e "$RED""Failed: [tar -xzf ${TARFILE} aotus-${AOTUSVERSION}]$NC"
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
  pwd
  ./waf configure --prefix="${AOTUSINSTALLDIR}" build

  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "${RED}Failed command: [./waf configure build].\nIf you get the error "/usr/bin/env: ‘python’: No such file or directory", either create a symbolic link from python to python 3 via\n\n    sudo ln -s /usr/bin/python3 /usr/bin/python\n\nor install the following package (Ubuntu)\n\n    sudo apt-get install python-is-python3\n$NC"
    exit
  else
    # Compile source files
    ./waf install
  fi

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed: [./waf install]$NC"
    exit
  fi

  # create modulefile if the installation seems successful (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -d "${AOTUSINSTALLDIR}/include" ]; then
    if [ ! -d "${AOTUSMODULEFILEDIR}" ]; then
      mkdir -p ${AOTUSMODULEFILEDIR}
    fi
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/aotusversion/'${AOTUSVERSION}'/gI' ${MODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${MODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${MODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${OPENMPIVERSION}'/gI' ${MODULEFILE}
    sed -i 's\AOTUSTOPDIR\'${AOTUSINSTALLDIR}'\gI' ${MODULEFILE}
  else
    echo -e "$RED""No module file created for AOTUS-${AOTUSVERSION}$NC"
    echo -e "$RED""no installation found in ${AOTUSINSTALLDIR}/include$NC"
  fi
else
  echo -e "$YELLOW""AOTUS-${AOTUSVERSION} already created: module file exists under ${MODULEFILE}$NC"
fi
