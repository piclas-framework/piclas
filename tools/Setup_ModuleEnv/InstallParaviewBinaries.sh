#!/bin/bash -i

#==============================================================================
# title       : InstallParaviewBindaries.sh
# description : This script installs paraview binaries (https://www.paraview.org/)
#               with a specified version defined by the download link below
# date        : Dec 04, 2024
# version     : 1.0
# usage       : sudo InstallParaviewBindaries.sh
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

INSTALLDIR=/opt
SOURCESDIR=/opt/sources
TEMPLATEPATH=$(echo `pwd`/moduletemplates/paraview/paraviewbinaries_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL Linaro (example ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64)
# For current releases, see: https://www.paraview.org/download/
#PARAVIEWVERSION='5.11.2-egl-binaries'

# ===== Adjust these three parameters when changing the version =====

PARAVIEWVERSION='5.11.2-osmesa-binaries'
DOWNLOADPATH="https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.11&type=binary&os=Linux&downloadFile=ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64.tar.gz"
PARAVIEWDIRNAME="ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64"

# ===================================================================
[ -z ${PARAVIEWVERSION} ] && echo "${RED}ERROR: Variable PARAVIEWVERSION cannot be empty${NC}" && exit 1
[ -z ${DOWNLOADPATH} ] && echo "${RED}ERROR: Variable DOWNLOADPATH cannot be empty${NC}" && exit 1
[ -z ${PARAVIEWDIRNAME} ] && echo "${RED}ERROR: Variable PARAVIEWDIRNAME cannot be empty${NC}" && exit 1

PARAVIEWINSTALLDIR=${INSTALLDIR}/paraview/${PARAVIEWVERSION}
MODULEFILE=${INSTALLDIR}/modules/modulefiles/paraview/paraview/${PARAVIEWVERSION}
BUILDDIR=${SOURCESDIR}/${PARAVIEWDIRNAME}
TARFILE="${SOURCESDIR}/${PARAVIEWDIRNAME}.tar.gz"

# Check if module file already exists under, e.g., /opt/modules/modulefiles/paraview/paraview/5.11.2-osmesa-binaries
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "This will install paraview binaries with version ${GREEN}${PARAVIEWVERSION}${NC} from ${GREEN}${DOWNLOADPATH}${NC}."
  echo -e "and place a module file under  ${GREEN}${MODULEFILE}${NC}"
  read -p "Press [Enter] to continue or [Crtl+c] to abort!"
fi

# Remove INSTALL module directory during re-run
if [[ -n ${1} ]]; then
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
    read -p "[Step 1 of 2: Remove module file] Are you sure that you want to delete ${MODULEFILE}?"
    if [ -n ${MODULEFILE} ]; then
      rm ${MODULEFILE}
    fi
  fi
fi

# Check if module file does NOT exist under, e.g., /opt/modules/modulefiles/paraview/paraview/5.11.2-osmesa-binaries
if [ ! -e "${MODULEFILE}" ]; then

  # Check if module is available (not required here, but for following libs)
  if [[ -n $(module purge 2>&1) ]]; then
    echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read #! /bin/bash -i)${NC}"
    exit
  fi

  # Switch to the sources directory
  cd ${SOURCESDIR}

  # Download tar.gz file from FTP server
  if [ ! -f ${TARFILE} ]; then
    wget --output-document=${TARFILE} ${DOWNLOADPATH}
  fi

  # Check if tar.gz file was correctly downloaded, abort script if non-existent
  if [ ! -f ${TARFILE} ]; then
    echo "no paraview install-file downloaded for ${PARAVIEWDIRNAME}"
    echo "check if ${DOWNLOADPATH} exists"
    exit
  fi

  # Extract tar.gz file
  echo "Extracting ${TARFILE} ..."
  tar -xzf ${TARFILE}

  # Check if extraction failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " " && ls -l ${TARFILE}
    echo "${RED} Failed to extract: [tar -xzf ${TARFILE}] (probably broken or failed download). Try removing ${TARFILE} before processing. Exit.${NC}"
    exit
  fi

  # Remove SOURCE ParaView-5.11.2-MPI-Linux-Python3.9-x86_64 directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] ; then
    read -p "[Step 2 of 2: Remove build directory/source files] Are you sure that you want to delete ${BUILDDIR} ?"
    if [ -n ${BUILDDIR} ]; then
      if [ -d ${BUILDDIR} ]; then
        # Delete the directory
        rm -r ${BUILDDIR}

        # Extract tar.gz file
        tar -xzf ${TARFILE}

        # Check if extraction failed
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
          echo " " && ls -l ${TARFILE}
          echo "${RED} Failed to extract: [tar -xzf ${TARFILE}] (probably broken or failed download). Try removing ${TARFILE} before processing. Exit.${NC}"
          exit
        fi
      fi
    fi
  fi

  # Install
  if [ -n ${BUILDDIR} ]; then
    # Check if the install directory exists
    if [ ! -d "${INSTALLDIR}/paraview" ]; then
      mkdir -p ${INSTALLDIR}/paraview
    fi
    if [ -d ${PARAVIEWINSTALLDIR} ]; then
      echo "${RED} The install directory [PARAVIEWINSTALLDIR=${PARAVIEWINSTALLDIR}] already exists. Exit.${NC}"
      exit
    fi
    if [ -d ${BUILDDIR} ]; then
      echo "mv ${BUILDDIR} ${PARAVIEWINSTALLDIR}"
      mv ${BUILDDIR} ${PARAVIEWINSTALLDIR}
    else
      echo "${RED} Failed to create [BUILDDIR=${BUILDDIR}]. Exit.${NC}"
      exit
    fi
  else
    echo "${RED} Failed to create [BUILDDIR=${BUILDDIR}]. Exit.${NC}"
    exit
  fi

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo "${RED}Failed: [mv ${BUILDDIR} ${PARAVIEWINSTALLDIR}]. Maybe the target directory already exists?${NC}"
    exit 1
  fi

  # Check if the binaries exist
  if [ -e "${PARAVIEWINSTALLDIR}/bin/pvserver" ] && [ -e "${PARAVIEWINSTALLDIR}/bin/pvbatch" ] ; then
    # Check if the module directory exists
    if [ ! -e "${INSTALLDIR}/modules/modulefiles/paraview/paraview" ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/paraview/paraview
    fi
    # Copy the module file
    cp ${TEMPLATEPATH} ${MODULEFILE}

    # Replace the temporary version in the module file by the current version
    sed -i 's/paraviewversion/'${PARAVIEWVERSION}'/g' ${MODULEFILE}

    # Remove SOURCE tar.gz file after successful installation
    # [ -f ${TARFILE} ] && rm ${TARFILE}
  else
    echo "${RED}ERROR in installation, no modulefile created${NC}"
    echo "checking  ${PARAVIEWINSTALLDIR}/bin/paraview"
    [ ! -e ${PARAVIEWINSTALLDIR}/bin/pvserver ] && echo "${RED}ERROR: pvserver not installed${NC}"
    [ ! -e ${PARAVIEWINSTALLDIR}/bin/pvbatch ] && echo "${RED}ERROR: pvbatch not installed${NC}"
    exit 1
  fi
else
  echo "${YELLOW}${PARAVIEWDIRNAME} already created (module file exists under ${MODULEFILE}). Run with -r to remove and re-install.${NC}"
fi

echo ""
echo "${GREEN}Successfully installed paraview.${NC}"