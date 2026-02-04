#!/bin/bash -i

#==============================================================================
# title       : InstallLinaro.sh
# description : This script installs linaro (https://www.linaroforge.com/)
#               with a specified version as given below via
#               below via LINAROVERSION='XX.XX.XX'
# date        : Nov 12, 2024
# version     : 1.0
# usage       : bash InstallLinaro.sh
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
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/linaro/linaro_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL Linaro (example linaro-24.0.6)
# For current releases, see: https://www.linaroforge.com/download-documentation
LINAROVERSION='24.0.6'

LINAROINSTALLDIR=${INSTALLDIR}/linaro/forge/${LINAROVERSION}
MODULEFILE=${INSTALLDIR}/modules/modulefiles/utilities/linaro/${LINAROVERSION}
LINARODIRNAME=linaro-forge-${LINAROVERSION}-linux-x86_64
BUILDDIR=${SOURCESDIR}/${LINARODIRNAME}
TARFILE=${SOURCESDIR}/${LINARODIRNAME}.tar

# Check if module file already exists under, e.g., /opt/modules/modulefiles/utilities/linaro/24.0.6
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "This will install linaro version ${GREEN}${LINAROVERSION}${NC}."
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

# Check if module file does NOT exist under, e.g., /opt/modules/modulefiles/utilities/linaro/24.0.6
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
    DOWNLOADPATH="https://downloads.linaroforge.com/${LINAROVERSION}/$(basename ${TARFILE})"
    wget ${DOWNLOADPATH}
  fi

  # Check if tar.gz file was correctly downloaded, abort script if non-existent
  if [ ! -f ${TARFILE} ]; then
    echo "no linaro install-file downloaded for ${LINARODIRNAME}"
    echo "check if ${DOWNLOADPATH} exists"
    exit
  fi

  # Extract tar.gz file
  tar -xf ${TARFILE}

  # Check if extraction failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " " && ls -l ${TARFILE}
    echo "${RED} Failed to extract: [tar -xzf ${TARFILE}] (probably broken or failed download). Try removing ${TARFILE} before processing. Exit.${NC}"
    exit
  fi

  # Remove SOURCE linaro-X.Y.Z/build/* directory during re-run
  if [[ ${1} =~ ^-r(erun)?$ ]] ; then
    read -p "[Step 2 of 2: Remove build directory/source files] Are you sure that you want to delete ${BUILDDIR} ?"
    if [ -n ${BUILDDIR} ]; then
      if [ -d ${BUILDDIR} ]; then
        # Delete the linaro directory
        rm -r ${BUILDDIR}

        # Extract tar.gz file
        tar -xf ${TARFILE}

        # Check if extraction failed
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
          echo " " && ls -l ${TARFILE}
          echo "${RED} Failed to extract: [tar -xzf ${TARFILE}] (probably broken or failed download). Try removing ${TARFILE} before processing. Exit.${NC}"
          exit
        fi
      fi
    fi
  fi

  # Change to build directory
  cd ${BUILDDIR}

  # Install linaro
  ./textinstall.sh --accept-licence ${LINAROINSTALLDIR}

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo "${RED}Failed: [./textinstall.sh --accept-licence]${NC}"
    exit 1
  fi

  # Check if the binaries exist
  if [ -e "${LINAROINSTALLDIR}/bin/ddt" ] && [ -e "${LINAROINSTALLDIR}/bin/forge" ] && [ -e "${LINAROINSTALLDIR}/bin/map" ] && [ -e "${LINAROINSTALLDIR}/bin/perf-report" ]; then
    # Check if the module directory for linaro exists
    if [ ! -e "${INSTALLDIR}/modules/modulefiles/utilities/linaro" ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities/linaro
    fi
    # Copy the module file
    cp ${TEMPLATEPATH} ${MODULEFILE}

    # Replace the temporary linaro version in the module file by the current version
    sed -i 's/linaroversion/'${LINAROVERSION}'/g' ${MODULEFILE}

    # Remove SOURCE tar.gz file after successful installation
    # [ -f ${TARFILE} ] && rm ${TARFILE}
  else
    echo "${RED}ERROR in linaro installation, no modulefile created${NC}"
    echo "checking  ${LINAROINSTALLDIR}/bin/linaro"
    [ ! -e ${LINAROINSTALLDIR}/bin/ddt ] && echo "${RED}ERROR: ddt not installed${NC}"
    [ ! -e ${LINAROINSTALLDIR}/bin/forge ] && echo "${RED}ERROR: forge not installed${NC}"
    [ ! -e ${LINAROINSTALLDIR}/bin/map ] && echo "${RED}ERROR: map not installed${NC}"
    [ ! -e ${LINAROINSTALLDIR}/bin/perf-report ] && echo "${RED}ERROR: perf-report not installed${NC}"
    exit 1
  fi
else
  echo "${YELLOW}${LINARODIRNAME} already created (module file exists under ${MODULEFILE}). Run with -r to remove and re-install.${NC}"
fi

echo ""
echo "${GREEN}Successfully installed linaro. For information on the usage, see the following:${NC}"
echo ""
echo "${GREEN}    RELEASE-NOTES:  ${BUILDDIR}/RELEASE-NOTES${NC}"
echo "${GREEN}    Documentation:  https://docs.linaroforge.com/24.0.4/html/forge/forge/index.html${NC}"
echo ""
echo "${GREEN}For information on how to install a licence, see the following:${NC}"
echo ""
echo "${GREEN}    Documentation:  https://docs.linaroforge.com/24.0.4/html/forge/forge/licensing/index.html${NC}"
echo ""
echo "${GREEN}    If you have a licence (workstation or evaluation licence), copy it into the directory ${LINAROINSTALLDIR}/licences${NC}"
echo ""