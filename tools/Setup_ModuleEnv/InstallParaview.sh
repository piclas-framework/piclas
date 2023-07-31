#!/bin/bash -i

#==============================================================================
# title       : InstallParaview.sh
# description : This script installs paraview with specific setting in a
#               pre-installed module env
# date        : Nov 27, 2019
# version     : 1.0
# usage       : bash InstallParaview.sh
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
# Function for checking modules
check_module () {
                 if [ -z "${2}" ]; then
                   echo "module for ${1} not found. Exit"
                   exit
                 else
                   echo "${1}: [${2}]"
                 fi
                }

# Function for loading modules and checking if they exist in a specific combination
load_module () {
                 ERROR=$(module load "$1" 2>&1 | grep -in "ERROR")
                 if [ ! -z "$ERROR" ]; then
                   echo " "
                   echo -e "$RED$ERROR$NC"
                   echo -e "$RED""Failed: [module load $1]. Exit$NC"
                   exit
                 else
                   module load "$1"
                 fi
                }

# --------------------------------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------------------------------
# Check command line arguments
RERUNMODE=0
LOADMODULES=1
for ARG in "$@"
do

  if [ ${ARG} == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "Input arguments:"
    echo "--help/-h            print help information"
    echo "--modules/-m         use modules defined in this script by the user."
    echo "                     Otherwise, find modules automatically."
    exit
  fi

  if [ ${ARG} == "--modules" ] || [ ${ARG} == "-m" ]; then
    LOADMODULES=0
    # Set desired versions
    #CMAKEVERSION=3.15.3-d
    #CMAKEVERSION=3.17.0-d
    #CMAKEVERSION=3.20.3
    #CMAKEVERSION=3.21.3
    CMAKEVERSION=3.24.2

    #GCCVERSION=9.2.0
    #GCCVERSION=9.3.0
    #GCCVERSION=10.1.0
    #GCCVERSION=10.2.0
    #GCCVERSION=11.2.0
    GCCVERSION=12.2.0

    #OPENMPIVERSION=3.1.4
    #OPENMPIVERSION=4.0.1
    #OPENMPIVERSION=4.0.2
    #OPENMPIVERSION=3.1.6
    #OPENMPIVERSION=4.1.1
    OPENMPIVERSION=4.1.4

    #HDF5VERSION=1.10.5
    #HDF5VERSION=1.10.6
    #HDF5VERSION=1.12.1
    HDF5VERSION=1.12.2
  fi

  # Check if re-run mode is selected by the user
  if [[ ${ARG} == "--rerun" ]] || [[ ${ARG} =~ ^-r(erun)?$ ]]; then
    RERUNMODE=1
  fi

done

# DOWNLOAD and INSTALL PARAVIEW (example Paraview-5.0.0)
#PARAVIEWVERSION=5.2.0
#PARAVIEWVERSION=5.3.0

# Version 5.8.0 on Ubuntu 20.04 requires QT5, which is installed by default, but not with all required development libs
#  sudo apt-get install libqt5x11extras5-dev
#  sudo apt-get install libqt5svg5-dev
#  sudo apt-get install qtxmlpatterns5-dev-tools
#  sudo apt-get install qttools5-dev qt5-default libxt-dev libgl1-mesa-dev
#  sudo apt-get install python3.8-dev
#PARAVIEWVERSION=5.8.0

# Version 5.9.0 on Ubuntu 21.04 requires the following (also maybe some of the above, if they are not installed)
# sudo apt-get install libqt5x11extras5-dev
# sudo apt-get install qtdeclarative5-dev
# sudo apt-get install qttools5-dev
PARAVIEWVERSION=5.9.1

# Fix for error ‘numeric_limits’ is not a member of ‘std’ that occurs for 5.9.1 is done automatically below

# --------------------------------------------------------------------------------------------------
# Check pre-requisites
# --------------------------------------------------------------------------------------------------

if [[ ${PARAVIEWVERSION} == '5.9.0' ]] || [[ ${PARAVIEWVERSION} == '5.9.1' ]]; then
  echo -e "${GREEN}Installing libqt5x11extras5-dev   qtdeclarative5-dev    qttools5-dev    for this version (${PARAVIEWVERSION}) of ParaView${NC}"
  sudo apt-get install libqt5x11extras5-dev -y
  sudo apt-get install qtdeclarative5-dev -y
  sudo apt-get install qttools5-dev -y
fi

# --------------------------------------------------------------------------------------------------
#  Settings
# --------------------------------------------------------------------------------------------------

NBROFCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')
#NBROFCORES=1 # set to fixed value when errors are encountered at random (not related to missing packages)
INSTALLDIR=/opt
SOURCESDIR=/opt/sources
MODULESDIR=/opt/modules/modulefiles
TEMPLATEPATH=$(echo `pwd`/moduletemplates/utilities/paraview/paraview_temp)
if [[ ! -f ${TEMPLATEPATH} ]]; then
  echo "${RED}ERROR: module template not found under ${TEMPLATEPATH}${NC}. Exit."
  exit
fi
TARFILE=${SOURCESDIR}/paraview-${PARAVIEWVERSION}-source.tar.gz

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p "${SOURCESDIR}"
fi

# Check if module is available (not required here, but for following libs)
if [[ -n $(module purge 2>&1) ]]; then
  echo -e "${RED}module: command not found.\nThis script must be run in an interactive shell (the first line must read '#! /bin/bash' -i)${NC}"
  exit
fi

# take the first gcc compiler installed with first compatible openmpi and hdf5
echo " "
if [[ $LOADMODULES -eq 1 ]]; then
  CMAKEVERSION=$(ls ${MODULESDIR}/utilities/cmake/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  OPENMPIVERSION=$(ls ${MODULESDIR}/MPI/openmpi/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  HDF5VERSION=$(ls ${MODULESDIR}/libraries/hdf5/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  echo -e "Modules found automatically.\n\nCMAKEVERSION=${CMAKEVERSION}\nGCCVERSION=${GCCVERSION}\nOPENMPIVERSION=${OPENMPIVERSION}\nHDF5VERSION=${HDF5VERSION}\n\nWARNING: The combination might not be possible!"
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Press [Enter] to continue or [Crtl+c] to abort!"
  fi
else
  echo "Modules defined by user. Check if the combination is possible!"
fi

check_module "cmake" "${CMAKEVERSION}"
check_module "gcc  " "${GCCVERSION}"
check_module "mpi  " "${OPENMPIVERSION}"
check_module "hdf5 " "${HDF5VERSION}"

PARAVIEWMODULEFILEDIR=${MODULESDIR}/utilities/paraview/${PARAVIEWVERSION}/gcc/${GCCVERSION}/openmpi/${OPENMPIVERSION}/hdf5
MODULEFILE=${PARAVIEWMODULEFILEDIR}/${HDF5VERSION}

# if no paraview module for this compiler found, install paraview and create module
if [ ! -e "${MODULEFILE}" ]; then
  echo -e "$GREEN""creating Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION} under$NC"
  echo -e "$GREEN""$MODULEFILE$NC"
  echo " "
  module purge
  load_module "cmake/${CMAKEVERSION}"
  load_module "gcc/${GCCVERSION}"
  load_module "openmpi/${OPENMPIVERSION}/gcc/${GCCVERSION}"
  load_module "hdf5/${HDF5VERSION}/gcc/${GCCVERSION}/openmpi/${OPENMPIVERSION}"
  module list
  echo " "
  echo -e "$GREEN""Important: If the compilation step fails, run the script again and if it still fails \n1) try compiling single, .i.e., remove -j from make -j or \n2) try make -j 2 (not all available threads)$NC"
  echo " "
  echo -e "This will install ParaView version ${GREEN}${PARAVIEWVERSION}${NC}.\nCompilation in parallel will be executed with ${GREEN}${NBROFCORES} threads${NC}."
  if [[ ${RERUNMODE} -eq 0 ]]; then
    read -p "Have the correct modules been loaded? If yes, press [Enter] to continue or [Crtl+c] to abort!"
  fi

  # Install destination
  PARAVIEWINSTALLDIR=/opt/paraview/${PARAVIEWVERSION}/gcc-${GCCVERSION}/openmpi-${OPENMPIVERSION}/hdf5-${HDF5VERSION}

  # Change to sources directors
  cd ${SOURCESDIR}

  # Download tar.gz file
  if [ ! -f ${TARFILE} ]; then
    wget --output-document=${TARFILE} "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSION%.*}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}.tar.gz"
  fi

  # Check if tar.gz file was correctly downloaded
  if [ ! -f ${TARFILE} ]; then
    echo -e "$RED""no source-file downloaded for Paraview-${PARAVIEWVERSION}$NC"
    echo -e "$RED""check if https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSION%.*}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}-source.tar.gz$NC"
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

  # Check if extraction failed
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed: [tar -xzf ${TARFILE} paraview-${PARAVIEWVERSION}]$NC"
    exit
  else
    # Check if decompressed directory exists
    if [ -d "${SOURCESDIR}/ParaView-v${PARAVIEWVERSION}" ]; then
      # Check if renamed directory exists and create backup of it
      if [ -d "${SOURCESDIR}/paraview-${PARAVIEWVERSION}" ]; then
        # Move directory, e.g., "paraview-5.8.0" to "paraview-5.8.0_bak"
        #echo -e "${SOURCESDIR}/paraview-${PARAVIEWVERSION} already exists."
        #rm -rf ${SOURCESDIR}/paraview-${PARAVIEWVERSION}_bak
        #mv ${SOURCESDIR}/paraview-${PARAVIEWVERSION} ${SOURCESDIR}/paraview-${PARAVIEWVERSION}_bak

        # Inquiry: Continue the installation with the existing files OR remove them all and start fresh
        while true; do
          echo " "
          echo "${YELLOW}${SOURCESDIR}/paraview-${PARAVIEWVERSION} already exists.${NC}"
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
              [Nn]* ) rm -rf ${SOURCESDIR}/paraview-${PARAVIEWVERSION} ; break;;
              * ) echo "Please answer yes or no.";;
          esac
        done
      fi

      # Check if renamed directory was possibly deleted in the previous step OR this is a completely fresh installation
      if [ ! -d "${SOURCESDIR}/paraview-${PARAVIEWVERSION}" ]; then
        # The extracted directory is named, e.g., "ParaView-v5.8.0", which is renamed to "paraview-5.8.0" with small letters and no "v"
        mv ${SOURCESDIR}/ParaView-v${PARAVIEWVERSION} ${SOURCESDIR}/paraview-${PARAVIEWVERSION}
      fi
    fi
    #rm -rf paraview-${PARAVIEWVERSION}-source.tar.gz
  fi

  if [ ! -e "${SOURCESDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}" ]; then
    mkdir -p "${SOURCESDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}"
  fi

  # ------------------------------------------------------------------------------------------------------------------------------------------
  # Fix paraview files for 5.9.1 by adding #include <limits> in the include section due to the error ‘numeric_limits’ is not a member of ‘std’
  if [[ ${PARAVIEWVERSION} == '5.9.1' ]]; then

  # find /opt/sources/paraview-5.9.1/. -name "vtkGenericDataArrayLookupHelper.h"
  #   LINENBR=$(grep -n "#include" /opt/sources/paraview-5.9.1/VTK/Filters/HyperTree/vtkHyperTreeGridThreshold.cxx | tail -1 | cut -f1 -d:)
  #   sudo sed -i "$(echo $((LINENBR + 1)))i #include <limits>" /opt/sources/paraview-5.9.1/VTK/Rendering/Core/vtkColorTransferFunction.cxx
  #   sudo sed -i '29i #include <limits>' /opt/sources/paraview-5.9.1/VTK/Filters/HyperTree/vtkHyperTreeGridThreshold.cxx
  # sudo vim /opt/sources/paraview-5.9.1/VTK/Common/Core/vtkGenericDataArrayLookupHelper.h
  # #include <limits>

  # myarray=(vtkPiecewiseFunction.cxx
  #          vtkColorTransferFunction.cxx
  #          vtkHyperTreeGridThreshold.cxx
  #          vtkGenericDataArrayLookupHelper.h)
  #   for t in "${myarray[@]}"; do
  #     FILE=$(find /opt/sources/paraview-5.9.1/. -name "$t")
  #     if [[ -f ${FILE} ]]; then
  #       LINENBR=$(grep -n "#include" ${FILE} | tail -1 | cut -f1 -d:)
  #       if [[ -n ${LINENBR} ]]; then
  #         LIMITS=$(grep -n "#include" ${FILE} | grep limits)
  #         if [[ -z ${LIMITS} ]]; then
  #           echo "Adding #include <limits> to $FILE at line $LINENBR"
  #           sudo sed -i "$(echo $((LINENBR + 1)))i #include <limits>" ${FILE}
  #         else
  #           echo "${LIMITS} already found in ${FILE}"
  #         fi
  #       fi
  #     fi
  #   done

    cd "${SOURCESDIR}/paraview-${PARAVIEWVERSION}"
    PATCHFILE='vtk-gcc11.patch'
    wget -O ${PATCHFILE} https://gitlab.kitware.com/vtk/vtk/-/merge_requests/7554.patch

    # Check if .patch file was correctly downloaded
    if [ ! -f ${PATCHFILE} ]; then
      echo -e "$RED""no patch-file downloaded from https://gitlab.kitware.com/vtk/vtk/-/merge_requests/7554.patch ... Check the link or internet access$NC"
      exit
    fi

    # Apply the patch: patch accepts the --forward option to apply patches only once (in this case an error is returned)
    patch --forward -p1 -d VTK < ${PATCHFILE} || true

  fi
  # ------------------------------------------------------------------------------------------------------------------------------------------

  cd "${SOURCESDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}"

  # CMAKE COMPILE FLAGS DEPEND ON THE CHOSEN PARAVIEW VERSION!
  if [ "$PARAVIEWVERSION" == "5.2.0" ] || [ "$PARAVIEWVERSION" == "5.3.0" ] || [ "$PARAVIEWVERSION" == "5.4.0" ]; then
    cmake -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_TESTING=OFF \
      -DBUILD_EXAMPLES=OFF \
      -DPARAVIEW_ENABLE_PYTHON=ON \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_USE_VISITBRIDGE=OFF \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DCMAKE_INSTALL_PREFIX=${PARAVIEWINSTALLDIR} \
      ${SOURCESDIR}/paraview-${PARAVIEWVERSION}
  elif [ "$PARAVIEWVERSION" == "5.7.0" ] || [ "$PARAVIEWVERSION" == "5.8.0" ]; then
    cmake -DCMAKE_BUILD_TYPE=Release \
      -DPARAVIEW_USE_PYTHON=ON \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
      -DCMAKE_INSTALL_PREFIX=${PARAVIEWINSTALLDIR} \
      ${SOURCESDIR}/paraview-${PARAVIEWVERSION}
  elif [ "$PARAVIEWVERSION" == "5.9.1" ]; then
    cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_POLICY_DEFAULT_CMP0072=NEW \
      -DCMAKE_POLICY_DEFAULT_CMP0074=NEW \
      -DPARAVIEW_USE_PYTHON=ON \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DVTK_USE_SYSTEM_HDF5=ON \
      -DCMAKE_INSTALL_PREFIX=${PARAVIEWINSTALLDIR} \
      ${SOURCESDIR}/paraview-${PARAVIEWVERSION}
  else
    echo -e "$RED""ERROR: Set the correct cmake compile flags for the paraview version [$PARAVIEWVERSION] in InstallParaview.sh$NC"
    exit
  fi

  # Compile source files with NBROFCORES threads
  make 2>&1 | tee make.out

  # Check if compilation failed
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo -e "$RED""Failed: [make -j${NBROFCORES} 2>&1 | tee make.out]$NC"
    exit
  else
    make install 2>&1 | tee install.out
  fi

  # create modulefile if the installation seems successful (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -e "${PARAVIEWINSTALLDIR}/bin/paraview" ]; then
    if [ ! -d "${PARAVIEWMODULEFILEDIR}" ]; then
      mkdir -p ${PARAVIEWMODULEFILEDIR}
    fi
    cp ${TEMPLATEPATH} ${MODULEFILE}
    sed -i 's/paraviewversion/'${PARAVIEWVERSION}'/gI' ${MODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${MODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${MODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${OPENMPIVERSION}'/gI' ${MODULEFILE}
    sed -i 's/HDF5VERSIONFLAG/'${HDF5VERSION}'/gI' ${MODULEFILE}
  else
    echo -e "$RED""No module file created for Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION}$NC"
    echo -e "$RED""no installation found in ${PARAVIEWINSTALLDIR}/bin$NC"
  fi
else
  echo -e "$YELLOW""Paraview-${PARAVIEWVERSION} already created: module file exists under ${MODULEFILE}$NC"
fi
