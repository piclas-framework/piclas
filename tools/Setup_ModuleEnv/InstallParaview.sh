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

# Function for checken modules
check_module () {
                 if [ -z "${2}" ]; then
                   echo "module for ${1} not found. Exit"
                   exit
                 else
                   echo "${1}: ["${2}"]"
                 fi
                }

# Check command line arguments
LOADMODULES=1
for arg in "$@"
do
  if [ "$arg" == "--help" ] || [ "$arg" == "-h" ]; then
    echo "Input arguments:"
    echo "--help/-h            print help information"
    echo "--modules/-m         use modules defined in this script by the user."
    echo "                     Otherwise, find modules automatically."
    exit
  fi
  if [ "$arg" == "--modules" ] || [ "$arg" == "-m" ]; then
    LOADMODULES=0
    # Set desired versions
    #CMAKEVERSION=3.15.3-d
    CMAKEVERSION=3.17.0-d
    #GCCVERSION=9.2.0
    GCCVERSION=10.1.0
    #OPENMPIVERSION=3.1.4
    #OPENMPIVERSION=4.0.1
    OPENMPIVERSION=4.0.2
    HDF5VERSION=1.10.5
    break
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
PARAVIEWVERSION=5.8.0

INSTALLDIR=/opt
SOURCEDIR=/opt/Installsources
MODULESDIR=/opt/modules/modulefiles
MODULETEMPLATESDIR=/opt/Installsources/moduletemplates
MODULETEMPLATENAME=paraview_temp

if [ ! -d "${SOURCEDIR}" ]; then
  mkdir -p ${SOURCEDIR}
fi

# take the first gcc compiler installed with first compatible openmpi and hdf5
if [[ $LOADMODULES -eq 1 ]]; then
  CMAKEVERSION=$(ls ${MODULESDIR}/utilities/cmake/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  OPENMPIVERSION=$(ls ${MODULESDIR}/MPI/openmpi/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
  HDF5VERSION=$(ls ${MODULESDIR}/libraries/hdf5/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
fi

echo " "
echo "Modules found:"
check_module "cmake" ${CMAKEVERSION}
check_module "gcc  " ${GCCVERSION}
check_module "mpi  " ${OPENMPIVERSION}
check_module "hdf5 " ${HDF5VERSION}

PARAVIEWMODULEFILEDIR=${MODULESDIR}/utilities/paraview/${PARAVIEWVERSION}/gcc/${GCCVERSION}/openmpi/${OPENMPIVERSION}/hdf5
PARAVIEWMODULEFILE=${PARAVIEWMODULEFILEDIR}/${HDF5VERSION}

# if no paraview module for this compiler found, install paraview and create module
if [ ! -e "${PARAVIEWMODULEFILE}" ]; then
  echo "creating Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION} under"
  echo "$PARAVIEWMODULEFILE"
  echo " "
  echo "Important: If the compilation step fails, try compiling single, .i.e., remove -j from make -j"
  echo " "
  read -p "Press enter to continue"
  module purge
  module load cmake/${CMAKEVERSION}
  module load gcc/${GCCVERSION}
  module load openmpi/${OPENMPIVERSION}/gcc/${GCCVERSION}
  module load hdf5/${HDF5VERSION}/gcc/${GCCVERSION}/openmpi/${OPENMPIVERSION}
  module list
  echo " "

  # Install destination
  PARAVIEWINSTALLDIR=/opt/paraview/${PARAVIEWVERSION}/gcc-${GCCVERSION}/openmpi-${OPENMPIVERSION}/hdf5-${HDF5VERSION}

  # build and installation
  cd ${SOURCEDIR}
  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}-source.tar.gz" ]; then
    wget --output-document=paraview-${PARAVIEWVERSION}-source.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSION%.*}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}.tar.gz"
  fi
  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}-source.tar.gz" ]; then
    echo "no source-file downloaded for Paraview-${PARAVIEWVERSION}"
    echo "check if https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSION%.*}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}-source.tar.gz"
    exit
  fi
  tar -xzf paraview-${PARAVIEWVERSION}-source.tar.gz
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo " "
    echo "Failed: ["tar -xzf paraview-${PARAVIEWVERSION}-source.tar.gz paraview-${PARAVIEWVERSION}"]"
    exit
  else
    # Check if decompressed directory exists
    if [ -d "${SOURCEDIR}/ParaView-v${PARAVIEWVERSION}" ]; then
      # Check if renamed directory exists and create backup of it
      if [ -d "${SOURCEDIR}/paraview-${PARAVIEWVERSION}" ]; then
        # Move directory, e.g., "paraview-5.8.0" to "paraview-5.8.0_bak"
        mv ${SOURCEDIR}/paraview-${PARAVIEWVERSION} ${SOURCEDIR}/paraview-${PARAVIEWVERSION}_bak
      fi
      # The extracted directory is named, e.g., "ParaView-v5.8.0", which is renamed to "paraview-5.8.0" with small letters and no "v"
      mv ${SOURCEDIR}/ParaView-v${PARAVIEWVERSION} ${SOURCEDIR}/paraview-${PARAVIEWVERSION}
    fi
    #rm -rf paraview-${PARAVIEWVERSION}-source.tar.gz
  fi

  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}" ]; then
    mkdir -p ${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}
  fi
  cd ${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}

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
      ${SOURCEDIR}/paraview-${PARAVIEWVERSION}
  elif [ "$PARAVIEWVERSION" == "5.7.0" ] || [ "$PARAVIEWVERSION" == "5.8.0" ]; then
    cmake -DCMAKE_BUILD_TYPE=Release \
      -DPARAVIEW_USE_PYTHON=ON \
      -DPARAVIEW_USE_MPI=ON \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
      -DCMAKE_INSTALL_PREFIX=${PARAVIEWINSTALLDIR} \
      ${SOURCEDIR}/paraview-${PARAVIEWVERSION}
  else
    echo "ERROR: Set the correct cmake compile flags for the paraview version [$PARAVIEWVERSION]"
    exit
  fi
  make -j 2>&1 | tee make.out
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo "Failed: [make -j 2>&1 | tee make.out]"
    exit
  else
    make install 2>&1 | tee install.out
  fi


  # create modulefile if installation seems succesfull (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -e "${PARAVIEWINSTALLDIR}/bin/paraview" ]; then
    if [ ! -d "${PARAVIEWMODULEFILEDIR}" ]; then
      mkdir -p ${PARAVIEWMODULEFILEDIR}
    fi
    cp ${MODULETEMPLATESDIR}/utilities/paraview/${MODULETEMPLATENAME} ${PARAVIEWMODULEFILE}
    sed -i 's/paraviewversion/'${PARAVIEWVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${OPENMPIVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/HDF5VERSIONFLAG/'${HDF5VERSION}'/gI' ${PARAVIEWMODULEFILE}
  else
    echo "No module file created for Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION}"
    echo "no installation found in ${PARAVIEWINSTALLDIR}/bin"
  fi
else
  echo "Paraview-${PARAVIEWVERSION} already created: module file exists under ${PARAVIEWMODULEFILE}"
fi
