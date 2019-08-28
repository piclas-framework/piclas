#!/bin/bash

check_module () {
                 if [ -z "${2}" ]; then
                   echo "module for ${1} not found. Exit"
                   exit
                 else
                   echo "${1}: ["${2}"]"
                 fi
                }

# DOWNLOAD and INSTALL PARAVIEW (example Paraview-2.1.6)
PARAVIEWVERSION=5.2.0
PARAVIEWVERSIONTAG=5.2

INSTALLDIR=/opt
SOURCEDIR=/opt/Installsources
MODULESDIR=/opt/modules/modulefiles
MODULETEMPLATESDIR=/opt/Installsources/moduletemplates
MODULETEMPLATENAME=paraview_temp

if [ ! -d "${SOURCEDIR}" ]; then
  mkdir -p ${SOURCEDIR}
fi

# take the first gcc compiler installed with first compatible openmpi and hdf5
CMAKEVERSION=$(ls ${MODULESDIR}/utilities/cmake/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
GCCVERSION=$(ls ${MODULESDIR}/compilers/gcc/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
MPIVERSION=$(ls ${MODULESDIR}/MPI/openmpi/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
HDF5VERSION=$(ls ${MODULESDIR}/libraries/hdf5/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n 1 | tail -n 1)
PARAVIEWMODULEFILEDIR=${MODULESDIR}/utilities/paraview
PARAVIEWMODULEFILE=${PARAVIEWMODULEFILEDIR}/${PARAVIEWVERSION}

echo "Modules found:"
check_module "cmake" ${CMAKEVERSION}
check_module "gcc  " ${GCCVERSION}
check_module "mpi  " ${MPIVERSION}
check_module "hdf5 " ${HDF5VERSION}

# if no paraview module for this compiler found, install paraview and create module
if [ ! -e "${PARAVIEWMODULEFILE}" ]; then
  echo "creating Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION}"
  module purge
  module load cmake/${CMAKEVERSION}
  module load gcc/${GCCVERSION}
  module load openmpi/${MPIVERSION}/gcc/${GCCVERSION}
  module load hdf5/${HDF5VERSION}/gcc/${GCCVERSION}/openmpi/${MPIVERSION}

  PARAVIEWINSTALLDIR=/opt/paraview/${PARAVIEWVERSION}/gcc-${GCCVERSION}/openmpi-${MPIVERSION}/hdf5-${HDF5VERSION}

  # build and installation
  cd ${SOURCEDIR}
  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}-source.tar.gz" ]; then
    wget --output-document=paraview-${PARAVIEWVERSION}-source.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSIONTAG}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}.tar.gz"
  fi
  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}-source.tar.gz" ]; then
    echo "no source-file downloaded for Paraview-${PARAVIEWVERSION}"
    echo "check if https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PARAVIEWVERSIONTAG}&type=source&os=Sources&downloadFile=ParaView-v${PARAVIEWVERSION}-source.tar.gz"
    exit
  fi
  tar -xzf paraview-${PARAVIEWVERSION}-source.tar.gz
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo "Failed: ["tar -xzf paraview-${PARAVIEWVERSION}-source.tar.gz paraview-${PARAVIEWVERSION}"]"
  else
    if [ -d "${SOURCEDIR}/ParaView-v${PARAVIEWVERSION}" ]; then
      if [ ! -d "${SOURCEDIR}/paraview-${PARAVIEWVERSION}" ]; then
        mv ${SOURCEDIR}/ParaView-v${PARAVIEWVERSION} ${SOURCEDIR}/paraview-${PARAVIEWVERSION}
      fi
    fi
    #rm -rf paraview-${PARAVIEWVERSION}-source.tar.gz
  fi

  if [ ! -e "${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}" ]; then
    mkdir -p ${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}
  fi
  cd ${SOURCEDIR}/paraview-${PARAVIEWVERSION}/build_gcc/${GCCVERSION}

  cmake -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=OFF \
    -DBUILD_EXAMPLES=OFF \
    -DPARAVIEW_ENABLE_PYTHON=ON \
    -DPARAVIEW_USE_MPI=ON \
    -DPARAVIEW_USE_VISITBRIDGE=OFF \
    -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
    -DCMAKE_INSTALL_PREFIX=${PARAVIEWINSTALLDIR} \
    ${SOURCEDIR}/paraview-${PARAVIEWVERSION}
  make -j 2>&1 | tee make.out
  make install 2>&1 | tee install.out

  # create modulefile if installation seems succesfull (check if mpicc, mpicxx, mpifort exists in installdir)
  if [ -e "${PARAVIEWINSTALLDIR}/bin/paraview" ]; then
    if [ ! -d "${PARAVIEWMODULEFILEDIR}" ]; then
      mkdir -p ${PARAVIEWMODULEFILEDIR}
    fi
    cp ${MODULETEMPLATESDIR}/utilities/paraview/${MODULETEMPLATENAME} ${PARAVIEWMODULEFILE}
    sed -i 's/paraviewversion/'${PARAVIEWVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/CMAKEVERSIONFLAG/'${CMAKEVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/GCCVERSIONFLAG/'${GCCVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/MPIVERSIONFLAG/'${MPIVERSION}'/gI' ${PARAVIEWMODULEFILE}
    sed -i 's/HDF5VERSIONFLAG/'${HDF5VERSION}'/gI' ${PARAVIEWMODULEFILE}
  else
    echo "No module file created for Paraview-${PARAVIEWVERSION} for GCC-${GCCVERSION}"
    echo "no installation found in ${PARAVIEWINSTALLDIR}/bin"
  fi
else
  echo "Paraview-${PARAVIEWVERSION} already created (module file exists)"
  continue
fi
