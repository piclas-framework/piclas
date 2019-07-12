#!/bin/bash
INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
templatedir=/opt/Installsources/moduletemplates

if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL CMAKE (example cmake-3.4.3)
# CMAKEVERSION='3.4.3'
CMAKEVERSION='3.13.3'
CMAKEDIR=${INSTALLDIR}/cmake/${CMAKEVERSION}/standard
MODULEFILE=${INSTALLDIR}/modules/modulefiles/utilities/cmake/${CMAKEVERSION}-d

if [ ! -e ${MODULEFILE} ]; then
  echo "creating CMake-${CMAKEVERSION}"
  cd ${SOURCESDIR}
  if [ ! -e ${SOURCESDIR}/cmake-${CMAKEVERSION}.tar.gz ]; then
    wget "https://github.com/Kitware/CMake/releases/download/v${CMAKEVERSION}/cmake-${CMAKEVERSION}.tar.gz"
  fi
  tar -xzf cmake-${CMAKEVERSION}.tar.gz
  if [ ! -d ${SOURCESDIR}/cmake-${CMAKEVERSION}/build ]; then
    mkdir -p ${SOURCESDIR}/cmake-${CMAKEVERSION}/build
  fi
  cd ${SOURCESDIR}/cmake-${CMAKEVERSION}/build
  ../bootstrap --prefix=${CMAKEDIR}
  make -j 2>&1 | tee make.out
  make install 2>&1 | tee install.out

  if [ -e ${CMAKEDIR}/bin/cmake ] && [ -e ${CMAKEDIR}/bin/ccmake ]; then
    if [ ! -e ${INSTALLDIR}/modules/modulefiles/utilities/cmake ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities/cmake
    fi
    cp ${TEMPLATEDIR}/utilities/cmake/cmake_temp ${MODULEFILE}
    sed -i 's/cmakeversion/'${CMAKEVERSION}'/g' ${MODULEFILE}
    sed -i 's/compilerversion/standard/g' ${MODULEFILE}
    rm -rf cmake-${CMAKEVERSION}.tar.gz
  else
    echo "ERROR in cmake installation, no modulefile created"
    if [ -e ${CMAKEDIR}/bin/cmake ]; then
      echo "ERROR: cmake not installed"
    fi
    if [ -e ${CMAKEDIR}/bin/ccmake ]; then
      echo "ERROR: cmake-curses-gui not installed"
    fi
  fi
else
  echo "CMake-${CMAKEVERSION} already created (module file exists)"
fi
