#!/bin/bash

#==============================================================================
# title       : InstallGCC.sh
# description : This script installs the gcc compiler with a specified version 
#               as given below via GCCVERSION='X.X.X'
# date        : Nov 27, 2019
# version     : 1.0   
# usage       : bash InstallGCC.sh
# notes       : 
#==============================================================================

INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
MODULETEMPLATEDIR=/opt/Installsources/moduletemplates

cd $INSTALLDIR
if [ ! -e "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# NOTE:
#GCC depends on:
#
#    GMP: GNU Multiple Precision Arithmetic Library
#    MPFR: GNU Multiple-precision floating-point rounding library
#    MPC: GNU Multiple-precision C library
#    ELF: Executable and Linkable Format library
#    PPL: Parma Polyhedra Library (optional, for memory optimizations)


# DOWNLOAD and INSTALL GCC COMPILER (example gcc-7.4.0)
#GCCVERSION='7.4.0'
#GCCVERSION='8.3.0'

# 9.3.0: Building GCC requires: GMP 4.2+, MPFR 2.4.0+ and MPC 0.8.0+
# sudo apt-get install libmpfr-dev
# sudo apt-get install libmpc-dev
GCCVERSION='9.3.0'

# 10.1.0: Building GCC requires GMP 4.2+, MPFR 3.1.0+ and MPC 0.8.0+
# sudo apt-get install libmpfr-dev
# sudo apt-get install libmpc-dev
#GCCVERSION='10.1.0'

MODULEFILEDIR=${INSTALLDIR}/modules/modulefiles/compilers/gcc
MODULEFILE=${MODULEFILEDIR}/${GCCVERSION}

COMPILERDIR=${INSTALLDIR}'/compiler/gcc/'${GCCVERSION}

if [[ -n ${1} ]]; then
  if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
    rm ${MODULEFILE}
  fi
fi

if [ ! -e "${MODULEFILE}" ]; then
  echo "creating Compiler GCC-${GCCVERSION}"
  cd ${SOURCESDIR}
  if [ ! -e "${SOURCESDIR}/gcc-${GCCVERSION}.tar.gz" ]; then
    wget -O gcc-${GCCVERSION}.tar.gz "ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/gcc-${GCCVERSION}/gcc-${GCCVERSION}.tar.gz"
  fi
  if [ ! -e "${SOURCESDIR}/gcc-${GCCVERSION}.tar.gz" ]; then
    echo "no gcc install-file downloaded for GCC-${GCCVERSION}"
    echo "check if ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/gcc-${GCCVERSION}/gcc-${GCCVERSION}.tar.gz exists"
    exit
  fi
  tar -xzf gcc-${GCCVERSION}.tar.gz && rm -rf gcc-${GCCVERSION}.tar.gz
  if [ ! -d "${SOURCESDIR}/gcc-${GCCVERSION}/build" ]; then
    mkdir -p gcc-${GCCVERSION}/build
  fi
  if [[ ${1} =~ ^-r(erun)?$ ]] ; then
    rm gcc-${GCCVERSION}/build/* 
  fi
  cd gcc-${GCCVERSION}/build
  ../configure -v \
    --prefix=${COMPILERDIR} \
    --enable-languages=c,c++,objc,obj-c++,fortran \
    --enable-shared \
    --disable-multilib \
    --disable-bootstrap \
    --enable-checking=release \
    --with-sysroot=/ \
    --with-system-zlib
    # --enable-valgrind-annotations
  make -j 2>&1 | tee make.out
  if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo " "
    echo "Failed: [make -j 2>&1 | tee make.out]"
    exit
  else
    make install 2>&1 | tee install.out
  fi

  if [ ! -d "${MODULEFILEDIR}" ]; then
    mkdir -p ${MODULEFILEDIR}
  fi

  if [ -e "${COMPILERDIR}/bin/gcc" ] && [ -e "${COMPILERDIR}/bin/gfortran" ]; then
    cp ${MODULETEMPLATEDIR}/compilers/gcc/v_temp ${MODULEFILE}
    sed -i 's/versionflag/'${GCCVERSION}'/gI' ${MODULEFILE}
  else
    echo "compiler not installed, no modulefile created"
  fi
else
  echo "Compiler GCC-${GCCVERSION} already created (module file exists)"
fi
