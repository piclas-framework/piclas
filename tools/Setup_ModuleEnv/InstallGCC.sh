#!/bin/bash
INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
MODULETEMPLATEDIR=/opt/Installsources/moduletemplates

cd $INSTALLDIR
if [ ! -e "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL GCC COMPILER (example gcc-7.4.0)
GCCVERSION='7.4.0'
GCCVERSION='8.3.0'
#GCCVERSION='9.2.0'
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
    break
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
  make -j 2 2>&1 | tee make.out
  make install 2>&1 | tee install.out

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
