#!/bin/bash
INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
TEMPLATEDIR=/opt/Installsources/moduletemplates

if [ ! -d "${SOURCESDIR}" ]; then
  mkdir -p ${SOURCESDIR}
fi

# DOWNLOAD and INSTALL HDF5 for every compiler openmpi / mpich combination (example HDF5-1.8.18)
#HDF5VERSION=1.8.18
#HDF5VERSION=1.10.4
HDF5VERSION=1.10.5

HDF5DIR=${INSTALLDIR}'/hdf5/'${HDF5VERSION}

COMPILERNAMES='gcc intel'
for WHICHCOMPILER in ${COMPILERNAMES}; do
  NCOMPILERS=$(ls /opt/modules/modulefiles/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | wc -l)
  for i in $(seq 1 ${NCOMPILERS}); do
    COMPILERVERSION=$(ls /opt/modules/modulefiles/compilers/${WHICHCOMPILER}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n ${i} | tail -n 1)
    if [ ! -e "${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}" ]; then
      mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}
    fi

    # build hdf5 in single
    MODULEFILE=${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/single
    if [[ -n ${1} ]]; then
      if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
        rm ${MODULEFILE}
      fi
    fi
    if [ ! -e "${MODULEFILE}" ]; then
      echo "creating HDF5-${HDF5VERSION} library for ${WHICHCOMPILER}-${COMPILERVERSION} single"
      module purge
      if [[ -n $(module load ${WHICHCOMPILER}/${COMPILERVERSION}) ]]; then
        echo "module ${WHICHCOMPILER}/${COMPILERVERSION} not found "
        break
      fi
      module load ${WHICHCOMPILER}/${COMPILERVERSION}

      cd ${SOURCESDIR}
      if [ ! -e "${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz" ]; then
        wget "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz"
      fi
      if [ ! -e "${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz" ]; then
        echo "no hdf5 install-file downloaded for HDF5-${HDF5VERSION}"
        echo "check if https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5VERSION%.*}/hdf5-${HDF5VERSION}/src/hdf5-${HDF5VERSION}.tar.gz exists"
        break
      fi
      tar -xzf hdf5-${HDF5VERSION}.tar.gz hdf5-${HDF5VERSION}
      if [ ! -d "${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single" ]; then
        mkdir -p ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single
      fi
      if [[ ${1} =~ ^-r(erun)?$ ]] ; then
        rm ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single/* 
      fi
      cd ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/single/
      if [ "${WHICHCOMPILER}" == "gcc" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
      elif [ "${WHICHCOMPILER}" == "intel" ]; then
        ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/single --with-pic --enable-fortran --enable-fortran2003 --disable-shared CC=$(which icc) CXX=$(which icpc) FC=$(which ifort)
      fi
      make -j 2>&1 | tee make.out
      make install 2>&1 | tee install.out

      cp ${TEMPLATEDIR}/libraries/hdf5/single_template ${MODULEFILE}
      sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
      sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
      sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE} 
    else
      echo "HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} already created (module file exists)"
    fi

    MPINAMES='openmpi mpich'
    for WHICHMPI in ${MPINAMES}; do
      # build hdf5 with mpi
      if [ ! -d "${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}" ]; then
        mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}
      fi
      NMPI=$(ls /opt/modules/modulefiles/MPI/${WHICHMPI} | sed 's/ /\n/g' | grep -i "[0-9]\." | wc -l)
      for j in $(seq 1 ${NMPI}); do
        MPIVERSION=$(ls /opt/modules/modulefiles/MPI/${WHICHMPI}/ | sed 's/ /\n/g' | grep -i "[0-9]\." | head -n ${j} | tail -n 1)
        MODULEFILE=${INSTALLDIR}/modules/modulefiles/libraries/hdf5/${HDF5VERSION}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
        if [[ -n ${1} ]]; then
          if [[ ${1} =~ ^-r(erun)?$ ]] && [[ -f ${MODULEFILE} ]]; then
            rm ${MODULEFILE}
          fi
        fi
        if [ ! -e "${MODULEFILE}" ]; then
          module purge
          if [[ -n $(module load ${WHICHCOMPILER}/${COMPILERVERSION}) ]]; then
            echo "module ${WHICHCOMPILER}/${COMPILERVERSION} not found "
            break
          fi
          module load ${WHICHCOMPILER}/${COMPILERVERSION}
          if [[ -n $(module load ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION}) ]]; then
            echo "module ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION} not found "
            break
          fi
          module load ${WHICHMPI}/${MPIVERSION}/${WHICHCOMPILER}/${COMPILERVERSION}

          if [ ! -d "${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}" ]; then
            mkdir -p ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
          fi
          if [[ ${1} =~ ^-r(erun)?$ ]] ; then
            rm ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}/* 
          fi
          cd ${SOURCESDIR}/hdf5-${HDF5VERSION}/build_${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION}
          ${SOURCESDIR}/hdf5-${HDF5VERSION}/configure --prefix=${HDF5DIR}/${WHICHCOMPILER}/${COMPILERVERSION}/${WHICHMPI}/${MPIVERSION} --with-pic --enable-fortran --enable-fortran2003 --disable-shared --enable-parallel CC=$(which mpicc) CXX=$(which mpicxx) FC=$(which mpifort)
          make -j 2>&1 | tee make.out
          make install 2>&1 | tee install.out
          cp ${TEMPLATEDIR}/libraries/hdf5/mpi_template ${MODULEFILE}
          sed -i 's/whichcompiler/'${WHICHCOMPILER}'/gI' ${MODULEFILE}
          sed -i 's/compilerversion/'${COMPILERVERSION}'/gI' ${MODULEFILE}
          sed -i 's/hdf5version/'${HDF5VERSION}'/gI' ${MODULEFILE} 
          sed -i 's/whichmpi/'${WHICHMPI}'/gI' ${MODULEFILE}
          sed -i 's/mpiversion/'${MPIVERSION}'/gI' ${MODULEFILE}
        else
          echo "HDF5-${HDF5VERSION} for ${WHICHCOMPILER}-${COMPILERVERSION} and ${WHICHMPI}-${MPIVERSION} already created (module file exists)"
          continue
        fi
      done
    done

    rm -rf ${SOURCESDIR}/hdf5-${HDF5VERSION}.tar.gz
  done

done
