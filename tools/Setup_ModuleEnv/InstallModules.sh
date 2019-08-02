#!/bin/bash
INSTALLDIR=/opt
SOURCESDIR=/opt/Installsources
calcTrue() { awk 'BEGIN{printf "%d\n" , ('"$*"'?1:0)}';}

if [ ! -d ${SOURCESDIR} ]; then
  mkdir -p ${SOURCESDIR}
fi

# download and install modules framework if no modules are present
if [ ! -d ${MODULESHOME} ]; then
  if [ ! -d ${INSTALLDIR}/modules/3.2.10 ]; then
    echo "creating module evironment with modules-3.2.10"
    cd ${SOURCESDIR}
    if [ ! -e ${SOURCESDIR}/modules-3.2.10.tar.gz ]; then
      wget --output-document=modules-3.2.10.tar.gz "https://downloads.sourceforge.net/project/modules/Modules/modules-3.2.10/modules-3.2.10.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fmodules%2Ffiles%2FModules%2Fmodules-3.2.10%2Fmodules-3.2.10.tar.gz%2Fdownload%3Fuse_mirror%3Dkent&ts=1548854959"
    fi
    tar -xzf modules-3.2.10.tar.gz && rm -rf modules-3.2.10.tar.gz
    if [ ! -e ${SOURCESDIR}/modules-3.2.10.tar.gz ]; then
      mkdir -p ${SOURCESDIR}/modules-3.2.10/build
    fi
    cd ${SOURCESDIR}/modules-3.2.10/build

    if [ `calcTrue "$(echo 'puts $tcl_version;exit 0' | tclsh) < 8.5"` -eq 1 ]; then
      ../configure --prefix=${INSTALLDIR}/modules/3.2.10 --with-module-path=${INSTALLDIR}/modules/modulefiles
    else
      CPPFLAGS="-DUSE_INTERP_ERRORLINE" ../configure --prefix=${INSTALLDIR}/modules/3.2.10 --with-module-path=${INSTALLDIR}/modules/modulefiles
    fi
    make 2>&1 | tee make.out
    make install 2>&1 | tee install.out

    if [ -z "$(grep "if.*Modules.*3.2.10.*init.*bash.*then" /etc/profile)" ]; then
      sed -i '1 i\if [ -f /opt/modules/3.2.10/Modules/3.2.10/init/bash ]; then' /etc/profile
      sed -i '2 i\  . /opt/modules/3.2.10/Modules/3.2.10/init/bash' /etc/profile
      sed -i '3 i\fi' /etc/profile
    else
      echo "modules init already exists in /etc/profile"
      break
    fi
    if [ -z "$(grep "if.*Modules.*3.2.10.*init.*bash.*then" /etc/bash.bashrc)" ]; then
      sed -i '1 i\if [ -f /opt/modules/3.2.10/Modules/3.2.10/init/bash ]; then' /etc/bash.bashrc
      sed -i '2 i\  . /opt/modules/3.2.10/Modules/3.2.10/init/bash' /etc/bash.bashrc
      sed -i '3 i\fi' /etc/bash.bashrc
    else
      echo "modules init already exists in /etc/bash.bashrc"
      break
    fi
    source /etc/profile

    # Change modulefiles path in init -> ${MODULESHOME}/init/.modulespath
    # comment everything in .modulespath
    sed -i 's/^/\# /' ${MODULESHOME}/init/.modulespath
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/utilities" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/MPI" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/libraries" >> ${MODULESHOME}/init/.modulespath
    # echo "/home/\$\{USER\}/modulefiles" >> ${MODULESHOME}/init/.modulespath
    mkdir -p ${INSTALLDIR}/modules/modulefiles/
    mkdir -p ${INSTALLDIR}/modules/modulefiles/compilers
    mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities
    mkdir -p ${INSTALLDIR}/modules/modulefiles/MPI
    mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries
    # cd /opt/modules/modulefiles/compilers/gcc
  else
    echo "module evironment (modules-3.2.10) already created"
  fi
else
  echo "module evironment ($(module --version)) already existant"
  # Change modulefiles path in init -> ${MODULESHOME}/init/.modulespath
  # comment everything
  if [ ! -z $(grep "${INSTALLDIR}/modules/modulefiles/" ${MODULESHOME}/init/.modulespath) ]; then
    sed -i 's/^/\# /' ${MODULESHOME}/init/.modulespath
    # add:
    echo "/opt/modules/modulefiles/compilers" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/utilities" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/MPI" >> ${MODULESHOME}/init/.modulespath
    echo "/opt/modules/modulefiles/libraries" >> ${MODULESHOME}/init/.modulespath
    # echo "/home/\$\{USER\}/modulefiles" >> ${MODULESHOME}/init/.modulespath
    mkdir -p ${INSTALLDIR}/modules/modulefiles/
    mkdir -p ${INSTALLDIR}/modules/modulefiles/compilers
    mkdir -p ${INSTALLDIR}/modules/modulefiles/utilities
    mkdir -p ${INSTALLDIR}/modules/modulefiles/MPI
    mkdir -p ${INSTALLDIR}/modules/modulefiles/libraries
    # cd /opt/modules/modulefiles/compilers/gcc
  fi
fi

if [[ -z ${MODULESHOME} ]] || [[ -z $(grep -i "${INSTALLDIR}/modules/modulefiles" ${MODULESHOME}/init/.modulespath) ]]; then
  echo "restart operating system and this run script again"
fi
