\hypertarget{reggie}{}

# Regression Testing \label{chap:reggie}

The purpose of regression testing is summarized by the following quote:

> *Regression testing (rarely non-regression testing) is re-running functional and non-functional tests to ensure that previously developed and tested software still performs after a change.*

Wikipedia: https://en.wikipedia.org/wiki/Regression_testing

## Reggie2.0 Tool

PICLas is continuously tested by utilizing a Python based regression testing environment, which is
run by gitlab-runners. Therefore, the tool *reggie2.0* is used, which is found here: https://gitlab.com/reggie2.0/reggie2.0

Additionally, the different analysis methods that can be applied and the general structure of a 
regression test is described there.

Different tests are executed on check-in, during nightly or weekly testing. These tests are defined
in the file *.gitlab-ci.yml* that is located in the top level repository directory of PICLas. In
this file, various tests are defined, which are found under *regressioncheck/checks* and a summary
of the different tests is given under https://github.com/piclas-framework/piclas/blob/master/REGGIE.md

The automatic execution by *gitlab-runners* can be performed on any machine that is connected to the
internet and in the following section, the setup of such a machine is described

## Regression Server *Gitlab Runner* Setup
In a first step, the required software packages for PICLas and Reggie2.0 are installed on a new
system. In a second step, the *gitlab-runner* program is installed and the setup of runner is
described.

### Required Installation of Software on Clean Ubuntu Setup (18.04)
Latest tests on 

  * Ubuntu (18.04), 3 Jul 2019
  * Ubuntu server (18.04.3 LTS), 19 Nov 2019

The following packages can be installed automatically by using the script located at `./tools/Setup_ModuleEnv/InstallPackagesReggie.sh`.
The system inquiries can be skipped by forcing `yes` as input via

    yes | ./InstallPackagesRe

The script contains the following packages

```
# Check for updates
sudo apt-get update

# compiler
sudo apt-get install make cmake cmake-curses-gui gfortran g++ gcc 

# python
sudo apt-get install python python-numpy python3 python3-numpy python3-matplotlib python-matplotlib

# editors
sudo apt-get install gedit vim vim-runtime vim-common vim-tiny gdb vim-gui-common

# git && versions
sudo apt-get install git gitg qgit subversion

# paraview
sudo apt-get install libpython-dev libboost-dev libphonon-dev libphonon4 libxt-dev mesa-common-dev
#apt-get install qt4-default qt4-dev-tools libqt4-dev qt4-qmake libqt4-opengl-dev
sudo apt-get install qttools5-dev libqt5x11extras5-dev qt5-default libgl1-mesa-dev

# Tecplot
sudo apt-get install libstdc++5

# tools
sudo apt-get install gzip gimp htop meld gnuplot gnuplot-x11 vlc okular ddd gmsh unzip
sudo apt-get install openvpn openssl openssh-client

# for FLEXI/PICLas
sudo apt-get install liblapack3 liblapack-dev zlib1g-dev exuberant-ctags

# for documentation
sudo apt-get install pandoc pandoc-citeproc
sudo apt-get install texlive-full

# hdf5-file viewer
sudo apt-get install hdfview 

# Install libs for reggie
sudo apt-get install python-h5py

```

When no module environment is to be used on the server, the following packages are also required

```
# Further libs
sudo apt-get install hdf5-tools libhdf5-dev # this is maybe not required (do not install them if it works without these packages)


# openMPI
wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.3.tar.gz
tar -xvf openmpi-3.1.3.tar.gz
cd openmpi-3.1.3
mkdir -p build
cd build
../configure --prefix=/opt/openmpi/3.1.3
sudo make all install
export PATH="/opt/openmpi/3.1.3/bin:$PATH"
export LD_LIBRARY_PATH="/opt/openmpi/3.1.3/lib:$LD_LIBRARY_PATH"
cd ../..
rm openmpi-3.1.3.tar.gz


#HDF5
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.bz2
tar -xjf hdf5-1.10.5.tar.bz2
cd hdf5-1.10.5
mkdir -p build
cd build
cmake -DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5
sudo make && sudo make install
export HDF5_DIR=/opt/hdf5/1.10.5/share/cmake
cd ../..
rm hdf5-1.10.5.tar.bz2
```

otherwise a module environment can be installed at this point, see
`~/Flexi/piclas/tools/Setup_ModuleEnv/README.txt`, which is explained in detail in Chapter \ref{chap:tools} under Section \ref{sec:tools_module_env}.

When no module environment is to be used on the server, the following commands must be places in the
*.gitlab-ci.yml* file:

```
# Export the paths on new reggie2@reggie2 (no module env any more)
before_script:
  - ulimit -s unlimited
  - export PATH=/opt/openmpi/3.1.3/bin:$PATH
  - export LD_LIBRARY_PATH=/opt/openmpi/3.1.3/lib/:$LD_LIBRARY_PATH
  - export CMAKE_PREFIX_PATH=/opt/openmpi/3.1.3/share/cmake:$CMAKE_PREFIX_PATH
  - export CMAKE_LIBRARY_PATH=/opt/openmpi/3.1.3/lib:$CMAKE_LIBRARY_PATH
  - export HDF5_DIR=/opt/hdf5/1.10.5/share/cmake/
  - export PATH=/opt/hdf5/1.10.5/bin:$PATH
  - export LD_LIBRARY_PATH=/opt/hdf5/1.10.5/lib/:$LD_LIBRARY_PATH
  - export CMAKE_PREFIX_PATH=/opt/hdf5/1.10.5/:$CMAKE_PREFIX_PATH
  - export CMAKE_LIBRARY_PATH=/opt/hdf5/1.10.5/lib:$CMAKE_LIBRARY_PAT
```


otherwise, the correct environment must be loaded by adding the following in `/etc/profile`

    ```
    # Default modules
    module load gcc/9.2.0  cmake/3.15.3-d  openmpi/4.0.1/gcc/9.2.0  hdf5/1.10.5/gcc/9.2.0/openmpi/4.0.1
    ```

or by loading the modules directly in the gilab script file, e.g.,

    ```
    module load XX/XX
    ```

NOTE: The stack size limit has been removed here by `ulimit -s unlimited`, which might be required
by memory consuming programs

### Installation Steps for Gitlab Runners 
Latest tests on 

  * Ubuntu (18.04) with gitlab-runner 10.5.0 (10.5.0), 3 Jul 2019
  * Ubuntu server (18.04.3 LTS) with gitlab-runner 10.5.0 (10.5.0), 19 Nov 2019

1. Install gitlab-runner from ubuntu packages (choose old version to avoid problems https://gitlab.com/gitlab-org/gitlab-runner/issues/1379)
   This creates the user gitlab-runner and a home directory (for 10.5 in /var/lib/gitlab-runner/)

    ```
    sudo apt-get install gitlab-runner
    ```

2. Start the runner program

    ```
    sudo gitlab-runner start
    ```
3. Register a runner using a shell executor (follow the information in the official guideline as it 
   can vary slightly from version to version), on gitlab see `Settings` $\rightarrow$ `CI / CD Settings` 
   $\rightarrow$ `Runners` (registration token during setup)
    ```
    sudo gitlab-runner register
    ```
4. Restart the runner
    ```
    sudo gitlab-runner restart
    ```
5. create ssh keys for normal user and set up password free access to gitlab (iag) and gitlab.com (reggie)
    ```
    ssh-keygen -t ecdsa -b 521
    ```
    Add key to `Enabled deploy keys`. If multiple codes are on gitlab.com, add the key to one 
    repository and select the key on the other repositories via `Privately accessible deploy keys`.

    Clone a code from each platform to create known hosts then
    ```
    sudo cp .ssh/*  /var/lib/gitlab-runner/.ssh
    ```
    Check rights in this folder. If necessary make gitlab-runner own the files with 
    ```
    sudo chown -R gitlab-runner:gitlab-runner /var/lib/gitlab-runner/.ssh/
    ```
    If the runner is used to push to remote repositories, add the public key under *deploy keys*
    and execute, e.g., 
    ```
    sudo -u gitlab-runner git clone git@github.com:piclas-framework/piclas.git piclas_github
    ```
    to establish the first connection with the new repository and add the repo IP to the list of
    known hosts.
6. Start pipeline in gitlab or github for testing of reggie

NOTE: Interesting information is found in `/etc/systemd/system/gitlab-runner.service`.




### **Configuration files**

The runner services can be adjusted by changing the settings in the file

    /etc/systemd/system/gitlab-runner.service

in which the runner configuration file is specified:

```
[Unit]
Description=GitLab Runner
After=syslog.target network.target
ConditionFileIsExecutable=/usr/bin/gitlab-runner

[Service]
StartLimitInterval=5
StartLimitBurst=10
ExecStart=/usr/bin/gitlab-runner "run" "--working-directory" "/var/lib/gitlab-runner/" "--config" "/etc/gitlab-runner/config.toml" "--service" "gitlab-runner" "--syslog" "--user" "gitlab-runner"


Restart=always
RestartSec=120

[Install]
WantedBy=multi-user.target
```

The runner configuration settings can be edited by changing the settings in the file

    /etc/gitlab-runner/config.toml

where the number of runners, the concurrency level and runner limits are specified:

```
concurrent = 2
check_interval = 0

[[runners]]
  name = "myrunner1"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]

[[runners]]
  name = "myrunner2"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]

[[runners]]
  name = "myrunner3"
  url = "https://gitlab.com/"
  token = "XXXXXXXXXX"
  executor = "shell"
  limit = 1
  [runners.cache]
```

### Automatic Deployment to other platforms

1. Add the required ssh key to the deploy keys on the respective platform (e.g. github)
1. Clone a code from the platform to update the list of known hosts. Do not forget to copy the
    information to the correct location for the runner to have access to the platform
    ```
    sudo cp~/.ssh/.ssh/known_hosts /var/lib/gitlab-runner/.ssh/known_hosts
    ```
    This might have to be performed via the gitlab-runner user, which can be accomplished by
    executing the following command
    ```
    sudo -u gitlab-runner git clone git@github.com:piclas-framework/piclas.git piclas_github
    ```
1. PICLas deployment is performed by the gitlab runner in the *deployment stage*
    ```
    github:
      stage: deploy
      tags:
        - withmodules-concurrent
      script:
        - if [ -z "${DO_DEPLOY}" ]; then exit ; fi
        - rm -rf piclas_github || true ;
        - git clone -b master --single-branch git@gitlab.com:piclas/piclas.git piclas_github ;
        - cd piclas_github ;
        - git remote add piclas-framework git@github.com:piclas-framework/piclas.git ;
        - git push --force --follow-tags piclas-framework master ;
    ```

    This script clones the master branch of PICLas and deploys it on github.
