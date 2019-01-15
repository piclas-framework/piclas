\hypertarget{codeoptions}{}

# Code options \label{chap:codeoptions}
## Compiler options \label{sec:compileroptions}
This section describes the main configuration options which can be set when building **PICLas** using CMake. 
Some options are dependent on others being enabled (or disabled), so the available ones may change. 

The first set of options describe general CMake behaviour:

* ``CMAKE_BUILD_TYPE``:

    This statically specifies what build type (configuration) will be built in this build tree. Possible values are
    * **Release**
    
        "Normal" execution.
    
    * **Profile**
    
        Performance profiling using gprof.
    
    * **Debug**
    
        Debug compiler for detailed error messages during code development.
    
    * **SANI**
    
        Sanitizer compiler for even more detailed error messages during code development.
    
* ``CMAKE_HOSTNAME``:

    This will display the host name of the machine you are compiling on.

* ``CMAKE_INSTALL_PREFIX``:

    If “make install” is invoked or INSTALL is built, this directory is prepended onto all install directories. This variable defaults to /usr/local on UNIX.

For some external libraries and programs that **PICLas** uses, the following options apply:

* ``CTAGS_PATH``:

    This variable specifies the Ctags install directory, an optional program used to jump between tags in the source file.

* ``PICLas_BUILD_HDF5``: ON/OFF

    This will be set to ON if no prebuilt HDF5 installation was found on your machine. In this case a HDF5 version will be build and used instead.

* ``HDF5_DIR``:

    If you want to use a prebuilt HDF5 library that has been build using the CMake system, this directory should contain the CMake configuration file for HDF5 (optional).

