# Developer Tools

This chapter gives an overview over the tools and scripts for developers contained in the **PICLas** repository
It also provides references to the tutorials where their usage is explained.

## WarningsCheck

The maximum number of allowed warnings is limited to 10 and can be checked by using the script *test_max_warnings.sh*.
Navigate to the top-level directory of the repository and execute

```
./tools/test_max_warnings.sh
```

This creates a build directory *build_test_max_warnings* and launches *cmake*, which must be
configured by the user (hit *c* for configure). The user must then supply the desired compilation flags and must 
complete configuration by hitting *c* again and then continue the script by generating the make
files (hit *g* for generate). The
number of warnings is then supplied.

## Remove trailing white spaces

Script (RemoveTrailWhiteSpaces.sh) can be executed from any directory inside the project.
Searches all files (\*.f90, \*.h) in gitroot/src directory and removes trailing white spaces.
Before remove operation all files and number of changes are shown and user is asked wether action is to be performed.


## Module Environment

A set of scripts given in `./tools/Setup_ModuleEnv/` can be used for setting up a
module environment. The main steps are described in `./tools/Setup_ModuleEnv/README.txt`

Global settings to the appearance of the module env list can be changed in `/etc/profile`, e.g., by
adding (sudo required)

```
# -------------------------------------------------------------------------
# Display module information in list form
alias module='module -l'

# -------------------------------------------------------------------------
```

### FAQ: Common Problems

* After installing packages, cmake and GCC, the module environment is loaded by the script
  `./InstallMPIallCOMPILERS.sh` the first time. This can fail even though the environment can be
  loaded from a shell by hand (interactive shell).


