\hypertarget{tools}{}

# Tools Overview \label{chap:tools}

This section gives an overview over the tools and scripts for developers contained in the **PICLas** repository. It also provides references to the tutorials where their usage is explained.

## WarningsCheck

The maximum number of allowed warnings is limited to 10 and can be checked by using the script *test_max_warnings.sh*. Navigate to the top-level directory of the repository and execute

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
