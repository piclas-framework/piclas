#!/bin/bash

# Script for updating the piclas version number in ./src/globals/globals_vars.f90 and .github/workflows/cmake-ninja.yml
GLOBALS='./src/globals/globals_vars.f90'
WORKFLOW='.github/workflows/cmake-ninja.yml'

if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    NC="$(tput sgr0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
  fi
fi

# Check format of input argument
re='^([0-9]+\.){2}[0-9]+$'
# ^([0-9]+\.){2}$ is a POSIX ERE compliant pattern that matches:
#
# ^             - start of string
# ([0-9]+\.){2} - two occurrences of one or more digits and a dot
# [0-9]+        - one or more digits
# $             - end of string.
if ! [[ "$1" =~ $re ]] ; then
  echo "${RED}error: Incorrect version specified; must be in the format of #.#.#${NC}" >&2; exit 1
else
  MAJORVERSION=$(echo "$1" | cut -d "." -f1)
  MINORVERSION=$(echo "$1" | cut -d "." -f2)
  PATCHVERSION=$(echo "$1" | cut -d "." -f3)
fi

# Check if in top level directory
if [[ ! -f ${GLOBALS} ]]; then
  echo "${RED}Could not find ${GLOBALS}${NC}"
  echo "${RED}Please run this script from the piclas top level directory (where the src directory is located)${NC}"
  exit 1
fi

if [[ ! -f ${WORKFLOW} ]]; then
  echo "${RED}Could not find .github/workflows/cmake-ninja.yml${NC}"
  exit 1
fi

# Major version
CHECKMAJOR=$(grep -in "INTEGER,PARAMETER          :: MajorVersion" ${GLOBALS})
if [[ -z ${CHECKMAJOR} ]]; then
  echo "${RED}Could not 'INTEGER,PARAMETER          :: MajorVersion' in ${GLOBALS}${NC}"
  exit 1
else
  MAJORVERSIONSTR=$(printf '%-22s' "$MAJORVERSION")
  sed -i "s/INTEGER,PARAMETER          :: MajorVersion.*/INTEGER,PARAMETER          :: MajorVersion = ${MAJORVERSIONSTR}!> FileVersion number saved in each hdf5 file with hdf5 header/" ${GLOBALS}
fi

# Minor version
CHECKMINOR=$(grep -in "INTEGER,PARAMETER          :: MinorVersion" ${GLOBALS})
if [[ -z ${CHECKMINOR} ]]; then
  echo "${RED}Could not 'INTEGER,PARAMETER          :: MinorVersion' in ${GLOBALS}${NC}"
  exit 1
else
  MINORVERSIONSTR=$(printf '%-22s' "$MINORVERSION")
  sed -i "s/INTEGER,PARAMETER          :: MinorVersion.*/INTEGER,PARAMETER          :: MinorVersion = ${MINORVERSIONSTR}!> FileVersion number saved in each hdf5 file with hdf5 header/" ${GLOBALS}
fi

# Patch version
CHECKPATCH=$(grep -in "INTEGER,PARAMETER          :: PatchVersion" ${GLOBALS})
if [[ -z ${CHECKPATCH} ]]; then
  echo "${RED}Could not 'INTEGER,PARAMETER          :: PatchVersion' in ${GLOBALS}${NC}"
  exit 1
else
  PATCHVERSIONSTR=$(printf '%-22s' "$PATCHVERSION")
  sed -i "s/INTEGER,PARAMETER          :: PatchVersion.*/INTEGER,PARAMETER          :: PatchVersion = ${PATCHVERSIONSTR}!> FileVersion number saved in each hdf5 file with hdf5 header/" ${GLOBALS}
fi

# Complete version for AppImage container
CHECKWORKFLOW=$(grep -in "name: piclas-binaries-v" ${WORKFLOW})
if [[ -z ${CHECKWORKFLOW} ]]; then
  echo "${RED}Could not 'name: piclas-binaries-v' in ${WORKFLOW}${NC}"
  exit 1
else
  sed -i "s/.*name: piclas-binaries-v.*/        name: piclas-binaries-v${1}/" ${WORKFLOW}
fi
