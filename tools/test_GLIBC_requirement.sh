#!/bin/bash

# ----------------------------------------------------------------------------
# Colors
# ----------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
BOLDGREEN='\033[1;32m'
YELLOW='\033[0;33m'
CYAN='\e[37;44m'
NC='\033[0m'  # No Color

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
# Script for checking the GLIBC versions required by the AppImage container
APPIMAGE=./piclasLeapfrogHDG

# check that the appimage is present
if ! [ -f ${APPIMAGE} ]; then
  echo -e "Could not find '${APPIMAGE}', which is required for running this script. Exit."
  exit
fi

# Extract Appimage
./piclasLeapfrogHDG --appimage-extract  # creates ./squashfs-root

# Check piclas binary
printf "\n${BOLDGREEN}piclas${NC}\n"
objdump -T squashfs-root/usr/bin/piclas | grep -o 'GLIBC_[0-9.]*' | sort -Vu

# Find GLIBC for all libs
printf "\n${BOLDGREEN}all libraries${NC}\n"
find squashfs-root -type f -name '*.so*' | while read f; do
  v=$(objdump -T "$f" 2>/dev/null | grep -o 'GLIBC_[0-9.]*' | sort -Vu | tail -1)
  [ -n "$v" ] && echo "$v  $f"
done | sort -V | tail