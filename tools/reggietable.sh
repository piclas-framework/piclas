#!/bin/bash

# Script for automatically creating the reggie table for the content of the current folder that can be copied into REGGIE.md

if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    NC="$(tput sgr0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
  fi
fi

nDIRECTORY=0
nprocs=0
nfeature=0
ncomparing=0
nlink=0

# find max length of each element
for DIRECTORY in $(ls -d */) ; do
  procs=$(grep -in MPI ${DIRECTORY}command_line.ini | cut -d "=" -f2)
  feature=$(head -n1 ${DIRECTORY}readme.md | cut -d "#" -f2)
  comparing=$(grep "**Comparing**" ${DIRECTORY}readme.md | cut -d ":" -f2)
  string="$(pwd -P)${DIRECTORY}"
  link="[Link](regressioncheck${string#*regressioncheck}readme.md)"

  nDIRECTORY=$(( nDIRECTORY > ${#DIRECTORY} ? nDIRECTORY : ${#DIRECTORY} ))
  nprocs=$(( nprocs > ${#procs} ? nprocs : ${#procs} ))
  nfeature=$(( nfeature > ${#feature} ? nfeature : ${#feature} ))
  ncomparing=$(( ncomparing > ${#comparing} ? ncomparing : ${#comparing} ))
  nlink=$(( nlink > ${#link} ? nlink : ${#link} ))
done

count=0
# Display the lines
for DIRECTORY in $(ls -d */) ; do
  count=$((count + 1)) # must be placed after echo on HLRS, ForHLR1 (but not on local PC with zsh)
  procs=$(grep -in MPI ${DIRECTORY}command_line.ini | cut -d "=" -f2)
  feature=$(head -n1 ${DIRECTORY}readme.md | cut -d "#" -f2)
  comparing=$(grep "**Comparing**" ${DIRECTORY}readme.md | cut -d ":" -f2)
  string="$(pwd -P)${DIRECTORY}"
  link="[Link](regressioncheck${string#*regressioncheck}readme.md)"

  output="| %-4s | %-${nDIRECTORY}s |     | %-${nfeature}s | nProcs=%-${nprocs}s  | %-${ncomparing}s | %-${nlink}s |\n"
  printf "${output}" "${count}" "${DIRECTORY}" "${feature}" "${procs}" "${comparing}" "${link}"
done

