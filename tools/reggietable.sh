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
  header=$(head -n1 ${DIRECTORY}readme.md 2> /dev/null || echo "")
  if [[ -z "$header" ]];then
    feature=""
    comparing=""
  else
    feature=$(head -n1 ${DIRECTORY}readme.md | cut -d "#" -f2)
    comparing=$(grep "\*\*Comparing\*\*" ${DIRECTORY}readme.md | cut -d ":" -f2)
  fi
  if [[ -z "${comparing}" ]]; then
    found=$(cat ${DIRECTORY}analyze.ini | grep "PartAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, PartAnalyze"
      else
        comparing="PartAnalyze"
      fi
    fi
    found=$(cat ${DIRECTORY}analyze.ini | grep "FieldAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, FieldAnalyze"
      else
        comparing="FieldAnalyze"
      fi
    fi
    found=$(cat ${DIRECTORY}analyze.ini | grep "SurfaceAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, SurfaceAnalyze"
      else
        comparing="SurfaceAnalyze"
      fi
    fi
  fi
  string="$(pwd -P)/${DIRECTORY}"
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
  header=$(head -n1 ${DIRECTORY}readme.md 2> /dev/null || echo "")
  if [[ -z "$header" ]];then
    echo "Missing readme.md under ${DIRECTORY}readme.md"
    comparing=""
  else
    feature=$(head -n1 ${DIRECTORY}readme.md | cut -d "#" -f2)
    comparing=$(grep "\*\*Comparing\*\*" ${DIRECTORY}readme.md | cut -d ":" -f2)
  fi
  if [[ -z "${comparing}" ]]; then
    found=$(cat ${DIRECTORY}analyze.ini | grep "PartAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, PartAnalyze"
      else
        comparing="PartAnalyze"
      fi
    fi
    found=$(cat ${DIRECTORY}analyze.ini | grep "FieldAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, FieldAnalyze"
      else
        comparing="FieldAnalyze"
      fi
    fi
    found=$(cat ${DIRECTORY}analyze.ini | grep "SurfaceAnalyze.csv")
    if [[ -n "${found}" ]]; then
      if [[ -n ${comparing} ]]; then
        comparing="${comparing}, SurfaceAnalyze"
      else
        comparing="SurfaceAnalyze"
      fi
    fi
  fi
  string="$(pwd -P)/${DIRECTORY}"
  link="[Link](regressioncheck${string#*regressioncheck}readme.md)"

  output="| %-4s | %-${nDIRECTORY}s |     | %-${nfeature}s | %-${nprocs}s | %-${ncomparing}s | %-${nlink}s |\n"
  printf "${output}" "${count}" "${DIRECTORY}" "${feature}" "${procs}" "${comparing}" "${link}"
done