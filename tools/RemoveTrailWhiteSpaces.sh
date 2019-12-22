#!/bin/bash
# removes trailing white space in all fortran files found in [gitroot/src] directory
removetrailwhitespaces() {
  SEARCHSTRING="\s\+$"
  REPLACESTRING=""
  FILESTRING=$1 #"*.f90"
  SEARCHPATH="$(git rev-parse --show-toplevel)/src"
  FOUNDFILES=$(find "$SEARCHPATH" -iname "$FILESTRING")
  ACCEPTEDFILES=""
  for REPLACEFILE in ${FOUNDFILES}; do
    if [[ -n $(grep -i "${SEARCHSTRING}" "${REPLACEFILE}") ]]; then
      echo found $(grep -io "${SEARCHSTRING}" ${REPLACEFILE} | wc -l) trailing white space occurences in $REPLACEFILE
      if [[ -z ${ACCEPTEDFILES} ]]; then
        ACCEPTEDFILES+="${REPLACEFILE}"
      else
        ACCEPTEDFILES+="\n${REPLACEFILE}"
      fi
    fi
  done
  if [[ -n ${ACCEPTEDFILES} ]]; then
    echo 'Do you want to remove trailing white spaces in these files?'
    echo 'proceed: [y/n]'
    read -es -n 1 PROCEED
    if [[ "${PROCEED}" == 'y' ]]; then
      for REPLACEFILE in ${FOUNDFILES}; do
        if [[ -n $(grep -i "${SEARCHSTRING}" "${REPLACEFILE}") ]]; then
          echo changes in file ${REPLACEFILE}
          sed -i 's/'${SEARCHSTRING}'/'${REPLACESTRING}'/gI' ${REPLACEFILE}
        fi
      done
    else
      echo "no replaces done"
    fi
  else
    echo "no ${FILESTRING} files with trailing white spaces found in gitroot/src"
  fi
}

removetrailwhitespaces "*.f90"
removetrailwhitespaces "*.h"
