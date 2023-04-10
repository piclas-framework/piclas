#!/bin/bash

# Script for checking the file size of every changed/added/removed file compared with master.dev branch
# Depending on the file size, different colours indicate if the added files are too large

if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    NC="$(tput sgr0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
  fi
fi

# Get all changed files
CHANGED=$(git diff --name-only master.dev)

# Check if any changes are present
if [[ -n $CHANGED ]]; then
  # Sort found files by size (-S) in reverse ordering (-r)
  SORTED=$(ls -Shsr $CHANGED 2> /dev/null)
  # Loop over all changes
  for file in $SORTED; do
    # Check if path is a file that exists
    if [[ -f $file ]]; then
      # -b, --bytes           equivalent to '--apparent-size --block-size=1'
      LINE=$(du -h $file)
      FILESIZE=$(du -b $file | cut -d '	' -f1) # this is a tab, not a white space
      if [[ $FILESIZE -gt 1000000 ]]; then
        printf "${RED}$LINE${NC}\n"
      else
        if [[ $FILESIZE -gt 100000 ]]; then
          printf "${YELLOW}$LINE${NC}\n"
        else
          printf "${GREEN}$LINE${NC}\n"
        fi
      fi
    fi
  done
  echo ""
  FILESIZE=$(du -h $file | cut -d '	' -f1) # this is a tab, not a white space
  printf "The largest file is [$FILESIZE]\n"
else
  echo "no changes to master.dev found. Exit."
fi


