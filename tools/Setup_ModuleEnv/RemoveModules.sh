#!/bin/bash -i

#==============================================================================
# title       : InstallPETSc.sh
# description : This script installs PETSc with specific setting in a
#               pre-installed module env
# date        : Apr 19, 2022
# version     : 1.0
# usage       : bash InstallPETSc.sh
# notes       : Bash in run interactively via "-i" to use "module load/purge"
#               commands
#==============================================================================

# Check privilege
if [[ "$EUID" -ne 0 ]]; then
  echo "Please run as root"
  exit 1
fi

# --------------------------------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------------------------------
dashes='............................................................................................'
whites='     '

if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    TERMCOLS=$(tput cols)
    BOLD="$(tput bold)"
    UNDERLINE="$(tput smul)"
    STANDOUT="$(tput smso)"
    NORMAL="$(tput sgr0)"
    NC="$(tput sgr0)"
    BLACK="$(tput setaf 0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
    BLUE="$(tput setaf 4)"
    MAGENTA="$(tput setaf 5)"
    CYAN="$(tput setaf 6)"
    WHITE="$(tput setaf 7)"
  fi
fi

# --------------------------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------------------------
read_char () {
if [ -n "$ZSH_VERSION" ]; then
  read -r -q "answer?"
else
  read -r -n 1 answer
fi
echo "$answer"
}

user_inquiry () { #echo "inquiry: [Y/n]"
printf "$*: [y/n]"
answer=$(read_char)
printf "\n"

if [[ ${answer} == 'y' ]]; then
  return 0
elif [[ ${answer} == 'n' ]]; then
  return 1
else
  return 1
fi
}

display_line () {
  if [[ ${1} == ${RED} ]]; then
    string='delete'
  else
    string=' keep '
  fi
  printf "[${1}${string}${NC}] "
  if [[ ${1} == ${RED} ]]; then
    printf "${1}${i}${NC}"
  else
    printf "$i"
  fi
  printf "\n"
}

remove_elements () { # Check list content: use * to prevent error: Arrays implicitly concatenate in [[ ]]. Use a loop (or explicit * instead of @).
  if [[ -z ${deletelist[*]} ]]; then
    return
  else
    for i in "${deletelist[@]}" ; do
      display_line ${RED}
    done
    echo ""
    user_inquiry "Remove all ${1} labelled [${RED}"delete"${NC}]?"
    if [[ $? -eq 0 ]]; then
      for i in "${deletelist[@]}" ; do
        if [[ "$i" == "/opt/"* ]]; then
          sudo rm -rf "$i"
          #echo "sudo rm -rf $i"
        fi
      done
      echo ""
    fi
  fi
}

# --------------------------------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------------------------------

FILES=$(find /opt/modules/modulefiles -name "*" -type f)
FIRSTCOMPILER=$(ls /opt/compiler/gcc/ | head -n1)
FIRSTCOMPILERGCC=$(echo gcc/$FIRSTCOMPILER)

echo "Non-compiler stuff"
for i in $FILES; do # Not recommended, will break on whitespace
  MODULEFILE=$(echo "$i")
  OPTDIR=$(grep -E "set.* topdir.* /opt" "$i" | sed 's/  */:/g' | cut -d ":" -f3)
  if [[ -n ${OPTDIR} ]]; then
    OPTSIZ=$(du -hs ${OPTDIR} 2>/dev/null | sed 's/\t/  /g' | cut -d " " -f1)
    if [[ -n ${OPTSIZ} ]]; then
      GCCFOUND=$(echo $OPTDIR | grep -i "gcc")
      if [[ -z ${GCCFOUND} ]]; then
        printf '%s %s [%s] %s %s\n' "$MODULEFILE" "${dashes:${#MODULEFILE}}" "$OPTSIZ" "${whites:${#OPTSIZ}}" "$OPTDIR"
      fi
    fi
  fi
done

printf "\n"
COMPILERS=$(ls /opt/compiler/gcc/)
for COMPILER in $COMPILERS; do # Not recommended, will break on whitespace
  if [[ -n ${COMPILER} ]]; then
    unset deletelist
    FIRSTCOMPILERGCC=$(echo gcc/$COMPILER)
    echo "Checking $FIRSTCOMPILERGCC"
    for i in $FILES; do # Not recommended, will break on whitespace
      MODULEFILE=$(echo "$i")
      if [[ "$MODULEFILE" == *"$FIRSTCOMPILERGCC"* ]]; then
        OPTDIR=$(grep -E "set.* topdir.* /opt" "$i" | sed 's/  */:/g' | cut -d ":" -f3)
        if [[ -n ${OPTDIR} ]]; then
          OPTSIZ=$(du -hs ${OPTDIR} 2>/dev/null | sed 's/\t/  /g' | cut -d " " -f1)
          if [[ -n ${OPTSIZ} ]]; then
            printf '%s %s [%s] %s %s\n' "$MODULEFILE" "${dashes:${#MODULEFILE}}" "$OPTSIZ" "${whites:${#OPTSIZ}}" "$OPTDIR"
            if [ ${#deletelist[@]} -eq 0 ]; then
              declare -a deletelist=("$MODULEFILE") # declare array for first entry
            else
              deletelist+=("$MODULEFILE")           # add element to array for following entries
            fi
              deletelist+=("$OPTDIR")           # add element to array for following entries
          else
            OPTSIZ="-"
            printf ${RED}'%s %s [%s] %s %s\n'${NC} "$MODULEFILE" "${dashes:${#MODULEFILE}}" "$OPTSIZ" "${whites:${#OPTSIZ}}" "Dangling module: Could not find the corresponding files in the /opt/... path!"
          fi
        fi
      fi
    done
    remove_elements "files/directories"
  fi
  printf "\n"
done