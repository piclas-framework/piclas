#!/bin/bash
dots='..........................................................................................................................................................'

# Set colours
if test -t 1; then # if terminal
  NbrOfColors=$(which tput > /dev/null && tput colors) # supports color
  if test -n "$NbrOfColors" && test $NbrOfColors -ge 8; then
    NC="$(tput sgr0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
    YELLOW="$(tput setaf 3)"
  fi
fi

# Test for python 3
[[ ! -x "$(command -v python3)" ]] && echo "${RED}This tool requires python3 to run! Exit.${NC}" && exit 1

# Check command line arguments
for ARG in "$@"; do
  if [ "${ARG}" == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "Input arguments:"
    echo ""
    echo "  --help/-h            Print this help information."
    echo "  \$1                   Name of the .zip file, e.g., piclas-binaries-v3.4.0.zip or path to the file"
    echo "                       if it is located elsewhere. If the zip file is already extracted, supply"
    echo "                       the path to the binares, if they are in the current directory, supply a '.' (a dot)."
    echo "  \$2                   Path to the regressioncheck directory in the piclas repository, e.g.,"
    echo "                       ~/piclas/regressioncheck (where all the regression check directories are)."
    echo "                       Do not supply a specific reggie test case."
    echo "  \$3 (optional)        Path to the regressioncheck tool reggie.py, e.g., ~/reggie2.0/reggie.py"
    echo "                       and is only required if the reggie tool is not installed via pip."
    echo "                       Supply the path if unsure if reggie has been installed via pip or not."
    echo ""
    echo "Usage example:"
    echo ""
    echo "  ~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck [~/reggie2.0/reggie.py]"
    echo ""
    exit 0
  fi
done

if [[ -z "${2}" ]]; then
    echo -e "${RED}\$2 is empty. Supply the path to the regressioncheck directory in the second command line argument \$2, e.g., by running\n\n~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck${NC}" && exit 1
else
  if ! [[ -d "${2}" ]]; then
      echo -e "${RED}The path [${2}] does not exist. Supply the path to the regressioncheck directory in the second command line argument \$2, e.g., by running\n\n~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck${NC}" && exit 1
  else
    # Test for realpath command
    [[ ! -x "$(command -v realpath)" ]] && echo "${RED}This tool requires realpath to run! Exit.${NC}" && exit 1
    # Convert realtive path to absolute path
    REGRESSIONCHECKS=$(realpath ${2})
    # Test reggie directories
    [[ ! -d "${REGRESSIONCHECKS}/CHE_DSMC" ]]       && echo "${RED}CHE_DSMC not found under ${REGRESSIONCHECKS}.${NC}" && exit 1
    [[ ! -d "${REGRESSIONCHECKS}/CHE_poisson" ]]    && echo "${RED}CHE_poisson not found under ${REGRESSIONCHECKS}.${NC}" && exit 1
    [[ ! -d "${REGRESSIONCHECKS}/NIG_piclas2vtk" ]] && echo "${RED}NIG_piclas2vtk not found under ${REGRESSIONCHECKS}.${NC}" && exit 1
    [[ ! -d "${REGRESSIONCHECKS}/NIG_SuperB" ]]     && echo "${RED}NIG_SuperB not found under ${REGRESSIONCHECKS}.${NC}" && exit 1
    [[ ! -d "${REGRESSIONCHECKS}/CHE_BGK" ]]        && echo "${RED}CHE_BGK not found under ${REGRESSIONCHECKS}.${NC}" && exit 1
    [[ ! -d "${REGRESSIONCHECKS}/CHE_FPFlow" ]]     && echo "${RED}CHE_FPFlow not found under ${REGRESSIONCHECKS}.${NC}" && exit 1

    # Check $3
    if [[ -n "${3}" ]]; then
      [[ ! -f "${3}" ]] && echo -e "${RED}\$3 must point to reggie.py. Run again, e.g.,\n\n~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck ~/reggie2.0/reggie.py'${NC}" && exit 1
      reggieX () {
        python3 $1
      }

      # Run reggieX
      reggieX "$3 --help" > reggieX.log

      # Check if execution failed
      ERRORCODE=$?
      if [ ${ERRORCODE} -ne 0 ]; then
        echo " " && echo -e "${RED}Failed: [reggieX $3 --help > reggieX.log]${NC}" && echo -e "${RED}$(cat reggieX.log)${NC}"&& exit 1
      fi
    else
      # Test for reggie command
      [[ ! -x "$(command -v reggie)" ]] && echo -e "${RED}This tool requires 'reggie' as binary or the path to reggie.py supplied by \$3. Run again, e.g.,\n\n~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck ~/reggie2.0/reggie.py'${NC}" && exit 1
    fi
  fi
fi

# Check if dir is supplied or tar file
if [[ -d ${1} ]]; then
  DIR=${1}
  ZIPFILE=''

  # Get list of AppImages
  FILES=$(find ${DIR} -maxdepth 1 -type f -executable)
  [[ -z ${FILES} ]] && echo -e "${RED}error: no binaries found under ${DIR}${NC}" >&2 && exit 1
else
  # Check format of input argument
  re='^piclas.+\.zip$'
  # ^          asserts position at start of a line
  # piclas     matches the characters piclas literally (case sensitive)
  # .          matches any character (except for line terminators)
  # +          matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy)
  # \.         matches the character . with index 4610 (2E16 or 568) literally (case sensitive)
  # zip        matches the characters zip literally (case sensitive)
  # $          asserts position at the end of a line
  ! [[ "$(basename ${1})" =~ $re ]] && echo -e "${RED}error: Incorrect file specified; must be in the format of piclas-binaries-v3.4.0.zip or piclas-linux64.zip${NC}" >&2 && exit 1
  ZIPFILE=${1}
  DIR=''

  # Check if .zip file was correctly supplied
  if [ ! -f ${ZIPFILE} ]; then
    echo -e "${RED}no zip-file found for ${ZIPFILE}${NC}" && exit 1
  fi

  # Extract
  unzip "$ZIPFILE" | tee unzip.log

  # Check if extraction failed
  ERRORCODE=$?
  if [ ${ERRORCODE} -ne 0 ]; then
    echo " " && echo -e "${RED}Failed: [unzip ${ZIPFILE} | tee unzip.log]${NC}" && exit 1
  fi

  # Get list of AppImages
  FILES=$(grep inflating unzip.log | cut -d ":" -f2 | grep -v ".txt")
  [[ -z ${FILES} ]] && echo -e "${RED}error: no binaries found in unzip.log${NC}" >&2 && exit 1
fi

# Loop over files. Important that no "" is used in the for/in satement to
# access the elements (binaries) one by one and not everything as a single string
GLOBALERROR=0
for FILE in ${FILES}; do
  # Just to be sure, make executable
  chmod +x "$FILE"

  # Run
  s="./${FILE} --help > ${FILE}.log"
  printf 'Running [%s] %s ' "$s" "${dots:${#s}}"
  ./${FILE} --help > ${FILE}.log
  ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC}"; else echo "${GREEN}OK${NC}"; fi;

  # Run reggie
  FOF='-p'
  FOF=''
  if [[ -n "${3}" ]]; then
    if [[ $(basename ${FILE}) == "piclasDSMC" ]]       ; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_DSMC/BC_PorousBC";                             printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_DSMC/BC_PorousBC"                             > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasLeapfrogHDG" ]]; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_poisson/SurfFlux_ThermionicEmission_Schottky"; printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_poisson/SurfFlux_ThermionicEmission_Schottky" > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasDSMC" ]]       ; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_piclas2vtk";                                   printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_piclas2vtk"                                > piclas2vtk.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "superB" ]]           ; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_SuperB";                                       printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_SuperB"                                       > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasBGK" ]]        ; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_BGK/2D_VTS_Insert_CellLocal";                  printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_BGK/2D_VTS_Insert_CellLocal"                  > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasFP" ]]         ; then s="python3 $3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_FPFlow/2D_VTS_Insert_CellLocal";               printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggieX "$3 $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_FPFlow/2D_VTS_Insert_CellLocal"               > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
  else
    if [[ $(basename ${FILE}) == "piclasDSMC" ]]       ; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_DSMC/BC_PorousBC";                                 printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/CHE_DSMC/BC_PorousBC"                             > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasLeapfrogHDG" ]]; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_poisson/SurfFlux_ThermionicEmission_Schottky";     printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/CHE_poisson/SurfFlux_ThermionicEmission_Schottky" > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasDSMC" ]]       ; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_piclas2vtk";                                       printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/NIG_piclas2vtk"                                   > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "superB" ]]           ; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/NIG_SuperB";                                           printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/NIG_SuperB"                                       > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasBGK" ]]        ; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_BGK/2D_VTS_Insert_CellLocal";                      printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/CHE_BGK/2D_VTS_Insert_CellLocal"                  > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
    if [[ $(basename ${FILE}) == "piclasFP" ]]         ; then s="reggie $FOF -e ${FILE} ${REGRESSIONCHECKS}/CHE_FPFlow/2D_VTS_Insert_CellLocal";                   printf 'Running [%s] %s ' "$s" "${dots:${#s}}"; reggie $FOF -e "${FILE}" "${REGRESSIONCHECKS}/CHE_FPFlow/2D_VTS_Insert_CellLocal"               > ${FILE}.log; ERRORCODE=$?; if [[ ${ERRORCODE} -ne 0  ]]; then echo "${RED}FAIL${NC} (see ${FILE}.log)"; else echo "${GREEN}OK${NC}"; fi; fi;
  fi

  [[ ${ERRORCODE} -ne 0 ]] && GLOBALERROR=1
done

if [ ${GLOBALERROR} -ne 0 ]; then
  echo " " && echo -e "${RED}Failed one or more tests.${NC}" && exit 1
fi