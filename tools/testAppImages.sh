#!/bin/bash
#==================================================================================================
# testAppImages.sh: Run each PICLas AppImage with --help and selected reggie regression checks
#==================================================================================================
set -o pipefail

dots='.................................................................................................................................................................................................'

# Set colours
if [[ -t 1 ]] && command -v tput > /dev/null; then # if terminal and tput available
  NbrOfColors=$(tput colors 2> /dev/null)
  if [[ -n "${NbrOfColors}" ]] && (( NbrOfColors >= 8 )); then
    NC="$(tput sgr0)"
    RED="$(tput setaf 1)"
    GREEN="$(tput setaf 2)"
  fi
fi

die() { echo -e "${RED}$*${NC}" >&2; exit 1; }

# Regression checks per binary (paths relative to the regressioncheck directory $2)
# Note that NIG_piclas2vtk is not run with the piclas2vtk executable but with the piclas executable
# (built for DSMC in this case)
declare -A TESTS=(
  [piclasDSMC]="CHE_DSMC/BC_PorousBC NIG_piclas2vtk"
  [piclasLeapfrogHDG]="CHE_poisson/SurfFlux_ThermionicEmission_Schottky"
  [superB]="NIG_SuperB"
  [piclasBGK]="CHE_BGK/2D_VTS_Insert_CellLocal"
  [piclasFP]="CHE_FPFlow/2D_VTS_Insert_CellLocal"
)

usage() {
  cat << EOF
Input arguments:

  --help/-h            Print this help information.
  \$1                   Name of the .zip file, e.g., piclas-binaries-v3.4.0.zip or path to the file
                       if it is located elsewhere. If the zip file is already extracted, supply
                       the path to the binaries, if they are in the current directory, supply a '.' (a dot).
  \$2                   Path to the regressioncheck directory in the piclas repository, e.g.,
                       ~/piclas/regressioncheck (where all the regression check directories are).
                       Do not supply a specific reggie test case.
  \$3 (optional)        Path to the regressioncheck tool reggie.py, e.g., ~/reggie/reggie.py
                       and is only required if the reggie tool is not installed via pip.
                       Supply the path if unsure if reggie has been installed via pip or not.

Usage example:

  ~/piclas/tools/testAppImages.sh piclas-binaries-v3.4.0.zip ~/piclas/regressioncheck [~/reggie/reggie.py]

EOF
}

EXAMPLE="Run, e.g.,\n\n~/piclas/tools/testAppImages.sh piclas-binaries-v1.2.3.zip ~/piclas/regressioncheck [~/reggie/reggie.py]"

#==================================================================================================
# Check command line arguments and requirements
#==================================================================================================
for ARG in "$@"; do
  if [[ "${ARG}" == "--help" || "${ARG}" == "-h" ]]; then usage; exit 0; fi
done

command -v python3  > /dev/null || die "This tool requires python3 to run! Exit."
command -v realpath > /dev/null || die "This tool requires realpath to run! Exit."

[[ -n "${1:-}" ]] || die "\$1 is empty. Supply the .zip file or the directory containing the binaries. ${EXAMPLE}"
[[ -n "${2:-}" ]] || die "\$2 is empty. Supply the path to the regressioncheck directory. ${EXAMPLE}"
[[ -d "${2}"   ]] || die "The path [${2}] does not exist. Supply the path to the regressioncheck directory. ${EXAMPLE}"

# Convert relative path to absolute path
REGRESSIONCHECKS=$(realpath "${2}")

# Test that all required regression checks exist
for t in ${TESTS[@]}; do
  [[ -d "${REGRESSIONCHECKS}/${t}" ]] || die "${t} not found under ${REGRESSIONCHECKS}."
done

# Select reggie (pip-installed binary or reggie.py supplied via $3)
if [[ -n "${3:-}" ]]; then
  [[ -f "${3}" ]] || die "\$3 must point to reggie.py. ${EXAMPLE}"
  REGGIE=(python3 "${3}")
else
  command -v reggie > /dev/null || die "This tool requires 'reggie' as binary or the path to reggie.py supplied by \$3. ${EXAMPLE}"
  REGGIE=(reggie)
fi
if ! "${REGGIE[@]}" --help > reggie.log 2>&1; then
  echo -e "${RED}Failed: [${REGGIE[*]} --help > reggie.log]${NC}" >&2
  cat reggie.log >&2
  exit 1
fi

# reggie runs every case via mpirun, even with a single rank
command -v mpirun > /dev/null || die "error: mpirun not on PATH (load an MPI module, and activate any venv BEFORE loading modules)."

#==================================================================================================
# Get list of binaries: either from a directory or by extracting a .zip file
#==================================================================================================
if [[ -d "${1}" ]]; then
  mapfile -t FILES < <(find "${1}" -maxdepth 1 -type f -executable | sort)
  (( ${#FILES[@]} )) || die "error: no binaries found under ${1}"
else
  # Must be in the format piclas*.zip
  re='^piclas.+\.zip$'
  # ^          asserts position at start of a line
  # piclas     matches the characters piclas literally (case sensitive)
  # .          matches any character (except for line terminators)
  # +          matches the previous token between one and unlimited times, as many times as possible, giving back as needed (greedy)
  # \.         matches the character . with index 4610 (2E16 or 568) literally (case sensitive)
  # zip        matches the characters zip literally (case sensitive)
  # $          asserts position at the end of a line
  [[ "$(basename "${1}")" =~ ${re} ]] || die "error: Incorrect file specified; must be in the format of piclas-binaries-v3.4.0.zip or piclas-linux64.zip"
  ZIPFILE="${1}"
  [[ -f "${ZIPFILE}" ]] || die "no zip-file found for ${ZIPFILE}"

  # Extract (-o: overwrite files without prompting)
  unzip -o "${ZIPFILE}" | tee unzip.log || die "Failed: [unzip -o ${ZIPFILE} | tee unzip.log]"

  # Binaries only: skip .txt files and directory entries
  mapfile -t FILES < <(unzip -Z1 "${ZIPFILE}" | grep -v -e '\.txt$' -e '/$')
  (( ${#FILES[@]} )) || die "error: no binaries found in ${ZIPFILE}"
fi

#==================================================================================================
# Run tests
#==================================================================================================
GLOBALERROR=0
run() { # usage: run <logfile> <cmd...>
  local log=$1; shift
  local s="$* > ${log}"
  printf 'Running [%s] %s ' "${s}" "${dots:${#s}}"
  if "$@" > "${log}" 2>&1; then
    echo -e "${GREEN}OK${NC}"
  else
    echo -e "${RED}FAIL${NC} (see ${log})"
    GLOBALERROR=1
  fi
}

for FILE in "${FILES[@]}"; do
  [[ "${FILE}" == */* ]] || FILE="./${FILE}"
  chmod +x "${FILE}"
  name=$(basename "${FILE}")
  run "${name}.help.log" "${FILE}" --help
  for t in ${TESTS[${name}]:-}; do
    run "${name}_${t//\//_}.log" "${REGGIE[@]}" -e "${FILE}" "${REGRESSIONCHECKS}/${t}"
  done
done

if (( GLOBALERROR )); then
  echo -e "\n${RED}Failed one or more tests.${NC}" >&2
  exit 1
fi
echo -e "\n${GREEN}All tests passed.${NC}"
