#!/bin/bash
# ./tools/cmake-builds-ini.sh regressioncheck/WEK_PIC_maxwell/builds.ini

# Check command line arguments
for ARG in "$@"; do
  if [ "${ARG}" == "--help" ] || [ ${ARG} == "-h" ]; then
    echo "Input arguments:"
    echo ""
    echo "  --help/-h            Print this help information."
    echo "  \$1                   Path to a reggie directory that contains a builds.ini file"
    echo ""
    echo "Usage example:"
    echo ""
    echo "  ~/piclas/tools/cmake-builds-ini.sh ~/piclas/regressioncheck/NIG_PIC_maxwell_RK4_p_adaption"
    echo ""
    exit 0
  fi
done
filepath="$1"
if [[ ! -f "$filepath" ]]; then
  filepath="$1"/builds.ini
fi

if [[ -f "$filepath" ]]; then
  echo "Testing $filepath"
  padding='........................................'
  test=$(grep -v "EXCLUDE" "$filepath" | grep -vi "nocrosscombination" | grep -vE "^\s*\!" | grep -vi bInArY | sed s/'\s'//g | sed '/^[[:space:]]*$/d' | sed 's/\!.*//')
  cmakeopntions=$(while IFS= read -r line; do name="$(echo "$line" | cut -d "=" -f1)";value="$(echo "$line" | cut -d "=" -f2)";printf ' %s %s %s\n' "$name" "${padding:${#name}}" "$value"; done <<< "$test")
  echo -e "\n${cmakeopntions}"
  test=$(grep -v "EXCLUDE" "$filepath" | grep -vi "nocrosscombination" | grep -vE "^\s*\!" | grep -vi bInArY | sed s/'\s'//g | sed '/^[[:space:]]*$/d' | sed 's/\!.*//' | sed 's/,.*//')
  cmakeinput=$(while IFS= read -r line; do printf ' -D%s' "$line"; done <<< "$test")
  echo -e "\nSelect the first set of options via\n\ncmake .. ${cmakeinput}"
else
  printf 'File/Folder [%s] does not exist.\nRun "./cmake-builds-ini.sh --help" for more information.\nExit.\n' "$1"
fi