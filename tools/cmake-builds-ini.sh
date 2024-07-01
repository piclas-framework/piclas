#!/bin/bash
# ./tools/cmake-builds-ini.sh regressioncheck/WEK_PIC_maxwell/builds.ini

if [[ -f "$1" ]]; then
  padding='........................................'
  test=$(grep -v "EXCLUDE" "$1" | grep -vi "nocrosscombination" | grep -vE "^\s*\!" | grep -vi bInArY | sed s/'\s'//g | sed '/^[[:space:]]*$/d' | sed 's/\!.*//')
  cmakeopntions=$(while IFS= read -r line; do name="$(echo "$line" | cut -d "=" -f1)";value="$(echo "$line" | cut -d "=" -f2)";printf ' %s %s %s\n' "$name" "${padding:${#name}}" "$value"; done <<< "$test")
  echo -e "\n${cmakeopntions}"
  test=$(grep -v "EXCLUDE" "$1" | grep -vi "nocrosscombination" | grep -vE "^\s*\!" | grep -vi bInArY | sed s/'\s'//g | sed '/^[[:space:]]*$/d' | sed 's/\!.*//' | sed 's/,.*//')
  cmakeinput=$(while IFS= read -r line; do printf ' -D%s' "$line"; done <<< "$test")
  echo -e "\nSelect the first set of options via\n\ncmake .. ${cmakeinput}"
else
  printf 'File [%s] does not exist. Exit.\n' "$1"
fi
