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
  # print and read in keys and values
  padding='........................................'
  test=$(grep -v "EXCLUDE" "$filepath" | grep -vi "nocrosscombination" | grep -vE "^\s*\!" | grep -vi bInArY | sed s/'\s'//g | sed '/^[[:space:]]*$/d' | sed 's/\!.*//')
  cmakeopntions=$(while IFS= read -r line; do name="$(echo "$line" | cut -d "=" -f1)";value="$(echo "$line" | cut -d "=" -f2)";printf ' %s %s %s\n' "$name" "${padding:${#name}}" "$value"; done <<< "$test")
  echo -e "\n${cmakeopntions}"
  declare -a keys=()
  declare -A values_map
  while IFS= read -r line; do name="$(echo "$line" | cut -d "=" -f1)"; value="$(echo "$line" | cut -d "=" -f2)"; keys+=("$name"); values_map["$name"]="$value"; done <<< "$test"

  # read in EXCLUDE and nocross patterns
  exclude_patterns=$(grep -i "^EXCLUDE:" "$filepath" | sed 's/^EXCLUDE://i')
  # Improvement: also exclude nocrosscombination
  # nocross_patterns=$(grep -i "^nocrosscombination:" "$filepath" | sed 's/^nocrosscombination://i')

  should_exclude() {
    # read in combination
    local combo="$1"
    while IFS=',' read -ra conditions; do
      local all_conditions_met=1
      for condition in "${conditions[@]}"; do
        local var="${condition%%=*}"
        local val="${condition#*=}"
        # check if current condition (e.g. PICLAS_EQNSYSNAME=poisson) is not in the combination
        if ! [[ "$combo" == *"-D${var}=${val}"* ]]; then
          # if its not in the combination, build is allowed and other conditions must not be checked (e.g. PICLAS_EQNSYSNAME=maxwell, but exclude is for poisson and any other option)
          all_conditions_met=0
          break
        fi
      done
      # all conditions of current exclusion pattern were met, so build is excluded
      if [[ $all_conditions_met -eq 1 ]]; then
        return 0
      fi
    done <<< "$exclude_patterns"
    return 1
  }

  generate_combinations() {
    local prefix="$1"
    # use level to check how "deep" the combination is
    local -i level="$2"
    # combination contains all keys
    if [[ $level -eq ${#keys[@]} ]]; then
      if ! should_exclude "$prefix"; then
        echo "$prefix"
      fi
      return
    fi
    local key="${keys[$level]}"
    # key=$(echo "$key" | xargs)
    local value_string="${values_map[$key]}"
    IFS=',' read -ra values <<< "$value_string"
    # build combination recursively, creating a new 'path' for each value of current key and increase level
    for val in "${values[@]}"; do
      local new_prefix
      if [[ -z "$prefix" ]]; then
        # first key, no prefix
        new_prefix="-D${key}=${val}"
      else
        # not first key, add new key/val
        new_prefix="${prefix} -D${key}=${val}"
      fi
      generate_combinations "$new_prefix" $((level + 1))
    done
  }

  combinations=()
  while IFS= read -r combo; do
    combinations+=("$combo")
  done < <(generate_combinations "" 0)

  echo -e "\nAll possible combinations:"
  # Find max index width for proper alignment
  max_idx=${#combinations[@]}
  idx_width=${#max_idx}
  for i in "${!combinations[@]}"; do
    printf "[%${idx_width}d] cmake .. %s\n" "$((i+1))" "${combinations[$i]}"
  done

  echo -e "\nTotal number of combinations: ${#combinations[@]}"
  echo "Enter the combination:"
  read -r selection
  selection=$(echo "$selection" | tr ',' ' ')
  selected_numbers=($selection)
  for num in "${selected_numbers[@]}"; do
    if [[ $num =~ ^[0-9]+$ ]] && ((num > 0)) && ((num <= ${#combinations[@]})); then
      index=$((num-1))
      eval "cmake .. ${combinations[$index]} && make -j"
    else
      echo "Invalid selection: $num"
    fi
  done
else
  printf 'File [%s] does not exist. Exit.\n' "$1"
fi