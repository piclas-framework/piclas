#!/usr/bin/env bash
# Strips one or more build directories down to only what's needed to keep as an artifact/cache:
#   - the "bin" folder and shared objects in "lib"
#   - with --keep-coverage: also the .gcno files throughout the build tree
#
# Usage:
#   ./clean_build_dir.sh [--keep-coverage] <build_directory> [<build_directory> ...]

set -euo pipefail

KEEP_COVERAGE=false
dirs=()

for arg in "$@"; do
  case "$arg" in
    --keep-coverage)
      KEEP_COVERAGE=true
      ;;
    *)
      dirs+=("$arg")
      ;;
  esac
done

if [ "${#dirs[@]}" -lt 1 ]; then
  echo "Usage: $0 [--keep-coverage] <build_directory> [<build_directory> ...]" >&2
  exit 1
fi

for dir in "${dirs[@]}"; do
  # Sanity check to make sure this is actually (probably) a PICLas build directory
  if [ ! -d "$dir/bin" ]; then
    echo "\"$dir\" has no \"bin\" folder, refusing to strip it."
    continue
  fi

  echo "Stripping build dir \"$dir\""
  if [ "$KEEP_COVERAGE" = true ]; then
    echo "  --keep-coverage set: preserving .gcno files"
  fi

  for entry in "$dir"/*; do
    name="$(basename "$entry")"

      # Keep bin and shared objects in lib
    case "$name" in
      bin)
        ;;
      lib)
        rm -f $entry/lib*.a
        ;;
      *)
        if [ "$KEEP_COVERAGE" = true ]; then
          if [ -d "$entry" ]; then
            # Delete everything except .gcno files, then prune any directories
            # left empty as a result (depth-first, so nested empties are caught
            # before their parents, and the entry itself is removed if it ends
            # up holding nothing but empty directories).
            find "$entry" -type f ! -name '*.gcno' -delete
            find "$entry" -depth -type d -empty -delete
          elif [ -f "$entry" ] && [[ "$name" != *.gcno ]]; then
            # Plain files (e.g. CMakeCache.txt, Makefile) can't be walked by
            # find, so remove them directly.
            rm -f -- "$entry"
          fi
          # Other entry types (symlinks, etc.) are left alone.
        else
          rm -rf -- "$entry"
        fi
        ;;
    esac
  done
done
