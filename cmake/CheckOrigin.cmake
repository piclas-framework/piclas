# =========================================================================
# Git configuration
# =========================================================================
# Determine whether the source tree itself is the root of a git repository. Checking
# `git rev-parse` alone is not enough: git searches parent directories, so a manually
# downloaded source tree (e.g. an archive without a .git directory) placed inside another
# git repository would otherwise pick up that outer .git. Compare the repository top-level
# against the source dir and skip the git hook setup unless they match.
EXECUTE_PROCESS(COMMAND git rev-parse --show-toplevel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GIT_TOPLEVEL OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE GIT_COMMON_ERR ERROR_QUIET)
GET_FILENAME_COMPONENT(GIT_TOPLEVEL "${GIT_TOPLEVEL}" REALPATH)
GET_FILENAME_COMPONENT(SOURCE_REALPATH "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
IF(NOT GIT_COMMON_ERR EQUAL 0 OR NOT "${GIT_TOPLEVEL}" STREQUAL "${SOURCE_REALPATH}")
  MESSAGE(STATUS "Not a git repository, skipping git hook setup")
ELSE()
  # Check where the code originates
  EXECUTE_PROCESS(COMMAND git ls-remote --get-url WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GIT_ORIGIN OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Checking git origin: " ${GIT_ORIGIN})

  # Perform checks only if origin points to boltzplatz.eu, other origins can't commit
  IF("${GIT_ORIGIN}" MATCHES "piclas.boltzplatz.eu")
    # Install the tracked git hooks by pointing core.hooksPath at the version-controlled
    # .githooks directory. A relative path is resolved against each working tree's root, so
    # this works for regular clones and linked worktrees alike (each runs its own hooks) and
    # avoids a shared symlink into a single working tree that would dangle once it is removed.
    SET(HOOKS_DIR ".githooks")
    EXECUTE_PROCESS(COMMAND git config --get core.hooksPath OUTPUT_VARIABLE HOOKSPATH OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    IF ("${HOOKSPATH}" STREQUAL "" OR "${HOOKSPATH}" STREQUAL ".git/hooks")
      # Not configured (or still git's default): point it at the tracked hooks
      EXECUTE_PROCESS(COMMAND git config core.hooksPath ${HOOKS_DIR} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
      MESSAGE(STATUS "Set git core.hooksPath to ${HOOKS_DIR}")
    ELSEIF (NOT "${HOOKSPATH}" STREQUAL "${HOOKS_DIR}")
      # A different, user-defined hooks path is set: do not overwrite it, just warn
      MESSAGE (WARNING "Custom git hooks path '${HOOKSPATH}' detected. Please ensure to call ${HOOKS_DIR}/pre-commit manually.")
    ENDIF()
  ENDIF()
ENDIF()
