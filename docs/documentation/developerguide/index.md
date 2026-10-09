# Developer Guide

The developer guide is intended to be for people who come in more close contact with PICLas, i.e., code developers and performance
analysts as well as people who are tasked with working or extending the documentation of PICLas.

This guide is organized to guide the first implementation steps as well as provide a complete overview of
the simulation code's features from a developer's point of view.

* Chapter {ref}`developerguide/git_workflow:Development Workflow` describes the development process with issues, milestones and
  merge requests on GitLab.
* Chapter {ref}`developerguide/styleguide:Style Guide` describes the rules for writing code, such as the file header and the header
  of functions and subroutines.
* Chapter {ref}`developerguide/linting:Linting` describes the Fortran linter Fortitude and the Python linter Ruff, and how to fix
  their violations.
* Chapter {ref}`developerguide/code_extension:Code Extension` describes how to extend standardized code blocks, e.g. add a new
  variable to the output, and lists useful functions and subroutines.
* Chapter {ref}`developerguide/documentation:Writing Documentation` describes how to build the documentation locally and how to
  prepare figures.
* Chapter {ref}`developerguide/release:Release & Deployment` describes how to create a release, deploy it to GitHub and build the
  AppImage executable.
* Chapter {ref}`developerguide/mpi:MPI Implementation` gives the rules for MPI and shared memory windows and describes how PICLas is
  parallelized.
* Chapter {ref}`developerguide/reggie:Regression Testing` describes the continuous integration through regression testing.
* Chapter {ref}`developerguide/unittest:Unit tests` shows which unit tests check individual key components of the source code.
* Chapter {ref}`developerguide/troubleshooting:Troubleshooting` collects known problems from past debugging sessions, compiler
  options and the warnings check.
* Chapter {ref}`developerguide/performance:Performance Analysis` describes tools for measuring the computational performance.
