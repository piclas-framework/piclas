# Developer Guide

The developer guide is intended to be for people who come in more close contact with PICLas, i.e., code developers and performance
analysts as well as people who are tasked with working or extending the documentation of PICLas.

```{toctree}
---
maxdepth: 1
caption: Table of Contents
numbered:
---
git_workflow.md
building_guide.md
styleguide.md
bestpractices.md
troubleshooting.md
code_extension.md
functions_and_subroutines.md
mpi.md
reggie.md
unittest.md
compiler.md
tools.md
performance.md
appimage.md
examples.md
```

This guide is organized to guide the first implementation steps as well as provide a complete overview of
the simulation code's features from a developer's point of view.

* The first Chapter {ref}`developerguide/git_workflow:GitLab Workflow` shall give an overview over the development workflow within
  the Gitlab environment, and the necessary steps to create a release, deploy updates to the Collaborative Numerics Group and GitHub.
* The second Chapter {ref}`developerguide/styleguide:Style Guide` describes the rules and guidelines regarding code development
  such as how the header of functions and subroutines look like.
* Chapter {ref}`developerguide/building_guide:Documentation` describes how to build the html/pdf files
  locally before committing changes to the repository and pointers on writing documentation.
* Chapter {ref}`developerguide/code_extension:Code Extension` describes how to extend standardized code blocks, e.g. add a new variable to the output.
* Chapter {ref}`developerguide/functions_and_subroutines:Useful Functions` contains a summary of useful functions
  and subroutines that might be re-used in the future.
* Chapter {ref}`developerguide/mpi:MPI Implementation` describes how PICLas subroutines and functions are parallelized.
* Chapter {ref}`developerguide/reggie:Regression Testing` summarizes the PICLas' continuous integration through regression testing.
* Chapter {ref}`developerguide/unittest:Unit Tests` shows which unit tests are used to test individual key components of the source code.
* Chapter {ref}`developerguide/compiler:Compiler Options` gives an overview of compiler options that are used in PICLas and their
  purpose.
* Chapter {ref}`developerguide/tools:Developer Tools` gives an overview over the tools and scripts for developers.
* Chapter {ref}`developerguide/performance:Performance Analysis` describes different tools that can be utilized for measuring the
  computational performance
* Chapter {ref}`developerguide/appimage:Building the AppImage Executable` described how an AppImage executable of HOPR is created.
* Chapter {ref}`developerguide/examples:Markdown Examples` gives a short overview of how to include code, equations, figures, tables
  etc. in the user and developer guides in Markdown.
