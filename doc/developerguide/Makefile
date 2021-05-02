# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_developerguide.md 010_git_workflow.md 020_styleguide.md 025_functions_and_subroutines.md 030_mpi.md 040_reggie.md 050_compiler.md 060_tools.md 070_performance.md
# SRC +=  050_tutorials.md tutorials/*/*.md 070_unittest.md 080_develop_guide.md 099_parameter.md 
ifneq ($(BUILDALL),)
  SRC += missing.md
endif
SRC += references.md

all: developerguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

developerguide:
	pandoc $(SRC)   -o developerguide.pdf -s -N --filter pandoc-citeproc --from markdown+implicit_figures --listings --template=../template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N --filter pandoc-citeproc --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: developerguide

