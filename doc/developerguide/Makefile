# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_developerguide.md 010_git_workflow.md 020_styleguide.md 025_functions_and_subroutines.md 030_mpi.md 040_reggie.md 050_compiler.md 060_tools.md 070_performance.md
# SRC +=  050_tutorials.md tutorials/*/*.md 070_unittest.md 080_develop_guide.md 099_parameter.md
ifneq ($(BUILDALL),)
  SRC += missing.md
endif
SRC += references.md

# GNU Make can only perform numerical comparison, split and compare by parts
PANDOCVERSION := $(shell pandoc --version | grep ^pandoc | sed 's/^.* //g')
PANDOCMAJOR   := $(shell echo $(PANDOCVERSION) | cut -f1 -d.)
PANDOCMINOR   := $(shell echo $(PANDOCVERSION) | cut -f2 -d.)

# Pandoc changed citeproc in 2.11
PANDOCFILTER:=--filter pandoc-citeproc
ifeq ($(shell expr $(PANDOCMAJOR) \= 2), 1)
ifeq ($(shell expr $(PANDOCMINOR) \>= 11), 1)
PANDOCFILTER:=--citeproc
endif
endif

# Update template to current version if possible
TEMPLATE_INFO := $(shell pandoc -D latex > /dev/null; echo $$?)
ifeq ($(shell expr $(TEMPLATE_INFO) \= 0), 1)
  $(info Found pandoc latex template, using system version)
  $(shell pandoc -D latex > ../template.tex)
else
  $(info pandoc latex template not found on system, using git version)
endif

all: developerguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

developerguide:
	pandoc $(SRC)   -o developerguide.pdf -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf           -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: developerguide

