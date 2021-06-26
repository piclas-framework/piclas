# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_userguide.md 010_installation.md 020_workflow.md 030_features_models.md 040_visu_output.md 050_tools.md
SRC +=  060_tutorials.md 080_cluster_guide.md 099_parameter.md
# SRC +=  tutorials/*/*.md 070_unittest.md 080_cluster_guide.md 099_parameter.md
ifneq ($(BUILDALL),)
  SRC += missing.md
endif
SRC += 100_appendix.md
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

# Copy header file to this directory
$(shell cp ../header.tex .)

# Get piclas version from globals_vars.f90
GLOBALSVARS := "../../src/globals/globals_vars.f90"
PICLAS_MAJOR_VERSION := $(shell grep "INTEGER.*PARAMETER.*MajorVersion.*\=" "$(GLOBALSVARS)" | cut -d "=" -f2 | cut -f1 -d! | sed 's/[[:space:]]//g')
PICLAS_MINOR_VERSION := $(shell grep "INTEGER.*PARAMETER.*MinorVersion.*\=" "$(GLOBALSVARS)" | cut -d "=" -f2 | cut -f1 -d! | sed 's/[[:space:]]//g')
PICLAS_PATCH_VERSION := $(shell grep "INTEGER.*PARAMETER.*PatchVersion.*\=" "$(GLOBALSVARS)" | cut -d "=" -f2 | cut -f1 -d! | sed 's/[[:space:]]//g')
PICLAS_VERSION := $(shell echo $(PICLAS_MAJOR_VERSION).$(PICLAS_MINOR_VERSION).$(PICLAS_PATCH_VERSION))

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
	@echo "Extracted PICLas version from ../../src/globals/globals_vars.f90: v$(PICLAS_VERSION)"
	sed -i 's/.*PICLas v.*/  PICLas v$(PICLAS_VERSION)/g' 000_userguide.md
	pandoc $(SRC)   -o userguide.pdf -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex
	rm header.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: userguide

