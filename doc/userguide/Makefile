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

# Pandoc changed citeproc in 2.11
PANDOCVERSION := $(shell pandoc --version | grep ^pandoc | sed 's/^.* //g')
PANDOCMAJOR   := $(shell echo $(PANDOCVERSION) | cut -f1 -d.)
PANDOCMINOR   := $(shell echo $(PANDOCVERSION) | cut -f2 -d.)
ifeq ($(shell expr $(PANDOCMAJOR) \= 2), 0)
  $(error No or unkown version of pandoc. Please use pandoc v2.x to build this documentation.)
endif

ifeq ($(shell expr $(PANDOCMINOR) \>= 11), 1)
  PANDOCFILTER=--citeproc
else
  PANDOCFILTER=--filter pandoc-citeproc
endif

# Update template to current version if possible
TEMPLATE_INFO := $(shell pandoc -D latex > /dev/null; echo $$?)
ifeq ($(shell expr $(TEMPLATE_INFO) \= 0), 1)
  $(info Found pandoc latex template, using system version)
  $(shell pandoc -D latex > template.tex)
else
  $(info pandoc latex template not found on system, using git version)
endif

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
	pandoc $(SRC)   -o userguide.pdf -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N $(PANDOCFILTER) --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: userguide

