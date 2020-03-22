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

all: userguide
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: DOCUMENTATION GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

userguide:
	pandoc $(SRC)   -o userguide.pdf -s -N --filter pandoc-citeproc --from markdown+implicit_figures --listings --template=../template.tex

# utility targets

tests:
	pandoc tests.md -o test.pdf      -s -N --filter pandoc-citeproc --from markdown+implicit_figures --listings --template=../template.tex

.PHONY: userguide

