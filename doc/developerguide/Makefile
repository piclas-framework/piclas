# "make" builds all

# Switch to build all document parts - or only ready parts
BUILDALL=

SRC=000_developerguide.md 
SRC +=  080_guidelines.md 150_pic.md
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

