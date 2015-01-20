# set paths and compiler flags
include Makefile.defs

# List of executables to be built within the package
#PROGRAMS = boltzplatz

# "make" builds all
all: boltzplatz
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: ALL EXECUTABLES GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

boltzplatz: shared
	cd src && touch deplist.mk && $(MAKE) -j #create deplist.mk for builddebs to prevent errors if not existing

shared:
	cd share && $(MAKE) -j

shared_teclib:
	cd share && $(MAKE) -j teclib
	
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' SUCCESS: ALL TESTS GENERATED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

# utility targets
.PHONY: clean veryclean cleanshare

clean:
	cd src   && $(MAKE) clean
	rm -f lib/$(BOLTZPLATZ_LIB)

veryclean:
	cd src   && $(MAKE) veryclean
	rm -f lib/$(BOLTZPLATZ_LIB)
	rm -f *~ */*~ */*/*~

cleanshare:
	cd share && $(MAKE) clean

new: 
	make veryclean
	make
