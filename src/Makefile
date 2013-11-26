#******************************************************************************
#
# B O L T Z P L A T Z    M A K E F I L E
#
#******************************************************************************

MAINDIR = ../
include ../Makefile.defs

SRCFILE = srcfiles.mk
ifneq ($(PARTICLES),)
  SRCFILE = particles/srcfiles_particles.mk
endif
ifeq ($(EQNSYS),electrostatic_pois)
  SRCFILE = particles/srcfiles_particles_pois.mk
endif
ifeq ($(EQNSYS),maxwell_pois)
  SRCFILE = particles/srcfiles_particles_pois.mk
endif
SRC_RAW = $(shell cat $(SRCFILE))
SRC     = $(subst $$(EQNSYS),$(EQNSYS),$(SRC_RAW))

OBJ = $(SRC:.f90=.o)

all: builddeps
	$(MAKE) $(BOLTZPLATZ_LIB)
	$(FC) $(FCFLAGS) $(INCDIRS) -c boltzplatz.f90 -o boltzplatz.o
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' BOLTZPLATZ COMPILED, NOW LINKING...'
	@echo '----------------------------------------------------------------------------------------'
	$(FC) $(FLFLAGS) boltzplatz.o ../lib/$(BOLTZPLATZ_LIB) -o ../bin/boltzplatz \
	$(LIBDIRS) \
	$(EXTERNAL_LIBS)
	@echo '----------------------------------------------------------------------------------------'
	@echo ' BOLTZPLATZ LINKED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

builddeps:
	@python ../tools/builddeps --target=$(EQNSYS) --source=$(SRCFILE) deplist.mk

$(BOLTZPLATZ_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(BOLTZPLATZ_LIB) $(OBJ)

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCDIRS) -c $< -o $@

#------------------------------------------------------------------------------
# DEPENDENCIES
#------------------------------------------------------------------------------

include deplist.mk

#------------------------------------------------------------------------------
# UTILITY TARGETS
#------------------------------------------------------------------------------
.PHONY: clean veryclean

clean:
	rm -f $(OBJ) *.mod boltzplatz.o *.i

veryclean: clean
	rm -f *~ */*~
