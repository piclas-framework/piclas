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
DEPEX = $(shell ls `cat depexceptions.mk`)

SRCPP = $(SRC:.f90=.i90)
OBJ   = $(SRC:.f90=.o)

all: preproc builddeps
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

preproc:
	@echo 'MODULE MOD_PreProcFlags' > preprocflags.f90
	@echo 'IMPLICIT NONE' >> preprocflags.f90
	@echo 'CHARACTER(LEN=1000) :: PREPROC_FLAGS = &' >> preprocflags.f90
	@echo \"$(FC) $(FCFLAGS) $(INCDIRS)\" | sed -r 's/(.{80})/\1"\/\/\&\n\"/g' >> preprocflags.f90
	@echo 'END MODULE MOD_PreProcFlags' >> preprocflags.f90

builddeps:
	$(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) boltzplatz.f90 $(PREOUT) boltzplatz.i90 ; \
	for file in $(SRC); \
	do \
	  outfile=$$(echo $$file | sed s/f90/i90/) ; \
	  $(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) $$file $(PREOUT) $$outfile ; \
	done
	@echo $(SRCPP) | sed -e 's/ /\n/g' > srcfiles_pp.mk
	@python ../tools/builddeps --source=srcfiles_pp.mk deplist.mk --exceptions=depexclude.mk

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
	rm -f $(OBJ) $(SRCPP) *.mod boltzplatz.o *.i

veryclean: clean
	rm -f *~ */*~
