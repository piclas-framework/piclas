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
SRCS  = $(SRC:.f90=.s)

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
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' CHECK PREPROCESSOR ...'
	@echo '----------------------------------------------------------------------------------------'
	@echo '----------------------------------------------------------------------------------------'
	@echo 'MODULE MOD_PreProcFlags' > preprocflags_tmp.f90 ; \
	echo 'IMPLICIT NONE' >> preprocflags_tmp.f90 ; \
	echo 'CHARACTER(LEN=1000) :: PREPROC_FLAGS = &' >> preprocflags_tmp.f90 ; \
	echo \"$(FC) $(FCFLAGS) $(INCDIRS)\" | sed -r 's/(.{120})/\1"\/\/\&\n\"/g' >> preprocflags_tmp.f90 ; \
	echo 'END MODULE MOD_PreProcFlags' >> preprocflags_tmp.f90 ; \
	if [ ! -f preprocflags.f90 ] ; then \
	  touch preprocflags.f90 ; \
	fi ; \
	if  cmp -s preprocflags.f90 preprocflags_tmp.f90   ; then  \
	  echo '   -> NO CHANGE IN PREPROCESSOR DIRECTIVES' ; \
	else \
	  echo '   -> PREPROCESSOR DIRECTIVES CHANGED (GENERATE NEW DEPLIST)'  ; \
	  rm -f srcfiles_preproc.mk ; \
	fi ; \
	mv preprocflags_tmp.f90 preprocflags.f90
	@echo '----------------------------------------------------------------------------------------'
	@echo ' CHECK PREPROCESSOR FINISHED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

builddeps:
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' BUILD DEPLIST...'
	@echo '----------------------------------------------------------------------------------------'
	@ if [ ! -f srcfiles_preproc.mk ] || [ ! -f deplist.mk ] ; then \
	  if echo "$(PREPROC)" | grep -q "DCRAY" ; then \
	    $(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) boltzplatz.f90 ; \
	  else \
	    $(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) boltzplatz.f90 $(PREOUT) boltzplatz.i90 ; \
	  fi ; \
	  echo '   APPLY PREPROCESSOR DIRECTIVES...' ; \
	  for file in $(SRC) ; \
	  do \
	    outfile=$$(echo $$file | sed s/f90/i90/) ; \
	    if echo "$(PREPROC)" | grep -q "DCRAY" ; then \
	      $(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) $$file ; \
        locfile=$$(echo $$(basename $$file) | sed s/f90/i/) ; \
	  		mv $$locfile $$outfile ; \
	  	else \
	      $(FC) $(FPRE) $(PREPROC) -I . $(INCDIRS) $$file $(PREOUT) $$outfile ; \
	  	fi ; \
	  done ; \
	  echo '   CREATE DEPLIST...' ; \
	  echo $(SRCPP) | sed -e 's/ /\n/g' > srcfiles_preproc.mk ; \
	  python ../tools/builddeps --source=srcfiles_preproc.mk deplist.mk --exceptions=depexclude.mk ; \
	else \
	  echo '   DEPLIST IS LEFT UNCHANGED' ; \
	fi
	@echo '----------------------------------------------------------------------------------------'
	@echo ' BUILD DEPLIST FINISHED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

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
	rm -f $(OBJ) $(SRCS) *.mod boltzplatz.o boltzplatz.i90
	rm -f $(SRCPP) boltzplatz.i90

veryclean: clean
	rm deplist.mk srcfiles_preproc.mk
	rm -f $(SRCPP) boltzplatz.i90
	rm -f *~ */*~ */*/*~
