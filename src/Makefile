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
ifneq ($(IMEX),)
  SRCFILE = particles/srcfiles_particles_imex.mk
endif
ifeq ($(EQNSYS),poisson)
  SRCFILE = particles/srcfiles_particles_hdg.mk
endif
ifeq ($(EQNSYS),magnetostatic)
  SRCFILE = particles/srcfiles_particles_hdg.mk
endif
ifeq ($(EQNSYS),magnetostatic_poisson)
  SRCFILE = particles/srcfiles_particles_hdg.mk
endif

SRC_RAW = $(shell cat $(SRCFILE))
SRC     = $(subst $$(EQNSYS),$(EQNSYS),$(SRC_RAW))

OBJ = $(SRC:.f90=.o)

all: preproc builddeps
	$(MAKE) $(BOLTZPLATZ_LIB)
	$(FC) $(FCFLAGS) $(INCDIRS) -c boltzplatz.f90 -o boltzplatz.o
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' BOLTZPLATZ COMPILED, NOW LINKING...'
	@echo '----------------------------------------------------------------------------------------'
	mkdir -p ../bin
	touch ../bin/userblock.txt
	rm ../bin/userblock.txt
	$(FC) $(FLFLAGS) boltzplatz.o ../lib/$(BOLTZPLATZ_LIB) -o ../bin/boltzplatz \
	$(LIBDIRS) \
	$(EXTERNAL_LIBS)
	@echo '----------------------------------------------------------------------------------------'
	@echo ' BOLTZPLATZ LINKED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	$(MAKE) userblock

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
	@python ../tools/builddeps --target=$(EQNSYS) --source=$(SRCFILE) deplist.mk
	@echo '----------------------------------------------------------------------------------------'
	@echo ' BUILD DEPLIST FINISHED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

$(BOLTZPLATZ_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(BOLTZPLATZ_LIB) $(OBJ)

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCDIRS) -c $< -o $@

userblock:
	@echo $(PWD)
	@echo '{[( BUILD_CMD )]}' > ../bin/userblock.txt
	@echo '$(FC) $(FCFLAGS) $(INCDIRS) -c boltzplatz.f90 -o boltzplatz.o' >> ../bin/userblock.txt
	@echo '{[( LINK_CMD )]}' >> ../bin/userblock.txt
	@echo '$(FC) $(FLFLAGS) boltzplatz.o ../lib/$(BOLTZPLATZ_LIB) -o ../bin/boltzplatz $(LIBDIRS) $(EXTERNAL_LIBS)' >> ../bin/userblock.txt
	@cd .. ;\
	 echo '{[( GIT REVISIONS )]}'        >> ./bin/userblock.txt ;\
	 git rev-parse --abbrev-ref HEAD     >> ./bin/userblock.txt ;\
	 git log --oneline -1                >> ./bin/userblock.txt ;\
	 echo 'master: '                     >> ./bin/userblock.txt ;\
	 git log origin/master --oneline -1  >> ./bin/userblock.txt ;\
	 echo '{[( GIT DIFF )]}'             >> ./bin/userblock.txt ;\
	 git diff -p origin/master..HEAD     >> ./bin/userblock.txt ;\
	 git diff -p                         >> ./bin/userblock.txt ;\
	 echo '{[( GIT URL )]}'       			 >> ./bin/userblock.txt ;\
   git config --get remote.origin.url	 >> ./bin/userblock.txt ;\
	 cd src

#------------------------------------------------------------------------------
# DEPENDENCIES
#------------------------------------------------------------------------------

include deplist.mk

#------------------------------------------------------------------------------
# UTILITY TARGETS
#------------------------------------------------------------------------------
.PHONY: clean veryclean

clean:
	rm -f $(OBJ) *.mod boltzplatz.o *.o

veryclean: clean
	rm -f *~ */*~
