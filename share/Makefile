# set paths and compiler flags
MAINDIR = ../
include ../Makefile.defs

ifeq ($(EQNSYS),maxwell)
all: hdf5 recipes teclib  
else ifeq ($(EQNSYS),electrostatic)
all: hdf5 recipes teclib  
else ifeq ($(EQNSYS),electrostatic_pois)
all: hdf5 recipes teclib  
else ifeq ($(EQNSYS),maxwell_pois)
all: hdf5 recipes teclib  
else
all: hdf5 teclib recipes
endif

dir:
	@# Create the directory for the share programs 
	@ if [ ! -d $(SHARECOMP)/ ] ; then \
		mkdir $(SHARECOMP)/ ; \
	fi

teclib: dir
	@# Tecplot Binary library                   -> currently not in use
	@ if [ -n "$(TECLIB_VERSION)" ] && [ ! -f $(SHARECOMP)/tecio-$(TECLIB_VERSION)/tecio.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   T E C P L O T                                     =";\
		echo "===============================================================";\
		if [ ! -d $(SHARECOMP)/tecio-$(TECLIB_VERSION) ] ; then \
			cd $(SHARECOMP)/ && tar -xvzf ../tecio-$(TECLIB_VERSION).tar.gz 1>/dev/null && cd .. ; \
		fi ; \
		cd $(SHARECOMP)/tecio-$(TECLIB_VERSION) && ./Runmake linuxg27x64.24 -tecio && rm pltview pltview.o && ln -s tecio.a libtecio.a; \
		echo "===============================================================";\
		echo "= T E C P L O T   D O N E                                     =";\
		echo "===============================================================";\
	fi

hdf5: dir
	@# HDF5 libraries
	@ if [ -n "$(HDF5_VERSION)" ] && [ ! -f $(SHARECOMP)/hdf5-$(HDF5_VERSION)/hdf5/lib/libhdf5.a ] ; then \
		echo "===============================================================";\
		echo "= M A K E   H D F 5                                           =";\
		echo "===============================================================";\
		if [ ! -f hdf5-$(HDF5_VERSION).tar.bz2 ] ; then \
		  wget -N http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2 ; \
		  if [ "$?" != "0" ] ; then \
			  echo "wget failed. Trying normal releases..." && wget -N http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2 ; \
		  fi ; \
		fi ; \
		if [ ! -d $(SHARECOMP)/hdf5-$(HDF5_VERSION) ] ; then \
			cd $(SHARECOMP)/ && tar -xvjf ../hdf5-$(HDF5_VERSION).tar.bz2 1>/dev/null && cd .. ; \
		fi ; \
		cd $(SHARECOMP)/hdf5-$(HDF5_VERSION) && CC='$(CC)' FC='$(FC)' ./configure --enable-fortran --enable-fortran2003 --disable-shared && $(MAKE) && $(MAKE) install ; \
		echo "===============================================================";\
		echo "= H D F 5   D O N E                                           =";\
		echo "===============================================================";\
	fi

recipes: dir
	@# NUMERICAL RECIPES
	@ if [ ! -f $(MAINDIR)lib/librecipes.a ] || [ ! -f $(MAINDIR)src/nr.mod ] ; then \
		if [ ! -d $(SHARECOMP)/RECIPES ] ; then \
			cd $(SHARECOMP) && tar -xvzf ../recipes.tar.gz 1>/dev/null && cd .. ; \
		fi ; \
		cd  $(SHARECOMP)/RECIPES && $(MAKE) ; \
	fi		
	

# utility targets
.PHONY: clean
clean:  
	rm -f *~
	@# TECLIB library
	@ if [ -d $(SHARECOMP)/tecio-$(TECLIB_VERSION) ] ; then \
		cd $(SHARECOMP)/tecio-$(TECLIB_VERSION) && rm -f tecio.a && rm -f libtecio.a; \
	fi
	@# HDF5 libraries
	@ if [ -d $(SHARECOMP)/hdf5-$(HDF5_VERSION) ] ; then \
		cd $(SHARECOMP)/hdf5-$(HDF5_VERSION) && $(MAKE) clean && $(MAKE) uninstall ; \
	fi
	@# Numerical Recipes
	@ if [ -d $(SHARECOMP)/RECIPES ] ; then \
		cd $(SHARECOMP)/RECIPES && $(MAKE) clean ; \
	fi
	@ rm -f $(MAINDIR)lib/librecipes.a
