#
# $Id: GNUmakefile,v 1.62 1998-11-20 21:34:19 lijewski Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT  = $(TOP)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = FALSE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 2
DIM	      = 3
WHICH_HG      =
COMP          = KCC
USE_WINDOWS   = FALSE
USE_MPI       = FALSE
USE_MPI       = TRUE
USE_NETCDF    = FALSE
USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE
#
# The only stencils in which we're interested when using USE_HGPROJ_SERIAL.
#
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
ifeq ($(DIM),2)
PRVERSION = v9
else
PRVERSION = v7
endif
endif
#
# Base name of the executable.
#
EBASE = amr

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(USE_HGPROJ_SERIAL),FALSE)
INCLUDE_LOCATIONS += . $(TOP)/include
LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
LIBRARIES         += -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d -lbox$(DIM)d
else
INCLUDE_LOCATIONS += . $(TOP)/hgproj-serial
INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
INCLUDE_LOCATIONS += $(TOP)/include
LIBRARY_LOCATIONS += $(TOP)/hgproj-serial/lib/$(machineSuffix) $(TOP)/lib/$(machineSuffix)
LIBRARIES         += -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d.$(PRVERSION) -lbox$(DIM)d
endif

ifeq ($(USE_MPI),TRUE)
LIBRARIES += -lmpi
ifeq ($(MACHINE),OSF1)
INCLUDE_LOCATIONS += /usr/local/mpi/include
LIBRARY_LOCATIONS += /usr/local/mpi/lib/alpha/ch_p4
endif
ifeq ($(MACHINE),AIX)
INCLUDE_LOCATIONS += /usr/lpp/ppe.poe/include
LIBRARY_LOCATIONS += /usr/lpp/ppe.poe/lib
endif
endif

ifeq ($(USE_WINDOWS),TRUE)
LIBRARIES += -lgraph
LIBRARIES += -lX11 
DEFINES   += -DBL_USE_WINDOWS
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES         += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif

ifeq ($(MACHINE),OSF1)
#
# Some additional stuff for our preferred development/debugging environment.
#
ifeq ($(PRECISION),DOUBLE)
FFLAGS += -real_size 64
endif
FDEBF += -C
FDEBF += -warn argument_checking
FDEBF += -warn declarations
ifneq ($(FC), f90)
FDEBF += -warn truncated_source
FDEBF += -warn unused
endif
endif
#
# For Running 3rd Only
#
3RD = 1
3RD =
ifdef 3RD
LDFLAGS += --link_command_prefix 3rd
LDFLAGS += -non_shared -v
endif

vpath %.a $(LIBRARY_LOCATIONS)

all: $(executable)

$(executable): $(LIBRARIES)
#
# Build and install all libraries needed by IAMRAll in an appropriate order.
#
libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/amrlib;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
	cd $(TOP)/hgproj;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) LBASE=proj EBASE= install
else
	cd $(TOP)/hgproj-serial; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) PRVERSION=$(PRVERSION) USE_MPI=$(USE_MPI) LBASE=proj EBASE=
endif
	cd $(TOP)/mglib;     $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(TOP)/mk/Make.rules
