#
# $Id: GNUmakefile,v 1.52 1998-11-02 20:50:05 lijewski Exp $
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
# Base name of the executable.
#
EBASE = amr

PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)

include ../mk/Make.defs ./Make.package

INCLUDE_LOCATIONS += . $(TOP)/include

LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)

LIBRARIES += -lmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d -lbox$(DIM)d

ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
LIBRARIES += -lmpi
ifeq ($(MACHINE), OSF1)
INCLUDE_LOCATIONS += /usr/local/mpi/include
LIBRARY_LOCATIONS += /usr/local/mpi/lib/alpha/ch_p4
endif
ifeq ($(MACHINE), AIX)
INCLUDE_LOCATIONS += /usr/lpp/ppe.poe/include
LIBRARY_LOCATIONS += /usr/lpp/ppe.poe/lib
endif
endif

ifeq ($(USE_WINDOWS),TRUE)
LIBRARIES += -lgraph
LIBRARIES += -lX11 
DEFINES += -DBL_USE_WINDOWS
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
INCLUDE_LOCATIONS += .
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif

ifeq ($(MACHINE), OSF1)
#
# Some additional stuff for our preferred development/debugging environment.
#
ifeq ($(PRECISION), DOUBLE)
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

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include ../mk/Make.rules
