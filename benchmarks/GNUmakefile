#
# $Id: GNUmakefile,v 1.3 1999-03-01 18:50:34 sstanley Exp $
#
PBOXLIB_HOME = ../../

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 3
COMP          = KCC
USE_MPI       = FALSE
#
# Base name of the executable.
#
EBASE = ViscBench

include $(TOP)/mk/Make.defs ./Make.package

CEXE_sources += $(EBASE)$(DIM)d.cpp

DEFINES += -DBL_PARALLEL_IO

#
# Local INCLUDES
#
INCLUDE_LOCATIONS += $(TOP)/pAmrvis
vpath %.H   . $(TOP)/pAmrvis
vpath %.cpp . $(TOP)/pAmrvis
vpath %.F   . $(TOP)/pAmrvis


INCLUDE_LOCATIONS += . .. $(TOP)/include
LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
LIBRARIES         += -lbox$(DIM)d -lbndry$(DIM)d

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

ifeq ($(USE_ARRAYVIEW),TRUE)
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
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

vpath %.H   .. 
vpath %.cpp .. 
vpath %.a   $(LIBRARY_LOCATIONS)

all: $(executable)

ifeq ($(MACHINE),T3E)
#
# We don't explicity set the path of the mpi directory on the T3E.
#
$(executable): $(filter-out -lmpi, $(LIBRARIES))
else
$(executable): $(LIBRARIES)
endif
#
# Build and install all libraries needed in an appropriate order.
#
libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean

include $(TOP)/mk/Make.rules
