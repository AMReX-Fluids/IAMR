
PBOXLIB_HOME = ../..

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
# Control Make Procedure
#
BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE


#
# Base name of the executable.
#
EBASE = nsSlabStat
CEXE_sources += $(EBASE).cpp

DEFINES += -DBL_PARALLEL_IO
include $(TOP)/mk/Make.defs ./Make.package

#
### Library and Include Search Paths ###
#######################################
ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/pBoxLib_2/Make.package

  INCLUDE_LOCATIONS += . 
  INCLUDE_LOCATIONS += $(TOP)/SlabStatTools
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2
else
  INCLUDE_LOCATIONS += . 
  INCLUDE_LOCATIONS += $(TOP)/include
  LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
  LIBRARIES         += -lbox$(DIM)d
endif

#
#
#
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
# Set Search Paths
#
ifeq ($(BUILD_INPLACE), TRUE)
  vpath %.H   . $(TOP)/SlabStatTools $(TOP)/pBoxLib_2
  vpath %.cpp . $(TOP)/SlabStatTools $(TOP)/pBoxLib_2
  vpath %.F   . $(TOP)/SlabStatTools $(TOP)/pBoxLib_2 $(TOP)/amrlib
else
  vpath %.H   . $(TOP)/SlabStatTools $(TOP)/include
  vpath %.cpp . $(TOP)/SlabStatTools  
  vpath %.F   . $(TOP)/SlabStatTools $(TOP)/amrlib
  vpath %.a   $(LIBRARY_LOCATIONS)
endif

#
### Make Rules ###
##################

all: $(executable)

ifneq ($(BUILD_INPLACE), TRUE)
#
# Build and install all libraries needed in an appropriate order.
#
libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install

#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
endif

include $(TOP)/mk/Make.rules
