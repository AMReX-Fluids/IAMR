
BOXLIB_HOME = /path/to/BoxLib

TOP = ../..

#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
DIM	  = 2
COMP	  = g++
USE_MPI   = FALSE

DEFINES += -DBL_PARALLEL_IO
EBASE = ViscBench

include $(TOP)/mk/Make.defs 
include ./Make.package

#include $(TOP)/amrlib/Make.package
include $(TOP)/bndrylib/Make.package
include $(TOP)/BoxLib/Make.package

INCLUDE_LOCATIONS += .
#INCLUDE_LOCATIONS += $(TOP)/amrlib
INCLUDE_LOCATIONS += $(TOP)/bndrylib
INCLUDE_LOCATIONS += $(TOP)/BoxLib
INCLUDE_LOCATIONS += $(TOP)/pAmrvis

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

vpath_cpp := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath_H   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath_FH  := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath_h   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath_F   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath_f   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/BoxLib
vpath %.cpp $(vpath_cpp)
vpath %.H   $(vpath_H)
vpath %.FH  $(vpath_FH)
vpath %.h   $(vpath_h)
vpath %.F   $(vpath_F)
vpath %.f   $(vpath_f)

all: $(executable)

#
# Create MS-IDE dsp file to build this project
#
DspOlevel = 1
DspFileRoot = $(EBASE)$(DIM)d
ifeq ($(DEBUG), TRUE)
  DspOlevel := 0
endif
dsp: $(DspFileRoot).dsp
.dspDepends: $(FEXE_sources) $(FEXE_headers) $(fEXE_sources) $(CEXE_headers) $(CEXE_sources) $(INCLUDE_LOCATIONS)
	@echo $^ > .dspDepends

$(DspFileRoot).dsp: .dspDepends
	@(echo Building DSP file = $(DspFileRoot).dsp;\
	$(TOP)/scripts/dsp.mak -p $(DspFileRoot) -t $(TOP) -d $(DIM) -o $(DspFileRoot).dsp -O $(DspOlevel) -f $^)

.PHONY: dsp .dspDepends libs cleanlibs nodata 

include $(TOP)/mk/Make.rules
