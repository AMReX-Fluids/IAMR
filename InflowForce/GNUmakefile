TOP            = ../../../Parallel
PBOXLIB_HOME   = $(TOP)
PRECISION      = DOUBLE
DEBUG	       = TRUE
PROFILE        = FALSE
FCOMP          = Intel
DIM    	       = 3
BUILD_IN_PLACE = TRUE
EBASE          = turbMerge
USE_MPI        = FALSE
COMP           = g++

include $(TOP)/mk/Make.defs
include $(TOP)/BoxLib/Make.package

MySrcDirs = . $(TOP)/BoxLib

INCLUDE_LOCATIONS += $(MySrcDirs)

CEXE_sources += main.cpp
FEXE_sources += main_F.F
FEXE_headers += main_F.H

vpath %.cpp $(MySrcDirs)
vpath %.F   $(MySrcDirs)
vpath %.H   $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
vpath %.a   $(LIBRARY_LOCATIONS)

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

.PHONY: dsp .dspDepends

include $(TOP)/mk/Make.rules
