TOP            = ../../../Parallel
PBOXLIB_HOME   = $(TOP)
PRECISION      = DOUBLE
DEBUG	       = TRUE
PROFILE        = FALSE
COMP           = Intel
FCOMP          = Intel
DIM    	       = 3
BUILD_IN_PLACE = TRUE
EBASE          = turbMerge
USE_MPI        = FALSE

include $(TOP)/mk/Make.defs
include $(TOP)/BoxLib/Make.package

MySrcDirs = . $(TOP)/BoxLib

INCLUDE_LOCATIONS += $(MySrcDirs)

CEXE_sources += main.cpp

vpath %.cpp $(MySrcDirs)
vpath %.F   $(MySrcDirs)
vpath %.H   $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)
vpath %.a   $(LIBRARY_LOCATIONS)

all: $(executable)

include $(TOP)/mk/Make.rules
