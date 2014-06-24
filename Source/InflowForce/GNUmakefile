
PRECISION      = DOUBLE
DEBUG	       = TRUE
COMP           = g++
FCOMP          = gfortran
DIM    	       = 3
BUILD_IN_PLACE = TRUE
EBASE          = turbMerge
USE_MPI        = FALSE

BOXLIB_HOME ?= ../../../BoxLib

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

Bpack   := ./Make.package
Blocs   := .

include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package

MySrcDirs = . $(BOXLIB_HOME)/Src/C_BaseLib

INCLUDE_LOCATIONS += $(MySrcDirs)

CEXE_sources += main.cpp

vpath %.cpp $(MySrcDirs)
vpath %.F   $(MySrcDirs)
vpath %.H   $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules


