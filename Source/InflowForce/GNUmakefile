
PRECISION      = DOUBLE
DEBUG	       = TRUE
COMP           = g++
DIM    	       = 3
BUILD_IN_PLACE = TRUE
EBASE          = turbMerge
USE_MPI        = FALSE

AMREX_HOME ?= ../../../amrex

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Bpack   := ./Make.package
Blocs   := .

include $(AMREX_HOME)/Src/Base/Make.package

MySrcDirs = . $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(MySrcDirs)

CEXE_sources += main.cpp

vpath %.cpp $(MySrcDirs)
vpath %.F   $(MySrcDirs)
vpath %.H   $(MySrcDirs)
vpath %.h   $(MySrcDirs)
vpath %.f   $(MySrcDirs)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


