
AMREX_HOME ?= ../../../../amrex

TOP = $(AMREX_HOME)

#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
DIM	  = 3
COMP	  = g++
USE_MPI   = FALSE

DEFINES += -DBL_PARALLEL_IO
EBASE = ViscBench

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package

vpathdir += ..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata


vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
