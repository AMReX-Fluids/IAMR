
AMREX_HOME ?= ../../../amrex

TOP = $(AMREX_HOME)

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

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
#include $(AMREX_HOME)/Tools/C_util/Make.package

vpathdir += ..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
vpathdir += $(AMREX_HOME)/Src/Boundary

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
vpathdir += $(AMREX_HOME)/Src/Amr

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
vpathdir += $(AMREX_HOME)/Tools/C_util


vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
