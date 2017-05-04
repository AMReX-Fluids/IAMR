
AMREX_HOME = /path/to/amrex

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
include $(TOP)/amrex/Make.package

INCLUDE_LOCATIONS += .
#INCLUDE_LOCATIONS += $(TOP)/amrlib
INCLUDE_LOCATIONS += $(TOP)/bndrylib
INCLUDE_LOCATIONS += $(TOP)/amrex
INCLUDE_LOCATIONS += $(TOP)/pAmrvis

vpath_cpp := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath_H   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath_FH  := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath_h   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath_F   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath_f   := . $(TOP)/pAmrvis $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/amrex
vpath %.cpp $(vpath_cpp)
vpath %.H   $(vpath_H)
vpath %.FH  $(vpath_FH)
vpath %.h   $(vpath_h)
vpath %.F   $(vpath_F)
vpath %.f   $(vpath_f)

all: $(executable)

include $(TOP)/mk/Make.rules
