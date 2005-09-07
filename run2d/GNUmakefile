#
# $Id: GNUmakefile,v 1.30 2005-09-07 19:23:37 car Exp $
#
PBOXLIB_HOME = ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
PROFILE   = FALSE
DIM	  = 2
COMP	  = g++
FCOMP     = g77
COMP	  = Intel
FCOMP	  = Intel
USE_MPI   = TRUE
USE_MPI   = FALSE

USE_THREADS	= TRUE
USE_THREADS	= FALSE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE
#
# What stencil do you want to use?
#
PRVERSION = v5

ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  DEFINES += -DBL_USE_HGPROJ_SERIAL
  EBASE = amr-serial
else
  EBASE = amr
endif

ifeq ($(DIM),2)
  #
  # SYNCREG_2D.F builds RHS differently for 5 and 9 point stencils.
  #
  # Note that parallel hgproj currently only implements v5.
  #
  ifeq ($(findstring 9, $(PRVERSION)), 9)
    DEFINES += -DBL_PRVERSION=9
    # DEFINES += -DBL_PRVERSION=5
  endif
endif

include $(TOP)/mk/Make.defs 
include ./Make.package

fincludes=$(includes)

Bdirs := BoxLib iamrlib amrlib tensorMG mglib bndrylib
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  Bdirs += hgproj-serial
else
  Bdirs += hgproj
endif
ifeq ($(USE_HYPRE),TRUE)
  Bdirs   += HypreABec
  HYPRE_DIR=$(HOME)/amr/hypre
  INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
  LIBRARIES += -lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_blas -lHYPRE_utilities
  DEFINES += -DMG_USE_HYPRE
endif
ifeq ($(USE_METIS),TRUE)
  Bdirs += metis-4.0/Lib
  DEFINES += -DBL_USE_METIS
endif

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

# If We are using the f90 solvers then must define Flocs, and include

ifdef FBOXLIB_HOME
  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  DEFINES += -DMG_USE_FBOXLIB
  Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
endif

VPATH_LOCATIONS   += . $(Blocs) $(Flocs)
INCLUDE_LOCATIONS += . $(Blocs) $(Flocs)

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

$(executable): $(LIBRARIES)

include $(TOP)/mk/Make.rules
