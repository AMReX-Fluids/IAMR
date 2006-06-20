#
# $Id: GNUmakefile,v 1.34 2006-06-20 16:42:36 car Exp $
#
PBOXLIB_HOME = ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = FALSE
PROFILE   = FALSE
DIM	  = 2
COMP	  = g++
FCOMP     = g77
COMP	  = mpiCC
FCOMP	  = g95
COMP	  = Intel
FCOMP	  = Intel
USE_MPI   = TRUE
USE_MPI   = FALSE

USE_FLCTS = FALSE
USE_FLCTS = TRUE

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

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

include $(TOP)/mk/Make.defs 
include ./Make.package

# CXX = mpiCC

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
ifeq ($(USE_FLCTS),TRUE)
  Bdirs += iamrlib/InflowForce
endif

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

# If We are using the f90 solvers then must define Flocs, and include

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  VPATH_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver

  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/
  Fdirs   := boxlib mg extern/SPARSKIT extern/LAPACK
  Flocs   := $(foreach dir, $(Fdirs), $(FBOXLIB_HOME)/$(dir))
  VPATH_LOCATIONS += $(Flocs)
endif

VPATH_LOCATIONS   += . $(Blocs)
INCLUDE_LOCATIONS += . $(Blocs)

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
