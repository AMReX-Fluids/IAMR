#
# $Id: GNUmakefile,v 1.29 2005-09-07 19:03:12 car Exp $
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
DEBUG	  = FALSE
DEBUG	  = TRUE
PROFILE   = FALSE
DIM       = 3
COMP	  = g++
COMP	  = Intel
FCOMP	  = Intel
USE_MPI   = TRUE
USE_MPI   = FALSE
USE_FLCTS = TRUE
USE_FLCTS = FALSE

#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE

USE_VAMPIRE = TRUE
USE_VAMPIRE = FALSE

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  #
  # What stencil do you want to use?
  #
  PRVERSION = v7
  EBASE = amr-serial
  DEFINES += -DBL_USE_HGPROJ_SERIAL
else
  EBASE = amr
endif

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(USE_ARRAYVIEW),TRUE)
  DEFINES += -DBL_USE_ARRAYVIEW
  DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

Bdirs := BoxLib iamrlib amrlib tensorMG mglib bndrylib 
ifeq ($(USE_FLCTS),TRUE)
  Bdirs += iamrlib/InflowForce
endif
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  Bdirs += hgproj-serial
else
  Bdirs += hgproj
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

ifeq ($(MACHINE),T3E)

  # We don't explicity set the path of the mpi directory on the T3E.

  $(executable): $(filter-out $(HYPRE_LIBS) -lmpi, $(LIBRARIES))
else
  $(executable): $(filter-out $(HYPRE_LIBS), $(LIBRARIES))
endif

include $(TOP)/mk/Make.rules
