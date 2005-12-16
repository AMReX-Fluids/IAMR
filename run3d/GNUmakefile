#
# $Id: GNUmakefile,v 1.32 2005-12-16 23:08:18 lijewski Exp $
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
COMP	  = blrts_xlC
FCOMP	  = blrts_xlf
USE_MPI   = FALSE
USE_MPI   = TRUE
USE_FLCTS = TRUE
USE_FLCTS = FALSE

USE_HYPRE = FALSE
USE_METIS = FALSE
#
# A "trick" to make the BGL frontend compile CGL code.
#
USE_BGL = FALSE

#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE

USE_VAMPIRE = TRUE
USE_VAMPIRE = FALSE

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

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

ifeq ($(USE_BGL),TRUE)
  NEEDS_FLUSH_F := TRUE
endif

include $(TOP)/mk/Make.defs

ifeq ($(USE_BGL),TRUE)
  MACHINE=BGL
  CXX=blrts_xlC
  FC=blrts_xlf
  fC=blrts_xlf
  CPPFLAGS+=-I/bgl/BlueLight/ppcfloor/bglsys/include
  CPPFLAGS+= -DMPICH_SKIP_MPICXX
  LDFLAGS+=-L/bgl/BlueLight/ppcfloor/bglsys/lib
  LDFLAGS+=-L/opt/ibmcmp/xlf/9.1/blrts_lib
  BL_MPI_LIBS=-lxlf90 -lxlfmath -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts
  FORTLINK:=LOWERCASE
endif

include ./Make.package

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

include $(TOP)/mk/Make.rules
