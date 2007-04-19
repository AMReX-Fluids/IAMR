#
# $Id: GNUmakefile,v 1.43 2007-04-19 22:20:51 aaspden Exp $
#
PBOXLIB_HOME = ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION   = DOUBLE
DEBUG       = FALSE
PROFILE     = FALSE
DIM         = 2
COMP        = Intel
FCOMP       = Intel
USE_MPI     = FALSE
USE_THREADS = FALSE
USE_XT3     = FALSE

# This allows GetForce to call a Fortran routine
DEFINES += -DGENGETFORCE

EBASE = amr

include $(TOP)/mk/Make.defs 
include ./Make.package

fincludes=$(includes)

Bdirs := BoxLib iamrlib amrlib mglib bndrylib tensorMG hgproj

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

ifeq ($(USE_XT3),TRUE)
  LIBRARIES += -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lpghpf
  DEFINES += -DBL_PGI -DMPICH_SKIP_MPICXX -DBL_XT3
ifeq ($(USE_MPI),TRUE)
  CXX := CC -target=catamount
  FC  := ftn -target=catamount
  fC  := ftn -target=catamount
else
  CXX := pgCC
  FC  := pgf90
  fC  := pgf90
endif
  CXXOPTF += -O -Minline=levels:10 --no_exceptions
  FOPTF   += -O -fastsse -Minline=levels:10
  fOPTF   += -O -fastsse -Minline=levels:10
endif

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
