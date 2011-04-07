#
# $Id: GNUmakefile,v 1.50 2011-04-07 21:21:53 lijewski Exp $
#
PBOXLIB_HOME = ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION   = DOUBLE
DEBUG       = TRUE
PROFILE     = FALSE
DIM         = 2

COMP        = g++
FCOMP       = gfortran

USE_MPI     = TRUE
USE_MPI     = FALSE
USE_OMP     = TRUE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE

EBASE = amr

include $(TOP)/mk/Make.defs 
include ./Make.package

fincludes=$(includes)

Bdirs := BoxLib iamrlib amrlib mglib bndrylib tensorMG hgproj

Bpack	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	:= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)

ifdef FBOXLIB_HOME
  DEFINES += -DMG_USE_FBOXLIB

  include $(PBOXLIB_HOME)/MGT_Solver/Make.package
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)/MGT_Solver
  VPATH_LOCATIONS   += $(PBOXLIB_HOME)/MGT_Solver
 
  include $(FBOXLIB_HOME)/mg/FParallelMG.mak
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/mg/
  VPATH_LOCATIONS   += $(FBOXLIB_HOME)/mg
 
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/SPARSKIT
  VPATH_LOCATIONS   += $(FBOXLIB_HOME)/extern/SPARSKIT
 
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/LAPACK
  VPATH_LOCATIONS   += $(FBOXLIB_HOME)/extern/LAPACK
 
  include $(FBOXLIB_HOME)/extern/BLAS/Make.package
  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/extern/BLAS
  VPATH_LOCATIONS   += $(FBOXLIB_HOME)/extern/BLAS

  INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/boxlib
  VPATH_LOCATIONS   += $(FBOXLIB_HOME)/boxlib
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

$(executable):

include $(TOP)/mk/Make.rules
