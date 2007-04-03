#
# $Id: GNUmakefile,v 1.44 2007-04-03 17:39:53 lijewski Exp $
#
PBOXLIB_HOME = ../..
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME = ../../../fParallel
FBOXLIB_HOME =

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
DEBUG	  = FALSE
PROFILE   = FALSE
PROFILE   = TRUE
DIM       = 3
COMP	  = Intel
FCOMP     = Intel
USE_MPI   = FALSE
USE_MPI   = TRUE
USE_FLCTS = FALSE
USE_FLCTS = TRUE
USE_HYPRE = FALSE
USE_METIS = FALSE
USE_XT3   = FALSE
USE_BGL   = FALSE
EBASE     = amr

ifeq ($(USE_BGL),TRUE)
  CXXOPTF += -O2 -qmaxmem=-1
  FOPTF   += -O2 -qmaxmem=-1
  fOPTF   += -O2 -qmaxmem=-1
endif

USE_VAMPIRE = TRUE
USE_VAMPIRE = FALSE

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

ifeq ($(USE_BGL),TRUE)
  NEEDS_FLUSH_F := TRUE
endif

include $(TOP)/mk/Make.defs

ifeq ($(USE_BGL),TRUE)
  MACHINE=BGL
  CXX=mpicxx.ibm
  FC=blrts_xlf
  fC=blrts_xlf
  CPPFLAGS+= -DMPICH_SKIP_MPICXX
  LDFLAGS+=-L/bgl/BlueLight/ppcfloor/bglsys/lib
  LDFLAGS+=-L/opt/ibmcmp/xlf/bg/10.1/blrts_lib
  BL_MPI_LIBS=-lxlf90 -lxlfmath
  FORTLINK:=LOWERCASE
endif

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
  CXXOPTF += -O -fastsse -Minline=levels:10 --no_exceptions
  FOPTF   += -O -fastsse -Minline=levels:10
  fOPTF   += -O -fastsse -Minline=levels:10
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

INCLUDE_LOCATIONS += . $(Blocs)
VPATH_LOCATIONS   += . $(Blocs)

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

ifeq ($(WHICHLINUX), COLUMBIA)
o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.PROF.MPI.EXE -O3 -hlo0 -ip  -I. -I.. -I. -I../../LMC/ -I../../iamrlib/ -I../../tensorMG/ -I../../hgproj/ -I../../mglib/ -I../../amrlib/ -I../../bndrylib/ -I../../BoxLib/ -I../../iamrlib/InflowForce -I../../pAmrvis -I../../iamrlib/InflowForce -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o
o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.MPI.EXE -O3 -hlo0 -ip  -I. -I.. -I. -I../../LMC/ -I../../iamrlib/ -I../../tensorMG/ -I../../hgproj/ -I../../mglib/ -I../../amrlib/ -I../../bndrylib/ -I../../BoxLib/ -I../../iamrlib/InflowForce -I../../pAmrvis -I../../iamrlib/InflowForce -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o
o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.EXE -O3 -hlo0 -ip -I. -I. -I../../BoxLib -I../../iamrlib -I../../amrlib -I../../tensorMG -I../../mglib -I../../bndrylib -I../../iamrlib/InflowForce -I../../hgproj -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o
endif

include $(TOP)/mk/Make.rules
