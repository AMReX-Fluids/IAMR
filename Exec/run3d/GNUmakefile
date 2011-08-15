#PBOXLIB_HOME defines the directory in which we will find the BoxLib directory
PBOXLIB_HOME = $(HOME)/GIT/BoxLib
 
#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION = DOUBLE

DEBUG	  = FALSE
DEBUG	  = TRUE

PROFILE   = TRUE
PROFILE   = FALSE

DIM       = 3

COMP	  = gcc
FCOMP     = gfortran

USE_MPI   = FALSE
USE_MPI   = TRUE

USE_FLCTS = TRUE
USE_FLCTS = FALSE

USE_HYPRE = FALSE
USE_METIS = FALSE

USE_BGL   = FALSE

USE_VELOCITY = TRUE
USE_VELOCITY = FALSE

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

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE

# This is for the time-dependent turbulence forcing term
#DEFINES += -DDO_IAMR_FORCE

Bpack   := ./Make.package
Blocs   := .

ifeq ($(WHICHLINUX), COLUMBIA)
o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.PROF.MPI.EXE -O3 -hlo0 -ip  -I. -I.. -I. -I../../LMC/ -I../../iamrlib/ -I../../tensorMG/ -I../../hgproj/ -I../../mglib/ -I../../amrlib/ -I../../bndrylib/ -I../../BoxLib/ -I../../iamrlib/InflowForce -I../../pAmrvis -I../../iamrlib/InflowForce -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.PROF.MPI.EXE/amr_real3d.2.o
o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.MPI.EXE -O3 -hlo0 -ip  -I. -I.. -I. -I../../LMC/ -I../../iamrlib/ -I../../tensorMG/ -I../../hgproj/ -I../../mglib/ -I../../amrlib/ -I../../bndrylib/ -I../../BoxLib/ -I../../iamrlib/InflowForce -I../../pAmrvis -I../../iamrlib/InflowForce -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.MPI.EXE/amr_real3d.2.o
o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o:
	ifort -module o/3d.Linux.Intel.Intel.EXE -O3 -hlo0 -ip -I. -I. -I../../BoxLib -I../../iamrlib -I../../amrlib -I../../tensorMG -I../../mglib -I../../bndrylib -I../../iamrlib/InflowForce -I../../hgproj -c  ../../hgproj/amr_real3d.2.f -o o/3d.Linux.Intel.Intel.EXE/amr_real3d.2.o
endif

include $(TOP)/Exec/Make.IAMR
