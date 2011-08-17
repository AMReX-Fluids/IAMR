#BOXLIB_HOME defines the directory in which we will find the BoxLib directory
BOXLIB_HOME = $(HOME)/BoxLib

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION = DOUBLE

DEBUG	  = TRUE
DEBUG	  = FALSE

PROFILE   = TRUE
PROFILE   = FALSE

DIM       = 3

COMP	  = g++
FCOMP     = gfortran

USE_MPI   = TRUE
USE_MPI   = FALSE

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

ifeq (${USE_VELOCITY}, TRUE)
  #AMRVIS_DIR defines the directory in which we will find pAmrvis (e.g. DataServices, AmrData and FABUTIL)
  AMRVIS_DIR = $(HOME)/Parallel/pAmrvis
endif

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE

# This is for the time-dependent turbulence forcing term
#DEFINES += -DDO_IAMR_FORCE

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
