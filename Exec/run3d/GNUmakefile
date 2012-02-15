#BOXLIB_HOME defines the directory in which we will find the BoxLib directory
BOXLIB_HOME ?= ../../../BoxLib

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION = DOUBLE

DEBUG	  = FALSE
DEBUG	  = TRUE

DIM       = 3

COMP	  = g++
FCOMP     = gfortran

USE_MPI   = TRUE
USE_MPI   = FALSE

USE_FLCTS = TRUE
USE_FLCTS = FALSE

USE_HYPRE = FALSE
USE_METIS = FALSE

USE_VELOCITY = FALSE
USE_VELOCITY = TRUE

USE_F90_SOLVERS = TRUE
USE_F90_SOLVERS = FALSE

EBASE     = amr

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

ifeq (${USE_VELOCITY}, TRUE)
  #AMRVIS_DIR defines the directory in which we will find pAmrvis (e.g. DataServices, AmrData and FABUTIL)
  #AMRVIS_DIR = $(BOXLIB_HOME)/Src/Extern/amrdata
  AMRVIS_DIR ?= ../../../Amrvis
endif

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE

# This is for the time-dependent turbulence forcing term
#DEFINES += -DDO_IAMR_FORCE

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
