#AMREX_HOME defines the directory in which we will find the AMReX directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

DIM        = 3
COMP	   = gnu
DEBUG	   = FALSE
USE_MPI    = TRUE
USE_OPENMP = FALSE

PRECISION = DOUBLE

TINY_PROFILE = TRUE

USE_FLCTS = TRUE
USE_FLCTS = FALSE

USE_HYPRE = FALSE
USE_METIS = FALSE

USE_VELOCITY = TRUE
USE_VELOCITY = FALSE

USE_F90_SOLVERS = FALSE
USE_F90_SOLVERS = TRUE

USE_XBLAS = TRUE
USE_XBLAS = FALSE

EBASE     = amr

ifeq (${USE_FLCTS}, TRUE)
  DEFINES += -DBL_DO_FLCT
endif

ifeq (${USE_XBLAS}, TRUE)
   XTRADEFS   += -DXBLAS
   XTRAINCLOC += $(HOME)/tmp/xblas-1.0.248/src
   XTRALIBLOC += $(HOME)/tmp/xblas-1.0.248/src
   XTRALIBS   += -lxblas
endif

ifeq (${USE_VELOCITY}, TRUE)
  #AMRVIS_DIR defines the directory in which we will find pAmrvis (e.g. DataServices, AmrData and FABUTIL)
  AMRVIS_DIR = $(AMREX_HOME)/Src/Extern/amrdata
endif

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
DEFINES += -DMOREGENGETFORCE

# This is for the time-dependent turbulence forcing term
DEFINES += -DDO_IAMR_FORCE

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
