#AMREX_HOME defines the directory in which we will find the AMReX directory
AMREX_HOME       ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../../../AMReX-Hydro

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

DIM        = 3
COMP	   = gcc
DEBUG	   = FALSE
USE_MPI    = TRUE
USE_OMP    = FALSE

USE_CUDA = FALSE

USE_TURBULENT_FORCING = TRUE
USE_FAST_FORCE        = FALSE

#TEST=TRUE
#USE_ASSERTION=TRUE

PRECISION = DOUBLE

USE_HYPRE = FALSE
USE_METIS = FALSE

USE_VELOCITY = TRUE
USE_VELOCITY = FALSE

USE_XBLAS = TRUE
USE_XBLAS = FALSE

USE_SENSEI_INSITU = FALSE

EBASE     = amr

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

Blocs   := .

include ./Make.package
include $(TOP)/Exec/Make.IAMR
