#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 3
COMP        = gnu

DEBUG       = TRUE
USE_MPI     = FALSE
USE_OMP     = FALSE
PROFILE     = FALSE


# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE
#DEFINES += -DBOUSSINESQ

USE_VELOCITY = FALSE

USE_EB = TRUE
ifeq ($(DIM),3)
USE_ALGOIM = TRUE
else
USE_ALGOIM = FALSE
endif

EBASE = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
