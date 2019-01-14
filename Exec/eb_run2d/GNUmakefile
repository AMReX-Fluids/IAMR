#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 2
COMP        = gnu

DEBUG       = TRUE
USE_MPI     = TRUE
USE_OMP     = TRUE
PROFILE     = FALSE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE
#DEFINES += -DMOREGENGETFORCE
#DEFINES += -DBOUSSINESQ

USE_VELOCITY = TRUE

EBASE = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
