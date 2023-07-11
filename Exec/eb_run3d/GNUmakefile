#AMREX_HOME defines the directory in which we will find the AMReX directory
AMREX_HOME       ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../../../AMReX-Hydro

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 3
COMP        = gnu

USE_MPI     = TRUE
USE_OMP     = FALSE
PROFILE     = FALSE

USE_VELOCITY = FALSE

DEBUG       = FALSE

TEST          = TRUE
USE_ASSERTION = TRUE

USE_EB = TRUE

EBASE = amr

Blocs   := .

include $(TOP)/Exec/Make.IAMR
