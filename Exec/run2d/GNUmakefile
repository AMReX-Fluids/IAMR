#AMREX_HOME defines the directory in which we will find the AMReX directory
AMREX_HOME       ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../../../AMReX-Hydro

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 2
COMP        = gnu

USE_MPI     = FALSE
USE_OMP     = FALSE
USE_CUDA    = FALSE
PROFILE     = FALSE

DEBUG         = FALSE
TEST          = TRUE
USE_ASSERTION = TRUE

USE_SENSEI_INSITU = FALSE

EBASE = amr

Blocs   := .

include $(TOP)/Exec/Make.IAMR
