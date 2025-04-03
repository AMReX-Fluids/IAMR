#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex
AMREX_HYDRO_HOME ?= ../../../AMReX-Hydro

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DIM         = 2 # DIM
COMP        = gnu

DEBUG       = FALSE
USE_MPI     = TRUE # MPI
USE_OMP     = TRUE # OMP
PROFILE     = TRUE

USE_CUDA = FALSE # CUDA

USE_SENSEI_INSITU = FALSE

EBASE = amr

Blocs   := .

include $(TOP)/Exec/Make.IAMR
