#BOXLIB_HOME defines the directory in which we will find the BoxLib directory
BOXLIB_HOME ?= ../../../BoxLib

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION = DOUBLE

DEBUG	  = TRUE
DEBUG	  = FALSE

DIM       = 3

COMP	  = g++
FCOMP     = gfortran

USE_MPI   = TRUE

USE_F90_SOLVERS = FALSE
USE_F90_SOLVERS = TRUE

# This allows GetForce to call a Fortran routine
DEFINES += -DBOUSSINESQ

EBASE     = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
