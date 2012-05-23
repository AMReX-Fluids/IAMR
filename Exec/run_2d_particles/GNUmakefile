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

DIM       = 2

COMP	  = g++
FCOMP     = gfortran

USE_PARTICLES   = TRUE

USE_MPI   = TRUE

USE_F90_SOLVERS = TRUE

EBASE     = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
