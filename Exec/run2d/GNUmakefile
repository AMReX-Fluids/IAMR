#BOXLIB_HOME defines the directory in which we will find the BoxLib directory
BOXLIB_HOME = $(HOME)/GIT/BoxLib

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DEBUG       = TRUE

DIM         = 2

COMP        = g++
FCOMP       = gfortran

USE_MPI     = TRUE
USE_MPI     = FALSE

USE_F90_SOLVERS = TRUE
USE_F90_SOLVERS = FALSE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE

EBASE = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
