#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

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

USE_PARTICLES   = TRUE

USE_MPI   = TRUE

# This allows GetForce to call a Fortran routine
DEFINES += -DBOUSSINESQ

EBASE     = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
