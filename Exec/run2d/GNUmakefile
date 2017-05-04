#AMREX_HOME defines the directory in which we will find the BoxLib directory
AMREX_HOME ?= ../../../amrex

#TOP defines the directory in which we will find Source, Exec, etc.
TOP = ../..

#
# Variables for the user to set ...
#

PRECISION   = DOUBLE

DEBUG       = TRUE

DIM         = 2

COMP        = g++

USE_MPI     = TRUE
USE_MPI     = FALSE

# This allows GetForce to call a Fortran routine
#DEFINES += -DGENGETFORCE

EBASE = amr

Bpack   := ./Make.package
Blocs   := .

include $(TOP)/Exec/Make.IAMR
