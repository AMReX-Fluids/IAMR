#
# $Id: GNUmakefile,v 1.65 1998-12-01 23:52:34 lijewski Exp $
#
PBOXLIB_HOME = ..

TOP = $(PBOXLIB_HOME)
#
# Where libraries and include files will be installed.
#
INSTALL_ROOT = $(TOP)
#
# Variables for the user to set ...
#
PRECISION     = DOUBLE
DEBUG	      = TRUE
PROFILE       = FALSE
DIM	      = 3
COMP          = KCC
USE_WINDOWS   = FALSE
USE_MPI       = TRUE
USE_ARRAYVIEW = FALSE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE
#
# The only stencils in which we're interested when using USE_HGPROJ_SERIAL.
#
ifeq ($(USE_HGPROJ_SERIAL),TRUE)
DEFINES += -DBL_USE_HGPROJ_SERIAL
ifeq ($(DIM),2)
PRVERSION = v9
else
PRVERSION = v7
endif
endif
#
# The base name of the library we're building.
#
LBASE = iamr

include $(TOP)/mk/Make.defs ./Make.package

ifeq ($(USE_HGPROJ_SERIAL),FALSE)
INCLUDE_LOCATIONS += . $(TOP)/include
else
INCLUDE_LOCATIONS += . $(TOP)/hgproj-serial
INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
INCLUDE_LOCATIONS += $(TOP)/include
endif

ifeq ($(USE_MPI),TRUE)
ifeq ($(MACHINE),OSF1)
INCLUDE_LOCATIONS += /usr/local/mpi/include
endif
ifeq ($(MACHINE),AIX)
INCLUDE_LOCATIONS += /usr/lpp/ppe.poe/include
endif
endif

ifeq ($(KCC_VERSION),3.3)
CXXFLAGS += --one_instantiation_per_object
endif
#
# Libraries to close against.
#
ifeq ($(COMP),KCC)
LibsToCloseAgainst := $(TOP)/lib/$(machineSuffix)/libmg$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libamr$(DIM)d.a
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbndry$(DIM)d.a
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libproj$(DIM)d.a
else
LibsToCloseAgainst += $(TOP)/hgproj-serial/lib/$(machineSuffix)/libproj$(DIM)d.$(PRVERSION).a
endif
LibsToCloseAgainst += $(TOP)/lib/$(machineSuffix)/libbox$(DIM)d.a
endif

all: $(optionsLib)

include $(TOP)/mk/Make.rules
