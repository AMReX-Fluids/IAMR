#
# $Id: GNUmakefile,v 1.20 2001-03-28 18:51:09 car Exp $
#
PBOXLIB_HOME = ../..

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
PROFILE   = FALSE
DIM	  = 2
COMP	  = g++
COMP      = KCC
USE_MPI   = FALSE
NAMESPACE = TRUE
NAMESPACE = FALSE
USE_THREADS	  = TRUE
USE_THREADS	  = FALSE

#
# BoxLib3 Support
#
USE_BL3 = TRUE
USE_BL3 = FALSE
#
# Hypre support
#
USE_HYPRE=TRUE
USE_HYPRE=FALSE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE
#
# What stencil do you want to use?
#
PRVERSION = v5

ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  DEFINES += -DBL_USE_HGPROJ_SERIAL
  EBASE = amr-serial
else
  EBASE = amr
endif

ifeq ($(DIM),2)
#
# SYNCREG_2D.F builds RHS differently for 5 and 9 point stencils.
#
# Note that parallel hgproj currently only implements v5.
#
ifeq ($(findstring 9, $(PRVERSION)), 9)
DEFINES += -DBL_PRVERSION=9
else
DEFINES += -DBL_PRVERSION=5
endif
endif

include $(TOP)/mk/Make.defs 
include ./Make.package

BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE

ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/iamrlib/Make.package
  include $(TOP)/mglib/Make.package
  include $(TOP)/tensorMG/Make.package
  include $(TOP)/amrlib/Make.package
  include $(TOP)/bndrylib/Make.package
  ifeq ($(USE_HGPROJ_SERIAL),FALSE)
    include $(TOP)/hgproj/Make.package
  else
    include $(TOP)/hgproj-serial/Make.package
  endif
  include $(TOP)/pBoxLib_2/Make.package
  INCLUDE_LOCATIONS += .
  INCLUDE_LOCATIONS += $(TOP)/iamrlib
  INCLUDE_LOCATIONS += $(TOP)/mglib
  INCLUDE_LOCATIONS += $(TOP)/tensorMG
  INCLUDE_LOCATIONS += $(TOP)/amrlib
  INCLUDE_LOCATIONS += $(TOP)/bndrylib
  ifeq ($(USE_HGPROJ_SERIAL),FALSE)
    INCLUDE_LOCATIONS += $(TOP)/hgproj
  else
    INCLUDE_LOCATIONS += $(TOP)/hgproj-serial
    INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
  endif
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2
  ifeq ($(USE_HYPRE),TRUE)
    HYPRE_DIR := $(PBOXLIB_HOME)/hypre
    LIBRARY_LOCATIONS += $(HYPRE_DIR)/lib
    INCLUDE_LOCATIONS += $(HYPRE_DIR)/include
    INCLUDE_LOCATIONS += $(HYPRE_DIR)/HypreABec
    LIBRARIES +=  -lHYPRE_ls			\
	 -lHYPRE_mv				\
	 -lHYPRE_parcsr_ls			\
	 -lHYPRE_DistributedMatrixPilutSolver 	\
	 -lHYPRE_ParaSails 			\
	 -lHYPRE_MatrixMatrix 			\
	 -lHYPRE_DistributedMatrix 		\
	 -lHYPRE_IJ_mv				\
	 -lHYPRE_parcsr_mv			\
	 -lHYPRE_seq_mv				\
	 -lHYPRE_utilities
    CPPFLAGS += -DMG_USE_HYPRE
    include $(HYPRE_DIR)/HypreABec/Make.package
  endif
else
  ifeq ($(USE_HGPROJ_SERIAL),FALSE)
    INCLUDE_LOCATIONS += . .. $(TOP)/include
    LIBRARY_LOCATIONS += $(TOP)/lib/$(machineSuffix)
    LIBRARIES         += -liamr$(DIM)d -lmg$(DIM)d -lmcmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d -lbox$(DIM)d
  else
    INCLUDE_LOCATIONS += . .. $(TOP)/hgproj-serial
    INCLUDE_LOCATIONS += $(TOP)/hgproj-serial/include/$(DIM)d.$(PRVERSION)
    INCLUDE_LOCATIONS += $(TOP)/include
    LIBRARY_LOCATIONS += $(TOP)/hgproj-serial/lib/$(machineSuffix) $(TOP)/lib/$(machineSuffix)
    LIBRARIES         += -liamr$(DIM)d -lmg$(DIM)d -lmcmg$(DIM)d -lamr$(DIM)d -lbndry$(DIM)d -lproj$(DIM)d.$(PRVERSION) -lbox$(DIM)d
  endif
  ifeq ($(USE_HYPRE),TRUE)
    error hypre not here yet
  endif
endif

ifeq ($(USE_BL3),TRUE)
  DEFINES += -DBL3_PROFILING
  ifeq ($(USE_MPI),TRUE)
    DEFINES += -DBL3_USE_MPI
    CEXE_sources += BL3Mpi.cpp
  endif
  ifeq ($(USE_THREADS),TRUE)
    DEFINES += -DBL3_PTHREADS
  endif
  CEXE_headers += Thread.H Profiler.H BoxLib3.H Utility.H Assert.H Error.H
  CEXE_sources += BL3Thread.cpp BL3Profiler.cpp BL3Utility.cpp 
  CEXE_sources += BL3Error.cpp
  CEXE_sources += BL3Parallel.cpp  BL3BoxLib3.cpp
  INCLUDE_LOCATIONS += $(PBOXLIB_HOME)
  MySrcDirs += $(PBOXLIB_HOME)/BoxLib3
endif

ifeq ($(MACHINE),OSF1)
#
# Some additional stuff for our preferred development/debugging environment.
#
  ifeq ($(PRECISION),DOUBLE)
    FFLAGS += -real_size 64
  endif
  FDEBF += -C
  FDEBF += -warn argument_checking
  FDEBF += -warn declarations
  ifneq ($(FC), f90)
    FDEBF += -warn truncated_source
    FDEBF += -warn unused
  endif
endif
#
# For Running 3rd Only
#
3RD = 1
3RD =
ifdef 3RD
LDFLAGS += --link_command_prefix 3rd
LDFLAGS += -non_shared -v
endif

ifeq ($(COMP),KCC)
  CXXFLAGS += --strict
  CXXFLAGS += --display_error_number
  CXXFLAGS += --diag_suppress 450
endif

ifeq ($(BUILD_INPLACE), TRUE)
  vpath_cpp := . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib
  vpath_H   := . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib
  vpath_FH  := . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib
  vpath_h   := . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib
  vpath_F   := . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib
  ifeq ($(USE_HGPROJ_SERIAL),FALSE)
    vpath_cpp += $(TOP)/hgproj
    vpath_H   += $(TOP)/hgproj
    vpath_FH  += $(TOP)/hgproj
    vpath_h   += $(TOP)/hgproj
    vpath_F   += $(TOP)/hgproj
  else
    vpath_cpp += $(TOP)/hgproj-serial
    vpath_H   += $(TOP)/hgproj-serial
    vpath_FH  += $(TOP)/hgproj-serial
    vpath_h   += $(TOP)/hgproj-serial
    vpath_F   += $(TOP)/hgproj-serial
  endif
  ifeq ($(USE_HYPRE),TRUE)
    vpath_cpp += $(HYPRE_DIR)/HypreABec
    vpath_H   += $(HYPRE_DIR)/HypreABec
    vpath_F   += $(HYPRE_DIR)/HypreABec
  endif
  ifeq ($(USE_BL3),TRUE)
    vpath_cpp += $(MySrcDirs)
    vpath_H   += $(MySrcDirs)
  endif
  vpath_cpp += $(TOP)/pBoxLib_2
  vpath_H   += $(TOP)/pBoxLib_2
  vpath_FH  += $(TOP)/pBoxLib_2
  vpath_h   += $(TOP)/pBoxLib_2
  vpath_F   += $(TOP)/pBoxLib_2
  vpath %.cpp $(vpath_cpp)
  vpath %.H   $(vpath_H)
  vpath %.FH  $(vpath_FH)
  vpath %.h   $(vpath_h)
  vpath %.F   $(vpath_F)
else
  vpath %.H   ..
  vpath %.cpp ..
  vpath %.a   $(LIBRARY_LOCATIONS)
endif

all: $(executable)

$(executable): $(LIBRARIES)

ifneq ($(BUILD_INPLACE), TRUE)

#
# Build and install all libraries needed in an appropriate order.
#

MyMakeLine = $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) \
		DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) \
		USE_HGPROJ_SERIAL=$(USE_HGPROJ_SERIAL) PRVERSION=$(PRVERSION) \
		DEFINES="$(DEFINES)"

libs:
	cd $(TOP)/pBoxLib_2; $(MyMakeLine) install
	cd $(TOP)/bndrylib;  $(MyMakeLine) install
	cd $(TOP)/amrlib;    $(MyMakeLine) install
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
	cd $(TOP)/hgproj;    $(MyMakeLine) LBASE=proj EBASE= install
else
	cd $(TOP)/hgproj-serial; $(MyMakeLine) LBASE=proj EBASE=
endif
ifeq ($(USE_HYPRE),TRUE)
	cd $(HYPRE_DIR)/HypreABec; $(MyMakeLine) install
endif
	cd $(TOP)/mglib;     $(MyMakeLine) install
	cd $(TOP)/tensorMG;  $(MyMakeLine) install
	cd $(TOP)/iamrlib;   $(MyMakeLine) install
#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MyMakeLine) clean
	cd $(TOP)/bndrylib;  $(MyMakeLine) clean
	cd $(TOP)/amrlib;    $(MyMakeLine) clean
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
	cd $(TOP)/hgproj;    $(MyMakeLine) LBASE=proj EBASE= clean
else
	cd $(TOP)/hgproj-serial; $(MyMakeLine) LBASE=proj EBASE= clean
endif
ifeq ($(USE_HYPRE),TRUE)
	cd $(HYPRE_DIR)/HypreABec; $(MyMakeLine) clean
endif
	cd $(TOP)/mglib;     $(MyMakeLine) clean
	cd $(TOP)/tensorMG;  $(MyMakeLine) clean
	cd $(TOP)/iamrlib;   $(MyMakeLine) clean
endif

include $(TOP)/mk/Make.rules
