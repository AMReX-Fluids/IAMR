#
# $Id: GNUmakefile,v 1.19 2000-04-05 18:29:29 sstanley Exp $
#
PBOXLIB_HOME = ../..

TOP = $(PBOXLIB_HOME)
#
# Variables for the user to set ...
#
PRECISION = DOUBLE
DEBUG	  = TRUE
DEBUG	  = FALSE
PROFILE   = FALSE
DIM       = 3
COMP      = KCC
USE_MPI   = TRUE
USE_MPI   = FALSE
#
# Use hgproj-serial -- only for testing.
#
# Touch Projection.cpp if you want to change the state of USE_HGPROJ_SERIAL.
#
USE_HGPROJ_SERIAL = TRUE
USE_HGPROJ_SERIAL = FALSE


USE_VAMPIRE = TRUE
USE_VAMPIRE = FALSE


ifeq ($(USE_HGPROJ_SERIAL),TRUE)
  #
  # What stencil do you want to use?
  #
  PRVERSION = v7
  EBASE = amr-serial
  DEFINES += -DBL_USE_HGPROJ_SERIAL
else
  EBASE = amr
endif

include $(TOP)/mk/Make.defs ./Make.package

BUILD_INPLACE:=FALSE
BUILD_INPLACE:=TRUE

#CXXFLAGS += --strict --display_error_number
#CXXFLAGS +=  --diag_suppress 450

ifeq ($(BUILD_INPLACE), TRUE)
  include $(TOP)/iamrlib/Make.package
  include $(TOP)/mglib/Make.package
  include $(TOP)/tensorMG/Make.package
  include $(TOP)/amrlib/Make.package
  include $(TOP)/bndrylib/Make.package
  include $(TOP)/hgproj/Make.package
  include $(TOP)/pBoxLib_2/Make.package
  INCLUDE_LOCATIONS += .
  INCLUDE_LOCATIONS += $(TOP)/iamrlib
  INCLUDE_LOCATIONS += $(TOP)/mglib
  INCLUDE_LOCATIONS += $(TOP)/tensorMG
  INCLUDE_LOCATIONS += $(TOP)/amrlib
  INCLUDE_LOCATIONS += $(TOP)/bndrylib
  INCLUDE_LOCATIONS += $(TOP)/hgproj
  INCLUDE_LOCATIONS += $(TOP)/pBoxLib_2
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
endif

ifeq ($(USE_MPI),TRUE)
  ifeq ($(MACHINE),OSF1)
    LIBRARIES += -lmpi
    INCLUDE_LOCATIONS += /usr/local/mpi/include
    LIBRARY_LOCATIONS += /usr/local/mpi/lib/alpha/ch_p4
  endif
  ifeq ($(MACHINE),AIX)
    LIBRARIES += -lmpi
    INCLUDE_LOCATIONS += /usr/lpp/ppe.poe/include
    LIBRARY_LOCATIONS += /usr/lpp/ppe.poe/lib
  endif
  ifeq ($(MACHINE), Linux)
    LIBRARIES += -lmpich
  endif
endif

ifeq ($(MACHINE),T3E)
  ifeq ($(USE_VAMPIRE),TRUE)
    #  use vampir trace
    LIBRARY_LOCATIONS += /usr/local/pkg/acts/VT151/lib
    LIBRARY_LOCATIONS += /opt/ctl/mpt/mpt/lib
    XTRALIBS += -lVT -lpmpi
  endif
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
  DEFINES += -DBL_USE_ARRAYVIEW
  DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

ifeq ($(MACHINE),IRIX64)
  ifeq ($(WHICHIRIX64),ESCHER)
    XTRALIBS += -IPA
  endif
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


XTRALIBS += --backend -lpat --backend pat.cld


ifeq ($(BUILD_INPLACE), TRUE)
  vpath %.cpp . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.H   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.h   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
  vpath %.F   . $(TOP)/iamrlib $(TOP)/mglib $(TOP)/tensorMG $(TOP)/amrlib $(TOP)/bndrylib $(TOP)/hgproj $(TOP)/pBoxLib_2
else
  vpath %.H   ..
  vpath %.cpp ..
  vpath %.a   $(LIBRARY_LOCATIONS)
endif

all: $(executable)

ifeq ($(MACHINE),T3E)
#
# We don't explicity set the path of the mpi directory on the T3E.
#
  $(executable): $(filter-out -lmpi, $(LIBRARIES))
else
  $(executable): $(LIBRARIES)
endif

ifneq ($(BUILD_INPLACE), TRUE)
#
# Build and install all libraries needed in an appropriate order.
#
libs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/amrlib;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
	cd $(TOP)/hgproj;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) LBASE=proj EBASE= install
else
	cd $(TOP)/hgproj-serial; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) PRVERSION=$(PRVERSION) USE_MPI=$(USE_MPI) LBASE=proj EBASE=
endif
	cd $(TOP)/mglib;     $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/tensorMG;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) install
	cd $(TOP)/iamrlib;   $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) PRVERSION=$(PRVERSION) USE_HGPROJ_SERIAL=$(USE_HGPROJ_SERIAL) install
#
# Cleanup libraries.
#
cleanlibs:
	cd $(TOP)/pBoxLib_2; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/bndrylib;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/amrlib;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
ifeq ($(USE_HGPROJ_SERIAL),FALSE)
	cd $(TOP)/hgproj;    $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) LBASE=proj EBASE= clean
else
	cd $(TOP)/hgproj-serial; $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) PRVERSION=$(PRVERSION) USE_MPI=$(USE_MPI) LBASE=proj EBASE= clean
endif
	cd $(TOP)/mglib;     $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/tensorMG;  $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) clean
	cd $(TOP)/iamrlib;   $(MAKE) PRECISION=$(PRECISION) PROFILE=$(PROFILE) COMP=$(COMP) DEBUG=$(DEBUG) DIM=$(DIM) USE_MPI=$(USE_MPI) USE_HGPROJ_SERIAL=$(USE_HGPROJ_SERIAL) clean
endif

include $(TOP)/mk/Make.rules
