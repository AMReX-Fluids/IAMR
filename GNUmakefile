PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
PROFILE = FALSE

DIM	= 2
DIM	= 3

WHICH_HG=_new
WHICH_HG=

COMP = KCC

USE_WINDOWS=FALSE
USE_BSP=TRUE
USE_BSP=FALSE
USE_MPI=FALSE
USE_MPI=FALSE
USE_NETCDF=FALSE
USE_ARRAYVIEW = TRUE
USE_ARRAYVIEW = FALSE

#
# What is the type of BSP device?
#
ifeq ($(USE_BSP),TRUE)
BSP_DEVICE = SHMEM_SYSV
endif

EBASE = amr
LBASE = 

HERE = .

PBOXLIB_HOME = ..
include ../mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

INCLUDE_LOCATIONS += . ../pBoxLib_2 ../amrlib ../bndrylib ../mglib ../hgproj$(WHICH_HG)

ifeq ($(USE_BSP), TRUE)
DEFINES += -DBL_USE_BSP
ifeq ($(BSP_MACHINE), OSF1)
BSP_HOME = /usr/people/vince/Parallel/BSP/BSP
endif
ifeq ($(MACHINE), T3E)
ifeq ($(WHICHT3E), NERSC)
DEFINES += -DBL_T3E_NERSC
BSP_HOME = /u1/vince/BSP
endif
ifeq ($(WHICHT3E), NAVO)
DEFINES += -DBL_T3E_NAVO
BSP_HOME = /home/Cvince/BSP
endif
endif
endif

ifeq ($(USE_MPI), TRUE)
DEFINES += -DBL_USE_MPI
MPI_HOME = /usr/local/mpi
endif



# FillPatch switches
#DEFINES += -DUSEUNRAVELEDFILLPATCH=1
DEFINES += -DUSEUNRAVELEDFILLPATCH=0
#DEFINES += -DUSEOLDFILLPATCH=1
DEFINES += -DUSEOLDFILLPATCH=0
#DEFINES += -DNEWFPMINBOX=1
DEFINES += -DNEWFPMINBOX=0

#
# Uncomment this is you want to use Parallel I/O.
#
DEFINES += -DBL_PARALLEL_IO

ifeq ($(USE_WINDOWS),TRUE)
LIBRARIES += -lgraph
LIBRARIES += -lX11 
DEFINES += -DBL_USE_WINDOWS
endif

ifeq ($(USE_ARRAYVIEW),TRUE)
#ARRAYVIEWDIR  = /usr/local/ccse/ArrayView
#ARRAYVIEWDIR  = /usr/people/vince/Visualization/ArrayView
ARRAYVIEWDIR  = .
INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
#LIBRARY_LOCATIONS += $(ARRAYVIEWDIR) 
#LIBRARIES += -larrayview$(DIM)d.$(machineSuffix) 
DEFINES += -DBL_USE_ARRAYVIEW
DEFINES += -DBL_ARRAYVIEW_TAGBOX
endif

ifeq ($(USE_BSP), TRUE)
INCLUDE_LOCATIONS += $(BSP_HOME)/include
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/$(BSP_MACHINE)
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/$(BSP_MACHINE)/$(BSP_DEVICE)
endif

ifeq ($(USE_MPI), TRUE)
INCLUDE_LOCATIONS += $(MPI_HOME)/include
LIBRARY_LOCATIONS += $(MPI_HOME)/lib/alpha/ch_p4
LIBRARIES += -lmpi
endif

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif

ifeq ($(USE_BSP), TRUE)
ifeq ($(DEBUG), TRUE)
LIBRARIES += -lbspcore_O0 -lbsplevel1_O0
else
ifeq ($(BSP_MACHINE), OSF1)
LIBRARIES += -lbspcore_O2 -lbsplevel1_O0
else
LIBRARIES += -lbspcore_O2 -lbsplevel1_O2
endif
endif
#
# exception library (for newest bsplib)
#
ifeq ($(BSP_MACHINE), OSF1)
LIBRARY_LOCATIONS += /usr/ccs/lib/cmplrs/cc
LIBRARIES += -lexc 
endif
endif

ifeq ($(BSP_MACHINE), OSF1)
#
# Some additional stuff for our preferred development/debugging environment.
#
ifeq ($(PRECISION), DOUBLE)
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

3RD = 1
3RD=
ifdef 3RD
# FOR RUNNING 3RD ONLY
LDFLAGS += --link_command_prefix 3rd
#CXXDEBF = +K0 --link_command_prefix 3rd -non_shared
LDFLAGS += -non_shared -v
#LIBRARIES += -ldnet_stub
#FDEBF += -automatic
# FOR RUNNING 3RD ONLY
endif

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

CFLAGS +=
COPTF +=
CDEBF +=

include $(HERE)/Make.package 

#INCLUDE_LOCATIONS += ../hgproj/include/3d.v7

FOPTF = -fast

vpath %.cpp : . ../pBoxLib_2 ../amrlib ../bndrylib ../mglib ../hgproj$(WHICH_HG)
vpath %.F   : . ../amrlib ../bndrylib ../mglib ../hgproj$(WHICH_HG)

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include ../mk/Make.rules
