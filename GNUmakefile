PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
PROFILE = FALSE

DIM	= 3
DIM	= 2

#
# Holy Grail stuff ...
#
ifeq ($(DIM),2)
#
# c5: -DHG_CONSTANT -DHG_CROSS_STENCIL
# v5: -DHG_CROSS_STENCIL
# c9: -DHG_CONSTANT
# v9:
#
CPPFLAGS +=
else
#
# c7: -DHG_CONSTANT -DHG_CROSS_STENCIL
# v7: -DHG_CROSS_STENCIL
#
CPPFLAGS += -DHG_CROSS_STENCIL
endif

COMP = KCC

USE_WINDOWS=FALSE
USE_BSP=FALSE
USE_BSP=TRUE
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

include ../mk/Make.defs

CPPFLAGS += -DBL_USE_NEW_HFILES

INCLUDE_LOCATIONS += . ../pBoxLib_2 ../amrlib

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
#DEFINES += -DBL_PARALLEL_IO

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
endif

###### exception library (for newest bsplib)
ifeq ($(BSP_MACHINE), OSF1)
LIBRARY_LOCATIONS += /usr/ccs/lib/cmplrs/cc
LIBRARIES += -lexc 
endif

#CXXFLAGS = -g --diag_suppress 177
#CXXFLAGS = --strict_warnings
ifeq ($(BSP_MACHINE), OSF1)
FFLAGS += -real_size 64
FDEBF += -C
FDEBF += - -warn argument_checking
FDEBF += - -warn declarations
FDEBF += - -warn truncated_source
FDEBF += - -warn unused
FOPTF  = -fast -O5 -tune ev5
endif

DEFINES += -DGUTHAMR 

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

FOPTF = -fast

vpath %.cpp : . ../pBoxLib_2 ../amrlib
vpath %.F   : . ../amrlib

all: $(executable)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include ../mk/Make.rules
