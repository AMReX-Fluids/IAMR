PRECISION = DOUBLE
DEBUG	= FALSE
DEBUG	= TRUE
PROFILE = FALSE

DIM	= 2
DIM	= 3

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
USE_BSP=TRUE
USE_BSP=FALSE
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

INCLUDE_LOCATIONS += $(HERE) ../pBoxLib_2

# bsp parallel locations
ifeq ($(USE_BSP),TRUE)
DEFINES += -DBL_USE_BSP
BSP_HOME = /usr/people/vince/Parallel/BSP/BSP
INCLUDE_LOCATIONS += $(BSP_HOME)/include
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/OSF1
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/OSF1/$(BSP_DEVICE)
#LIBRARIES += -lbspcore_O2 -lbsplevel1_O0
LIBRARIES += -lbspcore_O0 -lbsplevel1_O0
###### exception library (for newest bsplib)
# end bsp parallel locations
LIBRARY_LOCATIONS += /usr/ccs/lib/cmplrs/cc
LIBRARIES += -lexc
endif


# FillPatch switches
#DEFINES += -DUSEUNRAVELEDFILLPATCH=1
DEFINES += -DUSEUNRAVELEDFILLPATCH=0
#DEFINES += -DUSEOLDFILLPATCH=1
DEFINES += -DUSEOLDFILLPATCH=0
#DEFINES += -DNEWFPMINBOX=1
DEFINES += -DNEWFPMINBOX=0


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
endif

#--------------- netcdf library

ifeq ($(USE_NETCDF),TRUE)
LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin
endif

#DEFINES += -DNEWUTIL
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

#FC = f77  -warn declarations -extend_source -check_bounds
#FC = f90  -warn declarations -extend_source -check_bounds
FC = f77  -warn declarations -extend_source
FC = f90  -warn declarations -extend_source

CXXFLAGS +=
CXXOPTF +=
CXXDEBF +=

FFLAGS += 
#FOPTF += -fpe2
#FDEBF += -fpe2

CFLAGS +=
COPTF +=
CDEBF +=

include $(HERE)/Make.package 

vpath %.cpp :$(HERE) ../pBoxLib_2
vpath %.F :$(HERE)

all: $(executable)
#all: $(objForExecs) godzillaLink

#$(executable): $(optionsLib)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include ../mk/Make.rules
