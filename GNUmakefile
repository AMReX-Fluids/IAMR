PRECISION = DOUBLE
DEBUG	= TRUE
PROFILE = FALSE
DIM	= 2
PRVERSION = v9
COMP = KCC

EBASE = amr
LBASE = 

GUTHAMR_HOME=/usr/people/vince/Parallel/guthamr
ARRAYVIEWDIR  = /usr/local/ccse/ArrayView

include $(GUTHAMR_HOME)/mk/Make.defs

HERE = .

INCLUDE_LOCATIONS += $(HERE)
INCLUDE_LOCATIONS += $(GUTHAMR_HOME)/include
INCLUDE_LOCATIONS += ./include/$(DIM)d.$(PRVERSION)
LIBRARY_LOCATIONS += ./graphtools

# bsp parallel locations
BSP_HOME = /usr/people/vince/Parallel/BSP/BSP
INCLUDE_LOCATIONS += $(BSP_HOME)/include
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/OSF1
LIBRARY_LOCATIONS += $(BSP_HOME)/lib/OSF1/$(BSP_DEVICE)
LIBRARIES += -lbspcore_O2 -lbsplevel1_O0
###### exception library (for newest bsplib)
LIBRARY_LOCATIONS += /usr/ccs/lib/cmplrs/cc
LIBRARIES += -lexc
# end bsp parallel locations


# FillPatch switches
#DEFINES += -DUSEUNRAVELEDFILLPATCH=1
DEFINES += -DUSEUNRAVELEDFILLPATCH=0
DEFINES += -DUSEOLDFILLPATCH=1
#DEFINES += -DUSEOLDFILLPATCH=0
#DEFINES += -DNEWFPMINBOX=1
DEFINES += -DNEWFPMINBOX=0


LIBRARIES += -lgraph
LIBRARIES += -lX11 

#ifeq ($(DEBUG),TRUE)
INCLUDE_LOCATIONS += $(ARRAYVIEWDIR)
LIBRARY_LOCATIONS += $(ARRAYVIEWDIR) 
LIBRARIES += -larrayview$(DIM)d.$(machineSuffix) 
#endif

#--------------- netcdf library

LIBRARIES += /usr/people/stevens/bin/libnetcdf.a
INCLUDE_LOCATIONS += /usr/people/stevens/bin

DEFINES += -DGUTHAMR -DNEWUTIL -DHAS_WINDOWS

#3RD = 1
ifdef 3RD
# FOR RUNNING 3RD ONLY
CXXFLAGS += --link_command_prefix 3rd
CXXDEBF = +K0 --link_command_prefix 3rd
LIBRARIES += -ldnet_stub
FDEBF += -automatic
# FOR RUNNING 3RD ONLY
endif

CXXFLAGS += --diag_suppress 177
CXXOPTF +=
CXXDEBF +=

FFLAGS += 
FOPTF += 
FDEBF += -C

CFLAGS +=
COPTF +=
CDEBF +=

include $(HERE)/Make.package 

vpath %.C :$(HERE)
vpath %.F :$(HERE)
vpath %.a :$(GUTHAMR_HOME)/lib/$(machineSuffix)

all: $(executable)
#all: $(objForExecs) godzillaLink

#$(executable): $(optionsLib)

godzillaLink:
	rsh godzilla 'cd $(PWD); make $(executable)

include $(GUTHAMR_HOME)/mk/Make.rules
