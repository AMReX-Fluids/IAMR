#
# $Id: Make.package.mg,v 1.1 1997-07-08 23:08:04 vince Exp $
#

MGLIB_BASE=EXE

ifeq ($(LBASE),mg)
    MGLIB_BASE=LIB
endif

#
# Sources ...
#
C$(MGLIB_BASE)_headers += ABecLaplacian.H \
                          BndryData.H \
                          BoundCond.H \
                          CGSolver.H \
                          InterpBndryData.H \
                          LO_BCTYPES.H \
                          Laplacian.H \
                          LinOp.H \
                          Mask.H \
                          MultiGrid.H \
                          WriteMultiFab.H

C$(MGLIB_BASE)_sources += ABecLaplacian.C \
                          BndryData.C \
                          CGSolver.C \
                          InterpBndryData.C \
                          Laplacian.C \
                          LinOp.C \
                          Mask.C \
                          MultiGrid.C \
                          WriteMultiFab.C

F$(MGLIB_BASE)_headers += ABec_F.H \
                          CG_F.H \
                          INTERPBNDRYDATA_F.H \
                          LO_F.H \
                          LP_F.H \
                          MG_F.H

F$(MGLIB_BASE)_sources += ABec_$(DIM)D.F \
                          CG_$(DIM)D.F \
                          INTERPBNDRYDATA_$(DIM)D.F \
                          LO_$(DIM)D.F \
                          LO_UTIL.F \
                          LP_$(DIM)D.F \
                          MG_$(DIM)D.F
