.SUFFIXES: .f .for

DIM=3

!IF $(DIM) == 2
HGDEF=
!ELSE
HGDEF=/DHG_CROSS_STENCIL
!ENDIF

FPP_FLAGS=/ansi /nologo \
	/S. /S..\pboxlib_2 $(HGDEF) \
	/DBL_NEED_RAND \
	/DBL_NO_FORT_FLUSH \
	/DBL_SPACEDIM=$(DIM) \
	/DBL_USE_DOUBLE \
	/DBL_LANG_FORT

.f.for:
	fpp $(FPP_FLAGS) $*.f | perl ..\scripts\strip72 -c > $*.for

FORT_FILES= \
	ABec_$(DIM)D.For \
	amr_real$(DIM)d.For \
	CG_$(DIM)D.For \
	COORDSYS_$(DIM)D.For \
	DERIVE_$(DIM)D.For \
	DIFFUSION_$(DIM)D.For \
	FILCC_$(DIM)D.For \
	FLUXREG_$(DIM)D.For \
	GODUNOV_$(DIM)D.For \
	hg_avg$(DIM)d.For \
	hg_multi$(DIM)d.For \
	hg_proj$(DIM)d.For \
	INTERPBNDRYDATA_$(DIM)D.For \
	INTERP_$(DIM)D.For \
	LO_$(DIM)D.For \
	LP_$(DIM)D.For \
	MACOPERATOR_$(DIM)D.For \
	MACPROJ_$(DIM)D.For \
	MG_$(DIM)D.For \
	NAVIERSTOKES_$(DIM)D.For \
	PROB_$(DIM)D.For \
	PROJECTION_$(DIM)D.For \
	SYNCREG_$(DIM)D.For \
	VISCOPERATOR_$(DIM)D.For

FORT_FILES_ND= \
	GODUNOV_F.For \
	LO_UTIL.For

all:	$(FORT_FILES) $(FORT_FILES_ND)
