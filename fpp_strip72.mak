.SUFFIXES: .f .for

FPP_FLAGS=/ansi /nologo /S. /S.\include\2d.v9 /DBL_NEED_RAND /DBL_NO_FORT_FLUSH /DBL_SPACEDIM=2 /DBL_USE_DOUBLE /DBL_LANG_FORT

.f.for:
	fpp $(FPP_FLAGS) $*.f | perl strip72 -c > $*.for

FORT_FILES= \
	ABec_2D.For \
	ABec_3D.For \
	amr_grav2d.For \
	amr_grav3d.For \
	amr_real2d.For \
	amr_real3d.For \
	CG_2D.For \
	CG_3D.For \
	cont2d.For \
	cont3d.For \
	COORDSYS_2D.For \
	COORDSYS_3D.For \
	DERIVE_2D.For \
	DERIVE_3D.For \
	DIFFUSION_2D.For \
	DIFFUSION_3D.For \
	FILCC_2D.For \
	FILCC_3D.For \
	FLUXREG_2D.For \
	FLUXREG_3D.For \
	GODUNOV_2D.For \
	GODUNOV_3D.For \
	GODUNOV_F.For \
	hg_avg2d.For \
	hg_avg3d.For \
	hg_multi2d.For \
	hg_multi3d.For \
	hg_proj2d.For \
	hg_proj3d.For \
	INTERPBNDRYDATA_2D.For \
	INTERPBNDRYDATA_3D.For \
	INTERP_2D.For \
	INTERP_3D.For \
	LO_2D.For \
	LO_3D.For \
	LO_UTIL.For \
	LP_2D.For \
	LP_3D.For \
	MACOPERATOR_2D.For \
	MACOPERATOR_3D.For \
	MACPROJ_2D.For \
	MACPROJ_3D.For \
	mainIBD_2D.For \
	mainIBD_3D.For \
	MAKESLICE_3D.For \
	MG_2D.For \
	MG_3D.For \
	NAVIERSTOKES_2D.For \
	NAVIERSTOKES_3D.For \
	PROB_2D.For \
	PROJECTION_2D.For \
	PROJECTION_3D.For \
	SYNCREG_2D.For \
	SYNCREG_3D.For \
	VISCOPERATOR_2D.For \
	VISCOPERATOR_3D.For

all:	$(FORT_FILES)
