#ifndef _BuildAlias_H_
#define _BuildAlias_H_

#include <BoxArray.H>
#include <MultiFab.H>

#include <hg_projector.H>

extern level_mesh make_level_mesh(const Box& domain, const BoxArray& bs);
extern amr_mesh make_amr_mesh(Box domain, int nlev, const int rat[],
		              const BoxArray *bs[]);
extern level_real make_level_real(const level_mesh& mesh,
                  MultiFab& mf, int icomp = -1);
extern amr_real make_amr_real(const amr_mesh& mesh, MultiFab mf[]);
extern amr_real make_amr_real(const amr_mesh& mesh, MultiFab& mf);

extern void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab& mf,
		       int level, int comp);

extern void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab* mf[],
		       int levc, int levf, int comp);

extern void bldAmrReal(amr_real& areal, amr_mesh mesh, MultiFab * mf_crse,
		       MultiFab * mf_fine, int levc, int levf, int comp);

#endif

