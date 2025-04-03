
#ifndef NS_LS_H
#define NS_LS_H

#include <AMReX_MultiFab.H>

using namespace amrex;

void phi_to_heavi(const Geometry& geom, int epsilon, MultiFab& phi, MultiFab& heaviside);
void heavi_to_rhoormu(MultiFab& heaviside, Real var1, Real var2, MultiFab& outmf);

Real calculate_eps (const Geometry& geom, int epsilon);
Real calculate_eps_one (const Geometry& geom, int reinit_levelset);
Real calculate_eps_two (const Geometry& geom, int reinit_levelset);
void levelset_diffcomp (Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& phi_cc_grad,
                        MultiFab& phi_ctime,
                        MultiFab& phi1,
                        MultiFab& phi2,
                        int epsG1,
                        int epsG2);

#endif
