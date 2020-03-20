
//
// Godunov is the object which calculates advective terms for iamr.
//

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <Godunov.H>
#include <GODUNOV_F.H>

#include <AMReX_BCRec.H>

#include <algorithm>

#ifdef AMREX_USE_EB
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB_utils.H>
#endif

#define GEOM_GROW 1
#define XVEL 0
#define YVEL 1
#define ZVEL 2

using namespace amrex;

namespace
{
    bool initialized = false;

    const int use_unlimited_slopes_DEF = 0;
}
//
// Set default values for these in Initialize()!!!
//
static int hyp_grow;

int  Godunov::verbose;
int  Godunov::ppm_type;
int  Godunov::slope_order;
int  Godunov::corner_couple;
int  Godunov::use_forces_in_trans;

void
Godunov::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    hyp_grow = 3;

    Godunov::verbose             = 0;
    Godunov::ppm_type            = 0;
    Godunov::slope_order         = 4;
    Godunov::corner_couple       = 1;
    Godunov::use_forces_in_trans = 0;

    int use_unlimited_slopes = use_unlimited_slopes_DEF;

    ParmParse pp("godunov");

    pp.query("v",                    verbose);
    pp.query("ppm_type",             ppm_type);
    pp.query("slope_order",          slope_order);
    pp.query("corner_couple",        corner_couple);
    pp.query("use_forces_in_trans",  use_forces_in_trans);
    pp.query("use_unlimited_slopes", use_unlimited_slopes);

    if (ppm_type == 2)
    {
	hyp_grow = 4;
    }

#ifdef AMREX_USE_EB

//     //nghost from incflo.

    hyp_grow = 4;
#endif

#if (BL_SPACEDIM==2)
    BL_ASSERT(slope_order==1 || slope_order==2 || slope_order==4);
#else
    BL_ASSERT(slope_order==1 || slope_order==4);
#endif

    set_params(slope_order, use_unlimited_slopes);

    amrex::ExecOnFinalize(Godunov::Finalize);

    initialized = true;
}

void
Godunov::Finalize ()
{
    initialized = false;
}

//
// Construct the Godunov Object.
//
Godunov::Godunov (int max_size)
{
  Initialize();
}

Godunov::~Godunov ()
{
    ;
}


//
// Advection functions follow.
//

//
// Compute upwinded FC velocities by extrapolating CC values in space and time
void
Godunov::ExtrapVelToFaces (const amrex::Box&  box,
                           const amrex_real*  dx,
                           amrex_real         dt,
                           D_DECL(      FArrayBox&  umac,       FArrayBox&  vmac,       FArrayBox&  wmac),
                           D_DECL(const Vector<int>& ubc, const Vector<int>& vbc, const Vector<int>& wbc),
                           const amrex::FArrayBox&  U,
                           amrex::FArrayBox&  tforces)
{
  extrap_vel_to_faces(box.loVect(),box.hiVect(),
                      BL_TO_FORTRAN_ANYD(U),
                      ubc.dataPtr(),BL_TO_FORTRAN_N_ANYD(tforces,0),BL_TO_FORTRAN_ANYD(umac),
                      vbc.dataPtr(),BL_TO_FORTRAN_N_ANYD(tforces,1),BL_TO_FORTRAN_ANYD(vmac),
#if (AMREX_SPACEDIM == 3)
                      wbc.dataPtr(),BL_TO_FORTRAN_N_ANYD(tforces,2),BL_TO_FORTRAN_ANYD(wmac),
                      &corner_couple,
#endif
                      &dt, dx, &use_forces_in_trans, &ppm_type);
}

void
Godunov::AdvectScalars(const Box&  box,
                       const Real* dx,
                       Real        dt,
                       D_DECL(const FArrayBox&   Ax, const FArrayBox&   Ay, const FArrayBox&   Az),
                       D_DECL(const FArrayBox& umac, const FArrayBox& vmac, const FArrayBox& wmac),   int vcomp,
                       D_DECL(      FArrayBox& xflx,       FArrayBox& yflx,       FArrayBox& zflx),   int flxcomp,
                       D_DECL(      FArrayBox& xstate,     FArrayBox& ystate,     FArrayBox& zstate), int ecomp,
                       const FArrayBox& Sfab,   int first_scalar, int num_scalars,
                       const FArrayBox& Forces, int fcomp,
                       const FArrayBox& Divu,   int ducomp,
                       FArrayBox& aofs,         int state_ind,
                       const amrex::Vector<AdvectionForm>& advectionType, const amrex::Vector<int>& state_bc,
                       AdvectionScheme adv_scheme, const amrex::FArrayBox& V)
{
    AMREX_ASSERT(Sfab.nComp() >= first_scalar + num_scalars);

    Vector<int> use_conserv_diff(num_scalars);
    for (int i=0; i<num_scalars; ++i) {
        use_conserv_diff[i] = (advectionType[state_ind+i] == Conservative) ? 1 : 0;
    }

    // Extrapolate to cell faces (store result in flux container)
    const int state_fidx = first_scalar + 1;
    extrap_state_to_faces(box.loVect(),box.hiVect(),
                          BL_TO_FORTRAN_N_ANYD(Sfab,first_scalar), &num_scalars,
                          BL_TO_FORTRAN_N_ANYD(Forces,fcomp),
                          BL_TO_FORTRAN_N_ANYD(Divu,ducomp),
                          BL_TO_FORTRAN_N_ANYD(umac,vcomp),     BL_TO_FORTRAN_N_ANYD(xstate,ecomp),
                          BL_TO_FORTRAN_N_ANYD(vmac,vcomp),     BL_TO_FORTRAN_N_ANYD(ystate,ecomp),
#if (AMREX_SPACEDIM == 3)
                          BL_TO_FORTRAN_N_ANYD(wmac,vcomp),     BL_TO_FORTRAN_N_ANYD(zstate,ecomp),
                          &corner_couple,
#endif
                          &dt, dx, &(state_bc[0]), &state_fidx,
                          &use_forces_in_trans, &ppm_type, &(use_conserv_diff[0]));

    //amrex::Print() << ystate;

    // ComputeAofs erase the edge state values to write the fluxes
    // So here we make a copy to keep separated fluxes and edge state
    xflx.copy<RunOn::Host>(xstate,ecomp,flxcomp,num_scalars);
    yflx.copy<RunOn::Host>(ystate,ecomp,flxcomp,num_scalars);
#if (AMREX_SPACEDIM == 3)
    zflx.copy<RunOn::Host>(zstate,ecomp,flxcomp,num_scalars);
#endif
    //
    // C component indices starts from 0, Fortran from 1
    //
    //int fort_ind = state_ind+1;

    // Convert face states to face fluxes (return in place) and compute flux divergence
    for (int i=0; i<num_scalars; ++i) { // FIXME: Loop required because conserv_diff flag only scalar
        ComputeAofs (box,
                     D_DECL(Ax,  Ay,  Az),  D_DECL(0,0,0),
                     D_DECL(umac,vmac,wmac),D_DECL(vcomp,vcomp,vcomp),
                     D_DECL(xflx,yflx,zflx),D_DECL(flxcomp+i,flxcomp+i,flxcomp+i),
                     V,0,aofs,state_ind+i,use_conserv_diff[i]);
    }
}


//
// Advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
// FArrayBox work sized as in edge_states.
//

void
Godunov::AdvectState (const Box&  box,
                      const Real* dx,
                      Real        dt,
                      FArrayBox&  areax,
                      FArrayBox&  uedge,
                      FArrayBox&  xflux,
                      FArrayBox&  areay,
                      FArrayBox&  vedge,
                      FArrayBox&  yflux,
#if (BL_SPACEDIM == 3)
                      FArrayBox&  areaz,
                      FArrayBox&  wedge,
                      FArrayBox&  zflux,
#endif
                      FArrayBox&  U,
                      FArrayBox&  S,
                      FArrayBox&  tforces,
                      FArrayBox&  divu,
                      int         fab_ind,
                      FArrayBox&  aofs,
                      int         aofs_ind,
                      int         iconserv,
                      int         state_ind,
                      const int*  bc,
                      AdvectionScheme scheme,
                      FArrayBox&  vol)
{
     //Compute edge states for an advected quantity.

    const int state_fidx = state_ind+1;
    const int nc = 1;

    extrap_state_to_faces(box.loVect(),box.hiVect(),
                          BL_TO_FORTRAN_N_ANYD(S,fab_ind), &nc,
                          BL_TO_FORTRAN_N_ANYD(tforces,fab_ind),
                          BL_TO_FORTRAN_ANYD(divu),
                          BL_TO_FORTRAN_ANYD(uedge),     BL_TO_FORTRAN_ANYD(xflux),
                          BL_TO_FORTRAN_ANYD(vedge),     BL_TO_FORTRAN_ANYD(yflux),
#if (AMREX_SPACEDIM == 3)
                          BL_TO_FORTRAN_ANYD(wedge),     BL_TO_FORTRAN_ANYD(zflux),
                          &corner_couple,
#endif
                          &dt, dx, bc, &state_fidx,
                          &use_forces_in_trans, &ppm_type, &iconserv);

    ComputeAofs (box,
                 D_DECL(areax,areay,areaz),D_DECL(0,0,0),
                 D_DECL(uedge,vedge,wedge),D_DECL(0,0,0),
                 D_DECL(xflux,yflux,zflux),D_DECL(0,0,0),
                 vol,0,aofs,aofs_ind,iconserv);

}

//
// Compute the advective derivative from fluxes.
//
void
Godunov::ComputeAofs (const Box& grd,
                      const FArrayBox& areax, const FArrayBox& uedge, const FArrayBox& xflux,
                      const FArrayBox& areay, const FArrayBox& vedge, const FArrayBox& yflux,
#if (BL_SPACEDIM == 3 )
                      const FArrayBox& areaz, const FArrayBox& wedge, const FArrayBox& zflux,
#endif
                      const FArrayBox& vol,
                      FArrayBox& aofs, int aofs_ind, int iconserv) const
{
    ComputeAofs(grd,D_DECL(areax,areay,areaz), D_DECL(0,0,0),
                D_DECL(uedge,vedge,wedge), D_DECL(0,0,0),
                D_DECL(xflux,yflux,zflux), D_DECL(0,0,0),
                vol, 0, aofs, aofs_ind, iconserv);
}

void
Godunov::ComputeAofs (const Box& grd,
                      D_DECL(const FArrayBox& areax,
                             const FArrayBox& areay,
                             const FArrayBox& areaz),
                      D_DECL(int axcomp,int aycomp,int azcomp),
                      D_DECL(const FArrayBox& uedge,
                             const FArrayBox& vedge,
                             const FArrayBox& wedge),
                      D_DECL(int ucomp, int vcomp, int wcomp),
                      D_DECL(const FArrayBox& xflux,
                             const FArrayBox& yflux,
                             const FArrayBox& zflux),
                      D_DECL(int fxcomp,int fycomp,int fzcomp),
                      const FArrayBox& vol, int volcomp,
                      FArrayBox& aofs,int acomp, int iconserv ) const
{

    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();

    adv_forcing( aofs.dataPtr(acomp),ARLIM(aofs.loVect()), ARLIM(aofs.hiVect()),

                      xflux.dataPtr(fxcomp), ARLIM(xflux.loVect()), ARLIM(xflux.hiVect()),
                      uedge.dataPtr(ucomp),  ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
                      areax.dataPtr(axcomp), ARLIM(areax.loVect()), ARLIM(areax.hiVect()),

                      yflux.dataPtr(fycomp), ARLIM(yflux.loVect()), ARLIM(yflux.hiVect()),
                      vedge.dataPtr(vcomp),  ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
                      areay.dataPtr(aycomp), ARLIM(areay.loVect()), ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)
                     zflux.dataPtr(fzcomp), ARLIM(zflux.loVect()), ARLIM(zflux.hiVect()),
                     wedge.dataPtr(wcomp),  ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
                     areaz.dataPtr(azcomp), ARLIM(areaz.loVect()), ARLIM(areaz.hiVect()),
#endif
                     vol.dataPtr(volcomp), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
                      lo, hi, &iconserv);
}

//
// Sync advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
//

void
Godunov::SyncAdvect (const Box&  box,
                     const Real* dx,
                     Real        dt,
                     int         level,
                     const FArrayBox& areax,
                     FArrayBox& uedge,
                     FArrayBox& ucorr,
                     FArrayBox& xflux,
                     const FArrayBox& areay,
                     FArrayBox& vedge,
                     FArrayBox& vcorr,
                     FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                     const FArrayBox& areaz,
                     FArrayBox& wedge,
                     FArrayBox& wcorr,
                     FArrayBox& zflux,
#endif
                     FArrayBox& U,
                     FArrayBox& S,
                     FArrayBox& tforces,
                     FArrayBox& divu,
                     int        fab_ind,
                     FArrayBox& sync,
                     int        sync_ind,
                     int        iconserv,
                     int        state_ind,
                     const int* bc,
                     AdvectionScheme scheme,
                     const FArrayBox& vol)
{
    //
    // Error block.
    //
    BL_ASSERT(S.nComp()       >= BL_SPACEDIM);
    BL_ASSERT(S.nComp()       >= fab_ind    );
    BL_ASSERT(tforces.nComp() >= fab_ind    );
    BL_ASSERT(sync.nComp()    >= sync_ind   );

    BL_ASSERT((ucorr.box()).contains(xflux.box()));
    // below Assert fails while above passes, so the two are not equivalent
    //BL_ASSERT(ucorr.box()     >= xflux.box());
    BL_ASSERT(ucorr.nComp()   >= 1          );

    BL_ASSERT((vcorr.box()).contains(yflux.box()));
    BL_ASSERT(vcorr.nComp()   >= 1          );
#if (BL_SPACEDIM == 3)
    BL_ASSERT((wcorr.box()).contains(zflux.box()));
    BL_ASSERT(wcorr.nComp()   >= 1          );
#endif


    //
    // Compute the edge states.
    //

    const int state_fidx = state_ind+1;
    const int nc = 1;

    extrap_state_to_faces(box.loVect(),box.hiVect(),
                          BL_TO_FORTRAN_N_ANYD(S,fab_ind), &nc,
                          BL_TO_FORTRAN_N_ANYD(tforces,fab_ind),
                          BL_TO_FORTRAN_N_ANYD(divu,0),
                          BL_TO_FORTRAN_ANYD(uedge),     BL_TO_FORTRAN_ANYD(xflux),
                          BL_TO_FORTRAN_ANYD(vedge),     BL_TO_FORTRAN_ANYD(yflux),
#if (AMREX_SPACEDIM == 3)
                          BL_TO_FORTRAN_ANYD(wedge),     BL_TO_FORTRAN_ANYD(zflux),
                          &corner_couple,
#endif
                          &dt, dx, bc, &state_fidx,
                          &use_forces_in_trans, &ppm_type, &iconserv);

    //
    // Compute the advective tendency for the mac sync.
    //
    ComputeSyncAofs(box,
                    areax, ucorr, xflux,
                    areay, vcorr, yflux,
#if (BL_SPACEDIM == 3)
                    areaz, wcorr, zflux,
#endif
                    vol, sync, sync_ind, iconserv);
}

//
// Compute the advective derivative of corrective fluxes for the mac sync.
//

void
Godunov::ComputeSyncAofs (const Box& grd,
                          const FArrayBox& areax,
                          FArrayBox& ucorr,
                          FArrayBox& xflux,
                          const FArrayBox& areay,
                          FArrayBox& vcorr,
                          FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                          const FArrayBox& areaz,
                          FArrayBox& wcorr,
                          FArrayBox& zflux,
#endif
                          const FArrayBox& vol,
                          FArrayBox& sync,
                          int        sync_ind,
                          int        iconserv)
{
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    sync_adv_forcing(sync.dataPtr(sync_ind), ARLIM(sync.loVect()), ARLIM(sync.hiVect()),

                          xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
                          ucorr.dataPtr(),ARLIM(ucorr.loVect()),ARLIM(ucorr.hiVect()),
                          areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),

                          yflux.dataPtr(),ARLIM(yflux.loVect()),ARLIM(yflux.hiVect()),
                          vcorr.dataPtr(),ARLIM(vcorr.loVect()),ARLIM(vcorr.hiVect()),
                          areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
#if (BL_SPACEDIM == 3)
                          zflux.dataPtr(),ARLIM(zflux.loVect()),ARLIM(zflux.hiVect()),
                          wcorr.dataPtr(),ARLIM(wcorr.loVect()),ARLIM(wcorr.hiVect()),
                          areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
#endif
                          vol.dataPtr(), ARLIM(vol.loVect()), ARLIM(vol.hiVect()),
                          lo, hi);
}

//
// Correct a conservatively-advected scalar for under-over shoots.
//
void
Godunov::ConservativeScalMinMax (FArrayBox& Sold,
                                 FArrayBox& Snew,
                                 int        ind_old_s,
                                 int        ind_old_rho,
                                 int        ind_new_s,
                                 int        ind_new_rho,
                                 const int* bc,
                                 const Box& grd)
{
    const int *solo       = Sold.loVect();
    const int *sohi       = Sold.hiVect();
    const int *snlo       = Snew.loVect();
    const int *snhi       = Snew.hiVect();
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const Real *Sold_dat  = Sold.dataPtr(ind_old_s);
    const Real *Rho_dat   = Sold.dataPtr(ind_old_rho);
    const Real *Snew_dat  = Snew.dataPtr(ind_new_s);
    const Real *Rhon_dat  = Snew.dataPtr(ind_new_rho);

#if (BL_SPACEDIM == 3)
    Box flatbox(grd);
    int zlen = flatbox.length(BL_SPACEDIM-1);
    flatbox.growHi(BL_SPACEDIM-1,3-zlen);
    FArrayBox smin(flatbox,1);
    FArrayBox smax(flatbox,1);
    const Real *smin_dat = smin.dataPtr();
    const Real *smax_dat = smax.dataPtr();
#endif

    consscalminmax (Sold_dat, Rho_dat, ARLIM(solo), ARLIM(sohi),
			 Snew_dat, Rhon_dat,ARLIM(snlo), ARLIM(snhi),
#if (BL_SPACEDIM == 3)
                         smin_dat, smax_dat,
                         ARLIM(lo), ARLIM(hi),
#endif
                         lo, hi, bc);
}

//
// Correct a convectively-advected scalar for under-over shoots.
//
void
Godunov::ConvectiveScalMinMax (FArrayBox& Sold,
                               FArrayBox& Snew,
                               int        ind_old,
                               int        ind_new,
                               const int* bc,
                               const Box& grd)
{
    const int *slo        = Sold.loVect();
    const int *shi        = Sold.hiVect();
    const int *snlo       = Snew.loVect();
    const int *snhi       = Snew.hiVect();
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const Real *Sold_dat  = Sold.dataPtr(ind_old);
    const Real *Snew_dat  = Snew.dataPtr(ind_new);

#if (BL_SPACEDIM == 3)
    Box flatbox(grd);
    int zlen = flatbox.length(BL_SPACEDIM-1);
    flatbox.growHi(BL_SPACEDIM-1,3-zlen);
    FArrayBox smin(flatbox,1);
    FArrayBox smax(flatbox,1);
    const Real *smin_dat = smin.dataPtr();
    const Real *smax_dat = smax.dataPtr();
#endif

    convscalminmax (Sold_dat,
                         ARLIM(slo), ARLIM(shi),
                         Snew_dat,
                         ARLIM(snlo), ARLIM(snhi),
#if (BL_SPACEDIM == 3)
                         smin_dat, smax_dat,
                         ARLIM(lo), ARLIM(hi),
#endif
                         lo, hi, bc);
}


//
// Estimate the maximum change in velocity magnitude since previous iteration.
//
Real
Godunov::maxchng_velmag (FArrayBox&  U_old,
			 FArrayBox&  U_new,
			 const Box&  grd)
{
    BL_ASSERT( U_old.nComp()   >= BL_SPACEDIM );
    BL_ASSERT( U_new.nComp()   >= BL_SPACEDIM );

    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const int *uo_lo  = U_old.loVect();
    const int *uo_hi  = U_old.hiVect();
    const int *un_lo  = U_new.loVect();
    const int *un_hi  = U_new.hiVect();
    const Real *Uodat = U_old.dataPtr();
    const Real *Undat = U_new.dataPtr();

	Real max_change = 0.0;
    fort_maxchng_velmag(Uodat, ARLIM(uo_lo), ARLIM(uo_hi),
               	   Undat, ARLIM(un_lo), ARLIM(un_hi),
               	   lo, hi, &max_change);
    return max_change;
}

//
// Estimate the extrema of umac edge velocities and rho.
//

Real
Godunov::test_umac_rho (FArrayBox&  umac,
                        FArrayBox&  vmac,
#if (BL_SPACEDIM == 3)
                        FArrayBox&  wmac,
#endif
                        FArrayBox&  rho,
                        const Box&  grd,
                        const Real* dx,
                        const Real  dt,
                        const Real* u_max)
{
    //
    // Test block.
    //
    D_TERM(BL_ASSERT(umac.nComp() == 1);,
           BL_ASSERT(vmac.nComp() == 1);,
           BL_ASSERT(wmac.nComp() == 1););

    BL_ASSERT(rho.nComp()  == 1);

    const int *lo  = grd.loVect();
    const int *hi  = grd.hiVect();
    const int *ulo = umac.loVect();
    const int *uhi = umac.hiVect();
    const int *vlo = vmac.loVect();
    const int *vhi = vmac.hiVect();
    const int *rlo = rho.loVect();
    const int *rhi = rho.hiVect();
    const Real *um = umac.dataPtr();
    const Real *vm = vmac.dataPtr();
    const Real *rh = rho.dataPtr();

#if (BL_SPACEDIM == 3)
    const int *wlo = wmac.loVect();
    const int *whi = wmac.hiVect();
    const Real *wm = wmac.dataPtr();
#endif

    Real cfl;
    fort_test_umac_rho(um, ARLIM(ulo), ARLIM(uhi),
                       vm, ARLIM(vlo), ARLIM(vhi),
#if (BL_SPACEDIM == 3)
                       wm, ARLIM(wlo), ARLIM(whi),
#endif
                       rh, ARLIM(rlo), ARLIM(rhi),
                       lo, hi, &dt, dx, &cfl, u_max);
    return cfl;
}

//
// Source term functions follow
//

//
// Compute the update rule, this is useful for 1st order RK.
//
// psi^n+1 = psi^n + dt*tf^n
//

void
Godunov::Add_tf (const FArrayBox& Sold,
                 FArrayBox& Snew,
                 int        start_ind,
                 int        num_comp,
                 const FArrayBox& tforces,
                 int        tf_ind,
                 const Box& grd,
                 Real       dt) const
{
    BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
    BL_ASSERT(tforces.nComp() >= tf_ind    + num_comp);

    const int *solo   = Sold.loVect();
    const int *sohi   = Sold.hiVect();
    const int *snlo   = Snew.loVect();
    const int *snhi   = Snew.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SOdat = Sold.dataPtr(start_ind);
    Real *SNdat = Snew.dataPtr(start_ind);
    const Real *TFdat = tforces.dataPtr(tf_ind);

    update_tf(SOdat, ARLIM(solo), ARLIM(sohi),
                   SNdat, ARLIM(snlo), ARLIM(snhi),
                   TFdat, ARLIM(tlo), ARLIM(thi),
                   lo, hi, &dt, &num_comp);
}

//
// Compute the update rule
//
// psi^n+1 = psi^n - dt*aofs + dt*tforces
//

void
Godunov::Add_aofs_tf (const FArrayBox& Sold,
                      FArrayBox& Snew,
                      int        start_ind,
                      int        num_comp,
                      const FArrayBox& Aofs,
                      int        aofs_ind,
                      const FArrayBox& tforces,
                      int        tf_ind,
                      const Box& grd,
                      Real       dt) const
{
    BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Aofs.nComp()    >= aofs_ind  + num_comp);
    BL_ASSERT(tforces.nComp() >= tf_ind    + num_comp);

    const int *solo   = Sold.loVect();
    const int *sohi   = Sold.hiVect();
    const int *snlo   = Snew.loVect();
    const int *snhi   = Snew.hiVect();
    const int *alo    = Aofs.loVect();
    const int *ahi    = Aofs.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SOdat = Sold.dataPtr(start_ind);
    Real *SNdat = Snew.dataPtr(start_ind);
    const Real *AOdat = Aofs.dataPtr(aofs_ind);
    const Real *TFdat = tforces.dataPtr(tf_ind);

    update_aofs_tf(SOdat, ARLIM(solo), ARLIM(sohi),
                        SNdat, ARLIM(snlo), ARLIM(snhi),
                        AOdat, ARLIM(alo), ARLIM(ahi),
                        TFdat, ARLIM(tlo), ARLIM(thi),
                        lo, hi, &dt, &num_comp);
}

//
// Compute the update rule for velocities
//
// psi^n+1 = psi^n - dt*aofs - dt*gp/rho + dt*tforces
//

void
Godunov::Add_aofs_tf_gp (const FArrayBox& Uold,
                         FArrayBox& Unew,
                         const FArrayBox& Aofs,
                         const FArrayBox& tforces,
                         const FArrayBox& gp,
                         const FArrayBox& rho,
                         const Box& grd,
                         Real       dt) const
{
    BL_ASSERT(Unew.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(Aofs.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(tforces.nComp() >= BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      == BL_SPACEDIM);
    BL_ASSERT(rho.nComp()     == 1          );

    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const int *uolo   = Uold.loVect();
    const int *uohi   = Uold.hiVect();
    const int *unlo   = Unew.loVect();
    const int *unhi   = Unew.hiVect();
    const int *alo    = Aofs.loVect();
    const int *ahi    = Aofs.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *UOdat = Uold.dataPtr();
    Real *UNdat = Unew.dataPtr();
    const Real *AOdat = Aofs.dataPtr();
    const Real *TFdat = tforces.dataPtr();
    const Real *GPdat = gp.dataPtr();
    const Real *RHdat = rho.dataPtr();

    update_aofs_tf_gp(UOdat, ARLIM(uolo), ARLIM(uohi),
                           UNdat, ARLIM(unlo), ARLIM(unhi),
                           AOdat, ARLIM(alo), ARLIM(ahi),
                           TFdat, ARLIM(tlo), ARLIM(thi),
                           GPdat, ARLIM(glo), ARLIM(ghi),
                           RHdat, ARLIM(rlo), ARLIM(rhi),
                           lo, hi, &dt);
}

//
// Compute total source term for velocities, weighted by rho.
//
// tforces = (tforces - gp)/rho
//

void
Godunov::Sum_tf_gp (FArrayBox& tforces, int Tcomp,
                    const FArrayBox& gp, int Gcomp,
                    const FArrayBox& rho, int Rcomp) const
{
    BL_ASSERT(rho.nComp()     > Rcomp);
    BL_ASSERT(tforces.nComp() > Tcomp + BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      > Gcomp + BL_SPACEDIM);

    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    Real *TFdat = tforces.dataPtr(Tcomp);
    const Real *GPdat = gp.dataPtr(Gcomp);
    const Real *RHdat = rho.dataPtr(Rcomp);

    fort_sum_tf_gp(TFdat, ARLIM(tlo), ARLIM(thi),
                   GPdat, ARLIM(glo), ARLIM(ghi),
                   RHdat, ARLIM(rlo), ARLIM(rhi),
                   tlo, thi);
}

//
// Compute total source term for velocities, weighted by rho.
//
// tforces = (tforces + visc - gp)/rho
//

void
Godunov::Sum_tf_gp_visc (FArrayBox&       tforces,
                         const FArrayBox& visc,
                         const FArrayBox& gp,
                         const FArrayBox& Rho) const
{
    Sum_tf_gp_visc (tforces, 0, visc, 0, gp, 0, Rho, 0);
}

void
Godunov::Sum_tf_gp_visc (FArrayBox&       tforces,
                         int              Tcomp,
                         const FArrayBox& visc,
                         int              Vcomp,
                         const FArrayBox& gp,
                         int              Gcomp,
                         const FArrayBox& rho,
                         int              Rcomp) const
{
    BL_ASSERT(rho.nComp()     > Rcomp);
    BL_ASSERT(tforces.nComp() >= Tcomp+BL_SPACEDIM);
    BL_ASSERT(visc.nComp()    >= Vcomp+BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      == Gcomp+BL_SPACEDIM);

    const int *vlo    = visc.loVect();
    const int *vhi    = visc.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    Real *TFdat = tforces.dataPtr(Tcomp);
    const Real *VIdat = visc.dataPtr(Vcomp);
    const Real *GPdat = gp.dataPtr(Gcomp);
    const Real *RHdat = rho.dataPtr(Rcomp);

    fort_sum_tf_gp_visc(TFdat, ARLIM(tlo), ARLIM(thi),
                        VIdat, ARLIM(vlo), ARLIM(vhi),
                        GPdat, ARLIM(glo), ARLIM(ghi),
                        RHdat, ARLIM(rlo), ARLIM(rhi),
                        tlo, thi);
}

//
// Compute total source term for scalars.  Note for compatibility
// The switch iconserv, determines the form of the total source term
//
// iconserv==1   => tforces = tforces - divU*S
//
// iconserv==0   => tforces = (tforces)/rho
//

void
Godunov::Sum_tf_divu (const FArrayBox& S,
                      int        s_ind,
                      FArrayBox& tforces,
                      int        t_ind,
                      int        num_comp,
                      const FArrayBox& divu,
                      int        d_ind,
                      const FArrayBox& rho,
                      int        r_ind,
                      int        iconserv) const
{
    BL_ASSERT(S.nComp()       >= s_ind+num_comp);
    BL_ASSERT(tforces.nComp() >= t_ind+num_comp);
    BL_ASSERT(divu.nComp()    > d_ind          );
    BL_ASSERT(rho.nComp()     > r_ind          );

    const int *slo    = S.loVect();
    const int *shi    = S.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *dlo    = divu.loVect();
    const int *dhi    = divu.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *Sdat  = S.dataPtr(s_ind);
    Real *TFdat = tforces.dataPtr(t_ind);
    const Real *DUdat = divu.dataPtr(d_ind);
    const Real *RHdat = rho.dataPtr(r_ind);

    fort_sum_tf_divu(Sdat,  ARLIM(slo), ARLIM(shi),
                     TFdat, ARLIM(tlo), ARLIM(thi),
                     DUdat, ARLIM(dlo), ARLIM(dhi),
                     RHdat, ARLIM(rlo), ARLIM(rhi),
                     tlo, thi, &num_comp, &iconserv);
}

//
// Compute total source term for scalars.  Note for compatibility
// The switch iconserv, determines the form of the total source term
//
// iconserv==1   => tforces = tforces + visc - divU*S
//
// iconserv==0   => tforces = (tforces+ visc)/rho
//
void
Godunov::Sum_tf_divu_visc (const FArrayBox& S,
                           FArrayBox& tforces,
                           int s_ind, int num_comp,
                           const FArrayBox& visc, int v_ind,
                           const FArrayBox& divu,
                           const FArrayBox& rho,
                           int iconserv) const
{
    Sum_tf_divu_visc(S, s_ind, tforces, s_ind, num_comp,
                     visc, v_ind, divu, 0, rho, 0, iconserv);
}


void
Godunov::Sum_tf_divu_visc (const FArrayBox& S,
                           int        s_ind,
                           FArrayBox& tforces,
                           int        t_ind,
                           int        num_comp,
                           const FArrayBox& visc,
                           int        v_ind,
                           const FArrayBox& divu,
                           int        d_ind,
                           const FArrayBox& rho,
                           int        r_ind,
                           int        iconserv) const
{
    BL_ASSERT(S.nComp()       >= s_ind+num_comp);
    BL_ASSERT(tforces.nComp() >= t_ind+num_comp);
    BL_ASSERT(divu.nComp()    >  d_ind);
    BL_ASSERT(visc.nComp()    >= v_ind+num_comp);
    BL_ASSERT(rho.nComp()     >  r_ind);

    const int *slo    = S.loVect();
    const int *shi    = S.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *dlo    = divu.loVect();
    const int *dhi    = divu.hiVect();
    const int *vlo    = visc.loVect();
    const int *vhi    = visc.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *Sdat  = S.dataPtr(s_ind);
    Real *TFdat = tforces.dataPtr(t_ind);
    const Real *DUdat = divu.dataPtr(d_ind);
    const Real *VIdat = visc.dataPtr(v_ind);
    const Real *RHdat = rho.dataPtr(r_ind);

    fort_sum_tf_divu_visc(Sdat,  ARLIM(slo), ARLIM(shi),
                          TFdat, ARLIM(tlo), ARLIM(thi),
                          DUdat, ARLIM(dlo), ARLIM(dhi),
                          VIdat, ARLIM(vlo), ARLIM(vhi),
                          RHdat, ARLIM(rlo), ARLIM(rhi),
                          tlo, thi, &num_comp, &iconserv);
}


bool
Godunov::are_any(const Vector<AdvectionForm>& advectionType,
                 const AdvectionForm         testForm,
                 const int                   sComp,
                 const int                   nComp)
{
    for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
        if (advectionType[comp] == testForm)
            return true;
    }

    return false;
}

int
Godunov::hypgrow ()
{
    return hyp_grow;
}

int
Godunov::how_many(const Vector<AdvectionForm>& advectionType,
                  const AdvectionForm         testForm,
                  const int                   sComp,
                  const int                   nComp)
{
    int counter = 0;

    for (int comp = sComp; comp < sComp + nComp; ++comp)
    {
        if (advectionType[comp] == testForm)
            ++counter;
    }

    return counter;
}

#ifdef AMREX_USE_EB


void Godunov::ComputeSlopes( MultiFab&  a_state,
                             int a_comp,
                             D_DECL(MultiFab& a_xsl,
                                    MultiFab& a_ysl,
                                    MultiFab& a_zsl),
                             int a_sl_comp,
                             int a_ncomp,
                             const Vector<BCRec>& a_bcs,
                             const Box& a_domain)
{
    BL_PROFILE("Godunov::ComputeSlopes");
    AMREX_ALWAYS_ASSERT(a_state.nComp() >= a_comp    + a_ncomp - 1);
    D_TERM(AMREX_ALWAYS_ASSERT(a_xsl.nComp()   >= a_sl_comp + a_ncomp - 1);,
	   AMREX_ALWAYS_ASSERT(a_ysl.nComp()   >= a_sl_comp + a_ncomp - 1);,
	   AMREX_ALWAYS_ASSERT(a_zsl.nComp()   >= a_sl_comp + a_ncomp - 1); )
    AMREX_ALWAYS_ASSERT(a_bcs.size()    == a_ncomp);

    EB_set_covered(a_state, a_comp, a_ncomp, 1, COVERED_VAL);

    // We initialize slopes to zero in the grown domain ... this is essential
    //    to handle the up-winding at outflow faces (see inclfo)
    D_TERM( a_xsl.setVal(0.0, a_sl_comp, a_ncomp, a_xsl.nGrow());,
            a_ysl.setVal(0.0, a_sl_comp, a_ncomp, a_ysl.nGrow());,
            a_zsl.setVal(0.0, a_sl_comp, a_ncomp, a_zsl.nGrow()););

    // ... then set them to this large number in the interior in order to be sure
    //     that no "bad" values go unnoticed
    D_TERM( a_xsl.setVal( COVERED_VAL, a_sl_comp, a_ncomp,0);,
            a_ysl.setVal( COVERED_VAL, a_sl_comp, a_ncomp,0);,
            a_zsl.setVal( COVERED_VAL, a_sl_comp, a_ncomp,0););

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box bx = mfi.tilebox ();

        // This is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  Sborder_fab = static_cast<EBFArrayBox const&>(a_state[mfi]);
        const EBCellFlagFab&      flags = Sborder_fab.getEBCellFlagFab();

        if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
        {
            const auto& state_fab = a_state.array(mfi);
            D_TERM( const auto& xs_fab = a_xsl.array(mfi);,
                    const auto& ys_fab = a_ysl.array(mfi);,
                    const auto& zs_fab = a_zsl.array(mfi););

            // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
            if ( flags.getType(amrex::grow(bx,1)) == FabType::regular )
            {
                AMREX_FOR_4D(bx, a_ncomp, i, j, k, n,
                {
                    // X direction
                    Real du_xl = 2.0*(state_fab(i  ,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));
                    Real du_xr = 2.0*(state_fab(i+1,j,k,a_comp+n) - state_fab(i  ,j,k,a_comp+n));
                    Real du_xc = 0.5*(state_fab(i+1,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,a_sl_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;

                    // Y direction
                    Real du_yl = 2.0*(state_fab(i,j  ,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));
                    Real du_yr = 2.0*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j  ,k,a_comp+n));
                    Real du_yc = 0.5*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,a_sl_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;

#if (AMREX_SPACEDIM == 3)
                    // Z direction
                    Real du_zl = 2.0*(state_fab(i,j,k  ,a_comp+n) - state_fab(i,j,k-1,a_comp+n));
                    Real du_zr = 2.0*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k  ,a_comp+n));
                    Real du_zc = 0.5*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k-1,a_comp+n));

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,a_sl_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
#endif
                });

                Gpu::synchronize();
            }
            else
            {
                const auto& flag_fab = flags.array();

                AMREX_FOR_4D(bx, a_ncomp, i, j, k, n,
                {
                    if (flag_fab(i,j,k).isCovered())
                    {
                        D_TERM( xs_fab(i,j,k,a_sl_comp+n) = 0.0;,
                                ys_fab(i,j,k,a_sl_comp+n) = 0.0;,
                                zs_fab(i,j,k,a_sl_comp+n) = 0.0;);
                    }
                    else
                    {
                        // X direction
                        Real du_xl = (!flag_fab(i,j,k).isConnected(-1,0,0)) ? 0.0 :
                            2.0*(state_fab(i  ,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));
                        Real du_xr = (!flag_fab(i,j,k).isConnected( 1,0,0)) ? 0.0 :
                            2.0*(state_fab(i+1,j,k,a_comp+n) - state_fab(i  ,j,k,a_comp+n));
                        Real du_xc = 0.5*(state_fab(i+1,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));

                        Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                        xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                        xs_fab(i,j,k,a_sl_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;

                        // Y direction
                        Real du_yl = (!flag_fab(i,j,k).isConnected(0,-1,0)) ? 0.0 :
                            2.0*(state_fab(i,j  ,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));
                        Real du_yr = (!flag_fab(i,j,k).isConnected(0, 1,0)) ? 0.0 :
                            2.0*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j  ,k,a_comp+n));
                        Real du_yc = 0.5*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));

                        Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                        yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                        ys_fab(i,j,k,a_sl_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
#if (AMREX_SPACEDIM == 3)
                        // Z direction
                        Real du_zl = (!flag_fab(i,j,k).isConnected(0,0,-1)) ? 0.0 :
                            2.0*(state_fab(i,j,k  ,a_comp+n) - state_fab(i,j,k-1,a_comp+n));
                        Real du_zr = (!flag_fab(i,j,k).isConnected(0,0, 1)) ? 0.0 :
                            2.0*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k  ,a_comp+n));
                        Real du_zc = 0.5*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k-1,a_comp+n));

                        Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                        zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                        zs_fab(i,j,k,a_sl_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
#endif
                    }
                });

                Gpu::synchronize();
            } // end of cut cell region

            // MR: CHECK
            // Deal with BCs -- this part must be double checked
            // For now, assume that BCs are in place, i.e. ghost nodes are filled with the correct
            // Dirichlet's value.
            // Also, for the time being we assume that ext_dir is the only BC type which require
            // a special treatment
            //
            const auto& flag_fab  = flags.array();
            const auto& bc = a_bcs.dataPtr();

            AMREX_FOR_4D(bx, a_ncomp, i, j, k, n,
            {
                if ( (i == a_domain.smallEnd(0)) && !flag_fab(i,j,k).isCovered() && (bc[n].lo(0) == BCType::ext_dir) )
                {
                    Real du_xl = 2.0*(state_fab(i  ,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));
                    Real du_xr = 2.0*(state_fab(i+1,j,k,a_comp+n) - state_fab(i  ,j,k,a_comp+n));
                    Real du_xc = (state_fab(i+1,j,k,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i-1,j,k,a_comp+n))/3.0;

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,a_sl_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
                }
                if ( (i == a_domain.bigEnd(0)) && !flag_fab(i,j,k).isCovered() && (bc[n].hi(0) == BCType::ext_dir) )
                {
                    Real du_xl = 2.0*(state_fab(i  ,j,k,a_comp+n) - state_fab(i-1,j,k,a_comp+n));
                    Real du_xr = 2.0*(state_fab(i+1,j,k,a_comp+n) - state_fab(i  ,j,k,a_comp+n));
                    Real du_xc = -(state_fab(i-1,j,k,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i+1,j,k,a_comp+n))/3.0;

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,a_sl_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
                }

                if ( (j == a_domain.smallEnd(1)) && !flag_fab(i,j,k).isCovered() && (bc[n].lo(1) == BCType::ext_dir) )
                {
                    Real du_yl = 2.0*(state_fab(i,j  ,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));
                    Real du_yr = 2.0*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j  ,k,a_comp+n));
                    Real du_yc = (state_fab(i,j+1,k,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i,j-1,k,a_comp+n))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,a_sl_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
                }
                if ( (j == a_domain.bigEnd(1)) && !flag_fab(i,j,k).isCovered() && (bc[n].hi(1) == BCType::ext_dir) )
                {
                    Real du_yl = 2.0*(state_fab(i,j  ,k,a_comp+n) - state_fab(i,j-1,k,a_comp+n));
                    Real du_yr = 2.0*(state_fab(i,j+1,k,a_comp+n) - state_fab(i,j  ,k,a_comp+n));
                    Real du_yc = -(state_fab(i,j-1,k,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i,j+1,k,a_comp+n))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,a_sl_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
                }

#if (AMREX_SPACEDIM == 3)
                if ( (k == a_domain.smallEnd(2)) && !flag_fab(i,j,k).isCovered() && (bc[n].lo(2) == BCType::ext_dir) )
                {
                    Real du_zl = 2.0*(state_fab(i,j,k  ,a_comp+n) - state_fab(i,j,k-1,a_comp+n));
                    Real du_zr = 2.0*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k  ,a_comp+n));
                    Real du_zc = (state_fab(i,j,k+1,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i,j,k-1,a_comp+n))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,a_sl_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
                }
                if ( (k == a_domain.bigEnd(2)) && !flag_fab(i,j,k).isCovered() && (bc[n].hi(2) == BCType::ext_dir) )
                {
                    Real du_zl = 2.0*(state_fab(i,j,k  ,a_comp+n) - state_fab(i,j,k-1,a_comp+n));
                    Real du_zr = 2.0*(state_fab(i,j,k+1,a_comp+n) - state_fab(i,j,k  ,a_comp+n));
                    Real du_zc = -(state_fab(i,j,k-1,a_comp+n)+3.0*state_fab(i,j,k,a_comp+n)-4.0*state_fab(i,j,k+1,a_comp+n))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,a_sl_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
                }
#endif
            });

        }
    }
}

//
// Compute
//
void
Godunov::ComputeConvectiveTerm (MultiFab& a_state,
                                const int a_state_comp,
                                MultiFab& a_conv,
                                const int a_conv_comp,
                                const int a_ncomp,
                                D_DECL(MultiFab& a_fx,
                                       MultiFab& a_fy,
                                       MultiFab& a_fz),
                                D_DECL(MultiFab& edgst_x,
                                       MultiFab& edgst_y,
                                       MultiFab& edgst_z),
                                D_DECL(const MultiFab& a_umac,
                                       const MultiFab& a_vmac,
                                       const MultiFab& a_wmac),
                                D_DECL(const MultiFab& a_xsl,
                                       const MultiFab& a_ysl,
                                       const MultiFab& a_zsl),
                                const int a_sl_comp,
                                const Vector<BCRec>&  a_bcs,
                                const Geometry& a_geom,
                                int known_edgestate)
{
    AMREX_ALWAYS_ASSERT(a_state.hasEBFabFactory());
    AMREX_ALWAYS_ASSERT(a_state.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(a_bcs.size()==a_ncomp);

    // Compute fluxes
    ComputeFluxes( D_DECL(a_fx, a_fy, a_fz),
                   D_DECL(edgst_x, edgst_y, edgst_z),
                   a_state, a_state_comp, a_ncomp,
                   D_DECL(a_xsl, a_ysl, a_zsl), a_sl_comp,
                   D_DECL(a_umac, a_vmac, a_wmac),
                   a_geom, a_bcs, known_edgestate);

    // Compute divergence
    bool   already_on_centroids(true);
    int    nghost(2);
    MultiFab conv_tmp( a_conv.boxArray(), a_conv.DistributionMap(), a_ncomp, nghost,
                       MFInfo(), a_state.Factory() );

    conv_tmp.setVal(0.0);

    Array<MultiFab*,AMREX_SPACEDIM> fluxes;
    D_TERM(fluxes[0] = &a_fx;,
	   fluxes[1] = &a_fy;,
	   fluxes[2] = &a_fz;);

    EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes),
                         a_geom, already_on_centroids);

    // Redistribute divergence
    const amrex::MultiFab* weights;
    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_state.Factory());
    weights = &(ebfactory.getVolFrac());
    amrex::single_level_weighted_redistribute( 0, {conv_tmp}, {a_conv}, {*weights}, a_conv_comp, a_ncomp, {a_geom} );

    Gpu::synchronize();

    areaWeightFluxes(D_DECL(a_fx,a_fy,a_fz), a_ncomp, a_geom);
}

//
// Compute
//
void
Godunov::ComputeSyncConvectiveTerm (MultiFab& a_state,
                                    const int a_state_comp,
                                    MultiFab& a_conv,
                                    const int a_conv_comp,
                                    const int a_ncomp,
                                    D_DECL(MultiFab& a_fx,
                                           MultiFab& a_fy,
                                           MultiFab& a_fz),
                                    D_DECL(MultiFab& edgst_x,
                                           MultiFab& edgst_y,
                                           MultiFab& edgst_z),
                                    D_DECL(const MultiFab& a_umac,
                                           const MultiFab& a_vmac,
                                           const MultiFab& a_wmac),
                                    D_DECL(const MultiFab& a_ucorr,
                                           const MultiFab& a_vcorr,
                                           const MultiFab& a_wcorr),
                                    D_DECL(const MultiFab& a_xsl,
                                           const MultiFab& a_ysl,
                                           const MultiFab& a_zsl),
                                    const int a_sl_comp,
                                    const Vector<BCRec>&  a_bcs,
                                    const Geometry& a_geom,
                                    int known_edgestate)
{
    AMREX_ALWAYS_ASSERT(a_state.hasEBFabFactory());
    AMREX_ALWAYS_ASSERT(a_state.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(a_bcs.size()==a_ncomp);

    // Compute fluxes
    ComputeSyncFluxes( D_DECL(a_fx, a_fy, a_fz),
                       D_DECL(edgst_x, edgst_y, edgst_z),
                       a_state, a_state_comp, a_ncomp,
                       D_DECL(a_xsl, a_ysl, a_zsl), a_sl_comp,
                       D_DECL(a_umac, a_vmac, a_wmac),
                       D_DECL(a_ucorr, a_vcorr, a_wcorr),
                       a_geom, a_bcs, known_edgestate);

    // Compute divergence
    bool   already_on_centroids(true);
    int    nghost(2);
    MultiFab conv_tmp( a_conv.boxArray(), a_conv.DistributionMap(), a_ncomp, nghost,
                       MFInfo(), a_conv.Factory() );
    //
    // a_conv (== Vsync or Ssync) may already contain corrections, so need 2nd temp
    // for resistribution
    MultiFab conv_tmp_redist( a_conv.boxArray(), a_conv.DistributionMap(), a_ncomp, nghost,
                       MFInfo(), a_conv.Factory() );

    
    conv_tmp.setVal(0.0);

    Array<MultiFab*,AMREX_SPACEDIM> fluxes;
    D_TERM( fluxes[0] = &a_fx;,
	    fluxes[1] = &a_fy;,
	    fluxes[2] = &a_fz;);

    EB_computeDivergence(conv_tmp, GetArrOfConstPtrs(fluxes),
                         a_geom, already_on_centroids);

    // Redistribute divergence
    const amrex::MultiFab* weights;
    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_state.Factory());
    weights = &(ebfactory.getVolFrac());
    amrex::single_level_weighted_redistribute( 0, {conv_tmp}, {conv_tmp_redist}, {*weights}, 0, a_ncomp, {a_geom} );

    // Now add this update to a_conv
    // Note that redistribute() does not fill ghost cells, so if desired, would
    // need to fill ghosts here.
    MultiFab::Add(a_conv,conv_tmp_redist,0,a_conv_comp,a_ncomp,0);
    
    Gpu::synchronize();

    areaWeightFluxes(D_DECL(a_fx,a_fy,a_fz), a_ncomp, a_geom);
}

//
// create area-weighted fluxes, i.e. extensive fluxes
//
void
Godunov::areaWeightFluxes( D_DECL(MultiFab& a_fx,
				  MultiFab& a_fy,
				  MultiFab& a_fz),
			   const int        a_ncomp,
			   const Geometry&  a_geom)
{
    const Real*  dx = a_geom.CellSize();
    Real area[AMREX_SPACEDIM];
#if ( AMREX_SPACEDIM == 3 )
    area[0] = dx[1]*dx[2];
    area[1] = dx[0]*dx[2];
    area[2] = dx[0]*dx[1];
#else
    area[0] = dx[1];
    area[1] = dx[0];
#endif

    const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_fx.Factory());
    const FabArray<EBCellFlagFab>& flagMF = ebfactory.getMultiEBCellFlagFab();
    
    Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
    areafrac  = ebfactory.getAreaFrac();
    
    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_fx,mfi_info); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox(IntVect(AMREX_D_DECL(0,0,0)));
	
	D_TERM( const auto& fx = a_fx.array(mfi);,
		const auto& fy = a_fy.array(mfi);,
		const auto& fz = a_fz.array(mfi););

	D_TERM( const Box ubx = mfi.tilebox(IntVect(AMREX_D_DECL(1,0,0)));,
		const Box vbx = mfi.tilebox(IntVect(AMREX_D_DECL(0,1,0)));,
		const Box wbx = mfi.tilebox(IntVect(AMREX_D_DECL(0,0,1))););
	
	AMREX_FOR_4D(ubx, a_ncomp, i, j, k, n, {fx(i,j,k,n) *= area[0];});
	AMREX_FOR_4D(vbx, a_ncomp, i, j, k, n, {fy(i,j,k,n) *= area[1];});
#if (AMREX_SPACEDIM==3)
	AMREX_FOR_4D(wbx, a_ncomp, i, j, k, n, {fz(i,j,k,n) *= area[2];});
#endif
	
	//
	// Deal with irregular cells 
	//
        auto fabtyp = flagMF[mfi].getType(bx);
	
	if (fabtyp == FabType::regular ) continue;
	
 	D_TERM( const auto& afrac_x = areafrac[0]->array(mfi);,
		const auto& afrac_y = areafrac[1]->array(mfi);,
		const auto& afrac_z = areafrac[2]->array(mfi););
	  
	if ( fabtyp != FabType::covered )
        {
	  //
	  // Account for "effective areas" for cut cells
	  //
	  AMREX_FOR_4D(ubx, a_ncomp, i, j, k, n, {fx(i,j,k,n) *= afrac_x(i,j,k);});
	  AMREX_FOR_4D(vbx, a_ncomp, i, j, k, n, {fy(i,j,k,n) *= afrac_y(i,j,k);});
#if (AMREX_SPACEDIM==3)
	  AMREX_FOR_4D(wbx, a_ncomp, i, j, k, n, {fz(i,j,k,n) *= afrac_z(i,j,k);});
#endif
        }
    }

    D_TERM( a_fx.FillBoundary(a_geom.periodicity());,
	    a_fy.FillBoundary(a_geom.periodicity());,
	    a_fz.FillBoundary(a_geom.periodicity()); );
}

#endif
