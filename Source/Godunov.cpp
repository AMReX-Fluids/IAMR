
//
// Godunov is the object which calculates advective terms for iamr.
//

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <Godunov.H>
#include <GODUNOV_F.H>

#include <algorithm>

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
Godunov::maxchng_velmag ( MultiFab const&  u_old,
			  MultiFab const&  u_new )
{
    Real max_change = 0.0;

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#if !defined(AMREX_USE_GPU) && defined(_OPENMP)
#pragma omp parallel  if (!system::regtest_reduction) reduction(max:max_change)
#endif
    for (MFIter mfi(u_old,true); mfi.isValid(); ++mfi)
    {
        const auto& bx   = mfi.tilebox();
        const auto& uold = u_old[mfi].array();
        const auto& unew = u_new[mfi].array();

        reduce_op.eval(bx, reduce_data, [uold, unew]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real uold_mag = 0.0;
            Real unew_mag = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                uold_mag += uold(i,j,k,d)*uold(i,j,k,d);
                unew_mag += unew(i,j,k,d)*unew(i,j,k,d);
            }

            uold_mag = std::sqrt(uold_mag);
            unew_mag = std::sqrt(unew_mag);

            return std::abs(unew_mag-uold_mag);
        });

        max_change = std::max(amrex::get<0>(reduce_data.value()),
                              max_change);
    }

    ParallelDescriptor::ReduceRealMax(max_change);

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
// Compute total source term for velocities, weighted by rho.
//
// tforces = (tforces + visc - gp)/rho
//
void
Godunov::Sum_tf_gp_visc ( Box const&                 bx,
                          Array4<Real>        const& tforces,
                          Array4<Real const>  const& visc,
                          Array4<Real const>  const& gp,
                          Array4<Real const>  const& rho ) const
{
    amrex::ParallelFor(bx, AMREX_SPACEDIM, [tforces, visc, gp, rho]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        tforces(i,j,k,n) = ( tforces(i,j,k,n) + visc(i,j,k,n)
                             - gp(i,j,k,n) ) / rho(i,j,k);
    });
}


void
Godunov::Sum_tf_divu_visc ( amrex::Box                        const& bx,
                            amrex::Array4<amrex::Real>        const& tforces,
                            amrex::Array4<amrex::Real const>  const& visc,
                            amrex::Array4<amrex::Real const>  const& divu,
                            amrex::Array4<amrex::Real const>  const& S,
                            amrex::Array4<amrex::Real const>  const& rho,
                            int ncomp, int iconserv ) const
{

    if (iconserv == 1)
    {
        amrex::ParallelFor(bx, ncomp, [tforces, visc, S, divu]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            tforces(i,j,k,n) = tforces(i,j,k,n) + visc(i,j,k,n)
                                 - S(i,j,k,n) * divu(i,j,k);
        });
    }
    else
    {
        amrex::ParallelFor(bx, ncomp, [tforces, visc, rho]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            tforces(i,j,k,n) = ( tforces(i,j,k,n) + visc(i,j,k,n) ) / rho(i,j,k);
        });
    }

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
