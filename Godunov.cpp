
//
// $Id: Godunov.cpp,v 1.26 2001-08-09 22:42:00 marc Exp $
//

//
// Godunov is the object which calculates advective terms for iamr.
//
#include <winstd.H>

#include <LO_BCTYPES.H>
#include <CoordSys.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <FArrayBox.H>
#include <Godunov.H>
#include <GODUNOV_F.H>

#include <algorithm>

#define GEOM_GROW 1
#define HYP_GROW  3
#define XVEL 0
#define YVEL 1
#define ZVEL 2

int Godunov::verbose               = 0;
int Godunov::slope_order           = 4;
int Godunov::use_forces_in_trans   = 0;
const int use_unlimited_slopes_DEF = 0;

//
// Construct the Godunov Object.
//

Godunov::Godunov ()
    :
    max_1d(0)
{
    read_params();
    //
    // Set to 512 initiallly.
    //
    ZeroScratch();
    SetScratch(512);
}

//
// Size the 1D workspace explicitly.
//

Godunov::Godunov (int max_size)
    :
    max_1d(max_size)
{
    read_params();
    ZeroScratch();
    SetScratch(max_size);
}

//
// Read parameters from input file and command line.
//

void
Godunov::read_params ()
{
    //
    // Read parameters from input file and command line.
    //
    ParmParse pp("godunov");

    pp.query("v",verbose);

    pp.query("slope_order",slope_order);
#if (BL_SPACEDIM==2)
    BL_ASSERT(slope_order==1 || slope_order==2 || slope_order==4);
#else
    BL_ASSERT(slope_order==1 || slope_order==4);
#endif
    pp.query("use_forces_in_trans",use_forces_in_trans);
    int use_unlimited_slopes=use_unlimited_slopes_DEF;
    pp.query("use_unlimited_slopes",use_unlimited_slopes);

    FORT_SET_PARAMS(slope_order,use_unlimited_slopes);
}

//
// Set 1d scratch space as empty.
//

void
Godunov::ZeroScratch ()
{
    D_TERM(stxlo=0;,  stylo=0;,  stzlo=0;);
    D_TERM(stxhi=0;,  styhi=0;,  stzhi=0;);
    D_TERM(slxscr=0;, slyscr=0;, slzscr=0;);
}

//
// Initialize 1d scratch space to a bogus value.
//

void
Godunov::SetBogusScratch ()
{
#ifndef NDEBUG
    const Real bogus_value = 1.e20;

    for (int i = 0 ; i < scr_size ; i++)
    {
        D_TERM(stxlo[i]=bogus_value;,
               stylo[i]=bogus_value;,
               stzlo[i]=bogus_value;);

        D_TERM(stxhi[i]=bogus_value;,
               styhi[i]=bogus_value;,
               stzhi[i]=bogus_value;);

        D_TERM(slxscr[i]=bogus_value;,
               slyscr[i]=bogus_value;,
               slzscr[i]=bogus_value;);
    }
#endif /*NDEBUG*/
}

//
// Set 1d scratch space.
//

void
Godunov::SetScratch (int max_size)
{
    //
    // Set sizing parameters.
    //
    if (max_size <= max_1d)
        return;
    else
        max_1d = std::max(max_1d,max_size);
    scr_size = (max_size+2*HYP_GROW)*4;
    //
    // Get rid of the old scratch space.
    //
    RemScratch();
    //
    // Construct arrays.
    //
    D_TERM(stxlo  = new Real[scr_size];,
           stylo  = new Real[scr_size];,
           stzlo  = new Real[scr_size];);

    D_TERM(stxhi  = new Real[scr_size];,
           styhi  = new Real[scr_size];,
           stzhi  = new Real[scr_size];);

    D_TERM(slxscr = new Real[scr_size];,
           slyscr = new Real[scr_size];,
           slzscr = new Real[scr_size];);
}

//
// Remove 1D scratch space.
//

void
Godunov::RemScratch ()
{
    D_TERM(delete [] stxlo;,  delete [] stylo;,  delete [] stzlo;);
    D_TERM(delete [] stxhi;,  delete [] styhi;,  delete [] stzhi;);
    D_TERM(delete [] slxscr;, delete [] slyscr;, delete [] slzscr;);
}

//
// Destructor destroys work arrays.
//

Godunov::~Godunov ()
{
    RemScratch();
}

//
// Set up the farrayboxes for computing edge states, This function
// returns, the Farrayboxes where the fluxes are stored as these
// are used in NavierStokes.
//
// It also computes the transverse advective velocities.
//
// The amount of workspace needed in FArrayBox work is currently 2*SDIM+1.
//

void
Godunov::Setup (const Box&       grd,
                const Real*      dx,
                Real             dt,
                int              velpred,
                FArrayBox&       xflux,
                const int*       ubc,
                FArrayBox&       yflux,
                const int*       vbc,
#if (BL_SPACEDIM == 3)
                FArrayBox&       zflux,
                const int*       wbc,
#endif
                FArrayBox&       U,
                FArrayBox&       rho, 
                const FArrayBox& tforces)
{
    BL_ASSERT(rho.nComp() == 1);
    BL_ASSERT(U.nComp() >= BL_SPACEDIM);
    //
    // Compute the edge boxes.
    //
    D_TERM(xflux_bx = grd; xflux_bx.surroundingNodes(0);,
           yflux_bx = grd; yflux_bx.surroundingNodes(1);,
           zflux_bx = grd; zflux_bx.surroundingNodes(2););
    //
    // Create storage for fluxes.
    //
    if (!velpred)
    {
        D_TERM(xflux.resize(xflux_bx,1);,
               yflux.resize(yflux_bx,1);,
               zflux.resize(zflux_bx,1););
    }
    //
    // Ensure 1D scratch space is large enough.
    //
    SetScratch(BoxLib::grow(grd,HYP_GROW).longside());
    //
    // Create the advective velocities and FAB workspace for GODUNOV Box.
    //
    work_bx = BoxLib::grow(grd,1);
    work.resize(work_bx,2*BL_SPACEDIM+1);
    D_TERM(uad.resize(work_bx,1);,
           vad.resize(work_bx,1);,
           wad.resize(work_bx,1););

    SetBogusScratch();
    //
    // Test the cell-centered velocities.
    //
    // Real u_max[3];
    // test_u_rho( S, rho, grd, dx, dt, u_max );
    //
    // Create the bounds and pointers.
    //
    const int *lo       = grd.loVect();
    const int *hi       = grd.hiVect();
    const int *u_lo     = U.loVect();
    const int *u_hi     = U.hiVect();
    const int *w_lo     = work.loVect();
    const int *w_hi     = work.hiVect();
    const Real *u_dat   = U.dataPtr(XVEL);
    const Real *uad_dat = uad.dataPtr();
    const Real *v_dat   = U.dataPtr(YVEL);
    const Real *vad_dat = vad.dataPtr();
#if (BL_SPACEDIM == 3)
    const Real *w_dat   = U.dataPtr(ZVEL);
    const Real *wad_dat = wad.dataPtr();
    const Real *xhi_dat = work.dataPtr(1);
    const Real *yhi_dat = work.dataPtr(2);
    const Real *zhi_dat = work.dataPtr(3);
    const Real *slx_dat = work.dataPtr(4);
    const Real *sly_dat = work.dataPtr(5);
    const Real *slz_dat = work.dataPtr(6);
#else
    const Real *xhi_dat = work.dataPtr(1);
    const Real *yhi_dat = work.dataPtr(2);
    const Real *slx_dat = work.dataPtr(3);
    const Real *sly_dat = work.dataPtr(4);
#endif 
    const Real* tforcedat = 0;

    if (use_forces_in_trans)
    {
        tforcedat = tforces.dataPtr();
    }
    //
    // Compute the transverse velocities.
    //
    FORT_TRANSVEL(u_dat, uad_dat, xhi_dat, slx_dat, ubc, slxscr, 
                  v_dat, vad_dat, yhi_dat, sly_dat, vbc, slyscr, 
#if (BL_SPACEDIM == 3)
                  w_dat, wad_dat, zhi_dat, slz_dat, wbc, slzscr,
#endif    
                  ARLIM(u_lo), ARLIM(u_hi),
                  ARLIM(w_lo), ARLIM(w_hi), 
                  lo, hi, &dt, dx, &use_forces_in_trans, tforcedat);
}

//
// Advection functions follow.
//

//
// Compute the edge states using the advective transverse velocities
// The amount of workspace needed in FArrayBox work is currently 2*SDIM+1.
//

void
Godunov::edge_states (const Box&  grd,
                      const Real* dx,
                      Real        dt,
                      int         velpred,
                      FArrayBox&  uedge,
                      FArrayBox&  stx,
                      FArrayBox&  vedge,
                      FArrayBox&  sty,
#if (BL_SPACEDIM == 3)               
                      FArrayBox&  wedge,
                      FArrayBox&  stz,
#endif
                      FArrayBox&  U,
                      FArrayBox&  S,
                      FArrayBox&  tforces,
                      int         fab_ind,
                      int         state_ind,
                      const int*  bc)
{
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(U.nComp()       >= BL_SPACEDIM);
    BL_ASSERT(S.nComp()       >= fab_ind    );
    BL_ASSERT(tforces.nComp() >= fab_ind    );

    BL_ASSERT(uedge.box()     == xflux_bx   );
    BL_ASSERT(uedge.nComp()   >= 1          );
    BL_ASSERT(stx.nComp()     >= 1          );

    BL_ASSERT(vedge.box()     == yflux_bx   );
    BL_ASSERT(vedge.nComp()   >= 1          );
    BL_ASSERT(sty.nComp()     >= 1          );
#if (BL_SPACEDIM == 3)
    BL_ASSERT(wedge.box()     == zflux_bx   );
    BL_ASSERT(wedge.nComp()   >= 1          );
    BL_ASSERT(stz.nComp()     >= 1          );
#endif    
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *s_lo       = S.loVect();
    const int *s_hi       = S.hiVect();
    const int *u_lo       = uedge.loVect();
    const int *u_hi       = uedge.hiVect();
    const int *v_lo       = vedge.loVect();
    const int *v_hi       = vedge.hiVect();
    const int *ww_lo      = work.loVect();
    const int *ww_hi      = work.hiVect();
    const Real *s_dat     = S.dataPtr(fab_ind);
    const Real *u_dat     = U.dataPtr(XVEL);
    const Real *v_dat     = U.dataPtr(YVEL);
    const Real *tfr_dat   = tforces.dataPtr(fab_ind);
    const Real *uad_dat   = uad.dataPtr();
    const Real *vad_dat   = vad.dataPtr();
    const Real *uedge_dat = uedge.dataPtr();
    const Real *vedge_dat = vedge.dataPtr();
    const Real *stx_dat   = stx.dataPtr();
    const Real *sty_dat   = sty.dataPtr();
    //
    // Set work space to bogus values.
    //
    SetBogusScratch();

#if (BL_SPACEDIM == 3)
    const int *w_lo       = wedge.loVect();
    const int *w_hi       = wedge.hiVect();
    const Real *w_dat     = U.dataPtr(ZVEL);
    const Real *wad_dat   = wad.dataPtr();
    const Real *wedge_dat = wedge.dataPtr();
    const Real *stz_dat   = stz.dataPtr();
    const Real *xhi_dat   = work.dataPtr(0);
    const Real *yhi_dat   = work.dataPtr(0);
    const Real *zhi_dat   = work.dataPtr(0);
    const Real *xlo_dat   = work.dataPtr(1);
    const Real *ylo_dat   = work.dataPtr(2);
    const Real *zlo_dat   = work.dataPtr(3);
    const Real *slx_dat   = work.dataPtr(4);
    const Real *sly_dat   = work.dataPtr(5);
    const Real *slz_dat   = work.dataPtr(6);
#else
    const Real *xhi_dat   = work.dataPtr(0);
    const Real *yhi_dat   = work.dataPtr(0);
    const Real *xlo_dat   = work.dataPtr(1);
    const Real *ylo_dat   = work.dataPtr(2);
    const Real *slx_dat   = work.dataPtr(3);
    const Real *sly_dat   = work.dataPtr(4);
#endif
    //
    // C component indices starts from 0, Fortran from 1
    //
    int fort_ind = state_ind+1;  

    FORT_ESTATE(s_dat, tfr_dat, ARLIM(s_lo), ARLIM(s_hi),

                u_dat, xlo_dat, xhi_dat, slx_dat, uad_dat,
                slxscr, stxlo, stxhi,
                uedge_dat, stx_dat, ARLIM(u_lo), ARLIM(u_hi), 

                v_dat, ylo_dat, yhi_dat, sly_dat, vad_dat,
                slyscr, stylo, styhi,
                vedge_dat, sty_dat, ARLIM(v_lo), ARLIM(v_hi), 
#if (BL_SPACEDIM == 3)
                w_dat, zlo_dat, zhi_dat, slz_dat, wad_dat,
                slzscr, stzlo, stzhi,
                wedge_dat, stz_dat, ARLIM(w_lo), ARLIM(w_hi), 
#endif
                ARLIM(ww_lo), ARLIM(ww_hi),
                bc, lo, hi, &dt, dx, &fort_ind, &velpred, 
                &use_forces_in_trans);
}

//
// Compute the edge states for The Mac projection.
// FArrayBox work sized as in edge_states.
//

void
Godunov::ComputeUmac (const Box&  grd,
                      const Real* dx,
                      Real        dt, 
                      FArrayBox&  umac,
                      const int*  ubc, 
                      FArrayBox&  vmac,
                      const int*  vbc, 
#if (BL_SPACEDIM == 3)
                      FArrayBox&  wmac,
                      const int*  wbc, 
#endif
                      FArrayBox&  U,
                      FArrayBox&  tforces)
{
    int velpred = 1;
    //
    // 2D calls.
    //
#if (BL_SPACEDIM == 2)                  
    edge_states(grd, dx, dt, velpred,
                umac, umac,
                vmac, vmac,
                U, U, tforces, XVEL, XVEL, ubc);
    edge_states(grd, dx, dt, velpred,
                umac, umac,
                vmac, vmac,
                U, U, tforces, YVEL, YVEL, vbc);
#endif
    //
    // 3D calls.
    //
#if (BL_SPACEDIM == 3)                  
    edge_states(grd, dx, dt, velpred,
                umac, umac,
                vmac, vmac,
                wmac, wmac,
                U, U, tforces, XVEL, XVEL, ubc);
    edge_states(grd, dx, dt, velpred,
                umac, umac,
                vmac, vmac,
                wmac, wmac,
                U, U, tforces, YVEL, YVEL, vbc);
    edge_states(grd, dx, dt, velpred,
                umac, umac,
                vmac, vmac,
                wmac, wmac,
                U, U, tforces, ZVEL, ZVEL, wbc);
#endif
}

//
// Advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
// FArrayBox work sized as in edge_states.
//

void
Godunov::AdvectState (const Box&  grd,
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
                      int         fab_ind,
                      FArrayBox&  aofs,
                      int         aofs_ind,
                      int         iconserv,
                      int         state_ind,
                      const int*  bc,
                      FArrayBox&  vol)
{
    int velpred = 0;
    //
    // Compute edge states for an advected quantity.
    //
    edge_states(grd, dx, dt, velpred,
                uedge, xflux,
                vedge, yflux,
#if (BL_SPACEDIM == 3)             
                wedge, zflux,
#endif
                U, S, tforces, fab_ind, state_ind, bc);
    //
    // Compute the advective tendency.
    //
    ComputeAofs( grd,
                 areax, uedge, xflux,  
                 areay, vedge, yflux,  
#if (BL_SPACEDIM == 3)                             
                 areaz, wedge, zflux,
#endif                     
                 vol, aofs, aofs_ind, iconserv);
}

//
// Compute the advective derivative from fluxes.
//

void
Godunov::ComputeAofs (const Box& grd, 
                      FArrayBox& areax,
                      FArrayBox& uedge,
                      FArrayBox& xflux,  
                      FArrayBox& areay,
                      FArrayBox& vedge,
                      FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                               
                      FArrayBox& areaz,
                      FArrayBox& wedge,
                      FArrayBox& zflux,
#endif
                      FArrayBox& vol,
                      FArrayBox& aofs,
                      int        aofs_ind,
                      int        iconserv )
{
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *a_lo       = aofs.loVect();
    const int *a_hi       = aofs.hiVect();
    const Real *aofs_dat  = aofs.dataPtr(aofs_ind);
    const int *ax_lo      = areax.loVect();
    const int *ax_hi      = areax.hiVect();
    const Real *areax_dat = areax.dataPtr();
    const int *ay_lo      = areay.loVect();
    const int *ay_hi      = areay.hiVect();
    const Real *areay_dat = areay.dataPtr();
    const int *v_lo       = vol.loVect();
    const int *v_hi       = vol.hiVect();
    const Real *vol_dat   = vol.dataPtr();
    const int *xflux_lo   = xflux.loVect();
    const int *xflux_hi   = xflux.hiVect();
    const Real *xflux_dat = xflux.dataPtr();
    const Real *uedge_dat = uedge.dataPtr();
    const int *yflux_lo   = yflux.loVect();
    const int *yflux_hi   = yflux.hiVect();
    const Real *yflux_dat = yflux.dataPtr();
    const Real *vedge_dat = vedge.dataPtr();
#if (BL_SPACEDIM == 3)
    const int *az_lo      = areaz.loVect();
    const int *az_hi      = areaz.hiVect();
    const Real *areaz_dat = areaz.dataPtr();
    const int *zflux_lo   = zflux.loVect();
    const int *zflux_hi   = zflux.hiVect();
    const Real *zflux_dat = zflux.dataPtr();
    const Real *wedge_dat = wedge.dataPtr();
#endif
    //
    // Compute the advective tendency.
    //
    FORT_ADV_FORCING(aofs_dat, ARLIM(a_lo), ARLIM(a_hi),
                     xflux_dat, uedge_dat, ARLIM(xflux_lo), ARLIM(xflux_hi),
                     areax_dat, ARLIM(ax_lo), ARLIM(ax_hi),
                     yflux_dat, vedge_dat, ARLIM(yflux_lo), ARLIM(yflux_hi),
                     areay_dat, ARLIM(ay_lo), ARLIM(ay_hi),
#if (BL_SPACEDIM == 3)                                                    
                     zflux_dat, wedge_dat, ARLIM(zflux_lo), ARLIM(zflux_hi),
                     areaz_dat, ARLIM(az_lo), ARLIM(az_hi),
#endif
                     vol_dat, ARLIM(v_lo), ARLIM(v_hi),
                     lo, hi, &iconserv);
}

//
// Sync advect a state component.
// This routine assumes uad,vad,wad have been precomputed.
//

void
Godunov::SyncAdvect (const Box&  grd,
                     const Real* dx,
                     Real        dt,
                     int         level,
                     FArrayBox& areax,
                     FArrayBox& uedge,
                     FArrayBox& ucorr,
                     FArrayBox& xflux,
                     FArrayBox& areay,
                     FArrayBox& vedge,
                     FArrayBox& vcorr,
                     FArrayBox& yflux,
#if (BL_SPACEDIM == 3)
                     FArrayBox& areaz,
                     FArrayBox& wedge,
                     FArrayBox& wcorr,
                     FArrayBox& zflux,
#endif
                     FArrayBox& U,
                     FArrayBox& S,
                     FArrayBox& tforces,
                     int        fab_ind,
                     FArrayBox& sync,
                     int        sync_ind,
                     int        iconserv,
                     int        state_ind,
                     const int* bc,
                     FArrayBox& vol)
{
    int velpred = 0;
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(S.nComp()       >= BL_SPACEDIM);
    BL_ASSERT(S.nComp()       >= fab_ind    );
    BL_ASSERT(tforces.nComp() >= fab_ind    );
    BL_ASSERT(sync.nComp()    >= sync_ind   );

    BL_ASSERT(ucorr.box()     == xflux_bx   );
    BL_ASSERT(ucorr.nComp()   >= 1          );

    BL_ASSERT(vcorr.box()     == yflux_bx   );
    BL_ASSERT(vcorr.nComp()   >= 1          );
#if (BL_SPACEDIM == 3)
    BL_ASSERT(wcorr.box()     == zflux_bx   );
    BL_ASSERT(wcorr.nComp()   >= 1          );
#endif    
    //
    // Compute the edge states.
    //
    edge_states(grd, dx, dt, velpred,
                uedge, xflux,
                vedge, yflux,
#if (BL_SPACEDIM == 3)     
                wedge, zflux,
#endif
                U, S, tforces, fab_ind, state_ind, bc);

    //
    // Compute the advective tendency for the mac sync.
    //
    ComputeSyncAofs(grd,
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
                          FArrayBox& areax,
                          FArrayBox& ucorr,
                          FArrayBox& xflux,  
                          FArrayBox& areay,
                          FArrayBox& vcorr,
                          FArrayBox& yflux,  
#if (BL_SPACEDIM == 3)                             
                          FArrayBox& areaz,
                          FArrayBox& wcorr,
                          FArrayBox& zflux,
#endif                     
                          FArrayBox& vol,
                          FArrayBox& sync,
                          int        sync_ind,
                          int        iconserv)
{
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *s_lo       = sync.loVect();
    const int *s_hi       = sync.hiVect();
    const Real *sync_dat  = sync.dataPtr(sync_ind);
    const int *ax_lo      = areax.loVect();
    const int *ax_hi      = areax.hiVect();
    const Real *areax_dat = areax.dataPtr();
    const int *ay_lo      = areay.loVect();
    const int *ay_hi      = areay.hiVect();
    const Real *areay_dat = areay.dataPtr();
    const int *v_lo       = vol.loVect();
    const int *v_hi       = vol.hiVect();
    const Real *vol_dat   = vol.dataPtr();
    const int *xflux_lo   = xflux.loVect();
    const int *xflux_hi   = xflux.hiVect();
    const Real *xflux_dat = xflux.dataPtr();
    const Real *ucorr_dat = ucorr.dataPtr();
    const int *yflux_lo   = yflux.loVect();
    const int *yflux_hi   = yflux.hiVect();
    const Real *yflux_dat = yflux.dataPtr();
    const Real *vcorr_dat = vcorr.dataPtr();
#if (BL_SPACEDIM == 3)
    const int *az_lo      = areaz.loVect();
    const int *az_hi      = areaz.hiVect();
    const Real *areaz_dat = areaz.dataPtr();
    const int *zflux_lo   = zflux.loVect();
    const int *zflux_hi   = zflux.hiVect();
    const Real *zflux_dat = zflux.dataPtr();
    const Real *wcorr_dat = wcorr.dataPtr();
#endif
    //
    // Compute the corrective tendency.
    //
    FORT_SYNC_ADV_FORCING(sync_dat, ARLIM(s_lo), ARLIM(s_hi),
                           
                          xflux_dat, ucorr_dat, ARLIM(xflux_lo), ARLIM(xflux_hi),
                          areax_dat, ARLIM(ax_lo), ARLIM(ax_hi),
                           
                          yflux_dat, vcorr_dat, ARLIM(yflux_lo), ARLIM(yflux_hi),
                          areay_dat, ARLIM(ay_lo), ARLIM(ay_hi),
#if (BL_SPACEDIM == 3)                                             
                          zflux_dat, wcorr_dat, ARLIM(zflux_lo), ARLIM(zflux_hi),
                          areaz_dat, ARLIM(az_lo), ARLIM(az_hi),
#endif
                          vol_dat, ARLIM(v_lo), ARLIM(v_hi),
                          lo, hi, &iconserv);
}    

//
// Correct a scalar for under-over shoots.
//

void
Godunov::ScalMinMax (FArrayBox& Sold,
                     FArrayBox& Snew,
                     int        ind, 
                     const int* bc,
                     const Box& grd)
{
    const int *slo       = Sold.loVect();
    const int *shi       = Sold.hiVect();
    const int *lo        = grd.loVect();
    const int *hi        = grd.hiVect();
    const Real *Sold_dat = Sold.dataPtr(ind);
    const Real *Snew_dat = Snew.dataPtr(ind);

#if (BL_SPACEDIM == 3)
    Box flatbox(grd);
    int zlen = flatbox.length()[BL_SPACEDIM-1];
    flatbox.growHi(BL_SPACEDIM-1,3-zlen);
    FArrayBox smin(flatbox,1);
    FArrayBox smax(flatbox,1);
    const Real *smin_dat = smin.dataPtr();
    const Real *smax_dat = smax.dataPtr(); 
#endif

    FORT_SCALMINMAX (Sold_dat, ARLIM(slo), ARLIM(shi),
                     Snew_dat, ARLIM(slo), ARLIM(shi),
#if (BL_SPACEDIM == 3)
                     smin_dat, smax_dat,
                     ARLIM(lo), ARLIM(hi),
#endif
                     lo, hi, bc);
}

//
// Diagnostic functions follow
//

//
// Estimate the maximum allowable timestep at a cell center.
//

Real
Godunov::estdt (FArrayBox&  U,
                FArrayBox&  tforces,
                FArrayBox&  rho,
                const Box&  grd,
                const Real* dx,
                Real        cfl,
                Real*       u_max)
{
    BL_ASSERT( U.nComp()       >= BL_SPACEDIM );
    BL_ASSERT( tforces.nComp() >= BL_SPACEDIM );
    BL_ASSERT( rho.nComp()     == 1           );

    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const int *vlo    = U.loVect();
    const int *vhi    = U.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *Udat  = U.dataPtr();
    const Real *tfdat = tforces.dataPtr();
    const Real *rdat  = rho.dataPtr();

    Real dt;
    FORT_ESTDT(Udat,  ARLIM(vlo), ARLIM(vhi),
               tfdat, ARLIM(tlo), ARLIM(thi),
               rdat,  ARLIM(rlo), ARLIM(rhi),
               lo, hi, &dt, dx, &cfl, u_max);
    return dt;
}

//
// Estimate the extrema of cell-centered us and rho.
//

Real
Godunov::test_u_rho (FArrayBox&  U,
                     FArrayBox&  rho, 
                     const Box&  grd,
                     const Real* dx,
                     const Real  dt,
                     const Real* u_max)
{
    BL_ASSERT(U.nComp()   >= BL_SPACEDIM);
    BL_ASSERT(rho.nComp() == 1          );
    
    const int *lo  = grd.loVect();
    const int *hi  = grd.hiVect();
    const int *vlo = U.loVect();
    const int *vhi = U.hiVect();
    const int *rlo = rho.loVect();
    const int *rhi = rho.hiVect();
    const Real *rh = rho.dataPtr();
    D_TERM(const Real *u  = U.dataPtr(XVEL);,
           const Real *v  = U.dataPtr(YVEL);,
           const Real *w  = U.dataPtr(ZVEL););

    Real cflmax = 0;
    FORT_TEST_U_RHO(u,  ARLIM(vlo), ARLIM(vhi),
                    v,  ARLIM(vlo), ARLIM(vhi),
#if (BL_SPACEDIM == 3)                          
                    w,  ARLIM(vlo), ARLIM(vhi),
#endif
                    rh, ARLIM(rlo), ARLIM(rhi),
                    lo, hi, &dt, dx, &cflmax, u_max, &verbose);
    return cflmax;
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
    FORT_TEST_UMAC_RHO(um, ARLIM(ulo), ARLIM(uhi),
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
Godunov::Add_tf (FArrayBox& Sold,
                 FArrayBox& Snew,
                 int        start_ind,
                 int        num_comp, 
                 FArrayBox& tforces,
                 int        tf_ind,
                 const Box& grd,
                 Real       dt)
{
    BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
    BL_ASSERT(tforces.nComp() >= tf_ind    + num_comp);

    const int *slo    = Sold.loVect();
    const int *shi    = Sold.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SOdat = Sold.dataPtr(start_ind);
    const Real *SNdat = Snew.dataPtr(start_ind);
    const Real *TFdat = tforces.dataPtr(tf_ind);
    
    FORT_UPDATE_TF(SOdat, ARLIM(slo), ARLIM(shi), 
                   SNdat, ARLIM(slo), ARLIM(shi),
                   TFdat, ARLIM(tlo), ARLIM(thi),
                   lo, hi, &dt, &num_comp);
}

//
// Correct the 1st order RK to 2nd order via
//
// psi^n+1 = psi^* + (dt/2)*(tf^* - tf^n)
//

void
Godunov::Correct_tf (FArrayBox& Sstar,
                     FArrayBox& Snp1,
                     int        start_ind,
                     int        num_comp, 
                     FArrayBox& tfstar,
                     FArrayBox& tfn,
                     int        tf_ind,
                     const Box& grd,
                     Real       dt)
{
    BL_ASSERT(Snp1.nComp()   >= start_ind + num_comp);
    BL_ASSERT(Sstar.nComp()  >= start_ind + num_comp);
    BL_ASSERT(tfstar.nComp() >= tf_ind    + num_comp);
    BL_ASSERT(tfn.nComp()    >= tf_ind    + num_comp);

    const int *slo    = Sstar.loVect();
    const int *shi    = Sstar.hiVect();
    const int *tlo    = tfstar.loVect();
    const int *thi    = tfstar.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SSdat = Sstar.dataPtr(start_ind);
    const Real *SPdat = Snp1.dataPtr(start_ind);
    const Real *TSdat = tfstar.dataPtr(tf_ind);
    const Real *TNdat = tfn.dataPtr(tf_ind);
    
    FORT_CORRECT_TF(SSdat, SPdat, ARLIM(slo), ARLIM(shi),
                    TSdat, TNdat, ARLIM(tlo), ARLIM(thi),
                    lo, hi, &dt, &num_comp);
}

//
// Compute the update rule
//
// psi^n+1 = psi^n - dt*aofs + dt*tforces
//

void
Godunov::Add_aofs_tf (FArrayBox& Sold,
                      FArrayBox& Snew,
                      int        start_ind,
                      int        num_comp,
                      FArrayBox& Aofs,
                      int        aofs_ind,
                      FArrayBox& tforces,
                      int        tf_ind,
                      const Box& grd,
                      Real       dt)
{
    BL_ASSERT(Snew.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Sold.nComp()    >= start_ind + num_comp);
    BL_ASSERT(Aofs.nComp()    >= aofs_ind  + num_comp);
    BL_ASSERT(tforces.nComp() >= tf_ind    + num_comp);

    const int *slo    = Sold.loVect();
    const int *shi    = Sold.hiVect();
    const int *alo    = Aofs.loVect();
    const int *ahi    = Aofs.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SOdat = Sold.dataPtr(start_ind);
    const Real *SNdat = Snew.dataPtr(start_ind);
    const Real *AOdat = Aofs.dataPtr(aofs_ind);
    const Real *TFdat = tforces.dataPtr(tf_ind);
    
    FORT_UPDATE_AOFS_TF(SOdat, ARLIM(slo), ARLIM(shi), 
                        SNdat, ARLIM(slo), ARLIM(shi),
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
Godunov::Add_aofs_tf_gp (FArrayBox& Uold,
                         FArrayBox& Unew,
                         FArrayBox& Aofs,
                         FArrayBox& tforces,
                         FArrayBox& gp,
                         FArrayBox& rho, 
                         const Box& grd,
                         Real       dt)
{
    BL_ASSERT(Unew.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(Uold.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(Aofs.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(tforces.nComp() >= BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      == BL_SPACEDIM);
    BL_ASSERT(rho.nComp()     == 1          );
    
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const int *ulo    = Uold.loVect();
    const int *uhi    = Uold.hiVect();
    const int *alo    = Aofs.loVect();
    const int *ahi    = Aofs.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *UOdat = Uold.dataPtr();
    const Real *UNdat = Unew.dataPtr();
    const Real *AOdat = Aofs.dataPtr();
    const Real *TFdat = tforces.dataPtr();
    const Real *GPdat = gp.dataPtr();
    const Real *RHdat = rho.dataPtr();
    
    FORT_UPDATE_AOFS_TF_GP(UOdat, ARLIM(ulo), ARLIM(uhi),
                           UNdat, ARLIM(ulo), ARLIM(uhi),
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
Godunov::Sum_tf_gp (FArrayBox& tforces, 
                    FArrayBox& gp,
                    FArrayBox& rho)
{
    BL_ASSERT(rho.nComp()     == 1);
    BL_ASSERT(tforces.nComp() >= BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      == BL_SPACEDIM);
    
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *TFdat = tforces.dataPtr();
    const Real *GPdat = gp.dataPtr();
    const Real *RHdat = rho.dataPtr();
     
    FORT_SUM_TF_GP(TFdat, ARLIM(tlo), ARLIM(thi),
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
Godunov::Sum_tf_gp_visc (FArrayBox& tforces,
                         FArrayBox& visc, 
                         FArrayBox& gp,
                         FArrayBox& rho)
{
    BL_ASSERT(rho.nComp()     == 1);
    BL_ASSERT(tforces.nComp() >= BL_SPACEDIM);
    BL_ASSERT(visc.nComp()    >= BL_SPACEDIM);
    BL_ASSERT(gp.nComp()      == BL_SPACEDIM);
    
    const int *vlo    = visc.loVect();  
    const int *vhi    = visc.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *glo    = gp.loVect();
    const int *ghi    = gp.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *TFdat = tforces.dataPtr();
    const Real *VIdat = visc.dataPtr();
    const Real *GPdat = gp.dataPtr();
    const Real *RHdat = rho.dataPtr();
     
    FORT_SUM_TF_GP_VISC(TFdat, ARLIM(tlo), ARLIM(thi),
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
Godunov::Sum_tf_divu (FArrayBox& S,
                      FArrayBox& tforces,
                      int        s_ind,
                      int        num_comp,
                      FArrayBox& divu,
                      FArrayBox& rho,
                      int        iconserv)
{
    BL_ASSERT(S.nComp()       >= s_ind+num_comp);
    BL_ASSERT(tforces.nComp() >= s_ind+num_comp);
    BL_ASSERT(divu.nComp()    == 1             );
    BL_ASSERT(rho.nComp()     == 1             );
    
    const int *slo    = S.loVect();
    const int *shi    = S.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *dlo    = divu.loVect();
    const int *dhi    = divu.hiVect();
    const int *rlo    = rho.loVect();
    const int *rhi    = rho.hiVect();
    const Real *Sdat  = S.dataPtr(s_ind);
    const Real *TFdat = tforces.dataPtr(s_ind);
    const Real *DUdat = divu.dataPtr();
    const Real *RHdat = rho.dataPtr();
     
    FORT_SUM_TF_DIVU(Sdat,  ARLIM(slo), ARLIM(shi),
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
Godunov::Sum_tf_divu_visc (FArrayBox& S,
                           FArrayBox& tforces,
                           int        s_ind,
                           int        num_comp,
                           FArrayBox& visc,
                           int        v_ind,
                           FArrayBox& divu,
                           FArrayBox& rho,
                           int        iconserv)
{
    BL_ASSERT(S.nComp()       >= s_ind+num_comp);
    BL_ASSERT(tforces.nComp() >= s_ind+num_comp);
    BL_ASSERT(divu.nComp()    == 1             );
    BL_ASSERT(visc.nComp()    >= v_ind+num_comp);
    BL_ASSERT(rho.nComp()     == 1             );
    
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
    const Real *TFdat = tforces.dataPtr(s_ind);
    const Real *DUdat = divu.dataPtr();
    const Real *VIdat = visc.dataPtr(v_ind);
    const Real *RHdat = rho.dataPtr();
     
    FORT_SUM_TF_DIVU_VISC(Sdat,  ARLIM(slo), ARLIM(shi),
                          TFdat, ARLIM(tlo), ARLIM(thi),
                          DUdat, ARLIM(dlo), ARLIM(dhi),
                          VIdat, ARLIM(vlo), ARLIM(vhi),
                          RHdat, ARLIM(rlo), ARLIM(rhi),
                          tlo, thi, &num_comp, &iconserv);
}


bool
Godunov::are_any(const Array<AdvectionForm>& advectionType,
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
Godunov::how_many(const Array<AdvectionForm>& advectionType,
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
