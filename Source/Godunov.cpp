
//
// Godunov is the object which calculates advective terms for iamr.
//

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
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

#if (BL_SPACEDIM==2)
    BL_ASSERT(slope_order==1 || slope_order==2 || slope_order==4);
#else
    BL_ASSERT(slope_order==1 || slope_order==4);
#endif

    FORT_SET_PARAMS(slope_order, use_unlimited_slopes);

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
    SetScratch(max_size);
}

//
// Initialize 1d scratch space to a bogus value.
//

void
Godunov::SetBogusScratch ()
{
#ifndef NDEBUG
    const Real bogus_value = 1.e200;

    for (int i = 0 ; i < stxlo.size() ; ++i)
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
    int scr_size = (max_size+2*hyp_grow)*4;
    //
    // Construct arrays.
    //
    D_TERM(stxlo.resize(scr_size);,
           stylo.resize(scr_size);,
           stzlo.resize(scr_size););

    D_TERM(stxhi.resize(scr_size);,
           styhi.resize(scr_size);,
           stzhi.resize(scr_size););

    D_TERM(slxscr.resize(scr_size);,
           slyscr.resize(scr_size);,
           slzscr.resize(scr_size););
}

Godunov::~Godunov ()
{
    ;
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
Godunov::BuildWorkSpace (const Box& grd, const Real* dx, Real dt)
{
    //
    // Ensure 1D scratch space is large enough.
    //
    SetScratch(amrex::grow(grd,hyp_grow).longside());

    work_bx = amrex::grow(grd,1);
    D_TERM(;,
           work.resize(work_bx,5);,
           work.resize(work_bx,19););

    D_TERM(uad.resize(work_bx,1);,
           vad.resize(work_bx,1);,
           wad.resize(work_bx,1););

    SetBogusScratch();

    Box g1box = Box(grd).grow(1);
    smp.resize(g1box,2);
    I.resize(g1box,2*BL_SPACEDIM);

    Box g2box = Box(grd).grow(2);
    dsvl.resize(g2box,1);
    
    if (ppm_type == 2) g1box = Box(grd).grow(2);

    D_TERM(sedgex.resize(amrex::surroundingNodes(g1box,0),1);,
           sedgey.resize(amrex::surroundingNodes(g1box,1),1);,
           sedgez.resize(amrex::surroundingNodes(g1box,2),1););
}

void
Godunov::AllocEdgeBoxes (const Box& grd,
                         D_DECL(FArrayBox& xflux,FArrayBox& yflux,FArrayBox& zflux))
{
    //
    // Create storage for fluxes to be used by calling routine
    //
    D_TERM(xflux_bx = amrex::surroundingNodes(grd,0);,
           yflux_bx = amrex::surroundingNodes(grd,1);,
           zflux_bx = amrex::surroundingNodes(grd,2););

    D_TERM(xflux.resize(xflux_bx,1);,
           yflux.resize(yflux_bx,1);,
           zflux.resize(zflux_bx,1););
}

void
Godunov::ComputeTransverVelocities (const Box& grd, const Real* dx, Real dt,
                                    D_DECL(const int* ubc,const int* vbc,const int* wbc),
                                    D_DECL(const FArrayBox& U,const FArrayBox& V,const FArrayBox& W),
                                    D_DECL(int Ucomp,int Vcomp,int Wcomp),
                                    const FArrayBox& tforces, int Tcomp)
{
    D_TERM(BL_ASSERT(U.nComp() > Ucomp);,
           BL_ASSERT(V.nComp() > Vcomp);,
           BL_ASSERT(W.nComp() > Wcomp););

    BL_ASSERT(work_bx.contains(amrex::grow(grd,1)));
    D_TERM(;,
           BL_ASSERT(work.nComp() >= 5);,
           BL_ASSERT(work.nComp() >= 19););    

    D_TERM(BL_ASSERT(uad.box().contains(work_bx));,
           BL_ASSERT(vad.box().contains(work_bx));,
           BL_ASSERT(wad.box().contains(work_bx)););

    BL_ASSERT(smp.box().contains(amrex::grow(grd,1)));
    BL_ASSERT(I.box().contains(amrex::grow(grd,1)));
    BL_ASSERT(I.nComp() >= 2*BL_SPACEDIM);

    BL_ASSERT(dsvl.box().contains(amrex::grow(grd,2)));

    int nGrow = ppm_type==2 ? 2 : 1;
    const Box& gbox = amrex::grow(grd,nGrow);
    D_TERM(BL_ASSERT(sedgex.box().contains(amrex::surroundingNodes(gbox,0)));,
           BL_ASSERT(sedgey.box().contains(amrex::surroundingNodes(gbox,1)));,
           BL_ASSERT(sedgez.box().contains(amrex::surroundingNodes(gbox,2))););

    //
    // Create the bounds and pointers.
    //
    const int *lo       = grd.loVect();
    const int *hi       = grd.hiVect();
    const int *u_lo     = U.loVect();
    const int *u_hi     = U.hiVect();
    const int *w_lo     = work.loVect();
    const int *w_hi     = work.hiVect();
    const Real *u_dat   = U.dataPtr(Ucomp);
    const Real *uad_dat = uad.dataPtr();
    const Real *v_dat   = V.dataPtr(Vcomp);
    const Real *vad_dat = vad.dataPtr();
    Real* sm = smp.dataPtr(0);
    Real* sp = smp.dataPtr(1);
    Real* Imx = I.dataPtr(0);
    Real* Ipx = I.dataPtr(1);
    Real* Imy = I.dataPtr(2);
    Real* Ipy = I.dataPtr(3);
#if (BL_SPACEDIM == 3)
    Real* Imz = I.dataPtr(4);
    Real* Ipz = I.dataPtr(5);
#endif

#if (BL_SPACEDIM == 3)
    const Real *w_dat   = W.dataPtr(Wcomp);
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
        tforcedat = tforces.dataPtr(Tcomp);
    }
    //
    // Compute the transverse velocities.
    //
    FORT_TRANSVEL(u_dat, uad_dat, xhi_dat, slx_dat, ubc, slxscr.dataPtr(), Imx, Ipx,
                  sedgex.dataPtr(), ARLIM(sedgex.loVect()),  ARLIM(sedgex.hiVect()),
                  v_dat, vad_dat, yhi_dat, sly_dat, vbc, slyscr.dataPtr(), Imy, Ipy, 
                  sedgey.dataPtr(), ARLIM(sedgey.loVect()),  ARLIM(sedgey.hiVect()),
#if (BL_SPACEDIM == 3)
                  w_dat, wad_dat, zhi_dat, slz_dat, wbc, slzscr.dataPtr(), Imz, Ipz,
                  sedgez.dataPtr(), ARLIM(sedgez.loVect()),  ARLIM(sedgez.hiVect()),
#endif    
                  ARLIM(u_lo), ARLIM(u_hi),
                  ARLIM(w_lo), ARLIM(w_hi), 
                  ARLIM(I.loVect()), ARLIM(I.hiVect()),
                  dsvl.dataPtr(), ARLIM(dsvl.loVect()), ARLIM(dsvl.hiVect()),
                  sm, sp, ARLIM(smp.loVect()), ARLIM(smp.hiVect()),
                  lo, hi, &dt, dx, &use_forces_in_trans, tforcedat, &ppm_type);
}

void
Godunov::Setup (const Box& grd, const Real* dx, Real dt, int velpred, 
                FArrayBox& xflux, const int* ubc,
                FArrayBox& yflux, const int* vbc,
#if (BL_SPACEDIM == 3 )
                FArrayBox& zflux, const int* wbc,
#endif
                const FArrayBox& U, const FArrayBox& /*Rho*/, 
                const FArrayBox& tforces)
{
    Setup(grd,dx,dt,velpred,
          D_DECL(xflux,yflux,zflux), D_DECL(ubc,vbc,wbc),
          D_DECL(U,U,U), D_DECL(0,1,2), tforces,0);
}

void
Godunov::Setup (const Box& grd, const Real* dx, Real dt, int velpred,
                D_DECL(FArrayBox& xflux,FArrayBox& yflux,FArrayBox& zflux),
                D_DECL(const int* ubc,const int* vbc,const int* wbc),
                D_DECL(const FArrayBox& U,const FArrayBox& V,const FArrayBox& W),
                D_DECL(int Ucomp,int Vcomp,int Wcomp),
                const FArrayBox& tforces, int Tcomp)
{
    BuildWorkSpace(grd,dx,dt);

    if (!velpred)
        AllocEdgeBoxes(grd,D_DECL(xflux,yflux,zflux));

    ComputeTransverVelocities(grd,dx,dt,D_DECL(ubc,vbc,wbc),D_DECL(U,V,W),
                              D_DECL(Ucomp,Vcomp,Wcomp),tforces,Tcomp);
}

//
// Advection functions follow.
//

//
// Compute the edge states using the advective transverse velocities
// The amount of workspace needed in FArrayBox work is currently 2*SDIM+1.
//
void
Godunov::edge_states (const Box& grd, const Real* dx, Real dt, int velpred,
                      const FArrayBox& uedge, FArrayBox& stx,
                      const FArrayBox& vedge, FArrayBox& sty,
#if (BL_SPACEDIM == 3 )
                      const FArrayBox& wedge, FArrayBox& stz,
#endif
                      const FArrayBox& U, const FArrayBox& S, 
                      const FArrayBox& tforces, const FArrayBox& divu,
                      int fab_ind, int state_ind, const int *bc,
                      int use_conserv, AdvectionScheme which_scheme)
{
    edge_states(grd,dx,dt,velpred,
                D_DECL(uedge,vedge,wedge), D_DECL(0,0,0),
                D_DECL(stx,sty,stz), D_DECL(0,0,0),
                D_DECL(U,U,U), D_DECL(0,1,2),
                S, fab_ind, tforces, fab_ind, divu, 0,
                state_ind, bc, use_conserv, which_scheme);
}

void
Godunov::edge_states( const Box &grd, const Real *dx, Real dt, int velpred,
                      D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                      D_DECL(int mCompX,       int mCompY,       int mCompZ),
                      D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                      D_DECL(int eCompX,       int eCompY,       int eCompZ),
                      D_DECL(const FArrayBox& U, const FArrayBox &V, const FArrayBox &W),
                      D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                      const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, const FArrayBox& divu, int Dcomp,
                      int state_ind, const int *bc, int use_conserv, AdvectionScheme advection_scheme)
{ 

#if (BL_SPACEDIM == 3)
  BL_ASSERT(advection_scheme != BDS);
#endif

  if (advection_scheme == PRE_MAC) 
  {
      edge_states_orig(grd,dx,dt,velpred, D_DECL(uedge,vedge,wedge), D_DECL(mCompX,mCompY,mCompZ), D_DECL(stx,sty,stz),
                       D_DECL(eCompX,eCompY,eCompZ), D_DECL(U,V,W), D_DECL(Ucomp,Vcomp,Wcomp),
                       S,Scomp,tforces,Tcomp,state_ind,bc);
  }
  else if (advection_scheme == FPU) 
  {
      edge_states_fpu(grd,dx,dt,D_DECL(uedge,vedge,wedge), D_DECL(mCompX,mCompY,mCompZ), D_DECL(stx,sty,stz),
                      D_DECL(eCompX,eCompY,eCompZ), S,Scomp,tforces,Tcomp,divu,Dcomp,state_ind,bc,use_conserv);
  }
#if (BL_SPACEDIM == 2)
  else if (advection_scheme == BDS) 
  {
      edge_states_bds(grd,dx,dt,D_DECL(uedge,vedge,wedge), D_DECL(mCompX,mCompY,mCompZ), D_DECL(stx,sty,stz),
                      D_DECL(eCompX,eCompY,eCompZ), S,Scomp,tforces,Tcomp,divu,Dcomp,state_ind,bc,use_conserv);
  }
#endif
} 

void
Godunov::edge_states_orig( const Box &grd, const Real *dx, Real dt, int velpred,
                           D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                           D_DECL(int mCompX,       int mCompY,       int mCompZ),
                           D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                           D_DECL(int eCompX,       int eCompY,       int eCompZ),
                           D_DECL(const FArrayBox& U, const FArrayBox &V, const FArrayBox &W),
                           D_DECL(int Ucomp,        int Vcomp,        int Wcomp),
                           const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp,
                           int state_ind, const int *bc)
{
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(S.nComp()       > Scomp      );
    BL_ASSERT(tforces.nComp() > Tcomp      );

    BL_ASSERT(uedge.nComp()   > mCompX     );
    BL_ASSERT(stx.nComp()     > eCompX     );
    BL_ASSERT(U.nComp()       > Ucomp);

    BL_ASSERT(vedge.nComp()   > mCompY     );
    BL_ASSERT(sty.nComp()     > eCompY     );
    BL_ASSERT(V.nComp()       > Vcomp);

#if (BL_SPACEDIM == 3)
    BL_ASSERT(wedge.nComp()   > mCompZ     );
    BL_ASSERT(stz.nComp()     > eCompZ     );
    BL_ASSERT(W.nComp()       > Wcomp);
#endif    
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *s_lo       = S.loVect();
    const int *s_hi       = S.hiVect();
    const int *ww_lo      = work.loVect();
    const int *ww_hi      = work.hiVect();
    const Real *s_dat     = S.dataPtr(Scomp);
    const Real *u_dat     = U.dataPtr(Ucomp);
    const Real *v_dat     = V.dataPtr(Vcomp);
    const Real *tfr_dat   = tforces.dataPtr(Tcomp);
    const Real *uad_dat   = uad.dataPtr();
    const Real *vad_dat   = vad.dataPtr();
    //
    // Set work space to bogus values.
    //
    SetBogusScratch();

    // Build additional work space
    Box g1box = Box(grd).grow(1);
    smp.resize(g1box,2);
    I.resize(g1box,2*BL_SPACEDIM);

    Real* sm = smp.dataPtr(0);
    Real* sp = smp.dataPtr(1);
    Real* Imx = I.dataPtr(0);
    Real* Ipx = I.dataPtr(1);
    Real* Imy = I.dataPtr(2);
    Real* Ipy = I.dataPtr(3);
#if (BL_SPACEDIM == 3)
    Real* Imz = I.dataPtr(4);
    Real* Ipz = I.dataPtr(5);
#endif

    Box g2box = Box(grd).grow(2);
    dsvl.resize(g2box,1);
    
    if (ppm_type == 2) g1box = Box(grd).grow(2);

    D_TERM(sedgex.resize(g1box.surroundingNodes(0),1);,
           sedgey.resize(g1box.surroundingNodes(1),1);,
           sedgez.resize(g1box.surroundingNodes(2),1););

#if (BL_SPACEDIM == 3)
    const Real *w_dat     = W.dataPtr(Wcomp);
    const Real *wad_dat   = wad.dataPtr();
    const Real *stz_dat   = stz.dataPtr(eCompZ);
    const Real *xhi_dat   = work.dataPtr(0);
    const Real *yhi_dat   = work.dataPtr(1);
    const Real *zhi_dat   = work.dataPtr(2);
    const Real *xlo_dat   = work.dataPtr(3);
    const Real *ylo_dat   = work.dataPtr(4);
    const Real *zlo_dat   = work.dataPtr(5);
    const Real *slx_dat   = work.dataPtr(6);
    const Real *sly_dat   = work.dataPtr(7);
    const Real *slz_dat   = work.dataPtr(8);
    const Real *xedge_dat = work.dataPtr(9);
    const Real *yedge_dat = work.dataPtr(10);
    const Real *zedge_dat = work.dataPtr(11);
    const Real *xylo_dat  = work.dataPtr(12);
    const Real *xzlo_dat  = work.dataPtr(13);
    const Real *yxlo_dat  = work.dataPtr(14);
    const Real *yzlo_dat  = work.dataPtr(15);
    const Real *zxlo_dat  = work.dataPtr(16);
    const Real *zylo_dat  = work.dataPtr(17);
    const Real *xyhi_dat  = work.dataPtr(18);
    const Real *xzhi_dat  = work.dataPtr(18);
    const Real *yxhi_dat  = work.dataPtr(18);
    const Real *yzhi_dat  = work.dataPtr(18);
    const Real *zxhi_dat  = work.dataPtr(18);
    const Real *zyhi_dat  = work.dataPtr(18);
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
                slxscr.dataPtr(), stxlo.dataPtr(), stxhi.dataPtr(),
                uedge.dataPtr(mCompX), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
                stx.dataPtr(eCompX),  ARLIM(  stx.loVect()), ARLIM(  stx.hiVect()), Imx, Ipx,
                sedgex.dataPtr(), ARLIM(sedgex.loVect()), ARLIM(sedgex.hiVect()),

                v_dat, ylo_dat, yhi_dat, sly_dat, vad_dat,
                slyscr.dataPtr(), stylo.dataPtr(), styhi.dataPtr(),
                vedge.dataPtr(mCompY), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
                sty.dataPtr(eCompY),  ARLIM(  sty.loVect()), ARLIM(  sty.hiVect()), Imy, Ipy, 
                sedgey.dataPtr(), ARLIM(sedgey.loVect()), ARLIM(sedgey.hiVect()),

#if (BL_SPACEDIM == 3)
                w_dat, zlo_dat, zhi_dat, slz_dat, wad_dat,
                slzscr.dataPtr(), stzlo.dataPtr(), stzhi.dataPtr(),
                wedge.dataPtr(mCompZ), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
                stz.dataPtr(eCompZ),  ARLIM(  stz.loVect()), ARLIM(  stz.hiVect()), Imz, Ipz,
                sedgez.dataPtr(), ARLIM(sedgez.loVect()), ARLIM(sedgez.hiVect()),

		xedge_dat, yedge_dat, zedge_dat,
		xylo_dat, xzlo_dat, yxlo_dat, yzlo_dat, zxlo_dat, zylo_dat,
		xyhi_dat, xzhi_dat, yxhi_dat, yzhi_dat, zxhi_dat, zyhi_dat,

		&corner_couple,
#endif
                ARLIM(ww_lo), ARLIM(ww_hi),
                ARLIM(I.loVect()), ARLIM(I.hiVect()),
                dsvl.dataPtr(), ARLIM(dsvl.loVect()), ARLIM(dsvl.hiVect()),
                sm, sp, ARLIM(smp.loVect()), ARLIM(smp.hiVect()),
                bc, lo, hi, &dt, dx, &fort_ind, &velpred, 
                &use_forces_in_trans, &ppm_type);
}

void
Godunov::edge_states_fpu( const Box &grd, const Real *dx, Real dt,
                          D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, 
                          const FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int iconserv)
{
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(S.nComp()       > Scomp      );
    BL_ASSERT(tforces.nComp() > Tcomp      );
    BL_ASSERT(divu.nComp()    > Dcomp      );

    BL_ASSERT(uedge.nComp()   > mCompX     );
    BL_ASSERT(stx.nComp()     > eCompX     );

    BL_ASSERT(vedge.nComp()   > mCompY     );
    BL_ASSERT(sty.nComp()     > eCompY     );

#if (BL_SPACEDIM == 3)
    BL_ASSERT(wedge.nComp()   > mCompZ     );
    BL_ASSERT(stz.nComp()     > eCompZ     );
#endif    
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *s_lo       = S.loVect();
    const int *s_hi       = S.hiVect();
    const int *ww_lo      = work.loVect();
    const int *ww_hi      = work.hiVect();
    const Real *s_dat     = S.dataPtr(Scomp);
    const Real *tfr_dat   = tforces.dataPtr(Tcomp);
    const int *t_lo       = tforces.loVect();
    const int *t_hi       = tforces.hiVect();
    const Real *divu_dat  = divu.dataPtr(Dcomp);
    const int *d_lo       = divu.loVect();
    const int *d_hi       = divu.hiVect();
    //
    // Set work space to bogus values.
    //
    SetBogusScratch();

    // Build additional work space
    Box g1box = Box(grd).grow(1);
    smp.resize(g1box,2);
    I.resize(g1box,2*BL_SPACEDIM);

    Real* sm = smp.dataPtr(0);
    Real* sp = smp.dataPtr(1);
    Real* Imx = I.dataPtr(0);
    Real* Ipx = I.dataPtr(1);
    Real* Imy = I.dataPtr(2);
    Real* Ipy = I.dataPtr(3);
#if (BL_SPACEDIM == 3)
    Real* Imz = I.dataPtr(4);
    Real* Ipz = I.dataPtr(5);
#endif

    Box g2box = Box(grd).grow(2);
    dsvl.resize(g2box,1);
    
    if (ppm_type == 2) g1box = Box(grd).grow(2);

    D_TERM(sedgex.resize(g1box.surroundingNodes(0),1);,
           sedgey.resize(g1box.surroundingNodes(1),1);,
           sedgez.resize(g1box.surroundingNodes(2),1););

#if (BL_SPACEDIM == 3)
    const Real *xhi_dat   = work.dataPtr(0);
    const Real *yhi_dat   = work.dataPtr(1);
    const Real *zhi_dat   = work.dataPtr(2);
    const Real *xlo_dat   = work.dataPtr(3);
    const Real *ylo_dat   = work.dataPtr(4);
    const Real *zlo_dat   = work.dataPtr(5);
    const Real *slx_dat   = work.dataPtr(6);
    const Real *sly_dat   = work.dataPtr(7);
    const Real *slz_dat   = work.dataPtr(8);
    const Real *xedge_dat = work.dataPtr(9);
    const Real *yedge_dat = work.dataPtr(10);
    const Real *zedge_dat = work.dataPtr(11);
    const Real *xylo_dat  = work.dataPtr(12);
    const Real *xzlo_dat  = work.dataPtr(13);
    const Real *yxlo_dat  = work.dataPtr(14);
    const Real *yzlo_dat  = work.dataPtr(15);
    const Real *zxlo_dat  = work.dataPtr(16);
    const Real *zylo_dat  = work.dataPtr(17);
    const Real *xyhi_dat  = work.dataPtr(18);
    const Real *xzhi_dat  = work.dataPtr(18);
    const Real *yxhi_dat  = work.dataPtr(18);
    const Real *yzhi_dat  = work.dataPtr(18);
    const Real *zxhi_dat  = work.dataPtr(18);
    const Real *zyhi_dat  = work.dataPtr(18);
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

    FORT_ESTATE_FPU(s_dat,    ARLIM(s_lo), ARLIM(s_hi),
                    tfr_dat,  ARLIM(t_lo), ARLIM(t_hi),
                    divu_dat, ARLIM(d_lo), ARLIM(d_hi),
                    
                    xlo_dat, xhi_dat, slx_dat,
                    slxscr.dataPtr(), stxlo.dataPtr(), stxhi.dataPtr(),
                    uedge.dataPtr(mCompX), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
                    stx.dataPtr(eCompX),   ARLIM(stx.loVect()),   ARLIM(stx.hiVect()), Imx, Ipx,
                    sedgex.dataPtr(), ARLIM(sedgex.loVect()), ARLIM(sedgex.hiVect()),
                    
                    ylo_dat, yhi_dat, sly_dat,
                    slyscr.dataPtr(), stylo.dataPtr(), styhi.dataPtr(),
                    vedge.dataPtr(mCompY), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
                    sty.dataPtr(eCompY),   ARLIM(sty.loVect()),   ARLIM(sty.hiVect()), Imy, Ipy,
                    sedgey.dataPtr(), ARLIM(sedgey.loVect()), ARLIM(sedgey.hiVect()),

#if (BL_SPACEDIM == 3)
                    zlo_dat, zhi_dat, slz_dat,
                    slzscr.dataPtr(), stzlo.dataPtr(), stzhi.dataPtr(),
                    wedge.dataPtr(mCompZ), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
                    stz.dataPtr(eCompZ),   ARLIM(stz.loVect()),   ARLIM(stz.hiVect()), Imz, Ipz,
                    sedgez.dataPtr(), ARLIM(sedgez.loVect()), ARLIM(sedgez.hiVect()),

		    xedge_dat, yedge_dat, zedge_dat,
		    xylo_dat, xzlo_dat, yxlo_dat, yzlo_dat, zxlo_dat, zylo_dat,
		    xyhi_dat, xzhi_dat, yxhi_dat, yzhi_dat, zxhi_dat, zyhi_dat,

		    &corner_couple,
#endif
                    ARLIM(ww_lo), ARLIM(ww_hi),
                    ARLIM(I.loVect()), ARLIM(I.hiVect()),
                    dsvl.dataPtr(), ARLIM(dsvl.loVect()), ARLIM(dsvl.hiVect()),
                    sm, sp, ARLIM(smp.loVect()), ARLIM(smp.hiVect()),
                    bc, lo, hi, &dt, dx, &fort_ind,
                    &use_forces_in_trans, &iconserv, &ppm_type);

}

void
Godunov::edge_states_bds( const Box &grd, const Real *dx, Real dt,
                          D_DECL(const FArrayBox &uedge, const FArrayBox &vedge, const FArrayBox &wedge),
                          D_DECL(int mCompX,       int mCompY,       int mCompZ),
                          D_DECL(FArrayBox &stx,   FArrayBox &sty,   FArrayBox &stz  ),
                          D_DECL(int eCompX,       int eCompY,       int eCompZ),
                          const FArrayBox &S, int Scomp, const FArrayBox &tforces, int Tcomp, const FArrayBox& divu, int Dcomp,
                          int state_ind, const int *bc, int iconserv)
{
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(S.nComp()       > Scomp      );
    BL_ASSERT(tforces.nComp() > Tcomp      );
    BL_ASSERT(divu.nComp()    > Dcomp      );

    BL_ASSERT(uedge.nComp()   > mCompX     );
    BL_ASSERT(stx.nComp()     > eCompX     );

    BL_ASSERT(vedge.nComp()   > mCompY     );
    BL_ASSERT(sty.nComp()     > eCompY     );

#if (BL_SPACEDIM == 3)
    BL_ASSERT(wedge.nComp()   > mCompZ     );
    BL_ASSERT(stz.nComp()     > eCompZ     );
#endif    
    //
    // Create the bounds and pointers.
    //
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const int *s_lo       = S.loVect();
    const int *s_hi       = S.hiVect();
    const int *ww_lo      = work.loVect();
    const int *ww_hi      = work.hiVect();
    const Real *s_dat     = S.dataPtr(Scomp);
    const Real *tfr_dat   = tforces.dataPtr(Tcomp);
    const Real *divu_dat  = divu.dataPtr(Dcomp);
    //
    // Set work space to bogus values.
    //
    SetBogusScratch();

#if (BL_SPACEDIM == 3)
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
      
#if (BL_SPACEDIM == 2)
    int fort_ind = state_ind+1;  
    FORT_ESTATE_BDS(s_dat, tfr_dat, divu_dat, ARLIM(s_lo), ARLIM(s_hi),
                    
                    xlo_dat, xhi_dat, slx_dat,
                    slxscr.dataPtr(), stxlo.dataPtr(), stxhi.dataPtr(),
                    uedge.dataPtr(mCompX), ARLIM(uedge.loVect()), ARLIM(uedge.hiVect()),
                    stx.dataPtr(eCompX),   ARLIM(stx.loVect()),   ARLIM(stx.hiVect()),
                    
                    ylo_dat, yhi_dat, sly_dat,
                    slyscr.dataPtr(), stylo.dataPtr(), styhi.dataPtr(),
                    vedge.dataPtr(mCompY), ARLIM(vedge.loVect()), ARLIM(vedge.hiVect()),
                    sty.dataPtr(eCompY),   ARLIM(sty.loVect()),   ARLIM(sty.hiVect()),
#if (BL_SPACEDIM == 3)
                    zlo_dat, zhi_dat, slz_dat,
                    slzscr.dataPtr(), stzlo.dataPtr(), stzhi.dataPtr(),
                    wedge.dataPtr(mCompZ), ARLIM(wedge.loVect()), ARLIM(wedge.hiVect()),
                    stz.dataPtr(eCompZ),   ARLIM(stz.loVect()),   ARLIM(stz.hiVect()),
#endif
                    ARLIM(ww_lo), ARLIM(ww_hi),
                    bc, lo, hi, &dt, dx, &fort_ind,
                    &use_forces_in_trans, &iconserv);
#endif
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
    int mCompX = 0;
    int mCompY = 0;
    //
    // 2D calls.
    //
#if (BL_SPACEDIM == 2)                  
    edge_states_orig(grd,dx,dt,velpred,umac,vmac,mCompX,mCompY,umac,vmac,mCompX,mCompY,
                     U,U,XVEL,YVEL,U,XVEL,tforces,XVEL,XVEL,ubc);

    edge_states_orig(grd,dx,dt,velpred,umac,vmac,mCompX,mCompY,umac,vmac,mCompX,mCompY,
                     U,U,XVEL,YVEL,U,YVEL,tforces,YVEL,YVEL,vbc);
#endif
    //
    // 3D calls.
    //
#if (BL_SPACEDIM == 3)                  
    int mCompZ = 0;

    edge_states_orig(grd,dx,dt,velpred,umac,vmac,wmac,mCompX,mCompY,mCompZ,umac,vmac,wmac,mCompX,mCompY,mCompZ,
                     U,U,U,XVEL,YVEL,ZVEL,U,XVEL,tforces,XVEL,XVEL,ubc);

    edge_states_orig(grd,dx,dt,velpred,umac,vmac,wmac,mCompX,mCompY,mCompZ,umac,vmac,wmac,mCompX,mCompY,mCompZ,
                     U,U,U,XVEL,YVEL,ZVEL,U,YVEL,tforces,YVEL,YVEL,vbc);

    edge_states_orig(grd,dx,dt,velpred,umac,vmac,wmac,mCompX,mCompY,mCompZ,umac,vmac,wmac,mCompX,mCompY,mCompZ,
                     U,U,U,XVEL,YVEL,ZVEL,U,ZVEL,tforces,ZVEL,ZVEL,wbc);

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
    int velpred = 0;
    //
    // Compute edge states for an advected quantity.
    //
    edge_states(grd, dx, dt, velpred,
                D_DECL(uedge,vedge,wedge), D_DECL(0,0,0),
                D_DECL(xflux,yflux,zflux), D_DECL(0,0,0),
                D_DECL(U,U,U), D_DECL(0,1,2),
                S, fab_ind, tforces, fab_ind, divu, 0, state_ind, bc, 
                iconserv, scheme);

    ComputeAofs (grd,
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

    FORT_ADV_FORCING( aofs.dataPtr(acomp),ARLIM(aofs.loVect()), ARLIM(aofs.hiVect()),

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
Godunov::SyncAdvect (const Box&  grd,
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
    int velpred = 0;
    //
    // Error block.
    //
    BL_ASSERT(S.box().contains(work_bx));

    BL_ASSERT(S.nComp()       >= BL_SPACEDIM);
    BL_ASSERT(S.nComp()       >= fab_ind    );
    BL_ASSERT(tforces.nComp() >= fab_ind    );
    BL_ASSERT(sync.nComp()    >= sync_ind   );

    BL_ASSERT(ucorr.box()     == xflux.box());
    BL_ASSERT(ucorr.nComp()   >= 1          );

    BL_ASSERT(vcorr.box()     == yflux.box());
    BL_ASSERT(vcorr.nComp()   >= 1          );
#if (BL_SPACEDIM == 3)
    BL_ASSERT(wcorr.box()     == zflux.box());
    BL_ASSERT(wcorr.nComp()   >= 1          );
#endif    
    //
    // Compute the edge states.
    //
    edge_states(grd, dx, dt, velpred,
                D_DECL(uedge,vedge,wedge), D_DECL(0,0,0),
                D_DECL(xflux,yflux,zflux), D_DECL(0,0,0),
                D_DECL(U,U,U), D_DECL(0,1,2),
                S, fab_ind, tforces, fab_ind, divu, 0, sync_ind, bc,
                iconserv, scheme);

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
    FORT_SYNC_ADV_FORCING(sync.dataPtr(sync_ind), ARLIM(sync.loVect()), ARLIM(sync.hiVect()),
                           
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
    const int *slo        = Sold.loVect();
    const int *shi        = Sold.hiVect();
    const int *lo         = grd.loVect();
    const int *hi         = grd.hiVect();
    const Real *Sold_dat  = Sold.dataPtr(ind_old_s);
    const Real *Snew_dat  = Snew.dataPtr(ind_new_s);
    const Real *Rho_dat   = Sold.dataPtr(ind_old_rho);
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

    FORT_CONSSCALMINMAX (Sold_dat, Snew_dat, Rho_dat, Rhon_dat,
                         ARLIM(slo), ARLIM(shi),
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

    FORT_CONVSCALMINMAX (Sold_dat, 
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
    FORT_MAXCHNG_VELMAG(Uodat, ARLIM(uo_lo), ARLIM(uo_hi),
               	   Undat, ARLIM(un_lo), ARLIM(un_hi),
               	   lo, hi, &max_change);
    return max_change;
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

    const int *slo    = Sold.loVect();
    const int *shi    = Sold.hiVect();
    const int *tlo    = tforces.loVect();
    const int *thi    = tforces.hiVect();
    const int *lo     = grd.loVect();
    const int *hi     = grd.hiVect();
    const Real *SOdat = Sold.dataPtr(start_ind);
    Real *SNdat = Snew.dataPtr(start_ind);
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
Godunov::Correct_tf (const FArrayBox& Sstar,
                     FArrayBox& Snp1,
                     int        start_ind,
                     int        num_comp, 
                     const FArrayBox& tfstar,
                     const FArrayBox& tfn,
                     int        tf_ind,
                     const Box& grd,
                     Real       dt) const
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
    Real *SPdat = Snp1.dataPtr(start_ind);
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

    const int *slo    = Sold.loVect();
    const int *shi    = Sold.hiVect();
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
    Real *UNdat = Unew.dataPtr();
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
     
    FORT_SUM_TF_DIVU_VISC(Sdat,  ARLIM(slo), ARLIM(shi),
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
