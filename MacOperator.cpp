//
// $Id: MacOperator.cpp,v 1.38 2006-09-08 21:21:02 almgren Exp $
//
#include <winstd.H>

#include <MacBndry.H>
#include <MacOperator.H>
#include <MacOpMacDrivers.H>
#include <MACOPERATOR_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

#ifdef MG_USE_HYPRE
#include <HypreABec.H>
#endif

#ifdef MG_USE_FBOXLIB
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#endif

#ifndef _PorousMedia_H_
enum StateType {State_Type=0, Press_Type};
#endif

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

MacOperator::MacOperator (const BndryData& mgb,
                          const Real*      h)
    :
    ABecLaplacian(mgb,h)
{}

MacOperator::~MacOperator () {}

//
// Define the meaning of gradient for the multigrid object.
//

void
MacOperator::setCoefficients (MultiFab*   area,
                              MultiFab&   mac_coef,
                              const Real* dx)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];
    BL_ASSERT(mac_coef.boxArray() == ba);
    //
    // First set scalar coeficients.
    //
    setScalars(0.0,1.0);
    //
    // Don't need to set a because alpha is set to zero.
    //
    const int n_grow = 0;

    D_TERM(MultiFab bxcoef(area[0].boxArray(),area[0].nComp(),n_grow);,
           MultiFab bycoef(area[1].boxArray(),area[1].nComp(),n_grow);,
           MultiFab bzcoef(area[2].boxArray(),area[2].nComp(),n_grow););
    D_TERM(bxcoef.setVal(0);,
           bycoef.setVal(0);,
           bzcoef.setVal(0););

    for (MFIter mfi(mac_coef); mfi.isValid(); ++mfi)
    {
        BL_ASSERT(ba[mfi.index()] == mfi.validbox());

        const Box& grd       = ba[mfi.index()];
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& bx        = bxcoef[mfi];
        FArrayBox& by        = bycoef[mfi];
        const FArrayBox& ax  = area[0][mfi];
        const FArrayBox& ay  = area[1][mfi];
        const FArrayBox& mcc = mac_coef[mfi];

        DEF_LIMITS(bx,bx_dat,bxlo,bxhi);
        DEF_LIMITS(by,by_dat,bylo,byhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);

        const int* mlo      = mcc.loVect();
        const int* mhi      = mcc.hiVect();
        const Real* mcc_dat = mcc.dataPtr();

#if (BL_SPACEDIM == 2)
        FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     mcc_dat,ARLIM(mlo),ARLIM(mhi),lo,hi,dx);
#endif
#if (BL_SPACEDIM == 3)
        FArrayBox& bz       = bzcoef[rhomfi];
        const FArrayBox& az = area[2][rhomfi];

        DEF_CLIMITS(az,az_dat,azlo,azhi);
        DEF_LIMITS(bz,bz_dat,bzlo,bzhi);

        FORT_MACCOEF(bx_dat,ARLIM(bxlo),ARLIM(bxhi),
                     by_dat,ARLIM(bylo),ARLIM(byhi),
                     bz_dat,ARLIM(bzlo),ARLIM(bzhi),
                     ax_dat,ARLIM(axlo),ARLIM(axhi),
                     ay_dat,ARLIM(aylo),ARLIM(ayhi),
                     az_dat,ARLIM(azlo),ARLIM(azhi),
                     mcc_dat,ARLIM(mlo),ARLIM(mhi),lo,hi,dx);
#endif
    }
  
    D_TERM(bCoefficients(bxcoef,0);,
           bCoefficients(bycoef,1);,
           bCoefficients(bzcoef,2););
}

//
// This function creates the initial rhs for use in the mac multgrid solve.
//

void
MacOperator::defRHS (MultiFab* area,
                     MultiFab& volume,
                     MultiFab& Rhs,
                     MultiFab* vel,
                     Real      scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];
    BL_ASSERT(Rhs.boxArray() == ba);

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(ba[Rhsmfi.index()] == Rhsmfi.validbox());

        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        const FArrayBox& ax  = area[0][Rhsmfi];
        const FArrayBox& ay  = area[1][Rhsmfi];
        const FArrayBox& vol = volume[Rhsmfi];
        const FArrayBox& ux  = vel[0][Rhsmfi];
        const FArrayBox& uy  = vel[1][Rhsmfi];
        FArrayBox& rhs       = Rhs[Rhsmfi];

        DEF_CLIMITS(ux,ux_dat,uxlo,uxhi);
        DEF_CLIMITS(uy,uy_dat,uylo,uyhi);
        DEF_CLIMITS(ax,ax_dat,axlo,axhi);
        DEF_CLIMITS(ay,ay_dat,aylo,ayhi);
        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);

#if (BL_SPACEDIM == 2)
        FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi), 
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,&scale);
#endif
#if (BL_SPACEDIM == 3)
        const FArrayBox& az = area[2][Rhsmfi];
        DEF_CLIMITS(az,az_dat,azlo,azhi);

        const FArrayBox& uz = vel[2][Rhsmfi];
        DEF_CLIMITS(uz,uz_dat,uzlo,uzhi);

        FORT_MACRHS(ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                    uy_dat,ARLIM(uylo),ARLIM(uyhi),
                    uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                    ax_dat,ARLIM(axlo),ARLIM(axhi),
                    ay_dat,ARLIM(aylo),ARLIM(ayhi),
                    az_dat,ARLIM(azlo),ARLIM(azhi),
                    vol_dat,ARLIM(vlo),ARLIM(vhi),
                    rhs_dat,ARLIM(rlo),ARLIM(rhi),
                    lo,hi,&scale);
#endif
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Apply the mac pressure gradient to a velocity field.
// init, means that velocities are initialized here.
//

void
mac_vel_update (int              init,
                D_DECL(FArrayBox& ux,
                       FArrayBox& uy,
                       FArrayBox& uz),
                const FArrayBox& phi,
                const FArrayBox* mac_coef_ptr,
                const Box&       grd,
                int              level,
                int              n,
                const Real*      dx,
                Real             scale)
{
    const int* lo        = grd.loVect();
    const int* hi        = grd.hiVect();

    const FArrayBox& mac_coef = *mac_coef_ptr;
    
    DEF_LIMITS(ux,ux_dat,uxlo,uxhi);
    DEF_LIMITS(uy,uy_dat,uylo,uyhi);
    DEF_CLIMITS(phi,phi_dat,p_lo,p_hi);

    const int* mlo      = mac_coef.loVect();
    const int* mhi      = mac_coef.hiVect();
    const Real* mcc_dat = mac_coef.dataPtr();
    
#if (BL_SPACEDIM == 2)
    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   mcc_dat,ARLIM(mlo),ARLIM(mhi),
                   lo,hi,dx,&scale);
#endif
#if (BL_SPACEDIM == 3)
    DEF_LIMITS(uz,uz_dat,uzlo,uzhi);
    
    FORT_MACUPDATE(&init,
                   ux_dat,ARLIM(uxlo),ARLIM(uxhi),
                   uy_dat,ARLIM(uylo),ARLIM(uyhi),
                   uz_dat,ARLIM(uzlo),ARLIM(uzhi),
                   phi_dat,ARLIM(p_lo),ARLIM(p_hi),
                   mcc_dat,ARLIM(mlo),ARLIM(mhi),
                   lo,hi,dx,&scale);
#endif
}

//
// Apply the mac pressure gradient to the divergent mac velocities.
// The resultant velocity field is nondivergent.
//

void
MacOperator::velUpdate (MultiFab*       Vel,
                        MultiFab&       Phi,
                        const MultiFab& mac_coef,
                        const Real*     dx,
                        Real            scale)
{
    //
    // Should check that all BoxArrays are consistant.
    //
    const BoxArray& ba = gbox[0];
    BL_ASSERT(mac_coef.boxArray() == ba);
    //
    // Set bndry data in ghost zones.
    //
    int apply_lev = 0;
    applyBC(Phi,0,1,apply_lev);

    for (MFIter Phimfi(Phi); Phimfi.isValid(); ++Phimfi)
    {
        BL_ASSERT(ba[Phimfi.index()] == Phimfi.validbox());

        const Box& grd = Phimfi.validbox();

        mac_vel_update(0, 
                       D_DECL(Vel[0][Phimfi],Vel[1][Phimfi],Vel[2][Phimfi]),
                       Phi[Phimfi],
                       &(mac_coef[Phimfi]), 
                       grd, 0, Phimfi.index(), dx, scale );
    }
}

//
// Multiply by volume*rhs_scale since reflux step (which computed rhs)
// divided by volume.
//

void
MacOperator::syncRhs (const MultiFab& Volume,
                      MultiFab&       Rhs,
                      Real            rhs_scale,
                      const Real*     dx)
{
    const BoxArray& ba = gbox[0];

    for (MFIter Rhsmfi(Rhs); Rhsmfi.isValid(); ++Rhsmfi)
    {
        BL_ASSERT(ba[Rhsmfi.index()] == Rhsmfi.validbox());

        const Box& grd       = Rhsmfi.validbox();
        const int* lo        = grd.loVect();
        const int* hi        = grd.hiVect();
        FArrayBox& rhs       = Rhs[Rhsmfi];
        const FArrayBox& vol = Volume[Rhsmfi];

        DEF_CLIMITS(vol,vol_dat,vlo,vhi);
        DEF_LIMITS(rhs,rhs_dat,rlo,rhi);
        FORT_MACSYNCRHS(rhs_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                        vol_dat,ARLIM(vlo),ARLIM(vhi),&rhs_scale);
    }
    Rhs.mult(-1.0,Rhs.nGrow());
}

//
// Driver functions follow.
//

//
// A driver function for computing a level MAC solve.
//

void
mac_level_driver (const MacBndry& mac_bndry,
		  const BCRec&    phys_bc,
                  const BoxArray& grids,
                  int             the_solver,
                  int             level,
                  const Real*     dx,
                  Real            dt,
                  Real            mac_tol,
                  Real            mac_abs_tol,
                  Real            rhs_scale,
                  MultiFab*       area,
                  MultiFab&       volume,
                  MultiFab&       mac_coef,
                  MultiFab&       Rhs,
                  MultiFab*       u_mac,
                  MultiFab*       mac_phi)
{
  BL_PROFILE("mac_level_driver");
  MacOperator mac_op(mac_bndry,dx);
  mac_op.setCoefficients(area,mac_coef,dx);
  mac_op.defRHS(area,volume,Rhs,u_mac,rhs_scale);
  mac_op.maxOrder(2);
  if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
      BoxLib::Error("Can't use CGSolver with maxorder > 2");
    }
  //
  // Construct MultiGrid or CGSolver object and solve system.
  //
  if (the_solver == 1)
    {
      bool use_mg_precond = true;
      CGSolver mac_cg(mac_op,use_mg_precond);
      mac_cg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
  else if (the_solver == 3 ) 
    {
#ifdef MG_USE_FBOXLIB
      std::vector<BoxArray> bav(1);
      bav[0] = mac_phi->boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = mac_bndry.getGeom();

      int mg_bc[2*BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      { 
        if ( geom[0].isPeriodic(i) )
        {
          mg_bc[i*2 + 0] = 0;
          mg_bc[i*2 + 1] = 0;
        }
      else
        {
          mg_bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
          mg_bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
        }
      }
      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);

      const MultiFab* aa_p[1]; 
      aa_p[0] = &(mac_op.aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      {
          bb_p[0][i] = &(mac_op.bCoefficients(i));
      }
      mgt_solver.set_mac_coefficients(aa_p, bb_p, mac_bndry);

      MultiFab* mac_phi_p[1];
      MultiFab* Rhs_p[1];
      mac_phi_p[0] = mac_phi;
      Rhs_p[0] = &Rhs;

      mgt_solver.solve(mac_phi_p, Rhs_p, mac_tol, mac_abs_tol, mac_bndry);
#else
      BoxLib::Error("mac_level_driver::mg_cpp not in this build");
#endif
    }
  else
    {
      MultiGrid mac_mg(mac_op);
      mac_mg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
  //
  // velUpdate will set bndry values for mac_phi.
  //
  mac_op.velUpdate(u_mac,*mac_phi,mac_coef,dx,-dt/2.0);
}

//
// A driver function for computing a sync MAC solve.
//

void
mac_sync_driver (const MacBndry& mac_bndry,
	         const BCRec&    phys_bc,
                 const BoxArray& grids,
                 int             the_solver,
                 int             level, 
                 const Real*     dx,
                 Real            dt,
                 Real            mac_sync_tol,
                 Real            mac_abs_tol,
                 Real            rhs_scale,
                 MultiFab*       area,
                 MultiFab&       volume,
                 MultiFab&       Rhs,
                 MultiFab*       mac_coef,
                 MultiFab*       mac_sync_phi)
{
  BL_PROFILE("mac_sync_driver");
  MacOperator mac_op(mac_bndry,dx);
  mac_op.maxOrder(2);
  mac_op.setCoefficients(area,*mac_coef, dx);
  mac_op.syncRhs(volume,Rhs,rhs_scale,dx);
  if (the_solver == 1 && mac_op.maxOrder() != 2)
    {
      BoxLib::Error("Can't use CGSolver with maxorder > 2");
    }
  //
  // Now construct MultiGrid or CGSolver object to solve system.
  //
  if (the_solver == 1)
    {
      bool use_mg_precond = true;
      CGSolver mac_cg(mac_op,use_mg_precond);
      mac_cg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
  else if ( the_solver == 2 )
    {
#ifdef MG_USE_HYPRE
      HypreABec hp(mac_sync_phi->boxArray(), mac_bndry, dx, 0, false);
      hp.setScalars(mac_op.get_alpha(), mac_op.get_beta());
      hp.aCoefficients(mac_op.aCoefficients());
      for ( int i = 0; i < BL_SPACEDIM; ++i )
        {
	  hp.bCoefficients(mac_op.bCoefficients(i), i);
        }
      hp.setup_solver(mac_sync_tol, mac_abs_tol, 50);
      hp.solve(*mac_sync_phi, Rhs, true);
      hp.clear_solver();
#else
      BoxLib::Error("mac_sync_driver: HypreABec not in this build");
#endif
    }
  else if (the_solver == 3 )
    {
#ifdef MG_USE_FBOXLIB
      std::vector<BoxArray> bav(1);
      bav[0] = mac_sync_phi->boxArray();
      std::vector<DistributionMapping> dmv(1);
      dmv[0] = Rhs.DistributionMap();
      bool nodal = false;
      std::vector<Geometry> geom(1);
      geom[0] = mac_bndry.getGeom();

      int mg_bc[2*BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      { 
        if ( geom[0].isPeriodic(i) )
        {
          mg_bc[i*2 + 0] = 0;
          mg_bc[i*2 + 1] = 0;
        }
      else
        {
          mg_bc[i*2 + 0] = phys_bc.lo(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
          mg_bc[i*2 + 1] = phys_bc.hi(i)==Outflow? MGT_BC_DIR : MGT_BC_NEU;
        }
      }

      MGT_Solver mgt_solver(geom, mg_bc, bav, dmv, nodal);

      const MultiFab* aa_p[1];
      aa_p[0] = &(mac_op.aCoefficients());
      const MultiFab* bb_p[1][BL_SPACEDIM];
      for ( int i = 0; i < BL_SPACEDIM; ++i )
      {
          bb_p[0][i] = &(mac_op.bCoefficients(i));
      }
      mgt_solver.set_mac_coefficients(aa_p, bb_p, mac_bndry);

      MultiFab* mac_phi_p[1];
      MultiFab* Rhs_p[1];
      mac_phi_p[0] = mac_sync_phi;
      Rhs_p[0] = &Rhs;

      mgt_solver.solve(mac_phi_p, Rhs_p, mac_sync_tol, mac_abs_tol, mac_bndry);
#else
      BoxLib::Error("mac_sync_driver::mg_cpp not in this build");
#endif
    }
  else
    {
      MultiGrid mac_mg(mac_op);
      mac_mg.solve(*mac_sync_phi,Rhs,mac_sync_tol,mac_abs_tol);
    }
    
  int mac_op_lev = 0;
  mac_op.applyBC(*mac_sync_phi,0,1,mac_op_lev);
}
