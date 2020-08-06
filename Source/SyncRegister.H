
#ifndef _SYNCREGISTER_H_
#define _SYNCREGISTER_H_

#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BndryRegister.H>
#include <AMReX_BCRec.H>

class SyncRegister
    :
    public amrex::BndryRegister
{
public:

    SyncRegister (const amrex::BoxArray& fine_boxes, 
                  const amrex::DistributionMapping& dmap,
                  const amrex::IntVect& ref_ratio);

    ~SyncRegister ();

    void CrseInit (amrex::MultiFab& Sync_resid_crse, const amrex::Geometry& crse_geom, amrex::Real mult);

    void FineAdd (amrex::MultiFab& Sync_resid_fine, const amrex::Geometry& crse_geom, amrex::Real mult);

    void CompAdd (amrex::MultiFab& Sync_resid_fine,
                  const amrex::Geometry& fine_geom,
                  const amrex::Geometry& crse_geom,
		  const amrex::BoxArray& Pgrids, amrex::Real mult);

    void InitRHS (amrex::MultiFab& rhs, const amrex::Geometry& geom, const amrex::BCRec& phys_bc);

private:

    amrex::FabSet  bndry_mask[2*AMREX_SPACEDIM];
    amrex::IntVect ratio;
};

#endif /*_SYNCREGISTER_H_*/
