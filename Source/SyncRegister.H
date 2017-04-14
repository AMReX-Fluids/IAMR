
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
    public BndryRegister
{
public:

    SyncRegister (const BoxArray& fine_boxes, const IntVect& ref_ratio);

    ~SyncRegister ();

    void CrseInit (MultiFab& Sync_resid_crse, const Geometry& crse_geom, Real mult);

    void FineAdd (MultiFab& Sync_resid_fine, const Geometry& crse_geom, Real mult);

    void CompAdd (MultiFab& Sync_resid_fine, const Geometry& fine_geom, const Geometry& crse_geom,
		  const BoxArray& Pgrids, Real mult);

    void InitRHS (MultiFab& rhs, const Geometry& geom, const BCRec& phys_bc);

private:

    FabSet  bndry_mask[2*BL_SPACEDIM];
    IntVect ratio;
};

#endif /*_SYNCREGISTER_H_*/
