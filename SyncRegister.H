
#ifndef _SYNCREGISTER_H_
#define _SYNCREGISTER_H_

//
// $Id: SyncRegister.H,v 1.17 2000-10-02 20:50:16 lijewski Exp $
//

#include <PArray.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <MultiFab.H>
#include <Geometry.H>
#include <BndryRegister.H>
#include <BCRec.H>

class SyncRegister
    :
    public BndryRegister
{
public:
    SyncRegister ();

    SyncRegister (const BoxArray& fine_boxes,
                  const IntVect&  ref_ratio,
                  int             fine_lev);

    ~SyncRegister ();

    void define (const BoxArray& fine_boxes,
                 const IntVect&  ref_ratio,
                 int             fine_lev);
    //
    // Access functions.
    //
    IntVect refRatio () const { return ratio; }

    int fineLevel () const { return fine_level; }

    int crseLevel () const { return fine_level-1; }

    const BoxArray& coarsenedBoxes () const { return grids; }

    void CrseInit (MultiFab*    Sync_resid_crse,
                   const Geometry&    crse_geom,
                   Real               mult);

    void FineAdd  (MultiFab*    Sync_resid_crse,
                   const Geometry&    fine_geom,
                   const Geometry&    crse_geom,
                   const BCRec*       phys_bc,
                   Real               mult);

    void CompAdd  (MultiFab*    Sync_resid_crse,
                   const Geometry&    fine_geom,
                   const Geometry&    crse_geom,
                   const BCRec*       phys_bc,
                   const BoxArray&    Pgrids,
                   Real               mult);

    void InitRHS (MultiFab&       rhs,
                  const Geometry& geom,
                  const BCRec*    phys_bc);
    //
    // Compute sum of values (overlap values only counted once).
    // WARNING: this is a diagnostic and is not efficient.
    //
    Real sum ();

protected:
    //
    // Helper member functions.
    //
    void incrementPeriodic (const Geometry& geom,
                            const Box&      domain,
                            const MultiFab& mf);

    void copyPeriodic (const Geometry& geom,
                       const Box&      domain,
                       MultiFab&       rhs) const;

    void multByBndryMask (MultiFab& rhs) const;
    //
    // The data.
    //
    FabSet  bndry_mask[2*BL_SPACEDIM];
    IntVect ratio;
    int     fine_level;
};

#endif /*_SYNCREGISTER_H_*/
