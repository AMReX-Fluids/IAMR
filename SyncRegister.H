//BL_COPYRIGHT_NOTICE

#ifndef _SYNCREGISTER_H_
#define _SYNCREGISTER_H_

//
// $Id: SyncRegister.H,v 1.13 1998-12-11 21:10:47 lijewski Exp $
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

    void CrseDVInit (const MultiFab& V,
                     const Geometry& geom, 
                     int             is_rz,
                     int**           crse_bc,
                     Real            mult = 1.0);

    void FineDVAdd (const MultiFab& V, 
                    const Real*     dx_fine,
                    const Geometry& crse_geom,
                    int             is_rz,
                    int**           fine_bc,
                    Real            mult);

    void CrseDsdtAdd (const MultiFab& dsdt,
                      const Geometry& geom,
                      int             is_rz,
                      int**           crse_bc, 
                      int             lowfix,
                      int             hifix,
                      Real            mult = 1.0);

    void FineDsdtAdd (const MultiFab& dsdt,
                      const Geometry& geom,
                      const Geometry& crse_geom,
                      int             is_rz,
                      int**           fine_bc, 
                      int             lowfix,
                      int             hifix,
                      Real             mult);

    void CompDVAdd (const MultiFab& V, 
                    const BoxArray& fine_boxes,
                    const Real*     dx_fine, 
                    const Geometry& fine_geom, 
                    const Geometry& crse_geom, 
                    int             is_rz,
                    int**           fine_bc,
                    Real            mult);

    void CrseLPhiAdd (const MultiFab& Phi,
                      const MultiFab& sigma,
                      const Geometry& geom,
                      int             is_rz,
                      Real            mult = 1.0);

    void FineLPhiAdd (const MultiFab& Phi,
                      const MultiFab& sigma,
                      const Real*     dx_fine,
                      const Geometry& geom, 
                      int             is_rz,
                      Real            mult);

    void CompLPhiAdd (const MultiFab& Phi,
                      const MultiFab& sigma,
                      const BoxArray& fine_boxes, 
                      const Real*     dx_fine, 
                      const Geometry& fine_geom, 
                      const Geometry& crse_geom, 
                      int             is_rz,
                      Real            mult);

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
