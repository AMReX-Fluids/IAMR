//BL_COPYRIGHT_NOTICE

//
// $Id: SyncRegister.cpp,v 1.49 1998-12-09 23:10:05 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <cstring>
#include <cstdio>
#else
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#endif

#include <BC_TYPES.H>
#include <SyncRegister.H>

#include <NAVIERSTOKES_F.H>
#ifndef FORT_HGC2N
#include <PROJECTION_F.H>
#endif

#include <SYNCREG_F.H>

SyncRegister::SyncRegister ()
{
    fine_level = -1;
    ratio      = IntVect::TheUnitVector();
    ratio.scale(-1);
}

SyncRegister::SyncRegister (const BoxArray& fine_boxes,
                            const IntVect&  ref_ratio,
                            int             fine_lev)
{
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
    define(fine_boxes,ref_ratio,fine_lev);
}

void
SyncRegister::define (const BoxArray& fine_boxes,
                      const IntVect&  ref_ratio,
                      int             fine_lev)
{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        assert(ratio[dir] == -1);
    assert(fine_boxes.isDisjoint());
    assert(!grids.ready());

    ratio      = ref_ratio;
    fine_level = fine_lev;
 
    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].resize(grids.length());
        bndry[face()].DefineGrids(grids);
        bndry[face()].DefineDistributionMap(grids);

        bndry_mask[face()].resize(grids.length());
        bndry_mask[face()].DefineGrids(grids);
        bndry_mask[face()].DefineDistributionMap(grids);
    }
    //
    // Construct disjoint "face" fabs that are node centered in all index
    // directions.  The "grow" function at the bottom of this loop shrinks
    // the domain in the index direction just allocated so that the fabs
    // in the other directions will not overlap.
    //
    const int myproc = ParallelDescriptor::MyProc();

    for (int k = 0; k < grids.length(); k++)
    {
        Box ndbox = ::surroundingNodes(grids[k]);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            const int* blo = ndbox.loVect();
            Box nd_lo(ndbox);
            nd_lo.setRange(dir,blo[dir],1);
            FabSet& lo = bndry[Orientation(dir,Orientation::low)];
            lo.setBox(k, nd_lo);
            if (lo.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!lo.defined(k));
                lo.clear(k);
                lo.setFab(k,new FArrayBox(nd_lo,1));
            }
            FabSet& lo_mask = bndry_mask[Orientation(dir,Orientation::low)];
            lo_mask.setBox(k, nd_lo);
            if (lo_mask.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!lo_mask.defined(k));
                lo_mask.clear(k);
                lo_mask.setFab(k,new FArrayBox(nd_lo,1));
            }
            const int* bhi = ndbox.hiVect();
            Box nd_hi(ndbox);
            nd_hi.setRange(dir,bhi[dir],1);
            FabSet& hi = bndry[Orientation(dir,Orientation::high)];
            hi.setBox(k, nd_hi);
            if (hi.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!hi.defined(k));
                hi.clear(k);
                hi.setFab(k,new FArrayBox(nd_hi,1));
            }
            FabSet& hi_mask = bndry_mask[Orientation(dir,Orientation::high)];
            hi_mask.setBox(k, nd_hi);
            if (hi_mask.DistributionMap()[k] == myproc)
            {
                //
                // Local.
                //
                assert(!hi_mask.defined(k));
                hi_mask.clear(k);
                hi_mask.setFab(k,new FArrayBox(nd_hi,1));
            }
            assert(ndbox.shortside() > 0);
        }
    }
}

SyncRegister::~SyncRegister () {}

Real
SyncRegister::sum ()
{
    //
    // Sum values in all registers.  Note that overlap is counted twice.
    //
    Box bb(grids[0]);
    for (int k = 1; k < grids.length(); k++)
    {
        bb.minBox(grids[k]);
    }
    bb.surroundingNodes();
    FArrayBox bfab(bb,1);
    bfab.setVal(0);

    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry[face()]); fsi.isValid(); ++fsi)
        {
            bfab.copy(fsi());
        }
    }
    Real tempSum = bfab.sum(0,1);
    ParallelDescriptor::ReduceRealSum(tempSum);
    return tempSum;
}

void
SyncRegister::increment (const FArrayBox& src)
{
    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry[face()]); fsi.isValid(); ++fsi)
        {
            fsi().plus(src);
        }
    }
}

void
SyncRegister::InitRHS (MultiFab&       rhs,
                       const Geometry& geom,
                       const BCRec*    phys_bc)
{
    rhs.setVal(0);

    const Box& domain = ::surroundingNodes(geom.Domain());

    const int MyProc = ParallelDescriptor::MyProc();

    FabSetCopyDescriptor fscd;

    Array<IntVect>    pshifts(27);
    vector<FillBoxId> fillBoxIDs;
    //
    // Fill Rhs From Bndry Registers.
    //
    // If periodic, copy the values from sync registers onto the nodes of the
    // rhs which are not covered by sync registers through periodic shifts.
    //
    if (geom.isAnyPeriodic())
    {
        vector<IntVect> shifts;

        for (OrientationIter face; face; ++face)
        {
            fscd.clear();
            shifts.clear();
            fillBoxIDs.clear();

            FabSetId faid = fscd.RegisterFabSet(&bndry[face()]);

            for (MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi)
            {
                for (int j = 0; j < grids.length(); j++)
                {
                    geom.periodicShift(domain,bndry[face()].fabbox(j),pshifts);

                    for (int iiv = 0; iiv < pshifts.length(); iiv++)
                    {
                        Box sbox = bndry[face()].fabbox(j) + pshifts[iiv];

                        if (sbox.intersects(mfi().box()))
                        {
                            sbox &= mfi().box();
                            sbox -= pshifts[iiv];

                            fillBoxIDs.push_back(fscd.AddBox(faid,
                                                             sbox,
                                                             0,
                                                             j,
                                                             0,
                                                             0,
                                                             mfi().nComp()));

                            assert(fillBoxIDs.back().box() == sbox);
                            //
                            // I need to save mfi.index() and pshifts[iiv].
                            //
                            fillBoxIDs.back().FabIndex(mfi.index());
                            //
                            // I'll maintain a parallel array for the IntVects.
                            //
                            shifts.push_back(pshifts[iiv]);
                        }
                    }
                }
            }

            fscd.CollectData();

            assert(fillBoxIDs.size() == shifts.size());

            for (int i = 0; i < fillBoxIDs.size(); i++)
            {
                assert(rhs.DistributionMap()[fillBoxIDs[i].FabIndex()] == MyProc);

                FArrayBox& fab = rhs[fillBoxIDs[i].FabIndex()];

                Box destbox = fillBoxIDs[i].box();

                destbox.shift(shifts[i]);

                fscd.FillFab(faid, fillBoxIDs[i], fab, destbox);
            }
        }
    }
    //
    // Overwrite above-set values on all nodes covered by a sync register.
    //
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyTo(rhs);
    }
    const int* dlo     = domain.loVect();
    const int* dhi     = domain.hiVect();
    const int* phys_lo = phys_bc->lo();
    const int* phys_hi = phys_bc->hi();

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            //
            // Any RHS point on the physical bndry must be multiplied by two
            // (only for ref-wall and inflow) and set to zero at outflow.
            //
            Box domlo(domain), domhi(domain);

            domlo.setRange(dir,dlo[dir],1);
            domhi.setRange(dir,dhi[dir],1);

            for (MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi)
            {
                if (domlo.intersects(mfi.validbox()))
                {
                    Box blo = mfi.validbox() & domlo;

                    if (phys_lo[dir] == Outflow)
                        mfi().setVal(0.0,blo,0,1);
                    else
                        mfi().mult(2.0,blo,0,1);
                }
                if (domhi.intersects(mfi.validbox()))
                {
                    Box bhi = mfi.validbox() & domhi;

                    if (phys_hi[dir] == Outflow)
                        mfi().setVal(0.0,bhi,0,1);
                    else
                        mfi().mult(2.0,bhi,0,1);
                }
            }
        } 
    }
    //
    // Set Up bndry_mask.
    //
    for (OrientationIter face; face; ++face)
    {
        bndry_mask[face()].setVal(0);
    }

    FArrayBox tmpfab;

    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry_mask[face()]); fsi.isValid(); ++fsi)
        {
            Box mask_cells = ::enclosedCells(::grow(fsi().box(),1));

            tmpfab.resize(mask_cells,1);
            tmpfab.setVal(0);

            for (int n = 0; n < grids.length(); n++)
            {
                if (mask_cells.intersects(grids[n]))
                {
                    tmpfab.setVal(1.0,(mask_cells & grids[n]),0,1);
                }
            }
 
            if (geom.isAnyPeriodic())
            {
                geom.periodicShift(geom.Domain(),mask_cells,pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    mask_cells += pshifts[iiv];

                    for (int n = 0; n < grids.length(); n++)
                    {
                        if (mask_cells.intersects(grids[n]))
                        {
                            Box intersect = mask_cells & grids[n];
                            intersect    -= pshifts[iiv];
                            tmpfab.setVal(1.0,intersect,0,1);
                        }
                    }

                    mask_cells -= pshifts[iiv];
                }
            }
            REAL* mask_dat = fsi().dataPtr();
            const int* mlo = fsi().loVect(); 
            const int* mhi = fsi().hiVect();
            REAL* cell_dat = tmpfab.dataPtr();
            const int* clo = tmpfab.loVect(); 
            const int* chi = tmpfab.hiVect();
        
            FORT_MAKEMASK(mask_dat,ARLIM(mlo),ARLIM(mhi),
                          cell_dat,ARLIM(clo),ARLIM(chi));
        }
    }

    Box node_domain = ::surroundingNodes(domain);
    //
    // Here double the cell contributions if at a non-periodic physical bdry.
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            Box domlo(node_domain), domhi(node_domain);

            domlo.setRange(dir,node_domain.loVect()[dir],1);
            domhi.setRange(dir,node_domain.hiVect()[dir],1);

            for (OrientationIter face; face; ++face)
            {
                for (FabSetIterator fsi(bndry_mask[face()]); fsi.isValid(); ++fsi)
                {
                    if (domlo.intersects(fsi().box()))
                    {
                        fsi().mult(2.0,fsi().box() & domlo,0,1);
                    }
                    if (domhi.intersects(fsi().box()))
                    {
                        fsi().mult(2.0,fsi().box() & domhi,0,1);
                    }
                }
            }
        }
    }
    //
    // Here convert from sum of cell contributions to 0 or 1.
    //
    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry_mask[face()]); fsi.isValid(); ++fsi)
        {
            REAL* mask_dat  = fsi().dataPtr();
            const int* mlo  = fsi().loVect(); 
            const int* mhi  = fsi().hiVect();

            FORT_CONVERTMASK(mask_dat,ARLIM(mlo),ARLIM(mhi));
        }
    }
    //
    // MULTIPLY RHS BY BNDRY_MASK
    //
    for (OrientationIter face; face; ++face)
    {
        fscd.clear();
        fillBoxIDs.clear();

        FabSetId faid = fscd.RegisterFabSet(&bndry_mask[face()]);

        for (MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi)
        {
            for (int j = 0; j < grids.length(); j++)
            {
                if (mfi().box().intersects(bndry_mask[face()].fabbox(j)))
                {
                    Box intersect = mfi().box() & bndry_mask[face()].fabbox(j);

                    fillBoxIDs.push_back(fscd.AddBox(faid,
                                                     intersect,
                                                     0,
                                                     j,
                                                     0,
                                                     0,
                                                     mfi().nComp()));

                    assert(fillBoxIDs.back().box() == intersect);
                    //
                    // Also save the index of our FAB needed filling.
                    //
                    fillBoxIDs.back().FabIndex(mfi.index());
                }
            }
        }

        fscd.CollectData();

        for (int i = 0; i < fillBoxIDs.size(); i++)
        {
            const FillBoxId& fbID = fillBoxIDs[i];

            assert(rhs.DistributionMap()[fbID.FabIndex()] == MyProc);

            FArrayBox& rhs_fab = rhs[fbID.FabIndex()];

            tmpfab.resize(fbID.box(), rhs_fab.nComp());

            fscd.FillFab(faid, fbID, tmpfab, fbID.box());

            rhs_fab.mult(tmpfab,fbID.box(),fbID.box(),0,0,rhs_fab.nComp());
        }
    }
}

static
void
Nullify (const BoxArray& grids,
         FArrayBox&      fab,
         const Box&      validbox,
         const Box&      subbox,
         int*            bc,
         int             ncomp,
         int             fine_index)
{
    fab.setVal(0,subbox,0,ncomp);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int bc_index = 2*BL_SPACEDIM*dir + dir;

        if (bc[bc_index] == EXT_DIR &&
            grids[fine_index].smallEnd(dir) == validbox.smallEnd(dir))
        {
            Box finesidelo(subbox);
            finesidelo.setRange(dir,finesidelo.smallEnd(dir)-1,1);
            fab.setVal(0,finesidelo,0,ncomp);
        }
        if (bc[bc_index+BL_SPACEDIM] == EXT_DIR &&
            grids[fine_index].bigEnd(dir) == validbox.bigEnd(dir))
        {
            Box finesidehi(subbox);
            finesidehi.setRange(dir,finesidehi.bigEnd(dir)+1,1);
            fab.setVal(0,finesidehi,0,ncomp);
        }
    }
}

void
SyncRegister::CrseDVInit (const MultiFab& U, 
                          const Geometry& geom, 
                          int             is_rz,
                          int**           crse_bc,
                          Real            mult)
{
    //
    // Zero all registers.
    //
    setVal(0);

    const int nghost = 1;

    MultiFab U_local(U.boxArray(),BL_SPACEDIM,nghost);
    //
    // First fill all the coarse cells, including ghost cells on periodic
    // and ext_dir edges, before worrying about zeroing out the ones under
    // fine grids.
    //
    for (MultiFabIterator mfi(U_local); mfi.isValid(); ++mfi)
    {
        //
        // U and U_local have same BoxArray and hence same DistributionMapping.
        //
        DependentMultiFabIterator dmfi_U(mfi, U);

        assert(dmfi_U.validbox() == mfi.validbox());

        mfi().setComplement(0,mfi.validbox(),0,BL_SPACEDIM);

        mfi().copy(dmfi_U(),dmfi_U.validbox(),0,dmfi_U.validbox(),0,BL_SPACEDIM);

        int* bc = crse_bc[mfi.index()];

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            int bc_index = 2*BL_SPACEDIM*dir + dir;
            //
            // Fill ghost cells outside of domain.
            // 
            if (bc[bc_index] == EXT_DIR &&
                dmfi_U.validbox().smallEnd(dir) == geom.Domain().smallEnd(dir))
            {
                Box sidelo(dmfi_U.validbox());
                sidelo.setRange(dir,sidelo.smallEnd(dir)-1,1);
                mfi().copy(dmfi_U(),sidelo,dir,sidelo,dir,1);
            }
            if (bc[bc_index+BL_SPACEDIM] == EXT_DIR &&
                dmfi_U.validbox().bigEnd(dir) == geom.Domain().bigEnd(dir))
            {
                Box sidehi(dmfi_U.validbox());
                sidehi.setRange(dir,sidehi.bigEnd(dir)+1,1);
                mfi().copy(dmfi_U(),sidehi,dir,sidehi,dir,1);
            }
        }
    }

    geom.FillPeriodicBoundary(U_local,true,false);
    //
    // Now do all the zeroing-out associated with the fine grids, on the
    // interior and on ghost cells on periodic and ext_dir edges.
    //
    Array<IntVect> pshifts(27);

    for (MultiFabIterator mfi(U_local); mfi.isValid(); ++mfi)
    {
        assert(mfi.validbox() == U.boxArray()[mfi.index()]);

        int* bc = crse_bc[mfi.index()];

        for (int fine = 0; fine < grids.length(); fine++)
        {
            if (grids[fine].intersects(mfi.validbox()))
            {
                Box subbox = mfi.validbox() & grids[fine];
                Nullify(grids,mfi(),mfi.validbox(),subbox,bc,BL_SPACEDIM,fine);
            }
            //
            // Now zero out under periodic translations of fine grids.
            //
            if (geom.isAnyPeriodic())
            {
                Box domain_plus(geom.Domain());

                for (int dir = 0; dir < BL_SPACEDIM; dir++) 
                    if (geom.isPeriodic(dir))
                        domain_plus.grow(dir,1);

                geom.periodicShift(domain_plus, grids[fine], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box fine_shifted = grids[fine] + pshifts[iiv];

                    if (fine_shifted.intersects(mfi().box()))
                    {
                        fine_shifted &= mfi().box();
                        Nullify(grids,mfi(),mfi.validbox(),fine_shifted,bc,BL_SPACEDIM,fine);
                    }
                }
            }
        }
    }
    //
    // Now compute node-centered divergence.
    //
    FArrayBox divu;

    for (MultiFabIterator mfi(U_local); mfi.isValid(); ++mfi)
    {
        const int* ulo  = mfi().loVect();
        const int* uhi  = mfi().hiVect();
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        divu.resize(ndbox,1);

        FORT_SRDIVU(mfi().dataPtr(),ARLIM(ulo),ARLIM(uhi),
                    divu.dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
                    ndlo,ndhi,geom.CellSize(),&mult,&is_rz);

        increment(divu);
    }
}

static
void
FillExtDir (FArrayBox&       dstfab,
            const FArrayBox& srcfab,
            const Box&       validbox,
            int*             bc)
{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int bc_index = 2*BL_SPACEDIM*dir + dir;

        if (bc[bc_index] == EXT_DIR)
        {
            Box sidelo(validbox);
            sidelo.setRange(dir,sidelo.smallEnd(dir)-1,1);
            dstfab.copy(srcfab,sidelo,dir,sidelo,dir,1);
        }
        if (bc[bc_index+BL_SPACEDIM] == EXT_DIR)
        {
            Box sidehi(validbox);
            sidehi.setRange(dir,sidehi.bigEnd(dir)+1,1);
            dstfab.copy(srcfab,sidehi,dir,sidehi,dir,1);
        }
    }
}

void
SyncRegister::incrementPeriodic (const Geometry& crse_geom,
                                 const Box&      crse_node_domain,
                                 const Box&      cboxlo,
                                 const Box&      cboxhi,
                                 FArrayBox&      cfablo,
                                 FArrayBox&      cfabhi,
                                 Array<IntVect>& pshifts)
{
    crse_geom.periodicShift(crse_node_domain, cboxlo, pshifts);

    for (int iiv = 0; iiv < pshifts.length(); iiv++)
    {
        cfablo.shift(pshifts[iiv]);
        increment(cfablo);
        cfablo.shift(-pshifts[iiv]);
    }

    crse_geom.periodicShift(crse_node_domain, cboxhi, pshifts);

    for (int iiv = 0; iiv < pshifts.length(); iiv++)
    {
        cfabhi.shift(pshifts[iiv]);
        increment(cfabhi);
        cfabhi.shift(-pshifts[iiv]);
    }
}

static
void
SizeTmpFabs (Box&           reglo,
             Box&           reghi,
             Box&           tboxlo,
             Box&           tboxhi,
             FArrayBox&     ffablo,
             FArrayBox&     ffabhi,
             Box&           cboxlo,
             Box&           cboxhi,
             FArrayBox&     cfablo,
             FArrayBox&     cfabhi,
             int            dir,
             const IntVect& ratio)
{
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
        if (idir < dir)
        {
            //
            // Previous direction have included these
            // points in the stencil already, shrink.
            //
            reglo.grow(idir,-1);
            reghi.grow(idir,-1);
        }
        if (idir != dir)
        {
            //
            // Need additional room for stencil calculation.
            //
            tboxlo.grow(idir,ratio[idir]-1);
            tboxhi.grow(idir,ratio[idir]-1);
        }
    }
    //
    // Define fine grid tmp fabs.
    //
    ffablo.resize(tboxlo,1);
    ffabhi.resize(tboxhi,1);
    ffablo.setVal(0);
    ffabhi.setVal(0);
    //
    // Define coarsened tmp fabs.
    //
    cboxlo.coarsen(ratio);
    cboxhi.coarsen(ratio);
    cfablo.resize(cboxlo,1);
    cfabhi.resize(cboxhi,1);
    cfablo.setVal(0);
    cfabhi.setVal(0);
}

void
SyncRegister::FineDVAdd (const MultiFab& U, 
                         const Real*     dx_fine, 
                         const Geometry& crse_geom, 
                         int             is_rz,
                         int**           fine_bc,
                         Real            mult)
{
    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    FArrayBox ufab, cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (ConstMultiFabIterator mfi(U); mfi.isValid(); ++mfi)
    {
        ufab.resize(::grow(mfi.validbox(),1),BL_SPACEDIM);
        ufab.setComplement(0,mfi.validbox(),0,BL_SPACEDIM);
        ufab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,BL_SPACEDIM);

        const int* ulo = ufab.box().loVect();
        const int* uhi = ufab.box().hiVect();
        int* bc        = fine_bc[mfi.index()];

        FillExtDir(ufab, mfi(), mfi.validbox(), bc);
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,
                        cboxlo,cboxhi,cfablo,cfabhi,dir,ratio);

            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            //
            // Compute divu on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                        reglo.loVect(),reglo.hiVect(),
                        dx_fine,&mult,&is_rz);
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                        reghi.loVect(),reghi.hiVect(),
                        dx_fine,&mult,&is_rz);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                incrementPeriodic(crse_geom,
                                  crse_node_domain,
                                  cboxlo,
                                  cboxhi,
                                  cfablo,
                                  cfabhi,
                                  pshifts);
            }
        }
    }
}

static
void
SetCenter (int             is_rz,
           Array<Real>&    rcen,
           const Geometry& geom,
           const Box&      dsdtbox)
{
    if (is_rz)
    {
        geom.GetCellLoc(rcen, dsdtbox, 0);
    }
    else
    {
        for (int i = 0; i < rcen.length(); i++)
            rcen[i] = 1.0;
    }
}

void
SyncRegister::CrseDsdtAdd (const MultiFab& dsdt,
                           const Geometry& geom,
                           int             is_rz,
                           int**           crse_bc, 
                           int             lowfix,
                           int             hifix,
                           Real            mult)
{
    FArrayBox divu;

    Array<IntVect> pshifts(27);

    const int nghost = 1;

    MultiFab dsdt_local(dsdt.boxArray(),1,nghost);
    //
    // Fill valid region of dsdt_local from dsdt and zero out ghost cells.
    //
    for (MultiFabIterator mfi(dsdt_local); mfi.isValid(); ++mfi)
    {
        mfi().setComplement(0,mfi.validbox(),0,1);
        mfi().copy(dsdt[mfi.index()], mfi.validbox());
    }

    geom.FillPeriodicBoundary(dsdt_local,true,false);

    for (MultiFabIterator mfi(dsdt_local); mfi.isValid(); ++mfi)
    {
        int* bc = crse_bc[mfi.index()];
        //
        // Zero out coarse contributions under fine grid.
        // Also, zero out contrib if dsdt box is on solid wall
        //
        for (int fine = 0; fine < grids.length(); fine++)
        {
            if (grids[fine].intersects(mfi().box()))
            {
                Box subbox = mfi().box() & grids[fine];

                Nullify(grids,mfi(),mfi.validbox(),subbox,bc,1,fine);
            }
        }
        //
        // Zero ghost cells under periodic fine grid.
        //
        geom.periodicShift(geom.Domain(), mfi().box(), pshifts);

        for (int iiv = 0; iiv < pshifts.length(); iiv++)
        {
            Box dsdtbox = mfi().box() + pshifts[iiv];

            for (int fine = 0; fine < grids.length(); fine++)
            {
                if (dsdtbox.intersects(grids[fine]))
                {
                    Box ovlp = dsdtbox & grids[fine];
                    ovlp    -= pshifts[iiv];
                    mfi().setVal(0,ovlp,0,1);
                }
            }
        }
        //
        // Average dsdt to nodes.
        //
        divu.resize(::surroundingNodes(mfi.validbox()),1);

        const int* ndlo   = divu.box().loVect();
        const int* ndhi   = divu.box().hiVect();
        const int* dsdtlo = mfi().box().loVect();
        const int* dsdthi = mfi().box().hiVect();
        const int* rlo    = mfi().box().loVect();
        const int* rhi    = mfi().box().hiVect();
        const int* domlo  = geom.Domain().loVect();
        const int* domhi  = geom.Domain().hiVect();

        Array<Real> rcen(mfi().box().length(0));

        SetCenter(is_rz, rcen, geom, mfi().box());

#if (BL_SPACEDIM==2)
        int nghost         = 0;
        Real hx            = geom.CellSize()[0];
        int extrap_edges   = 0;
        int extrap_corners = 0;
        FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                   mfi().dataPtr(),
                   rcen.dataPtr(), 
                   ARLIM(ndlo), ARLIM(ndhi), divu.dataPtr(), 
                   domlo, domhi, lowfix, hifix, &hx,
                   &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
        //
        // TODO -- make 3-D work !!!
        //
        divu.setVal(0);
#endif
        divu.negate();
        divu.mult(mult);
        increment(divu);
    }
}

void
SyncRegister::FineDsdtAdd (const MultiFab& dsdt,
                           const Geometry& geom,
                           const Geometry& crse_geom,
                           int             is_rz,
                           int**           fine_bc, 
                           int             lowfix,
                           int             hifix,
                           Real            mult)
{
    FArrayBox dsdtfab, cfablo, cfabhi, ffablo, ffabhi;
    FArrayBox ffablo_tmp, ffabhi_tmp;

    Array<IntVect> pshifts(27);

    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (ConstMultiFabIterator mfi(dsdt); mfi.isValid(); ++mfi)
    {
        dsdtfab.resize(::grow(mfi.validbox(),1),1);
        dsdtfab.setComplement(0,mfi.validbox(),0,1);
        dsdtfab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,1);

        const int* dsdtlo = dsdtfab.box().loVect();
        const int* dsdthi = dsdtfab.box().hiVect();
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,
                        cboxlo,cboxhi,cfablo,cfabhi,dir,ratio);

            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            //
            // Average dsdt to nodes on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            const int* dsdtlo = dsdtfab.box().loVect();
            const int* dsdthi = dsdtfab.box().hiVect();
            const int* rlo    = dsdtfab.box().loVect();
            const int* rhi    = dsdtfab.box().hiVect();
            const int* domlo  = geom.Domain().loVect();
            const int* domhi  = geom.Domain().hiVect();

            Array<Real> rcen(dsdtfab.box().length(0));

            SetCenter(is_rz, rcen, geom, dsdtfab.box());

            ffablo_tmp.resize(reglo,1);
#if (BL_SPACEDIM==2)
            int nghost         = 0;
            int hi_fix         = 0;
            Real hx            = geom.CellSize()[0];
            int extrap_edges   = 0;
            int extrap_corners = 0;
            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reglo.loVect()), ARLIM(reglo.hiVect()), 
                       ffablo_tmp.dataPtr(),
                       domlo, domhi, lowfix, hi_fix, &hx,
                       &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
            //
            // TODO -- make 3-D work !!!
            //
            ffablo_tmp.setVal(0);
#endif
            ffablo_tmp.negate();
            ffablo_tmp.mult(mult);
            ffablo.copy(ffablo_tmp);
            ffabhi_tmp.resize(reghi,1);
#if (BL_SPACEDIM==2)
            int low_fix = 0;
            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reghi.loVect()), ARLIM(reghi.hiVect()), 
                       ffabhi_tmp.dataPtr(), 
                       domlo, domhi, low_fix, hifix, &hx,
                       &extrap_edges, &extrap_corners, &is_rz);
#elif (BL_SPACEDIM==3)
            //
            // TODO -- make 3-D work !!!
            //
            ffabhi_tmp.setVal(0);
#endif
            ffabhi_tmp.negate();
            ffabhi_tmp.mult(mult);
            ffabhi.copy(ffabhi_tmp);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                incrementPeriodic(crse_geom,
                                  crse_node_domain,
                                  cboxlo,
                                  cboxhi,
                                  cfablo,
                                  cfabhi,
                                  pshifts);
            }
        }
    }
}

static
void
NullOverlap (const BoxArray& Pgrids,
             const Geometry& fine_geom,
             const Box&      reglo,
             FArrayBox&      ffablo,
             const Box&      reghi,
             FArrayBox&      ffabhi,
             Array<IntVect>& pshifts)
{
    const int n_comp   = 1;
    const int set_comp = 0;

    for (int i = 0; i < Pgrids.length(); i++)
    {
        if (reglo.intersects(Pgrids[i]))
        {
            ffablo.setVal(0,(reglo & Pgrids[i]),set_comp,n_comp);
        }

        fine_geom.periodicShift(reglo, Pgrids[i], pshifts);

        for (int iiv = 0; iiv < pshifts.length(); iiv++)
        {
            Box overlap_lo_per = Pgrids[i] + pshifts[iiv];
            overlap_lo_per    &= reglo;
            ffablo.setVal(0,overlap_lo_per,set_comp,n_comp);
        }

        if (reghi.intersects(Pgrids[i]))
        {
            ffabhi.setVal(0,(reghi & Pgrids[i]),set_comp,n_comp);
        }

        fine_geom.periodicShift(reghi, Pgrids[i], pshifts);

        for (int iiv = 0; iiv < pshifts.length(); iiv++)
        {
            Box overlap_hi_per = Pgrids[i] + pshifts[iiv];
            overlap_hi_per    &= reghi;
            ffabhi.setVal(0,overlap_hi_per,set_comp,n_comp);
        }
    }
}

void
SyncRegister::CompDVAdd (const MultiFab& U, 
                         const BoxArray& Pgrids,
                         const Real*     dx_fine, 
                         const Geometry& fine_geom, 
                         const Geometry& crse_geom, 
                         int             is_rz,
                         int**           fine_bc,
                         Real            mult)
{
    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    Array<IntVect> pshifts(27);

    FArrayBox ufab, cfablo, cfabhi, ffablo, ffabhi;

    for (ConstMultiFabIterator mfi(U); mfi.isValid(); ++mfi)
    {
        ufab.resize(::grow(mfi.validbox(),1),BL_SPACEDIM);
        ufab.setComplement(0,mfi.validbox(),0,BL_SPACEDIM);
        ufab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,BL_SPACEDIM);

        const int* ulo = ufab.box().loVect();
        const int* uhi = ufab.box().hiVect();
        int* bc        = fine_bc[mfi.index()];

        FillExtDir(ufab, mfi(), mfi.validbox(), bc);
        //
        // Now compute node centered surrounding box.
        //
        Box ndbox       = ::surroundingNodes(mfi.validbox());
        const int* ndlo = ndbox.loVect();
        const int* ndhi = ndbox.hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,
                        cboxlo,cboxhi,cfablo,cfabhi,dir,ratio);

            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            //
            // Compute divu on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                        reglo.loVect(),reglo.hiVect(),dx_fine,&mult,&is_rz);
            FORT_SRDIVU(ufab.dataPtr(),ARLIM(ulo),ARLIM(uhi),
                        ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                        reghi.loVect(),reghi.hiVect(),dx_fine,&mult,&is_rz);

            NullOverlap(Pgrids,fine_geom,reglo,ffablo,reghi,ffabhi,pshifts);
            //
            // Coarsen edge values.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                incrementPeriodic(crse_geom,
                                  crse_node_domain,
                                  cboxlo,
                                  cboxhi,
                                  cfablo,
                                  cfabhi,
                                  pshifts);
            }
        }
    }
}

void
SyncRegister::CrseLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Geometry& geom,
                           int             is_rz,
                           Real            mult)
{
    const int nghost = 1;
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    assert(Phi.boxArray().length() == Sigma.boxArray().length());

    MultiFab Sig_local(Sigma.boxArray(),1,nghost);
    MultiFab Phi_local(Phi.boxArray(),1,nghost);
    //
    // Copy valid region of Phi into Phi_local & zero out ghost cells.
    //
    for (MultiFabIterator mfi(Phi_local); mfi.isValid(); ++mfi)
    {
        mfi().setComplement(0,mfi.validbox(),0,1);
        mfi().copy(Phi[mfi.index()], mfi.validbox());
    }
    //
    // Copy valid region of Sigma into Sig_local & zero out ghost cells.
    // Also, zero out region covered by fine grid.
    //
    for (MultiFabIterator mfi(Sig_local); mfi.isValid(); ++mfi)
    {
        mfi().setComplement(0,mfi.validbox(),0,1);
        mfi().copy(Sigma[mfi.index()], mfi.validbox());

        for (int fine = 0; fine < grids.length(); fine++)
        {
            if (grids[fine].intersects(mfi.validbox()))
            {
                mfi().setVal(0,grids[fine] & mfi.validbox(),0,1);
            }
        }
    }

    geom.FillPeriodicBoundary(Sig_local, 0, 1, true, false);

    geom.FillPeriodicBoundary(Phi_local, 0, 1, false, false);
    //
    // Now compute node centered div(Sigma*grad(PHI)).
    //
    FArrayBox divgp;

    for (MultiFabIterator pmfi(Phi_local); pmfi.isValid(); ++pmfi)
    {
        DependentMultiFabIterator smfi(pmfi, Sig_local);

        const int* slo  = smfi().loVect();
        const int* shi  = smfi().hiVect();
        const int* p_lo = pmfi().loVect();
        const int* p_hi = pmfi().hiVect();
        divgp.resize(pmfi.validbox(),1);
        const int* glo  = divgp.loVect();
        const int* ghi  = divgp.hiVect();

        FORT_SRDGPHI(pmfi().dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                     smfi().dataPtr(),ARLIM(slo),ARLIM(shi),
                     divgp.dataPtr(),ARLIM(glo),ARLIM(ghi),
                     glo,ghi,geom.CellSize(),&mult,&is_rz);

        increment(divgp);
    }
}

void
SyncRegister::FineLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Real*     dx_fine,
                           const Geometry& crse_geom, 
                           int             is_rz,
                           Real            mult)
{
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    assert(Phi.boxArray().length() == Sigma.boxArray().length());

    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    FArrayBox pfab, sfab, cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (ConstMultiFabIterator pmfi(Phi); pmfi.isValid(); ++pmfi)
    {
        pfab.resize(::grow(pmfi.validbox(),1),1);
        pfab.setComplement(0,pmfi.validbox(),0,1);
        pfab.copy(pmfi(),pmfi.validbox());

        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();

        ConstDependentMultiFabIterator smfi(pmfi, Sigma);

        sfab.resize(::grow(smfi.validbox(),1),1);
        sfab.setComplement(0,smfi.validbox(),0,1);
        sfab.copy(smfi(),smfi.validbox());

        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();
        const int* ndlo    = pmfi.validbox().loVect();
        const int* ndhi    = pmfi.validbox().hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(pmfi.validbox()), tboxhi(pmfi.validbox());
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,
                        cboxlo,cboxhi,cfablo,cfabhi,dir,ratio);

            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            //
            // Compute divgp on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                         reglo.loVect(),reglo.hiVect(),
                         dx_fine,&mult,&is_rz);
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                         reghi.loVect(),reghi.hiVect(),
                         dx_fine,&mult,&is_rz);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);

            if (crse_geom.isAnyPeriodic())
            {
                incrementPeriodic(crse_geom,
                                  crse_node_domain,
                                  cboxlo,
                                  cboxhi,
                                  cfablo,
                                  cfabhi,
                                  pshifts);
            }
        }
    }
}

void
SyncRegister::CompLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const BoxArray& Pgrids, 
                           const Real*     dx_fine, 
                           const Geometry& fine_geom, 
                           const Geometry& crse_geom, 
                           int             is_rz,
                           Real            mult)
{
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    assert(Phi.boxArray().length() == Sigma.boxArray().length());

    const Box& crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    FArrayBox pfab, sfab, cfablo, cfabhi, ffablo, ffabhi;

    Array<IntVect> pshifts(27);

    for (ConstMultiFabIterator pmfi(Phi); pmfi.isValid(); ++pmfi)
    {
        ConstDependentMultiFabIterator smfi(pmfi, Sigma);

        pfab.resize(::grow(pmfi.validbox(),1),1);
        pfab.setComplement(0,pmfi.validbox(),0,1);
        pfab.copy(pmfi(),pmfi.validbox());
        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();

        sfab.resize(::grow(smfi.validbox(),1),1);
        sfab.setComplement(0,smfi.validbox(),0,1);
        sfab.copy(smfi(),smfi.validbox());
        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();
        const int* ndlo    = pmfi.validbox().loVect();
        const int* ndhi    = pmfi.validbox().hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            //
            // Determine region of interest, and size of tmp fabs.
            //
            Box tboxlo(pmfi.validbox()), tboxhi(pmfi.validbox());
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box cboxlo(tboxlo), cboxhi(tboxhi);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,
                        cboxlo,cboxhi,cfablo,cfabhi,dir,ratio);

            const int* flo_lo = tboxlo.loVect();
            const int* flo_hi = tboxlo.hiVect();
            const int* fhi_lo = tboxhi.loVect();
            const int* fhi_hi = tboxhi.hiVect();
            //
            // Compute divgp on fine grid edges in regions defined
            // by reglo and reghi.  Fabs are set to zero outside region.
            //
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                         reglo.loVect(),reglo.hiVect(),
                         dx_fine,&mult,&is_rz);
            FORT_SRDGPHI(pfab.dataPtr(),ARLIM(pfab_lo),ARLIM(pfab_hi),
                         sfab.dataPtr(),ARLIM(sfab_lo),ARLIM(sfab_hi),
                         ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                         reghi.loVect(),reghi.hiVect(),
                         dx_fine,&mult,&is_rz);

            NullOverlap(Pgrids,fine_geom,reglo,ffablo,reghi,ffabhi,pshifts);
            //
            // Coarsen edge value.
            //
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();
            FORT_SRCRSEREG(ffablo.dataPtr(),ARLIM(flo_lo),ARLIM(flo_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();
            FORT_SRCRSEREG(ffabhi.dataPtr(),ARLIM(fhi_lo),ARLIM(fhi_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
            //
            // Intersect and add to registers.
            //
            increment(cfablo);
            increment(cfabhi);
            //
            // Intersect and add to registers.
            //
            if (crse_geom.isAnyPeriodic())
            {
                incrementPeriodic(crse_geom,
                                  crse_node_domain,
                                  cboxlo,
                                  cboxhi,
                                  cfablo,
                                  cfabhi,
                                  pshifts);
            }
        }
    }
}
