//BL_COPYRIGHT_NOTICE

//
// $Id: SyncRegister.cpp,v 1.58 1999-06-30 22:40:04 almgren Exp $
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

//#define BL_SYNC_DEBUG 1

#ifdef BL_SYNC_DEBUG
static
Real
Norm_0 (const SyncRegister& br)
{
    Real r = 0;
    for (OrientationIter face; face; ++face)
    {
        for (ConstFabSetIterator cmfi(br[face()]); cmfi.isValid(); ++cmfi)
        {
            Real s = cmfi->norm(cmfi->box(), 0, 0, cmfi->nComp());
            r = (r > s) ? r : s;
        }
    }
    ParallelDescriptor::ReduceRealMax(r,ParallelDescriptor::IOProcessorNumber());
    return r;
}

static
Real
Norm_1 (const SyncRegister& br)
{
    Real r = 0;
    for (OrientationIter face; face; ++face)
    {
        for (ConstFabSetIterator cmfi(br[face()]); cmfi.isValid(); ++cmfi)
        {
            r += cmfi->norm(cmfi->box(), 1, 0, cmfi->nComp());
        }
    }
    ParallelDescriptor::ReduceRealSum(r,ParallelDescriptor::IOProcessorNumber());
    return r;
}

static
void
WriteNorms (const SyncRegister& sr,
            const char*         msg)
{
    Real n_0 = Norm_0(sr);
    Real n_1 = Norm_1(sr);
    if (ParallelDescriptor::IOProcessor())
    {
        cout << "*** " << msg << ": Norm_0: " << n_0 << endl;
        cout << "*** " << msg << ": Norm_1: " << n_1 << endl;
    }
}
#else
inline
void
WriteNorms (const SyncRegister&, const char*)
{}
#endif /*BL_SYNC_DEBUG*/

//
// Structure used by some SyncRegister functions.
//

struct SRRec
{
    SRRec (Orientation    face,
           int            index,
           const IntVect& shift)
        :
        m_shift(shift),
        m_idx(index),
        m_face(face)
    {}

    SRRec (Orientation face,
           int         index)
        :
        m_idx(index),
        m_face(face)
    {}

    FillBoxId   m_fbid;
    IntVect     m_shift;
    int         m_idx;
    Orientation m_face;
};

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
        BL_ASSERT(ratio[dir] == -1);

    BL_ASSERT(!grids.ready());
    BL_ASSERT(fine_boxes.isDisjoint());

    ratio      = ref_ratio;
    fine_level = fine_lev;
 
    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        //
        // Construct BoxArrays for the FabSets.
        //
        const Orientation lo = Orientation(dir,Orientation::low);
        const Orientation hi = Orientation(dir,Orientation::high);

        BoxArray loBA(grids.length());
        BoxArray hiBA(grids.length());

        for (int k = 0; k < grids.length(); k++)
        {
            Box ndbox = ::surroundingNodes(grids[k]);

            Box nd_lo = ndbox;
            nd_lo.setRange(dir,ndbox.smallEnd(dir),1);
            loBA.set(k,nd_lo);

            Box nd_hi = ndbox;
            nd_hi.setRange(dir,ndbox.bigEnd(dir),1);
            hiBA.set(k,nd_hi);
        }
        //
        // Define the FabSets.
        //
        bndry[lo].define(loBA,1);
        bndry_mask[lo].define(loBA,1);
        bndry[hi].define(hiBA,1);
        bndry_mask[hi].define(hiBA,1);
    }
}

SyncRegister::~SyncRegister () {}

Real
SyncRegister::sum ()
{
    //
    // Sum values in all registers.  Note that overlap is counted twice.
    //
    Box bb = grids[0];

    for (int k = 1; k < grids.length(); k++)
        bb.minBox(grids[k]);

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

//
// If periodic, copy the values from sync registers onto the nodes of the
// rhs which are not covered by sync registers through periodic shifts.
//

void
SyncRegister::copyPeriodic (const Geometry& geom,
                            const Box&      domain,
                            MultiFab&       rhs) const
{
    if (!geom.isAnyPeriodic()) return;

    const int            MyProc = ParallelDescriptor::MyProc();
    Array<IntVect>       pshifts(27);
    vector<SRRec>        srrec;
    FabSetCopyDescriptor fscd;
    FabSetId             fsid[2*BL_SPACEDIM];

    for (OrientationIter face; face; ++face)
    {
        fsid[face()] = fscd.RegisterFabSet((FabSet*) &bndry[face()]);

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

                        SRRec sr(face(), mfi.index(), pshifts[iiv]);

                        sr.m_fbid = fscd.AddBox(fsid[face()],
                                                sbox,
                                                0,
                                                j,
                                                0,
                                                0,
                                                mfi().nComp());

                        BL_ASSERT(sr.m_fbid.box() == sbox);

                        srrec.push_back(sr);
                    }
                }
            }
        }
    }

    fscd.CollectData();

    for (int i = 0; i < srrec.size(); i++)
    {
        BL_ASSERT(rhs.DistributionMap()[srrec[i].m_idx] == MyProc);

        FArrayBox& fab = rhs[srrec[i].m_idx];

        Box dstbox = srrec[i].m_fbid.box() + srrec[i].m_shift;

        fscd.FillFab(fsid[srrec[i].m_face], srrec[i].m_fbid, fab, dstbox);
    }
}

void
SyncRegister::multByBndryMask (MultiFab& rhs) const
{
    const int            MyProc = ParallelDescriptor::MyProc();
    FabSetCopyDescriptor fscd;
    FabSetId             fsid[2*BL_SPACEDIM];
    vector<SRRec>        srrec;
    FArrayBox            tmpfab;

    for (OrientationIter face; face; ++face)
    {
        fsid[face()] = fscd.RegisterFabSet((FabSet*) &bndry_mask[face()]);

        for (MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi)
        {
            for (int j = 0; j < grids.length(); j++)
            {
                if (mfi().box().intersects(bndry_mask[face()].fabbox(j)))
                {
                    Box intersect = mfi().box() & bndry_mask[face()].fabbox(j);

                    SRRec sr(face(), mfi.index());

                    sr.m_fbid = fscd.AddBox(fsid[face()],
                                            intersect,
                                            0,
                                            j,
                                            0,
                                            0,
                                            mfi().nComp());

                    BL_ASSERT(sr.m_fbid.box() == intersect);

                    srrec.push_back(sr);
                }
            }
        }
    }

    fscd.CollectData();

    for (int i = 0; i < srrec.size(); i++)
    {
        const Box& bx = srrec[i].m_fbid.box();

        BL_ASSERT(rhs.DistributionMap()[srrec[i].m_idx] == MyProc);

        FArrayBox& rhs_fab = rhs[srrec[i].m_idx];

        tmpfab.resize(bx, rhs_fab.nComp());

        fscd.FillFab(fsid[srrec[i].m_face], srrec[i].m_fbid, tmpfab, bx);

        rhs_fab.mult(tmpfab, bx, bx, 0, 0, rhs_fab.nComp());
    }
}

void
SyncRegister::InitRHS (MultiFab&       rhs,
                       const Geometry& geom,
                       const BCRec*    phys_bc)
{
    rhs.setVal(0);

    Box domain = ::surroundingNodes(geom.Domain());

    const int MyProc = ParallelDescriptor::MyProc();

    FabSetCopyDescriptor fscd;

    Array<IntVect> pshifts(27);
    //
    // Fill Rhs From Bndry Registers.
    //
    // If periodic, copy the values from sync registers onto the nodes of the
    // rhs which are not covered by sync registers through periodic shifts.
    //
    copyPeriodic(geom,domain,rhs);
    //
    // Overwrite above-set values on all nodes covered by a sync register.
    //
    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyTo(rhs);
    }
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

            domlo.setRange(dir,domain.smallEnd(dir),1);
            domhi.setRange(dir,domain.bigEnd(dir),1);

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

    ofstream os("old.fab");
    rhs[0].writeOn(os);
    cout << "RHS " << endl;
    cout << rhs[0] << endl;
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
                if (mask_cells.intersects(grids[n]))
                    tmpfab.setVal(1.0,(mask_cells & grids[n]),0,1);
 
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
            Real* mask_dat = fsi().dataPtr();
            const int* mlo = fsi().loVect(); 
            const int* mhi = fsi().hiVect();
            Real* cell_dat = tmpfab.dataPtr();
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

            domlo.setRange(dir,node_domain.smallEnd(dir),1);
            domhi.setRange(dir,node_domain.bigEnd(dir),1);

            for (OrientationIter face; face; ++face)
            {
                for (FabSetIterator fsi(bndry_mask[face()]); fsi.isValid(); ++fsi)
                {
                    if (domlo.intersects(fsi().box()))
                        fsi().mult(2.0,fsi().box() & domlo,0,1);

                    if (domhi.intersects(fsi().box()))
                        fsi().mult(2.0,fsi().box() & domhi,0,1);
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
            Real* mask_dat  = fsi().dataPtr();
            const int* mlo  = fsi().loVect(); 
            const int* mhi  = fsi().hiVect();

            FORT_CONVERTMASK(mask_dat,ARLIM(mlo),ARLIM(mhi));
        }
    }

    multByBndryMask(rhs);

    WriteNorms(*this,"SyncRegister::InitRHS(E)");
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
    setVal(0);

    MultiFab U_local(U.boxArray(),BL_SPACEDIM,1);
    //
    // First fill all the coarse cells, including ghost cells on periodic
    // and ext_dir edges, before worrying about zeroing out the ones under
    // fine grids.
    //
    for (MultiFabIterator mfi(U_local); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi_U(mfi, U);

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
        BL_ASSERT(mfi.validbox() == U.boxArray()[mfi.index()]);

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
    // divuMF and U_local will have same distribution mapping.
    //
    BoxArray divuBA = U_local.boxArray();

    divuBA.surroundingNodes();

    MultiFab divuMF(divuBA,1,0);

    for (MultiFabIterator mfi(U_local); mfi.isValid(); ++mfi)
    {
        const int* ulo   = mfi().loVect();
        const int* uhi   = mfi().hiVect();
        FArrayBox& divu  = divuMF[mfi.index()];
        const Box& ndbox = divu.box();
        const int* ndlo  = ndbox.loVect();
        const int* ndhi  = ndbox.hiVect();

        FORT_SRDIVU(mfi().dataPtr(),ARLIM(ulo),ARLIM(uhi),
                    divu.dataPtr(),ARLIM(ndlo),ARLIM(ndhi),
                    ndlo,ndhi,geom.CellSize(),&mult,&is_rz);
    }

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(divuMF,0,0,0,1);
    }

    WriteNorms(*this,"SyncRegister::CrseDVInit(E)");
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

static
void
SizeTmpFabs (Box&           reglo,
             Box&           reghi,
             Box&           tboxlo,
             Box&           tboxhi,
             FArrayBox&     ffablo,
             FArrayBox&     ffabhi,
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
}

enum HowToBuild { WITH_BOX, WITH_SURROUNDING_BOX };

static
void
BuildMFs (const MultiFab& mf,
          MultiFab        cloMF[BL_SPACEDIM],
          MultiFab        chiMF[BL_SPACEDIM],
          const IntVect&  ratio,
          HowToBuild      how)
{
    BL_ASSERT(how == WITH_BOX || how == WITH_SURROUNDING_BOX);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        BoxArray cloBA(mf.length());
        BoxArray chiBA(mf.length());

        for (int i = 0; i < mf.length(); i++)
        {
            Box bx = mf.boxArray()[i];

            if (how == WITH_SURROUNDING_BOX)
                bx = ::surroundingNodes(mf.boxArray()[i]);

            Box tboxlo(bx), tboxhi(bx);

            tboxlo.setRange(dir,bx.smallEnd(dir),1);
            tboxhi.setRange(dir,bx.bigEnd(dir),1);

            cloBA.set(i,::coarsen(tboxlo,ratio));
            chiBA.set(i,::coarsen(tboxhi,ratio));
        }

        cloMF[dir].define(cloBA,1,0,Fab_allocate);
        chiMF[dir].define(chiBA,1,0,Fab_allocate);

        cloMF[dir].setVal(0);
        chiMF[dir].setVal(0);
    }
}

void
SyncRegister::incrementPeriodic (const Geometry& geom,
                                 const Box&      domain,
                                 const MultiFab& mf)
{
    if (!geom.isAnyPeriodic()) return;

    const int              MyProc = ParallelDescriptor::MyProc();
    const BoxArray&        mfBA   = mf.boxArray();
    Array<IntVect>         pshifts(27);
    vector<SRRec>          srrec;
    FArrayBox              tmpfab;
    MultiFabCopyDescriptor mfcd;
    MultiFabId             mfid = mfcd.RegisterFabArray((MultiFab*) &mf);

    for (OrientationIter face; face; ++face)
    {
        for (FabSetIterator fsi(bndry[face()]); fsi.isValid(); ++fsi)
        {
            for (int j = 0; j < mfBA.length(); j++)
            {
                geom.periodicShift(domain, mfBA[j], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box sbox = mfBA[j] + pshifts[iiv];

                    if (sbox.intersects(fsi().box()))
                    {
                        sbox &= fsi().box();
                        sbox -= pshifts[iiv];

                        SRRec sr(face(), fsi.index(), pshifts[iiv]);

                        sr.m_fbid = mfcd.AddBox(mfid,
                                                sbox,
                                                0,
                                                j,
                                                0,
                                                0,
                                                mf.nComp());

                        BL_ASSERT(sr.m_fbid.box() == sbox);

                        srrec.push_back(sr);
                    }
                }
            }
        }
    }

    mfcd.CollectData();

    for (int i = 0; i < srrec.size(); i++)
    {
        FabSet& fabset = bndry[srrec[i].m_face];

        BL_ASSERT(fabset.DistributionMap()[srrec[i].m_idx] == MyProc);

        tmpfab.resize(srrec[i].m_fbid.box(), mf.nComp());

        mfcd.FillFab(mfid, srrec[i].m_fbid, tmpfab);

        tmpfab.shift(srrec[i].m_shift);

        fabset[srrec[i].m_idx].plus(tmpfab, 0, 0, mf.nComp());
    }
}

void
SyncRegister::FineDVAdd (const MultiFab& U, 
                         const Real*     dx_fine, 
                         const Geometry& crse_geom, 
                         int             is_rz,
                         int**           fine_bc,
                         Real            mult)
{
    WriteNorms(*this,"SyncRegister::FineDVAdd(B)");

    FArrayBox ufab, ffablo, ffabhi;

    MultiFab cloMF[BL_SPACEDIM], chiMF[BL_SPACEDIM];

    BuildMFs(U,cloMF,chiMF,ratio,WITH_SURROUNDING_BOX);
    //
    // Note that cloMF[dir] and chiMF[dir] have same distribution as U.
    //
    for (ConstMultiFabIterator mfi(U); mfi.isValid(); ++mfi)
    {
        ufab.resize(::grow(mfi.validbox(),1),BL_SPACEDIM);
        ufab.setComplement(0,mfi.validbox(),0,BL_SPACEDIM);
        ufab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,BL_SPACEDIM);

        const int* ulo = ufab.box().loVect();
        const int* uhi = ufab.box().hiVect();
        int* bc        = fine_bc[mfi.index()];

        FillExtDir(ufab, mfi(), mfi.validbox(), bc);

        Box ndbox = ::surroundingNodes(mfi.validbox());

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndbox.smallEnd(dir),1);
            tboxhi.setRange(dir,ndbox.bigEnd(dir),1);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,dir,ratio);

            FArrayBox& cfablo = cloMF[dir][mfi.index()];
            FArrayBox& cfabhi = chiMF[dir][mfi.index()];
            const Box& cboxlo = cfablo.box();
            const Box& cboxhi = cfabhi.box();
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
        }
    }
    //
    // Intersect and add to registers.
    //
    Box crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF[dir],0,0,0,1);
            bndry[face()].plusFrom(chiMF[dir],0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF[dir]);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF[dir]);
    }

    WriteNorms(*this,"SyncRegister::FineDVAdd(E)");
}

inline
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
                           Real            mult)
{
    WriteNorms(*this,"SyncRegister::CrseDsdtAdd(B)");

    Array<IntVect> pshifts(27);

    MultiFab dsdt_local(dsdt.boxArray(),1,1);
    //
    // This MultiFab will parallel dsdt_local.
    //
    BoxArray divuBA = dsdt_local.boxArray();

    divuBA.surroundingNodes();

    MultiFab divuMF(divuBA,1,0);
    //
    // Fill valid region of dsdt_local from dsdt and zero out ghost cells.
    //
    for (MultiFabIterator mfi(dsdt_local); mfi.isValid(); ++mfi)
    {
        mfi().setComplement(0,mfi.validbox(),0,1);
        mfi().copy(dsdt[mfi.index()], mfi.validbox());
    }
    geom.FillPeriodicBoundary(dsdt_local,true,false);

    int imax = geom.Domain().bigEnd()[0]+1;

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
        FArrayBox& divu   = divuMF[mfi.index()];
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

        int nghost         = 0;
        Real hx            = geom.CellSize()[0];
        FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                   mfi().dataPtr(),
                   rcen.dataPtr(), 
                   ARLIM(ndlo), ARLIM(ndhi), divu.dataPtr(), 
                   domlo, domhi, &hx, &is_rz, &imax);
        divu.negate();
        divu.mult(mult);
    }
 

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(divuMF,0,0,0,1);
    }

    WriteNorms(*this,"SyncRegister::CrseDsdtAdd(E)");
}

void
SyncRegister::FineDsdtAdd (const MultiFab& dsdt,
                           const Geometry& geom,
                           const Geometry& crse_geom,
                           int             is_rz,
                           int**           fine_bc, 
                           Real            mult)
{
    WriteNorms(*this,"SyncRegister::FineDsdtAdd(B)");

    FArrayBox dsdtfab, ffablo, ffabhi, ffablo_tmp, ffabhi_tmp;

    MultiFab cloMF[BL_SPACEDIM], chiMF[BL_SPACEDIM];

    BuildMFs(dsdt,cloMF,chiMF,ratio,WITH_SURROUNDING_BOX);

    int imax = geom.Domain().bigEnd()[0]+1;

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
        Box ndbox = ::surroundingNodes(mfi.validbox());

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndbox.smallEnd(dir),1);
            tboxhi.setRange(dir,ndbox.bigEnd(dir),1);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,dir,ratio);

            FArrayBox& cfablo = cloMF[dir][mfi.index()];
            FArrayBox& cfabhi = chiMF[dir][mfi.index()];
            const Box& cboxlo = cfablo.box();
            const Box& cboxhi = cfabhi.box();
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
            int nghost         = 0;
            Real hx            = geom.CellSize()[0];
            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reglo.loVect()), ARLIM(reglo.hiVect()), 
                       ffablo_tmp.dataPtr(),
                       domlo, domhi, &hx, &is_rz, &imax);
            ffablo_tmp.negate();
            ffablo_tmp.mult(mult);
            ffablo.copy(ffablo_tmp);
            ffabhi_tmp.resize(reghi,1);

            FORT_HGC2N(&nghost, ARLIM(dsdtlo), ARLIM(dsdthi), 
                       dsdtfab.dataPtr(),
                       rcen.dataPtr(), 
                       ARLIM(reghi.loVect()), ARLIM(reghi.hiVect()), 
                       ffabhi_tmp.dataPtr(), 
                       domlo, domhi, &hx, &is_rz, &imax);
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
        }
    }
    //
    // Intersect and add to registers.
    //
    Box crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF[dir],0,0,0,1);
            bndry[face()].plusFrom(chiMF[dir],0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF[dir]);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF[dir]);
    }

    WriteNorms(*this,"SyncRegister::FineDsdtAdd(E)");
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
    WriteNorms(*this,"SyncRegister::CompDVAdd(B)");

    Array<IntVect> pshifts(27);

    FArrayBox ufab, ffablo, ffabhi;

    MultiFab cloMF[BL_SPACEDIM], chiMF[BL_SPACEDIM];

    BuildMFs(U,cloMF,chiMF,ratio,WITH_SURROUNDING_BOX);
    //
    // Note that cloMF[dir] and chiMF[dir] have same distribution as U.
    //
    for (ConstMultiFabIterator mfi(U); mfi.isValid(); ++mfi)
    {
        ufab.resize(::grow(mfi.validbox(),1),BL_SPACEDIM);
        ufab.setComplement(0,mfi.validbox(),0,BL_SPACEDIM);
        ufab.copy(mfi(),mfi.validbox(),0,mfi.validbox(),0,BL_SPACEDIM);

        const int* ulo = ufab.box().loVect();
        const int* uhi = ufab.box().hiVect();
        int* bc        = fine_bc[mfi.index()];

        FillExtDir(ufab, mfi(), mfi.validbox(), bc);

        Box ndbox = ::surroundingNodes(mfi.validbox());

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            Box tboxlo(ndbox), tboxhi(ndbox);
            tboxlo.setRange(dir,ndbox.smallEnd(dir),1);
            tboxhi.setRange(dir,ndbox.bigEnd(dir),1);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,dir,ratio);

            FArrayBox& cfablo = cloMF[dir][mfi.index()];
            FArrayBox& cfabhi = chiMF[dir][mfi.index()];
            const Box& cboxlo = cfablo.box();
            const Box& cboxhi = cfabhi.box();
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
        }
    }
    //
    // Intersect and add to registers.
    //
    Box crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF[dir],0,0,0,1);
            bndry[face()].plusFrom(chiMF[dir],0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF[dir]);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF[dir]);
    }

    WriteNorms(*this,"SyncRegister::CompDVAdd(E)");
}

void
SyncRegister::CrseLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Geometry& geom,
                           int             is_rz,
                           Real            mult)
{
    WriteNorms(*this,"SyncRegister::CrseLPhiAdd(B)");
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    BL_ASSERT(Phi.boxArray().length() == Sigma.boxArray().length());

    MultiFab Sig_local(Sigma.boxArray(),1,1);
    MultiFab Phi_local(Phi.boxArray(),1,1);
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
    // divgpMF and Phi_local will have same distribution mapping.
    //
    MultiFab divgpMF(Phi_local.boxArray(),1,0);

    for (MultiFabIterator pmfi(Phi_local); pmfi.isValid(); ++pmfi)
    {
        DependentMultiFabIterator smfi(pmfi, Sig_local);

        const int* slo   = smfi().loVect();
        const int* shi   = smfi().hiVect();
        const int* p_lo  = pmfi().loVect();
        const int* p_hi  = pmfi().hiVect();
        FArrayBox& divgp = divgpMF[pmfi.index()];
        const int* glo   = divgp.loVect();
        const int* ghi   = divgp.hiVect();

        FORT_SRDGPHI(pmfi().dataPtr(),ARLIM(p_lo),ARLIM(p_hi),
                     smfi().dataPtr(),ARLIM(slo),ARLIM(shi),
                     divgp.dataPtr(),ARLIM(glo),ARLIM(ghi),
                     glo,ghi,geom.CellSize(),&mult,&is_rz);
    }

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(divgpMF,0,0,0,1);
    }

    WriteNorms(*this,"SyncRegister::CrseLPhiAdd(E)");
}

void
SyncRegister::FineLPhiAdd (const MultiFab& Phi,
                           const MultiFab& Sigma,
                           const Real*     dx_fine,
                           const Geometry& crse_geom, 
                           int             is_rz,
                           Real            mult)
{
    WriteNorms(*this,"SyncRegister::FineLPhiAdd(B)");
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    BL_ASSERT(Phi.boxArray().length() == Sigma.boxArray().length());

    FArrayBox pfab, sfab, ffablo, ffabhi;

    MultiFab cloMF[BL_SPACEDIM], chiMF[BL_SPACEDIM];

    BuildMFs(Phi,cloMF,chiMF,ratio,WITH_BOX);
    //
    // Note that cloMF[dir] and chiMF[dir] have same distribution as Phi.
    //
    for (ConstMultiFabIterator pmfi(Phi); pmfi.isValid(); ++pmfi)
    {
        ConstDependentMultiFabIterator smfi(pmfi, Sigma);

        pfab.resize(::grow(pmfi.validbox(),1),1);
        sfab.resize(::grow(smfi.validbox(),1),1);
        pfab.setComplement(0,pmfi.validbox(),0,1);
        sfab.setComplement(0,smfi.validbox(),0,1);
        pfab.copy(pmfi(),pmfi.validbox());
        sfab.copy(smfi(),smfi.validbox());

        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();
        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();
        const int* ndlo    = pmfi.validbox().loVect();
        const int* ndhi    = pmfi.validbox().hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            Box tboxlo(pmfi.validbox()), tboxhi(pmfi.validbox());
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,dir,ratio);

            FArrayBox& cfablo = cloMF[dir][pmfi.index()];
            FArrayBox& cfabhi = chiMF[dir][pmfi.index()];
            const Box& cboxlo = cfablo.box();
            const Box& cboxhi = cfabhi.box();
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
        }
    }
    //
    // Intersect and add to registers.
    //
    Box crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF[dir],0,0,0,1);
            bndry[face()].plusFrom(chiMF[dir],0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF[dir]);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF[dir]);
    }

    WriteNorms(*this,"SyncRegister::FineLPhiAdd(E)");
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
    WriteNorms(*this,"SyncRegister::CompLPhiAdd(B)");
    //
    // This code assumes Sigma and Phi have same processor distribution.
    //
    BL_ASSERT(Phi.boxArray().length() == Sigma.boxArray().length());

    Array<IntVect> pshifts(27);

    FArrayBox pfab, sfab, ffablo, ffabhi;

    MultiFab cloMF[BL_SPACEDIM], chiMF[BL_SPACEDIM];

    BuildMFs(Phi,cloMF,chiMF,ratio,WITH_BOX);

    for (ConstMultiFabIterator pmfi(Phi); pmfi.isValid(); ++pmfi)
    {
        ConstDependentMultiFabIterator smfi(pmfi, Sigma);

        pfab.resize(::grow(pmfi.validbox(),1),1);
        sfab.resize(::grow(smfi.validbox(),1),1);
        pfab.setComplement(0,pmfi.validbox(),0,1);
        sfab.setComplement(0,smfi.validbox(),0,1);
        pfab.copy(pmfi(),pmfi.validbox());
        sfab.copy(smfi(),smfi.validbox());

        const int* pfab_lo = pfab.loVect();
        const int* pfab_hi = pfab.hiVect();
        const int* sfab_lo = sfab.loVect();
        const int* sfab_hi = sfab.hiVect();
        const int* ndlo    = pmfi.validbox().loVect();
        const int* ndhi    = pmfi.validbox().hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            Box tboxlo(pmfi.validbox()), tboxhi(pmfi.validbox());
            tboxlo.setRange(dir,ndlo[dir],1);
            tboxhi.setRange(dir,ndhi[dir],1);
            Box reglo(tboxlo), reghi(tboxhi);

            SizeTmpFabs(reglo,reghi,tboxlo,tboxhi,ffablo,ffabhi,dir,ratio);

            FArrayBox& cfablo = cloMF[dir][pmfi.index()];
            FArrayBox& cfabhi = chiMF[dir][pmfi.index()];
            const Box& cboxlo = cfablo.box();
            const Box& cboxhi = cfabhi.box();
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
        }
    }
    //
    // Intersect and add to registers.
    //
    Box crse_node_domain = ::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF[dir],0,0,0,1);
            bndry[face()].plusFrom(chiMF[dir],0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF[dir]);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF[dir]);
    }

    WriteNorms(*this,"SyncRegister::CompLPhiAdd(E)");
}
