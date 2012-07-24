
#include <winstd.H>

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <deque>

#include <BC_TYPES.H>
#include <FluxRegister.H>
#include <SyncRegister.H>
#include <NAVIERSTOKES_F.H>
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
        BL_ASSERT(ratio[dir] == -1);

    BL_ASSERT(grids.size() == 0);
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
        BoxArray loBA(grids.size());
        BoxArray hiBA(grids.size());

        for (int k = 0, N = grids.size(); k < N; k++)
        {
            const Box ndbox = BoxLib::surroundingNodes(grids[k]);

            Box nd_lo = ndbox;
            Box nd_hi = ndbox;

            nd_lo.setRange(dir,ndbox.smallEnd(dir),1);
            nd_hi.setRange(dir,ndbox.bigEnd(dir),1);

            loBA.set(k,nd_lo);
            hiBA.set(k,nd_hi);
        }
        //
        // Define the FabSets.
        //
        const Orientation lo = Orientation(dir,Orientation::low);
        const Orientation hi = Orientation(dir,Orientation::high);

        bndry[lo].define(loBA,1);
        bndry_mask[lo].define(loBA,1);
        bndry[hi].define(hiBA,1);
        bndry_mask[hi].define(hiBA,1);
    }
}

SyncRegister::~SyncRegister () {}

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

    Array<IntVect>                pshifts(27);
    std::deque<FluxRegister::Rec> srrec;
    FabSetCopyDescriptor          fscd;
    FabSetId                      fsid[2*BL_SPACEDIM];

    for (OrientationIter face_it; face_it; ++face_it)
    {
        const Orientation face = face_it();

        fsid[face] = fscd.RegisterFabSet((FabSet*) &bndry[face]);

        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
            const int        index  = mfi.index();
            const FArrayBox& rhsfab = rhs[mfi];
            const int        ncomp  = rhsfab.nComp();

            for (int j = 0, N = grids.size(); j < N; j++)
            {
                const Box fabbox = bndry[face].fabbox(j);

                geom.periodicShift(domain,fabbox,pshifts);

                for (int i = 0, M = pshifts.size(); i < M; i++)
                {
                    Box sbox = fabbox + pshifts[i];

                    sbox &= rhsfab.box();

                    if (sbox.ok())
                    {
                        sbox -= pshifts[i];

                        FluxRegister::Rec sr(face, index, pshifts[i]);

                        sr.m_fbid = fscd.AddBox(fsid[face],
                                                sbox,
                                                0,
                                                j,
                                                0,
                                                0,
                                                ncomp);

                        BL_ASSERT(sr.m_fbid.box() == sbox);

                        srrec.push_back(sr);
                    }
                }
            }
        }
    }

    int nrecv = srrec.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no work to do.
        //
        return;

    fscd.CollectData();

    for (std::deque<FluxRegister::Rec>::const_iterator it = srrec.begin(),
             End = srrec.end();
         it != End;
         ++it)
    {
        BL_ASSERT(rhs.DistributionMap()[it->m_idx] == ParallelDescriptor::MyProc());

        FArrayBox& fab = rhs[it->m_idx];
        const Box  dbx = it->m_fbid.box() + it->m_shift;

        fscd.FillFab(fsid[it->m_face], it->m_fbid, fab, dbx);
    }
}

void
SyncRegister::multByBndryMask (MultiFab& rhs) const
{
    FabSetCopyDescriptor          fscd;
    FabSetId                      fsid[2*BL_SPACEDIM];
    std::deque<FluxRegister::Rec> srrec;
    FArrayBox                     tmpfab;

    std::vector< std::pair<int,Box> > isects;

    for (OrientationIter face_it; face_it; ++face_it)
    {
        const Orientation face = face_it();

        fsid[face] = fscd.RegisterFabSet((FabSet*) &bndry_mask[face]);

        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
            const int  index  = mfi.index();
            FArrayBox& rhsfab = rhs[mfi];
            const int  ncomp  = rhsfab.nComp();

            bndry_mask[face].boxArray().intersections(rhsfab.box(),isects);

            for (int i = 0, N = isects.size(); i < N; i++)
            {
                FluxRegister::Rec sr(face,index);

                sr.m_fbid = fscd.AddBox(fsid[face],
                                        isects[i].second,
                                        0,
                                        isects[i].first,
                                        0,
                                        0,
                                        ncomp);
                srrec.push_back(sr);
            }
        }
    }

    int nrecv = srrec.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no work to do.
        //
        return;

    fscd.CollectData();

    for (std::deque<FluxRegister::Rec>::const_iterator it = srrec.begin(),
             End = srrec.end();
         it != End;
         ++it)
    {
        const Box& bx = it->m_fbid.box();

        BL_ASSERT(rhs.DistributionMap()[it->m_idx] == ParallelDescriptor::MyProc());

        FArrayBox& rhs_fab = rhs[it->m_idx];

        tmpfab.resize(bx, rhs_fab.nComp());

        fscd.FillFab(fsid[it->m_face], it->m_fbid, tmpfab, bx);

        rhs_fab.mult(tmpfab, bx, bx, 0, 0, rhs_fab.nComp());
    }
}

void
SyncRegister::InitRHS (MultiFab&       rhs,
                       const Geometry& geom,
                       const BCRec*    phys_bc)
{
    rhs.setVal(0);

    const Box domain = BoxLib::surroundingNodes(geom.Domain());

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

    const Box node_domain = BoxLib::surroundingNodes(geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            Box domlo(node_domain), domhi(node_domain);

            domlo.setRange(dir,node_domain.smallEnd(dir),1);
            domhi.setRange(dir,node_domain.bigEnd(dir),1);

            for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
            {
                const Box blo = mfi.validbox() & domlo;
                const Box bhi = mfi.validbox() & domhi;

                if (blo.ok() && phys_lo[dir] == Outflow)
                    rhs[mfi].mult(0.0,blo,0,1);

                if (bhi.ok() && phys_hi[dir] == Outflow) 
                    rhs[mfi].mult(0.0,bhi,0,1);
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

    std::vector< std::pair<int,Box> > isects;

    for (OrientationIter face_it; face_it; ++face_it)
    {
        FabSet& fs = bndry_mask[face_it()];

        for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
        {
            FArrayBox& fab = fs[fsi];

            Box mask_cells = BoxLib::enclosedCells(BoxLib::grow(fab.box(),1));

            tmpfab.resize(mask_cells,1);
            tmpfab.setVal(0);

            grids.intersections(mask_cells,isects);

            for (int i = 0, N = isects.size(); i < N; i++)
            {
                tmpfab.setVal(1,isects[i].second,0,1);
            }
 
            if (geom.isAnyPeriodic())
            {
                geom.periodicShift(geom.Domain(),mask_cells,pshifts);

                for (int j = 0, M = pshifts.size(); j < M; j++)
                {
                    mask_cells += pshifts[j];

                    grids.intersections(mask_cells,isects);

                    for (int i = 0, N = isects.size(); i < N; i++)
                    {
                        Box& isect = isects[i].second;
                        isect     -= pshifts[j];
                        tmpfab.setVal(1,isect,0,1);
                    }

                    mask_cells -= pshifts[j];
                }
            }
            Real* mask_dat = fab.dataPtr();
            const int* mlo = fab.loVect(); 
            const int* mhi = fab.hiVect();
            Real* cell_dat = tmpfab.dataPtr();
            const int* clo = tmpfab.loVect(); 
            const int* chi = tmpfab.hiVect();
        
            FORT_MAKEMASK(mask_dat,ARLIM(mlo),ARLIM(mhi), cell_dat,ARLIM(clo),ARLIM(chi));
        }
    }
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

            for (OrientationIter face_it; face_it; ++face_it)
            {
                FabSet& fs = bndry_mask[face_it()];

                for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
                {
                    FArrayBox& fab = fs[fsi];

                    Box isect = fab.box() & domlo;

                    if (isect.ok())
                        fab.mult(2.0,isect,0,1);

                    isect = fab.box() & domhi;

                    if (isect.ok())
                        fab.mult(2.0,isect,0,1);
                }
            }
        }
    }
    //
    // Here convert from sum of cell contributions to 0 or 1.
    //
    for (OrientationIter face_it; face_it; ++face_it)
    {
        FabSet& fs = bndry_mask[face_it()];

        for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
        {
            FArrayBox& fab      = fs[fsi];
            Real*      mask_dat = fab.dataPtr();
            const int* mlo      = fab.loVect(); 
            const int* mhi      = fab.hiVect();

            FORT_CONVERTMASK(mask_dat,ARLIM(mlo),ARLIM(mhi));
        }
    }

    multByBndryMask(rhs);
}

static
void
BuildMFs (const MultiFab& mf,
          MultiFab&       cloMF,
          MultiFab&       chiMF,
          const IntVect&  ratio,
          int             dir)
{
    BoxArray cloBA(mf.size());
    BoxArray chiBA(mf.size());

    for (int i = 0, N = mf.size(); i < N; i++)
    {
        const Box& bx = mf.boxArray()[i];

        Box tboxlo(bx), tboxhi(bx);

        tboxlo.setRange(dir,bx.smallEnd(dir),1);
        tboxhi.setRange(dir,bx.bigEnd(dir),1);

        cloBA.set(i,BoxLib::coarsen(tboxlo,ratio));
        chiBA.set(i,BoxLib::coarsen(tboxhi,ratio));
    }

    cloMF.define(cloBA,1,0,Fab_allocate);
    chiMF.define(chiBA,1,0,Fab_allocate);

    cloMF.setVal(0);
    chiMF.setVal(0);
}

void
SyncRegister::incrementPeriodic (const Geometry& geom,
                                 const Box&      domain,
                                 const MultiFab& mf)
{
    if (!geom.isAnyPeriodic()) return;

    Array<IntVect>                pshifts(27);
    std::deque<FluxRegister::Rec> srrec;
    FArrayBox                     tmpfab;
    MultiFabCopyDescriptor        mfcd;
    const BoxArray&               mfba = mf.boxArray();
    MultiFabId                    mfid = mfcd.RegisterFabArray((MultiFab*) &mf);

    for (int j = 0, N = mfba.size(); j < N; j++)
    {
        const Box& bx = mfba[j];

        geom.periodicShift(domain, bx, pshifts);

        if (pshifts.empty()) continue;

        for (OrientationIter face_it; face_it; ++face_it)
        {
            const Orientation face = face_it();

            FabSet& fs = bndry[face];

            for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
            {
                const int        idx   = fsi.index();
                const int        ncomp = mf.nComp();
                const FArrayBox& fab   = fs[fsi];

                for (int i = 0, M = pshifts.size(); i < M; i++)
                {
                    Box sbox = bx + pshifts[i];

                    sbox &= fab.box();

                    if (sbox.ok())
                    {
                        sbox -= pshifts[i];

                        FluxRegister::Rec sr(face, idx, pshifts[i]);

                        sr.m_fbid = mfcd.AddBox(mfid,
                                                sbox,
                                                0,
                                                j,
                                                0,
                                                0,
                                                ncomp);

                        BL_ASSERT(sr.m_fbid.box() == sbox);

                        srrec.push_back(sr);
                    }
                }
            }
        }
    }

    int nrecv = srrec.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no work to do.
        //
        return;

    mfcd.CollectData();

    for (std::deque<FluxRegister::Rec>::const_iterator it = srrec.begin(),
             End = srrec.end();
         it != End;
         ++it)
    {
        FabSet& fabset = bndry[it->m_face];

        BL_ASSERT(fabset.DistributionMap()[it->m_idx] == ParallelDescriptor::MyProc());

        tmpfab.resize(it->m_fbid.box(), mf.nComp());

        mfcd.FillFab(mfid, it->m_fbid, tmpfab);

        tmpfab.shift(it->m_shift);

        fabset[it->m_idx].plus(tmpfab, 0, 0, mf.nComp());
    }
}

void
SyncRegister::CrseInit  (MultiFab*       Sync_resid_crse,
                         const Geometry& crse_geom, 
                         Real            mult)
{
    setVal(0);

    Sync_resid_crse->mult(mult);

    const Box crse_node_domain = BoxLib::surroundingNodes(crse_geom.Domain());

    incrementPeriodic(crse_geom, crse_node_domain, *Sync_resid_crse);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(*Sync_resid_crse,0,0,0,1);
    }
}

void
SyncRegister::CompAdd  (MultiFab*       Sync_resid_fine,
                        const Geometry& fine_geom, 
                        const Geometry& crse_geom, 
                        const BCRec*    phys_bc,
                        const BoxArray& Pgrids,
                        Real            mult)
{
    Array<IntVect> pshifts(27);

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(*Sync_resid_fine); mfi.isValid(); ++mfi)
    {
        const Box& sync_box = mfi.validbox();

        Pgrids.intersections(sync_box,isects);

        FArrayBox& syncfab = (*Sync_resid_fine)[mfi];

        for (int ii = 0, N = isects.size(); ii < N; ii++)
        {
            const int  i   = isects[ii].first;
            const Box& pbx = Pgrids[i];

            syncfab.setVal(0,isects[ii].second,0,1);

            fine_geom.periodicShift(sync_box, pbx, pshifts);

            for (int j = 0, M = pshifts.size(); j < M; j++)
            {
                Box isect = pbx + pshifts[j];
                isect    &= sync_box;
                syncfab.setVal(0,isect,0,1);
            }
        }
    }

    FineAdd(Sync_resid_fine,fine_geom,crse_geom,phys_bc,mult);
}

void
SyncRegister::FineAdd  (MultiFab* Sync_resid_fine,
                        const Geometry& fine_geom, 
                        const Geometry& crse_geom, 
                        const BCRec*    phys_bc,
                        Real            mult)
{
    Sync_resid_fine->mult(mult);

    const Box crse_node_domain = BoxLib::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        MultiFab cloMF, chiMF;

        BuildMFs(*Sync_resid_fine,cloMF,chiMF,ratio,dir);
        //
        // Coarsen edge values.
        //
        const int N = Sync_resid_fine->IndexMap().size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
        {
            const int k = Sync_resid_fine->IndexMap()[i];

            FArrayBox& finefab = (*Sync_resid_fine)[k];

            const int* resid_lo = finefab.box().loVect();
            const int* resid_hi = finefab.box().hiVect();

            FArrayBox& cfablo = cloMF[k];
            const Box& cboxlo = cfablo.box();
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();

            FORT_SRCRSEREG(finefab.dataPtr(),
                           ARLIM(resid_lo),ARLIM(resid_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());

            FArrayBox& cfabhi = chiMF[k];
            const Box& cboxhi = cfabhi.box();
            clo = cboxhi.loVect();
            chi = cboxhi.hiVect();

            FORT_SRCRSEREG(finefab.dataPtr(),
                           ARLIM(resid_lo),ARLIM(resid_hi),
                           cfabhi.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());
        }

        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            if (!crse_geom.isPeriodic(j))
            {
                //
                // Now points on the physical bndry must be doubled
                // for any boundary but outflow or periodic
                //
                Box domlo(crse_node_domain), domhi(crse_node_domain);

                domlo.setRange(j,crse_node_domain.smallEnd(j),1);
                domhi.setRange(j,crse_node_domain.bigEnd(j),1);

                for (MFIter mfi(*Sync_resid_fine); mfi.isValid(); ++mfi)
                {
                    FArrayBox& cfablo = cloMF[mfi.index()];
                    const Box& cboxlo = cfablo.box();

                    Box isectlo = domlo & cboxlo;
                    Box isecthi = domhi & cboxlo;

                    if (isectlo.ok())
                        cfablo.mult(2.0,isectlo,0,1);
                    if (isecthi.ok())
                        cfablo.mult(2.0,isecthi,0,1);

                    FArrayBox& cfabhi = chiMF[mfi.index()];
                    const Box& cboxhi = cfabhi.box();

                    isectlo = domlo & cboxhi;
                    isecthi = domhi & cboxhi;

                    if (isectlo.ok())
                        cfabhi.mult(2.0,isectlo,0,1);
                    if (isecthi.ok())
                        cfabhi.mult(2.0,isecthi,0,1);
                }
            }
        }
        //
        // Intersect and add to registers.
        //
        for (OrientationIter face; face; ++face)
        {
            bndry[face()].plusFrom(cloMF,0,0,0,1);
            bndry[face()].plusFrom(chiMF,0,0,0,1);
        }
        incrementPeriodic(crse_geom, crse_node_domain, cloMF);
        incrementPeriodic(crse_geom, crse_node_domain, chiMF);
    }
}

