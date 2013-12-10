
#include <winstd.H>

#include <cstdlib>
#include <cstring>
#include <cstdio>

#include <BC_TYPES.H>
#include <FluxRegister.H>
#include <SyncRegister.H>
#include <NAVIERSTOKES_F.H>
#include <SYNCREG_F.H>
#include <Profiler.H>

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
        const int N = grids.size();

        BoxArray loBA(N), hiBA(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < N; k++)
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
// Consolidate the various parallel operations in this file.
//

void
SyncRegister::SendRecvDoit (const MapOfCopyComTagContainers& m_SndTags,
                            const MapOfCopyComTagContainers& m_RcvTags,
                            const std::map<int,int>&         m_SndVols,
                            const std::map<int,int>&         m_RcvVols,
                            int                              ncomp,
                            Who                              who,
                            MultiFab*                        rhs,
                            const MultiFab*                  mf)
{
    BL_ASSERT(who == SyncRegister::CopyPeriodic || who ==  SyncRegister::IncrementPeriodic);

#ifdef BL_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (m_SndTags.empty() && m_RcvTags.empty())
        //
        // No parallel work for this MPI process to do.
        //
        return;

    FArrayBox          fab;
    Array<MPI_Status>  stats;
    Array<int>         recv_from, index;
    Array<double*>     recv_data, send_data;
    Array<MPI_Request> recv_reqs, send_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    double* the_recv_data = 0;

    FabArrayBase::PostRcvs(m_RcvTags,m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);
    //
    // Send the data.
    //
    for (MapOfCopyComTagContainers::const_iterator m_it = m_SndTags.begin(),
             m_End = m_SndTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = m_SndVols.find(m_it->first);

        BL_ASSERT(vol_it != m_SndVols.end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        double* data = static_cast<double*>(BoxLib::The_Arena()->alloc(N*sizeof(double)));
        double* dptr = data;

        for (CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                 End = m_it->second.end();
             it != End;
             ++it)
        {
            const Box& bx = it->box;
            fab.resize(bx,ncomp);

            if (who == SyncRegister::CopyPeriodic)
            {
                fab.copy(bndry[it->srcIndex][it->fabIndex],bx,0,bx,0,ncomp);
            }
            else
            {
                fab.copy((*mf)[it->fabIndex],bx,0,bx,0,ncomp);
            }

            const int Cnt = bx.numPts()*ncomp;
            memcpy(dptr,fab.dataPtr(),Cnt*sizeof(double));
            dptr += Cnt;
        }
        BL_ASSERT(data+N == dptr);

        if (FabArrayBase::do_async_sends)
        {
            send_data.push_back(data);
            send_reqs.push_back(ParallelDescriptor::Asend(data,N,m_it->first,SeqNum).req());
        }
        else
        {
            ParallelDescriptor::Send(data,N,m_it->first,SeqNum);
            BoxLib::The_Arena()->free(data);
        }
    }
    //
    // Now receive and unpack FAB data as it becomes available.
    //
    const int N_rcvs = m_RcvTags.size();

    index.resize(N_rcvs);
    stats.resize(N_rcvs);

    for (int NWaits = N_rcvs, completed; NWaits > 0; NWaits -= completed)
    {
        ParallelDescriptor::Waitsome(recv_reqs, completed, index, stats);

        for (int k = 0; k < completed; k++)
        {
            const double* dptr = recv_data[index[k]];

            BL_ASSERT(dptr != 0);

            MapOfCopyComTagContainers::const_iterator m_it = m_RcvTags.find(recv_from[index[k]]);

            BL_ASSERT(m_it != m_RcvTags.end());

            for (CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->box;
                fab.resize(bx,ncomp);
                const int Cnt = bx.numPts()*ncomp;
                memcpy(fab.dataPtr(),dptr,Cnt*sizeof(double));

                if (who == SyncRegister::CopyPeriodic)
                {
                    (*rhs)[it->fabIndex].copy(fab,bx,0,bx,0,ncomp);
                }
                else
                {
                    bndry[it->srcIndex][it->fabIndex].plus(fab,bx,bx,0,0,ncomp);
                }

                dptr += Cnt;
            }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !m_SndTags.empty())
        FabArrayBase::GrokAsyncSends(m_SndTags.size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

//
// If periodic, copy the values from sync registers onto the nodes of the
// rhs which are not covered by sync registers through periodic shifts.
//

void
SyncRegister::copyPeriodic (const Geometry& geom,
                            const Box&      domain,
                            MultiFab&       rhs)
{
    if (!geom.isAnyPeriodic()) return;

    BL_PROFILE("SyncRegister::copyPeriodic()");

    MapOfCopyComTagContainers  m_SndTags, m_RcvTags;
    std::map<int,int>          m_SndVols, m_RcvVols;
    Array<IntVect>             pshifts(27);
    const int                  MyProc  = ParallelDescriptor::MyProc();
    const int                  ncomp   = rhs.nComp();
    const DistributionMapping& dstDMap = rhs.DistributionMap();

    for (OrientationIter face_it; face_it; ++face_it)
    {
        const Orientation          face    = face_it();
        const FabSet&              fabset  = bndry[face];
        const DistributionMapping& srcDMap = fabset.DistributionMap();

        for (int i = 0; i < rhs.size(); i++)
        {
            const int dst_owner = dstDMap[i];

            for (int j = 0, N = grids.size(); j < N; j++)
            {
                const int src_owner = srcDMap[j];

                if (dst_owner != MyProc && src_owner != MyProc) continue;

                const Box rhsbox = rhs.fabbox(i);
                const Box fabbox = fabset.fabbox(j);

                if (domain.contains(BoxLib::grow(fabbox,1))) continue;

                geom.periodicShift(domain,fabbox,pshifts);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    Box dbx = fabbox + iv; dbx &= rhsbox;

                    if (!dbx.ok()) continue;

                    const Box sbx = dbx - iv;

                    FabArrayBase::CopyComTag tag;

                    if (dst_owner == MyProc)
                    {
                        if (src_owner == MyProc)
                        {
                            //
                            // Do the local work right here.
                            //
                            rhs[i].copy(fabset[j],sbx,0,dbx,0,ncomp);
                        }
                        else
                        {
                            tag.fabIndex = i;
                            tag.box      = dbx;

                            FabArrayBase::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,dbx);
                        }
                    }
                    else if (src_owner == MyProc)
                    {
                        tag.fabIndex = j;
                        tag.srcIndex = face;  // Store face in srcIndex!
                        tag.box      = sbx;

                        FabArrayBase::SetSendTag(m_SndTags,dst_owner,tag,m_SndVols,sbx);
                    }
                }
            }
        }
    }

    if (ParallelDescriptor::NProcs() == 1) return;

    SendRecvDoit(m_SndTags,m_RcvTags,m_SndVols,m_RcvVols,ncomp,SyncRegister::CopyPeriodic,&rhs,0);
}

void
SyncRegister::multByBndryMask (MultiFab& rhs) const
{
    BL_ASSERT(rhs.nComp() == 1);
 
    for (OrientationIter face; face; ++face)
    {
        BL_ASSERT(bndry_mask[face()].nComp() == 1);

        const int N = rhs.size();

        BoxArray ba(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
            ba.set(i,rhs.fabbox(i));

        FabSet fs(ba, 1);

        fs.setVal(1);

        fs.copyFrom(bndry_mask[face()]);

        for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
            rhs[mfi].mult(fs[mfi], ba[mfi.index()], ba[mfi.index()], 0, 0, 1);
    }
}

void
SyncRegister::InitRHS (MultiFab&       rhs,
                       const Geometry& geom,
                       const BCRec*    phys_bc)
{
    BL_PROFILE("SyncRegister::InitRHS()");

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
 
            if (geom.isAnyPeriodic() && !geom.Domain().contains(mask_cells))
            {
                geom.periodicShift(geom.Domain(),mask_cells,pshifts);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    mask_cells += iv;

                    grids.intersections(mask_cells,isects);

                    for (int i = 0, N = isects.size(); i < N; i++)
                    {
                        Box& isect = isects[i].second;
                        isect     -= iv;
                        tmpfab.setVal(1,isect,0,1);
                    }

                    mask_cells -= iv;
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

    grids.clear_hash_bin();
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
    const int N = mf.size();

    BoxArray cloBA(N), chiBA(N);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
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

    BL_PROFILE("SyncRegister::incrementPeriodic()");

    FArrayBox                   fab;
    FabArrayBase::CopyComTag    tag;
    MapOfCopyComTagContainers   m_SndTags, m_RcvTags;
    std::map<int,int>           m_SndVols, m_RcvVols;
    Array<IntVect>              pshifts(27);

    const int                  MyProc  = ParallelDescriptor::MyProc();
    const int                  ncomp   = mf.nComp();
    const BoxArray&            mfba    = mf.boxArray();
    const DistributionMapping& srcDMap = mf.DistributionMap();

    for (int j = 0, N = mfba.size(); j < N; j++)
    {
        const Box& bx = mfba[j];

        if (domain.contains(BoxLib::grow(bx,1))) continue;

        geom.periodicShift(domain, bx, pshifts);

        if (pshifts.empty()) continue;

        const int src_owner = srcDMap[j];

        for (OrientationIter face_it; face_it; ++face_it)
        {
            const Orientation          face    = face_it();
            FabSet&                    fabset  = bndry[face];
            const DistributionMapping& dstDMap = fabset.DistributionMap();

            for (int i = 0; i < fabset.size(); i++)
            {
                const int dst_owner = dstDMap[i];

                if (dst_owner != MyProc && src_owner != MyProc) continue;

                const Box fabbox = fabset.fabbox(i);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    Box dbx = bx + iv; dbx &= fabbox;

                    if (!dbx.ok()) continue;

                    const Box sbx = dbx - iv;

                    if (dst_owner == MyProc)
                    {
                        if (src_owner == MyProc)
                        {
                            //
                            // Do the local work right here.
                            //
                            fabset[i].plus(mf[j],sbx,dbx,0,0,ncomp);
                        }
                        else
                        {
                            tag.fabIndex = i;
                            tag.srcIndex = face;  // Store face in srcIndex!
                            tag.box      = dbx;

                            FabArrayBase::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,dbx);
                        }
                    }
                    else if (src_owner == MyProc)
                    {
                        tag.fabIndex = j;
                        tag.box      = sbx;

                        FabArrayBase::SetSendTag(m_SndTags,dst_owner,tag,m_SndVols,sbx);
                    }
                }
            }
        }
    }

    if (ParallelDescriptor::NProcs() == 1) return;

    SendRecvDoit(m_SndTags,m_RcvTags,m_SndVols,m_RcvVols,ncomp,SyncRegister::IncrementPeriodic,0,&mf);
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

            for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                 it != End;
                 ++it)
            {
                Box isect = pbx + *it;
                isect    &= sync_box;
                syncfab.setVal(0,isect,0,1);
            }
        }
    }

    Pgrids.clear_hash_bin();

    FineAdd(Sync_resid_fine,fine_geom,crse_geom,phys_bc,mult);
}

void
SyncRegister::FineAdd  (MultiFab* Sync_resid_fine,
                        const Geometry& fine_geom, 
                        const Geometry& crse_geom, 
                        const BCRec*    phys_bc,
                        Real            mult)
{
    BL_PROFILE("SyncRegister::FineAdd()");

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

