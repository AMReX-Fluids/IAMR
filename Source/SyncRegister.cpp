
#include <BC_TYPES.H>
#include <FluxRegister.H>
#include <SyncRegister.H>
#include <NAVIERSTOKES_F.H>
#include <SYNCREG_F.H>
#include <BLProfiler.H>

#ifdef _OPENMP
#include <omp.h>
#endif

SyncRegister::SyncRegister (const BoxArray& fine_boxes,
                            const IntVect&  ref_ratio)
    : ratio(ref_ratio)
{
    BL_ASSERT(grids.size() == 0);
    BL_ASSERT(fine_boxes.isDisjoint());

    grids = fine_boxes;
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	Orientation loface(dir, Orientation::low);
	Orientation hiface(dir, Orientation::high);

	BndryBATransformer lotrans(loface, IndexType::TheNodeType(), 0, 1, 0);
	BndryBATransformer hitrans(hiface, IndexType::TheNodeType(), 0, 1, 0);

	BoxArray loBA(grids, lotrans);
	BoxArray hiBA(grids, hitrans);

        bndry[loface].define(loBA,1);
        bndry_mask[loface].define(loBA,1);
        bndry[hiface].define(hiBA,1);
        bndry_mask[hiface].define(hiBA,1);
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

    FabArrayBase::PostRcvs(m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);
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
            const Box& bx = it->dbox;
            fab.resize(bx,ncomp);

            if (who == SyncRegister::CopyPeriodic)
            {
                fab.copy(bndry[it->srcIndex][it->dstIndex],bx,0,bx,0,ncomp);
            }
            else
            {
                fab.copy((*mf)[it->dstIndex],bx,0,bx,0,ncomp);
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
                const Box& bx = it->sbox;
                fab.resize(bx,ncomp);
                const int Cnt = bx.numPts()*ncomp;
                memcpy(fab.dataPtr(),dptr,Cnt*sizeof(double));

                if (who == SyncRegister::CopyPeriodic)
                {
                    (*rhs)[it->dstIndex].copy(fab,bx,0,bx,0,ncomp);
                }
                else
                {
                    bndry[it->srcIndex][it->dstIndex].plus(fab,bx,bx,0,0,ncomp);
                }

                dptr += Cnt;
            }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !m_SndTags.empty())
        FabArrayBase::WaitForAsyncSends(m_SndTags.size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

void
SyncRegister::multByBndryMask (MultiFab& rhs) const
{
    BL_ASSERT(rhs.nComp() == 1);
 
    int ngrow = rhs.nGrow();

    MultiFab tmp(rhs.boxArray(), 1, ngrow, rhs.DistributionMap());

    for (OrientationIter face; face; ++face)
    {
        BL_ASSERT(bndry_mask[face()].nComp() == 1);
	
	tmp.setVal(1.0);

	bndry_mask[face()].copyTo(tmp, ngrow, 0, 0, 1);

	MultiFab::Multiply(rhs, tmp, 0, 0, 1, ngrow);
    }
}

void /* note that rhs is on a different BoxArray */
SyncRegister::InitRHS (MultiFab& rhs, const Geometry& geom, const BCRec& phys_bc)
{
    BL_PROFILE("SyncRegister::InitRHS()");

    rhs.setVal(0);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].copyTo(rhs,0,0,0,bndry[face()].nComp(),geom.periodicity());
    }

    const int* phys_lo = phys_bc.lo();
    const int* phys_hi = phys_bc.hi();

    const Box& node_domain = BoxLib::surroundingNodes(geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!geom.isPeriodic(dir))
        {
            Box domlo(node_domain), domhi(node_domain);

            domlo.setRange(dir,node_domain.smallEnd(dir),1);
            domhi.setRange(dir,node_domain.bigEnd(dir),1);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
            {
                const Box& blo = mfi.validbox() & domlo;
                const Box& bhi = mfi.validbox() & domhi;

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

    for (OrientationIter face_it; face_it; ++face_it)
    {
        FabSet& fs = bndry_mask[face_it()];

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    FArrayBox tmpfab;
	    std::vector< std::pair<int,Box> > isects;	    
	    Array<IntVect> pshifts(26);

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
			
			grids.intersections(mask_cells+iv,isects);
			
			for (int i = 0, N = isects.size(); i < N; i++)
			{
			    Box& isect = isects[i].second;
			    isect     -= iv;
			    tmpfab.setVal(1,isect,0,1);
			}
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

#ifdef _OPENMP
#pragma omp parallel
#endif
                for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
                {
                    FArrayBox& fab = fs[fsi];

                    const Box& blo = fab.box() & domlo;

                    if (blo.ok())
                        fab.mult(2.0,blo,0,1);

                    const Box& bhi = fab.box() & domhi;

                    if (bhi.ok())
                        fab.mult(2.0,bhi,0,1);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
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

                const Box& fabbox = fabset.fabbox(i);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    Box dbx = bx + iv; dbx &= fabbox;

                    if (!dbx.ok()) continue;

                    const Box& sbx = dbx - iv;

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
                            tag.dstIndex = i;
                            tag.srcIndex = face;  // Store face in srcIndex!
                            tag.dbox     = dbx;
                            tag.sbox     = dbx;

                            FabArrayBase::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,dbx);
                        }
                    }
                    else if (src_owner == MyProc)
                    {
                        tag.dstIndex = j;
                        tag.sbox     = sbx;
                        tag.dbox     = sbx;

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
SyncRegister::CrseInit (MultiFab& Sync_resid_crse, const Geometry& crse_geom, Real mult)
{
    BL_PROFILE("SyncRegister::CrseInit()");

    setVal(0);

    Sync_resid_crse.mult(mult);

    for (OrientationIter face; face; ++face)
    {
        bndry[face()].plusFrom(Sync_resid_crse,0,0,0,1,crse_geom.periodicity());
    }
}

void
SyncRegister::CompAdd (MultiFab& Sync_resid_fine, const Geometry& fine_geom, const Geometry& crse_geom, 
		       const BoxArray& Pgrids, Real mult)
{
    BL_PROFILE("SyncRegister::CompAdd()");

    Array<IntVect> pshifts(27);

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(Sync_resid_fine); mfi.isValid(); ++mfi)
    {
        const Box& sync_box = mfi.validbox();

        Pgrids.intersections(sync_box,isects);

        FArrayBox& syncfab = Sync_resid_fine[mfi];

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

    FineAdd(Sync_resid_fine,crse_geom,mult);
}

void
SyncRegister::FineAdd (MultiFab& Sync_resid_fine, const Geometry& crse_geom, Real mult)
{
    BL_PROFILE("SyncRegister::FineAdd()");

    Sync_resid_fine.mult(mult);

    const Box& crse_node_domain = BoxLib::surroundingNodes(crse_geom.Domain());

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	Orientation loface(dir, Orientation::low);
	Orientation hiface(dir, Orientation::high);
	
        MultiFab cloMF(bndry_mask[loface].boxArray(),1,0,bndry_mask[loface].DistributionMap());
        MultiFab chiMF(bndry_mask[hiface].boxArray(),1,0,bndry_mask[hiface].DistributionMap());

	BL_ASSERT(cloMF.DistributionMap() == Sync_resid_fine.DistributionMap());
	BL_ASSERT(chiMF.DistributionMap() == Sync_resid_fine.DistributionMap());

        //
        // Coarsen edge values.
        //
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(Sync_resid_fine); mfi.isValid(); ++mfi)
        {
            FArrayBox& finefab = Sync_resid_fine[mfi];

            const int* resid_lo = finefab.box().loVect();
            const int* resid_hi = finefab.box().hiVect();

            FArrayBox& cfablo = cloMF[mfi];
            const Box& cboxlo = cfablo.box();
            const int* clo = cboxlo.loVect();
            const int* chi = cboxlo.hiVect();

            FORT_SRCRSEREG(finefab.dataPtr(),
                           ARLIM(resid_lo),ARLIM(resid_hi),
                           cfablo.dataPtr(),ARLIM(clo),ARLIM(chi),
                           clo,chi,&dir,ratio.getVect());

            FArrayBox& cfabhi = chiMF[mfi];
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

                for (MFIter mfi(Sync_resid_fine); mfi.isValid(); ++mfi)
                {
                    FArrayBox& cfablo = cloMF[mfi];
                    const Box& cboxlo = cfablo.box();

                    Box isectlo = domlo & cboxlo;
                    Box isecthi = domhi & cboxlo;

                    if (isectlo.ok())
                        cfablo.mult(2.0,isectlo,0,1);
                    if (isecthi.ok())
                        cfablo.mult(2.0,isecthi,0,1);

                    FArrayBox& cfabhi = chiMF[mfi];
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

