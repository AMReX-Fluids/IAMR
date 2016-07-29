
#ifndef _SYNCREGISTER_H_
#define _SYNCREGISTER_H_

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

    SyncRegister (const BoxArray& fine_boxes, const IntVect& ref_ratio);

    ~SyncRegister ();

    void CrseInit (MultiFab& Sync_resid_crse, const Geometry& crse_geom, Real mult);

    void FineAdd (MultiFab& Sync_resid_fine, const Geometry& crse_geom, Real mult);

    void CompAdd (MultiFab& Sync_resid_fine, const Geometry& fine_geom, const Geometry& crse_geom,
		  const BoxArray& Pgrids, Real mult);

    void InitRHS (MultiFab& rhs, const Geometry& geom, const BCRec& phys_bc);

protected:
    //
    // Helpers for various parallel routines.
    //
    typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

    typedef FabArrayBase::CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

    enum Who { CopyPeriodic, IncrementPeriodic };

    void SendRecvDoit (const MapOfCopyComTagContainers& m_SndTags,
                       const MapOfCopyComTagContainers& m_RcvTags,
                       const std::map<int,int>&         m_SndVols,
                       const std::map<int,int>&         m_RcvVols,
                       int                              ncomp,
                       Who                              who,
                       MultiFab*                        rhs,
                       const MultiFab*                  mf);
    //
    // Helper member functions.
    //
    void incrementPeriodic (const Geometry& geom,
                            const Box&      domain,
                            const MultiFab& mf);

    void copyPeriodic (const Geometry& geom,
                       const Box&      domain,
                       MultiFab&       rhs);

    void multByBndryMask (MultiFab& rhs) const;
    //
    // The data.
    //
    FabSet  bndry_mask[2*BL_SPACEDIM];
    IntVect ratio;
};

#endif /*_SYNCREGISTER_H_*/
