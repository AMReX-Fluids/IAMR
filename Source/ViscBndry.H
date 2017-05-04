
#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

#include <AMReX_InterpBndryData.H>
#include <AMReX_BC_TYPES.H>

class ViscBndry
    :
    public amrex::InterpBndryData
{
public:
    ViscBndry ();

    ViscBndry (const amrex::BoxArray& _grids,
               const amrex::DistributionMapping& _dmap,
               int             _ncomp,
               const amrex::Geometry& _geom);

    virtual ~ViscBndry () override;

    virtual void setBndryConds (const amrex::BCRec& phys_bc,
                                int          ratio) override;

    virtual void setBndryConds (const amrex::BCRec&   phys_bc,
                                const amrex::IntVect& ratio,
				int            comp=0) override;

    void setHomogValues (const amrex::BCRec&   bc,
                         const amrex::IntVect& ratio);
};

#endif /*_VISCBNDRYDATA_H_*/
