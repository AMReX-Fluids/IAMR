
#ifndef _VISCBNDRYTENSOR_H_
#define _VISCBNDRYTENSOR_H_

#include <AMReX_MCInterpBndryData.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MCLinOp.H>

class ViscBndryTensor
    :
    public amrex::MCInterpBndryData
{
public:
    ViscBndryTensor () {}

    ViscBndryTensor (const amrex::BoxArray& grids_,
                     const amrex::DistributionMapping& dmap_,
                     const amrex::Geometry& geom_)
        :
        amrex::MCInterpBndryData(grids_,dmap_,amrex::MCLinOp::bcComponentsNeeded(),geom_) {}

    virtual void setBndryConds (const amrex::BCRec& phys_bc,
                                int          ratio,
                                int          comp=0) override;

    void setHomogValues (const amrex::Vector<amrex::BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_VISCBNDRYTENSOR_H_*/
