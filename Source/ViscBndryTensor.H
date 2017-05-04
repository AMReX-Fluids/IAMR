
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

    ViscBndryTensor (const amrex::BoxArray& grids,
                     const amrex::DistributionMapping& dmap,
                     const amrex::Geometry& geom)
        :
        amrex::MCInterpBndryData(grids,dmap,amrex::MCLinOp::bcComponentsNeeded(),geom) {}

    virtual void setBndryConds (const amrex::BCRec& phys_bc,
                                int          ratio,
                                int          comp=0) override;

    void setHomogValues (const amrex::Array<amrex::BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_VISCBNDRYTENSOR_H_*/
