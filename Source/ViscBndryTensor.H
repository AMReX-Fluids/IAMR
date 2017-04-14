
#ifndef _VISCBNDRYTENSOR_H_
#define _VISCBNDRYTENSOR_H_

#include <AMReX_MCInterpBndryData.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MCLinOp.H>

class ViscBndryTensor
    :
    public MCInterpBndryData
{
public:
    ViscBndryTensor () {}

    ViscBndryTensor (const BoxArray& grids,
                     const Geometry& geom)
        :
        MCInterpBndryData(grids,MCLinOp::bcComponentsNeeded(),geom) {}

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0) override;

    void setHomogValues (const Array<BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_VISCBNDRYTENSOR_H_*/
