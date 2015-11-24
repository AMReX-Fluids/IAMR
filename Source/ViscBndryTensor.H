
#ifndef _VISCBNDRYTENSOR_H_
#define _VISCBNDRYTENSOR_H_

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>
#include <MCLinOp.H>
#include <BL_CXX11.H>

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
                                int          comp=0) BL_OVERRIDE;

    void setHomogValues (const Array<BCRec>& phys_bc,
                         int                 ratio);

};

#endif /*_VISCBNDRYTENSOR_H_*/
