#ifndef _VISCBNDRYTENSOR_H_
#define _VISCBNDRYTENSOR_H_

//
// $Id: ViscBndryTensor.H,v 1.1 1999-02-18 18:22:18 sstanley Exp $
//

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndryTensor
    :
    public MCInterpBndryData
{
public:
    ViscBndryTensor () {}

    ViscBndryTensor (const BoxArray& grids,
                     const Geometry& geom)
        :
#if BL_SPACEDIM == 2
        MCInterpBndryData(grids,4,geom) {}
#elif BL_SPACEDIM == 3
        MCInterpBndryData(grids,12,geom) {}
#endif


#if 0
    virtual void setBndryConds (const Array<BCRec>& phys_bc,
                                const Geometry&     geom,
                                IntVect&            ratio);

    void setHomogValues (const Array<BCRec>& phys_bc,
                         IntVect&            ratio);
#else
    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0);

    void setHomogValues (const Array<BCRec>& phys_bc,
                         int                 ratio);
#endif

};

#endif /*_VISCBNDRYTENSOR_H_*/
