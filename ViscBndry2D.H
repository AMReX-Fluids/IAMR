#ifndef _VISCBNDRY2D_H_
#define _VISCBNDRY2D_H_

//
// $Id: ViscBndry2D.H,v 1.6 1999-01-06 00:51:15 marc Exp $
//

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry2D
    :
    public MCInterpBndryData
{
public:
    ViscBndry2D () {}

    ViscBndry2D (const BoxArray& grids,
                 const Geometry& geom)
        :
        MCInterpBndryData(grids,4,geom) {}

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

#endif /*_VISCBNDRY2D_H_*/
