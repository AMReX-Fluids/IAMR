#ifndef _VISCBNDRY2D_H_
#define _VISCBNDRY2D_H_

//
// $Id: ViscBndry2D.H,v 1.4 1998-03-26 06:41:01 almgren Exp $
//

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry2D : public MCInterpBndryData
{
public:
    ViscBndry2D() : MCInterpBndryData()  {}

    ViscBndry2D(const BoxArray& _grids, const Geometry& geom) :
        MCInterpBndryData(_grids,4,geom) {};

#if 0
    virtual void setBndryConds(const Array<BCRec>& phys_bc, IntVect& ratio);

    void setHomogValues(const Array<BCRec>& phys_bc, IntVect& ratio);
#else
    virtual void setBndryConds(const Array<BCRec>& phys_bc, int ratio);

    void setHomogValues(const Array<BCRec>& phys_bc, int ratio);
#endif

};

#endif
