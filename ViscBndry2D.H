#ifndef _VISCBNDRY2D_H_
#define _VISCBNDRY2D_H_  "%W% %G%"

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry2D : public MCInterpBndryData
{
public:
    ViscBndry2D() : MCInterpBndryData()  {}

    ViscBndry2D(const BoxArray& _grids, const Geometry& geom) :
        MCInterpBndryData(_grids,4,geom) {};

#if 0
    virtual void setBndryConds(const Array<BCRec>& phys_bc,
                               const Geometry& geom, IntVect& ratio);

    void setHomogValues(const Array<BCRec>& phys_bc, IntVect& ratio);
#else
    virtual void setBndryConds(const Array<BCRec>& phys_bc,
                               const Geometry& geom, int ratio);

    void setHomogValues(const Array<BCRec>& phys_bc, int ratio);
#endif

};

#endif
