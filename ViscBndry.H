#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

//
// $Id: ViscBndry.H,v 1.3 1997-09-26 16:57:07 lijewski Exp $
//

#include <InterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry : public InterpBndryData
{
public:
    ViscBndry() : InterpBndryData()  {}

    ViscBndry(const BoxArray& _grids, int _ncomp,const Geometry& geom) :
        InterpBndryData(_grids,_ncomp,geom) {};

    virtual void setBndryConds(const BCRec& phys_bc,
                               const Geometry& geom, IntVect& ratio);

    void setHomogValues(const BCRec& bc, IntVect& ratio);
};

#endif
