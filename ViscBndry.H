#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

//
// $Id: ViscBndry.H,v 1.4 1998-03-26 06:41:01 almgren Exp $
//

#include <InterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry : public InterpBndryData
{
public:
    ViscBndry() : InterpBndryData()  {}

    ViscBndry(const BoxArray& _grids, int _ncomp,const Geometry& _geom) :
        InterpBndryData(_grids,_ncomp,_geom) {};

    virtual void setBndryConds(const BCRec& phys_bc, IntVect& ratio);

    void setHomogValues(const BCRec& bc, IntVect& ratio);
};

#endif
