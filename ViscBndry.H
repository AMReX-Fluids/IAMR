#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_  $Id: ViscBndry.H,v 1.1 1997-07-08 23:08:15 vince Exp $

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
