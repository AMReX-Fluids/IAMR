#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_ $Id: MacBndry.H,v 1.1 1997-07-08 23:08:14 vince Exp $

#include <InterpBndryData.H>

class MacBndry : public InterpBndryData
{
public:
    MacBndry() : InterpBndryData()  {}

    MacBndry(const BoxArray& _grids, int _ncomp, const Geometry& geom) :
        InterpBndryData(_grids,_ncomp,geom) {};

    virtual void setBndryConds(const BCRec& phys_bc,
                               const Geometry& geom, IntVect& ratio);
};

#endif
