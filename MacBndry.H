#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.3 1997-09-26 16:57:01 lijewski Exp $
//

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
