
#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.7 2000-10-02 20:50:14 lijewski Exp $
//

#include <InterpBndryData.H>

class MacBndry : public InterpBndryData
{
public:
    MacBndry() : InterpBndryData()  {}

    MacBndry(const BoxArray& _grids, int _ncomp, const Geometry& _geom) :
        InterpBndryData(_grids,_ncomp,_geom) {};

    virtual void setBndryConds(const BCRec& phys_bc, IntVect& ratio, int comp=0);
};

#endif
