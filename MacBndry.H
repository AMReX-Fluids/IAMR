//BL_COPYRIGHT_NOTICE

#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.5 1998-09-29 20:08:28 lijewski Exp $
//

#include <InterpBndryData.H>

class MacBndry : public InterpBndryData
{
public:
    MacBndry() : InterpBndryData()  {}

    MacBndry(const BoxArray& _grids, int _ncomp, const Geometry& _geom) :
        InterpBndryData(_grids,_ncomp,_geom) {};

    virtual void setBndryConds(const BCRec& phys_bc, IntVect& ratio);
};

#endif
