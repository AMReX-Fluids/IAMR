//BL_COPYRIGHT_NOTICE

#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.6 1999-01-04 18:15:17 marc Exp $
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
