
#ifndef _MACBNDRYDATA_H_
#define _MACBNDRYDATA_H_

//
// $Id: MacBndry.H,v 1.9 2006-09-13 18:03:02 lijewski Exp $
//

#include <InterpBndryData.H>

class MacBndry
    :
    public InterpBndryData
{
public:

    MacBndry ();

    MacBndry (const BoxArray& _grids,
              int             _ncomp,
              const Geometry& _geom);

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio);

    virtual void setBndryConds (const BCRec& phys_bc,
                                IntVect&     ratio,
                                int          comp=0);
};

#endif
