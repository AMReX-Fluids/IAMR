#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

//
// $Id: ViscBndry.H,v 1.5 1998-05-29 23:20:30 lijewski Exp $
//

#include <InterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry
    :
    public InterpBndryData
{
public:
    ViscBndry()
        :
        InterpBndryData()  {}

    ViscBndry (const BoxArray& _grids,
               int             _ncomp,
               const Geometry& _geom)
        :
        InterpBndryData(_grids,_ncomp,_geom)
    {}

    virtual void setBndryConds (const BCRec&   phys_bc,
                                /*const*/ IntVect& ratio);

    void setHomogValues (const BCRec& bc,
                         /*const*/ IntVect& ratio);
};

#endif /*_VISCBNDRYDATA_H_*/
