//BL_COPYRIGHT_NOTICE

#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

//
// $Id: ViscBndry.H,v 1.7 1999-01-04 18:15:17 marc Exp $
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
                                /*const*/ IntVect& ratio,
				int comp=0);

    void setHomogValues (const BCRec& bc,
                         /*const*/ IntVect& ratio);
};

#endif /*_VISCBNDRYDATA_H_*/
