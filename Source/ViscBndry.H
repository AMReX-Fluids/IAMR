
#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

#include <InterpBndryData.H>
#include <BC_TYPES.H>
#include <BL_CXX11.H>

class ViscBndry
    :
    public InterpBndryData
{
public:
    ViscBndry ();

    ViscBndry (const BoxArray& _grids,
               int             _ncomp,
               const Geometry& _geom);

    virtual ~ViscBndry () BL_OVERRIDE;

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio) BL_OVERRIDE;

    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
				int            comp=0) BL_OVERRIDE;

    void setHomogValues (const BCRec&   bc,
                         const IntVect& ratio);
};

#endif /*_VISCBNDRYDATA_H_*/
