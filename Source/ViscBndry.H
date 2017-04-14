
#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

#include <AMReX_InterpBndryData.H>
#include <AMReX_BC_TYPES.H>

class ViscBndry
    :
    public InterpBndryData
{
public:
    ViscBndry ();

    ViscBndry (const BoxArray& _grids,
               int             _ncomp,
               const Geometry& _geom);

    virtual ~ViscBndry () override;

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio) override;

    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
				int            comp=0) override;

    void setHomogValues (const BCRec&   bc,
                         const IntVect& ratio);
};

#endif /*_VISCBNDRYDATA_H_*/
