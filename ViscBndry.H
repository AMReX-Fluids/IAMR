#ifndef _VISCBNDRYDATA_H_
#define _VISCBNDRYDATA_H_

#include <InterpBndryData.H>
#include <BC_TYPES.H>

class ViscBndry : public InterpBndryData
{
public:
    ViscBndry() : InterpBndryData()  {}

    ViscBndry(const BoxArray& _grids, int _ncomp,const Geometry& geom) :
        InterpBndryData(_grids,_ncomp,geom) {};

    virtual void setBndryConds(const BCRec& phys_bc,
                               const Geometry& geom, IntVect& ratio);

    void setHomogValues(const BCRec& bc, IntVect& ratio);
};

#endif
