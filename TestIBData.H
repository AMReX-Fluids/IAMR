//
// $Id: TestIBData.H,v 1.1 1997-07-08 23:08:02 vince Exp $
//
#ifndef _TestIBData_H_
#define _TestIBData_H_

#include <InterpBndryData.H>
#include <BC_TYPES.H>

class TestIBData : public InterpBndryData
{
public:
    TestIBData() : InterpBndryData()  {}

    TestIBData(const BoxArray& _grids, int _ncomp,const Geometry& geom) :
        InterpBndryData(_grids,_ncomp,geom) {};

    void setBndryConds(const BCRec& bc,
		       const Geometry& geom, IntVect& ratio);
};

#endif
