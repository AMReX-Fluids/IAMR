#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

//
// $Id: LevelBld.H,v 1.3 1997-10-08 20:15:34 car Exp $
//

#include <REAL.H>
#include <Box.H>
#include <Geometry.H>
#include <Amr.H>
#include <AmrLevel.H>

#ifndef TESTINGMG
extern LevelBld* getLevelBld();
#endif

class LevelBld
{
public:
    virtual void variableSetUp() = 0;
    virtual void variableCleanUp() = 0;
    virtual AmrLevel* operator()()  = 0;
    virtual AmrLevel* operator()(Amr &papa, int lev,
				 const Geometry &geom_lev,
				 const BoxArray &ba, Real time) = 0;
};

#endif
