#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

//
// $Id: LevelBld.H,v 1.2 1997-09-26 16:57:01 lijewski Exp $
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
				 const BoxArray &ba, REAL time) = 0;
};

#endif
