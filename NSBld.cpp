#define _NS_BLD_C_  $Id: NSBld.cpp,v 1.2 1997-09-19 20:08:25 car Exp $

#include <NavierStokes.H>
#include <ParmParse.H>

// --------------------------------------------------------------------
// -----   NSBld class instantiation
// --------------------------------------------------------------------

NSBld nsbld;

LevelBld* getLevelBld()
{
    return &nsbld;
}

// --------------------------------------------------------------------
// -----   NSBld class implementation
// --------------------------------------------------------------------

void
NSBld::variableSetUp()
{
    NavierStokes::variableSetUp();
}

void
NSBld::variableCleanUp()
{
    NavierStokes::variableCleanUp();
}

AmrLevel*
NSBld::operator()()
{
    return new NavierStokes;
}

AmrLevel*
NSBld::operator()(Amr &papa, int lev, const Geometry &level_geom,
		  const BoxArray &ba, REAL time)
{
    return new NavierStokes(papa, lev, level_geom, ba, time);
}
