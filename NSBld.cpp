
//
// $Id: NSBld.cpp,v 1.8 2004-07-27 17:28:57 car Exp $
//

#include <winstd.H>

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
                  const BoxArray &ba, Real time)
{
    return new NavierStokes(papa, lev, level_geom, ba, time);
}
