
//
// $Id: NSBld.cpp,v 1.9 2011-08-08 19:46:32 lijewski Exp $
//

#include <winstd.H>

#include <NavierStokes.H>

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
