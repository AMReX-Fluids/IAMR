
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
