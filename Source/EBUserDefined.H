#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& /*geom*/,
              const int /*required_coarsening_level*/,
              const int /*max_coarsening_level*/)
{
    // ParmParse your geometry parameters

    // Build geometry pieces using EB2::* methods

    // Build your geometry shop using EB2::makeShop

    // Build geom using EB2::Build

    // We shoulnd't be here, copy this file in your run folder
    // and implement your geometry
    Abort("In default EBUserDefined function! Shouldn't be here. Copy and edit this file for your needs");
}
#endif
#endif
